/* app.js -- table rendering, filtering, and the 3Dmol mode viewer.
 *
 * The server sends one JSON payload with the scored modes and their
 * displacement vectors; everything below is client-side, which is why the
 * viewer behaves identically on localhost and on a serverless host.
 *
 * The viewer is fed OUR bond list, not 3Dmol's distance-based guess, so the
 * picture is drawn from exactly the connectivity s[V_S] was computed from.
 */
(function () {
  "use strict";

  var P = JSON.parse(document.getElementById("payload").textContent);
  var KEYS = ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "V_S"];
  var tbody = document.querySelector("#tbl tbody");
  var selected = null;
  var viewer = null, model = null;
  var framed = false;          // camera framed once, then left to the user

  // ---------------------------------------------------------------- table
  // The T/R toggle is only rendered when references were constructed (the
  // 3N-6 path). On a 3N upload the element does not exist -- nothing was
  // constructed -- so treat it as off rather than dereferencing null.
  function refsOn() {
    var el = document.getElementById("f-ref");
    return !!el && el.checked;
  }

  function rows() {
    var out = P.vibrations.slice();
    if (refsOn()) out = P.references.concat(out);
    return out;
  }

  function visible(r) {
    // Constructed T/R modes carry a synthetic frequency of 0, so the frequency
    // range says nothing about them -- applying it would silently hide them
    // whenever the low bound is above 0, which is every real molecule. Their
    // visibility is governed solely by the "T/R reference modes" checkbox.
    if (!r.is_reference && r.frequency !== null) {
      var lo = parseFloat(document.getElementById("f-lo").value);
      var hi = parseFloat(document.getElementById("f-hi").value);
      if (!isNaN(lo) && r.frequency < lo) return false;
      if (!isNaN(hi) && r.frequency > hi) return false;
    }
    if (["S", "B", "SB"].indexOf(r.label) !== -1) {
      return document.getElementById("f-" + r.label).checked;
    }
    return true;                       // external labels are not label-filtered
  }

  function fmt(v) {
    if (v === null || v === undefined) return "—";
    return (v >= 0 ? " " : "") + v.toFixed(3);
  }

  function render() {
    tbody.innerHTML = "";
    var shown = 0;
    rows().forEach(function (r) {
      if (!visible(r)) return;
      shown++;
      var tr = document.createElement("tr");
      tr.dataset.index = r.index;
      if (r.is_reference) tr.classList.add("ref");
      if (selected === r.index) tr.classList.add("sel");

      var cells = '<td class="c"><input type="checkbox" class="pick" checked></td>' +
                  '<td class="mode">' + r.name + "</td>" +
                  // A reference mode has no frequency -- it is a constructed
                  // ideal motion, not a solution of the Hessian -- so show a
                  // dash rather than a misleading 0.00 cm-1.
                  '<td class="n' + (r.is_reference ? " dim" : "") + '">' +
                  (r.is_reference || r.frequency === null
                     ? "—" : r.frequency.toFixed(2)) + "</td>";
      KEYS.forEach(function (k) {
        var cls = "n";
        // r.highlight is the dominant score on a constructed T/R reference and
        // V_S on a real mode -- see pipeline._row. On a real mode the six T/R
        // numbers are dimmed: they are diagnostics, not the mode's assigned
        // character, and only these six reference modes carry T/R character.
        // r.highlights follows the LABEL: one cell for a clean external or an
        // internal mode, two for a mixed external (V_S + the largest |T/R|),
        // which by definition is part vibration and part rigid-body motion.
        if (r.highlights.indexOf(k) !== -1) cls += " hot";
        else if (!r.is_reference && k !== "V_S") cls += " dim";
        cells += '<td class="' + cls + '">' + fmt(r.scores[k]) + "</td>";
      });
      // A starred external carries its vibrational character as
      // annotation="vibration=SB". Show that as a second chip rather than
      // spelling it out: [Tx*] [SB].
      var ann = "";
      if (r.annotation) {
        var vm = /^vibration=(\w+)$/.exec(r.annotation);
        ann = vm ? '<b class="lab ' + vm[1] + '">' + vm[1] + "</b>"
                 : '<span class="ann">' + r.annotation + "</span>";
      }
      cells += '<td><b class="lab ' + r.label.replace("*", "x") + '">' + r.label + "</b>" +
               ann + "</td>" +
               "<td>" + (r.irrep || "—") + "</td>";
      tr.innerHTML = cells;

      tr.addEventListener("click", function (e) {
        if (e.target.classList.contains("pick")) return;
        show(r);
      });
      tbody.appendChild(tr);
    });
    document.getElementById("count").textContent =
      shown + " of " + rows().length + " shown";
  }

  // --------------------------------------------------------------- viewer
  function elementColor() { return null; }   // let 3Dmol use its element palette

  function show(r) {
    selected = r.index;
    document.querySelectorAll("#tbl tbody tr").forEach(function (tr) {
      tr.classList.toggle("sel", parseInt(tr.dataset.index, 10) === r.index);
    });
    document.getElementById("vtitle").textContent =
      r.name + (r.frequency !== null && r.frequency !== 0
                ? "  ·  " + r.frequency.toFixed(1) + " cm⁻¹" : "");
    var lb = document.getElementById("vlabel");
    lb.textContent = r.label + " — " + r.label_text;
    lb.className = "lab " + r.label.replace("*", "x");

    // per-bond s_AB panel
    var bl = document.getElementById("bondlist");
    if (!r.bonds.length) { bl.innerHTML = '<span class="muted">—</span>'; }
    else {
      var max = Math.max.apply(null, r.bonds.map(function (b) { return Math.abs(b.s_AB); })) || 1;
      bl.innerHTML = r.bonds.map(function (b) {
        var w = (Math.abs(b.s_AB) / max * 100).toFixed(0);
        return '<div class="brow"><span class="bp">' + b.pair + "</span>" +
               '<span class="bbar"><i class="' + (b.s_AB < 0 ? "neg" : "pos") +
               '" style="width:' + w + '%"></i></span>' +
               '<span class="bv">' + (b.s_AB >= 0 ? "+" : "") + b.s_AB.toFixed(3) + "</span></div>";
      }).join("");
    }
    draw(r);
  }

  /* Extended XYZ: "elem x y z dx dy dz". 3Dmol's xyz parser reads columns
   * 4-6 as dx/dy/dz when a line has >= 7 fields, which is exactly what
   * vibrate() consumes -- so the displacement vectors go in with no
   * conversion. Building a model this way (rather than addModel() with no
   * data + addAtoms) is what the parser actually supports; the empty-model
   * route leaves this.atoms undefined and throws inside addMolData. */
  function extendedXyz(r) {
    var lines = [String(P.atoms.length), P.title + " " + r.name];
    for (var i = 0; i < P.atoms.length; i++) {
      var g = P.geometry[i], d = r.vector[i];
      lines.push(P.atoms[i] + " " +
                 g[0].toFixed(6) + " " + g[1].toFixed(6) + " " + g[2].toFixed(6) + " " +
                 d[0].toFixed(6) + " " + d[1].toFixed(6) + " " + d[2].toFixed(6));
    }
    return lines.join("\n") + "\n";
  }

  /* Principal axes of inertia. The molecule is rotated into this frame before
   * scoring, so Tx / Rz and friends are defined against exactly these
   * directions -- drawing them makes the score columns readable. Colours avoid
   * red, which is the displacement arrows. */
  var AXES = [
    { k: "x", v: [1, 0, 0], c: "#f59e0b" },
    { k: "y", v: [0, 1, 0], c: "#3b82f6" },
    { k: "z", v: [0, 0, 1], c: "#10b981" }
  ];

  function drawAxes() {
    viewer.removeAllLabels();
    if (!document.getElementById("axes").checked) return;
    var span = 1.0;
    P.geometry.forEach(function (g) {
      span = Math.max(span, Math.abs(g[0]), Math.abs(g[1]), Math.abs(g[2]));
    });
    var L = span * 1.35 + 0.6;
    AXES.forEach(function (a) {
      viewer.addArrow({
        start: { x: 0, y: 0, z: 0 },
        end: { x: a.v[0] * L, y: a.v[1] * L, z: a.v[2] * L },
        radius: 0.035, radiusRatio: 2.2, mid: 0.88, color: a.c
      });
      viewer.addLabel(a.k, {
        position: { x: a.v[0] * (L + 0.28), y: a.v[1] * (L + 0.28), z: a.v[2] * (L + 0.28) },
        fontSize: 12, fontColor: a.c, backgroundOpacity: 0.0,
        showBackground: false, inFront: true
      });
    });
  }

  function draw(r) {
    if (!viewer) return;

    // Stop the previous loop FIRST. animate() registers timers in
    // viewer.animationTimers and increments an internal counter; neither
    // removeAllModels() nor addModel() clears them. Without this, every row
    // click and every slider nudge stacks another loop, and each one calls
    // setFrame() on a different index into the same canvas -- which reads as
    // the atoms flickering.
    viewer.stopAnimate();
    viewer.removeAllModels();
    viewer.removeAllShapes();

    // assignBonds:false -- 3Dmol would otherwise guess bonds by distance, and
    // the picture could then show a different connectivity than the one
    // s[V_S] was computed from. We supply the real bond list below.
    model = viewer.addModel(extendedXyz(r), "xyz", { assignBonds: false });

    var atoms = model.selectedAtoms({});
    // drawBondSticks draws each bond once, from the lower-index atom to the
    // higher: `if (atom.index < partner.index)`. 3Dmol's xyz parser sets
    // `serial` but leaves `index` null, and `null < null` is false -- so every
    // bond was skipped and no sticks were ever drawn. Assign it explicitly.
    atoms.forEach(function (a, i) { a.index = i; });
    P.bonds.forEach(function (b) {
      atoms[b[0]].bonds.push(b[1]); atoms[b[0]].bondOrder.push(1);
      atoms[b[1]].bonds.push(b[0]); atoms[b[1]].bondOrder.push(1);
    });

    var amp = parseFloat(document.getElementById("amp").value);
    var frames = parseInt(document.getElementById("frm").value, 10);
    var arrows = document.getElementById("arrows").checked;
    // Red reads clearly against both the light and dark viewer backgrounds,
    // and against the CPK element colours (no CPK element is red except
    // oxygen, which is a darker #ff0d0d-family tone at sphere scale).
    var arrowSpec = arrows
      ? { color: "#e03131", radius: 0.06, radiusRatio: 2.0, mid: 0.75 }
      : undefined;

    // Builds the animation frames from dx/dy/dz. bothWays=true lays the frames
    // out as -amp -> 0 -> +amp, which is why the loop below is backAndForth:
    // a plain forward loop would snap from +amp straight back to -amp.
    model.vibrate(frames, amp, true, arrowSpec);

    // Style AFTER vibrate(): vibrate rebuilds the model's frames, and a style
    // set beforehand is applied to geometry it then replaces. Sticks are drawn
    // from atom.bonds, which is OUR connectivity (assignBonds:false above), so
    // the sticks and s[V_S] describe the same bonds.
    var showBonds = document.getElementById("bonds").checked;
    model.setStyle({}, showBonds
      ? { stick: { radius: 0.15 }, sphere: { scale: 0.25 } }
      : { sphere: { scale: 0.32 } });   // no sticks -> larger spheres to read

    drawAxes();

    // Frame the molecule once. Re-zooming on every mode change fights the
    // user's own rotate/zoom and looks like the view jumping.
    if (!framed) {
      viewer.zoomTo();
      framed = true;
    }

    if (document.getElementById("play").checked) {
      viewer.animate({ loop: "backAndForth", interval: 90 });
    } else {
      // setFrame is async (animate() itself chains off its promise), so render
      // only once the frame is actually in place.
      var done = viewer.setFrame(frames);     // frames = midpoint = equilibrium
      if (done && done.then) { done.then(function () { viewer.render(); }); }
      else { viewer.render(); }
    }
  }

  // -------------------------------------------------------------- exports
  var DL = JSON.parse(document.getElementById("downloads").textContent);

  function download(key, name, mime) {
    var blob = new Blob([DL[key]], { type: mime });
    var a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = name;
    document.body.appendChild(a); a.click();
    setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 0);
  }

  // ------------------------------------------------------------- wire-up
  ["f-S", "f-B", "f-SB", "f-ref"].forEach(function (id) {
    var el = document.getElementById(id);
    if (el) el.addEventListener("change", render);   // f-ref absent under 3N
  });
  ["f-lo", "f-hi"].forEach(function (id) {
    document.getElementById(id).addEventListener("input", render);
  });
  document.getElementById("chk-all").addEventListener("change", function () {
    var on = this.checked;
    document.querySelectorAll(".pick").forEach(function (c) { c.checked = on; });
  });
  document.getElementById("all").addEventListener("click", function () {
    ["f-S", "f-B", "f-SB"].forEach(function (i) { document.getElementById(i).checked = true; });
    document.getElementById("f-lo").value = Math.floor(P.freq_range[0]);
    document.getElementById("f-hi").value = Math.ceil(P.freq_range[1]);
    render();
  });
  document.getElementById("none").addEventListener("click", function () {
    ["f-S", "f-B", "f-SB"].forEach(function (i) { document.getElementById(i).checked = false; });
    render();
  });
  ["amp", "frm", "arrows", "play", "bonds", "axes"].forEach(function (id) {
    document.getElementById(id).addEventListener("change", function () {
      if (selected !== null) {
        var r = P.references.concat(P.vibrations).find(function (x) { return x.index === selected; });
        if (r) draw(r);
      }
    });
  });
  document.querySelectorAll("[data-dl]").forEach(function (b) {
    b.addEventListener("click", function () {
      var k = b.dataset.dl;
      if (k === "csv") download("csv", P.title + "_scores.csv", "text/csv");
      if (k === "csvall") download("csvall", P.title + "_scores_all.csv", "text/csv");
      if (k === "vsc") download("vsc", P.title + ".vsc", "text/plain");
    });
  });

  // keyboard: up/down move through visible rows
  document.addEventListener("keydown", function (e) {
    if (e.key !== "ArrowDown" && e.key !== "ArrowUp") return;
    if (/INPUT|TEXTAREA/.test(document.activeElement.tagName)) return;
    var trs = Array.prototype.slice.call(document.querySelectorAll("#tbl tbody tr"));
    if (!trs.length) return;
    var cur = trs.findIndex(function (t) { return t.classList.contains("sel"); });
    var next = e.key === "ArrowDown" ? Math.min(cur + 1, trs.length - 1) : Math.max(cur - 1, 0);
    if (cur === -1) next = 0;
    e.preventDefault();
    trs[next].click();
    trs[next].scrollIntoView({ block: "nearest" });
  });

  // ------------------------------------------------------------- startup
  render();
  try {
    viewer = $3Dmol.createViewer(document.getElementById("mol"), {
      backgroundColor: getComputedStyle(document.body)
        .getPropertyValue("--viewer-bg").trim() || "white"
    });
    if (P.vibrations.length) show(P.vibrations[0]);
  } catch (err) {
    document.getElementById("mol").innerHTML =
      '<p class="muted" style="padding:16px">3D viewer failed to start: ' + err + "</p>";
  }
})();
