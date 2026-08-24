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
  var picked = [];             // mode indices ticked for comparison, max 3
  var sortKey = null, sortDir = 0;   // 0 = file order, 1 = ascending, -1 = descending
  var selectedBond = null;           // index into P.bonds, highlighted both ways
  var BOND_HL = "#f2c14e";
  var cmpViews = [];           // one live 3Dmol viewer per compared mode
  var MAX_CMP = 3;
  var CMP_COLOURS = ["#0b5f52", "#8a4b1f", "#4a3f8f"];

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
    return sortRows(out);
  }

  /* Sort value for a column. Missing values are pushed to the end in BOTH
   * directions rather than sorting as zero, so a mode with no frequency does
   * not land among the low-frequency ones. */
  function sortVal(r, key) {
    if (key === "name") return r.index;          // "Vib 10" must follow "Vib 9"
    // A reference mode is displayed with a dash, not 0.00 -- it is constructed,
    // not a solution of the Hessian -- so it must sort as missing too, or it
    // lands among the lowest frequencies.
    if (key === "freq") return r.is_reference ? null : r.frequency;
    if (key === "label") return r.label;
    if (key === "irrep") return r.irrep;
    return r.scores[key];
  }

  function sortRows(list) {
    if (!sortKey || !sortDir) return list;
    var missing = [], present = [];
    list.forEach(function (r) {
      var v = sortVal(r, sortKey);
      (v === null || v === undefined || v === "" ? missing : present).push(r);
    });
    present.sort(function (a, b) {
      var x = sortVal(a, sortKey), y = sortVal(b, sortKey), c;
      if (typeof x === "string" || typeof y === "string") {
        c = String(x).localeCompare(String(y));
      } else {
        c = x - y;
      }
      // stable tie-break on file order, so equal scores keep a fixed order
      return (c || (a.index - b.index)) * sortDir;
    });
    return present.concat(missing);
  }

  function applySort(key) {
    if (sortKey !== key) { sortKey = key; sortDir = 1; }
    else if (sortDir === 1) { sortDir = -1; }
    else if (sortDir === -1) { sortKey = null; sortDir = 0; }   // back to file order
    else { sortDir = 1; }
    markHeaders();
    render();
  }

  function markHeaders() {
    document.querySelectorAll("#tbl thead th[data-sort]").forEach(function (th) {
      var k = th.dataset.sort;
      var on = (k === sortKey && sortDir !== 0);
      th.classList.toggle("sorted", on);
      var ind = th.querySelector(".ind");
      if (!ind) {
        ind = document.createElement("span");
        ind.className = "ind";
        th.appendChild(ind);
      }
      ind.textContent = on ? (sortDir === 1 ? " \u25b2" : " \u25bc") : "";
    });
  }

  /* ---- per-bond sorting -------------------------------------------------
   * ONE ordering drives both the side panel and the compare table. They list
   * the same bonds, and letting them disagree makes the two windows impossible
   * to read against each other -- which is the whole point of comparing.
   *
   * `mode` is the payload index of the row whose s_AB values order the list,
   * so clicking a mode column in the compare table orders the side panel by
   * that same mode rather than by whatever happens to be selected.
   */
  var bondSort = { key: null, dir: 0, mode: null };

  function naturalCmp(a, b) {
    // "C2-H3" must sort before "C10-H11", so compare digit runs numerically
    var x = String(a).match(/\d+|\D+/g) || [], y = String(b).match(/\d+|\D+/g) || [];
    for (var i = 0; i < Math.min(x.length, y.length); i++) {
      var p = x[i], q = y[i], d;
      if (/^\d/.test(p) && /^\d/.test(q)) {
        d = parseInt(p, 10) - parseInt(q, 10);
        if (d) return d;
      } else if (p !== q) { return p < q ? -1 : 1; }
    }
    return x.length - y.length;
  }

  function bondPair(row, i) {
    return (row && row.bonds && row.bonds[i] && row.bonds[i].pair) || P.bond_labels[i] || "";
  }

  /* Display order as ORIGINAL bond indices. Returning indices rather than
   * reordered bond objects is deliberate: data-b, selectedBond, P.bonds[] and
   * the viewer's click targets are all keyed on the original position, so a
   * sorted view must never renumber them. */
  function bondOrder(fallbackRow) {
    var ref = P.bonds.map(function (_, i) { return i; });
    if (!bondSort.key || !bondSort.dir) return ref;
    var vr = bondSort.mode === null ? fallbackRow : byIndex(bondSort.mode);
    if (bondSort.key === "val" && !(vr && vr.bonds && vr.bonds.length)) return ref;
    return ref.slice().sort(function (i, j) {
      var c = bondSort.key === "pair"
        ? naturalCmp(bondPair(fallbackRow, i), bondPair(fallbackRow, j))
        // signed, like the main table's score columns -- with the centred bars
        // this reads as a gradient from left-extending to right-extending
        : vr.bonds[i].s_AB - vr.bonds[j].s_AB;
      return (c || (i - j)) * bondSort.dir;   // stable tie-break on bond order
    });
  }

  function applyBondSort(key, mode) {
    mode = (mode === undefined ? null : mode);
    if (bondSort.key !== key || bondSort.mode !== mode) {
      bondSort.key = key; bondSort.mode = mode; bondSort.dir = 1;
    } else if (bondSort.dir === 1) {
      bondSort.dir = -1;
    } else {
      bondSort.key = null; bondSort.mode = null; bondSort.dir = 0;   // bond order
    }
    // Repaint the lists only. Going through show() would rebuild the model and
    // restart the animation, so sorting would visibly jolt the viewer.
    renderBondPanel(selected === null ? null : byIndex(selected));
    if (picked.length >= 2) renderCompare();
    markBondHeaders();
  }

  function markBondHeaders() {
    var on = bondSort.key && bondSort.dir;
    var arrow = bondSort.dir === 1 ? " \u25b2" : " \u25bc";
    document.querySelectorAll(".bsort").forEach(function (el) {
      var k = el.dataset.bsort;
      var m = el.dataset.mode === undefined || el.dataset.mode === ""
        ? null : parseInt(el.dataset.mode, 10);
      // The side panel's own s_AB header has no mode of its own: it is lit
      // whenever a value sort is active, whichever mode is driving it.
      var hit = on && k === bondSort.key &&
                (el.dataset.mode === undefined || el.dataset.mode === ""
                   ? true : m === bondSort.mode);
      el.classList.toggle("sorted", !!hit);
      var ind = el.querySelector(".ind");
      if (!ind) {
        ind = document.createElement("span");
        ind.className = "ind";
        el.appendChild(ind);
      }
      ind.textContent = hit ? arrow : "";
    });
    // Say which mode is driving the order when it is not the one on screen.
    var h = document.querySelector('.bhead .bsort[data-bsort="val"]');
    if (h) {
      var vr = bondSort.mode === null ? null : byIndex(bondSort.mode);
      h.title = (on && bondSort.key === "val" && vr && vr.index !== selected)
        ? "ordered by " + vr.name : "";
    }
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

      var on = picked.indexOf(r.index) !== -1;
      var cells = '<td class="c"><input type="checkbox" class="pick"' +
                  (on ? " checked" : "") + ' data-i="' + r.index + '"></td>' +
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
      var cb = tr.querySelector(".pick");
      if (cb) cb.addEventListener("change", function () { togglePick(r.index, cb); });
      tbody.appendChild(tr);
    });
    document.getElementById("count").textContent =
      shown + " of " + rows().length + " shown";
  }

  // ------------------------------------------------------------- compare
  function allRows() { return P.references.concat(P.vibrations); }

  function byIndex(i) {
    var all = allRows();
    for (var k = 0; k < all.length; k++) if (all[k].index === i) return all[k];
    return null;
  }

  function togglePick(i, cb) {
    var at = picked.indexOf(i);
    if (at !== -1) {
      picked.splice(at, 1);
    } else {
      if (picked.length >= MAX_CMP) {
        // At most three: refuse rather than silently dropping someone's
        // earlier choice, so the tick that did not take is visible.
        cb.checked = false;
        return;
      }
      picked.push(i);
    }
    renderCompare();
  }

  function swatch(k) {
    return '<span class="cswatch" style="background:' + CMP_COLOURS[k] + '"></span>';
  }

  function renderCompare() {
    var box = document.getElementById("compare");
    if (!box) return;
    document.getElementById("cmp-n").textContent = String(picked.length);
    if (picked.length < 2) {
      box.hidden = true;
      quietViews();          // stop the loops; keep the contexts for reuse
      syncTicks();
      return;
    }
    box.hidden = false;

    var rows = picked.map(byIndex).filter(Boolean);
    buildCompareViews(rows);

    // ---- the seven scores, side by side ----
    var h = "<table><tr><th>score</th>";
    rows.forEach(function (r, k) {
      h += "<th>" + swatch(k) + r.name + "</th>";
    });
    h += "</tr>";
    KEYS.forEach(function (key) {
      h += "<tr><td>" + key + "</td>";
      rows.forEach(function (r) {
        var hot = r.highlights.indexOf(key) !== -1 ? ' class="hot"' : "";
        h += "<td" + hot + ">" + fmt(r.scores[key]) + "</td>";
      });
      h += "</tr>";
    });
    h += "<tr><td>label</td>";
    rows.forEach(function (r) { h += "<td>" + r.label + "</td>"; });
    h += "</tr><tr><td>" + (P.is_emit ? "eigenvalue" : "ν / cm⁻¹") + "</td>";
    rows.forEach(function (r) {
      h += "<td>" + (r.frequency === null ? "—" : r.frequency.toFixed(2)) + "</td>";
    });
    h += "</tr></table>";
    document.getElementById("cmp-scores").innerHTML = h;

    // ---- per-bond s_AB ----
    // Every mode carries the same bond list in the same order (ModeScorer.bList),
    // so the rows line up without matching on names.
    var scale = 0;
    rows.forEach(function (r) {
      r.bonds.forEach(function (b) { scale = Math.max(scale, Math.abs(b.s_AB)); });
    });
    scale = scale || 1;

    // Each mode's header sorts the bonds by THAT mode's s_AB -- with three
    // modes side by side, "which bonds does this one move most" is the
    // question the panel exists to answer.
    var b = '<table><tr><th class="bsort" data-bsort="pair">bond</th>';
    rows.forEach(function (r, k) {
      b += '<th class="bsort" data-bsort="val" data-mode="' + r.index + '">' +
           swatch(k) + r.name + "</th>";
    });
    b += "</tr>";
    bondOrder(rows[0]).forEach(function (i) {
      b += '<tr class="cbrow' + (i === selectedBond ? " sel" : "") +
           '" data-b="' + i + '"><td>' + rows[0].bonds[i].pair + "</td>";
      rows.forEach(function (r) {
        var v = r.bonds[i] ? r.bonds[i].s_AB : 0;
        var w = (Math.abs(v) / scale * 50).toFixed(1);
        // The bar and number live in an inner span: making the <td> itself a
        // flex container takes it out of the table's column model and the mode
        // columns stack vertically instead of sitting side by side.
        b += '<td><span class="cv2"><span class="cbar"><i class="' +
             (v < 0 ? "neg" : "pos") + '" style="width:' + w + '%"></i></span>' +
             '<span class="num">' + (v >= 0 ? "+" : "") + v.toFixed(3) +
             "</span></span></td>";
      });
      b += "</tr>";
    });
    b += "</table>";
    var cb = document.getElementById("cmp-bonds");
    cb.innerHTML = b;
    cb.querySelectorAll("tr.cbrow").forEach(function (tr) {
      tr.addEventListener("click", function () {
        selectBond(parseInt(tr.dataset.b, 10));
      });
    });
    // headers are rebuilt with the table, so rebind every time
    cb.querySelectorAll("th.bsort").forEach(function (th) {
      th.addEventListener("click", function () {
        applyBondSort(th.dataset.bsort,
                      th.dataset.mode === undefined || th.dataset.mode === ""
                        ? null : parseInt(th.dataset.mode, 10));
      });
    });
    markBondHeaders();
    syncTicks();
  }

  /* One animated viewer per compared mode, running together -- seeing the modes
   * move side by side is the reason to compare; the tables are the detail.
   *
   * The viewers are created ONCE, up to MAX_CMP, and then reused. A browser
   * caps how many WebGL contexts may exist, and 3Dmol's GLViewer exposes no
   * teardown, so creating a fresh one on every tick would eventually exhaust
   * them and the viewers would quietly stop drawing. Spare slots are hidden,
   * not destroyed.
   */
  function ensureViews() {
    var host = document.getElementById("cmp-views");
    if (!host || typeof $3Dmol === "undefined") return false;
    if (cmpViews.length) return true;

    var bg = getComputedStyle(document.body)
               .getPropertyValue("--viewer-bg").trim() || "white";
    for (var k = 0; k < MAX_CMP; k++) {
      var card = document.createElement("div");
      card.className = "cview";
      card.innerHTML = '<div class="vh"></div><div class="cv"></div>';
      host.appendChild(card);
      var box = card.querySelector(".cv");
      try {
        cmpViews.push({ viewer: $3Dmol.createViewer(box, { backgroundColor: bg }),
                        card: card });
      } catch (e) {
        box.innerHTML = '<p class="muted" style="padding:14px">' + e + "</p>";
        return false;
      }
    }

    // Tie the cameras together: rotating, zooming or panning one panel does the
    // same to the others, so the modes stay in a common orientation.
    //
    // Deliberately NOT 3Dmol's linkViewer(). That propagates from show(), which
    // runs on every render -- including every animation frame -- so three
    // linked viewers push a full setView+render into each other 90ms, turning
    // 3 renders per tick into 9 and continuously overwriting each other's
    // cameras. Syncing on actual interaction costs nothing while idle.
    cmpViews.forEach(function (c, k) {
      var box = c.card.querySelector(".cv");
      if (!box) return;
      ["mousemove", "touchmove", "wheel", "mouseup", "touchend"].forEach(
        function (ev) {
          box.addEventListener(ev, function (e) {
            // a bare mousemove is not an interaction; only a drag is
            if (ev === "mousemove" && !e.buttons) return;
            scheduleSync(k);
          }, { passive: true });
        });
    });
    return true;
  }

  var syncPending = false;

  // rAF as a bare global, not off window: it is a global in every browser, and
  // reaching through window makes the lookup fail wherever window is not the
  // global object.
  var raf = (typeof requestAnimationFrame === "function")
    ? requestAnimationFrame
    : function (f) { return setTimeout(f, 16); };

  function scheduleSync(k) {
    if (syncPending) return;                       // coalesce to one per frame
    // Set the flag BEFORE scheduling: assigning the handle afterwards leaves it
    // stuck if the callback happens to run synchronously, and every later sync
    // is then dropped.
    syncPending = true;
    raf(function () {
      syncPending = false;
      syncFrom(k);
    });
  }

  function syncFrom(k) {
    var src = cmpViews[k];
    if (!src) return;
    var view;
    try { view = src.viewer.getView(); } catch (e) { return; }
    cmpViews.forEach(function (c, i) {
      if (i === k || c.card.hidden) return;
      // second argument suppresses re-propagation inside show()
      try { c.viewer.setView(view, true); } catch (e) {}
    });
  }

  function quietViews() {
    cmpViews.forEach(function (c) {
      try { c.viewer.stopAnimate(); } catch (e) {}
      c.card.hidden = true;
    });
  }

  function buildCompareViews(rows) {
    if (!ensureViews()) return;
    var host = document.getElementById("cmp-views");
    host.className = "cviews n" + rows.length;
    quietViews();

    var o = viewOpts();
    rows.forEach(function (r, k) {
      var c = cmpViews[k];
      if (!c) return;
      c.card.hidden = false;
      var freq = r.frequency === null ? "\u2014" : r.frequency.toFixed(1);
      c.card.querySelector(".vh").innerHTML =
        swatch(k) + "<b>" + r.name + "</b>" +
        '<span class="muted">' + freq + "</span>" +
        '<span class="sp"></span><b class="lab ' + r.label.replace("*", "x") +
        '">' + r.label + "</b>";
      // The canvas was sized while the card was hidden (display:none -> 0x0),
      // and the grid width changes with the number of modes shown, so the
      // viewer must be told its box again or it draws into a stale buffer --
      // which is why the third panel came up blank.
      try { c.viewer.resize(); } catch (e) {}
      c.model = renderMode(c.viewer, r, o);
      drawBondPicks(c.viewer);
      c.viewer.zoomTo();
      if (o.play) c.viewer.animate({ loop: "backAndForth", interval: 90 });
      else c.viewer.render();
    });
  }

  function syncTicks() {
    document.querySelectorAll("#tbl tbody .pick").forEach(function (cb) {
      cb.checked = picked.indexOf(parseInt(cb.dataset.i, 10)) !== -1;
    });
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

    renderBondPanel(r);
    draw(r);
  }

  /* The per-bond s_AB panel. Separate from show() so a sort can repaint the
   * list without touching the model or restarting the animation. */
  function renderBondPanel(r) {
    var bl = document.getElementById("bondlist");
    if (!bl) return;
    if (!r || !r.bonds.length) {
      bl.innerHTML = '<span class="muted">—</span>';
      return;
    }
    var max = Math.max.apply(null, r.bonds.map(function (b) { return Math.abs(b.s_AB); })) || 1;
    bl.innerHTML = bondOrder(r).map(function (bi) {
      var b = r.bonds[bi];
      // Half-width: the bar grows out from a centre line at zero, so the full
      // track spans -max..+max (matches the compare panel).
      var w = (Math.abs(b.s_AB) / max * 50).toFixed(1);
      return '<div class="brow' + (bi === selectedBond ? " sel" : "") +
             '" data-b="' + bi + '"><span class="bp">' + b.pair + "</span>" +
             '<span class="bbar"><i class="' + (b.s_AB < 0 ? "neg" : "pos") +
             '" style="width:' + w + '%"></i></span>' +
             '<span class="bv">' + (b.s_AB >= 0 ? "+" : "") + b.s_AB.toFixed(3) + "</span></div>";
    }).join("");
    bl.querySelectorAll(".brow").forEach(function (row) {
      row.addEventListener("click", function () {
        selectBond(parseInt(row.dataset.b, 10));
      });
    });
    scrollBondIntoView(bl);
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

  /* One invisible cylinder per bond: purely a click target, because 3Dmol's
   * sticks belong to the model and are not individually pickable. The highlight
   * itself is done in the model (see renderMode) so it follows the vibration --
   * a drawn cylinder would sit at the equilibrium geometry while the bond moved
   * through it. */
  function drawBondPicks(v) {
    v = v || viewer;
    if (!v) return;
    P.bonds.forEach(function (b, bi) {
      var a = P.geometry[b[0]], c = P.geometry[b[1]];
      v.addCylinder({
        start: { x: a[0], y: a[1], z: a[2] },
        end: { x: c[0], y: c[1], z: c[2] },
        radius: 0.26,
        color: "#ffffff",
        opacity: 0.0,          // invisible, but still picked
        clickable: true,
        callback: function () { selectBond(bi); }
      });
    });
  }

  /* Bring the selected row into view inside the bond list, which scrolls in
   * its own box.
   *
   * Measured with getBoundingClientRect rather than offsetTop: offsetTop is
   * relative to the nearest POSITIONED ancestor, and .viewer is
   * position:sticky, so it counted the panel header and the whole 3D canvas
   * too -- every row looked far below the window and the list slammed to the
   * bottom every time. Rects are independent of the positioning context.
   *
   * scrollTop is nudged by the overflow only, so a row already visible does
   * not move, and scrollIntoView() is avoided because it scrolls the page as
   * well -- jarring when the click came from the 3D view.
   */
  function scrollBondIntoView(bl) {
    if (selectedBond === null || !bl || !bl.getBoundingClientRect) return;
    var row = bl.querySelector('.brow[data-b="' + selectedBond + '"]');
    if (!row || !row.getBoundingClientRect) return;
    var r = row.getBoundingClientRect(), box = bl.getBoundingClientRect();
    if (r.top < box.top) {
      bl.scrollTop -= (box.top - r.top);
    } else if (r.bottom > box.bottom) {
      bl.scrollTop += (r.bottom - box.bottom);
    }
  }

  function selectBond(bi) {
    selectedBond = (selectedBond === bi) ? null : bi;   // click again to clear
    // Lists and colour only -- show()/renderCompare() would rebuild the models
    // and restart every animation, which is not what picking a bond asked for.
    renderBondPanel(selected === null ? null : byIndex(selected));
    markBondRows();
    recolorBond(viewer, model);
    cmpViews.forEach(function (c) {
      if (!c.card.hidden) recolorBond(c.viewer, c.model);
    });
  }

  /* The compare table's own selected-bond row. renderBondPanel() rebuilds the
   * side panel wholesale; this table is far bigger, and rebuilding it would
   * drop the scroll position, so only the class is moved. */
  function markBondRows() {
    document.querySelectorAll("#cmp-bonds tr.cbrow").forEach(function (tr) {
      tr.classList.toggle("sel", parseInt(tr.dataset.b, 10) === selectedBond);
    });
  }

  function drawOverlays() {
    // removeAllLabels() is global to the viewer, so axis letters and atom
    // numbers have to be (re)drawn together or one wipes the other.
    viewer.removeAllLabels();
    drawAxes();
    drawAtomNumbers();
  }

  /* Element + 1-based index, matching the bond labels (C1, H2, ...). Drawn at
   * the equilibrium geometry: labels are viewer-level, not part of the model,
   * so they cannot ride the vibration frames. The atom oscillates about this
   * position, so the number still reads as belonging to it. */
  function drawAtomNumbers() {
    var el = document.getElementById("labels");
    if (!el || !el.checked) return;
    P.atoms.forEach(function (sym, i) {
      var g = P.geometry[i];
      viewer.addLabel(sym + (i + 1), {
        position: { x: g[0], y: g[1], z: g[2] },
        fontSize: 11, fontColor: "#cfcfc8",
        backgroundOpacity: 0.45, backgroundColor: "#101114",
        borderThickness: 0, inFront: true
      });
    });
  }

  function drawAxes() {
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

  /* Paint the current selectedBond onto ONE atom array -- either a model's
   * live atoms or a single vibration frame's copy of them. Every existing
   * bondStyles entry is cleared first, so this is also the un-highlight path.
   */
  function paintBond(atoms) {
    atoms.forEach(function (a) {
      // A falsy entry is "no override" as far as drawBondSticks is concerned,
      // so blanking is enough -- the array itself can stay.
      if (a.bondStyles) {
        for (var i = 0; i < a.bondStyles.length; i++) a.bondStyles[i] = undefined;
      }
    });
    if (selectedBond === null || !P.bonds[selectedBond]) return;
    var sb = P.bonds[selectedBond];
    [[sb[0], sb[1]], [sb[1], sb[0]]].forEach(function (pair) {
      var a = atoms[pair[0]];
      if (!a) return;
      var t = a.bonds.indexOf(pair[1]);
      if (t < 0) return;
      a.bondStyles = a.bondStyles || [];
      a.bondStyles[t] = { color1: BOND_HL, color2: BOND_HL, radius: 0.24 };
    });
  }

  /* Move the highlight on an ALREADY BUILT model, without rebuilding it.
   *
   * Going through renderMode() would work, but it stops the loop, drops the
   * model and re-runs vibrate(), so the molecule snapped back to frame 0 and
   * the animation restarted every time a bond was clicked. The frames are the
   * only place the colour lives (vibrate() copies bondStyles into each one),
   * so repainting all of them in place moves the highlight while the loop
   * keeps running over them.
   */
  function recolorBond(v, m) {
    if (!v || !m || !m.selectedAtoms) return;
    paintBond(m.selectedAtoms({}));
    (m.frames || []).forEach(paintBond);
    // A running loop calls setFrame() on the next tick, which nulls molObj
    // itself and so picks the colour up; a paused viewer never does, and would
    // sit on stale geometry until something else invalidated it.
    m.molObj = null;
    v.render();
  }

  /* Build one mode into a viewer. Shared by the main viewer and the compare
   * viewers, so they cannot drift apart in how a mode is rendered. */
  function renderMode(v, r, o) {
    // Stop the previous loop FIRST. animate() registers timers in
    // viewer.animationTimers; neither removeAllModels() nor addModel() clears
    // them, and several loops calling setFrame() into one canvas is what the
    // atoms flickering looked like.
    v.stopAnimate();
    v.removeAllModels();
    v.removeAllShapes();

    // assignBonds:false -- 3Dmol would otherwise guess bonds by distance, and
    // the picture could then show a different connectivity than the one
    // s[V_S] was computed from. We supply the real bond list below.
    var m = v.addModel(extendedXyz(r), "xyz", { assignBonds: false });
    var atoms = m.selectedAtoms({});
    // drawBondSticks draws each bond from the lower atom.index to the higher,
    // and 3Dmol's xyz parser leaves index null -- null < null is false, so
    // every bond was skipped and no stick was ever drawn.
    atoms.forEach(function (a, i) { a.index = i; });
    P.bonds.forEach(function (b) {
      atoms[b[0]].bonds.push(b[1]); atoms[b[0]].bondOrder.push(1);
      atoms[b[1]].bonds.push(b[0]); atoms[b[1]].bondOrder.push(1);
    });

    // Colour the selected bond IN THE MODEL via atom.bondStyles, which
    // drawBondSticks reads per bond (indexed into atom.bonds). An overlay
    // cylinder cannot follow the vibration -- it sits at the equilibrium
    // geometry while the stick underneath swings through it -- whereas this is
    // part of the model, so vibrate() copies it into every frame. Set before
    // vibrate() for exactly that reason.
    paintBond(atoms);

    m.vibrate(o.frames, o.amp, true, o.arrow);
    // Style AFTER vibrate(): vibrate rebuilds the frames, and a style set
    // beforehand is applied to geometry it then replaces.
    m.setStyle({}, o.bonds
      ? { stick: { radius: 0.15 }, sphere: { scale: 0.25 } }
      : { sphere: { scale: 0.32 } });
    return m;
  }

  function viewOpts() {
    var arrows = document.getElementById("arrows").checked;
    return {
      amp: parseFloat(document.getElementById("amp").value),
      frames: parseInt(document.getElementById("frm").value, 10),
      bonds: document.getElementById("bonds").checked,
      // Red reads clearly against both viewer backgrounds and the CPK palette.
      arrow: arrows ? { color: "#e03131", radius: 0.06, radiusRatio: 2.0, mid: 0.75 }
                    : undefined,
      play: document.getElementById("play").checked
    };
  }

  function draw(r) {
    if (!viewer) return;
    var o = viewOpts();
    model = renderMode(viewer, r, o);
    drawOverlays();
    drawBondPicks();

    // Frame the molecule once. Re-zooming on every mode change fights the
    // user's own rotate/zoom and looks like the view jumping.
    if (!framed) { viewer.zoomTo(); framed = true; }

    if (o.play) {
      viewer.animate({ loop: "backAndForth", interval: 90 });
    } else {
      var done = viewer.setFrame(o.frames);
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
  document.querySelectorAll("#tbl thead th[data-sort]").forEach(function (th) {
    th.addEventListener("click", function () { applySort(th.dataset.sort); });
  });
  document.querySelectorAll(".bhead .bsort").forEach(function (el) {
    el.addEventListener("click", function () {
      // the side panel orders by the mode it is showing
      applyBondSort(el.dataset.bsort,
                    el.dataset.bsort === "val" ? selected : null);
    });
  });

  ["f-S", "f-B", "f-SB", "f-ref"].forEach(function (id) {
    var el = document.getElementById(id);
    if (el) el.addEventListener("change", render);   // f-ref absent under 3N
  });
  ["f-lo", "f-hi"].forEach(function (id) {
    document.getElementById(id).addEventListener("input", render);
  });
  var clearBtn = document.getElementById("cmp-clear");
  if (clearBtn) clearBtn.addEventListener("click", function () {
    picked = [];
    renderCompare();
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
  ["amp", "frm", "arrows", "play", "bonds", "axes", "labels"].forEach(function (id) {
    document.getElementById(id).addEventListener("change", function () {
      if (selected !== null) {
        var r = allRows().find(function (x) { return x.index === selected; });
        if (r) draw(r);
      }
      if (picked.length >= 2) renderCompare();   // keep the compare views in step
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

  // exposed for the headless test harness
  window.__cmp = {
    toggle: function (i) {
      var cb = { checked: picked.indexOf(i) === -1 };
      togglePick(i, cb);
      return cb.checked;
    },
    picked: function () { return picked.slice(); },
    sort: function (k) { applySort(k); },
    sortState: function () { return { key: sortKey, dir: sortDir }; },
    bondSort: function (k, m) { applyBondSort(k, m); },
    bondSortState: function () {
      return { key: bondSort.key, dir: bondSort.dir, mode: bondSort.mode };
    },
    bondOrder: function () {
      return bondOrder(selected === null ? null : byIndex(selected));
    },
    selectBond: function (b) { selectBond(b); },
    selectedBond: function () { return selectedBond; },
    clear: function () { picked = []; renderCompare(); }
  };

  // ------------------------------------------------------------- startup
  render();
  markHeaders();
  markBondHeaders();
  renderCompare();          // set the panel's state from the code, not markup
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
