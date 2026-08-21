/* builder.js -- type a molecule in field by field.
 *
 * Local-instance only (see config.IS_LOCAL). Builds a geometry grid, a bond
 * list and one displacement grid per mode, then assembles exactly the same
 * .vsc text the file upload would carry and posts it in the same field. The
 * server keeps a single parser; nothing here is a second definition of the
 * format.
 */
(function () {
  "use strict";

  var root = document.getElementById("builder");
  if (!root) return;

  var nAtoms = document.getElementById("b-atoms");
  var nModes = document.getElementById("b-modes");
  var geomBox = document.getElementById("b-geom");
  var bondBox = document.getElementById("b-bonds");
  var modeBox = document.getElementById("b-modes-box");
  var hidden = document.getElementById("pasted");

  var ELEMS = ["H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I",
               "B", "Si", "Se", "As", "Xe", "Al", "Ga", "Ge", "Sn", "Te"];

  function el(tag, cls, txt) {
    var e = document.createElement(tag);
    if (cls) e.className = cls;
    if (txt !== undefined) e.textContent = txt;
    return e;
  }

  function num(value, step) {
    var i = el("input");
    i.type = "number";
    i.step = step || "any";
    i.value = value === undefined ? "" : value;
    return i;
  }

  // ------------------------------------------------------------ geometry
  function buildGeometry(n) {
    geomBox.innerHTML = "";
    var head = el("div", "brow head");
    ["#", "element", "x / Å", "y / Å", "z / Å"].forEach(function (t) {
      head.appendChild(el("span", null, t));
    });
    geomBox.appendChild(head);

    for (var i = 1; i <= n; i++) {
      var r = el("div", "brow");
      r.appendChild(el("span", "idx", String(i)));
      var sel = el("select");
      sel.className = "b-el";
      ELEMS.forEach(function (s) {
        var o = el("option", null, s);
        o.value = s;
        sel.appendChild(o);
      });
      sel.value = "C";
      r.appendChild(sel);
      ["x", "y", "z"].forEach(function () { r.appendChild(num("")); });
      geomBox.appendChild(r);
    }
  }

  // -------------------------------------------------------- connectivity
  function bondRow(i, j, n) {
    var r = el("div", "brow bond");
    [i, j].forEach(function (v) {
      var s = el("select");
      for (var a = 1; a <= n; a++) {
        var o = el("option", null, String(a));
        o.value = String(a);
        s.appendChild(o);
      }
      s.value = String(v);
      r.appendChild(s);
    });
    var rm = el("button", "mini", "×");
    rm.type = "button";
    rm.title = "remove this bond";
    rm.addEventListener("click", function () { r.remove(); });
    r.appendChild(rm);
    return r;
  }

  function buildBonds(n) {
    bondBox.innerHTML = "";
    var head = el("div", "brow head bond");
    ["atom", "bonded to", ""].forEach(function (t) {
      head.appendChild(el("span", null, t));
    });
    bondBox.appendChild(head);
    bondBox.appendChild(bondRow(1, Math.min(2, n), n));

    var add = el("button", "mini add", "+ bond");
    add.type = "button";
    add.addEventListener("click", function () {
      bondBox.insertBefore(bondRow(1, Math.min(2, n), n), add);
    });
    bondBox.appendChild(add);
  }

  // --------------------------------------------------------------- modes
  function buildModes(nAt, nMo) {
    modeBox.innerHTML = "";
    for (var m = 1; m <= nMo; m++) {
      var box = el("div", "modeblock");
      var h = el("div", "modehead");
      h.appendChild(el("b", null, "mode " + m));
      var f = num("");
      f.className = "freq";
      f.placeholder = "freq / cm⁻¹ (optional)";
      f.step = "any";
      h.appendChild(f);
      box.appendChild(h);

      var head = el("div", "brow head disp");
      ["#", "dx", "dy", "dz"].forEach(function (t) {
        head.appendChild(el("span", null, t));
      });
      box.appendChild(head);

      for (var a = 1; a <= nAt; a++) {
        var r = el("div", "brow disp");
        r.appendChild(el("span", "idx", String(a)));
        ["dx", "dy", "dz"].forEach(function () { r.appendChild(num("")); });
        box.appendChild(r);
      }
      modeBox.appendChild(box);
    }
  }

  // ------------------------------------------------------------- preview
  var viewer = null, previewTimer = null;

  function readGeometry() {
    var atoms = [], bad = 0;
    Array.prototype.forEach.call(geomBox.querySelectorAll(".brow:not(.head)"),
      function (r) {
        var n = r.querySelectorAll("input");
        var xyz = [0, 1, 2].map(function (i) { return parseFloat(n[i].value); });
        if (xyz.some(isNaN)) { bad++; return; }
        atoms.push({ elem: r.querySelector("select").value, xyz: xyz });
      });
    return { atoms: atoms, incomplete: bad };
  }

  function readBonds(nAtoms) {
    var seen = {}, out = [];
    Array.prototype.forEach.call(
      bondBox.querySelectorAll(".brow.bond:not(.head)"), function (r) {
        var s = r.querySelectorAll("select");
        var i = parseInt(s[0].value, 10) - 1, j = parseInt(s[1].value, 10) - 1;
        if (i === j || i < 0 || j < 0 || i >= nAtoms || j >= nAtoms) return;
        var k = Math.min(i, j) + "-" + Math.max(i, j);
        if (seen[k]) return;
        seen[k] = 1;
        out.push([i, j]);
      });
    return out;
  }

  function status(msg) {
    var el = document.getElementById("b-vstat");
    if (el) el.textContent = msg;
  }

  function preview() {
    var box = document.getElementById("b-view");
    if (!box || typeof $3Dmol === "undefined") return;
    var g = readGeometry();
    if (g.atoms.length < 1) { status("enter coordinates"); return; }

    if (!viewer) {
      try {
        viewer = $3Dmol.createViewer(box, {
          backgroundColor: getComputedStyle(document.body)
            .getPropertyValue("--viewer-bg").trim() || "white"
        });
      } catch (e) { status("viewer unavailable"); return; }
    }

    var bonds = readBonds(g.atoms.length);
    var lines = [String(g.atoms.length), "preview"];
    g.atoms.forEach(function (a) {
      lines.push(a.elem + " " + a.xyz[0] + " " + a.xyz[1] + " " + a.xyz[2]);
    });

    viewer.removeAllModels();
    // assignBonds:false, then wire OUR bonds -- the preview must show the
    // connectivity that will be scored, not a distance-based guess.
    var model = viewer.addModel(lines.join("\n") + "\n", "xyz",
                                { assignBonds: false });
    var at = model.selectedAtoms({});
    // drawBondSticks needs a comparable index on each atom or it draws nothing.
    at.forEach(function (a, i) { a.index = i; });
    bonds.forEach(function (b) {
      at[b[0]].bonds.push(b[1]); at[b[0]].bondOrder.push(1);
      at[b[1]].bonds.push(b[0]); at[b[1]].bondOrder.push(1);
    });
    model.setStyle({}, { stick: { radius: 0.15 }, sphere: { scale: 0.25 } });
    viewer.zoomTo();
    viewer.render();

    status(g.atoms.length + " atom" + (g.atoms.length === 1 ? "" : "s") + ", " +
           bonds.length + " bond" + (bonds.length === 1 ? "" : "s") +
           (g.incomplete ? "  \u00b7 " + g.incomplete + " incomplete" : ""));
  }

  function schedulePreview() {
    clearTimeout(previewTimer);
    previewTimer = setTimeout(preview, 250);   // debounce while typing
  }

  function build() {
    var a = Math.max(2, Math.min(200, parseInt(nAtoms.value, 10) || 3));
    var m = Math.max(1, Math.min(400, parseInt(nModes.value, 10) || 1));
    nAtoms.value = a;
    nModes.value = m;
    buildGeometry(a);
    buildBonds(a);
    buildModes(a, m);
    root.classList.add("ready");
    [geomBox, bondBox].forEach(function (box) {
      box.addEventListener("input", schedulePreview);
      box.addEventListener("change", schedulePreview);
      box.addEventListener("click", schedulePreview);   // + / x buttons
    });
    schedulePreview();
  }

  // ------------------------------------------------------------ assemble
  function val(input) {
    // String(): a real <input>.value is always a string, but coercing keeps
    // this honest if the field is ever populated programmatically.
    var v = String(input.value === undefined || input.value === null ? "" : input.value).trim();
    return v === "" ? "0.000000" : v;      // a blank number field reads as zero
  }

  function toVsc() {
    var lines = ["#VSCORE 1.0", "#TITLE   typed in", "", "[GEOMETRY] Angstrom"];

    var rows = geomBox.querySelectorAll(".brow:not(.head)");
    Array.prototype.forEach.call(rows, function (r, i) {
      var sel = r.querySelector("select");
      var n = r.querySelectorAll("input");
      lines.push([i + 1, sel.value, val(n[0]), val(n[1]), val(n[2])].join("  "));
    });

    lines.push("", "[CONNECTIVITY]");
    var seen = {}, partners = {};
    Array.prototype.forEach.call(bondBox.querySelectorAll(".brow.bond:not(.head)"),
      function (r) {
        var s = r.querySelectorAll("select");
        var i = parseInt(s[0].value, 10), j = parseInt(s[1].value, 10);
        if (i === j) return;                       // a self-bond is not a bond
        var key = Math.min(i, j) + "-" + Math.max(i, j);
        if (seen[key]) return;                     // de-duplicate
        seen[key] = 1;
        (partners[i] = partners[i] || []).push(j);
      });
    Object.keys(partners).forEach(function (i) {
      lines.push(i + "  " + partners[i].join("  "));
    });

    var blocks = modeBox.querySelectorAll(".modeblock");
    lines.push("", "[MODES] " + blocks.length);
    Array.prototype.forEach.call(blocks, function (b, mi) {
      var f = b.querySelector("input.freq");
      var fv = String(f.value == null ? "" : f.value).trim();
      lines.push("  mode " + (mi + 1) + (fv ? "   freq=" + fv : ""));
      Array.prototype.forEach.call(b.querySelectorAll(".brow.disp:not(.head)"),
        function (r, ai) {
          var n = r.querySelectorAll("input");
          lines.push([ai + 1, val(n[0]), val(n[1]), val(n[2])].join("  "));
        });
    });
    return lines.join("\n") + "\n";
  }

  // --------------------------------------------------------------- wire
  document.getElementById("b-build").addEventListener("click", build);
  var form = document.getElementById("form");
  if (form) {
    form.addEventListener("submit", function () {
      if (document.getElementById("mode").value === "paste") {
        hidden.value = toVsc();
      }
    });
  }
  // expose for the headless test harness
  window.__builder = { build: build, toVsc: toVsc, preview: preview,
                       readGeometry: readGeometry, readBonds: readBonds };
})();
