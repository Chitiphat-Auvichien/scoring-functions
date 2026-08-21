"""Runs app.js against the REAL 3Dmol library under JavaScriptCore.

Why this exists: an earlier harness stubbed out $3Dmol entirely, so it happily
"passed" while the viewer threw `undefined is not an object (evaluating
'this.atoms.length')` in a real browser -- addModel() with no data leaves
this.atoms undefined inside addMolData. These tests load the actual
3Dmol-min.js and only stub WebGL rendering, so the parsing, bonding and
vibrate() paths are genuinely exercised.

Skipped automatically if jsc is unavailable (i.e. off macOS).
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
JS = ROOT / "tests" / "js"
JSC = Path("/System/Library/Frameworks/JavaScriptCore.framework/"
           "Versions/A/Helpers/jsc")
sys.path.insert(0, str(ROOT))

from app.core.parsers import parse_connectivity, parse_gaussian_log  # noqa: E402
from app.core.pipeline import analyse  # noqa: E402
from app.main import _embed  # noqa: E402

def _find_reference_data():
    """Locate the reference data whichever way the app is laid out.

    The webapp lives INSIDE scoring-functions/ (so the data is ../data), but it
    was developed as a sibling directory (../scoring-functions/data). Hardcoding
    one of those made pytestmark's skipif silently skip the whole suite when the
    app moved -- 42 skipped, reported as success. Try both, and say which was
    found so a miss is visible.
    """
    for cand in (ROOT.parent / "data",
                 ROOT.parent / "scoring-functions" / "data"):
        if (cand / "logs").is_dir() and (cand / "results").is_dir():
            return cand
    return ROOT.parent / "data"          # nonexistent -> skipif reports it


REF = _find_reference_data()

pytestmark = [
    pytest.mark.skipif(not JSC.exists(), reason="jsc (JavaScriptCore) not available"),
    pytest.mark.skipif(not REF.exists(), reason="reference data not present"),
]

CASES = ["HOCl", "CO2", "C6H6", "C10H8", "XeF2Cl2"]


def _run_js(stem, tmp_path, probe, payload=None, stub="dom_stub.js"):
    """Drive app.js over a payload. Defaults to the .log path for `stem`;
    pass `payload` to exercise a case the .log path cannot produce."""
    if payload is None:
        log = REF / "logs" / f"{stem}.log"
        com = REF / "gjf" / f"{stem}.com"
        g = parse_gaussian_log(log.read_text(errors="replace"))
        bonds = parse_connectivity(com.read_text(errors="replace"), len(g["atoms"]))
        payload = analyse(g["atoms"], g["coords"], bonds, g["modes"], title=stem)

    data = tmp_path / "data.js"
    data.write_text("var PAYLOAD_TXT=%s; var DL_TXT=%s;" % (
        json.dumps(_embed(payload)),
        json.dumps(_embed({"csv": "", "csvall": "", "vsc": ""}))))

    out = subprocess.run(
        [str(JSC), str(JS / "dom.js"), str(ROOT / "app/static/3Dmol-min.js"),
         str(data), str(JS / stub), str(JS / "real_viewer.js"),
         str(ROOT / "app/static/app.js"), "-e", probe],
        capture_output=True, text=True, cwd=ROOT)
    assert out.returncode == 0, out.stderr or out.stdout
    return payload, out.stdout.strip()


@pytest.mark.parametrize("stem", CASES)
def test_viewer_starts_without_error(stem, tmp_path):
    """The startup try/catch must not have written a failure into #mol."""
    _, got = _run_js(stem, tmp_path,
                     "print(document.getElementById('mol').innerHTML || 'NONE');")
    assert got == "NONE", f"{stem}: viewer failed to start -> {got}"


@pytest.mark.parametrize("stem", CASES)
def test_model_atoms_bonds_and_frames(stem, tmp_path):
    """Real 3Dmol must parse every atom, take OUR bonds, and build frames."""
    payload, got = _run_js(stem, tmp_path, """
        var m = CALLS.viewer.models[0], a = m.selectedAtoms({});
        var nb = 0; a.forEach(function (x) { nb += x.bonds.length; });
        print(JSON.stringify({atoms:a.length, bonds:nb/2,
                              frames:m.frames ? m.frames.length : 0,
                              dx:a[0].dx, dy:a[0].dy, dz:a[0].dz}));
    """)
    r = json.loads(got)
    assert r["atoms"] == payload["n_atoms"]
    # Bond count must match the connectivity we supplied -- if 3Dmol were
    # guessing by distance this would drift.
    assert r["bonds"] == len(payload["bonds"])
    assert r["frames"] == 20                      # 10 frames, bothWays
    first = payload["vibrations"][0]["vector"][0]
    assert r["dx"] == pytest.approx(first[0], abs=1e-5)
    assert r["dy"] == pytest.approx(first[1], abs=1e-5)
    assert r["dz"] == pytest.approx(first[2], abs=1e-5)


def test_vibrate_called_with_expected_arguments(tmp_path):
    _, got = _run_js("HOCl", tmp_path, "print(JSON.stringify(CALLS.vibrate));")
    v = json.loads(got)
    assert v["frames"] == 10 and v["amp"] == 1 and v["both"] is True
    assert v["arrow"]["color"] == "#e03131"      # red displacement arrows


# ----------------------------------------------------------------------
# Regression guard for the flicker bug.
# ----------------------------------------------------------------------
def test_animation_loops_do_not_stack(tmp_path):
    """Exactly one animation loop may be live, however much the user clicks.

    viewer.animate() registers timers in viewer.animationTimers and bumps an
    internal counter; neither removeAllModels() nor addModel() clears them. Any
    draw() path that forgets stopAnimate() therefore leaves the previous loop
    running, and several loops calling setFrame() on different indices into one
    canvas is what "the atoms flicker a lot" looks like.
    """
    _, got = _run_js("C6H6", tmp_path, """
        var out = {start: CALLS.live};
        var trs = document.querySelectorAll('#tbl tbody tr');
        for (var i = 0; i < 8; i++) trs[i % trs.length].click();
        out.afterClicks = CALLS.live;
        var a = document.getElementById('amp'); a.value = '2';
        for (var j = 0; j < 4; j++) a.fire('change');
        out.afterSliders = CALLS.live;
        out.started = CALLS.animateCount; out.stopped = CALLS.stopCount;
        out.zooms = CALLS.zooms;
        print(JSON.stringify(out));
    """)
    r = json.loads(got)
    assert r["start"] == 1
    assert r["afterClicks"] == 1, f"{r['afterClicks']} loops live after 8 clicks"
    assert r["afterSliders"] == 1, f"{r['afterSliders']} loops live after sliders"
    # every start must be matched by a stop
    assert r["started"] == r["stopped"]
    # and the camera is framed once, not re-zoomed on every mode change
    assert r["zooms"] == 1, f"camera re-zoomed {r['zooms']} times"


def test_pause_holds_equilibrium_frame(tmp_path):
    """Unchecking 'animate' should park on the undisplaced geometry."""
    _, got = _run_js("HOCl", tmp_path, """
        var pl = document.getElementById('play');
        pl.checked = false; pl.fire('change');
        print(JSON.stringify({live: CALLS.live, frame: CALLS.frame,
                              total: CALLS.vibrate.frames * 2}));
    """)
    r = json.loads(got)
    assert r["live"] == 0
    assert r["frame"] == r["total"] // 2      # midpoint == equilibrium


def test_frequency_filter_does_not_hide_reference_modes(tmp_path):
    """The constructed T/R modes must survive the frequency range filter.

    They carry a synthetic frequency of 0, while the range defaults to the real
    modes' own min/max -- so applying the filter to them hid all six on every
    real molecule (ethylene showed "12 of 18"). Their visibility belongs to the
    'T/R reference modes' checkbox alone.
    """
    payload, got = _run_js("C2H4", tmp_path, """
        var P = JSON.parse(PAYLOAD_TXT);
        // reproduce the real page default: low bound = lowest real frequency
        var lo = document.getElementById('f-lo');
        lo.value = String(Math.floor(P.freq_range[0]));
        lo.fire('input');
        var refs = 0, vibs = 0;
        document.querySelectorAll('#tbl tbody tr').forEach(function (tr) {
            (tr.classList.contains('ref') ? refs++ : vibs++);
        });
        print(JSON.stringify({refs: refs, vibs: vibs, lo: lo.value}));
    """)
    r = json.loads(got)
    assert float(r["lo"]) > 0, "test setup: low bound should be above 0"
    assert r["refs"] == 6, f"{r['refs']} reference rows visible with lo={r['lo']}"
    assert r["vibs"] == payload["n_modes"]


def test_reference_rows_show_no_frequency(tmp_path):
    """A constructed mode is not a Hessian solution -- show a dash, not 0.00."""
    _, got = _run_js("C2H4", tmp_path, """
        var tr = document.querySelectorAll('#tbl tbody tr')[0];
        var m = tr.innerHTML.match(/<td class="n[^"]*">\\s*([^<]*?)\\s*</);
        print(JSON.stringify({isRef: tr.classList.contains('ref'), freq: m[1]}));
    """)
    r = json.loads(got)
    assert r["isRef"] is True
    assert r["freq"] == "—", f"reference frequency rendered as {r['freq']!r}"


def test_dom_stub_defaults_match_the_template():
    """The JS harness must start from the same control defaults as the page.

    tests/js/dom_stub.js hardcodes initial checkbox/input values. When the
    template's defaults changed (f-ref became checked) the stub kept the old
    value, and the reference-mode tests passed against a state the real page
    never has. Compare the two directly rather than trusting they stay in sync.
    """
    tpl = (ROOT / "app/templates/result.html").read_text()
    stub = (JS / "dom_stub.js").read_text()
    for control in ("f-ref", "f-S", "f-B", "f-SB", "arrows", "play"):
        checked_in_tpl = bool(re.search(
            r'id="%s"[^>]*\bchecked\b' % re.escape(control), tpl))
        m = re.search(r'get\("%s"\)\.checked\s*=\s*(true|false)'
                      % re.escape(control), stub)
        if m is None:
            continue                     # stub leaves it at its own default
        checked_in_stub = m.group(1) == "true"
        assert checked_in_stub == checked_in_tpl, (
            f"{control}: template checked={checked_in_tpl} but "
            f"dom_stub.js sets {checked_in_stub}")


def test_frequency_bounds_never_exclude_the_extreme_modes(tmp_path):
    """The default range must include every mode, at both ends.

    The bounds were rendered with Jinja's round(), which can move a bound PAST
    the mode it should include: HCOOH's lowest frequency is 610.9187, rounds to
    611, and the `frequency < lo` test then hid it ("20 of 21 shown"). The same
    happens at the top when the highest frequency rounds down. floor/ceil fix
    both ends.
    """
    import math
    for stem in ("HCOOH", "C2H4", "C6H6", "C10H8", "HOCl"):
        payload, got = _run_js(stem, tmp_path, """
            var P = JSON.parse(PAYLOAD_TXT);
            var lo = document.getElementById('f-lo'), hi = document.getElementById('f-hi');
            lo.value = String(Math.floor(P.freq_range[0]));
            hi.value = String(Math.ceil(P.freq_range[1]));
            lo.fire('input');
            var n = 0;
            document.querySelectorAll('#tbl tbody tr').forEach(function () { n++; });
            print(JSON.stringify({rows: n}));
        """)
        want = payload["n_modes"] + len(payload["references"])
        assert json.loads(got)["rows"] == want, \
            f"{stem}: {json.loads(got)['rows']} of {want} rows shown"


def test_template_uses_floor_and_ceil_for_the_bounds():
    """Guard the fix at its source, not just its effect."""
    tpl = (ROOT / "app/templates/result.html").read_text()
    assert "round(0, 'floor')" in tpl, "f-lo must floor, not round"
    assert "round(0, 'ceil')" in tpl, "f-hi must ceil, not round"
    assert "freq_range[0]|round(0)|int" not in tpl
    assert "freq_range[1]|round(0)|int" not in tpl


# ----------------------------------------------------------------------
# Bonds must actually render, and the toggle must switch them.
# ----------------------------------------------------------------------
_STICK_PROBE = """
    var m = CALLS.viewer.models[0];
    function verts(style) {
      var atoms = m.selectedAtoms({}).map(function (a) {
        var c = {}; for (var k in a) c[k] = a[k]; c.style = style; return c;
      });
      var obj = m.createMolObj(atoms, {}), t = 0;
      (obj.children || []).forEach(function (c) {
        if (c.geometry && c.geometry.geometryGroups)
          c.geometry.geometryGroups.forEach(function (g) { t += g.vertices || 0; });
      });
      return t;
    }
    print(JSON.stringify({
      indices: m.selectedAtoms({}).map(function (x) { return x.index; }),
      stick: verts({stick: {radius: 0.15}}),
      sphere: verts({sphere: {scale: 0.32}})
    }));
"""


@pytest.mark.parametrize("stem", ["HOCl", "C6H6", "C10H8", "CO2"])
def test_bond_sticks_are_actually_generated(stem, tmp_path):
    """Sticks must produce geometry, not merely be requested in the style.

    3Dmol's drawBondSticks draws each bond once, from the lower-index atom to
    the higher: `if (atom.index < partner.index)`. Its xyz parser sets `serial`
    but leaves `index` null, and `null < null` is false -- so every bond was
    silently skipped and the viewer drew spheres floating unconnected. The
    style said `stick`, the atoms carried the right bonds, and nothing was
    rendered. Only the vertex count catches that.
    """
    _, got = _run_js(stem, tmp_path, _STICK_PROBE)
    r = json.loads(got)
    assert r["indices"] == list(range(len(r["indices"]))), \
        f"atom.index must be 0..N-1, got {r['indices'][:6]}"
    assert r["stick"] > 0, "stick style produced no geometry -- bonds are not drawn"
    assert r["sphere"] > 0


def test_bonds_toggle_switches_the_style(tmp_path):
    _, got = _run_js("C6H6", tmp_path, """
        function fire(id) { document.getElementById(id).fire('change'); }
        var out = {on: CALLS.style};
        var b = document.getElementById('bonds');
        b.checked = false; fire('bonds'); out.off = CALLS.style;
        b.checked = true;  fire('bonds'); out.back = CALLS.style;
        out.loops = CALLS.live;
        print(JSON.stringify(out));
    """)
    r = json.loads(got)
    assert "stick" in r["on"] and "sphere" in r["on"]
    assert "stick" not in r["off"], "bonds off must drop the stick style"
    assert "sphere" in r["off"], "spheres must remain when bonds are hidden"
    assert "stick" in r["back"]
    assert r["loops"] == 1, "toggling bonds must not stack animation loops"


# ----------------------------------------------------------------------
# Highlighting follows the label; axes are drawn in the principal frame.
# ----------------------------------------------------------------------
def test_mixed_external_highlights_both_scores(tmp_path):
    """A starred label is part external and part vibration, so one number
    cannot describe it: highlight V_S AND the largest |T/R|."""
    log = REF / "logs" / "C6H6.log"
    payload, got = _run_js("C6H6", tmp_path, """
        var K = ['Tx','Ty','Tz','Rx','Ry','Rz','V_S'], out = [];
        document.querySelectorAll('#tbl tbody tr').forEach(function (tr) {
            var h = tr.innerHTML, hot = [], ci = 0, re = /<td class="(n[^"]*)">/g, m;
            while ((m = re.exec(h)) !== null) {
                if (ci > 0 && ci < 8 && m[1].indexOf('hot') > -1) hot.push(K[ci - 1]);
                ci++;
            }
            var lab = (h.match(/class="lab ([^"]*)">([^<]*)</) || [])[2];
            out.push({label: lab, hot: hot});
        });
        print(JSON.stringify(out));
    """)
    rows = json.loads(got)
    assert rows, "no rows rendered"
    for r in rows:
        lab = r["label"]
        if lab and lab[0] in "TR" and lab.endswith("*"):
            assert len(r["hot"]) == 2, f"{lab}: {r['hot']}"
            assert "V_S" in r["hot"]
            assert any(k != "V_S" for k in r["hot"])
        elif lab and lab[0] in "TR":
            assert r["hot"] == [lab], f"clean external {lab} highlighted {r['hot']}"
        else:
            assert r["hot"] == ["V_S"], f"internal {lab} highlighted {r['hot']}"


@pytest.mark.parametrize("stem", ["HOCl", "C6H6", "CO2"])
def test_axes_are_drawn_and_toggle(stem, tmp_path):
    """Three labelled arrows along the principal axes, and a toggle."""
    payload, got = _run_js(stem, tmp_path, """
        function fire(id) { document.getElementById(id).fire('change'); }
        var on = {n: CALLS.arrows.length,
                  labels: CALLS.labels.map(function (l) { return l.text; }),
                  colors: CALLS.arrows.map(function (a) { return a.color; }),
                  len: CALLS.arrows.length ? CALLS.arrows[0].end.x : 0};
        var a = document.getElementById('axes');
        a.checked = false; fire('axes');
        var off = {n: CALLS.arrows.length, labels: CALLS.labels.length};
        print(JSON.stringify({on: on, off: off, loops: CALLS.live}));
    """)
    r = json.loads(got)
    assert r["on"]["n"] == 3, f"expected 3 axis arrows, got {r['on']['n']}"
    assert r["on"]["labels"] == ["x", "y", "z"]
    assert len(set(r["on"]["colors"])) == 3, "axes must be distinguishable"
    assert "#e03131" not in r["on"]["colors"], \
        "axis colour clashes with the red displacement arrows"
    # axis length must scale past the molecule, not sit inside it
    span = max(abs(c) for g in payload["geometry"] for c in g)
    assert r["on"]["len"] > span, f"axis {r['on']['len']} shorter than span {span}"
    assert r["off"]["n"] == 0 and r["off"]["labels"] == 0, "toggle must clear axes"
    assert r["loops"] == 1, "toggling axes must not stack animation loops"


def _payload_with_mixed_externals():
    """A 3N pool whose T/R modes are built about NON-principal axes, so the
    external slots come back starred. The .log path cannot produce this: it
    constructs clean principal-frame references."""
    from app.core.scoring import ModeScorer
    log = REF / "logs" / "HCOOH.log"
    com = REF / "gjf" / "HCOOH.com"
    g = parse_gaussian_log(log.read_text(errors="replace"))
    bonds = parse_connectivity(com.read_text(errors="replace"), len(g["atoms"]))
    sc = ModeScorer(g["atoms"], g["coords"], bonds)      # COM only, not aligned
    modes = [{"frequency": None, "vector": m["vector"]}
             for m in sc.construct_T() + sc.construct_R()]
    modes += [{"frequency": m["frequency"], "vector": m["vector"]} for m in g["modes"]]
    return analyse(g["atoms"], g["coords"], bonds, modes, mode_set="3n")


def test_starred_external_shows_a_second_label_chip(tmp_path):
    """A mixed external renders as two chips, e.g. [Tx*][SB] -- not as the raw
    annotation text "vibration=SB"."""
    pay = _payload_with_mixed_externals()
    assert any(r["label"].endswith("*") for r in pay["vibrations"]), \
        "test setup produced no mixed external"
    _, got = _run_js(None, tmp_path, """
        var out = [];
        document.querySelectorAll('#tbl tbody tr').forEach(function (tr) {
            var cell = (tr.innerHTML.match(/<td>(<b class="lab[\\s\\S]*?)<\\/td>/) || [])[1] || '';
            var chips = (cell.match(/<b class="lab [^"]*">([^<]*)<\\/b>/g) || [])
                        .map(function (c) { return c.replace(/<[^>]*>/g, ''); });
            out.push({chips: chips, raw: cell});
        });
        print(JSON.stringify(out));
    """, payload=pay)
    rows = json.loads(got)
    assert rows
    starred = [r for r in rows if r["chips"] and r["chips"][0].endswith("*")]
    assert starred, "no mixed external rendered"
    for r in starred:
        assert len(r["chips"]) == 2, f"expected two chips, got {r['chips']}"
        assert r["chips"][1] in ("S", "B", "SB"), r["chips"]
        assert "vibration=" not in r["raw"], "raw annotation text must not be shown"
    for r in rows:
        if r["chips"] and not r["chips"][0].endswith("*"):
            assert len(r["chips"]) == 1, f"unstarred row got {r['chips']}"


def test_3n_hides_the_reference_toggle_and_js_survives_it(tmp_path):
    """Under 3N nothing is constructed, so the T/R toggle is not rendered.

    app.js must then not dereference the missing element: rows() previously did
    document.getElementById("f-ref").checked unconditionally, which throws and
    takes the whole table with it.
    """
    from app.core.parsers import parse_vsc
    vsc = ROOT / "examples" / "C6H6_EMIT.vsc"
    if not vsc.exists():
        pytest.skip("EMIT example not present")
    v = parse_vsc(vsc.read_text())
    pay = analyse(v["atoms"], v["coords"], v["bonds"], v["modes"], mode_set="3n")
    assert pay["references"] == []

    _, got = _run_js(None, tmp_path, """
        // mirror the page defaults; the stub's are hardcoded and EMIT
        // eigenvalues run negative
        var P = JSON.parse(PAYLOAD_TXT);
        var lo = document.getElementById('f-lo'), hi = document.getElementById('f-hi');
        lo.value = String(Math.floor(P.freq_range[0]));
        hi.value = String(Math.ceil(P.freq_range[1]));
        lo.fire('input');
        var n = 0;
        document.querySelectorAll('#tbl tbody tr').forEach(function () { n++; });
        print(JSON.stringify({
            err: document.getElementById('mol').innerHTML || '',
            rows: n,
            toggle: document.getElementById('f-ref') === null,
            loops: CALLS.live
        }));
    """, payload=pay, stub="dom_stub_noref.js")
    r = json.loads(got)
    assert r["err"] == "", f"app.js failed without the toggle: {r['err']}"
    assert r["toggle"] is True, "test setup: f-ref should be absent"
    assert r["rows"] == pay["n_modes"]
    assert r["loops"] == 1


def test_template_only_renders_the_toggle_with_references():
    tpl = (ROOT / "app/templates/result.html").read_text()
    i = tpl.index('id="f-ref"')
    before = tpl[:i]
    assert "{% if p.references %}" in before[-400:], \
        "the T/R toggle must be guarded by {% if p.references %}"


# ----------------------------------------------------------------------
# Mode comparison: up to three at once, with their per-bond s_AB.
# ----------------------------------------------------------------------
def test_compare_panel_appears_only_with_two_or_more(tmp_path):
    _, got = _run_js("HOCl", tmp_path, """
        var C = window.__cmp, box = document.getElementById('compare');
        var out = {start: box.hidden};
        C.toggle(6); out.one = box.hidden;
        C.toggle(7); out.two = box.hidden;
        out.count = document.getElementById('cmp-n').textContent;
        C.clear();   out.cleared = box.hidden;
        print(JSON.stringify(out));
    """)
    r = json.loads(got)
    assert r["start"] is True, "must start hidden, from the code not the markup"
    assert r["one"] is True, "one mode is not a comparison"
    assert r["two"] is False
    assert r["count"] == "2"
    assert r["cleared"] is True


def test_at_most_three_modes_can_be_compared(tmp_path):
    """The fourth tick is refused rather than silently evicting an earlier one."""
    _, got = _run_js("C6H6", tmp_path, """
        var C = window.__cmp;
        [6, 7, 8].forEach(function (i) { C.toggle(i); });
        var before = C.picked();
        var accepted = C.toggle(9);
        print(JSON.stringify({before: before, accepted: accepted, after: C.picked()}));
    """)
    r = json.loads(got)
    assert len(r["before"]) == 3
    assert r["accepted"] is False, "the fourth tick must be refused"
    assert r["after"] == r["before"], "an earlier choice must not be dropped"


def test_compare_shows_scores_and_per_bond_for_each_mode(tmp_path):
    payload, got = _run_js("HOCl", tmp_path, """
        var C = window.__cmp;
        C.toggle(6); C.toggle(8);
        function txt(id) {
          return document.getElementById(id).innerHTML
                 .replace(/<\\/tr>/g, '\\n').replace(/<[^>]*>/g, ' ')
                 .replace(/[ \\t]+/g, ' ');
        }
        print(JSON.stringify({scores: txt('cmp-scores'), bonds: txt('cmp-bonds')}));
    """)
    r = json.loads(got)
    # every score, the label and the frequency row
    for key in ("Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "V_S", "label"):
        assert key in r["scores"], key
    assert "Vib 1" in r["scores"] and "Vib 3" in r["scores"]

    # one row per bond, values matching the payload
    vib = {v["name"]: v for v in payload["vibrations"]}
    for b in vib["Vib 1"]["bonds"]:
        assert b["pair"] in r["bonds"], b["pair"]
    a = vib["Vib 1"]["bonds"][1]["s_AB"]          # O1-Cl3 on the O-Cl stretch
    assert ("%+.3f" % a) in r["bonds"], r["bonds"]


def test_untick_removes_only_that_mode(tmp_path):
    _, got = _run_js("C6H6", tmp_path, """
        var C = window.__cmp;
        [6, 7, 8].forEach(function (i) { C.toggle(i); });
        C.toggle(7);
        print(JSON.stringify(C.picked()));
    """)
    assert json.loads(got) == [6, 8]


def test_compare_shows_one_animated_viewer_per_mode(tmp_path):
    """Side-by-side animation is the point of comparing."""
    _, got = _run_js("C6H6", tmp_path, """
        var C = window.__cmp, host = document.getElementById('cmp-views');
        function vis(){ return host.children.filter(function(c){return !c.hidden;}).length; }
        C.toggle(6); C.toggle(7);
        var two = {visible: vis(), cls: host.className};
        C.toggle(8);
        var three = {visible: vis(), cls: host.className};
        C.clear();
        print(JSON.stringify({two: two, three: three, cleared: vis(),
                              loops: CALLS.live}));
    """)
    r = json.loads(got)
    assert r["two"]["visible"] == 2 and r["two"]["cls"] == "cviews n2"
    assert r["three"]["visible"] == 3 and r["three"]["cls"] == "cviews n3"
    assert r["cleared"] == 0


def test_compare_viewers_do_not_leak_webgl_contexts(tmp_path):
    """Created once and reused, never per tick.

    A browser caps concurrent WebGL contexts and 3Dmol's GLViewer exposes no
    teardown, so making a fresh viewer on every toggle would exhaust them and
    the panels would quietly stop drawing.
    """
    _, got = _run_js("C6H6", tmp_path, """
        var C = window.__cmp;
        C.toggle(6); C.toggle(7);
        var after2 = CALLS.viewers.length;
        for (var i = 0; i < 8; i++) { C.toggle(8); C.toggle(8); }
        var churned = CALLS.viewers.length;
        C.toggle(8);
        print(JSON.stringify({after2: after2, churned: churned,
                              atThree: CALLS.viewers.length}));
    """)
    r = json.loads(got)
    assert r["after2"] == 4, f"expected 1 main + 3 compare, got {r['after2']}"
    assert r["churned"] == r["after2"], "contexts grew while toggling"
    assert r["atThree"] == r["after2"]


def test_every_compared_viewer_is_resized_when_shown(tmp_path):
    """A canvas sized while its card was display:none draws into a stale buffer.

    The third panel came up blank for exactly this reason: the viewers are
    created once with all three cards laid out, then hidden, and the grid width
    changes with how many are shown. Each must be told its box again.
    """
    _, got = _run_js("C6H6", tmp_path, """
        var C = window.__cmp;
        CALLS.resizes = 0;
        C.toggle(6); C.toggle(7); C.toggle(8);
        var host = document.getElementById('cmp-views');
        print(JSON.stringify({
          visible: host.children.filter(function (c) { return !c.hidden; }).length,
          resizes: CALLS.resizes
        }));
    """)
    r = json.loads(got)
    assert r["visible"] == 3
    assert r["resizes"] >= 3, "each shown viewer must be resized"


def test_bond_table_puts_each_mode_in_its_own_column(tmp_path):
    """One row per bond, one cell per mode.

    A <td> made into a flex container drops out of the table's column model,
    so the mode cells stacked vertically instead of sitting side by side.
    """
    payload, got = _run_js("HOCl", tmp_path, """
        var C = window.__cmp;
        C.toggle(6); C.toggle(7); C.toggle(8);
        var html = document.getElementById('cmp-bonds').innerHTML;
        var body = html.split('</tr>').slice(1).filter(function (r) {
            return r.indexOf('<td') !== -1;
        });
        print(JSON.stringify({
          rows: body.length,
          cellsPerRow: body.map(function (r) { return (r.match(/<td/g) || []).length; }),
          tdIsFlex: /td class="cv2"/.test(html)
        }));
    """)
    r = json.loads(got)
    nbonds = len(payload["vibrations"][0]["bonds"])
    assert r["rows"] == nbonds, f"{r['rows']} rows for {nbonds} bonds"
    # 1 bond-name cell + one per compared mode
    assert all(c == 4 for c in r["cellsPerRow"]), r["cellsPerRow"]
    assert r["tdIsFlex"] is False, "the <td> must stay a table cell"


def test_compare_viewers_share_one_camera(tmp_path):
    """Interacting with any panel moves the others by the same amount.

    Deliberately NOT 3Dmol's linkViewer(): that propagates from show(), which
    runs on every render including every animation frame, so three linked
    viewers push a full setView+render into each other every 90 ms -- 9 renders
    per tick instead of 3, each overwriting the others' cameras. That is what
    broke with three modes selected. Syncing on interaction costs nothing while
    idle.
    """
    _, got = _run_js("C6H6", tmp_path, """
        var C = window.__cmp, host = document.getElementById('cmp-views');
        C.toggle(6); C.toggle(7); C.toggle(8);
        var vs = CALLS.viewers.slice(1);
        var out = { links: vs.map(function (v) { return v.linkedViewers.length; }) };

        vs[0]._view = [9,8,7,6,5,4,3,2];
        host.children[0].querySelector('.cv').fire('mousemove', {buttons: 1});
        out.fromFirst = [vs[1].getView(), vs[2].getView()];

        vs[2]._view = [2,2,2,2,2,2,2,2];
        host.children[2].querySelector('.cv').fire('wheel', {});
        out.fromLast = [vs[0].getView(), vs[1].getView()];

        CALLS.setViews = 0;
        host.children[0].querySelector('.cv').fire('mousemove', {buttons: 0});
        out.idleSyncs = CALLS.setViews;
        print(JSON.stringify(out));
    """)
    r = json.loads(got)
    assert r["links"] == [0, 0, 0], \
        "linkViewer must not be used -- it propagates on every animation frame"
    assert r["fromFirst"] == [[9,8,7,6,5,4,3,2]] * 2, "drag did not reach the others"
    assert r["fromLast"] == [[2,2,2,2,2,2,2,2]] * 2, "sync must work from any panel"
    assert r["idleSyncs"] == 0, "a mousemove with no button held is not an interaction"


# ----------------------------------------------------------------------
# Sorting the results table by column.
# ----------------------------------------------------------------------
_NAMES = """
        function names() {
          return document.querySelectorAll('#tbl tbody tr').map(function (tr) {
            return (tr.innerHTML.match(/class="mode">([^<]*)</) || [])[1];
          });
        }
"""


def test_sort_cycles_ascending_descending_then_file_order(tmp_path):
    _, got = _run_js("C10H8", tmp_path, _NAMES + """
        var C = window.__cmp, out = {};
        out.original = names();
        C.sort('freq'); out.asc = names();  out.ascState = C.sortState();
        C.sort('freq'); out.desc = names(); out.descState = C.sortState();
        C.sort('freq'); out.back = names(); out.backState = C.sortState();
        print(JSON.stringify(out));
    """)
    r = json.loads(got)
    assert r["ascState"] == {"key": "freq", "dir": 1}
    assert r["descState"] == {"key": "freq", "dir": -1}
    assert r["backState"] == {"key": None, "dir": 0}
    assert r["back"] == r["original"], "a third click must restore file order"
    assert r["asc"] != r["original"]
    # descending is the ascending order reversed, apart from the no-frequency rows
    vib_asc = [n for n in r["asc"] if n.startswith("Vib")]
    vib_desc = [n for n in r["desc"] if n.startswith("Vib")]
    assert vib_desc == list(reversed(vib_asc))


def test_sort_by_score_matches_the_payload(tmp_path):
    payload, got = _run_js("C10H8", tmp_path, _NAMES + """
        var C = window.__cmp;
        C.sort('V_S');
        print(JSON.stringify(names()));
    """)
    order = [n for n in json.loads(got) if n.startswith("Vib")]
    by_score = sorted(payload["vibrations"], key=lambda r: (r["scores"]["V_S"], r["index"]))
    assert order == [r["name"] for r in by_score]


def test_rows_without_a_frequency_sort_last_in_both_directions(tmp_path):
    """A reference mode shows a dash, not 0.00, so it must not sort as zero."""
    _, got = _run_js("C10H8", tmp_path, _NAMES + """
        var C = window.__cmp, out = {};
        C.sort('freq'); out.asc = names().slice(-6);
        C.sort('freq'); out.desc = names().slice(-6);
        print(JSON.stringify(out));
    """)
    r = json.loads(got)
    refs = {"Tx", "Ty", "Tz", "Rx", "Ry", "Rz"}
    assert set(r["asc"]) == refs, r["asc"]
    assert set(r["desc"]) == refs, r["desc"]


def test_sort_indicator_marks_only_the_active_column(tmp_path):
    _, got = _run_js("C10H8", tmp_path, """
        var C = window.__cmp;
        function heads() {
          return document.querySelectorAll('#tbl thead th').map(function (th) {
            var ind = th.querySelector('.ind');
            return {k: th.dataset.sort, on: th.classList.contains('sorted'),
                    ind: ind ? ind.textContent.trim() : ''};
          });
        }
        C.sort('Tx');
        print(JSON.stringify(heads()));
    """)
    hs = json.loads(got)
    active = [h for h in hs if h["on"]]
    assert len(active) == 1 and active[0]["k"] == "Tx"
    assert active[0]["ind"] == "\u25b2"
    assert all(h["ind"] == "" for h in hs if not h["on"])


# ----------------------------------------------------------------------
# Bond selection, both directions.
# ----------------------------------------------------------------------
def test_every_bond_gets_an_invisible_click_target(tmp_path):
    """3Dmol's sticks belong to the model and are not individually pickable, so
    a cylinder is laid over each bond to catch the click. It is never drawn --
    the highlight itself lives in the model."""
    payload, got = _run_js("C10H8", tmp_path, """
        print(JSON.stringify({
          n: CALLS.cylinders.length,
          clickable: CALLS.cylinders.every(function (c) { return c.clickable; }),
          visible: CALLS.cylinders.filter(function (c) { return c.opacity > 0; }).length
        }));
    """)
    r = json.loads(got)
    assert r["n"] == len(payload["bonds"]), "one target per bond"
    assert r["clickable"] is True
    assert r["visible"] == 0, "the click targets must never be drawn"


_BOND_STATE = """
    function state(bi) {
      var P = JSON.parse(PAYLOAD_TXT);
      var m = CALLS.viewer.models[0], a = m.selectedAtoms({});
      var b = P.bonds[bi];
      var t0 = a[b[0]].bonds.indexOf(b[1]);
      var neighbour = null;
      a[b[0]].bonds.forEach(function (x, k) {
        if (x !== b[1] && neighbour === null) neighbour = k;
      });
      return {
        styled: a[b[0]].bondStyles ? a[b[0]].bondStyles[t0] || null : null,
        other: (neighbour !== null && a[b[0]].bondStyles)
                 ? a[b[0]].bondStyles[neighbour] || null : null,
        inFrames: (m.frames && m.frames[3] && m.frames[3][b[0]].bondStyles)
                    ? m.frames[3][b[0]].bondStyles[t0] || null : null
      };
    }
"""


def test_clicking_a_bond_in_3d_highlights_its_score_row(tmp_path):
    _, got = _run_js("C10H8", tmp_path, _BOND_STATE + """
        var C = window.__cmp, bl = document.getElementById('bondlist');
        CALLS.cylinders[3].callback();                  // as 3Dmol would
        print(JSON.stringify({
          selected: C.selectedBond(),
          row: (bl._html.match(/brow sel" data-b="(\\d+)"/) || [])[1],
          st: state(3)
        }));
    """)
    r = json.loads(got)
    assert r["selected"] == 3
    assert r["row"] == "3", "the matching score row must be marked"
    assert r["st"]["styled"], "the bond itself must be coloured"


def test_the_highlight_is_in_the_model_so_it_moves_with_the_bond(tmp_path):
    """Colouring the bond via atom.bondStyles rather than drawing a cylinder.

    An overlay sits at the equilibrium geometry while the stick swings through
    it. bondStyles is part of the model, so vibrate() copies it into every
    frame and the colour travels with the bond.
    """
    _, got = _run_js("C10H8", tmp_path, _BOND_STATE + """
        window.__cmp.selectBond(7);
        print(JSON.stringify(state(7)));
    """)
    r = json.loads(got)
    assert r["styled"] is not None, "the selected bond must carry a style"
    assert r["styled"]["color1"] == r["styled"]["color2"] == "#f2c14e"
    assert r["inFrames"] == r["styled"], \
        "the style must reach the animation frames, or it will not move"
    assert r["other"] is None, "a neighbouring bond must not be tinted"


def test_clicking_a_score_row_colours_that_bond(tmp_path):
    _, got = _run_js("C10H8", tmp_path, _BOND_STATE + """
        var C = window.__cmp, bl = document.getElementById('bondlist');
        bl.querySelectorAll('.brow')[7].fire('click');
        print(JSON.stringify({selected: C.selectedBond(), st: state(7)}));
    """)
    r = json.loads(got)
    assert r["selected"] == 7
    assert r["st"]["styled"] is not None
    assert r["st"]["inFrames"] is not None


def test_clicking_the_same_bond_again_clears_it(tmp_path):
    _, got = _run_js("C10H8", tmp_path, _BOND_STATE + """
        var C = window.__cmp;
        C.selectBond(4);
        var on = state(4).styled;
        C.selectBond(4);
        print(JSON.stringify({on: on, off: state(4).styled,
                              selected: C.selectedBond()}));
    """)
    r = json.loads(got)
    assert r["on"] is not None
    assert r["off"] is None, "the colour must be removed"
    assert r["selected"] is None


def test_selecting_a_bond_scrolls_its_row_into_view(tmp_path):
    """Scroll the list by the overflow only, never to the end.

    The first version used row.offsetTop, which is measured from the nearest
    POSITIONED ancestor -- and .viewer is position:sticky, so it included the
    panel header and the whole 3D canvas. Every row read as far below the
    window and the list slammed to the bottom on every selection. Rects are
    independent of the positioning context.
    """
    payload, got = _run_js("C10H8", tmp_path, """
        var C = window.__cmp, bl = document.getElementById('bondlist');
        function rect(i) { return bl.querySelectorAll('.brow')[i].getBoundingClientRect(); }
        function visible(i) {
          var r = rect(i), b = bl.getBoundingClientRect();
          return r.top >= b.top && r.bottom <= b.bottom;
        }
        var n = JSON.parse(PAYLOAD_TXT).bonds.length;
        var out = {maxScroll: n * 18 - bl.clientHeight,
                   offsetTopIsAncestorRelative: bl.querySelectorAll('.brow')[0].offsetTop};
        bl.scrollTop = 0;
        C.selectBond(3);
        out.near = {top: bl.scrollTop, visible: visible(3)};
        C.selectBond(3); C.selectBond(n - 8);
        out.mid = {top: bl.scrollTop, visible: visible(n - 8)};
        C.selectBond(n - 8); C.selectBond(0);
        out.back = {top: bl.scrollTop, visible: visible(0)};
        print(JSON.stringify(out));
    """)
    r = json.loads(got)
    assert r["offsetTopIsAncestorRelative"] > r["maxScroll"], \
        "test setup: offsetTop must exceed the scroll range, as it does under a sticky ancestor"
    assert r["near"]["top"] == 0 and r["near"]["visible"], \
        "a row already in view must not move the list"
    assert r["mid"]["visible"], "a distant row must be scrolled into view"
    assert 0 < r["mid"]["top"] < r["maxScroll"], \
        f"scrolled to {r['mid']['top']} of {r['maxScroll']} -- must be the overflow, not the end"
    assert r["back"]["top"] == 0 and r["back"]["visible"], "and scroll back up"


def test_atom_numbers_toggle_and_match_the_bond_labels(tmp_path):
    """Numbering must read the same as the bond labels (C1, H2 ...), or the
    per-bond panel and the picture disagree about which atom is which."""
    payload, got = _run_js("C10H8", tmp_path, """
        function fire(id) { document.getElementById(id).fire('change'); }
        var out = {off: CALLS.labels.map(function (l) { return l.text; })};
        document.getElementById('labels').checked = true; fire('labels');
        out.on = CALLS.labels.map(function (l) { return l.text; });
        document.getElementById('axes').checked = false; fire('axes');
        out.numbersOnly = CALLS.labels.map(function (l) { return l.text; });
        document.getElementById('labels').checked = false; fire('labels');
        out.bothOff = CALLS.labels.length;
        print(JSON.stringify(out));
    """)
    r = json.loads(got)
    n = len(payload["atoms"])
    assert r["off"] == ["x", "y", "z"], "numbers are off by default"
    assert len(r["on"]) == 3 + n, "axis letters plus one label per atom"
    assert len(r["numbersOnly"]) == n, "turning axes off must not drop the numbers"
    assert r["bothOff"] == 0

    expected = [s + str(i + 1) for i, s in enumerate(payload["atoms"])]
    assert r["numbersOnly"] == expected
    # and they agree with how bonds are named
    first_atom = payload["bond_labels"][0].split("-")[0]
    assert first_atom in expected


def test_viewer_controls_are_laid_out_in_two_rows():
    """Sliders on one line, toggles on the next -- flex wrapping put the first
    toggles up beside the sliders and the rest on their own line."""
    import re
    tpl = (ROOT / "app/templates/result.html").read_text()
    block = tpl[tpl.index('class="vctl"'):tpl.index('class="bondtab"')]
    rows = re.findall(r'<div class="vrow">(.*?)</div>', block, re.S)
    assert len(rows) == 2, f"expected two control rows, found {len(rows)}"
    for control in ("amp", "frm"):
        assert f'id="{control}"' in rows[0], f"{control} belongs on the slider row"
    for control in ("bonds", "axes", "arrows", "labels", "play"):
        assert f'id="{control}"' in rows[1], f"{control} belongs on the toggle row"


# ----------------------------------------------------------------------
# Bond selection reaches the comparison panel too.
# ----------------------------------------------------------------------
_CMP_BOND = """
    function colouredIn(bi) {
      var P = JSON.parse(PAYLOAD_TXT), b = P.bonds[bi];
      return CALLS.viewers.slice(1).map(function (v) {
        var m = v.models[v.models.length - 1];
        if (!m) return null;                       // an unused (hidden) slot
        var a = m.selectedAtoms({}), t = a[b[0]].bonds.indexOf(b[1]);
        return !!(a[b[0]].bondStyles && a[b[0]].bondStyles[t]);
      });
    }
"""


def test_compare_bond_table_rows_select_the_bond(tmp_path):
    payload, got = _run_js("C10H8", tmp_path, _CMP_BOND + """
        var C = window.__cmp;
        C.toggle(6); C.toggle(7);
        var cb = document.getElementById('cmp-bonds');
        var bl = document.getElementById('bondlist');
        var rows = cb.querySelectorAll('tr.cbrow');
        rows[5].fire('click');
        print(JSON.stringify({
          rows: rows.length,
          selected: C.selectedBond(),
          markedInCompare: /class="cbrow sel" data-b="5"/.test(cb._html),
          markedInPanel: /class="brow sel" data-b="5"/.test(bl._html),
          coloured: colouredIn(5)
        }));
    """)
    r = json.loads(got)
    assert r["rows"] == len(payload["bonds"]), "one row per bond"
    assert r["selected"] == 5
    assert r["markedInCompare"], "the compared table must mark the row"
    assert r["markedInPanel"], "and so must the single-mode panel -- one selection"
    active = [c for c in r["coloured"] if c is not None]
    assert active and all(active), "the bond must be coloured in every compared viewer"


def test_bonds_are_clickable_inside_a_compare_viewer(tmp_path):
    _, got = _run_js("C10H8", tmp_path, _CMP_BOND + """
        var C = window.__cmp;
        C.toggle(6); C.toggle(7);
        CALLS.cylinders[2].callback();              // a target in a compare view
        print(JSON.stringify({selected: C.selectedBond(), coloured: colouredIn(2)}));
    """)
    r = json.loads(got)
    assert r["selected"] == 2
    active = [c for c in r["coloured"] if c is not None]
    assert active and all(active)


def test_selecting_a_bond_before_comparing_survives_the_compare(tmp_path):
    """A bond chosen in the single view must still be coloured once a
    comparison opens -- the viewers are built after the selection."""
    _, got = _run_js("C10H8", tmp_path, _CMP_BOND + """
        var C = window.__cmp;
        C.selectBond(9);
        C.toggle(6); C.toggle(7);
        print(JSON.stringify({selected: C.selectedBond(), coloured: colouredIn(9)}));
    """)
    r = json.loads(got)
    assert r["selected"] == 9
    active = [c for c in r["coloured"] if c is not None]
    assert active and all(active)


def test_compare_rows_do_not_reuse_a_grid_class():
    """The compare table's rows must not carry a class that app.css declares
    `display:grid`.

    `.brow` is the side panel's bond row -- a div laid out as
    `grid-template-columns:74px 1fr 56px`. Put that class on a <tr> and the
    browser silently drops table layout: the cells become grid items in those
    three tracks, so mode columns stop lining up with their headers and a
    third mode's cell wraps onto a second grid row. Nothing throws; the table
    just renders wrong, which is why this is a static check.
    """
    css = (ROOT / "app" / "static" / "app.css").read_text()
    js = (ROOT / "app" / "static" / "app.js").read_text()

    # every class selector whose block sets display:grid
    grid = set()
    for sel, body in re.findall(r"([^{}]+)\{([^{}]*)\}", css):
        if re.search(r"display\s*:\s*grid", body):
            grid.update(re.findall(r"\.([A-Za-z0-9_-]+)", sel))
    assert "brow" in grid, "guard is only meaningful while .brow is a grid"

    # classes app.js puts on a <tr>
    on_tr = set()
    for frag in re.findall(r"<tr class=\\?[\"']([A-Za-z0-9_ -]+)", js):
        on_tr.update(frag.split())

    clash = on_tr & grid
    assert not clash, (
        f"class(es) {sorted(clash)} are used on a <tr> but declared display:grid "
        "in app.css -- this destroys table layout silently")


# ----------------------------------------------------------------------
# Per-bond sorting, and the centred bar geometry it is read against.
# ----------------------------------------------------------------------
_ORDER = """
    function order(el) {
      return (el._html.match(/data-b="(\\d+)"/g) || [])
             .map(function (m) { return +m.replace(/\\D/g, ''); });
    }
    var BL = document.getElementById('bondlist');
    var CB = document.getElementById('cmp-bonds');
    var H = document.querySelectorAll('.bhead .bsort');   // [pair, val]
"""


def _adamantane():
    """A molecule with enough bonds for ordering to be meaningful, and with
    labels that separate natural from lexicographic sorting (C6-C8 vs C6-C22)."""
    from app.core.parsers import parse_vsc
    d = parse_vsc((ROOT / "examples" / "C10H16.vsc").read_text())
    d.pop("title", None)
    d.pop("source", None)
    return analyse(**d, title="C10H16")


def test_bond_bars_grow_from_a_centred_zero(tmp_path):
    """Positive right, negative left, matching the compare panel -- so the bar
    can never exceed half the track."""
    _, got = _run_js(None, tmp_path,
                     "print(document.getElementById('bondlist')._html);",
                     payload=_adamantane())
    widths = [float(w) for w in re.findall(r'width:([0-9.]+)%', got)]
    kinds = re.findall(r'<i class="(pos|neg)"', got)
    assert widths, "no bars rendered"
    assert max(widths) <= 50.0 + 1e-9, f"bar exceeds half the track: {max(widths)}"
    assert abs(max(widths) - 50.0) < 1e-9, "the largest |s_AB| should fill its half"
    assert set(kinds) == {"pos", "neg"}, "both signs must be represented here"


def test_bond_sort_cycles_and_drives_both_windows(tmp_path):
    payload, got = _run_js(None, tmp_path, _ORDER + """
        var C = window.__cmp;
        C.toggle(6); C.toggle(7);
        var out = {def: order(BL), cmpDef: order(CB)};
        H[1].fire('click');  out.asc = order(BL); out.ascCmp = order(CB);
        out.stateAsc = C.bondSortState();
        H[1].fire('click');  out.desc = order(BL);
        H[1].fire('click');  out.off = order(BL); out.stateOff = C.bondSortState();
        print(JSON.stringify(out));
    """, payload=_adamantane())
    r = json.loads(got)
    n = len(payload["bonds"])
    assert r["def"] == list(range(n)), "default is bond order"
    assert r["cmpDef"] == r["def"]
    assert r["stateAsc"]["key"] == "val" and r["stateAsc"]["dir"] == 1
    assert r["asc"] != r["def"], "sorting must actually reorder"
    assert sorted(r["asc"]) == list(range(n)), "a sort is a permutation, not a rebuild"
    assert r["desc"] == r["asc"][::-1], "second click reverses exactly"
    assert r["off"] == r["def"], "third click returns to bond order"
    assert r["stateOff"]["key"] is None
    # the two windows must never disagree about order
    assert r["ascCmp"] == r["asc"], "compare table must follow the same ordering"


def test_bond_sort_by_name_is_natural_not_lexicographic(tmp_path):
    """C6-C8 before C6-C22. Plain string order puts "C22" first because '2'<'8',
    which scatters a numbered atom list."""
    payload, got = _run_js(None, tmp_path, _ORDER + """
        H[0].fire('click');
        print(JSON.stringify(order(BL)));
    """, payload=_adamantane())
    labels = [payload["bond_labels"][i] for i in json.loads(got)]
    assert labels.index("C6-C8") < labels.index("C6-C22"), labels[:12]
    assert labels.index("C4-C11") < labels.index("C4-C19")


def test_compare_column_sorts_by_that_mode(tmp_path):
    """Each mode's header orders the bonds by ITS values -- with three modes up,
    'which bonds does this one move' is the question being asked."""
    payload, got = _run_js(None, tmp_path, _ORDER + """
        var C = window.__cmp;
        C.toggle(6); C.toggle(7);
        var th = CB.querySelectorAll('th.bsort');
        th[2].fire('click');                       // the SECOND mode's column
        print(JSON.stringify({state: C.bondSortState(),
                              panel: order(BL), cmp: order(CB),
                              nth: th.length}));
    """, payload=_adamantane())
    r = json.loads(got)
    assert r["nth"] == 3, "bond column + one per compared mode"
    assert r["state"]["key"] == "val"
    assert r["state"]["mode"] == 7, "ordered by the mode whose column was clicked"
    assert r["panel"] == r["cmp"], "one ordering across both windows"
    # and it really is that mode's values, ascending
    row = [x for x in payload["vibrations"] if x["index"] == 7][0]
    vals = [row["bonds"][i]["s_AB"] for i in r["panel"]]
    assert vals == sorted(vals), "not ordered by the clicked mode's s_AB"


def test_sorting_keeps_the_selected_bond(tmp_path):
    """data-b carries the ORIGINAL bond index, so a sort must not renumber the
    selection out from under the viewer."""
    _, got = _run_js(None, tmp_path, _ORDER + """
        var C = window.__cmp;
        C.selectBond(5);
        H[1].fire('click');
        print(JSON.stringify({
          selected: C.selectedBond(),
          markedRow: /class="brow sel" data-b="5"/.test(BL._html),
          rows: order(BL).length
        }));
    """, payload=_adamantane())
    r = json.loads(got)
    assert r["selected"] == 5, "sorting changed which bond is selected"
    assert r["markedRow"], "the selected bond must stay marked after sorting"
