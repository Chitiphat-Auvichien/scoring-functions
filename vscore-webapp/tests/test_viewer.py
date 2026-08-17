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


def _run_js(stem, tmp_path, probe):
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
         str(data), str(JS / "dom_stub.js"), str(JS / "real_viewer.js"),
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
        for (var j = 0; j < 4; j++) a._listeners.change.call(a);
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
        pl.checked = false; pl._listeners.change.call(pl);
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
        lo._listeners.input.call(lo);
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
            lo._listeners.input.call(lo);
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
