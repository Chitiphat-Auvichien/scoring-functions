"""The type-it-in builder, driven headlessly against the real builder.js.

It assembles the same .vsc text the upload path carries, so the server keeps a
single parser. These tests fill the form field by field and check the result
scores identically to the same modes read from the Gaussian log.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
JS = ROOT / "tests" / "js"
JSC = Path("/System/Library/Frameworks/JavaScriptCore.framework/"
           "Versions/A/Helpers/jsc")
sys.path.insert(0, str(ROOT))

from app.core.parsers import parse_gaussian_log, parse_connectivity  # noqa: E402
from app.core.pipeline import analyse  # noqa: E402


def _find_reference_data():
    for cand in (ROOT.parent / "data", ROOT.parent / "scoring-functions" / "data"):
        if (cand / "logs").is_dir():
            return cand
    return ROOT.parent / "data"


REF = _find_reference_data()

pytestmark = [
    pytest.mark.skipif(not JSC.exists(), reason="jsc not available"),
    pytest.mark.skipif(not REF.exists(), reason="reference data not present"),
]

_FILL = """
REG['b-atoms'].value='3'; REG['b-modes'].value='2'; __builder.build();
var G=[['O','0.038069','1.197522','0.000000'],
       ['H','-0.951724','1.378354','0.000000'],
       ['Cl','0.038069','-0.644619','0.000000']];
REG['b-geom'].querySelectorAll('.brow:not(.head)').forEach(function(r,i){
  r.querySelector('select').value=G[i][0];
  var n=r.querySelectorAll('input');
  n[0].value=G[i][1]; n[1].value=G[i][2]; n[2].value=G[i][3];
});
REG['b-bonds'].querySelectorAll('.add')[0].click();
var br=REG['b-bonds'].querySelectorAll('.brow.bond:not(.head)');
br[0].querySelectorAll('select')[1].value='2';
br[1].querySelectorAll('select')[1].value='3';
var D=[[['0.00088','0.85015','-0.00000'],['-0.16190','0.30478','-0.00000'],
        ['0.00426','-0.39765','0.00000']],
       [['-0.06138','0.01478','-0.00000'],['0.97981','-0.18968','0.00000'],
        ['-0.00017','-0.00130','0.00000']]];
REG['b-modes-box'].querySelectorAll('.modeblock').forEach(function(b,mi){
  b.querySelector('input.freq').value = mi===0?'667.6406':'3368.5426';
  b.querySelectorAll('.brow.disp:not(.head)').forEach(function(r,ai){
    var n=r.querySelectorAll('input');
    n[0].value=D[mi][ai][0]; n[1].value=D[mi][ai][1]; n[2].value=D[mi][ai][2];
  });
});
"""


def _run(probe):
    out = subprocess.run(
        [str(JSC), str(JS / "builder_dom.js"),
         str(ROOT / "app/static/builder.js"), "-e", probe],
        capture_output=True, text=True, cwd=ROOT)
    assert out.returncode == 0, out.stderr or out.stdout
    return out.stdout


def test_builder_constructs_the_requested_grid():
    got = _run("""
        REG['b-atoms'].value='5'; REG['b-modes'].value='3'; __builder.build();
        print(JSON.stringify({
          geom: REG['b-geom'].querySelectorAll('.brow:not(.head)').length,
          bonds: REG['b-bonds'].querySelectorAll('.brow.bond:not(.head)').length,
          blocks: REG['b-modes-box'].querySelectorAll('.modeblock').length,
          disp: REG['b-modes-box'].querySelectorAll('.brow.disp:not(.head)').length
        }));
    """)
    import json
    r = json.loads(got)
    assert r["geom"] == 5
    assert r["bonds"] == 1                 # one starter row
    assert r["blocks"] == 3
    assert r["disp"] == 15                 # 3 modes x 5 atoms


def test_typed_in_modes_score_exactly_like_the_gaussian_log():
    """Field-by-field entry must reproduce the source calculation."""
    vsc = _run(_FILL + "print(__builder.toVsc());")

    from app.core.parsers import parse_vsc
    v = parse_vsc(vsc)
    typed = analyse(v["atoms"], v["coords"], v["bonds"], v["modes"])

    g = parse_gaussian_log((REF / "logs/HOCl.log").read_text(errors="replace"))
    b = parse_connectivity((REF / "gjf/HOCl.com").read_text(errors="replace"),
                           len(g["atoms"]))
    direct = analyse(g["atoms"], g["coords"], b, [g["modes"][0], g["modes"][2]])

    assert typed["mode_set"] == "arbitrary"
    assert typed["atoms"] == ["O", "H", "Cl"]
    assert typed["bond_labels"] == direct["bond_labels"]
    for a, c in zip(typed["vibrations"], direct["vibrations"]):
        assert a["label"] == c["label"]
        assert a["frequency"] == pytest.approx(c["frequency"])
        for k in c["scores"]:
            assert a["scores"][k] == pytest.approx(c["scores"][k], abs=1e-9)


def test_builder_deduplicates_and_drops_self_bonds():
    got = _run("""
        REG['b-atoms'].value='3'; REG['b-modes'].value='1'; __builder.build();
        var add = REG['b-bonds'].querySelectorAll('.add')[0];
        add.click(); add.click();
        var br = REG['b-bonds'].querySelectorAll('.brow.bond:not(.head)');
        // 1-2, then 2-1 (same bond), then 3-3 (self)
        br[0].querySelectorAll('select')[0].value='1'; br[0].querySelectorAll('select')[1].value='2';
        br[1].querySelectorAll('select')[0].value='2'; br[1].querySelectorAll('select')[1].value='1';
        br[2].querySelectorAll('select')[0].value='3'; br[2].querySelectorAll('select')[1].value='3';
        print(__builder.toVsc());
    """)
    conn = got.split("[CONNECTIVITY]")[1].split("[MODES]")[0].strip().splitlines()
    assert conn == ["1  2"], conn


def test_blank_fields_become_zero():
    got = _run("""
        REG['b-atoms'].value='3'; REG['b-modes'].value='1'; __builder.build();
        print(__builder.toVsc());
    """)
    assert "0.000000" in got
    from app.core.parsers import parse_vsc
    v = parse_vsc(got)                     # must still be a legal file
    assert len(v["modes"]) == 1


# ----------------------------------------------------------------------
# Live preview of the structure being typed.
# ----------------------------------------------------------------------
_PREVIEW_SETUP = """
REG['b-atoms'].value='3'; REG['b-modes'].value='1'; __builder.build();
var G=[['O','0.038069','1.197522','0.000000'],
       ['H','-0.951724','1.378354','0.000000'],
       ['Cl','0.038069','-0.644619','0.000000']];
REG['b-geom'].querySelectorAll('.brow:not(.head)').forEach(function(r,i){
  r.querySelector('select').value=G[i][0];
  var n=r.querySelectorAll('input');
  n[0].value=G[i][1]; n[1].value=G[i][2]; n[2].value=G[i][3];
});
REG['b-bonds'].querySelectorAll('.add')[0].click();
var br=REG['b-bonds'].querySelectorAll('.brow.bond:not(.head)');
br[0].querySelectorAll('select')[1].value='2';
br[1].querySelectorAll('select')[1].value='3';
__builder.preview();
"""


def test_preview_renders_the_typed_structure():
    import json
    got = _run(_PREVIEW_SETUP + """
        var m = CALLS.models[0], a = m.selectedAtoms({});
        print(JSON.stringify({
          elems: a.map(function(x){ return x.elem; }),
          opts: m.opts,
          indices: a.map(function(x){ return x.index; }),
          bonds: a.map(function(x){ return x.bonds; }),
          style: CALLS.style,
          status: REG['b-vstat'].textContent
        }));
    """)
    r = json.loads(got)
    assert r["elems"] == ["O", "H", "Cl"]
    # the preview must show the bonds that will be scored, not a distance guess
    assert r["opts"]["assignBonds"] is False
    assert r["bonds"] == [[1, 2], [0], [0]]
    # drawBondSticks compares atom.index; without it no stick is ever drawn
    assert r["indices"] == [0, 1, 2]
    assert "stick" in r["style"] and "sphere" in r["style"]
    assert r["status"] == "3 atoms, 2 bonds"


def test_preview_reflects_edits_immediately():
    """A duplicated bond collapses, and the count says so."""
    got = _run(_PREVIEW_SETUP + """
        var br = REG['b-bonds'].querySelectorAll('.brow.bond:not(.head)');
        br[1].querySelectorAll('select')[1].value='2';     // 1-2 twice now
        __builder.preview();
        print(REG['b-vstat'].textContent);
    """)
    assert got.strip() == "3 atoms, 1 bond"


def test_preview_reports_incomplete_rows_rather_than_guessing():
    got = _run("""
        REG['b-atoms'].value='3'; REG['b-modes'].value='1'; __builder.build();
        var rows = REG['b-geom'].querySelectorAll('.brow:not(.head)');
        var n = rows[0].querySelectorAll('input');
        n[0].value='1.0'; n[1].value='0.0'; n[2].value='0.0';
        __builder.preview();
        print(REG['b-vstat'].textContent);
    """)
    assert "1 atom" in got and "incomplete" in got, got
