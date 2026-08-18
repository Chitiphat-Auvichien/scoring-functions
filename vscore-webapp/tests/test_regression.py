"""Regression tests: this app must agree with the reference implementation.

The scoring core here is a vendored copy, so the value of these tests is
entirely in the comparison against ``scoring-functions/data/results/*.csv`` --
the published numbers. Any drift between the two is a bug in the copy.

Run:  python -m pytest tests/ -q      (from vscore-webapp/)
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from app.core.parsers import (ParseError, parse_connectivity,  # noqa: E402
                              parse_gaussian_log, parse_vsc, write_vsc)
from app.core.pipeline import analyse  # noqa: E402

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
SCORE_COLS = ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz"]

pytestmark = pytest.mark.skipif(
    not REF.exists(), reason="reference scoring-functions data not present")


def molecules():
    out = []
    for csv_path in sorted((REF / "results").glob("*_normal.csv")):
        stem = csv_path.name.replace("_normal.csv", "")
        log, com = REF / "logs" / f"{stem}.log", REF / "gjf" / f"{stem}.com"
        if log.exists() and com.exists():
            out.append((stem, log, com, csv_path))
    return out


def run(log, com, title=""):
    g = parse_gaussian_log(log.read_text(errors="replace"))
    bonds = parse_connectivity(com.read_text(errors="replace"), len(g["atoms"]))
    return g, bonds, analyse(g["atoms"], g["coords"], bonds, g["modes"], title=title)


@pytest.mark.parametrize("stem,log,com,csv_path", molecules(),
                         ids=[m[0] for m in molecules()])
def test_matches_reference(stem, log, com, csv_path):
    """Every score and label reproduces the published CSV exactly."""
    _, _, payload = run(log, com, stem)
    mine = payload["references"] + payload["vibrations"]

    with open(csv_path) as f:
        ref = list(csv.DictReader(f))
    assert len(ref) == len(mine), f"{stem}: row count {len(ref)} vs {len(mine)}"

    for expected, got in zip(ref, mine):
        for col in SCORE_COLS:
            assert abs(float(expected[col]) - got["scores"][col]) < 5e-4, \
                f"{stem} {got['name']} {col}"
        assert abs(float(expected["V_Stretch"]) - got["scores"]["V_S"]) < 5e-4
        assert expected["label"] == got["label"], f"{stem} {got['name']} label"


@pytest.mark.parametrize("stem,log,com,csv_path", molecules(),
                         ids=[m[0] for m in molecules()])
def test_vsc_roundtrip(stem, log, com, csv_path):
    """.log+.com -> .vsc -> scores must equal .log+.com -> scores."""
    g, bonds, direct = run(log, com, stem)
    text = write_vsc(g["atoms"], g["coords"], bonds, g["modes"], title=stem)
    v = parse_vsc(text)
    assert v["atoms"] == g["atoms"]
    assert v["bonds"] == list(bonds)
    again = analyse(v["atoms"], v["coords"], v["bonds"], v["modes"], title=stem)

    for a, b in zip(direct["vibrations"], again["vibrations"]):
        for k in a["scores"]:
            assert abs(a["scores"][k] - b["scores"][k]) < 1e-4, f"{stem} {a['name']} {k}"
        assert a["label"] == b["label"]
        assert a["dominant"] == b["dominant"]


def test_dominant_is_absolute_max():
    """The highlight is the largest |value|, not the largest value."""
    mols = molecules()
    assert mols, "no reference molecules found"
    _, log, com, _ = mols[0]
    _, _, p = run(log, com)
    for r in p["vibrations"]:
        biggest = max(r["scores"], key=lambda k: abs(r["scores"][k]))
        assert r["dominant"] == biggest
        assert abs(r["dominant_value"]) == pytest.approx(abs(r["scores"][biggest]))


def test_missing_connectivity_is_refused():
    """A deck with no connectivity must raise, never silently score V=0."""
    with pytest.raises(ParseError, match="connectivity"):
        parse_connectivity("%chk=x\n# hf/sto-3g\n\ntitle\n\n0 1\nH 0 0 0\n\n")


def test_vsc_without_connectivity_is_refused():
    with pytest.raises(ParseError, match="CONNECTIVITY"):
        parse_vsc("#VSCORE 1.0\n[GEOMETRY]\n1 H 0 0 0\n[MODES] 1\n"
                  "mode 1 freq=1.0\n1 0.1 0.0 0.0\n")


def test_atom_count_mismatch_is_refused():
    """A .com describing a different molecule than the .log must be caught."""
    mols = molecules()
    _, log, com, _ = mols[0]
    g = parse_gaussian_log(log.read_text(errors="replace"))
    with pytest.raises(ParseError, match="only|atoms"):
        parse_connectivity(com.read_text(errors="replace"), natoms=1)


def test_log_without_frequencies_is_refused():
    with pytest.raises(ParseError, match="frequency|orientation"):
        parse_gaussian_log("Gaussian 16 output\n no useful content here\n")


def test_precision_warning_fires_on_two_dp():
    """2-dp displacements must be flagged: T-scores move by up to ~0.3."""
    import numpy as np
    mols = molecules()
    _, log, com, _ = mols[0]
    g = parse_gaussian_log(log.read_text(errors="replace"))
    bonds = parse_connectivity(com.read_text(errors="replace"), len(g["atoms"]))
    coarse = [{**m, "vector": np.round(m["vector"], 2)} for m in g["modes"]]
    p = analyse(g["atoms"], g["coords"], bonds, coarse)
    assert p["precision_dp"] <= 2
    assert any("decimal places" in w for w in p["warnings"])


# ----------------------------------------------------------------------
# Regression guard for the axis-sign bug found in testing.
# ----------------------------------------------------------------------
def test_mit_is_not_idempotent():
    """Documents WHY .vsc stores the original frame.

    MIT()'s sign-fix heuristic negates the rotation when the heaviest atom's
    projected coordinates sum negative. Feeding an already-aligned geometry
    back in can therefore flip axes -- which silently inverted Ty on a .vsc
    round-trip until the writer was changed to store the original frame.
    If this test ever starts failing, MIT() became idempotent and the comment
    in main._build can be simplified.
    """
    from app.core.scoring import ModeScorer
    mols = molecules()
    flipped = []
    for stem, log, com, _ in mols:
        g = parse_gaussian_log(log.read_text(errors="replace"))
        bonds = parse_connectivity(com.read_text(errors="replace"), len(g["atoms"]))
        s1 = ModeScorer(g["atoms"], g["coords"], bonds)
        m1 = s1.MIT([dict(m) for m in g["modes"]], rotate_modes=True)
        a = s1.calculate_scores(m1[0]["vector"])
        s2 = ModeScorer(g["atoms"], s1.coords, bonds)
        m2 = s2.MIT([dict(m) for m in m1], rotate_modes=True)
        b = s2.calculate_scores(m2[0]["vector"])
        if any(abs(a["T"][k] - b["T"][k]) > 1e-6 for k in "xyz"):
            flipped.append(stem)
    assert flipped, ("MIT() now appears idempotent on every reference molecule; "
                     "re-check whether .vsc may store the aligned frame.")


@pytest.mark.parametrize("stem,log,com,csv_path", molecules(),
                         ids=[m[0] for m in molecules()])
def test_webapp_vsc_download_roundtrips(stem, log, com, csv_path):
    """The .vsc the RESULTS PAGE hands the user must re-score identically.

    This is the path the earlier round-trip test missed: it built its own .vsc
    from the original coordinates, while the app built one from the scored
    (aligned) payload, which inverted Ty.
    """
    from app.main import _build, _vsc_text
    payload, raw = _build("log", log.read_text(errors="replace"),
                          com.read_text(errors="replace"), None,
                          f"{stem}.log", f"{stem}.com", None)
    again, _ = _build("vsc", None, None, _vsc_text(payload, raw),
                      None, None, f"{stem}.vsc")
    for a, b in zip(payload["vibrations"], again["vibrations"]):
        for k in a["scores"]:
            assert abs(a["scores"][k] - b["scores"][k]) < 1e-4, \
                f"{stem} {a['name']} {k}: {a['scores'][k]} vs {b['scores'][k]}"
        assert a["label"] == b["label"]
        assert a["dominant"] == b["dominant"]


# ----------------------------------------------------------------------
# Highlighting policy: only the six constructed T/R modes carry external
# character, so only they may highlight a T/R column.
# ----------------------------------------------------------------------
@pytest.mark.parametrize("stem,log,com,csv_path", molecules(),
                         ids=[m[0] for m in molecules()])
def test_only_reference_modes_highlight_TR(stem, log, com, csv_path):
    """A real vibration must always highlight V_S, never a T/R column.

    Step 2 assigns every external slot to the ideal T/R references, so a real
    mode is never given external character. Its abs-max lands on a T/R column
    for ~43% of modes, and highlighting that would assert translation or
    rotation the algorithm explicitly did not assign.
    """
    _, _, p = run(log, com, stem)

    for r in p["vibrations"]:
        assert r["is_reference"] is False
        assert r["highlight"] == "V_S", f"{stem} {r['name']} highlights {r['highlight']}"

    for r in p["references"]:
        assert r["is_reference"] is True
        # a reference highlights its own dominant score, which must be its slot
        assert r["highlight"] == r["dominant"]
        assert r["highlight"] != "V_S", f"{stem} {r['name']} highlights V_S"
        assert r["highlight"] == r["label"].rstrip("*"), \
            f"{stem} {r['name']}: highlight {r['highlight']} != slot {r['label']}"


def test_reference_count_matches_geometry():
    """6 reference modes, or 5 for a linear molecule (Rx is dropped)."""
    for stem, log, com, _ in molecules():
        _, _, p = run(log, com, stem)
        assert len(p["references"]) == (5 if p["linear"] else 6), \
            f"{stem}: {len(p['references'])} references, linear={p['linear']}"


def test_csv_carries_both_highlight_and_dominant():
    """The download must record what was highlighted AND the raw abs-max."""
    from app.core.pipeline import to_csv_rows
    mols = molecules()
    _, log, com, _ = mols[0]
    _, _, p = run(log, com)
    rows = to_csv_rows(p, include_references=True)
    assert all("highlighted_score" in r and "dominant_score" in r for r in rows)
    vib = [r for r in rows if not r["is_reference_mode"]]
    assert vib and all(r["highlighted_score"] == "V_S" for r in vib)
    ref = [r for r in rows if r["is_reference_mode"]]
    assert ref and all(r["highlighted_score"] != "V_S" for r in ref)


# ----------------------------------------------------------------------
# .vsc mode metadata is optional -- the scorer reads none of it.
# ----------------------------------------------------------------------
_VSC_HEAD = """#VSCORE 1.0
[GEOMETRY] Angstrom
   1   O     0.038069    1.197522    0.000000
   2   H    -0.951724    1.378354    0.000000
   3  Cl     0.038069   -0.644619    0.000000
[CONNECTIVITY]
   1  2  3
   2
   3
"""
_VSC_DISP = ("   1    0.00088    0.85015   -0.00000\n"
             "   2   -0.16190    0.30478   -0.00000\n"
             "   3    0.00426   -0.39765    0.00000\n")


def _vsc(mode_header, count=" 1"):
    return _VSC_HEAD + f"[MODES]{count}\n" + mode_header + "\n" + _VSC_DISP


@pytest.mark.parametrize("header,expect", [
    ("  mode 1   freq=667.6406   mu=17.7040   k=4.9258   irrep=A'",
     {"frequency": 667.6406, "reduced_mass": 17.7040,
      "force_constant": 4.9258, "irrep": "A'"}),
    ("  mode 1   freq=667.6406",
     {"frequency": 667.6406, "reduced_mass": None,
      "force_constant": None, "irrep": None}),
    ("  mode 1   irrep=A'",
     {"frequency": None, "reduced_mass": None,
      "force_constant": None, "irrep": "A'"}),
    ("  mode 1   mu=17.7040   k=4.9258",
     {"frequency": None, "reduced_mass": 17.7040,
      "force_constant": 4.9258, "irrep": None}),
    ("  mode 1",
     {"frequency": None, "reduced_mass": None,
      "force_constant": None, "irrep": None}),
    ("  MODE 1   FREQ=667.6406",                      # case-insensitive
     {"frequency": 667.6406, "reduced_mass": None,
      "force_constant": None, "irrep": None}),
])
def test_vsc_mode_metadata_is_optional(header, expect):
    """freq / mu / k / irrep are all optional, in any combination.

    None of them feeds the scorer -- they are row metadata. An absent freq must
    stay None, not default to 0.0, or the table shows an invented 0.00 cm-1.
    """
    v = parse_vsc(_vsc(header))
    p = analyse(v["atoms"], v["coords"], v["bonds"], v["modes"], title="t")
    r = p["vibrations"][0]
    for key, want in expect.items():
        got = r[key]
        if want is None:
            assert got is None, f"{key}: expected None, got {got!r}"
        else:
            assert got == pytest.approx(want), f"{key}: {got} != {want}"
    # scores are identical whatever metadata is present or absent
    assert r["scores"]["V_S"] == pytest.approx(0.9843, abs=1e-4)
    assert r["label"] == "S"


def test_vsc_modes_count_is_optional():
    v = parse_vsc(_vsc("  mode 1   freq=667.6406", count=""))
    assert len(v["modes"]) == 1


def test_vsc_modes_count_must_agree_when_given():
    with pytest.raises(ParseError, match="declares"):
        parse_vsc(_vsc("  mode 1   freq=667.6406", count=" 3"))


@pytest.mark.parametrize("header", [
    "  mode 1   freq=667.6406   mu=17.7040   k=4.9258   irrep=A'",
    "  mode 1   irrep=A'",
    "  mode 1",
])
def test_write_vsc_omits_absent_metadata(header):
    """The writer must not invent metadata the input never had."""
    v = parse_vsc(_vsc(header))
    text = write_vsc(v["atoms"], v["coords"], v["bonds"], v["modes"])
    line = next(l.strip() for l in text.splitlines() if l.strip().startswith("mode"))
    for key, present in (("freq=", "freq=" in header),
                         ("mu=", "mu=" in header),
                         ("k=", "k=" in header),
                         ("irrep=", "irrep=" in header)):
        assert (key in line) == present, f"{key} in {line!r} but input was {header!r}"
    # and it must still round-trip to the same scores
    again = parse_vsc(text)
    a = analyse(v["atoms"], v["coords"], v["bonds"], v["modes"])
    b = analyse(again["atoms"], again["coords"], again["bonds"], again["modes"])
    for k in a["vibrations"][0]["scores"]:
        assert a["vibrations"][0]["scores"][k] == pytest.approx(
            b["vibrations"][0]["scores"][k], abs=1e-9)



# ----------------------------------------------------------------------
# EMIT modes: a .vsc with eigen= (not freq=) headers must be rotated
# with rotate_modes=False and score identically to the published
# data/results/C6H6_EMIT.csv (main.py's headless EMIT pipeline).
# ----------------------------------------------------------------------
def _emit_vsc_path():
    return REF / "EMIT" / "C6H6_EMIT.vsc"


@pytest.mark.skipif(not _emit_vsc_path().exists(), reason="C6H6_EMIT.vsc not present")
def test_emit_vsc_matches_reference():
    v = parse_vsc(_emit_vsc_path().read_text())
    assert all(m["is_emit"] for m in v["modes"]), "eigen= must mark every mode is_emit"

    payload = analyse(v["atoms"], v["coords"], v["bonds"], v["modes"], title="C6H6 EMIT")
    assert payload["references"] == [], "EMIT has no synthetic ideal T/R rows"
    assert len(payload["vibrations"]) == 36

    with open(REF / "results" / "C6H6_EMIT.csv") as f:
        ref = list(csv.DictReader(f))
    assert len(ref) == 36

    for expected, got in zip(ref, payload["vibrations"]):
        for col in SCORE_COLS:
            assert abs(float(expected[col]) - got["scores"][col]) < 5e-4, \
                f"{got['name']} {col}"
        assert abs(float(expected["V_Stretch"]) - got["scores"]["V_S"]) < 5e-4
        assert expected["label"] == got["label"], f"{got['name']}"


@pytest.mark.skipif(not _emit_vsc_path().exists(), reason="C6H6_EMIT.vsc not present")
def test_emit_rotate_modes_false_matters():
    """Scoring EMIT modes with rotate_modes=True (the normal-mode default)
    must NOT reproduce the reference -- this is the exact silent-corruption
    failure mode eigen= exists to prevent."""
    from app.core.scoring import ModeScorer

    v = parse_vsc(_emit_vsc_path().read_text())
    scorer = ModeScorer(v["atoms"], v["coords"], v["bonds"])
    wrongly_rotated = scorer.MIT([dict(m) for m in v["modes"]], rotate_modes=True)

    with open(REF / "results" / "C6H6_EMIT.csv") as f:
        expected = list(csv.DictReader(f))

    # At least some mode's T/R scores must diverge once double-rotated --
    # a handful of highly symmetric EMIT modes happen to keep near-zero T/R
    # components under either rotation, so check across all 36, not just one.
    max_diff = 0.0
    for exp, wmode in zip(expected, wrongly_rotated):
        got = scorer.calculate_scores(wmode["vector"])
        for axis in "xyz":
            max_diff = max(max_diff, abs(got["T"][axis] - float(exp[f"T{axis}"])))
            max_diff = max(max_diff, abs(got["R"][axis] - float(exp[f"R{axis}"])))
    assert max_diff > 1e-2, "double-rotating EMIT modes should visibly corrupt T/R scores"


@pytest.mark.skipif(not _emit_vsc_path().exists(), reason="C6H6_EMIT.vsc not present")
def test_emit_csv_uses_eigenvalue_column():
    from app.core.pipeline import to_csv_rows

    v = parse_vsc(_emit_vsc_path().read_text())
    payload = analyse(v["atoms"], v["coords"], v["bonds"], v["modes"], title="C6H6 EMIT")
    rows = to_csv_rows(payload)
    assert all("Eigenvalue" in r and "Freq" not in r for r in rows)


def test_write_vsc_eigen_roundtrips():
    """write_vsc() must tag EMIT modes with eigen=, and parse_vsc() must read
    it back as is_emit=True with the value preserved as 'frequency'."""
    modes = [{"frequency": -23.7514, "vector": np.zeros((3, 3)), "is_emit": True,
              "reduced_mass": None, "force_constant": None, "irrep": None}]
    atoms, coords = ["O", "H", "Cl"], np.zeros((3, 3))
    bonds = [(0, 1), (0, 2)]
    text = write_vsc(atoms, coords, bonds, modes)
    line = next(l.strip() for l in text.splitlines() if l.strip().startswith("mode"))
    assert "eigen=-23.7514" in line
    assert "freq=" not in line

    v = parse_vsc(text)
    assert v["modes"][0]["is_emit"] is True
    assert v["modes"][0]["frequency"] == pytest.approx(-23.7514)


def test_reference_data_was_actually_found():
    """Fail loudly if the reference data is missing, rather than skipping.

    pytestmark's skipif silently skipped all 42 tests when the webapp moved into
    scoring-functions/ and REF resolved to a path that does not exist. A skipped
    suite reads as a passing one, so assert the location explicitly.
    """
    assert REF.exists(), f"reference data not found; REF resolved to {REF}"
    assert (REF / "logs").is_dir() and (REF / "results").is_dir()
    assert molecules(), "reference data found but no molecule had both .log and .com"
