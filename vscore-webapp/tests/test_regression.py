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
# HOCl is 3 atoms, so a valid file must carry 3N-6 = 3 modes (or 3N = 9).
# The metadata under test varies on mode 1; modes 2 and 3 just make the count
# legal, since a partial set is refused rather than guessed at.
_VSC_FILLER = ("  mode 2\n"
               "   1   -0.04924   -0.02717    0.00000\n"
               "   2    0.17600    0.97748   -0.00000\n"
               "   3    0.01906   -0.02189    0.00000\n"
               "  mode 3\n"
               "   1   -0.05672    0.02347   -0.00000\n"
               "   2    0.97331   -0.21988    0.00000\n"
               "   3   -0.00108   -0.00035    0.00000\n")


def _vsc(mode_header, count=" 3"):
    return (_VSC_HEAD + f"[MODES]{count}\n" + mode_header + "\n"
            + _VSC_DISP + _VSC_FILLER)


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
    assert len(v["modes"]) == 3


def test_vsc_modes_count_must_agree_when_given():
    with pytest.raises(ParseError, match="declares"):
        parse_vsc(_vsc("  mode 1   freq=667.6406", count=" 7"))


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



def test_reference_data_was_actually_found():
    """Fail loudly if the reference data is missing, rather than skipping.

    pytestmark's skipif silently skipped all 42 tests when the webapp moved into
    scoring-functions/ and REF resolved to a path that does not exist. A skipped
    suite reads as a passing one, so assert the location explicitly.
    """
    assert REF.exists(), f"reference data not found; REF resolved to {REF}"
    assert (REF / "logs").is_dir() and (REF / "results").is_dir()
    assert molecules(), "reference data found but no molecule had both .log and .com"


# ----------------------------------------------------------------------
# .vsc mode_set: 3N (T/R supplied) vs 3N-6 (vibrations only)
# ----------------------------------------------------------------------
def _molecule(stem="HCOOH"):
    for s, log, com, _ in molecules():
        if s == stem:
            return parse_gaussian_log(log.read_text(errors="replace")), \
                   parse_connectivity(com.read_text(errors="replace"),
                                      len(parse_gaussian_log(
                                          log.read_text(errors="replace"))["atoms"]))
    pytest.skip(f"{stem} not in the reference set")


def _complete_3n(g, bonds):
    """The full 3N set: ideal T/R built in the PRINCIPAL frame, then the
    vibrations, all in that same frame. Constructing T/R about non-principal
    axes instead makes each one a mixture of Rx/Ry/Rz, which no single external
    slot can claim -- see test_3n_recovery_degrades_off_principal_axes."""
    from app.core.scoring import ModeScorer
    sc = ModeScorer(g["atoms"], g["coords"], bonds)
    rot = sc.MIT([dict(m) for m in g["modes"]], rotate_modes=True)
    refs = sc.construct_T() + sc.construct_R()
    modes = [{"frequency": None, "vector": m["vector"]} for m in refs]
    modes += [{"frequency": m["frequency"], "vector": m["vector"],
               "irrep": m.get("irrep")} for m in rot]
    return sc.coords, modes


def test_3n_finds_the_six_TR_modes_without_constructing_any():
    g, bonds = _molecule()
    coords, modes = _complete_3n(g, bonds)
    p = analyse(g["atoms"], coords, bonds, modes, mode_set="3n")

    assert p["references"] == [], "3N must not construct reference modes"
    assert p["n_modes"] == 3 * len(g["atoms"])
    ext = [r for r in p["vibrations"] if r["label"][0] in "TR"]
    assert len(ext) == 6, f"{len(ext)} external slots filled, expected 6"
    # the six supplied T/R modes are the ones claimed, and cleanly
    assert [r["name"] for r in ext] == [f"Vib {i}" for i in range(1, 7)]
    assert all(not r["label"].endswith("*") for r in ext), \
        "a clean 3N set should not need starred externals"
    assert {r["label"] for r in ext} == {"Tx", "Ty", "Tz", "Rx", "Ry", "Rz"}
    assert p["warnings"] == []


def test_3n_and_3n6_agree_on_the_real_vibrations():
    """Both paths must classify the 3N-6 physical modes identically."""
    g, bonds = _molecule()
    a = analyse(g["atoms"], g["coords"], bonds, g["modes"], mode_set="3n-6")
    coords, modes = _complete_3n(g, bonds)
    c = analyse(g["atoms"], coords, bonds, modes, mode_set="3n")
    vib = [r for r in c["vibrations"] if r["label"][0] not in "TR"]

    assert len(vib) == len(a["vibrations"])
    for x, y in zip(a["vibrations"], vib):
        assert x["label"] == y["label"]
        assert x["frequency"] == pytest.approx(y["frequency"], abs=1e-3)
        # magnitudes only: re-aligning an already-principal geometry can flip
        # axis signs (MIT is not idempotent), which is a frame convention and
        # does not affect the labels -- the assignment uses |score|.
        for k in x["scores"]:
            assert abs(x["scores"][k]) == pytest.approx(abs(y["scores"][k]), abs=2e-3)


def test_3n_recovery_degrades_off_principal_axes():
    """Documents the limitation the 3N warning exists for.

    T/R constructed about non-principal axes become mixtures of Rx/Ry/Rz once
    MIT rotates into the principal frame. No single slot score is then high
    enough to win, so a real vibration can take the slot and a genuine rotation
    falls through to Step 4.
    """
    from app.core.scoring import ModeScorer
    g, bonds = _molecule()
    sc = ModeScorer(g["atoms"], g["coords"], bonds)      # COM only, NOT aligned
    refs = sc.construct_T() + sc.construct_R()
    modes = [{"frequency": None, "vector": m["vector"]} for m in refs]
    modes += [{"frequency": m["frequency"], "vector": m["vector"]} for m in g["modes"]]

    p = analyse(g["atoms"], g["coords"], bonds, modes, mode_set="3n")
    ext = [r["name"] for r in p["vibrations"] if r["label"][0] in "TR"]
    assert ext != [f"Vib {i}" for i in range(1, 7)], (
        "off-axis T/R now recovers perfectly; the 3N warning text and this "
        "test's rationale should be revisited")
    assert any("mixed external" in n for n in p["notes"]), \
        "the user must be told the recovery was not clean"
    assert not any("mixed external" in w for w in p["warnings"]), \
        "a starred label is information, not a fault -- it belongs in notes"


def test_mode_set_is_detected_from_the_counts():
    """3N vs 3N-6 is never ambiguous, so the user does not declare it."""
    g, bonds = _molecule()
    vib = analyse(g["atoms"], g["coords"], bonds, g["modes"])
    assert vib["mode_set"] == "3n-6"
    assert len(vib["references"]) == 6

    coords, modes = _complete_3n(g, bonds)
    full = analyse(g["atoms"], coords, bonds, modes)
    assert full["mode_set"] == "3n"
    assert full["references"] == []


def test_a_count_matching_neither_is_refused():
    """A partial set cannot be classified: the external slots would be filled
    from whatever happened to be supplied."""
    g, bonds = _molecule()                      # 5 atoms -> 3N=15, 3N-6=9
    for n in (1, 5, 8, 10, 14, 16):
        with pytest.raises(ParseError, match="neither 3N"):
            analyse(g["atoms"], g["coords"], bonds, g["modes"][:1] * n)


def test_declared_mode_set_must_match_the_counts():
    g, bonds = _molecule()
    with pytest.raises(ParseError, match="but 3n was requested"):
        analyse(g["atoms"], g["coords"], bonds, g["modes"], mode_set="3n")


def test_linear_molecule_uses_3n_minus_5():
    g, bonds = _molecule("CO2")
    p = analyse(g["atoms"], g["coords"], bonds, g["modes"], mode_set="3n-6")
    assert p["linear"] is True
    assert len(p["references"]) == 5, "linear molecules get 5 references, not 6"
    assert p["n_modes"] == 3 * len(g["atoms"]) - 5


def test_unknown_mode_set_raises():
    g, bonds = _molecule()
    with pytest.raises(ValueError, match="mode_set"):
        analyse(g["atoms"], g["coords"], bonds, g["modes"], mode_set="every")


def test_viewer_frame_is_the_score_frame():
    """The geometry and vectors the viewer draws must be the frame the scores
    are defined in -- otherwise the drawn axes would not mean Tx/Ty/Tz.

    Checked by recomputing s[T] straight from the payload the browser receives
    and comparing with the reported score, and by confirming the ideal Tx
    reference displaces purely along +x.
    """
    import numpy as np
    for stem, log, com, _ in molecules()[:6]:
        _, _, p = run(log, com, stem)
        for r in p["vibrations"][:3]:
            V = np.array(r["vector"])
            L = np.sqrt((V ** 2).sum(axis=1))
            mask = L > 1e-6
            U = np.zeros_like(V)
            U[mask] = V[mask] / L[mask, None]
            tx, ty, tz = U.sum(axis=0) / len(V)
            for got, key in ((tx, "Tx"), (ty, "Ty"), (tz, "Tz")):
                assert got == pytest.approx(r["scores"][key], abs=2e-4), (
                    f"{stem} {r['name']} {key}: payload gives {got}, "
                    f"score says {r['scores'][key]} -- viewer frame != score frame")

        by_name = {r["name"]: r for r in p["references"]}
        for name, axis in (("Tx", 0), ("Ty", 1), ("Tz", 2)):
            if name not in by_name:
                continue
            V = np.array(by_name[name]["vector"])
            assert np.allclose(V[:, axis], V[0, axis]), f"{stem} {name} not uniform"
            for other in (0, 1, 2):
                if other != axis:
                    assert np.allclose(V[:, other], 0.0, atol=1e-6), \
                        f"{stem} {name} has a component off the {name[-1]} axis"


def test_frame_rotation_is_reported():
    """The angle between the input orientation and the score frame is surfaced,
    so a molecule that looks turned in the viewer is explained rather than
    mistaken for a mismatch between the axes and the scores."""
    for stem, log, com, _ in molecules()[:5]:
        _, _, p = run(log, com, stem)
        assert "frame_rotation_deg" in p
        assert 0.0 <= p["frame_rotation_deg"] <= 180.0


def test_rotation_labels_match_the_physical_axis_in_the_principal_frame():
    """A mode labelled Rx must actually rotate about x, for T/R built in the
    principal frame.

    s[R_Q] measures whether the tangential SENSE is consistent about Q, not how
    much of the rotation axis lies along Q. For a planar molecule an
    off-principal rotation can therefore saturate at |1.000| on the wrong axis:
    HCOOH's rotation about [0.956, -0.294, 0] scored s[Ry] = -1.000 (all five
    sign terms agreed) against s[Rx] = +0.600 (one disagreed), and was labelled
    Ry while visibly rotating about x. Built in the principal frame the labels
    and the physical axes agree exactly.
    """
    import numpy as np
    from app.core.scoring import ModeScorer
    g, bonds = _molecule()
    coords, modes = _complete_3n(g, bonds)
    p = analyse(g["atoms"], coords, bonds, modes, mode_set="3n")
    G = np.array(p["geometry"])

    for r in p["vibrations"]:
        if not (r["label"].startswith("R") and len(r["label"]) == 2):
            continue
        D = np.array(r["vector"])
        A, b = [], []
        for rr, d in zip(G, D):
            A.append([[0, rr[2], -rr[1]], [-rr[2], 0, rr[0]], [rr[1], -rr[0], 0]])
            b.append(d)
        w, *_ = np.linalg.lstsq(np.vstack(A), np.concatenate(b), rcond=None)
        n = np.linalg.norm(w)
        assert n > 1e-9, f"{r['name']} labelled {r['label']} is not a rotation"
        u = w / n
        axis = "xyz"[int(np.argmax(np.abs(u)))]
        assert r["label"] == "R" + axis, (
            f"{r['name']} labelled {r['label']} but rotates about {axis} "
            f"(axis {np.round(u, 3).tolist()})")
        assert abs(u["xyz".index(axis)]) > 0.99, \
            f"{r['name']} axis {np.round(u, 3).tolist()} is not a clean {axis} rotation"
