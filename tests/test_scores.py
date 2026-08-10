"""Golden-reference regression tests for the scoring engine.

Run from ``Github/scoring-functions/``:
    py -m pytest tests/            (with pytest)
    py tests/test_scores.py       (standalone; no pytest needed)

Golden values: the published water table (tab:water, consensus s[R]) and the
benzene-EMIT score-level targets. Invariants: Sum s_AB == s[V_S] and score ranges.
These freeze the engine's numerical behaviour so regressions are caught immediately.

**2026-07-07 basename update:** the old `data/logs/water.log`/`data/gjf/water.gjf`
and `data/logs/benzene.log`/`data/gjf/benzene.com` were replaced (not just
renamed) by the author when `data/mol_list_method.csv`'s 72-molecule roster
was finalized -- `H2O` and `C6H6` are the new,
finalized on-disk basenames (`data/EMIT/*_EMIT.txt` renamed to match, same
session). Benzene's EMIT-derived numbers below are numerically UNCHANGED
(re-verified against a live run) -- the rename there was content-preserving.
Water's numbers DID change: `H2O.log` is a genuinely different
(corrected) calculation from the old `water.log`, not the same water calc
under a new name (confirmed: old `water.log`'s engine-parsed frequencies
were 1628.029/3887.192/4005.506 cm-1, mismatching even the literature/
data_score.csv values; the new file's 1722.457/3501.507/3660.797 cm-1 match
data_score.csv). `test_water_tab_water`'s pinned values are updated to the
new file's real output (re-derived from a live run, not fudged) -- these no
longer necessarily match the manuscript's currently-typeset `tab:water`
(that reconciliation is separate, downstream manuscript work).
"""
import os
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from src.parser import GaussianParser, EMITParser   # noqa: E402
from src.scoring import ModeScorer                  # noqa: E402

TOL = 5e-4  # 3 decimal places


def _normal_table(mol):
    """Mirror main.py for normal modes: parse -> MIT(rotate) -> +ideal T/R -> score."""
    gp = GaussianParser(os.path.join(ROOT, "data", "logs", f"{mol}.log"))
    data = gp.parse(parse_modes=True)
    scorer = ModeScorer(data["atoms"], data["coords"], data["bonds"])
    rotated = scorer.MIT(data["modes"], rotate_modes=True)
    labels = ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz"] + [f"Vib{i+1}" for i in range(len(rotated))]
    modes = scorer.construct_T() + scorer.construct_R() + rotated
    table = {}
    for lbl, m in zip(labels, modes):
        sc = scorer.calculate_scores(m["vector"])
        table[lbl] = {"Tx": sc["T"]["x"], "Ty": sc["T"]["y"], "Tz": sc["T"]["z"],
                      "Rx": sc["R"]["x"], "Ry": sc["R"]["y"], "Rz": sc["R"]["z"], "V": sc["V"]}
    return scorer, table


def _emit_table(mol):
    """Mirror main.py for EMIT modes: parse geometry -> EMIT modes -> MIT(no mode rotation) -> score."""
    gp = GaussianParser(os.path.join(ROOT, "data", "logs", f"{mol}.log"))
    data = gp.parse(parse_modes=False)
    ep = EMITParser(os.path.join(ROOT, "data", "EMIT", f"{mol}_EMIT.txt"), len(data["atoms"]))
    emit = ep.parse()
    scorer = ModeScorer(data["atoms"], data["coords"], data["bonds"])
    rotated = scorer.MIT(emit, rotate_modes=False)  # EMIT modes already in principal axes
    table = {}
    for i, m in enumerate(rotated):
        sc = scorer.calculate_scores(m["vector"])
        table[f"EMIT {i+1}"] = {"Tx": sc["T"]["x"], "Ty": sc["T"]["y"], "Tz": sc["T"]["z"],
                                "Rx": sc["R"]["x"], "Ry": sc["R"]["y"], "Rz": sc["R"]["z"], "V": sc["V"]}
    return scorer, table


def test_water_tab_water():
    """Score-level water golden values from the finalized roster basename
    H2O (consensus s[R]). Re-derived from a live run against the
    new file -- these are NOT the manuscript's currently-typeset tab:water
    numbers (see module docstring: the old water.log was a different,
    mismatched calculation); manuscript reconciliation is separate,
    downstream work."""
    _, t = _normal_table("H2O")
    for ext in ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz"]:
        assert abs(t[ext][ext] - 1.0) < TOL, f"external {ext} should be 1.000"
    assert abs(t["Tx"]["Rz"] - 0.0411) < TOL
    assert abs(t["Tz"]["Rx"] - (-0.3333)) < TOL
    assert abs(t["Vib1"]["V"] - 0.1343) < TOL         # bend (sigma)
    assert abs(t["Vib2"]["V"] - 0.9966) < TOL         # nu_s
    assert abs(t["Vib3"]["V"] - 0.9984) < TOL         # nu_as
    # 2026-08-10 G16 promotion: Vib3's Tx/Rz pair flipped sign together
    # (old G09 H2O.log: Tx=-0.1963, Rz=-0.2959; new G16 H2O.log: Tx=+0.1963,
    # Rz=+0.2959) -- an arbitrary principal-axis/eigensolver sign convention
    # (CLAUDE.md: "sign-of-frame is arbitrary; classification uses |score|"),
    # not a logic bug -- both components flipped together and magnitudes are
    # bit-identical to 4 dp, confirming a coherent frame flip.
    assert abs(t["Vib3"]["Tx"] - 0.1963) < TOL
    assert abs(t["Vib3"]["Rz"] - 0.2959) < TOL


def test_co2_linear():
    """Linear molecule: molecular-axis rotation = 0 (n_R=2), no crash; 2 stretches, 2 bends."""
    _, t = _normal_table("CO2")
    diag = sorted(round(t[ax][ax], 3) for ax in ["Rx", "Ry", "Rz"])
    assert diag == [0.0, 1.0, 1.0], f"expected one zero R-external (axis), got {diag}"
    vibs = sorted(round(t[k]["V"], 3) for k in t if k.startswith("Vib"))
    assert vibs == [0.0, 0.0, 1.0, 1.0], f"expected 2 bends (0) + 2 stretches (1), got {vibs}"


def test_benzene_emit_targets():
    """Score-level benzene-EMIT targets from the manuscript. Numbers
    unchanged vs. the pre-2026-07-07 basename (benzene -> C6H6
    was a content-preserving rename; re-verified against a live run)."""
    _, t = _emit_table("C6H6")
    assert abs(t["EMIT 34"]["V"] - 0.6667) < TOL
    assert abs(t["EMIT 35"]["V"] - 0.5774) < TOL
    assert abs(t["EMIT 36"]["V"] - 0.0) < TOL
    assert abs(abs(t["EMIT 2"]["Ry"]) - 0.1427) < TOL
    assert abs(abs(t["EMIT 9"]["Ry"]) - 0.2151) < TOL
    # the s[R] non-monotonicity: EMIT 2 has more Ry character yet smaller |s[Ry]| than EMIT 9
    assert abs(t["EMIT 2"]["Ry"]) < abs(t["EMIT 9"]["Ry"])
    # EMIT 3 sits in the -23.7514 degenerate eigenvalue block: several H atoms carry
    # ~4e-7 A numerical noise (symmetry requires exactly 0) in the raw Gaussian EMIT
    # file. Tscore must reject that noise (EPS_DENOM floor), not promote it to a
    # full-weight unit-vector contribution (regression for the EPS_DISP/EPS_DENOM bug).
    assert abs(t["EMIT 3"]["Tz"] - (-0.1667)) < TOL
    assert abs(t["EMIT 3"]["Tx"] - 0.0) < TOL
    assert abs(t["EMIT 3"]["Ty"] - 0.0) < TOL


def test_tscore_ignores_subthreshold_noise():
    """Tscore must use the same noise floor as Rscore/Vscore (EPS_DENOM), not the looser EPS_DISP.

    Reproduces the mechanism found in degenerate benzene EMIT eigenvectors: a symmetry-required-
    zero atom carries a tiny (~1e-7) numerical noise residual. Before the fix, EPS_DISP=1e-8 let
    this leak through as a full-weight unit-vector contribution to s[T]; it must now be excluded.
    """
    atoms = ["O", "H", "H"]
    coords = [[0.0, 0.0, 0.117], [0.0, 0.757, -0.470], [0.0, -0.757, -0.470]]
    bonds = [(0, 1), (0, 2)]
    scorer = ModeScorer(atoms, coords, bonds)

    # Atom 0 is symmetry-required to be exactly zero in this eigenvector but carries
    # a ~1e-7 noise residual (its *entire* displacement, not a small component riding
    # on top of real motion) -- exactly the pattern seen in the raw benzene EMIT file.
    mode_vec = np.zeros((3, 3))
    mode_vec[0] = [1e-7, 0.0, 0.0]
    mode_vec[1] = [0.0, 0.0, 1.0]
    mode_vec[2] = [0.0, 0.0, 1.0]

    sc = scorer.calculate_scores(mode_vec)
    assert abs(sc["T"]["x"]) < 1e-9, f"noise-level x-displacement leaked into Tx: {sc['T']['x']}"
    assert abs(sc["T"]["z"] - 2.0 / 3.0) < TOL


def test_score_ranges():
    """s[T],s[R] in [-1,1]; s[V_S] in [0,1] across all modes of water and benzene."""
    for mol in ["H2O", "C6H6"]:
        _, t = _normal_table(mol)
        for lbl, row in t.items():
            for k in ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz"]:
                assert -1.0 - TOL <= row[k] <= 1.0 + TOL, f"{mol} {lbl} {k}={row[k]}"
            assert -TOL <= row["V"] <= 1.0 + TOL, f"{mol} {lbl} V={row['V']}"


def test_bond_decomposition_sums_to_vscore():
    """Sum of per-bond |s_AB| equals s[V_S] for every mode (water + benzene).
    s_AB is signed (positive = stretching, negative = compressing); s[V_S]
    itself sums magnitudes, so the comparison must too."""
    for mol in ["H2O", "C6H6"]:
        gp = GaussianParser(os.path.join(ROOT, "data", "logs", f"{mol}.log"))
        data = gp.parse(parse_modes=True)
        scorer = ModeScorer(data["atoms"], data["coords"], data["bonds"])
        for m in scorer.MIT(data["modes"], rotate_modes=True):
            scorer.calculate_scores(m["vector"])
            total = sum(abs(b["s_AB"]) for b in scorer.score_bonds())
            assert abs(total - scorer.Vscore()) < 1e-6, f"{mol}: Sum |s_AB| != s[V_S]"


def test_parser_fails_loud_on_bad_emit():
    """EMITParser raises on a malformed (too-short) matrix rather than scoring garbage."""
    gp = GaussianParser(os.path.join(ROOT, "data", "logs", "H2O.log"))
    natoms = len(gp.parse(parse_modes=False)["atoms"])
    bad = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False)
    bad.write("CART EMIT modes\n0.1 0.2 0.3\nEigenvalues:\n1.0\n")
    bad.close()
    raised = False
    try:
        EMITParser(bad.name, natoms).parse()
    except ValueError:
        raised = True
    finally:
        os.unlink(bad.name)
    assert raised, "EMITParser should raise on a 3-value matrix when 3N x 3N expected"


if __name__ == "__main__":
    # Standalone runner so the suite works even without pytest installed.
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    passed = 0
    for fn in tests:
        try:
            fn()
            print(f"PASS  {fn.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"FAIL  {fn.__name__}: {e}")
        except Exception as e:  # noqa: BLE001
            print(f"ERROR {fn.__name__}: {type(e).__name__}: {e}")
    print(f"\n{passed}/{len(tests)} passed")
    sys.exit(0 if passed == len(tests) else 1)
