"""Golden-reference regression tests for the scoring engine.

Run from ``Github/scoring-functions/``:
    py -m pytest tests/            (with pytest)
    py tests/test_scores.py       (standalone; no pytest needed)

Golden values: the published water table (tab:water, consensus s[R]) and the
benzene-EMIT score-level targets. Invariants: Sum s_AB == s[V_S] and score ranges.
These freeze the engine's numerical behaviour so regressions are caught immediately.
"""
import os
import sys
import tempfile

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
    """Reproduce tab:water to 3 dp (consensus s[R])."""
    _, t = _normal_table("water")
    for ext in ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz"]:
        assert abs(t[ext][ext] - 1.0) < TOL, f"external {ext} should be 1.000"
    assert abs(t["Tx"]["Rz"] - 0.049) < TOL
    assert abs(t["Tz"]["Rx"] - (-0.333)) < TOL
    assert abs(t["Vib1"]["V"] - 0.061) < TOL          # bend (sigma)
    assert abs(t["Vib2"]["V"] - 1.000) < TOL          # nu_s
    assert abs(t["Vib3"]["V"] - 0.998) < TOL          # nu_as
    assert abs(t["Vib3"]["Tx"] - (-0.190)) < TOL
    assert abs(t["Vib3"]["Rz"] - (-0.295)) < TOL


def test_co2_linear():
    """Linear molecule: molecular-axis rotation = 0 (n_R=2), no crash; 2 stretches, 2 bends."""
    _, t = _normal_table("co2_mp2_3-21g")
    diag = sorted(round(t[ax][ax], 3) for ax in ["Rx", "Ry", "Rz"])
    assert diag == [0.0, 1.0, 1.0], f"expected one zero R-external (axis), got {diag}"
    vibs = sorted(round(t[k]["V"], 3) for k in t if k.startswith("Vib"))
    assert vibs == [0.0, 0.0, 1.0, 1.0], f"expected 2 bends (0) + 2 stretches (1), got {vibs}"


def test_benzene_emit_targets():
    """Score-level benzene-EMIT targets from the manuscript."""
    _, t = _emit_table("benzene")
    assert abs(t["EMIT 34"]["V"] - 0.6667) < TOL
    assert abs(t["EMIT 35"]["V"] - 0.5774) < TOL
    assert abs(t["EMIT 36"]["V"] - 0.0) < TOL
    assert abs(abs(t["EMIT 2"]["Ry"]) - 0.1427) < TOL
    assert abs(abs(t["EMIT 9"]["Ry"]) - 0.2151) < TOL
    # the s[R] non-monotonicity: EMIT 2 has more Ry character yet smaller |s[Ry]| than EMIT 9
    assert abs(t["EMIT 2"]["Ry"]) < abs(t["EMIT 9"]["Ry"])


def test_score_ranges():
    """s[T],s[R] in [-1,1]; s[V_S] in [0,1] across all modes of water and benzene."""
    for mol in ["water", "benzene"]:
        _, t = _normal_table(mol)
        for lbl, row in t.items():
            for k in ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz"]:
                assert -1.0 - TOL <= row[k] <= 1.0 + TOL, f"{mol} {lbl} {k}={row[k]}"
            assert -TOL <= row["V"] <= 1.0 + TOL, f"{mol} {lbl} V={row['V']}"


def test_bond_decomposition_sums_to_vscore():
    """Sum of per-bond s_AB equals s[V_S] for every mode (water + benzene)."""
    for mol in ["water", "benzene"]:
        gp = GaussianParser(os.path.join(ROOT, "data", "logs", f"{mol}.log"))
        data = gp.parse(parse_modes=True)
        scorer = ModeScorer(data["atoms"], data["coords"], data["bonds"])
        for m in scorer.MIT(data["modes"], rotate_modes=True):
            scorer.calculate_scores(m["vector"])
            total = sum(b["s_AB"] for b in scorer.score_bonds())
            assert abs(total - scorer.Vscore()) < 1e-6, f"{mol}: Sum s_AB != s[V_S]"


def test_parser_fails_loud_on_bad_emit():
    """EMITParser raises on a malformed (too-short) matrix rather than scoring garbage."""
    gp = GaussianParser(os.path.join(ROOT, "data", "logs", "water.log"))
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
