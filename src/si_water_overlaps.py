"""SI figure for referee point R1.8: water T/R scores next to Cartesian and mass-weighted overlaps.

The referee notes that s[T_Q] and s[R_Q] (eq:tscore, eq:rscore) give nonzero cross-scores for modes
that are orthogonal (water: -0.33 between T_z and R_x; up to 0.30 for nu_as). The T/R definitions are
kept; the SI shows instead that "orthogonal" depends on the metric. For each of the nine H2O modes
(the six ideal T/R references, then delta, nu_s, nu_as) this module computes, against the six ideal
references T_x..R_z:

  MS      the manuscript's scores s[T_Q], s[R_Q], from the program's own ModeScorer.
  o_Cart  plain Cartesian overlap (cosine), metric w^A = 1.
  o_MW    mass-weighted overlap, metric w^A = m^A.

Overlap in metric w (R1.8 side report, eq. general and its m = 2 identity):
    <a, b>_w = sum_A w^A a^A . b^A,     o_w(d, e) = <d, e>_w / (||d||_w ||e||_w),
with ideal vectors e^A = Q (T_Q) and e^A = Q x r^A (R_Q), r^A from the centre of mass. The side
report's general score s_{w,m} with moment m = 2 equals o_w exactly, so the o_MW panel is also the
"m = 2, w^A = m^A" member of that family. Normal modes are orthonormal only in the mass-weighted
metric: the ideal T/R block of o_MW is the identity and the vibrations have o_MW = 0 (Eckart
conditions), whereas in the Cartesian metric T_z and R_x overlap by -0.76 (T is centred at the
centroid, R at the CM).

Frame: the program frame of main.build_scorer_and_final (src/frame.py, R1.5 rule): for an
asymmetric top such as water, the principal axes of inertia in increasing-moment order
(x = a, y = b, z = c), right-handed, origin at the CM. Overlaps (and, up to the sign set by the
axis naming, the scores) do not depend on this choice beyond the axis labels.

    python -m src.si_water_overlaps          # writes the CSV and the figure

Outputs: data/results/si_water_overlaps.csv (one row per mode; columns <quantity>_<slot>) and
data/figures/fig_si_water_overlaps.{pdf,png}.
"""

import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from src.figures import _style, _savefig

SLOTS = ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz"]
QUANTITIES = ["MS", "oCart", "oMW"]
FRAME = "program frame (R1.5 rule): principal axes x=a, y=b, z=c, right-handed, origin at CM"
FIG_LABEL = "fig_si_water_overlaps"
CSV_NAME = "si_water_overlaps.csv"
# H2O vibrations in increasing frequency: bend, symmetric stretch, asymmetric stretch.
WATER_VIB_NAMES = ["delta", "nu_s", "nu_as"]
_TEX = {"Tx": r"T$_\mathrm{x}$", "Ty": r"T$_\mathrm{y}$", "Tz": r"T$_\mathrm{z}$",
        "Rx": r"R$_\mathrm{x}$", "Ry": r"R$_\mathrm{y}$", "Rz": r"R$_\mathrm{z}$",
        "delta": r"$\delta$", "nu_s": r"$\nu_\mathrm{s}$", "nu_as": r"$\nu_\mathrm{as}$"}
_TITLE = {"MS": r"(a) $s[\mathrm{T}_Q]$, $s[\mathrm{R}_Q]$", "oCart": r"(b) $o_\mathrm{Cart}$",
          "oMW": r"(c) $o_\mathrm{MW}$"}
_NORM_TOL = 1e-12


def weighted_overlap(d, e, w):
    """Cosine between displacement fields d and e (both (N, 3)) in the metric with atom weights w:
    <d, e>_w / (||d||_w ||e||_w). Returns 0 if either vector has zero norm (unit(0) := 0)."""
    d = np.asarray(d, float)
    e = np.asarray(e, float)
    w = np.asarray(w, float)
    de = np.sum(w * np.einsum("ij,ij->i", d, e))
    nd = np.sqrt(np.sum(w * np.einsum("ij,ij->i", d, d)))
    ne = np.sqrt(np.sum(w * np.einsum("ij,ij->i", e, e)))
    if nd < _NORM_TOL or ne < _NORM_TOL:
        return 0.0
    return float(de / (nd * ne))


def ideal_references(scorer):
    """The six ideal T/R displacement fields (T_x..R_z) in the scorer's current frame."""
    refs = scorer.construct_T() + scorer.construct_R()
    return [np.asarray(r["vector"], float) for r in refs]


def ms_scores(scorer, d):
    """The manuscript's s[T_x..T_z], s[R_x..R_z] for one displacement field."""
    s = scorer.calculate_scores(np.asarray(d, float))
    return np.array([s["T"][q] for q in "xyz"] + [s["R"][q] for q in "xyz"])


def compute(scorer, modes):
    """MS, o_Cart and o_MW (each len(modes) x 6) for a list of mode dicts."""
    refs = ideal_references(scorer)
    m = np.asarray(scorer.masses, float)
    one = np.ones_like(m)
    out = {"MS": [], "oCart": [], "oMW": []}
    for md in modes:
        d = np.asarray(md["vector"], float)
        out["MS"].append(ms_scores(scorer, d))
        out["oCart"].append([weighted_overlap(d, e, one) for e in refs])
        out["oMW"].append([weighted_overlap(d, e, m) for e in refs])
    return {k: np.array(v, float) for k, v in out.items()}


def mode_names(modes, molecule):
    """Display keys: ideal slots keep their labels; H2O vibrations get delta/nu_s/nu_as."""
    names = []
    vib = [i for i, md in enumerate(modes) if md.get("label") not in SLOTS]
    order = sorted(vib, key=lambda i: modes[i]["frequency"])
    rank = {i: k for k, i in enumerate(order)}
    for i, md in enumerate(modes):
        if md.get("label") in SLOTS:
            names.append(md["label"])
        elif molecule == "H2O" and len(order) == len(WATER_VIB_NAMES):
            names.append(WATER_VIB_NAMES[rank[i]])
        else:
            names.append(md.get("label", f"Vib {i + 1}"))
    return names


def build_table(scorer, modes, molecule):
    """One wide DataFrame: rows = modes, columns = <quantity>_<slot> (+ mode info and frame)."""
    vals = compute(scorer, modes)
    names = mode_names(modes, molecule)
    rows = []
    for i, md in enumerate(modes):
        row = {"molecule": molecule, "mode": names[i], "label": md.get("label", ""),
               "frequency_cm-1": float(md.get("frequency", 0.0))}
        for q in QUANTITIES:
            for k, s in enumerate(SLOTS):
                row[f"{q}_{s}"] = vals[q][i, k]
        row["frame"] = FRAME
        rows.append(row)
    return pd.DataFrame(rows)


def _cell_text(v):
    """'.04' / '−.33' / '1.00' (two decimals, no leading zero, true minus)."""
    s = f"{v:.2f}"
    if s.startswith("0."):
        s = s[1:]
    elif s.startswith("-0."):
        s = "−" + s[2:]
    elif s.startswith("-"):
        s = "−" + s[1:]
    return s


def plot(df, out_dir, label=FIG_LABEL, blank_below=0.005):
    """Three 9 x 6 heatmaps (MS, o_Cart, o_MW) with one shared diverging colour bar centred at 0."""
    _style()
    rows = list(df["mode"])
    n_ext = sum(r in SLOTS for r in rows)
    cmap = plt.get_cmap("PuOr")  # positive = purple, as on the R1.8 slide
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 3.5), sharey=True,
                            gridspec_kw=dict(wspace=0.07))
    for ax, q in zip(axs, QUANTITIES):
        S = df[[f"{q}_{s}" for s in SLOTS]].to_numpy(float)
        im = ax.imshow(S, cmap=cmap, vmin=-1.0, vmax=1.0, aspect="auto")
        for i in range(S.shape[0]):
            for k in range(S.shape[1]):
                v = S[i, k]
                if abs(v) < blank_below:
                    continue
                ax.text(k, i, _cell_text(v), ha="center", va="center", fontsize=7.5,
                        color="white" if abs(v) > 0.55 else "black")
        if 0 < n_ext < len(rows):
            ax.axhline(n_ext - 0.5, color="0.15", lw=0.9)
        ax.set_xticks(range(len(SLOTS)))
        ax.set_xticklabels([_TEX[s] for s in SLOTS], fontsize=9)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([_TEX.get(r, r) for r in rows], fontsize=9)
        ax.tick_params(length=2.5)
        ax.set_title(_TITLE[q], fontsize=11, loc="left")
        if ax is not axs[0]:
            ax.tick_params(axis="y", left=False)
        for sp in ax.spines.values():
            sp.set_visible(True)
            sp.set_linewidth(0.6)
    cb = fig.colorbar(im, ax=list(axs), fraction=0.025, pad=0.02, ticks=[-1, -0.5, 0, 0.5, 1])
    cb.ax.tick_params(labelsize=8, length=2.5)
    cb.set_label("score / overlap", fontsize=9)
    cb.outline.set_linewidth(0.6)
    paths = _savefig(fig, out_dir, label)
    plt.close(fig)
    return paths


def run(data_dir="data", out_dir="data/figures", results_dir="data/results", molecule="H2O"):
    """Score and overlap the molecule's ideal T/R references and normal modes; write the wide CSV
    and the figure. Returns dict(table=DataFrame, csv=path, pdf=path, png=path, frame=str)."""
    from main import load_inputs, build_scorer_and_final

    raw, _ = load_inputs(molecule, "normal", data_dir=data_dir)
    scorer, final = build_scorer_and_final(raw, "normal")
    df = build_table(scorer, final, molecule)
    os.makedirs(results_dir, exist_ok=True)
    csv_path = os.path.join(results_dir, CSV_NAME)
    df.to_csv(csv_path, index=False, float_format="%.6f")
    pdf, png = plot(df, out_dir)
    return dict(table=df, csv=csv_path, pdf=pdf, png=png, frame=FRAME,
                frame_info=getattr(scorer, "frame_info", None))


if __name__ == "__main__":
    res = run()
    pd.set_option("display.width", 200)
    print(res["table"].drop(columns=["frame", "molecule"]).round(3).to_string(index=False))
    print("frame:", res["frame"], res["frame_info"])
    print("wrote", res["csv"], res["pdf"], res["png"])
