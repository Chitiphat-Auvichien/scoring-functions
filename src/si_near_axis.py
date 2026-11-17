"""SI figure for referee point R1.6: stability of s[R_Q] for an atom near the axis, illustrated with water.

The referee notes that eq:rscore normalizes the radius vector r_Q^A, so an atom slightly off the
Q-axis contributes with full weight, and asks how stable s[R_Q] is for such atoms. This module
moves the O atom of H2O off its C2 axis and scores the modes again:

  geometry  O is moved perpendicular to the C2 axis, in direction psi (psi = 0: within the molecular
            plane; psi = 90 deg: out of it). The axis DIRECTIONS are held fixed (the program frame of
            the undistorted molecule); the origin follows the centre of mass, as ModeScorer does for
            every input. O then lies |r_Q^O| = Delta (1 - m_O/M) from the axis; Delta is chosen so that
            |r_Q^O| runs over R_MIN..R_MAX on a log grid. The figure thus shows how eq:rscore responds
            to an atom's distance from a given axis: the case of an atom that lies off its axis by a
            coordinate error (R1.6 side report: rounding residues).
  modes     the three normal modes of the undistorted molecule (delta, nu_s, nu_as), unchanged, and the
            ideal rotation about the C2 axis rebuilt for each moved geometry as the program builds its
            references (ModeScorer.construct_R: d^A = Q x r^A from the CM, normalized).
  scores    s[R_C2] with the program's cut-offs (|r_Q^A| <= ON_AXIS_TOL = 1e-5 A: on the axis, not
            counted; |d^A| <= EPS_DENOM = 1e-6: term 0), and without any tolerance (only exact zeros),
            plus O's own contribution w^O c^O / n to the score.

Each term c^A = (r_Q^A x d^A).Q / (|r_Q^A| |d^A|) depends only on directions, so O's term does not
shrink as O approaches the axis: s[R_C2] of nu_as steps from 0 to a value set by psi when O crosses
the cut-off, and then stays constant. In delta and nu_s, O moves along the axis, so its term is 0 in
every direction psi. In the ideal rotation, O moves round the axis like every other atom, so its term
is +1 and s[R_C2] stays 1; it would be (N - 1)/N = 2/3 if O were counted while on the axis (no N_Q).

Why the axis directions are held fixed: if the principal axes are recomputed after O is moved, water
hides the effect. Three atoms always lie in a plane, so O moved out of the plane leaves a C2v
molecule whose recomputed C2 axis passes through O again (recomputed_axis_distance()), and an
in-plane offset is parallel to O's nu_as displacement (term 0).

Frame: the program frame of main.build_scorer_and_final (src/frame.py, R1.5 rule): for water the
principal axes x = a, y = b, z = c, right-handed, origin at the CM; y is the C2 axis, x lies in the
molecular plane, z is perpendicular to it.

Reference values marked on the figure (R1.6 side report; DECISIONS.md, 2026-11-09): the on-axis
atoms of the 85 roster molecules lie at most ROUNDING_MAX from their axis for any rounding of the
six printed decimals; the closest atom off an axis by symmetry lies CLOSEST_OFF from it (OBr4 O).

    python -m src.si_near_axis              # writes the CSV and the figure

Outputs: data/results/si_near_axis_water.csv (one row per (psi, |r_Q^O|) grid point; columns
<quantity>_<mode>) and data/figures/fig_si_near_axis_water.{pdf,png}.
"""

import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, FancyArrowPatch

from src.figures import COLORS, _savefig, _style
from src.scoring import EPS_DENOM, ON_AXIS_TOL, ModeScorer

FIG_LABEL = "fig_si_near_axis_water"
CSV_NAME = "si_near_axis_water.csv"
FRAME = ("program frame (R1.5 rule) of the undistorted molecule, axis directions held fixed: x=a (in plane), "
         "y=b (C2), z=c (out of plane), right-handed; origin at the CM of each moved geometry")
FRAME_SHORT = "program frame (R1.5), axis directions fixed, origin at CM"   # per-row CSV label (FRAME in full above)
MODES = ["delta", "nu_s", "nu_as", "R_C2"]
R_MIN, R_MAX, N_R = 1e-9, 1e-2, 141          # O's distance from the C2 axis (A), log grid (O moved <= 0.09 A)
N_PSI = 145                                  # psi = 0, 2.5, ..., 360 deg
ROUNDING_MAX, CLOSEST_OFF = 3.5e-6, 1.33e-5  # A, R1.6 side report
PSI_CUT = 90.0                               # direction of the line cuts (O moved out of the plane)
_TEX = {"delta": r"$\delta$", "nu_s": r"$\nu_\mathrm{s}$", "nu_as": r"$\nu_\mathrm{as}$",
        "R_C2": r"ideal rotation R$_{C_2}$"}
C_O = "#332288"                              # the near-axis atom (indigo, as in the R1.6 side report)


# ---------------------------------------------------------------------------------------- geometry
def water(data_dir="data"):
    """Undistorted water in the program frame: coordinates, symbols, bonds, masses, the three
    vibrations (sorted by frequency), and the indices of O, the C2 axis, the in-plane axis and the
    out-of-plane axis."""
    from main import load_inputs, build_scorer_and_final

    raw, _ = load_inputs("H2O", "normal", data_dir=data_dir)
    scorer, final = build_scorer_and_final(raw, "normal")
    X0 = np.asarray(scorer.coords, float)
    sym = list(raw["atoms"])
    vib = sorted([m for m in final if str(m.get("label", "")).startswith("Vib")], key=lambda m: m["frequency"])
    o = [i for i, s in enumerate(sym) if s.upper() == "O"][0]
    r_axis = [np.linalg.norm(np.delete(X0[o], k)) for k in range(3)]
    k_c2 = int(np.argmin(r_axis))
    k_out = int(np.argmin(np.abs(X0).max(axis=0)))           # every atom has 0 on the out-of-plane axis
    k_in = ({0, 1, 2} - {k_c2, k_out}).pop()
    return dict(raw=raw, X0=X0, sym=sym, bonds=raw["bonds"], masses=np.asarray(scorer.masses, float),
                vib={n: np.asarray(m["vector"], float) for n, m in zip(MODES[:3], vib)},
                freq={n: float(m["frequency"]) for n, m in zip(MODES[:3], vib)},
                o=o, k_c2=k_c2, k_in=k_in, k_out=k_out, r_O_axis=float(r_axis[k_c2]),
                out_of_plane_max=float(np.abs(X0[:, k_out]).max()))


def direction(w, psi_deg):
    u = np.zeros(3)
    u[w["k_in"]] = np.cos(np.radians(psi_deg)); u[w["k_out"]] = np.sin(np.radians(psi_deg))
    return u


def moved(w, r, psi_deg):
    """Coordinates, at their CM, with O moved so that it lies r from the C2 axis in direction psi
    (axis directions fixed)."""
    m = w["masses"]; o = w["o"]
    X = w["X0"].copy()
    X[o] += r / (1.0 - m[o] / m.sum()) * direction(w, psi_deg)
    return X - (m[:, None] * X).sum(axis=0) / m.sum()


def ideal_rotation(w, X):
    """The program's ideal rotation about the C2 axis (ModeScorer.construct_R) for coordinates X at their CM."""
    sc = ModeScorer(w["sym"], X, w["bonds"])
    return np.asarray(sc.construct_R()[w["k_c2"]]["vector"], float)


# ---------------------------------------------------------------------------------------- scores
def rscore_pieces(X, d, k, r_tol, d_tol):
    """eq:rscore about axis k with cut-offs r_tol (on the axis) and d_tol (displacement).
    Returns (score, n, w, c): counting weights w^A and terms c^A (0 where not computed)."""
    Q = np.eye(3)[k]
    rp = X - np.outer(X @ Q, Q)
    lr = np.linalg.norm(rp, axis=1); ld = np.linalg.norm(d, axis=1)
    w = lr > r_tol; ok = w & (ld > d_tol)
    c = np.zeros(len(X))
    c[ok] = (np.cross(rp[ok], d[ok]) @ Q) / (lr[ok] * ld[ok])
    n = int(w.sum())
    return (c[w].sum() / n if n else 0.0), n, w, c


def score_point(w, X):
    """All quantities of one moved geometry."""
    k, o = w["k_c2"], w["o"]
    vecs = dict(w["vib"]); vecs["R_C2"] = ideal_rotation(w, X)
    row = {}
    for name, d in vecs.items():
        s, n, wt, c = rscore_pieces(X, d, k, ON_AXIS_TOL, EPS_DENOM)
        s0, n0, _, _ = rscore_pieces(X, d, k, 0.0, 0.0)
        row.update({f"s_{name}": s, f"n_{name}": n, f"cO_{name}": c[o] * wt[o],
                    f"contribO_{name}": c[o] * wt[o] / n if n else 0.0,
                    f"s_tol0_{name}": s0, f"dO_{name}": float(np.linalg.norm(d[o]))})
    return row, vecs


def scan(w, n_r=N_R, n_psi=N_PSI, validate_every=29):
    """The (psi, |r_Q^O|) grid. Every validate_every-th point is also scored with the program's own
    ModeScorer.calculate_scores (same fixed coordinates); the largest difference is returned."""
    rows, dev = [], 0.0
    radii = np.logspace(np.log10(R_MIN), np.log10(R_MAX), n_r)
    psis = np.linspace(0.0, 360.0, n_psi)
    i = 0
    for psi in psis:
        for r in radii:
            X = moved(w, r, psi)
            row, vecs = score_point(w, X)
            Q = np.eye(3)[w["k_c2"]]
            row.update(psi_deg=psi, r_target=r, r_Q_O=float(np.linalg.norm(X[w["o"]] - (X[w["o"]] @ Q) * Q)))
            if i % validate_every == 0:
                sc = ModeScorer(w["sym"], X, w["bonds"])
                for name, d in vecs.items():
                    dev = max(dev, abs(sc.calculate_scores(d)["R"]["xyz"[w["k_c2"]]] - row[f"s_{name}"]))
            rows.append(row); i += 1
    df = pd.DataFrame(rows)
    df["frame"] = FRAME_SHORT
    first = ["psi_deg", "r_target", "r_Q_O"]
    return df[first + [c for c in df.columns if c not in first]], dev


def no_nq_rotation(w):
    """s[R_C2] of the ideal rotation with O exactly on the axis but counted (no N_Q): O does not
    move, its term is 0, and the score is (N - 1)/N."""
    d = ideal_rotation(w, w["X0"])
    _, _, _, c = rscore_pieces(w["X0"], d, w["k_c2"], -1.0, EPS_DENOM)     # r_tol < 0: every atom counted
    return float(c.sum() / len(c))


def recomputed_axis_distance(w, r=1e-3):
    """O's distance from the recomputed C2 axis (principal axis nearest the old one, CM recomputed)
    after O is moved r out of the molecular plane, and after it is moved r within the plane."""
    m = w["masses"]; out = {}
    for name, psi in (("out_of_plane", 90.0), ("in_plane", 0.0)):
        X = moved(w, r, psi)
        Xc = X - (m[:, None] * X).sum(axis=0) / m.sum()
        I = np.eye(3) * np.sum(m * np.sum(Xc ** 2, axis=1)) - (m[:, None] * Xc).T @ Xc
        _, W = np.linalg.eigh(I)
        q = W[:, int(np.argmax(np.abs(W[w["k_c2"], :])))]
        p = Xc[w["o"]]
        out[name] = float(np.linalg.norm(p - (p @ q) * q))
    return out


# ---------------------------------------------------------------------------------------- figure
def _schematic(ax, w):
    """Water viewed down the C2 axis (offset exaggerated): O moved off the axis in direction psi."""
    H = [i for i in range(len(w["sym"])) if i != w["o"]]
    Hx = float(abs(w["X0"][H[0], w["k_in"]]))
    ax.set_aspect("equal"); ax.axis("off")
    ax.plot([-1.25, 1.25], [0, 0], color="0.75", lw=0.8, zorder=0)
    for x in (-Hx, Hx):
        ax.add_patch(plt.Circle((x, 0), 0.09, color="0.65", zorder=2))
        ax.text(x, -0.16, "H", ha="center", va="top", fontsize=7)
    ax.plot(0, 0, marker="o", ms=7, mfc="white", mec="black", mew=0.7, zorder=3)
    ax.plot(0, 0, marker=".", ms=2.5, color="black", zorder=4)
    ax.text(0, -0.16, f"$C_2$ axis ({'xyz'[w['k_c2']]})", fontsize=6.5, ha="center", va="top")
    psi, rr = np.radians(60.0), 0.62
    O = np.array([rr * np.cos(psi), rr * np.sin(psi)])
    ax.add_patch(FancyArrowPatch((0, 0), O, arrowstyle="-|>", mutation_scale=6, lw=0.8, color=C_O,
                                 shrinkA=3, shrinkB=6))
    ax.add_patch(plt.Circle(O, 0.1, color=C_O, zorder=3))
    ax.text(O[0], O[1] + 0.14, "O", fontsize=7, color=C_O, ha="center", va="bottom")
    ax.text(0.05, 0.36, r"$\boldsymbol{r}_Q^{\mathrm{O}}$", fontsize=7, color=C_O, ha="right")
    ax.add_patch(Arc((0, 0), 0.5, 0.5, theta1=0, theta2=np.degrees(psi), lw=0.6))
    ax.text(0.3, 0.1, r"$\psi$", fontsize=7)
    ax.add_patch(FancyArrowPatch(O, O + np.array([0.45, 0]), arrowstyle="-|>", mutation_scale=7, lw=1.1,
                                 color=COLORS["stretching"], shrinkA=5, shrinkB=0))
    ax.text(O[0] + 0.52, O[1], r"$\boldsymbol{d}^{\mathrm{O}}$", fontsize=7, color=COLORS["stretching"], va="center")
    ax.text(O[0] + 0.3, O[1] - 0.14, r"($\nu_\mathrm{as}$)", fontsize=6, color=COLORS["stretching"], ha="center", va="top")
    x0, y0 = -1.15, 0.62
    ax.annotate("", xy=(x0, y0 + 0.32), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", lw=0.6, color="0.3"))
    ax.annotate("", xy=(x0 + 0.32, y0), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", lw=0.6, color="0.3"))
    ax.text(x0, y0 + 0.35, "xyz"[w["k_out"]], fontsize=6.5, ha="center", va="bottom", color="0.3")
    ax.text(x0 + 0.35, y0, "xyz"[w["k_in"]], fontsize=6.5, ha="left", va="center", color="0.3")
    ax.text(0.15, -0.45, "viewed down the $C_2$ axis;\ngrey line: molecular plane;\noffset exaggerated", fontsize=6,
            color="0.45", ha="center", va="top")
    ax.set_xlim(-1.3, 1.6); ax.set_ylim(-0.95, 1.05)
    ax.set_title("(a) geometry", loc="left", fontsize=8)


def _guides(ax, horizontal):
    """Shaded: where the on-axis atoms of the roster lie (any rounding); dashed: the cut-off;
    dotted: the closest off-axis atom of the roster."""
    span = ax.axhspan if horizontal else ax.axvspan
    line = ax.axhline if horizontal else ax.axvline
    span(R_MIN / 10, ROUNDING_MAX, color="0.85", alpha=0.6, lw=0, zorder=0)
    line(ON_AXIS_TOL, color=COLORS["threshold"], ls="--", lw=0.8, zorder=3)
    line(CLOSEST_OFF, color=COLORS["threshold"], ls=":", lw=0.8, zorder=3)


def plot(df, w, s_no_nq, out_dir, label=FIG_LABEL):
    _style()
    plt.rcParams.update({"font.size": 7.5, "axes.labelsize": 7.5, "xtick.labelsize": 7, "ytick.labelsize": 7,
                         "legend.fontsize": 6.5, "axes.titlesize": 8})
    fig = plt.figure(figsize=(7.2, 2.55))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.15, 1.15, 0.05, 1.1], wspace=0.3)
    _schematic(fig.add_subplot(gs[0]), w)

    ax = fig.add_subplot(gs[1])
    piv = df.pivot(index="r_target", columns="psi_deg", values="contribO_nu_as").sort_index().sort_index(axis=1)
    third = 1.0 / 3.0
    im = ax.pcolormesh(piv.columns.to_numpy(), piv.index.to_numpy(), piv.to_numpy(), cmap="PuOr",
                       vmin=-third, vmax=third, shading="nearest", rasterized=True)
    ax.set_yscale("log"); ax.set_ylim(R_MIN, R_MAX); ax.set_xlim(0, 360)
    ax.axhline(ON_AXIS_TOL, color=COLORS["threshold"], ls="--", lw=0.8)
    ax.set_xticks([0, 90, 180, 270, 360]); ax.set_yticks([1e-9, 1e-7, 1e-5, 1e-3])
    ax.set_xlabel(r"direction of the offset, $\psi$ (deg)")
    ax.set_ylabel(r"$|\boldsymbol{r}_Q^{\mathrm{O}}|$ (Å)", labelpad=1)
    ax.set_title(r"(b) O's share of $s[\mathrm{R}_{C_2}]$ in $\nu_\mathrm{as}$", loc="left")
    ax.text(180, 1e-7, "O counted as on the axis:\nnot in $n$, share 0", fontsize=6.5, ha="center", va="center",
            color="0.3")
    cax = fig.add_subplot(gs[2])
    cb = fig.colorbar(im, cax=cax, ticks=[-third, 0, third])
    cb.ax.set_yticklabels([r"$-1/3$", "0", r"$+1/3$"], fontsize=7)
    cb.ax.set_title(r"$w^{\mathrm{O}}c^{\mathrm{O}}/n$", fontsize=7, pad=4); cb.outline.set_linewidth(0.6)

    ax = fig.add_subplot(gs[3])
    pos = ax.get_position(); ax.set_position([pos.x0 + 0.07, pos.y0, pos.width - 0.07, pos.height])
    cut = df[np.isclose(df.psi_deg, PSI_CUT)].sort_values("r_Q_O")
    style = {"R_C2": dict(color=COLORS["external"], ls="-", lw=1.6),
             "delta": dict(color=COLORS["bending"], ls="-", lw=1.0),
             "nu_s": dict(color=COLORS["stretching"], ls=(0, (1.5, 1.5)), lw=1.3),
             "nu_as": dict(color=COLORS["stretching"], ls="-", lw=1.3)}
    for name in ("R_C2", "delta", "nu_s", "nu_as"):
        lab = _TEX[name] + (f" ({w['freq'][name]:.0f} cm$^{{-1}}$)" if name in w["freq"] else "")
        ax.plot(cut.r_Q_O, cut[f"s_{name}"], label=lab, **style[name])
    ax.plot(cut.r_Q_O, cut["s_tol0_nu_as"], color=COLORS["stretching"], ls="--", lw=0.8, alpha=0.6,
            label=r"$\nu_\mathrm{as}$, no tolerance")
    ax.plot([R_MIN], [s_no_nq], marker="o", mfc="white", mec=COLORS["external"], mew=1.1, ms=5, ls="none",
            label=r"R$_{C_2}$, O on the axis but counted (no $N_Q$)", clip_on=False, zorder=5)
    _guides(ax, horizontal=False)
    ax.set_xscale("log"); ax.set_xlim(R_MIN, R_MAX); ax.set_ylim(-0.45, 1.1)
    ax.set_xticks([1e-9, 1e-7, 1e-5, 1e-3])
    ax.set_xlabel(r"$|\boldsymbol{r}_Q^{\mathrm{O}}|$ (Å), at $\psi=90°$")
    ax.set_ylabel(r"$s[\mathrm{R}_{C_2}]$", labelpad=1)
    ax.set_title("(c) scores", loc="left")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    h, l = ax.get_legend_handles_labels()
    h += [Line2D([], [], color=COLORS["threshold"], ls="--", lw=0.8), Line2D([], [], color=COLORS["threshold"], ls=":", lw=0.8),
          Patch(color="0.85", alpha=0.6, lw=0)]
    l += [r"on-axis cut-off, $10^{-5}$ Å", r"closest off-axis atom of the roster, $1.3\times10^{-5}$ Å",
          r"on-axis atoms of the roster, any rounding ($\leq3.5\times10^{-6}$ Å)"]
    fig.legend(h, l, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0),
               handlelength=2.0, columnspacing=1.0)
    paths = _savefig(fig, out_dir, label)
    plt.close(fig)
    return paths


# ---------------------------------------------------------------------------------------- run
def summary(df, w, s_no_nq, dev):
    """Numbers quoted in the caption/text, plus the checks (all should hold; see tests)."""
    above = df[df.r_Q_O > ON_AXIS_TOL]; below = df[df.r_Q_O <= ON_AXIS_TOL]
    sinp = np.sin(np.radians(df.psi_deg))
    rec = recomputed_axis_distance(w)
    return dict(
        frame=FRAME, c2_axis="xyz"[w["k_c2"]], in_plane_axis="xyz"[w["k_in"]], out_of_plane_axis="xyz"[w["k_out"]],
        frequencies=w["freq"], O_distance_undistorted=w["r_O_axis"], out_of_plane_max=w["out_of_plane_max"],
        max_dev_vs_ModeScorer=float(dev),
        nu_as_below_max_abs=float(below.s_nu_as.abs().max()),
        nu_as_above_sign_x_3_vs_sinpsi=float(np.abs(3 * above.s_nu_as.abs() - np.abs(sinp[above.index])).max()),
        nu_as_above_vs_sinpsi_upto_100fA=float(np.abs(3 * above[above.r_Q_O <= 1e-4].s_nu_as.abs()
                                                     - np.abs(sinp[above[above.r_Q_O <= 1e-4].index])).max()),
        nu_as_above_max_abs=float(above.s_nu_as.abs().max()),
        nu_as_tol0_max_dev_from_step=float(np.abs(3 * df.s_tol0_nu_as.abs() - np.abs(sinp)).max()),
        delta_nu_s_max_abs=float(df[["s_delta", "s_nu_s"]].abs().to_numpy().max()),
        R_C2_min=float(df.s_R_C2.min()), R_C2_max=float(df.s_R_C2.max()),
        R_C2_O_smallest_disp_when_counted=float(above.dO_R_C2.min()),
        R_C2_no_NQ=s_no_nq,
        recomputed_axis_O_distance_out_of_plane_1mA=rec["out_of_plane"],
        recomputed_axis_O_distance_in_plane_1mA=rec["in_plane"])


def run(data_dir="data", out_dir="data/figures", results_dir="data/results", n_r=N_R, n_psi=N_PSI):
    """Scan, write the wide CSV and the figure. Returns dict(csv, pdf, png, frame, summary)."""
    w = water(data_dir)
    df, dev = scan(w, n_r=n_r, n_psi=n_psi)
    s_no_nq = no_nq_rotation(w)
    os.makedirs(results_dir, exist_ok=True)
    csv_path = os.path.join(results_dir, CSV_NAME)
    df.to_csv(csv_path, index=False, float_format="%.8g")
    pdf, png = plot(df, w, s_no_nq, out_dir)
    return dict(csv=csv_path, pdf=pdf, png=png, frame=FRAME, summary=summary(df, w, s_no_nq, dev))


if __name__ == "__main__":
    res = run()
    for k, v in res["summary"].items():
        print(f"{k:48s} {v}")
    print("wrote", res["csv"], res["pdf"], res["png"])
