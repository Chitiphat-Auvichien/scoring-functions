"""SI figure for referee point R1.6: an atom close to the rotation axis in s[R_Q], illustrated with NH3.

The referee notes that eq:rscore normalizes the radius vector r_Q^A, so an atom slightly off the Q-axis
contributes with full weight. This module moves the N atom of NH3 over the plane perpendicular to its
C3 axis and scores mode 3 (e bend, 1796 cm-1) again, with and without the on-axis tolerance:

  geometry  N is first put exactly on the C3 axis (over the centroid of the three H atoms; the printed
            geometry has it ~7e-8 A off, a rounding residue). It is then placed at distance rho and
            direction theta from the axis, the H atoms fixed. The axis DIRECTIONS are those of the
            undistorted molecule (program frame, src/frame.py); the origin follows the centre of mass, as
            ModeScorer does for any input, so N itself is moved rho / (1 - m_N/M).
  grid      log-polar: rho on a log grid, theta every 2.5 deg, plus N exactly on the axis (rho = 0).
            With the tolerance: rho = 1e-9..1e-2 A; without it: 1e-20..1e-2 A (past the floating-point
            limit, ~1e-17 A, at which the centre of mass no longer resolves N's shift).
  scores    the contribution of N, i.e. its summand (r_Q/|r_Q| x d/|d|) . Q_hat in eq:rscore before division
            by N - N_Q (0 where N is counted as on the axis), and s[R_z] of mode 3, both
            with the program's tolerance (|r_Q^A| <= ON_AXIS_TOL = 1e-5 A: on the axis, in N_Q) and
            without any (N_Q = 0; a summand is 0 only where it is 0/0).

Because the summand depends only on the direction of r_Q^N, N's contribution keeps its full size, with a
sign set by the direction of its offset, at every distance; within the tolerance N is treated as exactly
on the axis and s[R_z] does not change (0.171, the H atoms alone).

    python -m src.si_near_axis              # writes the CSV and the figure

Outputs: data/results/si_near_axis_nh3.csv (one row per map point; column `tolerance` = "1e-5 A" or
"none"; `frame` = "program": program frame of the undistorted molecule, axis directions fixed, origin at
the CM of each moved geometry) and data/figures/fig_si_near_axis_nh3.{pdf,png}.
"""

import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.patches import FancyArrowPatch, Wedge

from src.figures import _savefig, _style
from src.scoring import EPS_DENOM, ON_AXIS_TOL, ModeScorer

FIG_LABEL = "fig_si_near_axis_nh3"
CSV_NAME = "si_near_axis_nh3.csv"
FRAME = ("program frame (src/frame.py) of the undistorted molecule, axis directions held fixed: z = C3 axis; "
         "origin at the CM of each moved geometry")
MODE = 3                                    # 1-based vibration number: 1796 cm-1, e bend
RANGE_TOL, RANGE_NOTOL = (1e-9, 1e-2), (1e-20, 1e-2)
PER_DEC_TOL, PER_DEC_NOTOL, N_THETA = 16, 12, 144
NO_TOL = (-1.0, 0.0)                        # (r_tol, d_tol): no N_Q, only 0/0 summands set to 0
R0 = 0.6                                    # radius of the central disk (N on the axis) in the log-polar maps
C_N = "#332288"                             # the near-axis atom (indigo; not a manuscript class colour)


# ---------------------------------------------------------------------------------------- geometry
def nh3(data_dir="data"):
    """NH3 in the program frame with N exactly on the C3 axis: coordinates (at the CM), symbols, bonds,
    masses, mode-3 vector, and the indices of N, the C3 axis and the two map axes."""
    from main import build_scorer_and_final, load_inputs
    raw, _ = load_inputs("NH3", "normal", data_dir)
    scorer, final = build_scorer_and_final(raw, "normal")
    sym = list(raw["atoms"]); m = np.asarray(scorer.masses, float)
    X = np.asarray(scorer.coords, float).copy()
    nN = [i for i, s in enumerate(sym) if s.upper() == "N"][0]
    k = int(np.argmin([np.linalg.norm(np.delete(X[nN], j)) for j in range(3)]))
    ka, kb = [j for j in range(3) if j != k]
    r_printed = float(np.linalg.norm(np.delete(X[nN], k)))
    H = [i for i in range(len(sym)) if i != nN]
    X[nN, [ka, kb]] = X[H][:, [ka, kb]].mean(axis=0)   # N over the H centroid: on the C3 axis through the CM
    X = X - (m[:, None] * X).sum(axis=0) / m.sum()
    vib = [f for f in final if str(f.get("label", "")).startswith("Vib")]
    return dict(X0=X, sym=sym, bonds=raw["bonds"], masses=m, n=nN, H=H, k=k, ka=ka, kb=kb, r_printed=r_printed,
                d=np.asarray(vib[MODE - 1]["vector"], float), freq=float(vib[MODE - 1]["frequency"]))


def place(w, rx, ry):
    """Coordinates (at their CM) with N at (rx, ry) from the C3 axis through the moved CM."""
    m = w["masses"]; X = w["X0"].copy()
    scale = 1.0 / (1.0 - m[w["n"]] / m.sum())
    X[w["n"], w["ka"]] += rx * scale; X[w["n"], w["kb"]] += ry * scale
    return X - (m[:, None] * X).sum(axis=0) / m.sum()


def rscore_pieces(X, d, k, r_tol, d_tol):
    """eq:rscore about axis k with the on-axis tolerance r_tol and the displacement floor d_tol.
    Returns (s, N - N_Q, counted, summands); a summand is 0 where it is not computed."""
    Q = np.eye(3)[k]
    rp = X - np.outer(X @ Q, Q)
    lr = np.linalg.norm(rp, axis=1); ld = np.linalg.norm(d, axis=1)
    counted = lr > r_tol; ok = counted & (ld > d_tol) & (lr > 0)
    c = np.zeros(len(X))
    c[ok] = (np.cross(rp[ok], d[ok]) @ Q) / (lr[ok] * ld[ok])
    n = int(counted.sum())
    return (c[counted].sum() / n if n else 0.0), n, counted, c


# ---------------------------------------------------------------------------------------- maps
def lp_grid(rmin, rmax, per_dec, n_theta=N_THETA):
    ndec = np.log10(rmax / rmin); nr = int(round(ndec * per_dec))
    ue = np.linspace(0.0, ndec, nr + 1); te = np.linspace(0.0, 2 * np.pi, n_theta + 1)
    return rmin * 10 ** (0.5 * (ue[1:] + ue[:-1])), 0.5 * (te[1:] + te[:-1])


def maps(w, per_dec_tol=PER_DEC_TOL, per_dec_notol=PER_DEC_NOTOL, n_theta=N_THETA):
    """Both log-polar maps in one table (column `tolerance`)."""
    rows = []
    for tol, rng, pd_, cut in (("1e-5 A", RANGE_TOL, per_dec_tol, (ON_AXIS_TOL, EPS_DENOM)),
                               ("none", RANGE_NOTOL, per_dec_notol, NO_TOL)):
        rho, th = lp_grid(*rng, pd_, n_theta)
        for r, t in [(0.0, 0.0)] + [(a, b) for b in th for a in rho]:
            X = place(w, r * np.cos(t), r * np.sin(t))
            s, n_, counted, c = rscore_pieces(X, w["d"], w["k"], *cut)
            rq = X[w["n"]]
            rows.append(dict(tolerance=tol, rho_A=r, theta_deg=np.degrees(t), r_x_A=rq[w["ka"]], r_y_A=rq[w["kb"]],
                             r_Q_N_A=float(np.hypot(rq[w["ka"]], rq[w["kb"]])),
                             contribution_N=c[w["n"]] * counted[w["n"]], s_R_z=s, N_minus_NQ=n_))
    return pd.DataFrame(rows)


def validate(w, df, n_sample=40, seed=1):
    """Own eq:rscore against ModeScorer.calculate_scores at sampled map points (with the tolerance)."""
    dev = 0.0
    for _, r in df[df.tolerance == "1e-5 A"].sample(n_sample, random_state=seed).iterrows():
        X = place(w, r.rho_A * np.cos(np.radians(r.theta_deg)), r.rho_A * np.sin(np.radians(r.theta_deg)))
        sc = ModeScorer(w["sym"], X, w["bonds"]).calculate_scores(w["d"])["R"]["xyz"[w["k"]]]
        dev = max(dev, abs(sc - r.s_R_z))
    return dev


# ---------------------------------------------------------------------------------------- figure
def centred_levels(step, lim):
    """Colour-bin edges of width `step` centred on the multiples of `step`, so 0 has a single colour."""
    k = int(round(lim / step))
    return np.round(step * (np.arange(-k, k + 2) - 0.5), 6)


LEVELS_C, LEVELS_S = centred_levels(0.1, 1.0), centred_levels(0.05, 0.4)


def _axes_marker(ax, x0, y0, right, up, L):
    kw = dict(arrowstyle="-|>", lw=0.6, color="0.3", mutation_scale=6, shrinkA=0, shrinkB=0)
    ax.annotate("", xy=(x0 + L, y0), xytext=(x0, y0), arrowprops=kw)
    ax.annotate("", xy=(x0, y0 + L), xytext=(x0, y0), arrowprops=kw)
    ax.text(x0 + L * 1.12, y0, right, va="center", ha="left", color="0.3")
    ax.text(x0, y0 + L * 1.12, up, va="bottom", ha="center", color="0.3")


def end_view(ax, w):
    """(a) NH3 viewed along the C3 axis: N on the axis with its displacement, the H atoms (not to scale)
    with their displacement projections, and a ring giving N's contribution for each direction of offset."""
    E = np.eye(3); Q = E[w["k"]]; ka, kb, n, d = w["ka"], w["kb"], w["n"], w["d"]
    dh = d[n] / np.linalg.norm(d[n])
    cm = plt.get_cmap("PuOr"); norm = BoundaryNorm(LEVELS_C, cm.N)
    ax.set_aspect("equal"); ax.axis("off")
    R_H = 1.08
    for h in w["H"]:
        v = w["X0"][h, [ka, kb]] / np.linalg.norm(w["X0"][h, [ka, kb]])
        ax.plot([0, R_H * v[0]], [0, R_H * v[1]], color="0.75", lw=0.8, zorder=0)
        ax.add_patch(plt.Circle(R_H * v, 0.08, color="0.65", zorder=2))
        dp = d[h, [ka, kb]] / np.linalg.norm(d[h])
        if np.linalg.norm(dp) > 1e-3:
            ax.add_patch(FancyArrowPatch(R_H * v, R_H * v + 0.34 * dp, arrowstyle="-|>", mutation_scale=6, lw=0.8,
                                         color="0.3", shrinkA=3, shrinkB=0, zorder=3))
    r1, r2 = 0.5, 0.72
    for t in np.arange(0.0, 360.0, 2.0):
        u = np.cos(np.radians(t + 1)) * E[ka] + np.sin(np.radians(t + 1)) * E[kb]
        ax.add_patch(Wedge((0, 0), r2, t, t + 2.05, width=r2 - r1, color=cm(norm(float(np.cross(u, dh) @ Q))),
                           lw=0, zorder=1))
    for r in (r1, r2):
        ax.add_patch(plt.Circle((0, 0), r, fill=False, lw=0.4, color="0.5", zorder=1))
    for t in (0, 90, 180, 270):
        u = np.cos(np.radians(t)) * E[ka] + np.sin(np.radians(t)) * E[kb]; c = float(np.cross(u, dh) @ Q)
        p = 0.61 * np.array([np.cos(np.radians(t)), np.sin(np.radians(t))])
        ax.text(*p, f"{c:+.0f}".replace("+0", "0").replace("-0", "0"), ha="center", va="center", zorder=4,
                color="white" if abs(c) > 0.6 else "black")
    ax.plot(0, 0, marker="o", ms=7, mfc=C_N, mec="black", mew=0.6, zorder=3)
    du = dh[[ka, kb]] / np.linalg.norm(dh[[ka, kb]])
    ax.add_patch(FancyArrowPatch((0, 0), 0.44 * du, arrowstyle="-|>", mutation_scale=7, lw=1.1, color="black",
                                 shrinkA=4, shrinkB=0, zorder=5))
    pp = np.array([-du[1], du[0]]); pp = pp if pp[1] <= 0 else -pp
    ax.text(*(0.26 * du + 0.17 * pp), r"$\boldsymbol{d}^{\mathrm{N}}$", ha="center", va="center", zorder=5)
    ax.text(0.1, 0.12, "N", color=C_N, ha="left", va="bottom", zorder=5)
    _axes_marker(ax, -1.35, 0.85, "xyz"[ka], "xyz"[kb], L=0.28)
    ax.set_xlim(-1.45, 1.45); ax.set_ylim(-1.45, 1.45)


def lp_panel(ax, M, value, rmin, rmax, cmap, levels, cut, label_every, h_angles, axes_lab):
    """One log-polar map: radius R0 + log10(rho/rmin) at the true direction; central disk = N on the axis;
    white decade rings (every decade for short ranges, else every `label_every`; labelled every `label_every`); dashed circle = the tolerance; grey lines = H directions."""
    cm = plt.get_cmap(cmap); norm = BoundaryNorm(levels, cm.N)
    G = M[M.rho_A > 0]; z0 = M[M.rho_A == 0][value].to_numpy()
    piv = G.pivot_table(index="theta_deg", columns="rho_A", values=value)
    rho = piv.columns.to_numpy(); th = np.radians(piv.index.to_numpy())
    uc = np.log10(rho / rmin); du = uc[1] - uc[0]; ue = np.r_[uc - du / 2, uc[-1] + du / 2]
    dt = th[1] - th[0]; te = np.r_[th - dt / 2, th[-1] + dt / 2]
    Dg, Tg = np.meshgrid(R0 + ue, te)
    im = ax.pcolormesh(Dg * np.cos(Tg), Dg * np.sin(Tg), piv.to_numpy(), cmap=cm, norm=norm, shading="flat",
                       rasterized=True, zorder=1)
    if len(z0):
        ax.add_patch(plt.Circle((0, 0), R0, color=cm(norm(z0[0])), lw=0, zorder=2))
        ax.add_patch(plt.Circle((0, 0), R0, fill=False, lw=0.4, color="0.35", zorder=3))
    Dmax = R0 + np.log10(rmax / rmin)
    for ang in h_angles:
        a = np.radians(ang)
        ax.plot([R0 * np.cos(a), Dmax * np.cos(a)], [R0 * np.sin(a), Dmax * np.sin(a)], color="0.45", lw=0.5, zorder=5)
    e0, e1 = int(np.ceil(np.log10(rmin) - 1e-9)), int(np.floor(np.log10(rmax) + 1e-9))
    ring_step = 1 if e1 - e0 <= 10 else label_every          # long ranges: rings only where labelled
    for e in range(e1, e0 - 1, -ring_step):
        D = R0 + e - np.log10(rmin)
        ax.add_patch(plt.Circle((0, 0), D, fill=False, lw=0.35, color="white", alpha=0.8, zorder=5))
        if (e1 - e) % label_every == 0:
            a = np.radians(112.5)
            ax.text(D * np.cos(a), D * np.sin(a), f"$10^{{{e}}}$", fontsize=6.5, ha="center", va="center", zorder=7,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.8))
    if cut:
        ax.add_patch(plt.Circle((0, 0), R0 + np.log10(cut / rmin), fill=False, ls="--", lw=0.9, color="#555555", zorder=6))
    ax.add_patch(plt.Circle((0, 0), Dmax, fill=False, lw=0.6, color="0.2", zorder=6))
    _axes_marker(ax, 0.8 * Dmax, 0.8 * Dmax, axes_lab[0], axes_lab[1], L=0.16 * Dmax)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-1.03 * Dmax, 1.03 * Dmax); ax.set_ylim(-1.03 * Dmax, 1.03 * Dmax)
    return im


def plot(df, w, out_dir, label=FIG_LABEL):
    """(a) end view; (b, c) contribution of N and s[R_z] without tolerance; (d, e) the same with 1e-5 A."""
    _style()
    plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
                         "axes.titlesize": 8.5})
    X0, ka, kb = w["X0"], w["ka"], w["kb"]
    h_ang = [float(np.degrees(np.arctan2(X0[h, kb], X0[h, ka]))) for h in w["H"]]
    fig = plt.figure(figsize=(7.2, 4.6))
    gs = fig.add_gridspec(2, 6, width_ratios=[1.0, 1, 0.045, 0.22, 1, 0.045], wspace=0.06, hspace=0.16)
    ax = fig.add_subplot(gs[:, 0]); end_view(ax, w); ax.set_title("(a)", loc="left")
    D_no, D_tol = df[df.tolerance == "none"], df[df.tolerance == "1e-5 A"]
    panels = [(0, 1, D_no, "contribution_N", "PuOr", LEVELS_C, RANGE_NOTOL, None, 3, "(b)", "contribution of N", 0.5),
              (0, 4, D_no, "s_R_z", "BrBG", LEVELS_S, RANGE_NOTOL, None, 3, "(c)", r"$s[\mathrm{R}_z]$", 0.2),
              (1, 1, D_tol, "contribution_N", "PuOr", LEVELS_C, RANGE_TOL, ON_AXIS_TOL, 3, "(d)", "contribution of N", 0.5),
              (1, 4, D_tol, "s_R_z", "BrBG", LEVELS_S, RANGE_TOL, ON_AXIS_TOL, 3, "(e)", r"$s[\mathrm{R}_z]$", 0.2)]
    for row, col, D, val, cmap, lev, rng, cut, every, tag, cbl, tstep in panels:
        ax = fig.add_subplot(gs[row, col])
        im = lp_panel(ax, D, val, rng[0], rng[1], cmap, lev, cut, every, h_ang, ("xyz"[ka], "xyz"[kb]))
        ax.set_title(tag, loc="left")
        lim = 1.0 if tstep == 0.5 else 0.4
        cb = fig.colorbar(im, cax=fig.add_subplot(gs[row, col + 1]), ticks=np.round(np.arange(-lim, lim + 1e-9, tstep), 1))
        cb.set_label(cbl, labelpad=2); cb.outline.set_linewidth(0.6)
    paths = _savefig(fig, out_dir, label)
    plt.close(fig)
    return paths


# ---------------------------------------------------------------------------------------- run
def summary(df, w, dev):
    tol, no = df[df.tolerance == "1e-5 A"], df[df.tolerance == "none"]
    inside, outside = tol[tol.r_Q_N_A <= ON_AXIS_TOL], tol[tol.r_Q_N_A > ON_AXIS_TOL]
    off = no[no.rho_A > 0]
    return dict(
        mode=f"NH3 mode {MODE}, {w['freq']:.0f} cm-1, about the C3 axis ({'xyz'[w['k']]})",
        N_printed_offset_A=w["r_printed"], max_dev_vs_ModeScorer=dev,
        inside_tolerance_s=(float(inside.s_R_z.min()), float(inside.s_R_z.max())),
        outside_tolerance_s=(float(outside.s_R_z.min()), float(outside.s_R_z.max())),
        no_tolerance_contribution=(float(off.contribution_N.min()), float(off.contribution_N.max())),
        no_tolerance_s=(float(off.s_R_z.min()), float(off.s_R_z.max())),
        points=(len(tol), len(no)))


def run(data_dir="data", out_dir="data/figures", results_dir="data/results"):
    w = nh3(data_dir)
    df = maps(w)
    dev = validate(w, df)
    os.makedirs(results_dir, exist_ok=True)
    csv_path = os.path.join(results_dir, CSV_NAME)
    df.assign(frame="program").to_csv(csv_path, index=False, float_format="%.8g")   # frame: see FRAME
    pdf, png = plot(df, w, out_dir)
    return dict(csv=csv_path, pdf=pdf, png=png, frame=FRAME, summary=summary(df, w, dev))


if __name__ == "__main__":
    out = run()
    for k, v in out.items():
        print(f"{k}: {v}")
