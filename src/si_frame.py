"""SI analyses of the principal-axis frame (JCC revision, referee comments R1.5 and R1.6).

    py -m src.si_frame            # from Github/scoring-functions/

Inputs: the tracked normal-mode intermediates of every roster molecule
(data/mol_list_method.csv, mol_type in ideal / non-ideal / test / multi-centre),
read with main.load_inputs(mol, "normal"). The frame is the program's own
(src/frame.py: select_frame, with the same inputs as main.build_scorer_and_final:
COM-shifted input coordinates and the input normal vibrational modes, no ideal
modes); "allowed frames" are src.frame.candidate_frames (symmetric and spherical
tops) or the single frame (asymmetric tops: eigenvectors; linear molecules: the
distances from the axes do not depend on the turn about the molecular axis).

R1.5 -- frame spread (fig_si_frame_spread, si_frame_spread.csv)
    For a normal vibrational mode k, K in {T, R} and a set of frames F,
        Delta_k[K] = max over the rotatable axes Q of ( max_F |s_k[K_Q]| - min_F |s_k[K_Q]| ).
    |s_k[K_Q]| depends only on the direction of the axis Q (eq:tscore, eq:rscore), so
    the frames are sampled through their axis directions u:
      asymmetric  Q = x, y, z; the 8 signed principal frames (u = +-axis): Delta = 0
      linear      Q = y, z;    turns about the molecular axis, 1-deg steps
      symmetric   Q = x, y;    turns about the C_n axis z, 1-deg steps
      spherical   Q = x, y, z; N_SPHERE Fibonacci directions + the directions of the
                  atoms + rings of directions 0.01-5 deg around them + the axes of the
                  allowed frames (the side report sampled 4000 directions without the
                  rings; its spherical-top spreads are lower by up to ~0.01)
    spread_allowed: the same over the allowed frames only (tops). Values whose mode has
    max |s_k[K_Q]| < MIN_ABS in every frame are kept in the CSV but not plotted.

R1.6 -- on-axis cut-off (fig_si_onaxis_gap, fig_si_onaxis_displacement,
    si_frame_onaxis_pairs.csv, si_frame_onaxis_displacement.csv)
    Atom A is ON-AXIS (BY SYMMETRY) for the Q-axis (through the CM) if
      (i)   a rotation C_n (n = 2..8) about Q maps the molecule onto itself and A onto itself, or
      (ii)  A is at the CM by symmetry: fixed by inversion, by rotations about two
            non-parallel axes, or by a C_n(u) and the mirror plane perpendicular to u, or
      (iii) the molecule is linear and Q is the molecular axis;
    otherwise OFF-AXIS (BY SYMMETRY). The class never uses the cut-off itself; the
    figures test that the cut-off ON_AXIS_TOL (src.scoring) separates the two classes.
    The displacement figure adds |d^A| of every atom in every normal vibrational mode
    against the floor EPS_DENOM: on-axis pairs are not counted (N_Q); an off-axis atom
    with |d^A| <= EPS_DENOM counts in N - N_Q but adds a zero term (unit(0) := 0).

Ported from the side-report scripts r15_frame_spread.py / r15_cases_spread.py and
r16_symmetric_on_axis.py / r16_cutoff_gap.py / r16_selected_frames.py /
r16_displacement_floor.py; the frame functions are src/frame.py's.
"""

import os
import time

import numpy as np
import pandas as pd

from . import frame as FR
from .scoring import EPS_DENOM, ON_AXIS_TOL, ModeScorer

MOL_TYPES = ("ideal", "non-ideal", "test", "multi-centre")
# data set of each mol_type (the manuscript's single-centre set = ideal + non-ideal)
GROUP = {"ideal": "single-centre", "non-ideal": "single-centre", "test": "test", "multi-centre": "benzene"}
CASES = ("asymmetric", "linear", "symmetric", "spherical")
TURN_STEP_DEG = 1.0      # symmetric tops and linear molecules: turn grid
N_SPHERE = 20000         # spherical tops: Fibonacci directions ...
RING_RADII_DEG = (0.01, 0.03, 0.1, 0.3, 1.0, 2.0, 3.0, 5.0)   # ... plus rings around each atom direction
RING_AZIMUTHS = 72
MIN_ABS = 0.01           # plotted: values whose mode reaches |s_k[K_Q]| >= MIN_ABS in some frame
FREQ_TOL = 0.05          # cm-1: modes this close form one degenerate set (CSV column set_size)
ZERO_TOL = 1e-12         # r_Q at or below this counts as exactly 0 (drawn on the "0" tick)


# ---------------------------------------------------------------- inputs
def roster(data_dir="data"):
    r = pd.read_csv(os.path.join(data_dir, "mol_list_method.csv"))
    return r[r["mol_type"].isin(MOL_TYPES)].reset_index(drop=True)


def molecule_inputs(mol, data_dir="data"):
    """COM-shifted input coordinates, symbols, masses, bonds, normal-mode vectors and frequencies,
    the top kind and the program's frame (src/frame.py, inputs as in main.build_scorer_and_final)."""
    from main import load_inputs
    raw, _ = load_inputs(mol, "normal", data_dir=data_dir)
    sc = ModeScorer(raw["atoms"], raw["coords"], raw["bonds"])
    X, sym, m, bonds = sc.coords.copy(), list(raw["atoms"]), sc.masses, raw["bonds"]
    vecs = [np.asarray(md["vector"], dtype=float) for md in raw["modes"]]
    freqs = np.array([md["frequency"] for md in raw["modes"]], dtype=float)
    kind, z, n = FR.detect_top(X, sym, m)
    F, info = FR.select_frame(X, sym, m, bonds, vecs, keep_input_frame=False)
    return dict(molecule=mol, X=X, sym=sym, m=m, bonds=bonds, vecs=vecs, freqs=freqs,
                kind=kind, z=z, n=n, F=F, info=info)


def allowed_frames(d):
    """The allowed frames and the index of the selected one among them."""
    if d["kind"] in ("symmetric", "spherical"):
        frames = FR.candidate_frames(d["X"], d["sym"], d["m"], d["kind"], d["z"])
        idx = [i for i, f in enumerate(frames) if np.allclose(f, d["F"], atol=1e-8)]
        if len(idx) != 1:
            raise RuntimeError(f"{d['molecule']}: selected frame matches {len(idx)} allowed frames")
        return frames, idx[0]
    return [d["F"]], 0


def set_sizes(freqs, tol=FREQ_TOL):
    return np.array([int(np.sum(np.abs(freqs - f) < tol)) for f in freqs])


# ---------------------------------------------------------------- per-axis scores along any direction
def axis_scores(X, vecs, U):
    """s_k[T_u] and s_k[R_u] (eq:tscore, eq:rscore as src.scoring.ModeScorer evaluates them)
    of every mode k for every unit direction u (rows of U), all in the coordinates of X.
    For u = an axis of a frame f this equals ModeScorer on X @ f.T and d @ f.T.
    Returns (T, R), each of shape (n_modes, n_directions)."""
    U = np.atleast_2d(np.asarray(U, dtype=float))
    n = len(X)
    rp = X[:, None, :] - (X @ U.T)[:, :, None] * U[None, :, :]      # r_perp, N x M x 3
    lr = np.linalg.norm(rp, axis=2)                                  # N x M
    off = lr > ON_AXIS_TOL                                           # N_Q: on the axis -> excluded
    n_off = off.sum(axis=0)
    T = np.zeros((len(vecs), len(U)))
    R = np.zeros((len(vecs), len(U)))
    for k, d in enumerate(vecs):
        ld = np.linalg.norm(d, axis=1)
        move = ld > EPS_DENOM                                        # unit(0) := 0
        unit_d = np.zeros_like(d)
        unit_d[move] = d[move] / ld[move, None]
        T[k] = U @ unit_d.sum(axis=0) / n
        ok = off & move[:, None]
        num = np.einsum("nmc,mc->nm", np.cross(rp, d[:, None, :]), U)
        term = np.where(ok, num / np.where(ok, lr * ld[:, None], 1.0), 0.0)
        R[k] = np.where(n_off > 0, term.sum(axis=0) / np.maximum(n_off, 1), 0.0)
    return T, R


def fibonacci_sphere(n):
    g = np.arange(n) + 0.5
    ph = np.arccos(1 - 2 * g / n)
    th = np.pi * (1 + 5 ** 0.5) * g
    return np.c_[np.sin(ph) * np.cos(th), np.sin(ph) * np.sin(th), np.cos(ph)]


def _turn_dirs(a, b, step_deg=TURN_STEP_DEG):
    phi = np.radians(np.arange(0.0, 360.0, step_deg))
    return np.outer(np.cos(phi), a) + np.outer(np.sin(phi), b)


def spread_directions(d, n_sphere=N_SPHERE):
    """Axis directions that sample every principal frame of the molecule, grouped per rotatable axis
    (a list of direction arrays; the spread is the largest over the groups of max - min |s|)."""
    F, kind, X = d["F"], d["kind"], d["X"]
    if kind == "asymmetric":
        return [np.array([F[q], -F[q]]) for q in range(3)]           # the 8 signed frames, per axis
    if kind == "linear":
        return [_turn_dirs(F[1], F[2])]                               # y and z sweep the same circle
    if kind == "symmetric":
        return [_turn_dirs(F[0], F[1])]                               # x and y sweep the same circle
    frames, _ = allowed_frames(d)
    atoms = [x / np.linalg.norm(x) for x in X if np.linalg.norm(x) > 1e-2]
    # Near an atom direction the atom's term of s[R_u] depends on the azimuth of approach, so the
    # extremes sit just beside it: rings around +-each atom direction (the lattice alone converges
    # to ~0.01 only; with the rings 2e-4 between 2e4 and 1e5 lattice points).
    caps = [_rings(s * u) for u in atoms for s in (1, -1)]
    U = np.vstack([fibonacci_sphere(n_sphere)] + atoms + caps + [f for f in frames])
    return [U]                                                        # every axis sweeps every direction


def _rings(u, radii_deg=RING_RADII_DEG, n_az=RING_AZIMUTHS):
    a = np.eye(3)[int(np.argmin(np.abs(u)))]
    e1 = np.cross(u, a)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(u, e1)
    az = np.linspace(0, 2 * np.pi, n_az, endpoint=False)
    circ = np.outer(np.cos(az), e1) + np.outer(np.sin(az), e2)
    return np.vstack([np.cos(r) * u + np.sin(r) * circ for r in np.radians(radii_deg)])


# ---------------------------------------------------------------- R1.5: frame spread
def frame_spread_rows(d, mol_type, point_group, n_sphere=N_SPHERE):
    groups = spread_directions(d, n_sphere)
    TR = [axis_scores(d["X"], d["vecs"], U) for U in groups]
    kind = d["kind"]
    if kind in ("symmetric", "spherical"):
        frames, _ = allowed_frames(d)
        qs = (0, 1) if kind == "symmetric" else (0, 1, 2)
        TA = [axis_scores(d["X"], d["vecs"], np.array([f[q] for f in frames])) for q in qs]
    sizes = set_sizes(d["freqs"])
    rows = []
    for k in range(len(d["vecs"])):
        for j, K in enumerate("TR"):
            A = [np.abs(tr[j][k]) for tr in TR]
            spread = max(a.max() - a.min() for a in A)
            row = dict(molecule=d["molecule"], mol_type=mol_type, group=GROUP[mol_type],
                       point_group=point_group, case=kind, mode=k + 1, freq=d["freqs"][k],
                       set_size=int(sizes[k]), kind=K, spread=float(spread),
                       max_abs=float(max(a.max() for a in A)), spread_allowed=np.nan,
                       n_directions=int(sum(len(a) for a in A)), n_allowed=np.nan)
            if kind in ("symmetric", "spherical"):
                AA = [np.abs(ta[j][k]) for ta in TA]
                row.update(spread_allowed=float(max(a.max() - a.min() for a in AA)), n_allowed=len(frames))
            rows.append(row)
    return rows


# ---------------------------------------------------------------- R1.6: on-axis by symmetry
def perm_of(X, sym, A, tol=None):
    """Atom permutation induced by the operation A (3x3), or None if A is not a symmetry operation within tol."""
    tol = FR.SYM_TOL if tol is None else tol
    Y = X @ A.T
    other = np.array(sym, dtype=object)
    out = []
    for i in range(len(X)):
        dist = np.linalg.norm(X - Y[i], axis=1)
        dist[other != sym[i]] = np.inf
        j = int(np.argmin(dist))
        if dist[j] > tol or j in out:
            return None
        out.append(j)
    return out


def _unit(v):
    return v / np.linalg.norm(v)


def rotation_axes(X, sym, m, frames):
    """Every C_n axis (n = 8..2, the largest n per line) among the candidate lines: the frame rows,
    the C2 axes, the atom directions and the inertia eigenvectors. Returns [(u, n, permutation)]."""
    cand = [f[k] for f in frames for k in range(3)] + FR.c2_axes(X, sym, FR.SYM_TOL)
    cand += [_unit(x) for x in X if np.linalg.norm(x) > 1e-2]
    cand += list(np.linalg.eigh(FR.inertia(X, m))[1].T)
    axes = []
    for u in cand:
        u = _unit(u)
        if any(abs(abs(u @ a) - 1) < 1e-6 for a, _, _ in axes):
            continue
        for n in range(8, 1, -1):
            p = perm_of(X, sym, FR.rotation(u, 2 * np.pi / n))
            if p is not None:
                axes.append((u, n, p))
                break
    return axes


def at_cm_by_symmetry(X, sym, axes):
    N = len(X)
    fixed = {i: [u for u, _, p in axes if p[i] == i] for i in range(N)}
    inv = perm_of(X, sym, -np.eye(3))
    out = set()
    for i in range(N):
        if inv is not None and inv[i] == i:
            out.add(i)
            continue
        us = fixed[i]
        if any(abs(abs(a @ b) - 1) > 1e-6 for a in us for b in us):
            out.add(i)
            continue
        for u in us:                                                  # C_n(u) + sigma perpendicular to u
            s = perm_of(X, sym, np.eye(3) - 2 * np.outer(u, u))
            if s is not None and s[i] == i:
                out.add(i)
                break
    return out


def on_axis_by_symmetry(N, kind, axis_lin, axes, cm_set, q):
    """(atoms on the line through the CM along q by symmetry, reason)."""
    if kind == "linear" and abs(abs(q @ axis_lin) - 1) < 1e-6:
        return set(range(N)), "linear"
    on, why = set(cm_set), "at the CM"
    for u, n, p in axes:
        if abs(abs(u @ q) - 1) < 1e-6:
            on |= {i for i in range(N) if p[i] == i}
            why = f"C{n} about the axis"
    return on, why


def onaxis_pair_rows(d, mol_type, point_group):
    """One row per (allowed frame, axis, atom); `selected` marks the program's frame."""
    X, sym, m, kind = d["X"], d["sym"], d["m"], d["kind"]
    frames, isel = allowed_frames(d)
    axis_lin = d["z"] if kind == "linear" else None                   # detect_top: molecular axis
    axes = rotation_axes(X, sym, m, list(frames) + [np.eye(3)])       # Gaussian's axes are candidates too
    cm_set = at_cm_by_symmetry(X, sym, axes)
    rows = []
    for fi, f in enumerate(frames):
        Y = X @ f.T
        for k, q in enumerate("xyz"):
            on, why = on_axis_by_symmetry(len(X), kind, axis_lin, axes, cm_set, f[k])
            r_q = np.linalg.norm(np.delete(Y, k, axis=1), axis=1)
            for a in range(len(X)):
                rows.append(dict(molecule=d["molecule"], mol_type=mol_type, point_group=point_group,
                                 kind=kind, frame_index=fi, selected=fi == isel, n_allowed=len(frames),
                                 axis=q, atom=f"{sym[a]}{a + 1}", atom_index=a, r_Q=float(r_q[a]),
                                 r_par=float(abs(Y[a, k])), on_by_symmetry=a in on,
                                 at_cm_by_symmetry=a in cm_set,
                                 reason=("at the CM" if a in cm_set else why) if a in on else ""))
    return rows


def displacement_rows(d, pairs):
    """One row per (normal vibrational mode, atom, axis) in the selected frame."""
    sel = [p for p in pairs if p["selected"]]
    rows = []
    for k, v in enumerate(d["vecs"]):
        ld = np.linalg.norm(v, axis=1)
        for p in sel:
            a = p["atom_index"]
            region = ("not counted" if p["r_Q"] <= ON_AXIS_TOL
                      else "term 0" if ld[a] <= EPS_DENOM else "full term")
            rows.append(dict(molecule=p["molecule"], mol_type=p["mol_type"], point_group=p["point_group"],
                             kind=p["kind"], mode=k + 1, freq=d["freqs"][k], axis=p["axis"], atom=p["atom"],
                             atom_index=a, d=float(ld[a]), r_Q=p["r_Q"],
                             on_by_symmetry=p["on_by_symmetry"], region=region))
    return rows


# ---------------------------------------------------------------- figures
C_ON, C_OFF = "#E69F00", "#332288"        # on-/off-axis (by symmetry): colours the manuscript palette does not use
GROUP_COLORS = {"single-centre": "#E69F00", "test": "#0072B2", "benzene": "#009E73"}   # Okabe-Ito, as the side report
CASE_SHORT = {"asymmetric": "asym.", "linear": "linear", "symmetric": "symm.", "spherical": "spher."}
RQ = r"$|\boldsymbol{r}_{Q}^{A}|$"
DA = r"$|\boldsymbol{d}^{A}|$"
Z0 = 1e-8                                 # exact zeros are drawn here, on a tick labelled 0


def _figure_spread(S, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .figures import COLORS, _savefig, _style
    _style()
    V = S[S.max_abs >= MIN_ABS]
    tops = V[V.case.isin(["symmetric", "spherical"])]
    fig, axs = plt.subplots(1, 2, figsize=(10.0, 4.3), gridspec_kw=dict(width_ratios=[1.15, 1.0], wspace=0.32))

    ax = axs[0]
    bins = np.linspace(0, 0.75, 31)
    bottom = np.zeros(len(bins) - 1)
    for g, c in GROUP_COLORS.items():
        h, _ = np.histogram(tops[tops.group == g].spread, bins=bins)
        n_mol = tops[tops.group == g].molecule.nunique()
        ax.bar(bins[:-1], h, width=np.diff(bins), bottom=bottom, align="edge", color=c, edgecolor="white",
               lw=0.4, label=f"{g} ({n_mol})")
        bottom += h
    med, p90 = tops.spread.median(), tops.spread.quantile(0.9)
    top = bottom.max() * 1.25
    ax.axvline(med, color="black", lw=1.0)
    ax.axvline(p90, color="black", lw=1.0, ls="--")
    ax.text(med + 0.008, top * 0.97, f"median\n{med:.2f}", va="top", fontsize=9)
    ax.text(p90 + 0.008, top * 0.97, f"90th percentile\n{p90:.2f}", va="top", fontsize=9)
    ax.set_xlim(0, 0.75)
    ax.set_ylim(0, top)
    ax.set_xlabel(r"frame spread $\Delta_k[\mathrm{T}]$ or $\Delta_k[\mathrm{R}]$")
    ax.set_ylabel(r"number of values $\Delta_k[\mathrm{T}]$, $\Delta_k[\mathrm{R}]$")
    ax.legend(frameon=False, loc="center right", fontsize=9, bbox_to_anchor=(1.0, 0.55))
    ax.set_title("(a) symmetric and spherical tops; all principal-axis frames", loc="left", fontsize=10)

    ax = axs[1]
    rng = np.random.default_rng(0)
    grey = COLORS["threshold"]
    for i, c in enumerate(CASES):
        v = V[V.case == c].spread.to_numpy()
        bp = ax.boxplot(v, positions=[i], widths=0.5, patch_artist=True, showfliers=False,
                        medianprops=dict(color="black", lw=1.2))
        for bx in bp["boxes"]:
            bx.set(facecolor="white", edgecolor=grey)
        ax.scatter(i + rng.uniform(-0.17, 0.17, len(v)), v, s=3, color=grey, alpha=0.35, lw=0, zorder=3,
                   rasterized=True)
        ax.text(i, 0.79, f"n = {len(v)}", ha="center", fontsize=9)
    nmol = S.groupby("case").molecule.nunique()
    ax.set_xticks(range(len(CASES)))
    ax.set_xticklabels([f"{CASE_SHORT[c]}\n({nmol.get(c, 0)})" for c in CASES])
    ax.set_xlim(-0.6, len(CASES) - 0.4)
    ax.set_ylim(-0.02, 0.86)
    ax.set_xlabel("case (number of molecules)")
    ax.set_ylabel(r"frame spread $\Delta_k[\mathrm{T}]$, $\Delta_k[\mathrm{R}]$")
    ax.set_title(f"(b) all {S.molecule.nunique()} molecules; all principal-axis frames", loc="left", fontsize=10)
    return _savefig(fig, out_dir, "fig_si_frame_spread")


def _sci(v):
    e = int(np.floor(np.log10(v)))
    return rf"${v / 10 ** e:.1f}\times10^{{{e}}}$"


def _strip(ax, P, rng, title):
    from .figures import COLORS
    on = P[P.on_by_symmetry]
    off = P[~P.on_by_symmetry]
    on_nz = on[on.r_Q > ZERO_TOL].r_Q.to_numpy()
    n_zero = int((on.r_Q <= ZERO_TOL).sum())
    for yi, (v, col) in enumerate([(on_nz, C_ON), (off.r_Q.to_numpy(), C_OFF)]):
        ax.scatter(v, yi + rng.uniform(-0.27, 0.27, len(v)), s=6, color=col, lw=0, alpha=0.5, rasterized=True)
    lo, hi = on.r_Q.max(), off.r_Q.min()
    ax.axvspan(lo, hi, color=COLORS["background"], alpha=0.25, lw=0, zorder=0)
    ax.axvline(ON_AXIS_TOL, color=COLORS["threshold"], ls="--", lw=1.0, zorder=1)
    ax.set_xscale("log")
    ax.set_xlim(2e-8, 30)
    ax.set_ylim(-0.6, 1.6)
    ax.invert_yaxis()
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["on-axis\n(by symmetry)", "off-axis\n(by symmetry)"])
    ax.tick_params(axis="y", length=0)
    ax.text(20, 0, f"{len(on)} pairs\n({n_zero} with " + RQ + " = 0, not shown)", fontsize=9,
            ha="right", va="center")
    ax.text(2.6e-8, 1, f"{len(off)} pairs", fontsize=9, ha="left", va="center")
    ax.text(lo / 1.1, 0.36, _sci(lo) + " Å", fontsize=9, ha="right", va="top", color=C_ON)
    ax.text(hi * 1.1, 0.64, _sci(hi) + " Å", fontsize=9, ha="left", va="bottom", color=C_OFF)
    ax.set_title(title, loc="left", fontsize=11)


def _figure_gap(P, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from .figures import COLORS, _savefig, _style
    _style()
    fig, axs = plt.subplots(2, 1, figsize=(8.0, 5.4), sharex=True, gridspec_kw=dict(hspace=0.35))
    rng = np.random.default_rng(3)
    nm = P.molecule.nunique()
    sel = P[P.selected]
    nf = P.groupby(["molecule", "frame_index"]).ngroups
    _strip(axs[0], sel, rng, f"(a) frame chosen by the program; {nm} molecules")
    _strip(axs[1], P, rng, f"(b) every frame allowed by the symmetry rule; {nf} frames")
    axs[1].set_xlabel("distance from the $Q$-axis, " + RQ + " (Å)")
    fig.legend(handles=[Line2D([], [], color=COLORS["threshold"], ls="--", lw=1.0,
                                  label=r"cut-off $10^{-5}$ Å"),
                           Patch(facecolor=COLORS["background"], alpha=0.25, label="gap between the classes")],
               frameon=False, loc="upper center", bbox_to_anchor=(0.55, 0.0), fontsize=9, ncol=2)
    return _savefig(fig, out_dir, "fig_si_onaxis_gap")


def _zero_log_axis(ax, which, lim):
    from matplotlib.ticker import FixedLocator, NullLocator
    ticks = [Z0] + [10.0 ** k for k in range(-6, 1, 2)]
    labs = ["0"] + [f"$10^{{{k}}}$" for k in range(-6, 1, 2)]
    a = ax.xaxis if which == "x" else ax.yaxis
    a.set_major_locator(FixedLocator(ticks))
    a.set_minor_locator(NullLocator())
    (ax.set_xticklabels if which == "x" else ax.set_yticklabels)(labs)
    (ax.set_xlim if which == "x" else ax.set_ylim)(*lim)


def _axis_break(ax, which, at):
    """Two short slashes on the axis between the 0 tick and the log scale."""
    kw = dict(color="black", lw=0.8, clip_on=False, zorder=5)
    if which == "x":
        tr = ax.get_xaxis_transform()
        for dx in (0.88, 1.12):
            ax.plot([at * dx / 1.25, at * dx * 1.25], [-0.015, 0.015], transform=tr, **kw)
    else:
        tr = ax.get_yaxis_transform()
        for dy in (0.88, 1.12):
            ax.plot([-0.015, 0.015], [at * dy / 1.25, at * dy * 1.25], transform=tr, **kw)


def _figure_displacement(D, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .figures import COLORS, _savefig, _style
    _style()
    fig, ax = plt.subplots(figsize=(6.2, 5.6))
    rng = np.random.default_rng(5)
    xlim, ylim = (4e-9, 3.0), (4e-9, 30.0)
    for on, col in [(False, C_OFF), (True, C_ON)]:
        g = D[D.on_by_symmetry == on]
        x = np.where(g.d > 0, g.d, Z0 * 10 ** rng.uniform(-0.15, 0.15, len(g)))
        y = np.where(g.r_Q > ZERO_TOL, g.r_Q, Z0 * 10 ** rng.uniform(-0.15, 0.15, len(g)))
        ax.scatter(x, y, s=4, color=col, lw=0, alpha=0.35, rasterized=True,
                   label="on-axis (by symmetry)" if on else "off-axis (by symmetry)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    _zero_log_axis(ax, "x", xlim)
    _zero_log_axis(ax, "y", ylim)
    _axis_break(ax, "x", 3.5e-8)
    _axis_break(ax, "y", 3.5e-8)
    thr = COLORS["threshold"]
    ax.axhline(ON_AXIS_TOL, color=thr, ls="--", lw=1.0)
    # the displacement floor matters only for atoms that are counted (above the on-axis cut-off)
    ax.plot([EPS_DENOM, EPS_DENOM], [ON_AXIS_TOL, ylim[1]], color=thr, ls="--", lw=1.0)
    ax.text(2.0e-4, 2.2e-6, "not counted", ha="center", va="center", color=thr, fontsize=10)
    ax.text(1.3e-7, 3e-3, "term 0", ha="center", va="center", color=thr, fontsize=10)
    ax.text(3e-5, 3e-3, "full term", ha="center", va="center", color=thr, fontsize=10)
    ax.text(xlim[1] * 0.8, ON_AXIS_TOL / 1.3, r"$10^{-5}$ Å", ha="right", va="top", color=thr, fontsize=9)
    ax.text(EPS_DENOM * 1.2, ylim[1] * 0.7, r"$10^{-6}$", ha="left", va="top", color=thr, fontsize=9)
    ax.set_xlabel("displacement of the atom, " + DA)
    ax.set_ylabel("distance from the $Q$-axis, " + RQ + " (Å)")
    lg = ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2, markerscale=3,
                   fontsize=9, handletextpad=0.2, columnspacing=1.0)
    for h in lg.legend_handles:
        h.set_alpha(1)
    nm, nk = D.molecule.nunique(), D.groupby(["molecule", "mode"]).ngroups
    ax.set_title(f"frame chosen by the program; {nk} normal vibrational modes of {nm} molecules",
                 loc="left", fontsize=10, pad=24)
    return _savefig(fig, out_dir, "fig_si_onaxis_displacement")


# ---------------------------------------------------------------- driver
def compute(data_dir="data", molecules=None, n_sphere=N_SPHERE):
    """The three tables as DataFrames (spread, on-axis pairs, displacement)."""
    R = roster(data_dir)
    if molecules is not None:
        R = R[R.molecule.isin(molecules)]
    spread, pairs, disp = [], [], []
    for r in R.itertuples():
        d = molecule_inputs(r.molecule, data_dir)
        spread += frame_spread_rows(d, r.mol_type, r.point_group, n_sphere)
        p = onaxis_pair_rows(d, r.mol_type, r.point_group)
        pairs += p
        disp += displacement_rows(d, p)
    return pd.DataFrame(spread), pd.DataFrame(pairs), pd.DataFrame(disp)


def summary(S, P, D):
    V = S[S.max_abs >= MIN_ABS]
    tops = V[V.case.isin(["symmetric", "spherical"])]
    out = dict(spread_tops_n=len(tops), spread_tops_median=float(tops.spread.median()),
               spread_tops_p90=float(tops.spread.quantile(0.9)), spread_tops_max=float(tops.spread.max()))
    for c in CASES:
        g = V[V.case == c]
        out[f"spread_{c}_n"] = len(g)
        out[f"spread_{c}_molecules"] = int(S[S.case == c].molecule.nunique())
        out[f"spread_{c}_median"] = float(g.spread.median()) if len(g) else np.nan
    for tag, Q in (("selected", P[P.selected]), ("allowed", P)):
        on, off = Q[Q.on_by_symmetry], Q[~Q.on_by_symmetry]
        out[f"gap_{tag}_on_pairs"] = len(on)
        out[f"gap_{tag}_on_exact_zero"] = int((on.r_Q <= ZERO_TOL).sum())
        out[f"gap_{tag}_off_pairs"] = len(off)
        out[f"gap_{tag}_largest_on"] = float(on.r_Q.max())
        out[f"gap_{tag}_smallest_off"] = float(off.r_Q.min())
        out[f"gap_{tag}_on_above_cutoff"] = int((on.r_Q > ON_AXIS_TOL).sum())
        out[f"gap_{tag}_off_within_cutoff"] = int((off.r_Q <= ON_AXIS_TOL).sum())
    out["gap_allowed_frames"] = int(P.groupby(["molecule", "frame_index"]).ngroups)
    out["disp_points"] = len(D)
    for reg, g in D.groupby("region"):
        out[f"disp_{reg.replace(' ', '_')}"] = len(g)
    mv = D[D.d > 0]
    out["disp_smallest_nonzero_d"] = float(mv.d.min())
    out["disp_largest_d_le_floor"] = float(D[D.d <= EPS_DENOM].d.max())
    return out


def run(data_dir="data", out_dir="data/figures", results_dir="data/results"):
    t0 = time.time()
    S, P, D = compute(data_dir)
    os.makedirs(results_dir, exist_ok=True)
    paths = {}
    for name, df in (("spread", S), ("onaxis_pairs", P), ("onaxis_displacement", D)):
        p = os.path.join(results_dir, f"si_frame_{name}.csv")
        df.to_csv(p, index=False, float_format="%.6g")
        paths[f"csv_{name}"] = p
    for key, (pdf, png) in (("fig_spread", _figure_spread(S, out_dir)),
                            ("fig_gap", _figure_gap(P, out_dir)),
                            ("fig_displacement", _figure_displacement(D, out_dir))):
        paths[f"{key}_pdf"], paths[f"{key}_png"] = pdf, png
    out = summary(S, P, D)
    out.update(paths)                     # flat path strings (reproduce.py collects them)
    out["seconds"] = round(time.time() - t0, 1)
    return out


if __name__ == "__main__":
    res = run()
    for k, v in res.items():
        print(f"{k}: {v}")
