"""SI table and figures for referee point R1.4: basis dependence of s[V_S] inside degenerate sets.

Ported from the R1.4 side report (JCC/2nd submission/side_reports/R1.4_degenerate_modes/
r14_degenerate_modes.py, report section 2) with the revision's settings: tau_SB from
Thresholds.calibrated() (0.43) and the S/B ambiguity range src.figures.AMBIGUITY_RANGE (0.31-0.53).

Degenerate set (definition fixed by the author): modes of one molecule that
  (i) carry the SAME degenerate irrep (E, T or Pi) of the molecule's point group in
      data/mol_list_method.csv (`point_group`), with as many members as the irrep's dimension, and
 (ii) have equal frequencies (chained |dnu| < 0.5 cm-1).
Gaussian prints the irreps (mode key `irrep`). Close modes with different (1-D) irreps are scored as
single modes. Gaussian's '?A' labels (AlCl3) are resolved by the trace of the C3 operation on the
pair's mass-weighted span: -1 (E-type) makes the pair a set.

Inside a set any orthonormal recombination of the members is an equally valid set of normal modes.
Rotations mix the MASS-WEIGHTED vectors (the physical metric, in which normal modes are orthonormal);
each recombined vector is converted back to Cartesian displacements and scored with eq:vscore,
exactly as ModeScorer.Vscore (side report section 2.4). With B_k the Cartesian back-transforms of a
mass-weighted orthonormal basis of the set, a single member is d(c) = sum_k c_k B_k with |c| = 1 and a
basis is an orthogonal k x k matrix O whose rows are the members' coefficients. Per set:
  Vmin, Vmax, spread   lowest / highest s[V_S] of a single member over ALL directions c
                       (pairs: grid over the half circle + bounded refinement; triples: Gaussian's
                       vectors + 20 000 Fibonacci directions + Nelder-Mead from the 20 best)
  V_mean (s-bar)       subspace mean: s[V_S] averaged uniformly over directions
  member_avg_orig      s-bar_G, the mean of Gaussian's own members (the manuscript's Fig. 7 practice)
  w                    within-basis range of ONE basis: highest - lowest of its members
  member mean          mean of the members of ONE basis
Both w and the member mean are sampled over `n_rot` rotations (bases) per set: pairs on the grid
theta = 0, 90/n_rot, ... < 90 deg (member 2 at theta + 90); triples `n_rot` Haar-random bases. Their
lowest/highest values are refined by optimisation over the rotation angles (triples: three angles,
side report section 2.5). s[V_S] does not depend on the scoring frame, so the frame chosen by
main.build_scorer_and_final (src/frame.py) does not change any number here.

    python -m src.si_degenerate                 # full run (80 s on an idle machine; the 30 triples dominate)

Outputs
  data/results/si_degenerate_sets.csv          one row per set (152)
  data/results/si_degenerate_rotations.csv     per set, histogram counts over the n_rot rotations of
                                               w and of (member mean - s-bar) (long format)
  data/figures/fig_si_degenerate_all_sets      every set: position on 0-1 and range about s-bar
  data/figures/fig_si_degenerate_near_threshold sets whose range meets the ambiguity range
  data/figures/fig_si_degenerate_mean_comparison s-bar vs s-bar_G; s-bar_G - s-bar vs spread
  data/figures/fig_si_degenerate_rotations     distributions of w and member mean - s-bar over rotations
  data/figures/fig_si_degenerate_spread        histogram of the spreads
"""

import os
import re
import time
import zlib

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.optimize import minimize, minimize_scalar

from src.classifier import Thresholds
from src.figures import _style, _savefig, COLORS, AMBIGUITY_RANGE

TOL_CM = 0.5                      # cm-1, chained frequency tolerance of a set
N_DIR_TRIPLE = 20000              # Fibonacci directions for the subspace mean / extremes of triples
DEFAULT_N_ROT = 1800              # rotations (bases) per set
DEFAULT_SEED = 20261105
DEG_POINT_GROUPS = {"C3v", "D3h", "Td", "Oh", "C4v", "D4h", "D6h", "Dinfh"}
GROUP_OF_TYPE = {"ideal": "SC", "non-ideal": "SC", "test": "test", "multi-centre": "benzene"}
GROUP_ORDER = ("SC", "test", "benzene")
GROUP_LABEL = {"SC": "SC set", "test": "test set", "benzene": "benzene"}
GROUP_MARKER = {"SC": "o", "test": "s", "benzene": "^"}
GROUP_LINESTYLE = {"SC": "-", "test": "--", "benzene": ":"}
GROUP_GREY = {"SC": "#4D4D4D", "test": "#9A9A9A", "benzene": "#DADADA"}
SYM_CLASSES = ("E3", "T", "E4", "Pi")
# Symmetry-class colours: Paul Tol's muted set, none of them used by the manuscript's figures
# (S red / B blue stay reserved for labels, grey for external modes and the threshold).
SYM_CLASS_COLOR = {"E3": "#44AA99", "T": "#999933", "E4": "#882255", "Pi": "#DDCC77"}
SYM_CLASS_LABEL = {"E3": "E pair, 3- or 6-fold axis", "T": "T triple",
                   "E4": "E pair, 4-fold axis", "Pi": r"$\Pi$ pair (linear)"}
SV = r"$s[\mathrm{V_S}]$"
SBAR = r"$\bar{s}$"
SBAR_G = r"$\bar{s}_\mathrm{G}$"
SPREAD_LABEL = r"spread of $s[\mathrm{V_S}]$ (highest $-$ lowest)"
SETS_CSV = "si_degenerate_sets.csv"
ROTATIONS_CSV = "si_degenerate_rotations.csv"
W_BINS = np.round(np.arange(0.0, 0.1501, 0.0025), 6)
MM_BINS = np.round(np.arange(-0.15, 0.1501, 0.0025), 6)


# ---------------------------------------------------------------- set detection
def irrep_dim(label):
    """Dimension of a Gaussian irrep label (E/Pi/Delta -> 2, T/F -> 3, A/B/Sigma -> 1; else None)."""
    if not isinstance(label, str):
        return None
    u = label.upper().strip("'\"")
    if u.startswith(("PI", "DLT", "DELTA", "PHI")):
        return 2
    if u.startswith("E"):
        return 2
    if u.startswith(("T", "F")):
        return 3
    if u.startswith(("A", "B", "SG")):
        return 1
    return None


def sym_class(point_group, k):
    """How the point group acts on a set's coefficient space (side report section 2.5).
    E3: E pair under a 3- or 6-fold axis (C3v, D3h, D6h, E of Td/Oh); E4: E pair under a 4-fold axis
    (C4v, D4h); Pi: linear molecule (s[V_S] constant); T: three-fold set."""
    if k == 3:
        return "T"
    if point_group == "Dinfh":
        return "Pi"
    if point_group in ("C4v", "D4h"):
        return "E4"
    return "E3"


def _c3_character(scorer, V):
    """Trace of the C3 operation on the mass-weighted span of V (k, N, 3); the C3 axis is the
    principal axis of the non-degenerate moment. E-type pair -> -1; two A-type modes -> +2."""
    mom, ax = scorer.principal_axes()
    d = np.abs(np.diff(np.sort(mom)))
    u = ax[:, 2] if d[0] < d[1] else ax[:, 0]
    ang = 2 * np.pi / 3
    K = np.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]])
    R = np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K
    X = scorer.coords @ R.T
    perm = np.array([int(np.argmin(np.linalg.norm(scorer.coords - x, axis=1))) for x in X])
    if not (np.linalg.norm(scorer.coords[perm] - X, axis=1).max() < 1e-3 and len(set(perm)) == len(perm)):
        raise ValueError("no C3 axis along the non-degenerate principal axis")
    m_sqrt = np.sqrt(scorer.masses)[:, None]
    k = len(V)
    Q = np.linalg.qr((V * m_sqrt[None]).reshape(k, -1).T)[0].T.reshape(k, -1, 3)
    RQ = np.zeros_like(Q)
    RQ[:, perm, :] = Q @ R.T                    # displacement of atom A moves to atom perm[A]
    return float(np.trace(np.einsum("inc,jnc->ij", Q, RQ)))


def find_degenerate_sets(data_dir="data", molecules=None, tol=TOL_CM, keep_modes=False, diagnostics=None):
    """Degenerate sets of the SC (ideal + non-ideal), test and benzene (multi-centre) molecules.

    Scores come from the program's pipeline: main.load_inputs(mol, "normal") (tracked intermediates)
    and main.build_scorer_and_final (frame rule of src/frame.py). Returns a list of dicts, one per set,
    in roster order and increasing frequency, with keys molecule, group, mol_type, point_group, freq
    (mean member frequency, 4 dp), deg, irreps ('/'-joined Gaussian labels), status ('symmetry' or
    'symmetry-unlabelled' for '?A' pairs resolved by the C3 character), sym_class, member_index
    (0-based indices into the molecule's vibrational modes), member_freqs, freq_spread. With
    keep_modes=True each dict also holds 'scorer' (the molecule's ModeScorer, shared by its sets) and
    'vectors' ((k, N, 3) Cartesian member vectors in the scorer's frame). If `diagnostics` is a dict it
    receives 'close_pairs_not_degenerate' and 'c3_character' (lists of dicts). Order: roster, then
    frequency clusters in increasing frequency, then irrep label inside a cluster (as the side report;
    only C10H16 3051 T2 / 3052 E share a cluster).
    """
    import main as M          # deferred: main imports src modules at import time

    roster = pd.read_csv(os.path.join(data_dir, "mol_list_method.csv"))
    roster = roster[roster["mol_type"].isin(list(GROUP_OF_TYPE))]
    if molecules is not None:
        roster = roster[roster["molecule"].isin(list(molecules))]
    excluded, unlabelled, out = [], [], []
    for _, r in roster.iterrows():
        mol, pg = r["molecule"], r["point_group"]
        raw, _ = M.load_inputs(mol, "normal", data_dir)
        scorer, final = M.build_scorer_and_final(raw, "normal")
        modes = final[len(final) - len(raw["modes"]):]          # drop the prepended ideal T/R references
        fr = np.array([m["frequency"] for m in modes], float)
        irr = [m.get("irrep") for m in modes]
        order = np.argsort(fr, kind="stable")
        clusters, cur = [], [order[0]]
        for a, b in zip(order[:-1], order[1:]):
            if fr[b] - fr[a] < tol:
                cur.append(b)
            else:
                clusters.append(cur)
                cur = [b]
        clusters.append(cur)
        sets = []
        for grp in clusters:
            for j in grp:
                if (irrep_dim(irr[j]) or 1) > 1 and pg not in DEG_POINT_GROUPS:
                    raise ValueError(f"{mol}: degenerate irrep {irr[j]} in point group {pg}")
            rest = list(grp)
            for lab in sorted(set(x for x in (irr[j] for j in grp) if isinstance(x, str))):
                mem = [j for j in grp if irr[j] == lab]
                if (irrep_dim(lab) or 1) > 1 and len(mem) == irrep_dim(lab):
                    sets.append((mem, "symmetry"))
                    rest = [j for j in rest if j not in mem]
            unl = [j for j in rest if not isinstance(irr[j], str) or irr[j].startswith("?")]
            if len(unl) == 2 and pg in DEG_POINT_GROUPS and np.ptp(fr[unl]) < 0.01:
                chi = _c3_character(scorer, np.array([modes[j]["vector"] for j in unl], float))
                unlabelled.append(dict(molecule=mol, freq=float(fr[unl].mean()), c3_character=chi))
                if abs(chi + 1) < 1e-3:
                    sets.append((unl, "symmetry-unlabelled"))
                    rest = [j for j in rest if j not in unl]
            if len(rest) > 1:
                excluded.append(dict(molecule=mol, group=GROUP_OF_TYPE[r["mol_type"]], point_group=pg,
                                     freqs="/".join(f"{fr[j]:.3f}" for j in rest),
                                     irreps="/".join(str(irr[j]) for j in rest),
                                     freq_spread=float(np.ptp(fr[rest]))))
        # order as the side report: clusters by frequency, inside a cluster by irrep label
        for mem, status in sets:
            k = len(mem)
            d = dict(molecule=mol, group=GROUP_OF_TYPE[r["mol_type"]], mol_type=r["mol_type"], point_group=pg,
                     freq=round(float(fr[mem].mean()), 4), deg=k, irreps="/".join(str(irr[j]) for j in mem),
                     status=status, sym_class=sym_class(pg, k), member_index=[int(j) for j in mem],
                     member_freqs=[float(fr[j]) for j in mem], freq_spread=float(np.ptp(fr[mem])))
            if keep_modes:
                d["scorer"] = scorer
                d["vectors"] = np.array([modes[j]["vector"] for j in mem], float)
            out.append(d)
    if diagnostics is not None:
        diagnostics["close_pairs_not_degenerate"] = excluded
        diagnostics["c3_character"] = unlabelled
    return out


# ---------------------------------------------------------------- scoring inside a set
class VScan:
    """Vectorised s[V_S] (eq:vscore, with the scorer's bond weights, as ModeScorer.Vscore) of many
    members d(c) = sum_k c_k B_k of one set. B: (k, N, 3) Cartesian basis vectors."""

    def __init__(self, scorer, B):
        i = np.array([b[0] for b in scorer.bList])
        j = np.array([b[1] for b in scorer.bList])
        bvec = scorer.coords[j] - scorer.coords[i]
        self.bhat = bvec / np.linalg.norm(bvec, axis=1)[:, None]
        self.w = np.asarray(scorer.bond_weights, float)
        self.db = B[:, j, :] - B[:, i, :]                 # (k, n_bonds, 3)

    def __call__(self, C):
        """C: (M, k) coefficient rows -> (M,) s[V_S]."""
        db = np.einsum("mk,kbj->mbj", np.atleast_2d(C), self.db)
        L = np.linalg.norm(db, axis=2)
        P = np.abs(np.einsum("mbj,bj->mb", db, self.bhat))
        num = (self.w * L * P).sum(1)
        den = (self.w * L * L).sum(1)
        return np.where(den > 1e-30, num / np.maximum(den, 1e-300), 0.0)


def mass_weighted_basis(scorer, V):
    """Cartesian back-transforms B_k = Q_k / sqrt(m) of a mass-weighted orthonormal basis Q of span(V)."""
    k = len(V)
    m_sqrt = np.sqrt(scorer.masses)[:, None]
    Q = np.linalg.qr((V * m_sqrt[None]).reshape(k, -1).T)[0].T.reshape(k, -1, 3)
    return Q / m_sqrt[None]


def frame3(a, b, g):
    """Orthonormal coefficient rows of a three-fold basis from three angles (deg), side report 2.5:
    c1 = r-hat(a, b); c2 = cos g theta-hat + sin g phi-hat; c3 = c1 x c2."""
    a, b, g = np.radians([a, b, g])
    r = np.array([np.sin(a) * np.cos(b), np.sin(a) * np.sin(b), np.cos(a)])
    th = np.array([np.cos(a) * np.cos(b), np.cos(a) * np.sin(b), -np.sin(a)])
    ph = np.array([-np.sin(b), np.cos(b), 0.0])
    c2 = np.cos(g) * th + np.sin(g) * ph
    return np.vstack([r, c2, np.cross(r, c2)])


def fib_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]


def haar(k, rng):
    """Haar-random orthogonal k x k matrix (QR of a Gaussian matrix, sign fix of Mezzadri 2007)."""
    Q, R = np.linalg.qr(rng.normal(size=(k, k)))
    return Q * np.sign(np.diag(R))


def _angle_extremes3(fun, rng, n_start, n_refine):
    """Lowest and highest of fun(three angles) over all bases of a triple (Nelder-Mead from the best of
    n_start random angle triples)."""
    S0 = rng.uniform(0, 360, size=(n_start, 3))
    v0 = np.array([fun(x) for x in S0])
    out = []
    for sgn in (1, -1):
        best = (sgn * v0).min()
        for x0 in S0[np.argsort(sgn * v0)[:n_refine]]:
            res = minimize(lambda x: sgn * fun(x), x0, method="Nelder-Mead",
                           options=dict(xatol=1e-6, fatol=1e-10, maxiter=3000))
            best = min(best, res.fun)
        out.append(sgn * best)
    return out[0], out[1]


def _direction_extremes(f, k, C0, n_rot):
    """True min / max of f over unit vectors in R^k, from starting coefficients C0."""
    vals = f(C0)
    out = []
    if k == 2:
        th0 = np.arctan2(C0[:, 1], C0[:, 0])
        step = np.pi / n_rot                        # two grid spacings
        for sgn in (1, -1):
            best = (sgn * vals).min()
            for t in th0[np.argsort(sgn * vals)[:5]]:
                res = minimize_scalar(lambda t_: sgn * f(np.array([[np.cos(t_), np.sin(t_)]]))[0],
                                      bounds=(t - step, t + step), method="bounded", options=dict(xatol=1e-10))
                best = min(best, res.fun)
            out.append(sgn * best)
        return out[0], out[1]
    for sgn in (1, -1):
        best = (sgn * vals).min()
        for s0 in C0[np.argsort(sgn * vals)[:20]]:
            res = minimize(lambda c: sgn * f(c / np.linalg.norm(c))[0], s0, method="Nelder-Mead",
                           options=dict(xatol=1e-8, fatol=1e-11, maxiter=4000))
            best = min(best, res.fun)
        out.append(sgn * best)
    return out[0], out[1]


def _set_rng(seed, s):
    """Per-set random stream: a reduced run over a few sets reproduces the full run's numbers."""
    return np.random.default_rng([int(seed), zlib.crc32(f"{s['molecule']} {s['freq']:.4f}".encode())])


def analyse_set(s, tau, n_rot=DEFAULT_N_ROT, seed=DEFAULT_SEED, n_start=200, n_refine=12, n_dir=N_DIR_TRIPLE):
    """Scan one set (a find_degenerate_sets(..., keep_modes=True) entry). Returns (row, w, mm_dev):
    the CSV row, and the within-basis range and (member mean - s-bar) of each of the n_rot rotations."""
    sc, V, k = s["scorer"], s["vectors"], s["deg"]
    rng = _set_rng(seed, s)
    orig = np.array([sc.calculate_scores(v)["V"] for v in V])
    B = mass_weighted_basis(sc, V)
    f = VScan(sc, B)
    Corig = np.linalg.lstsq(B.reshape(k, -1).T, V.reshape(k, -1).T, rcond=None)[0].T
    Corig /= np.linalg.norm(Corig, axis=1)[:, None]
    if k == 2:
        th = np.linspace(0, np.pi, 2 * n_rot, endpoint=False)
        C = np.c_[np.cos(th), np.sin(th)]
        vals = f(C)
        mean, mean_se = float(vals.mean()), 0.0
        s1, s2 = vals[:n_rot], vals[n_rot:]                 # member 1 at theta, member 2 at theta + 90
        w = np.abs(s1 - s2)
        mm = 0.5 * (s1 + s2)
        gw = lambda t: -abs(f(np.array([[np.cos(t), np.sin(t)]]))[0] - f(np.array([[-np.sin(t), np.cos(t)]]))[0])
        w_hi = w.max()
        for t in th[:n_rot][np.argsort(-w)[:5]]:
            res = minimize_scalar(gw, bounds=(t - np.pi / n_rot, t + np.pi / n_rot), method="bounded",
                                  options=dict(xatol=1e-10))
            w_hi = max(w_hi, -res.fun)
        w_lo = float(w.min())
        mixed = (np.minimum(s1, s2) < tau) & (np.maximum(s1, s2) >= tau)
        mm_lo, mm_hi = float(mm.min()), float(mm.max())
        C0 = np.vstack([Corig, C])
    else:
        P = fib_sphere(n_dir)
        vals = f(P)
        mean, mean_se = float(vals.mean()), float(vals.std() / np.sqrt(len(vals)))
        SB = np.array([f(haar(3, rng)) for _ in range(n_rot)])     # (n_rot, 3) members of each basis
        w = np.ptp(SB, axis=1)
        mm = SB.mean(1)
        mixed = (SB.min(1) < tau) & (SB.max(1) >= tau)
        w_lo, w_hi = _angle_extremes3(lambda x: np.ptp(f(frame3(*x))), rng, n_start, n_refine)
        w_lo, w_hi = min(w_lo, w.min()), max(w_hi, w.max())
        mm_lo, mm_hi = _angle_extremes3(lambda x: f(frame3(*x)).mean(), rng, n_start, n_refine)
        mm_lo, mm_hi = min(mm_lo, mm.min(), orig.mean()), max(mm_hi, mm.max(), orig.mean())
        C0 = np.vstack([Corig, P])
    vmin, vmax = _direction_extremes(f, k, C0, n_rot)
    vmin, vmax = min(vmin, orig.min()), max(vmax, orig.max())
    lab = lambda v: "S" if v >= tau else "B"
    row = {key: s[key] for key in ("molecule", "group", "mol_type", "point_group", "freq", "deg", "irreps",
                                   "status", "sym_class")}
    row.update(
        member_freqs="/".join(f"{x:.4f}" for x in s["member_freqs"]), freq_spread=s["freq_spread"],
        V_orig="/".join(f"{v:.6f}" for v in orig), V_orig_min=float(orig.min()), V_orig_max=float(orig.max()),
        Vmin=float(vmin), Vmax=float(vmax), spread=float(vmax - vmin),
        V_mean=mean, V_mean_se=mean_se, V_sd=float(vals.std()),
        member_avg_orig=float(orig.mean()), member_minus_mean=float(orig.mean() - mean),
        member_avg_min=float(mm_lo), member_avg_max=float(mm_hi), mm_range=float(mm_hi - mm_lo),
        mm_p25=float(np.quantile(mm, 0.25)), mm_median=float(np.median(mm)), mm_p75=float(np.quantile(mm, 0.75)),
        mm_sd=float(mm.std()), mm_frac_S=float((mm >= tau).mean()),
        w_gauss=float(np.ptp(orig)), w_min=float(w_lo), w_p25=float(np.quantile(w, 0.25)),
        w_median=float(np.median(w)), w_p75=float(np.quantile(w, 0.75)), w_p90=float(np.quantile(w, 0.9)),
        w_mean=float(w.mean()), w_max=float(w_hi), w_frac_mixed_labels=float(mixed.mean()),
        tau_SB=float(tau), labels_orig="/".join(lab(v) for v in orig), label_mean=lab(mean),
        label_member_mean=lab(orig.mean()),
        orig_mixed_labels=bool((orig.min() < tau) and (orig.max() >= tau)),
        straddles_tau=bool(vmin < tau <= vmax),
        meets_ambiguity_range=bool(vmax >= AMBIGUITY_RANGE[0] and vmin <= AMBIGUITY_RANGE[1]),
        n_rot=int(n_rot), vscan_dev=float(np.abs(f(Corig) - orig).max()),
    )
    return row, w, mm - mean


def rotation_histograms(D, samples):
    """Long-format histogram counts per set over its rotations: quantity in {w, mm_minus_mean}."""
    rows = []
    for (_, r), (w, dm) in zip(D.iterrows(), samples):
        for q, x, bins in (("w", w, W_BINS), ("mm_minus_mean", dm, MM_BINS)):
            h, e = np.histogram(np.clip(x, bins[0], bins[-1]), bins=bins)
            for lo, hi, c in zip(e[:-1], e[1:], h):
                if c:
                    rows.append(dict(molecule=r.molecule, freq=r.freq, group=r.group, sym_class=r.sym_class,
                                     quantity=q, bin_lo=lo, bin_hi=hi, count=int(c)))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- figures
def _chem(mol):
    """Formula with subscript digits (mathtext)."""
    return re.sub(r"(\d+)", r"$_{\1}$", mol)


def _irrep_tex(irreps):
    """Gaussian's irrep string -> mathtext, e.g. PIU -> Pi_u, E1U -> E_1u, T2G -> T_2g; '?A' (AlCl3,
    D3h, resolved by the C3 character) -> E'."""
    lab = irreps.split("/")[0].replace("?A", "E'")
    prime = "".join(ch for ch in lab if ch in "'\"").replace('"', "''")
    core = "".join(ch for ch in lab if ch not in "'\"")
    head, sub = (r"\Pi", core[2:]) if core.upper().startswith("PI") else (core[0], core[1:])
    sub = sub.lower()
    return f"${head}{prime}" + (f"_{{\\mathrm{{{sub}}}}}$" if sub else "$")


def _members(row):
    return np.array([float(x) for x in str(row.V_orig).split("/")])


def _threshold_axes(ax, tau, vertical=True):
    lo, hi = AMBIGUITY_RANGE
    if vertical:
        ax.axvspan(lo, hi, color=COLORS["background"], alpha=0.35, lw=0, zorder=0)
        ax.axvline(tau, color=COLORS["threshold"], ls="--", lw=0.8, zorder=1)
    else:
        ax.axhspan(lo, hi, color=COLORS["background"], alpha=0.35, lw=0, zorder=0)
        ax.axhline(tau, color=COLORS["threshold"], ls="--", lw=0.8, zorder=1)


def _class_handles(lw=4):
    return [Line2D([], [], color=SYM_CLASS_COLOR[c], lw=lw, solid_capstyle="butt", label=SYM_CLASS_LABEL[c])
            for c in SYM_CLASSES]


def plot_all_sets(D, tau, out_dir, label="fig_si_degenerate_all_sets"):
    """Every set (replaces slide 8): per column, (left) the range of s[V_S] over all bases on the 0-1
    scale with a tick at s-bar; (right) the same range about s-bar, Gaussian's members (circles) and their
    mean s-bar_G (cross). Column 1: SC set; column 2: test set + benzene."""
    _style()
    cols = [("SC set", D[D.group == "SC"].sort_values("V_mean")),
            ("test set + benzene", pd.concat([D[D.group == "test"].sort_values("V_mean"),
                                              D[D.group == "benzene"].sort_values("V_mean")]))]
    nmax = max(len(x) for _, x in cols)
    dev = np.r_[D.Vmin - D.V_mean, D.Vmax - D.V_mean]
    lim = np.ceil(np.abs(dev).max() * 100 + 0.6) / 100
    fs_lab, fs_tick = 5.6, 6.5
    fig = plt.figure(figsize=(7.0, 9.6))
    top, bot = 0.965, 0.085
    row_h = (top - bot) / nmax
    for ci, (title, x) in enumerate(cols):
        x0 = 0.115 + ci * 0.5
        h = row_h * len(x)
        axA = fig.add_axes([x0, top - h, 0.165, h])
        axB = fig.add_axes([x0 + 0.185, top - h, 0.165, h])
        _threshold_axes(axA, tau)
        axB.axvline(0, color="k", lw=0.5, zorder=1)
        labels = []
        for y, (_, r) in enumerate(x.iterrows()):
            c = SYM_CLASS_COLOR[r.sym_class]
            axA.plot([r.Vmin, r.Vmax], [y, y], color=c, lw=2.6, solid_capstyle="butt", zorder=2)
            axA.plot(r.V_mean, y, "|", color="k", ms=4.5, mew=0.9, zorder=3)
            axB.plot([r.Vmin - r.V_mean, r.Vmax - r.V_mean], [y, y], color=c, lw=2.6, solid_capstyle="butt", zorder=2)
            axB.plot(_members(r) - r.V_mean, np.full(r.deg, y), "o", ms=2.5, mfc="white", mec="k", mew=0.5, zorder=3)
            axB.plot(r.member_avg_orig - r.V_mean, y, "x", ms=3.0, color="k", mew=0.7, zorder=4)
            star = "* " if r.straddles_tau else ""
            labels.append(f"{star}{_chem(r.molecule)} {r.freq:.0f} {_irrep_tex(r.irreps)}")
        nt = int((x.group == "test").sum())
        if 0 < nt < len(x):
            for a in (axA, axB):
                a.axhline(nt - 0.5, color="k", lw=0.6)
            for name, y0 in (("test set", -0.3), ("benzene", nt - 0.3)):
                axB.text(0.97, y0, name, transform=axB.get_yaxis_transform(), ha="right", va="top",
                         fontsize=fs_lab + 0.4, style="italic")
        for a in (axA, axB):
            a.set_ylim(len(x) - 0.5, -0.5)
            a.tick_params(axis="x", labelsize=fs_tick, length=2.5, pad=1.5)
            a.spines["top"].set_visible(True)
            a.spines["right"].set_visible(True)
            a.spines["top"].set_linewidth(0.6)
            a.spines["right"].set_linewidth(0.6)
        axA.set_yticks(range(len(x)))
        axA.set_yticklabels(labels, fontsize=fs_lab)
        axA.tick_params(axis="y", length=0, pad=2)
        axB.set_yticks([])
        axA.set_xlim(0, 1)
        axA.set_xticks([0, tau, 1])
        axA.set_xticklabels(["0", f"{tau:.2f}", "1"])
        axB.set_xlim(-lim, lim)
        major = 0.1 if lim >= 0.1 else 0.05
        tk = np.arange(-np.floor(lim / major) * major, lim + 1e-9, major)
        axB.set_xticks(tk)
        axB.set_xticklabels(["0" if abs(t) < 1e-9 else f"{t:.2f}".rstrip("0").replace("-", "−") for t in tk])
        axB.set_xticks(np.arange(-np.floor(lim / 0.05) * 0.05, lim + 1e-9, 0.05), minor=True)
        axB.tick_params(axis="x", which="minor", length=1.5)
        axA.set_xlabel(SV, fontsize=7.5, labelpad=1.5)
        axB.set_xlabel(SV + r"$\;-\;\bar{s}$", fontsize=7.5, labelpad=1.5)
        axA.set_title(f"({'ab'[ci]}) {title}", fontsize=8, loc="left", x=-0.62, pad=4, fontweight="bold")
    hand = _class_handles(lw=2.6)
    hand += [Line2D([], [], marker="|", ls="", color="k", ms=6, mew=1.0, label=r"subspace mean $\bar{s}$"),
             Line2D([], [], marker="o", ls="", mfc="white", mec="k", ms=3.5, mew=0.6, label="Gaussian's members"),
             Line2D([], [], marker="x", ls="", color="k", ms=4, mew=0.8, label=r"their mean $\bar{s}_\mathrm{G}$"),
             Patch(facecolor=COLORS["background"], alpha=0.35, lw=0, label="ambiguity range (0.31–0.53)"),
             Line2D([], [], color=COLORS["threshold"], ls="--", lw=0.8, label=r"$\tau_\mathrm{SB}$" + f" = {tau:.2f}")]
    fig.legend(handles=hand, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=6.8, frameon=False,
               columnspacing=1.4, handlelength=1.8)
    paths = _savefig(fig, out_dir, label)
    plt.close(fig)
    return paths


def near_threshold_sets(D, min_spread=0.01):
    """Sets whose range of s[V_S] over all bases meets the ambiguity range, with spread > min_spread."""
    lo, hi = AMBIGUITY_RANGE
    return D[(D.Vmax >= lo) & (D.Vmin <= hi) & (D.spread > min_spread)]


def plot_near_threshold(D, tau, out_dir, label="fig_si_degenerate_near_threshold"):
    """Sets near the threshold (replaces slide 9): range over all bases (bar, symmetry-class colour),
    Gaussian's members (circles filled by their label at tau_SB: S red, B blue) and s-bar (tick)."""
    _style()
    show = near_threshold_sets(D).copy()
    show["g"] = show.group.map({g: i for i, g in enumerate(GROUP_ORDER)})
    show = show.sort_values(["g", "V_mean"])
    n = len(show)
    fig, ax = plt.subplots(figsize=(4.6, 0.95 + 0.215 * n))
    _threshold_axes(ax, tau)
    yl = []
    for y, (_, r) in enumerate(show.iterrows()):
        ax.plot([r.Vmin, r.Vmax], [y, y], color=SYM_CLASS_COLOR[r.sym_class], lw=5, solid_capstyle="butt", zorder=2)
        ax.plot(r.V_mean, y, "|", color="k", ms=11, mew=1.8, zorder=3)
        for v in _members(r):
            c = COLORS["stretching"] if v >= tau else COLORS["bending"]
            ax.plot(v, y, "o", ms=5, mfc=c, mec="k", mew=0.6, zorder=4)
        yl.append(f"{'* ' if r.straddles_tau else ''}{_chem(r.molecule)} {r.freq:.0f} {_irrep_tex(r.irreps)}")
    bounds = np.flatnonzero(np.diff(show.g.to_numpy())) + 0.5
    for b in bounds:
        ax.axhline(b, color="k", lw=0.6)
    starts = np.r_[0, np.ceil(bounds)].astype(int)
    xmin = min(show.Vmin.min(), AMBIGUITY_RANGE[0]) - 0.02
    xmax = max(show.Vmax.max(), AMBIGUITY_RANGE[1]) + 0.02
    for st in starts:
        ax.text(xmax - 0.004, st - 0.25, GROUP_LABEL[show.group.iloc[st]], ha="right", va="top", fontsize=8,
                style="italic")
    ax.set_yticks(range(n))
    ax.set_yticklabels(yl, fontsize=8.5)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel(SV)
    present = [c for c in SYM_CLASSES if c in set(show.sym_class)]
    hand = [Line2D([], [], color=SYM_CLASS_COLOR[c], lw=5, solid_capstyle="butt", label=f"range, {SYM_CLASS_LABEL[c]}")
            for c in present]
    hand += [Line2D([], [], marker="|", ls="", color="k", ms=9, mew=1.6, label=r"subspace mean $\bar{s}$"),
             Patch(facecolor=COLORS["background"], alpha=0.35, lw=0, label="ambiguity range (0.31–0.53)")]
    hand2 = [Line2D([], [], marker="o", ls="", mfc=COLORS["stretching"], mec="k", mew=0.6, ms=5,
                    label="Gaussian's member, S"),
             Line2D([], [], marker="o", ls="", mfc=COLORS["bending"], mec="k", mew=0.6, ms=5,
                    label="Gaussian's member, B"),
             Line2D([], [], color=COLORS["threshold"], ls="--", lw=0.8, label=r"$\tau_\mathrm{SB}$" + f" = {tau:.2f}")]
    hand2 += [Line2D([], [], ls="", label=" ")] * max(0, len(hand) - len(hand2))
    ax.legend(handles=hand + hand2, loc="lower center", bbox_to_anchor=(0.42, 1.0), ncol=2, fontsize=7.5,
              frameon=False, columnspacing=1.2, handlelength=1.8)
    paths = _savefig(fig, out_dir, label)
    plt.close(fig)
    return paths


def plot_mean_comparison(D, tau, out_dir, label="fig_si_degenerate_mean_comparison", annotate=None):
    """(a) s-bar vs s-bar_G (replaces slide 29a); grey quadrants: the two would give different labels at
    tau_SB. (b) s-bar_G - s-bar vs spread (slide 29b); grey wedge |difference| <= spread."""
    _style()
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.3), gridspec_kw=dict(wspace=0.36))
    ax = axs[0]
    ax.add_patch(plt.Rectangle((0, tau), tau, 1 - tau, color=COLORS["background"], alpha=0.35, lw=0, zorder=0))
    ax.add_patch(plt.Rectangle((tau, 0), 1 - tau, tau, color=COLORS["background"], alpha=0.35, lw=0, zorder=0))
    ax.plot([0, 1], [0, 1], color="k", lw=0.6, zorder=1)
    ax.axvline(tau, color=COLORS["threshold"], ls="--", lw=0.8, zorder=1)
    ax.axhline(tau, color=COLORS["threshold"], ls="--", lw=0.8, zorder=1)
    for g in GROUP_ORDER:
        for c in SYM_CLASSES:
            d = D[(D.group == g) & (D.sym_class == c)]
            ax.plot(d.member_avg_orig, d.V_mean, GROUP_MARKER[g], ms=4, mfc=SYM_CLASS_COLOR[c], mec="k", mew=0.3,
                    alpha=0.95, zorder=3)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xlabel(r"mean of Gaussian's members $\bar{s}_\mathrm{G}$")
    ax.set_ylabel(r"subspace mean $\bar{s}$")
    ax = axs[1]
    diff = D.member_minus_mean
    xmax = np.ceil((D.spread.max() + 0.005) * 100) / 100
    ylo = np.floor((diff.min() - 0.006) * 100) / 100
    yhi = np.ceil((diff.max() + 0.006) * 100) / 100
    sp = np.array([0, xmax])
    ax.fill_between(sp, -sp, sp, color=COLORS["background"], alpha=0.35, lw=0, zorder=0)
    ax.axhline(0, color="k", lw=0.6, zorder=1)
    for g in GROUP_ORDER:
        for c in SYM_CLASSES:
            d = D[(D.group == g) & (D.sym_class == c)]
            ax.plot(d.spread, d.member_minus_mean, GROUP_MARKER[g], ms=4, mfc=SYM_CLASS_COLOR[c], mec="k", mew=0.3,
                    alpha=0.95, zorder=3)
    if annotate is None:
        annotate = D.loc[diff.abs().sort_values(ascending=False).index[:4]]
    for _, r in annotate.iterrows():
        ax.annotate(f"{_chem(r.molecule)} {r.freq:.0f}", (r.spread, r.member_minus_mean), xytext=(-5, 0),
                    textcoords="offset points", ha="right", va="center", fontsize=7.5)
    ax.set_xlim(0, xmax)
    ax.set_ylim(ylo, yhi)
    ax.set_box_aspect(1)
    ax.set_xlabel(SPREAD_LABEL)
    ax.set_ylabel(r"$\bar{s}_\mathrm{G} - \bar{s}$")
    for a, t in zip(axs, "ab"):
        a.text(-0.3, 1.03, f"({t})", transform=a.transAxes, fontweight="bold", va="bottom")
        a.xaxis.label.set_size(10.5)
        a.yaxis.label.set_size(10.5)
    hand = [Line2D([], [], marker="s", ls="", mfc=SYM_CLASS_COLOR[c], mec="k", mew=0.3, ms=5.5,
                   label=SYM_CLASS_LABEL[c]) for c in SYM_CLASSES]
    hand += [Line2D([], [], marker=GROUP_MARKER[g], ls="", mfc="white", mec="k", mew=0.6, ms=5, label=GROUP_LABEL[g])
             for g in GROUP_ORDER]
    fig.legend(handles=hand, loc="lower center", bbox_to_anchor=(0.5, 0.93), ncol=4, fontsize=8, frameon=False,
               columnspacing=1.2, handletextpad=0.3)
    paths = _savefig(fig, out_dir, label)
    plt.close(fig)
    return paths


def plot_rotations(D, H, out_dir, label="fig_si_degenerate_rotations"):
    """Distributions over rotations (replaces slide 30), from the histogram table H: (a) within-basis range
    w by data set; (b) w by symmetry class (E 4-fold and Pi sets: w = 0 in every basis, noted);
    (c) member mean - s-bar by data set. Counts on a log scale."""
    _style()
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.7), gridspec_kw=dict(wspace=0.36))
    nset = lambda x: x.groupby(["molecule", "freq"]).ngroups

    def stairs(ax, h, bins, q, **kw):
        cnt = h[h.quantity == q].groupby("bin_lo")["count"].sum()
        y = np.array([cnt.get(b, 0) for b in bins[:-1]], float)
        ax.stairs(np.where(y > 0, y, np.nan), bins, **kw)
        return y

    wx = H[H.quantity == "w"]
    wmax = np.ceil((wx.bin_hi.max()) * 100) / 100
    mx = H[H.quantity == "mm_minus_mean"]
    mlo, mhi = np.floor(mx.bin_lo.min() * 50) / 50, np.ceil(mx.bin_hi.max() * 50) / 50
    for g in GROUP_ORDER:
        h = H[H.group == g]
        lab = f"{GROUP_LABEL[g]} ({int((D.group == g).sum())})"
        stairs(axs[0], h, W_BINS, "w", color="k", ls=GROUP_LINESTYLE[g], lw=1.1, label=lab)
        stairs(axs[2], h, MM_BINS, "mm_minus_mean", color="k", ls=GROUP_LINESTYLE[g], lw=1.1, label=lab)
    for c in ("E3", "T"):
        h = H[H.sym_class == c]
        stairs(axs[1], h, W_BINS, "w", color=SYM_CLASS_COLOR[c], lw=1.3,
               label=f"{SYM_CLASS_LABEL[c]} ({int((D.sym_class == c).sum())})")
    n0 = {c: int((D.sym_class == c).sum()) for c in ("E4", "Pi")}
    axs[1].text(0.05, 0.04, f"E pair, 4-fold axis ({n0['E4']}) and\n" + rf"$\Pi$ pair ({n0['Pi']}): $w = 0$ in every basis",
                transform=axs[1].transAxes, ha="left", va="bottom", fontsize=6.5)
    for ax in axs:
        ax.set_yscale("log")
        ax.set_ylim(0.7, 3e6)
        ax.set_ylabel("number of rotations", fontsize=9)
        ax.tick_params(labelsize=8)
    for ax in axs[:2]:
        ax.set_xlim(0, wmax)
        ax.set_xlabel(r"within-basis range $w$")
    axs[2].set_xlim(mlo, mhi)
    axs[2].xaxis.set_major_locator(plt.MultipleLocator(0.05))
    axs[2].axvline(0, color="k", lw=0.5)
    axs[2].set_xlabel(r"member mean $-\;\bar{s}$")
    for ax in axs:
        ax.xaxis.label.set_size(9)
    for ax, loc in zip(axs, ("upper right", "upper right", "upper left")):
        ax.legend(fontsize=6.5, frameon=False, loc=loc, handlelength=2.0, borderaxespad=0.2)
    for a, t in zip(axs, "abc"):
        a.text(-0.28, 1.02, f"({t})", transform=a.transAxes, fontweight="bold", va="bottom")
    paths = _savefig(fig, out_dir, label)
    plt.close(fig)
    return paths


def plot_spread(D, out_dir, label="fig_si_degenerate_spread"):
    """Histogram of the spreads of all sets, stacked by data set; median (solid) and 90th percentile
    (dashed) marked (replaces the histogram on slide 8)."""
    _style()
    fig, ax = plt.subplots(figsize=(3.4, 2.7))
    xmax = np.ceil((D.spread.max() + 0.005) * 100) / 100
    bins = np.arange(0, xmax + 1e-9, 0.005)
    ax.hist([D[D.group == g].spread for g in GROUP_ORDER], bins=bins, stacked=True,
            color=[GROUP_GREY[g] for g in GROUP_ORDER], edgecolor="k", lw=0.3,
            label=[f"{GROUP_LABEL[g]} ({int((D.group == g).sum())})" for g in GROUP_ORDER])
    ymax = ax.get_ylim()[1]
    for q, ls, yf in ((0.5, "-", 0.97), (0.9, "--", 0.55)):
        v = D.spread.quantile(q)
        ax.axvline(v, color="k", lw=0.8, ls=ls)
        ax.text(v + 0.002, ymax * yf, f"{'median' if q == 0.5 else '90th percentile'}\n{v:.3f}", fontsize=7.5,
                va="top")
    ax.set_xlim(0, xmax)
    ax.set_xlabel(SPREAD_LABEL)
    ax.set_ylabel("number of sets")
    ax.legend(fontsize=7.5, frameon=False, loc="upper right")
    paths = _savefig(fig, out_dir, label)
    plt.close(fig)
    return paths


# ---------------------------------------------------------------- driver
def headline_numbers(D, samples, tau):
    """Numbers quoted in the SI text."""
    w_all = np.concatenate([s[0] for s in samples]) if samples else np.array([])
    dm_all = np.concatenate([s[1] for s in samples]) if samples else np.array([])
    imax = (D.member_minus_mean.abs()).idxmax()
    smax = D.spread.idxmax()
    lab = lambda v: v >= tau
    return dict(
        tau_SB=float(tau), ambiguity_range=AMBIGUITY_RANGE,
        n_sets=len(D), n_sets_by_group={g: int((D.group == g).sum()) for g in GROUP_ORDER},
        n_sets_by_class={c: int((D.sym_class == c).sum()) for c in SYM_CLASSES},
        n_2fold=int((D.deg == 2).sum()), n_3fold=int((D.deg == 3).sum()),
        spread_median=float(D.spread.median()), spread_p90=float(D.spread.quantile(0.9)),
        spread_max=float(D.spread.max()), spread_max_set=f"{D.molecule[smax]} {D.freq[smax]:.0f}",
        straddling_sets=[f"{m} {f:.0f}" for m, f in zip(D.molecule[D.straddles_tau], D.freq[D.straddles_tau])],
        sets_with_mixed_gaussian_labels=[f"{m} {f:.0f}" for m, f in
                                         zip(D.molecule[D.orig_mixed_labels], D.freq[D.orig_mixed_labels])],
        n_sets_meeting_ambiguity_range=int(D.meets_ambiguity_range.sum()),
        n_near_threshold_shown=len(near_threshold_sets(D)),
        r_member_mean_vs_subspace_mean=float(np.corrcoef(D.member_avg_orig, D.V_mean)[0, 1]),
        max_abs_member_minus_mean=float(D.member_minus_mean.abs().max()),
        max_abs_member_minus_mean_set=f"{D.molecule[imax]} {D.freq[imax]:.0f}",
        n_label_disagree_member_mean_vs_subspace_mean=int((lab(D.member_avg_orig) != lab(D.V_mean)).sum()),
        sets_member_mean_label_depends_on_basis=[f"{m} {f:.0f}" for m, f, a, b in
                                                 zip(D.molecule, D.freq, D.member_avg_min, D.member_avg_max)
                                                 if a < tau <= b],
        sets_mixed_labels_in_some_basis=[f"{m} {f:.0f}" for m, f, x in zip(D.molecule, D.freq, D.w_frac_mixed_labels)
                                         if x > 0],
        n_rotations=int(len(w_all)),
        rot_share_w_lt_0p01=float((w_all < 0.01).mean()) if len(w_all) else np.nan,
        rot_w_p90=float(np.quantile(w_all, 0.9)) if len(w_all) else np.nan,
        rot_w_max=float(w_all.max()) if len(w_all) else np.nan,
        rot_share_absmm_lt_0p01=float((np.abs(dm_all) < 0.01).mean()) if len(dm_all) else np.nan,
        rot_absmm_p99=float(np.quantile(np.abs(dm_all), 0.99)) if len(dm_all) else np.nan,
        rot_absmm_max=float(np.abs(dm_all).max()) if len(dm_all) else np.nan,
    )


def scan_sets(sets, tau, n_rot=DEFAULT_N_ROT, seed=DEFAULT_SEED, n_start=200, n_refine=12, n_dir=N_DIR_TRIPLE):
    """analyse_set over a list of sets -> (DataFrame of rows, list of (w, mm_dev) samples)."""
    rows, samples = [], []
    for s in sets:
        row, w, dm = analyse_set(s, tau, n_rot=n_rot, seed=seed, n_start=n_start, n_refine=n_refine, n_dir=n_dir)
        rows.append(row)
        samples.append((w, dm))
    return pd.DataFrame(rows), samples


def run(data_dir="data", out_dir="data/figures", results_dir="data/results", n_rot=DEFAULT_N_ROT, seed=DEFAULT_SEED,
        n_start=200, n_refine=12, molecules=None, figures=True):
    """Detect the sets, scan them, write the two CSVs and the five figures. Returns a dict with the
    output paths, the headline numbers and the runtime (s)."""
    t0 = time.time()
    tau = Thresholds.calibrated().tau_SB
    diag = {}
    sets = find_degenerate_sets(data_dir, molecules=molecules, keep_modes=True, diagnostics=diag)
    D, samples = scan_sets(sets, tau, n_rot=n_rot, seed=seed, n_start=n_start, n_refine=n_refine)
    if D["vscan_dev"].max() > 1e-10:
        raise RuntimeError(f"vectorised s[V_S] disagrees with ModeScorer.Vscore by {D['vscan_dev'].max():.1e}")
    os.makedirs(results_dir, exist_ok=True)
    sets_csv = os.path.join(results_dir, SETS_CSV)
    D.to_csv(sets_csv, index=False, float_format="%.6g")
    H = rotation_histograms(D, samples)
    rot_csv = os.path.join(results_dir, ROTATIONS_CSV)
    H.to_csv(rot_csv, index=False, float_format="%.6g")
    out = dict(sets_csv=sets_csv, rotations_csv=rot_csv, figures={},
               close_pairs_not_degenerate=diag["close_pairs_not_degenerate"], c3_character=diag["c3_character"])
    if figures:
        out["figures"]["all_sets"] = plot_all_sets(D, tau, out_dir)
        out["figures"]["near_threshold"] = plot_near_threshold(D, tau, out_dir)
        out["figures"]["mean_comparison"] = plot_mean_comparison(D, tau, out_dir)
        out["figures"]["rotations"] = plot_rotations(D, H, out_dir)
        out["figures"]["spread"] = plot_spread(D, out_dir)
    out.update(headline_numbers(D, samples, tau))
    out["runtime_s"] = time.time() - t0
    return out


if __name__ == "__main__":
    res = run()
    for k_, v_ in res.items():
        if k_ not in ("close_pairs_not_degenerate", "c3_character"):
            print(f"{k_}: {v_}")
