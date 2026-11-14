"""Choice of the principal-axis frame in which the T/R scores are evaluated.

The per-axis scores s[T_Q] and s[R_Q] depend on the frame. For an asymmetric top the
principal axes are unique up to sign; for a symmetric top the two perpendicular axes can
turn freely about the unique axis, for a spherical top all three can, and for a linear
molecule the frame can turn about the molecular axis. This module fixes the frame by a
rule that depends only on the geometry and the input modes, so the same input gives the
same frame in any orientation:

Step 1 (symmetry, symmetric and spherical tops). As in Gaussian's standard orientation:
    z is the C_n axis (the S_4 axis of D2d/S4 tops, which have no C_n with n >= 3); y passes
    through an atom of the smallest set of symmetry-equivalent off-axis atoms; x = y cross z.
    Spherical tops: every set of three mutually perpendicular C2 axes (C4 axes for O/Oh).
    Step 1 builds one frame per set of axis lines (candidate_lines): one per direction to an
    atom of that set, modulo 90 deg -- n of them (n/2 for even n, n/4 if 4 | n) when the set
    projects onto n directions; with an S_2n axis use 2n. Signs and names wait for step 3.
Step 2 (tie-break). Among the candidate sets of axis lines, the one
    with the largest value of the step-2 criterion over the input modes k (before any ideal
    references are added). The rotatable axes Q are x, y (symmetric top), x, y, z
    (spherical top) or y, z (linear molecule). Exact ties (TIE_TOL) go to step 3.
      "tmax" (default):  max over DIFFERENT modes of  sum_Q |s_k(Q)[T_Q]|, i.e. the largest
                         translational score on each rotatable axis, each axis from its own
                         mode (best assignment of the input modes to the in-plane T slots).
      "gap" (option):    sum_k (D_T,k + D_R,k), D = largest minus second-largest |s| of mode k
                         over the rotatable axes (the R1.5 report's sum D).
    Both give the same frame whenever the input contains translational and rotational modes
    (3N input) and for the benzene EMIT modes; with 3N-6 input they can differ.
Step 3 (naming; no physical effect). Among the allowed frames on the winning lines (named_frames;
    on every tied set of lines): +z towards the heaviest on-axis atom off the centre
    (else sum m z^3 > 0); the lowest-index atom in the yz plane on +y; remaining ties by
    the largest coordinates in input order. It fixes axis names and signs only: never a
    score magnitude, a label or the T/R assignment.
Asymmetric tops: eigenvectors of the inertia tensor in increasing moment order
    (x = a, y = b, z = c), right-handed, signs by the step-3 key.
Linear molecules: x = the molecular axis, +x towards the heaviest atom off the centre.
    The turn about x follows the mode vectors ("vector", default): each mode k votes for the
    direction of w_k = sum_A d_k^A projected onto the yz plane, with weight |w_k|. Directions
    count modulo 90 deg (y and z are both axes, and the two members of a degenerate pair
    point at 90 deg to each other). The turn is the direction with the highest density of
    votes (largest summed weight within LINEAR_WINDOW_DEG), averaged over the votes in that
    window; ties: the larger single vote, then the smallest turn. No vote (no mode moves the
    atoms' sum off the axis): +y along the first in-plane atomic displacement. Translations
    carry the largest possible weight (|w| = sqrt(N) for a normalised vector), so 3N input
    follows the T vectors. Option "scan": the turn with the largest step-2 criterion, found
    by rescoring every mode at turns of 0.1 deg (ties: the smallest turn). Either way y and z
    are then named by the modes (the first mode with an in-plane score has its larger
    component on +y).

All functions work in the *input* coordinates, shifted to the centre of mass. A frame F
is a 3x3 matrix whose rows are the new x, y, z axes; coordinates and mode vectors go into
it as X @ F.T and d @ F.T.
"""

import itertools

import numpy as np
from scipy.optimize import linear_sum_assignment

from .scoring import ModeScorer

# Atoms count as symmetry-equivalent within SYM_TOL (Angstrom). Roster tops are symmetric
# to <= 4.5e-6 A in the 6-decimal coordinates Gaussian prints; the ten C2v near-tops miss
# a C_n axis by >= 5e-5 A, so they stay asymmetric. Coordinates printed more coarsely than
# 1e-5 A need a larger value (main.py --sym-tol).
SYM_TOL = 1.5e-5
# sum-D ties: numerically exact ones only, so the choice depends on the scores alone.
TIE_TOL = 1e-6
# Linear molecule: smallest moment <= LIN_TOL * largest (= classifier.LINEAR_TOL).
LIN_TOL = 1e-6
# Spherical top: moments equal to this relative tolerance (and three perpendicular C2 axes).
SPHERICAL_TOL = 1e-4
# Step of the linear-molecule turn scan (degrees; linear rule "scan").
LINEAR_STEP_DEG = 0.1
# Linear rule "vector": votes within this angle (degrees, modulo 90) count as one direction.
LINEAR_WINDOW_DEG = 1.0
# Linear rule "vector": an in-plane |sum_A d^A| at or below this casts no vote.
LINEAR_VOTE_TOL = 1e-6
# Geometric tolerance for "on the axis", "in the plane", equivalent radii/heights (Angstrom).
GEOM_TOL = 1e-3

# Frame rule used by main.build_scorer_and_final: "symmetry" (this module) or "eigh"
# (the submitted version: eigenvectors of the inertia tensor with the heaviest-atom sign fix).
FRAME_RULES = ("symmetry", "eigh")
_frame_rule = "symmetry"
_keep_input_frame = False
# Step-2 criterion (main.py --tie-break): "tmax" (default) or "gap" (sum D).
TIE_BREAKS = ("tmax", "gap")
_tie_break = "tmax"
# Turn of a linear molecule about its axis (main.py --linear-turn): "vector" (default) or "scan".
LINEAR_RULES = ("vector", "scan")
_linear_rule = "vector"


def set_linear_rule(rule):
    global _linear_rule
    if rule not in LINEAR_RULES:
        raise ValueError(f"linear rule must be one of {LINEAR_RULES}, got {rule!r}")
    _linear_rule = rule


def get_linear_rule():
    return _linear_rule


def set_tie_break(rule):
    global _tie_break
    if rule not in TIE_BREAKS:
        raise ValueError(f"tie-break must be one of {TIE_BREAKS}, got {rule!r}")
    _tie_break = rule


def get_tie_break():
    return _tie_break


def set_frame_rule(rule):
    global _frame_rule
    if rule not in FRAME_RULES:
        raise ValueError(f"frame rule must be one of {FRAME_RULES}, got {rule!r}")
    _frame_rule = rule


def get_frame_rule():
    return _frame_rule


def set_keep_input_frame(flag):
    global _keep_input_frame
    _keep_input_frame = bool(flag)


def get_keep_input_frame():
    return _keep_input_frame


def set_sym_tol(tol):
    global SYM_TOL
    if not tol > 0:
        raise ValueError(f"sym_tol must be positive, got {tol!r}")
    SYM_TOL = float(tol)


# ---------------------------------------------------------------- geometry helpers
def inertia(X, m):
    return sum(mi * (r @ r * np.eye(3) - np.outer(r, r)) for mi, r in zip(m, X))


def rotation(u, t):
    """Rotation by angle t about the unit vector u (Rodrigues)."""
    u = u / np.linalg.norm(u)
    K = np.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]])
    return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * K @ K


def maps_onto_itself(X, sym, A, tol):
    """True if the operation A maps every atom onto a distinct atom of the same element within tol."""
    Y = X @ A.T
    used = set()
    other = np.array(sym, dtype=object)
    for i in range(len(X)):
        d = np.linalg.norm(X - Y[i], axis=1)
        d[other != sym[i]] = np.inf
        j = int(np.argmin(d))
        if d[j] > tol or j in used:
            return False
        used.add(j)
    return True


def c2_axes(X, sym, tol):
    """All C2 axes: candidates through an atom or through the midpoint of two like atoms."""
    cand = [x for x in X if np.linalg.norm(x) > 1e-2]
    cand += [(X[i] + X[j]) / 2 for i in range(len(X)) for j in range(i + 1, len(X)) if sym[i] == sym[j]]
    axes = []
    for c in cand:
        n = np.linalg.norm(c)
        if n < 1e-2:
            continue
        u = c / n
        if any(abs(abs(u @ a) - 1) < 1e-6 for a in axes):
            continue
        if maps_onto_itself(X, sym, 2 * np.outer(u, u) - np.eye(3), tol):
            axes.append(u)
    return axes


def detect_top(X, sym, m, tol=None):
    """('linear' | 'spherical' | 'symmetric' | 'asymmetric', unique axis or None, n of C_n).
    Symmetric: a proper C_n with n >= 3, or an S_4 axis (D2d, S4: highest proper axis C_2, returned n = 2)."""
    tol = SYM_TOL if tol is None else tol
    mom, W = np.linalg.eigh(inertia(X, m))
    if mom[0] <= LIN_TOL * mom.max():
        return "linear", W[:, 0], 0
    if (mom.max() - mom.min()) / mom.max() < SPHERICAL_TOL:
        ax = c2_axes(X, sym, tol)
        if any(abs(a @ b) < 1e-6 and abs(a @ c) < 1e-6 and abs(b @ c) < 1e-6
               for a in ax for b in ax for c in ax):
            return "spherical", None, 0
    for k in range(3):
        for n in range(8, 2, -1):
            if maps_onto_itself(X, sym, rotation(W[:, k], 2 * np.pi / n), tol):
                return "symmetric", W[:, k], n
    # S_4 without a proper C_n>=3 (allene, D2d; S4): two equal moments, so the eigenvectors alone would leave the
    # perpendicular axes arbitrary. (S_2n for n >= 3 contains C_n and is found above.)
    for k in range(3):
        u = W[:, k]
        if maps_onto_itself(X, sym, (np.eye(3) - 2 * np.outer(u, u)) @ rotation(u, np.pi / 2), tol):
            return "symmetric", u, 2
    return "asymmetric", None, 0


# the four right-handed sign patterns of a frame
SIGNS = [np.diag(d) for d in ((1, 1, 1), (-1, -1, 1), (-1, 1, -1), (1, -1, -1))]


def smallest_set(X, sym, m, z):
    """Atoms of the smallest set of equivalent off-axis atoms (same element, radius and
    |height|); ties: heaviest, then nearest the axis; the union if still tied."""
    zc = X @ z
    P = X - np.outer(zc, z)
    rad = np.linalg.norm(P, axis=1)
    groups = []
    for i in (i for i in range(len(X)) if rad[i] > GEOM_TOL):
        for g in groups:
            j = g[0]
            if sym[i] == sym[j] and abs(rad[i] - rad[j]) < GEOM_TOL and abs(abs(zc[i]) - abs(zc[j])) < GEOM_TOL:
                g.append(i)
                break
        else:
            groups.append([i])
    key = lambda g: (len(g), -m[g[0]], round(rad[g[0]], 3))
    k0 = min(key(g) for g in groups)
    return [i for g in groups if key(g) == k0 for i in g], P, rad


def _dedup(frames):
    uniq = []
    for f in frames:
        if not any(np.allclose(f, g, atol=1e-8) for g in uniq):
            uniq.append(f)
    return uniq


def candidate_lines(X, sym, m, kind, z, tol=None):
    """Step 1: one frame (rows = x, y, z) per set of candidate axis lines of a symmetric or spherical
    top. Signs and names are left to step 3 (named_frames); they never change an |s|.
      symmetric: z on the C_n axis, y through an atom of the smallest set, x = y cross z -- one set per
                 direction to those atoms modulo 90 deg (180 deg reverses y, 90 deg swaps x and y);
      spherical: one set per triple of perpendicular C2 axes (C4 for O/Oh): one for T..Oh, five for I/Ih."""
    tol = SYM_TOL if tol is None else tol
    if kind == "symmetric":
        atoms, P, rad = smallest_set(X, sym, m, z)
        cands = [np.array([np.cross(P[i] / rad[i], z), P[i] / rad[i], z]) for i in atoms]
    elif kind == "spherical":
        ax = c2_axes(X, sym, tol)
        c4 = [a for a in ax if maps_onto_itself(X, sym, rotation(a, np.pi / 2), tol)]
        pool = c4 if len(c4) >= 3 else ax
        cands = [np.array([a, b, c if np.linalg.det([a, b, c]) > 0 else -c])
                 for a, b, c in itertools.combinations(pool, 3)
                 if abs(a @ b) < 1e-6 and abs(a @ c) < 1e-6 and abs(b @ c) < 1e-6]
    else:
        return []
    reps, seen = [], set()
    for f in cands:
        k = line_set(f)
        if k not in seen:
            seen.add(k)
            reps.append(f)
    return reps


def named_frames(X, sym, m, kind, f):
    """Step 3's choices: every allowed right-handed frame on the axis lines of f.
      symmetric: z stays on the C_n axis; y along either in-plane line of f that passes through an atom
                 of the smallest set (both lines when 4 | n), with the 4 right-handed sign patterns;
      spherical: the 3! orderings x 4 right-handed sign patterns of the three lines."""
    if kind == "symmetric":
        z = f[2]
        atoms, P, rad = smallest_set(X, sym, m, z)
        dirs = P[atoms] / rad[atoms, None]
        out = []
        for y in (f[0], f[1]):
            if np.any(np.abs(dirs @ y) > 1 - 1e-6):
                base = np.array([np.cross(y, z), y, z])
                out += [D @ base for D in SIGNS]
    else:
        out = [g for g in (np.array([sg[k] * f[perm[k]] for k in range(3)])
                           for perm in itertools.permutations(range(3))
                           for sg in itertools.product((1, -1), repeat=3))
               if np.linalg.det(g) > 0]
    return _dedup(out)


def candidate_frames(X, sym, m, kind, z, tol=None):
    """All allowed right-handed frames of a symmetric or spherical top: the step-3 names of every step-1
    set of axis lines (for analyses; select_frame scores one frame per set)."""
    return [g for f in candidate_lines(X, sym, m, kind, z, tol) for g in named_frames(X, sym, m, kind, f)]


def _line_key(r):
    """One key per line (r and -r): the sign that makes the first clearly nonzero component positive.
    (Signing by the largest component fails for a line with two equal largest components, e.g. at 45 deg.)"""
    i = int(np.flatnonzero(np.abs(r) > 1e-4)[0])
    return tuple(np.round(r if r[i] > 0 else -r, 4) + 0.0)


def line_set(f):
    """The unordered set of axis lines of a frame (frames on the same lines give the same |s| up to names)."""
    return tuple(sorted(_line_key(r) for r in f))


def convention_key(X, sym, m, f):
    """Step 3 (naming). Smallest key wins:
      1) +z towards the heaviest on-axis atom OFF the centre (lowest index among equals); if none, sum m z^3 > 0;
      2) the lowest-index atom in the yz plane on +y;
      3) remaining ties: the largest coordinates in input order (x1, y1, z1, x2, ...) on a 0.01 A grid.
    Atom indices enter only where the alternatives are symmetry-equivalent."""
    Y = X @ f.T
    rad = np.hypot(Y[:, 0], Y[:, 1])
    on = [i for i in range(len(X)) if rad[i] < GEOM_TOL and abs(Y[i, 2]) > GEOM_TOL]
    if on:
        zkey = -float(np.sign(Y[max(on, key=lambda i: (m[i], -i)), 2]))
    else:
        mz3 = float(np.dot(m, Y[:, 2] ** 3))
        zkey = -float(np.sign(mz3)) if abs(mz3) > 1e-3 else 0.0
    ony = [i for i in range(len(X)) if abs(Y[i, 0]) < GEOM_TOL and Y[i, 1] > GEOM_TOL]
    lex = tuple(-np.round(Y / 1e-2).ravel())
    return (zkey, min(ony) if ony else len(X), lex)


# ---------------------------------------------------------------- scores in a trial frame
def _scores_in_frame(X, sym, bonds, vecs, f):
    sc = ModeScorer(sym, X @ f.T, bonds)
    return [sc.calculate_scores(v @ f.T) for v in vecs]


def sum_d(X, sym, bonds, vecs, f, axes="xy"):
    """sum over the modes of (D_T + D_R); D = largest minus second-largest |s| over `axes`."""
    tot = 0.0
    for r in _scores_in_frame(X, sym, bonds, vecs, f):
        for K in ("T", "R"):
            t = sorted((abs(r[K][q]) for q in axes), reverse=True)
            tot += t[0] - t[1]
    return float(tot)


def t_max(X, sym, bonds, vecs, f, axes="xy"):
    """largest sum over `axes` of |s[T_Q]|, each axis taken from a different mode
    (best assignment of the modes to the T slots of `axes`)."""
    A = np.array([[abs(r["T"][q]) for q in axes] for r in _scores_in_frame(X, sym, bonds, vecs, f)])
    if A.size == 0:
        return 0.0
    rows, cols = linear_sum_assignment(-A)
    return float(A[rows, cols].sum())


def step2_value(X, sym, bonds, vecs, f, axes="xy", tie_break=None):
    """The step-2 criterion of frame f (larger wins): t_max ("tmax") or sum_d ("gap")."""
    tie_break = _tie_break if tie_break is None else tie_break
    if tie_break not in TIE_BREAKS:
        raise ValueError(f"tie-break must be one of {TIE_BREAKS}, got {tie_break!r}")
    fn = t_max if tie_break == "tmax" else sum_d
    return fn(X, sym, bonds, vecs, f, axes)


# ---------------------------------------------------------------- asymmetric tops
def frame_asymmetric(X, sym, m):
    """Eigenvectors of I in increasing moment order, right-handed, signs by the naming key."""
    _, W = np.linalg.eigh(inertia(X, m))
    base = W.T.copy()
    if np.linalg.det(base) < 0:
        base[2] *= -1
    return min((D @ base for D in SIGNS), key=lambda f: convention_key(X, sym, m, f))


# ---------------------------------------------------------------- linear molecules
def _linear_start(X, m):
    """x = molecular axis towards the heaviest atom off the centre (lowest index among equals);
    u = the next eigenvector, made perpendicular to x (the start of the turn scan)."""
    _, W = np.linalg.eigh(inertia(X, m))
    ax = W[:, 0]
    off = [i for i in range(len(X)) if abs(X[i] @ ax) > GEOM_TOL]
    if off:
        i0 = max(off, key=lambda i: (m[i], -i))
        if X[i0] @ ax < 0:
            ax = -ax
    u = W[:, 1] - (W[:, 1] @ ax) * ax
    return ax, u / np.linalg.norm(u)


def _linear_frame_at(ax, u, phi):
    v = np.cross(ax, u)
    c, s_ = np.cos(phi), np.sin(phi)
    return np.array([ax, c * u + s_ * v, -s_ * u + c * v])          # right-handed, x = molecular axis


def linear_votes(vecs, ax, u, tol=None):
    """One vote per mode: (|w|, angle of w in [0, pi/2)) for w = sum_A d^A projected onto the
    plane perpendicular to ax, the angle measured from u towards ax x u. Modes with
    |w| <= tol cast no vote."""
    tol = LINEAR_VOTE_TOL if tol is None else tol
    v = np.cross(ax, u)
    out = []
    for d in vecs:
        w = np.asarray(d, dtype=float).sum(axis=0)
        a, b = float(w @ u), float(w @ v)
        r = float(np.hypot(a, b))
        if r > tol:
            out.append((r, float(np.arctan2(b, a) % (np.pi / 2))))
    return out


def vector_turn(vecs, ax, u, window_deg=None):
    """The turn (radians, in [0, pi/2)) with the highest density of mode votes, averaged over
    the votes within the window of the densest one. Returns (phi or None, density, n_votes)."""
    win = np.radians(LINEAR_WINDOW_DEG if window_deg is None else window_deg)
    votes = linear_votes(vecs, ax, u)
    if not votes:
        return None, 0.0, 0
    w = np.array([r for r, _ in votes])
    th = np.array([t for _, t in votes])
    q = np.pi / 2
    dist = np.abs((th[:, None] - th[None, :] + q / 2) % q - q / 2)       # angular distance modulo 90 deg
    near = dist <= win + 1e-12
    dens = near @ w
    best = min(range(len(votes)), key=lambda k: (-round(float(dens[k]), 9), -round(float(w[k]), 9), th[k]))
    sel = near[best]
    z = np.sum(w[sel] * np.exp(4j * th[sel]))                          # 4-fold mean of the peak's votes
    return float((np.angle(z) / 4) % q), float(dens[best]), len(votes)


def _linear_named(X, sym, bonds, vecs, ax, u, phi):
    """Of the four frames phi + k*90 deg, the one whose first mode with an in-plane score
    (|s_y| != |s_z|) has its larger component on +y, positive. If no mode has one (e.g. the
    bends of HCCH, whose in-plane T and R scores all vanish), the first in-plane atomic
    displacement (first mode in input order, first moving atom) decides the same way, so
    the names never depend on the input orientation."""
    cands = [_linear_frame_at(ax, u, phi + k * np.pi / 2) for k in range(4)]
    S = [_scores_in_frame(X, sym, bonds, vecs, f) for f in cands]
    for kmode in range(len(vecs)):
        for K in ("T", "R"):
            a0, b0 = S[0][kmode][K]["y"], S[0][kmode][K]["z"]
            if abs(abs(a0) - abs(b0)) > 1e-6:
                ok = [abs(s[kmode][K]["y"]) > abs(s[kmode][K]["z"]) and s[kmode][K]["y"] > 0 for s in S]
                return cands[int(np.argmax(ok))]
    for v in vecs:
        P = v - np.outer(v @ ax, ax)
        n_ = np.linalg.norm(P, axis=1)
        if n_.max() > 1e-6:
            p = P[int(np.argmax(n_ > 1e-6))]
            ok = [abs(p @ f[1]) > abs(p @ f[2]) + 1e-9 and p @ f[1] > 0 for f in cands]
            return cands[int(np.argmax(ok))]
    return cands[0]


def _linear_no_turn(ax, u, vecs, tie_break, scorings, why):
    """Nothing fixes the turn: +y along the perpendicular displacement of the first moving
    atom (first mode in input order)."""
    for v in vecs:
        P = v - np.outer(v @ ax, ax)
        n_ = np.linalg.norm(P, axis=1)
        if n_.max() > 1e-6:
            y = P[int(np.argmax(n_ > 1e-6))]
            y = y / np.linalg.norm(y)
            return (np.array([ax, y, np.cross(ax, y)]),
                    dict(phi=None, value=0.0, tie_break=tie_break, scorings=scorings, step=f"linear: {why}"))
    return (_linear_frame_at(ax, u, 0.0),
            dict(phi=None, value=0.0, tie_break=tie_break, scorings=scorings, step="linear: no in-plane displacement"))


def frame_linear(X, sym, m, bonds, vecs, step_deg=None, tie_break=None, linear_rule=None):
    """Turn about the molecular axis, from the mode vectors ("vector") or by a rescoring
    scan ("scan"). Returns (frame, info)."""
    tie_break = _tie_break if tie_break is None else tie_break
    linear_rule = _linear_rule if linear_rule is None else linear_rule
    if linear_rule not in LINEAR_RULES:
        raise ValueError(f"linear rule must be one of {LINEAR_RULES}, got {linear_rule!r}")
    ax, u = _linear_start(X, m)
    if linear_rule == "vector":
        phi, dens, nv = vector_turn(vecs, ax, u)
        if phi is None:
            return _linear_no_turn(ax, u, vecs, tie_break, 0, "no vote (vector)")
        info = dict(phi=float(np.degrees(phi)), value=dens, tie_break=tie_break, scorings=4,
                    step=f"linear: vector ({nv} votes)")
        return _linear_named(X, sym, bonds, vecs, ax, u, phi), info
    step = LINEAR_STEP_DEG if step_deg is None else step_deg
    turns = np.radians(np.arange(0.0, 90.0, step))
    sums = np.array([step2_value(X, sym, bonds, vecs, _linear_frame_at(ax, u, t), "yz", tie_break) for t in turns])
    if sums.max() < TIE_TOL:
        # No in-plane score (for "tmax": no in-plane T score) at all: the criterion is the
        # same at every turn.
        return _linear_no_turn(ax, u, vecs, tie_break, len(turns), "no in-plane score")
    i = int(np.argmax(sums >= sums.max() - TIE_TOL))
    phi = turns[i]
    info = dict(phi=float(np.degrees(phi)), value=float(sums[i]), tie_break=tie_break, scorings=len(turns) + 4,
                step=f"linear: scan ({tie_break})")
    return _linear_named(X, sym, bonds, vecs, ax, u, phi), info


# ---------------------------------------------------------------- the rule
def select_frame(X, sym, m, bonds, vecs, keep_input_frame=None, tol=None, tie_break=None, linear_rule=None):
    """Choose the frame for COM-shifted input coordinates X and the input mode vectors `vecs`
    (as given: no ideal references). Returns (F, info); info has kind, n, evaluations,
    tie_break (the step-2 criterion), values (its value per set of axis lines, descending)
    and the step that decided."""
    keep_input_frame = _keep_input_frame if keep_input_frame is None else keep_input_frame
    tie_break = _tie_break if tie_break is None else tie_break
    kind, z, n = detect_top(X, sym, m, tol)
    if kind == "asymmetric":
        return frame_asymmetric(X, sym, m), dict(kind=kind, n=0, evaluations=0, tie_break=tie_break,
                                                 values=[], step="asymmetric: eigenvectors + naming")
    if kind == "linear":
        f, info = frame_linear(X, sym, m, bonds, vecs, tie_break=tie_break, linear_rule=linear_rule)
        return f, dict(kind=kind, n=0, evaluations=info["scorings"], tie_break=tie_break, values=[info["value"]],
                       step=info["step"])
    axes = "xy" if kind == "symmetric" else "xyz"
    reps = candidate_lines(X, sym, m, kind, z, tol)                       # step 1: one frame per set of lines
    if len(reps) == 1:
        score, nev = [0.0], 0
    else:
        score = [step2_value(X, sym, bonds, vecs, f, axes, tie_break) for f in reps]   # step 2
        nev = len(reps)
    best = max(score)
    winners = [f for f, s in zip(reps, score) if s >= best - TIE_TOL]
    tied = len(winners)
    info = dict(kind=kind, n=n, evaluations=nev, tie_break=tie_break, values=sorted(score, reverse=True))
    if keep_input_frame:
        for f in reps:
            for g in named_frames(X, sym, m, kind, f):
                if np.allclose(g, np.eye(3), atol=1e-4):
                    return g, dict(info, step=f"{kind}: input frame kept")
    keep = [g for f in winners for g in named_frames(X, sym, m, kind, f)]   # step 3 names the winner
    if len(reps) == 1:
        step = f"{kind}: one set of axis lines + naming"
    elif tied == 1:
        step = f"{kind}: {tie_break} + naming"
    else:
        step = f"{kind}: {tie_break} tie + naming"
    return min(keep, key=lambda f: convention_key(X, sym, m, f)), dict(info, step=step)


def select_frame_for_scorer(scorer, vecs, keep_input_frame=None, tol=None, tie_break=None, linear_rule=None):
    """select_frame on a ModeScorer's current (COM-shifted) coordinates."""
    return select_frame(scorer.coords, scorer.symbols, scorer.masses, scorer.bList, vecs,
                        keep_input_frame=keep_input_frame, tol=tol, tie_break=tie_break, linear_rule=linear_rule)
