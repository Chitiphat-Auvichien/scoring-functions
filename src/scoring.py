import numpy as np
import math
from .utils import atomicMass

# --- Centralized numerical constants (JCC spec conventions) ---
EPS_DISP = 1e-8    # spec eps_disp: |d_A| <= this -> zero-motion (unit(0):=0)
EPS_NORM = 1e-9    # vector-normalization guard (unit(0):=0)
EPS_DENOM = 1e-6   # V-score denominator guard (Sigma|db|^2 near zero -> score 0)
RANGE_TOL = 1e-6   # tolerance for the range-invariant score asserts

# --- V-score bond-weighting variants (eq:vscore/eq:bondscore) ---
# 'none' is the original definition (weight every bond by |Delta b|^2 alone);
# 'mu' additionally weights each bond by its reduced mass mu_AB, so a bond's
# contribution tracks the kinetic energy of its relative stretching motion
# (1/2 mu |db/dt|^2) rather than raw displacement. mu appears in BOTH the
# numerator and the denominator, so s[V_S] stays a weighted mean of per-bond
# |cos| and both invariants survive: s[V_S] in [0,1] and Sum|s_AB| == s[V_S].
V_WEIGHTINGS = ("mu", "none")
DEFAULT_V_WEIGHTING = "mu"
_ACTIVE_V_WEIGHTING = DEFAULT_V_WEIGHTING


def resolve_v_weighting(name=None):
    """Validate a weighting name; None means 'use the process-active default'."""
    resolved = _ACTIVE_V_WEIGHTING if name is None else name
    if resolved not in V_WEIGHTINGS:
        raise ValueError(
            f"Unknown V-score weighting {resolved!r}; expected one of {V_WEIGHTINGS}")
    return resolved


def set_v_weighting(name):
    """Set the process-wide default weighting. Returns the previous value so
    callers (and tests) can restore it. main.py calls this once, right after
    parse_args(), so every downstream module -- library_ingest, calibrate,
    merge_ped_scores, compare_rerun -- picks the variant up without needing the
    option threaded through its own signature."""
    global _ACTIVE_V_WEIGHTING
    previous = _ACTIVE_V_WEIGHTING
    _ACTIVE_V_WEIGHTING = resolve_v_weighting(name)
    return previous


def get_v_weighting():
    """The process-active default, without constructing a ModeScorer."""
    return _ACTIVE_V_WEIGHTING


def size_vec(v):
    return math.sqrt(np.dot(v, v))


# --- Main Scorer Class ---

class ModeScorer:
    def __init__(self, atom_symbols, coords, bonds, v_weighting=None):
        """Initialize from parser output (atom symbols, coords, 0-based bond index pairs).

        v_weighting selects the eq:vscore bond weighting ('mu' or 'none');
        None means the process-active default (see set_v_weighting).
        """
        self.n = len(atom_symbols)
        self.symbols = list(atom_symbols)

        masses = []
        for sym in self.symbols:
            key = sym.lower()
            if key not in atomicMass:
                raise ValueError(f"Unrecognized element symbol '{sym}': no atomic mass on file")
            masses.append(atomicMass[key])
        self.masses = np.array(masses, dtype=float)

        self.coords = np.array(coords, dtype=float)  # (n, 3)
        self.dispVecs = np.zeros((self.n, 3))         # updated per mode, see calculate_scores()
        self.dispLengths = np.zeros(self.n)

        # Bond Setup
        self.nBond = len(bonds)
        self.bList = bonds
        self.bVec = []

        # Initial bond calculation
        self.update_bond_vectors()

        # eq:vscore bond weights -- needs bList and masses, both set above.
        self.v_weighting = resolve_v_weighting(v_weighting)
        self.bond_weights = self._bond_weights()

        # Move to Center of Mass
        self.COM()

    def _bond_weights(self):
        """Per-bond weight w_AB entering eq:vscore's numerator and denominator.

        'none' -> all 1.0, i.e. the original definition where a bond's weight is
        |Delta b_AB|^2 alone.
        'mu'   -> the reduced mass mu_AB = m_A*m_B/(m_A+m_B), rescaled by its
        maximum so w in (0, 1].

        The rescale is a mathematical no-op (the same constant divides both the
        numerator and the denominator of eq:vscore) but it earns its keep twice:
        it makes every bond of a homoleptic AB_n molecule weigh exactly 1.0, so
        such molecules score identically under both variants, and it keeps the
        summed denominator on the same scale as the unweighted one, so the
        EPS_DENOM zero-motion guard keeps its calibrated meaning.
        """
        if self.v_weighting == "none":
            return np.ones(self.nBond)
        mu = np.array([
            self.masses[i] * self.masses[j] / (self.masses[i] + self.masses[j])
            for (i, j) in self.bList
        ], dtype=float)
        if mu.size == 0:
            return mu
        return mu / mu.max()

    def update_bond_vectors(self):
        """Recalculate bond vectors based on current atom positions."""
        self.bVec = [self.coords[j] - self.coords[i] for (i, j) in self.bList]

    def COM(self):
        """Calculate Center of Mass and translate molecule."""
        total_mass = self.masses.sum()
        weighted = (self.masses[:, None] * self.coords).sum(axis=0)
        com = weighted / total_mass if total_mass > 0 else weighted
        self.coords = self.coords - com

        # Update bonds after translation (vectors shouldn't change, but good practice)
        self.update_bond_vectors()

    def _build_inertia_tensor(self):
        """Build the moment-of-inertia tensor at the current geometry (shared
        by MIT() and the principal_axes() accessor)."""
        x, y, z = self.coords[:, 0], self.coords[:, 1], self.coords[:, 2]
        m = self.masses

        XX = np.sum(m * (y**2 + z**2))
        YY = np.sum(m * (x**2 + z**2))
        ZZ = np.sum(m * (x**2 + y**2))
        XY = -np.sum(m * x * y)
        XZ = -np.sum(m * x * z)
        YZ = -np.sum(m * y * z)

        return np.array([
            [XX, XY, XZ],
            [XY, YY, YZ],
            [XZ, YZ, ZZ]
        ])

    def principal_axes(self):
        """Principal moments (ascending eigenvalues) and axes (eigenvector
        columns) of the inertia tensor at the current geometry. Moments are
        rotation-invariant, so this is consistent before or after MIT()."""
        tensor = self._build_inertia_tensor()
        moments, axes = np.linalg.eigh(tensor)
        return moments, axes

    def MIT(self, modes=None, rotate_modes=True):
        """Rotates the molecule (and, if rotate_modes, the mode displacement
        vectors) into the basis of principal axes of inertia."""
        tensor = self._build_inertia_tensor()
        eigVal, rot = np.linalg.eigh(tensor)  # eigenvectors as columns of rot

        # Sign-fix heuristic (ported from legacy atom.py, not independently
        # derived here): orient axes so the heaviest atom's projected
        # coordinates sum positive.
        heaviest_idx = int(np.argmax(self.masses))  # first max on ties, matches a plain scan
        h_new = self.coords[heaviest_idx] @ rot
        if h_new.sum() < 0.0:
            rot = -rot

        self.coords = self.coords @ rot
        self.update_bond_vectors()  # recompute rather than rotate existing vectors

        if modes is not None:
            if rotate_modes:
                rotated_modes = []
                for mode in modes:
                    new_mode = mode.copy()
                    new_mode['vector'] = mode['vector'] @ rot  # shape (N, 3) @ (3, 3)
                    rotated_modes.append(new_mode)
                return rotated_modes
            else:
                return modes
        return None

    def construct_T(self):
        """Constructs the 3 ideal translational reference modes (Tx, Ty, Tz)."""
        modes = []
        labels = ['Tx', 'Ty', 'Tz']

        for i in range(3):
            vec = np.zeros((self.n, 3))
            vec[:, i] = 1.0

            flat_norm = np.linalg.norm(vec)
            if flat_norm > EPS_NORM:
                vec = vec / flat_norm

            modes.append({
                "frequency": 0.0,
                "vector": vec,
                "label": labels[i]
            })
        return modes

    def construct_R(self):
        """Constructs the 3 ideal rotational reference modes (Rx, Ry, Rz);
        tangent direction per atom is (0,-z,y)/(z,0,-x)/(-y,x,0)."""
        self.COM()  # ensure we are at COM
        modes = []
        labels = ['Rx', 'Ry', 'Rz']

        x, y, z = self.coords[:, 0], self.coords[:, 1], self.coords[:, 2]
        zeros = np.zeros(self.n)
        rx_vecs = np.stack([zeros, -z, y], axis=1)
        ry_vecs = np.stack([z, zeros, -x], axis=1)
        rz_vecs = np.stack([-y, x, zeros], axis=1)

        for vecs, lbl in zip([rx_vecs, ry_vecs, rz_vecs], labels):
            flat_norm = np.linalg.norm(vecs)
            if flat_norm > EPS_NORM:
                vecs = vecs / flat_norm
            
            modes.append({
                "frequency": 0.0,
                "vector": vecs,
                "label": lbl
            })
            
        return modes

    def calculate_scores(self, mode_vector):
        """Load a mode's displacement vector (N_atoms, 3) and calculate all scores."""
        self.dispVecs = np.asarray(mode_vector, dtype=float)
        self.dispLengths = np.sqrt(np.sum(self.dispVecs**2, axis=1))

        scores = {
            "T": self.Tscore(),
            "R": self.Rscore(),
            "V": self.Vscore()
        }
        self._assert_score_ranges(scores)
        return scores

    @staticmethod
    def _assert_score_ranges(scores):
        """Range-invariant guards from the spec: s[T],s[R] in [-1,1]; s[V_S] in [0,1]."""
        for kind in ("T", "R"):
            for axis, val in scores[kind].items():
                if not (-1.0 - RANGE_TOL <= val <= 1.0 + RANGE_TOL):
                    raise ValueError(f"s[{kind}_{axis}]={val} out of [-1,1]")
        vs = scores["V"]
        if not (-RANGE_TOL <= vs <= 1.0 + RANGE_TOL):
            raise ValueError(f"s[V_S]={vs} out of [0,1]")

    def Tscore(self):
        """Calculates Translational Scores (Tx, Ty, Tz)."""
        n = self.n
        # EPS_DENOM (not the looser EPS_DISP) is the noise floor here:
        # EPS_DISP=1e-8 is too permissive against the ~1e-8-1e-6 noise
        # Gaussian prints for atoms symmetry-required to be exactly zero
        # in degenerate EMIT eigenvectors, which would otherwise be
        # promoted to a full-weight unit-vector contribution.
        mask = self.dispLengths > EPS_DENOM
        unit_disp = np.zeros((n, 3))
        unit_disp[mask] = self.dispVecs[mask] / self.dispLengths[mask, None]
        Tx, Ty, Tz = unit_disp.sum(axis=0)

        return {
            'x': Tx * (1.0/float(n)),
            'y': Ty * (1.0/float(n)),
            'z': Tz * (1.0/float(n))
        }

    def Rscore(self):
        """Rotational scores s[R_x], s[R_y], s[R_z] (eq:rscore).

        Per atom and axis Q, normalize r_perp = r-(r.Qhat)Qhat and the
        displacement d SEPARATELY, cross them, and take the Q-component:
            s[R_Q] = (1/(N-N_Q)) sum_offaxis (r_perp x d)_Q / (|r_perp| |d|)
        Normalizing r_perp and d separately (rather than by |r_perp x d|)
        keeps the sin(phi) factor, down-weighting motion that isn't purely
        tangential. N_Q = on-axis atoms (|r_perp| ~ 0), excluded.
        """
        axes = (np.array([1.0, 0.0, 0.0]),
                np.array([0.0, 1.0, 0.0]),
                np.array([0.0, 0.0, 1.0]))
        out = {}
        for key, Q in zip('xyz', axes):
            r_perp = self.coords - np.outer(self.coords @ Q, Q)
            lr = np.sqrt(np.sum(r_perp**2, axis=1))
            off_mask = lr > EPS_DENOM              # on the axis -> excluded (N_Q)
            n_off = int(np.sum(off_mask))

            ld = self.dispLengths
            d_mask = off_mask & (ld > EPS_DENOM)    # unit(0):=0 otherwise
            cross = np.cross(r_perp[d_mask], self.dispVecs[d_mask])
            total = np.sum((cross @ Q) / (lr[d_mask] * ld[d_mask]))

            out[key] = total / n_off if n_off > 0 else 0.0
        return out

    def _bond_contributions(self):
        """Per-bond pieces of the V-score (eq:vscore numerator/denominator),
        plus the diagnostic signed relative bond-length change. Every term
        carries the per-bond weight w_AB from _bond_weights() (1 under the
        'none' variant, the rescaled reduced mass mu_AB under 'mu'). Returns
        four parallel lists over self.bList:
          terms        : w_AB*|Δb_AB|^2 * |unit(Δb_AB).b_hat_AB|  (numerator
                         term, magnitude -- feeds s[V_S] via Vscore())
          signed_terms : w_AB*|Δb_AB|^2 * (unit(Δb_AB).b_hat_AB)  (signed
                         version of the same term, no outer abs() -- feeds the
                         signed per-bond s_AB in score_bonds(); positive =
                         stretching, negative = compressing)
          sqdisps      : w_AB*|Δb_AB|^2                        (denominator term)
          rel_db       : (|b_AB+Δd_AB| - |b_AB|) / |b_AB|. NOT part of
                         eq:vscore/eq:bondscore -- a purely geometric
                         diagnostic, so it stays UNWEIGHTED
                         (fig:bondscores' x-axis).
        Vscore() and score_bonds() both build on this, so the weight can never
        get out of step between them. s_AB is signed but s[V_S] is not:
        Σ|s_AB| == s[V_S] exactly (not Σ s_AB) -- w_AB divides out of that
        identity because it scales numerator and denominator alike.
        """
        terms = []
        signed_terms = []
        sqdisps = []
        rel_db = []
        for b in range(self.nBond):
            idx1, idx2 = self.bList[b]

            bVec_curr = self.bVec[b]  # equilibrium bond vector
            bLength = size_vec(bVec_curr)

            delDisp = self.dispVecs[idx2] - self.dispVecs[idx1]
            delDispLength = size_vec(delDisp)

            dot_val = np.dot(delDisp, bVec_curr)

            w = self.bond_weights[b]
            terms.append(w * abs(dot_val) * delDispLength / bLength)
            signed_terms.append(w * dot_val * delDispLength / bLength)
            sqdisps.append(w * delDispLength**2)

            if bLength > EPS_NORM:
                rel_db.append((size_vec(bVec_curr + delDisp) - bLength) / bLength)
            else:
                rel_db.append(0.0)
        return terms, signed_terms, sqdisps, rel_db

    def Vscore(self):
        """Calculates the Vibrational Score s[V_S] (eq:vscore):

            s[V_S] = Σ_b w_b*|Δb_b|^2*|unit(Δb_b).b_hat_b| / Σ_b w_b*|Δb_b|^2

        i.e. a weighted mean of the per-bond directional cosines, hence always
        in [0,1] whatever the (non-negative) weights. w_b is 1 under the 'none'
        variant and the reduced mass mu_AB under 'mu' -- see _bond_weights().
        """
        terms, _signed_terms, sqdisps, _ = self._bond_contributions()
        modeScr = sum(terms)
        denom = sum(sqdisps)

        if denom > EPS_DENOM:
            return modeScr / denom
        return 0.0

    def score_bonds(self):
        """Per-bond stretch contribution s_AB (eq:bondscore), signed:
        s_AB = w_AB*|Δb_AB|^2 * (unit(Δb_AB).b_hat_AB) / Σ_bonds w_b*|Δb_b|^2

        where w is the eq:vscore bond weight (1, or the reduced mass mu_AB
        under the 'mu' variant -- see _bond_weights()).

        Positive s_AB means the bond is stretching, negative means it is
        compressing; s[V_S] (Vscore()) itself remains non-negative and is
        computed from the magnitude internally, unaffected by this sign.

        Returns a list of dicts {'i','j','s_AB','rel_db','i_label','j_label'}
        aligned with self.bList ('i_label'/'j_label' are e.g. "C1", "H7" for
        readable bond identifiers). Asserts Σ|s_AB| == s[V_S] to 1e-6. Call
        calculate_scores() first so the atom dispVecs are populated.
        """
        _terms, signed_terms, sqdisps, rel_db = self._bond_contributions()
        denom = sum(sqdisps)

        bonds = []
        for b in range(self.nBond):
            idx1, idx2 = self.bList[b]
            s_AB = signed_terms[b] / denom if denom > EPS_DENOM else 0.0
            bonds.append({
                "i": idx1, "j": idx2, "s_AB": s_AB, "rel_db": rel_db[b],
                "i_label": f"{self.symbols[idx1]}{idx1 + 1}",
                "j_label": f"{self.symbols[idx2]}{idx2 + 1}",
            })

        total = sum(abs(bd["s_AB"]) for bd in bonds)
        vs = self.Vscore()
        if abs(total - vs) > 1e-6:
            raise ValueError(f"Sum of |per-bond s_AB| ({total}) != s[V_S] ({vs})")
        return bonds


def format_bond_map(bonds, value_key, dp=4):
    """Semicolon-joined 'iLabel-jLabel:value' string from a score_bonds()-
    style list of dicts, e.g. 'C1-C2:0.0342;C1-C6:0.0539'."""
    return ";".join(f"{b['i_label']}-{b['j_label']}:{b[value_key]:.{dp}f}" for b in bonds)


def parse_bond_string(s):
    """Inverse of format_bond_map: '"C1-C2:0.0342;C1-C6:0.0539;..."' ->
    {'C1-C2': 0.0342, ...}."""
    out = {}
    if not isinstance(s, str) or not s:
        return out
    for part in s.split(";"):
        atoms, val = part.split(":")
        out[atoms] = float(val)
    return out