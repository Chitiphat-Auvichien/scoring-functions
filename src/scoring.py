import numpy as np
import math
from .utils import atomicMass

# --- Centralized numerical constants (JCC spec conventions) ---
EPS_DISP = 1e-8    # spec eps_disp: |d_A| <= this -> zero-motion (unit(0):=0)
EPS_NORM = 1e-9    # vector-normalization guard (unit(0):=0)
EPS_DENOM = 1e-6   # V-score denominator guard (Sigma|db|^2 near zero -> score 0)
RANGE_TOL = 1e-6   # tolerance for the range-invariant score asserts

# --- Helper Classes to mimic atom.py structure ---

def sizeVec(v):
    return math.sqrt(np.dot(v, v))

class Coordinate:
    def __init__(self, x, y, z):
        self.X = x
        self.Y = y
        self.Z = z

class Atom:
    def __init__(self, element, x, y, z):
        self.symbol = element
        key = self.symbol.lower()
        if key not in atomicMass:
            raise ValueError(f"Unrecognized element symbol '{element}': no atomic mass on file")
        self.rMass = float(atomicMass[key])
        self.coord = Coordinate(x, y, z)
        # These will be updated for each mode
        self.dispVec = np.zeros(3)
        self.dispLength = 0.0

    # Helper accessors used in atom.py calculations
    def x(self): return self.coord.X
    def y(self): return self.coord.Y
    def z(self): return self.coord.Z

    def translation(self, XC, YC, ZC):
        self.coord.X -= XC
        self.coord.Y -= YC
        self.coord.Z -= ZC

# --- Main Scorer Class ---

class ModeScorer:
    def __init__(self, atom_symbols, coords, bonds):
        """Initialize from parser output (atom symbols, coords, 0-based bond index pairs)."""
        self.n = len(atom_symbols)
        self.atoms = []
        
        # Create Atom objects
        for i in range(self.n):
            sym = atom_symbols[i]
            x, y, z = coords[i]
            self.atoms.append(Atom(sym, x, y, z))
            
        # Bond Setup
        self.nBond = len(bonds)
        self.bList = bonds
        self.bVec = []
        
        # Initial bond calculation
        self.update_bond_vectors()

        # Move to Center of Mass
        self.COM()

    def update_bond_vectors(self):
        """Recalculate bond vectors based on current atom positions."""
        self.bVec = []
        for (i, j) in self.bList:
            vec = np.array([
                self.atoms[j].x() - self.atoms[i].x(),
                self.atoms[j].y() - self.atoms[i].y(),
                self.atoms[j].z() - self.atoms[i].z()
            ])
            self.bVec.append(vec)

    def COM(self):
        """Calculate Center of Mass and translate molecule."""
        totalMass = 0.0
        XMass = 0.0
        YMass = 0.0
        ZMass = 0.0
        
        for atom in self.atoms:
            rMass = atom.rMass
            totalMass += rMass
            XMass += rMass * atom.x()
            YMass += rMass * atom.y()
            ZMass += rMass * atom.z()
        
        if totalMass > 0:
            XMass /= totalMass
            YMass /= totalMass
            ZMass /= totalMass

        for atom in self.atoms:
            atom.translation(XMass, YMass, ZMass)
        
        # Update bonds after translation (vectors shouldn't change, but good practice)
        self.update_bond_vectors()

    def _build_inertia_tensor(self):
        """Build the moment-of-inertia tensor at the current geometry (shared
        by MIT() and the principal_axes() accessor)."""
        XX = YY = ZZ = 0.0
        XY = XZ = YZ = 0.0
        for atom in self.atoms:
            rMass = atom.rMass
            x, y, z = atom.x(), atom.y(), atom.z()

            XX += rMass * (y**2 + z**2)
            YY += rMass * (x**2 + z**2)
            ZZ += rMass * (x**2 + y**2)
            XY -= rMass * x * y
            XZ -= rMass * x * z
            YZ -= rMass * y * z

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
        heaviest_idx = 0
        max_mass = -1.0
        for i, atom in enumerate(self.atoms):
            if atom.rMass > max_mass:
                max_mass = atom.rMass
                heaviest_idx = i

        h_atom = self.atoms[heaviest_idx]
        h_x = h_atom.x()
        h_y = h_atom.y()
        h_z = h_atom.z()

        new_h_x = h_x * rot[0, 0] + h_y * rot[1, 0] + h_z * rot[2, 0]
        new_h_y = h_x * rot[0, 1] + h_y * rot[1, 1] + h_z * rot[2, 1]
        new_h_z = h_x * rot[0, 2] + h_y * rot[1, 2] + h_z * rot[2, 2]

        if (new_h_x + new_h_y + new_h_z) < 0.0:
            rot = -rot

        for atom in self.atoms:
            x, y, z = atom.x(), atom.y(), atom.z()
            atom.coord.X = x * rot[0, 0] + y * rot[1, 0] + z * rot[2, 0]
            atom.coord.Y = x * rot[0, 1] + y * rot[1, 1] + z * rot[2, 1]
            atom.coord.Z = x * rot[0, 2] + y * rot[1, 2] + z * rot[2, 2]

        self.update_bond_vectors()  # recompute rather than rotate existing vectors

        if modes is not None:
            if rotate_modes:
                rotated_modes = []
                for mode in modes:
                    vecs = mode['vector']  # shape (N, 3)
                    new_vecs = np.zeros_like(vecs)

                    for a in range(self.n):
                        x, y, z = vecs[a][0], vecs[a][1], vecs[a][2]
                        new_vecs[a][0] = x * rot[0, 0] + y * rot[1, 0] + z * rot[2, 0]
                        new_vecs[a][1] = x * rot[0, 1] + y * rot[1, 1] + z * rot[2, 1]
                        new_vecs[a][2] = x * rot[0, 2] + y * rot[1, 2] + z * rot[2, 2]

                    new_mode = mode.copy()
                    new_mode['vector'] = new_vecs
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

        rx_vecs = np.zeros((self.n, 3))
        ry_vecs = np.zeros((self.n, 3))
        rz_vecs = np.zeros((self.n, 3))

        for a in range(self.n):
            x = self.atoms[a].x()
            y = self.atoms[a].y()
            z = self.atoms[a].z()

            rx_vecs[a] = np.array([0.0, -z, y])
            ry_vecs[a] = np.array([z, 0.0, -x])
            rz_vecs[a] = np.array([-y, x, 0.0])

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
        for i in range(self.n):
            self.atoms[i].dispVec = mode_vector[i]
            self.atoms[i].dispLength = sizeVec(mode_vector[i])

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
        Tx, Ty, Tz = 0.0, 0.0, 0.0

        for atom in self.atoms:
            # EPS_DENOM (not the looser EPS_DISP) is the noise floor here:
            # EPS_DISP=1e-8 is too permissive against the ~1e-8-1e-6 noise
            # Gaussian prints for atoms symmetry-required to be exactly zero
            # in degenerate EMIT eigenvectors, which would otherwise be
            # promoted to a full-weight unit-vector contribution.
            if atom.dispLength > EPS_DENOM:
                Tx += atom.dispVec[0] / atom.dispLength
                Ty += atom.dispVec[1] / atom.dispLength
                Tz += atom.dispVec[2] / atom.dispLength
        
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
            total = 0.0
            n_off = 0  # off-axis atom count = N - N_Q
            for atom in self.atoms:
                r = np.array([atom.x(), atom.y(), atom.z()])
                r_perp = r - np.dot(r, Q) * Q
                lr = sizeVec(r_perp)
                if lr <= EPS_DENOM:                # on the axis -> excluded (N_Q)
                    continue
                n_off += 1
                ld = atom.dispLength
                if ld > EPS_DENOM:                 # unit(0):=0 otherwise
                    total += np.dot(np.cross(r_perp, atom.dispVec), Q) / (lr * ld)
            out[key] = total / n_off if n_off > 0 else 0.0
        return out

    def _bond_contributions(self):
        """Per-bond pieces of the V-score (eq:vscore numerator/denominator),
        plus the diagnostic signed relative bond-length change. Returns four
        parallel lists over self.bList:
          terms        : |Δb_AB|^2 * |unit(Δb_AB).b_hat_AB|  (numerator term,
                         magnitude -- feeds s[V_S] via Vscore(), unchanged)
          signed_terms : |Δb_AB|^2 * (unit(Δb_AB).b_hat_AB)  (signed version
                         of the same term, no outer abs() -- feeds the signed
                         per-bond s_AB in score_bonds(); positive = stretching,
                         negative = compressing)
          sqdisps      : |Δb_AB|^2                            (denominator term)
          rel_db       : (|b_AB+Δd_AB| - |b_AB|) / |b_AB|. NOT part of
                         eq:vscore/eq:bondscore -- a diagnostic only
                         (fig:bondscores' x-axis).
        Vscore() and score_bonds() both build on this. s_AB is signed but
        s[V_S] is not: Σ|s_AB| == s[V_S] exactly (not Σ s_AB).
        """
        terms = []
        signed_terms = []
        sqdisps = []
        rel_db = []
        for b in range(self.nBond):
            idx1, idx2 = self.bList[b]
            atom1 = self.atoms[idx1]
            atom2 = self.atoms[idx2]

            bVec_curr = self.bVec[b]  # equilibrium bond vector
            bLength = sizeVec(bVec_curr)

            delDisp = atom2.dispVec - atom1.dispVec
            delDispLength = sizeVec(delDisp)

            dot_val = np.dot(delDisp, bVec_curr)

            terms.append(abs(dot_val) * delDispLength / bLength)
            signed_terms.append(dot_val * delDispLength / bLength)
            sqdisps.append(delDispLength**2)

            if bLength > EPS_NORM:
                rel_db.append((sizeVec(bVec_curr + delDisp) - bLength) / bLength)
            else:
                rel_db.append(0.0)
        return terms, signed_terms, sqdisps, rel_db

    def Vscore(self):
        """Calculates the Vibrational Score s[V_S] (eq:vscore)."""
        terms, _signed_terms, sqdisps, _ = self._bond_contributions()
        modeScr = sum(terms)
        denom = sum(sqdisps)

        if denom > EPS_DENOM:
            return modeScr / denom
        return 0.0

    def score_bonds(self):
        """Per-bond stretch contribution s_AB (eq:bondscore), signed:
        s_AB = |Δb_AB|^2 * (unit(Δb_AB).b_hat_AB) / Σ_bonds |Δb|^2

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
                "i_label": f"{self.atoms[idx1].symbol}{idx1 + 1}",
                "j_label": f"{self.atoms[idx2].symbol}{idx2 + 1}",
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