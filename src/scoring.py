import numpy as np
import math
from .utils import atomicMass

# --- Centralized numerical constants (JCC spec conventions) ---
# Displacement cutoff ε_disp from the spec: atoms with |d_A| <= EPS_DISP are
# treated as zero-motion (unit(0):=0). Spec value is 1e-8.
EPS_DISP = 1e-8
# Vector-normalization guard (unit(0):=0); used when normalizing ideal T/R
# basis vectors and arbitrary direction vectors.
EPS_NORM = 1e-9
# Denominator guard for the V-score ratio (Σ|Δb|² near zero -> score 0).
EPS_DENOM = 1e-6
# Relative tolerance for grouping degenerate principal moments of inertia
# into axis blocks (symmetric/spherical tops). Provisional per the plan.
DEGEN_TOL = 1e-3
# Tolerance for the range-invariant score asserts.
RANGE_TOL = 1e-6

# --- Helper Classes to mimic atom.py structure ---

def sizeVec(v):
    return math.sqrt(np.dot(v, v))

def normalize(v):
    norm = sizeVec(v)
    if norm < EPS_NORM:
        return v
    return v / norm

class Coordinate:
    def __init__(self, x, y, z):
        self.X = x
        self.Y = y
        self.Z = z

class Atom:
    def __init__(self, element, x, y, z):
        self.symbol = element
        # Retrieve mass using lowercase symbol key
        self.rMass = float(atomicMass.get(self.symbol.lower(), 1.0))
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
        """
        Initialize with data from parser.
        atom_symbols: list of strings ['O', 'H', 'H']
        coords: list of lists [[x,y,z], ...]
        bonds: list of tuples [(0,1), (0,2)]
        """
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
        """Build the moment-of-inertia tensor at the current geometry.

        Factored out of MIT so the classifier accessors (principal_axes /
        axis_blocks) share one implementation instead of recomputing it.
        """
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
        """Principal moments and axes of the inertia tensor at the current geometry.

        Returns
        -------
        moments : np.ndarray, shape (3,)
            Principal moments of inertia in ascending order (eigenvalues).
        axes : np.ndarray, shape (3, 3)
            Principal-axis directions as columns (eigenvectors of the tensor).

        Note: moments are rotation-invariant, so this is consistent whether
        called before or after MIT(); after MIT the tensor is diagonal and the
        axes reduce to (a permutation/sign of) the identity.
        """
        tensor = self._build_inertia_tensor()
        moments, axes = np.linalg.eigh(tensor)
        return moments, axes

    def axis_blocks(self, rel_tol=DEGEN_TOL):
        """Group principal axes into degeneracy blocks by their moments.

        Axes whose principal moments are equal within a relative tolerance are
        collected into the same block (symmetric/spherical tops), where per-axis
        rotation assignment is ill-defined and the classifier must assign the
        block collectively.

        Returns a list of blocks, each a list of axis indices (0,1,2) referring
        to the ascending-moment ordering of principal_axes().
        """
        moments, _ = self.principal_axes()
        order = list(np.argsort(moments))
        scale = max(float(np.max(np.abs(moments))), 1e-12)
        blocks = []
        current = [order[0]]
        for k in range(1, len(order)):
            prev = moments[order[k - 1]]
            cur = moments[order[k]]
            if abs(cur - prev) <= rel_tol * scale:
                current.append(order[k])
            else:
                blocks.append(current)
                current = [order[k]]
        blocks.append(current)
        return blocks

    def MIT(self, modes=None, rotate_modes=True):
        """
        Rotates the molecule and displacement vectors into the basis of principal axes of rotation.
        Integrated from atom.py.
        """
        # 1. Compute Moment of Inertia Tensor
        tensor = self._build_inertia_tensor()

        # 2. Diagonalize (Principal Axes)
        # eigh returns eigenvalues and eigenvectors (columns of rot)
        eigVal, rot = np.linalg.eigh(tensor)
        
        # 3. Handle coordinate orientation (Heaviest atom check from atom.py)
        # Find heaviest atom
        heaviest_idx = 0
        max_mass = -1.0
        for i, atom in enumerate(self.atoms):
            if atom.rMass > max_mass:
                max_mass = atom.rMass
                heaviest_idx = i
        
        # Project heaviest atom coords onto new axes to check sign
        h_atom = self.atoms[heaviest_idx]
        h_x = h_atom.x()
        h_y = h_atom.y()
        h_z = h_atom.z()
        
        # Calculate new coordinates of heaviest atom temporarily
        new_h_x = h_x * rot[0, 0] + h_y * rot[1, 0] + h_z * rot[2, 0]
        new_h_y = h_x * rot[0, 1] + h_y * rot[1, 1] + h_z * rot[2, 1]
        new_h_z = h_x * rot[0, 2] + h_y * rot[1, 2] + h_z * rot[2, 2]
        
        if (new_h_x + new_h_y + new_h_z) < 0.0:
            # Invert rotation matrix (as per atom.py logic)
            rot = -rot

        # 4. Rotate Atoms
        for atom in self.atoms:
            x, y, z = atom.x(), atom.y(), atom.z()
            atom.coord.X = x * rot[0, 0] + y * rot[1, 0] + z * rot[2, 0]
            atom.coord.Y = x * rot[0, 1] + y * rot[1, 1] + z * rot[2, 1]
            atom.coord.Z = x * rot[0, 2] + y * rot[1, 2] + z * rot[2, 2]

        # 5. Rotate Bond Vectors
        # Recalculate is safer/easier than rotating existing vectors
        self.update_bond_vectors()

        # 6. Rotate Displacement Vectors (Modes)
        if modes is not None:
            if rotate_modes:
                rotated_modes = []
                for mode in modes:
                    # Mode vector shape: (N_atoms, 3)
                    vecs = mode['vector'] # shape (N, 3)
                    new_vecs = np.zeros_like(vecs)
                    
                    for a in range(self.n):
                        x, y, z = vecs[a][0], vecs[a][1], vecs[a][2]
                        # Apply same rotation
                        new_vecs[a][0] = x * rot[0, 0] + y * rot[1, 0] + z * rot[2, 0]
                        new_vecs[a][1] = x * rot[0, 1] + y * rot[1, 1] + z * rot[2, 1]
                        new_vecs[a][2] = x * rot[0, 2] + y * rot[1, 2] + z * rot[2, 2]
                    
                    # Store back
                    new_mode = mode.copy()
                    new_mode['vector'] = new_vecs
                    rotated_modes.append(new_mode)
                
                return rotated_modes
            else:
                return modes
        return None

    def construct_T(self):
        """
        Constructs 3 translational modes (Tx, Ty, Tz).
        Returns a list of 3 mode dictionaries.
        """
        modes = []
        labels = ['Tx', 'Ty', 'Tz']
        
        # Create vectors for x, y, z translation
        for i in range(3):
            # Shape (N, 3)
            vec = np.zeros((self.n, 3))
            
            # Set the i-th component to 1.0 for all atoms
            vec[:, i] = 1.0
            
            # Normalize the entire 3N vector
            # Flatten, calc norm, divide
            flat_norm = np.linalg.norm(vec)
            if flat_norm > 1e-9:
                vec = vec / flat_norm
            
            modes.append({
                "frequency": 0.0, # Placeholder
                "vector": vec,
                "label": labels[i] # Special tag
            })
        return modes

    def construct_R(self):
        """
        Constructs 3 rotational modes (Rx, Ry, Rz).
        Returns a list of 3 mode dictionaries.
        """
        self.COM() # Ensure we are at COM
        modes = []
        labels = ['Rx', 'Ry', 'Rz']
        
        # Arrays to accumulate displacements
        # Rx: cross(x-axis, r) -> vector along tangent
        # Tangent directions for rotation around axes:
        # Rx: (0, -z, y)
        # Ry: (z, 0, -x)
        # Rz: (-y, x, 0)
        
        rx_vecs = np.zeros((self.n, 3))
        ry_vecs = np.zeros((self.n, 3))
        rz_vecs = np.zeros((self.n, 3))
        
        for a in range(self.n):
            x = self.atoms[a].x()
            y = self.atoms[a].y()
            z = self.atoms[a].z()
            
            # Rx
            rx_vecs[a] = np.array([0.0, -z, y])
            # Ry
            ry_vecs[a] = np.array([z, 0.0, -x])
            # Rz
            rz_vecs[a] = np.array([-y, x, 0.0])
            
        # Normalize
        for vecs, lbl in zip([rx_vecs, ry_vecs, rz_vecs], labels):
            flat_norm = np.linalg.norm(vecs)
            if flat_norm > 1e-9:
                vecs = vecs / flat_norm
            
            modes.append({
                "frequency": 0.0,
                "vector": vecs,
                "label": lbl
            })
            
        return modes

    def calculate_scores(self, mode_vector):
        """
        Load a specific mode's displacement vector and calculate all scores.
        mode_vector: np.array of shape (N_atoms, 3)
        """
        # Load displacements into Atom objects
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
        for axis, val in scores["T"].items():
            assert -1.0 - RANGE_TOL <= val <= 1.0 + RANGE_TOL, \
                f"s[T_{axis}]={val} out of [-1,1]"
        for axis, val in scores["R"].items():
            assert -1.0 - RANGE_TOL <= val <= 1.0 + RANGE_TOL, \
                f"s[R_{axis}]={val} out of [-1,1]"
        vs = scores["V"]
        assert -RANGE_TOL <= vs <= 1.0 + RANGE_TOL, f"s[V_S]={vs} out of [0,1]"

    def Tscore(self):
        """Calculates Translational Scores (Tx, Ty, Tz). Adapted from atom.py."""
        n = self.n
        Tx, Ty, Tz = 0.0, 0.0, 0.0

        for atom in self.atoms:
            # Use EPS_DENOM (not the looser EPS_DISP) as the noise floor here,
            # matching Rscore/Vscore: EPS_DISP=1e-8 is too permissive relative
            # to the ~1e-8-1e-6 numerical noise Gaussian prints for atoms that
            # are symmetry-required to be exactly zero in degenerate EMIT
            # eigenvectors, which would otherwise be promoted to a full-weight
            # unit-vector contribution.
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

        Per atom and axis Q, normalize the radius vector r_perp = r-(r.Qhat)Qhat
        and the displacement d SEPARATELY, cross them, and take the Q-component:
            s[R_Q] = (1/(N-N_Q)) sum_offaxis (unit(r_perp) x unit(d)) . Qhat
                   = (1/(N-N_Q)) sum_offaxis (r_perp x d)_Q / (|r_perp| |d|).
        |unit(r_perp) x unit(d)| = sin(phi), phi = angle(r_perp, d): it is 1 only
        for a purely tangential (ideal-rotation) displacement and is reduced as d
        tilts toward radial, so non-rotational in-plane motion is down-weighted
        (unlike normalizing by |omega|=|r_perp x d|, which would discard sin phi).
        N_Q = atoms on the Q-axis (|r_perp| ~ 0), excluded; an atom with |d| ~ 0
        contributes a zero unit vector. A linear molecule's axis returns 0 (n_R=2).
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
        """Per-bond pieces of the V-score (eq:vscore numerator and denominator),
        plus the diagnostic signed relative bond-length change.

        Returns three parallel lists over self.bList:
          terms   : |(d_B - d_A) . b_AB| * |d_B - d_A| / |b_AB|
                    (== |Δb_AB|² * |unit(Δb_AB) · b̂_AB|, the numerator term)
          sqdisps : |d_B - d_A|²                         (the denominator term)
          rel_db  : (|b_AB + Δd_AB| - |b_AB|) / |b_AB|    (signed relative
                    bond-length change; Δd_AB = d_B - d_A, the SAME
                    displacement-difference convention as terms/sqdisps above.
                    This is NOT part of eq:vscore/eq:bondscore -- it is a
                    diagnostic quantity (fig:bondscores' x-axis) originally
                    reverse-engineered from the Excel workbook's own
                    "(|b'|-|b|)/|b|" per-bond column (see
                    src/excel_ingest.py's module docstring); computed here
                    from geometry + mode displacement using the engine's own
                    displacement convention, not read from the spreadsheet.)
        Both Vscore() and score_bonds() build on this so the per-bond s_AB
        sum back to s[V_S] exactly.
        """
        terms = []
        sqdisps = []
        rel_db = []
        for b in range(self.nBond):
            idx1, idx2 = self.bList[b]
            atom1 = self.atoms[idx1]
            atom2 = self.atoms[idx2]

            # Current (equilibrium) bond vector
            bVec_curr = self.bVec[b]
            bLength = sizeVec(bVec_curr)

            # Difference in displacement
            delDisp = atom2.dispVec - atom1.dispVec
            delDispLength = sizeVec(delDisp)

            # | (d2-d1) . bondVec | * |d2-d1| / |bondVec|
            dot_val = np.dot(delDisp, bVec_curr)

            terms.append(abs(dot_val) * delDispLength / bLength)
            sqdisps.append(delDispLength**2)

            if bLength > EPS_NORM:
                rel_db.append((sizeVec(bVec_curr + delDisp) - bLength) / bLength)
            else:
                rel_db.append(0.0)
        return terms, sqdisps, rel_db

    def Vscore(self):
        """Calculates Vibrational Score (V). Adapted from atom.py."""
        terms, sqdisps, _ = self._bond_contributions()
        modeScr = sum(terms)
        denom = sum(sqdisps)

        if denom > EPS_DENOM:
            return modeScr / denom
        return 0.0

    def score_bonds(self):
        """Per-bond stretch contribution s_AB (eq:bondscore), summing to s[V_S],
        plus the diagnostic signed relative bond-length change rel_db (see
        _bond_contributions()).

        s_AB = |Δb_AB|² * |unit(Δb_AB) · b̂_AB| / Σ_bonds |Δb|²  (global denominator)

        Returns a list of dicts {'i', 'j', 's_AB', 'rel_db', 'i_label',
        'j_label'} aligned with self.bList. 'i_label'/'j_label' are
        human-readable atom tags ("<symbol><1-based index>", e.g. "C1", "H7")
        for building interpretable bond identifiers like "C1-C2" instead of
        bare index pairs like "1-2". Asserts Σ s_AB == s[V_S] to 1e-6. Call
        calculate_scores()/load displacements first so the atom dispVecs are
        populated.
        """
        terms, sqdisps, rel_db = self._bond_contributions()
        denom = sum(sqdisps)

        bonds = []
        for b in range(self.nBond):
            idx1, idx2 = self.bList[b]
            s_AB = terms[b] / denom if denom > EPS_DENOM else 0.0
            bonds.append({
                "i": idx1, "j": idx2, "s_AB": s_AB, "rel_db": rel_db[b],
                "i_label": f"{self.atoms[idx1].symbol}{idx1 + 1}",
                "j_label": f"{self.atoms[idx2].symbol}{idx2 + 1}",
            })

        total = sum(bd["s_AB"] for bd in bonds)
        vs = self.Vscore()
        assert abs(total - vs) <= 1e-6, \
            f"Sum of per-bond s_AB ({total}) != s[V_S] ({vs})"
        return bonds