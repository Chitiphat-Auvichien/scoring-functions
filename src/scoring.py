import numpy as np
import math
from .utils import atomicMass

# --- Helper Classes to mimic atom.py structure ---

def sizeVec(v):
    return math.sqrt(np.dot(v, v))

def normalize(v):
    norm = sizeVec(v)
    if norm < 1e-9:
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

    def MIT(self, modes=None, rotate_modes=True):
        """
        Rotates the molecule and displacement vectors into the basis of principal axes of rotation.
        Integrated from atom.py.
        """
        # 1. Compute Moment of Inertia Tensor
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

        tensor = np.array([
            [XX, XY, XZ],
            [XY, YY, YZ],
            [XZ, YZ, ZZ]
        ])

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

        return {
            "T": self.Tscore(),
            "R": self.Rscore(),
            "V": self.Vscore()
        }

    def Tscore(self):
        """Calculates Translational Scores (Tx, Ty, Tz). Adapted from atom.py."""
        n = self.n
        Tx, Ty, Tz = 0.0, 0.0, 0.0
        
        for atom in self.atoms:
            if atom.dispLength > 1.0E-6:
                Tx += atom.dispVec[0] / atom.dispLength
                Ty += atom.dispVec[1] / atom.dispLength
                Tz += atom.dispVec[2] / atom.dispLength
        
        return {
            'x': Tx * (1.0/float(n)),
            'y': Ty * (1.0/float(n)),
            'z': Tz * (1.0/float(n))
        }

    def Rscore(self):
        """Calculates Rotational Scores (Rx, Ry, Rz). Adapted from atom.py."""
        tolerance = 1.0E-6
        Rx, Ry, Rz = 0.0, 0.0, 0.0
        Nx, Ny, Nz = 0, 0, 0

        for atom in self.atoms:
            x = atom.x()
            y = atom.y()
            z = atom.z()

            # Rx component
            radiusX = np.array([0.0, y, z])
            lengthX = sizeVec(radiusX)
            if lengthX > tolerance:
                Nx += 1
                if atom.dispLength > tolerance:
                    term = (radiusX[1]*atom.dispVec[2] - radiusX[2]*atom.dispVec[1])
                    Rx += term / (lengthX * atom.dispLength)

            # Ry component
            radiusY = np.array([x, 0.0, z])
            lengthY = sizeVec(radiusY)
            if lengthY > tolerance:
                Ny += 1
                if atom.dispLength > tolerance:
                    term = (radiusY[2]*atom.dispVec[0] - radiusY[0]*atom.dispVec[2])
                    Ry += term / (lengthY * atom.dispLength)

            # Rz component
            radiusZ = np.array([x, y, 0.0])
            lengthZ = sizeVec(radiusZ)
            if lengthZ > tolerance:
                Nz += 1
                if atom.dispLength > tolerance:
                    term = (radiusZ[0]*atom.dispVec[1] - radiusZ[1]*atom.dispVec[0])
                    Rz += term / (lengthZ * atom.dispLength)

        return {
            'x': Rx * (1.0/float(Nx)) if Nx > 0 else 0.0,
            'y': Ry * (1.0/float(Ny)) if Ny > 0 else 0.0,
            'z': Rz * (1.0/float(Nz)) if Nz > 0 else 0.0
        }

    def Vscore(self):
        """Calculates Vibrational Score (V). Adapted from atom.py."""
        modeScr = 0.0
        denom = 0.0
        
        for b in range(self.nBond):
            idx1, idx2 = self.bList[b]
            atom1 = self.atoms[idx1]
            atom2 = self.atoms[idx2]
            
            # Current bond vector
            bVec_curr = self.bVec[b]
            bLength = sizeVec(bVec_curr)

            # Difference in displacement
            delDisp = atom2.dispVec - atom1.dispVec
            delDispLength = sizeVec(delDisp)

            # | (d2-d1) . bondVec | * |d2-d1| / |bondVec|
            dot_val = np.dot(delDisp, bVec_curr)
            
            term = abs(dot_val) * delDispLength / bLength
            modeScr += term
            denom += delDispLength**2
        
        if denom > 1.0E-6:
            return modeScr / denom
        return 0.0