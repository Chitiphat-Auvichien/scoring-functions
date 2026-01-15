import numpy as np
import math
from .utils import atomicMass

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
        
        # Calculate initial bond vectors (Vector from i to j)
        for (i, j) in self.bList:
            vec = np.array([
                self.atoms[j].x() - self.atoms[i].x(),
                self.atoms[j].y() - self.atoms[i].y(),
                self.atoms[j].z() - self.atoms[i].z()
            ])
            self.bVec.append(vec)

        # Move to Center of Mass (Required for Rscore)
        self.COM()

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