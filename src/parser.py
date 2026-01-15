import re
import numpy as np
from .utils import get_symbol

class GaussianParser:
    def __init__(self, filepath):
        self.filepath = filepath
        with open(filepath, 'r') as f:
            self.lines = f.readlines()
        self.natoms = 0
        self.atom_symbols = []
        self.coordinates = []
        self.modes = []
        self.bonds = []

    def parse(self):
        self._parse_standard_orientation()
        self._parse_modes()
        self._parse_connectivity()
        
        # NOTE: Automatic bonding guess removed as requested.
        # If no bonds are found, the user will be prompted in main.py.

        return {
            "atoms": self.atom_symbols,
            "coords": np.array(self.coordinates),
            "bonds": self.bonds,
            "modes": self.modes
        }

    def _parse_standard_orientation(self):
        start_indices = [i for i, line in enumerate(self.lines) if "Standard orientation" in line]
        if not start_indices:
            start_indices = [i for i, line in enumerate(self.lines) if "Input orientation" in line]
        
        if not start_indices:
            raise ValueError("No orientation block found in log file.")
        
        start_idx = start_indices[-1] + 5
        coords = []
        symbols = []
        i = start_idx
        while "-----" not in self.lines[i]:
            parts = self.lines[i].split()
            if len(parts) < 6: break
            try:
                atomic_num = int(parts[1])
                x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                symbols.append(get_symbol(atomic_num))
                coords.append([x, y, z])
            except ValueError: pass
            i += 1
        
        self.natoms = len(coords)
        self.atom_symbols = symbols
        self.coordinates = coords

    def _parse_connectivity(self):
        """Attempts to parse Gaussian 'Geom=Connectivity' blocks."""
        # This is a basic parser for the 1 2 1.0 3 1.0 format
        # It looks for the block usually at the end or explicitly labeled
        # For robustness in this educational tool, we largely rely on manual input
        # if this fails, but we try a simple scan:
        
        # Look for a line that might be connectivity: integers followed by floats
        # This is complex to regex robustly across all Gaussian versions.
        # For now, we leave this empty. If Gaussian log doesn't have it explicitly
        # in a format we strictly know, we default to manual user input.
        pass

    def _parse_modes(self):
        hp_indices = [i for i, line in enumerate(self.lines) if "Coord Atom Element:" in line]
        if hp_indices:
            self._parse_hpmodes(hp_indices[-1])
        else:
            self._parse_standard_modes()

    def _parse_hpmodes(self, start_idx):
        freq_line_idx = -1
        for i in range(start_idx, 0, -1):
            if "Frequencies --" in self.lines[i]:
                freq_line_idx = i
                break
        if freq_line_idx == -1: return

        freqs = self.lines[freq_line_idx].split()[2:]
        num_modes = len(freqs)
        mode_data = np.zeros((num_modes, self.natoms, 3))
        
        i = start_idx + 1
        while i < len(self.lines):
            parts = self.lines[i].split()
            if len(parts) < 3: break
            try:
                coord_idx = int(parts[0]) - 1 
                atom_idx = int(parts[1]) - 1
                values = [float(v) for v in parts[3:]]
                for m in range(min(len(values), num_modes)):
                    mode_data[m, atom_idx, coord_idx] = values[m]
            except ValueError: break
            i += 1
            
        for m in range(num_modes):
            self.modes.append({"frequency": float(freqs[m]), "vector": mode_data[m]})

    def _parse_standard_modes(self):
        freq_lines = [i for i, line in enumerate(self.lines) if "Frequencies --" in line]
        for line_idx in freq_lines:
            freqs = self.lines[line_idx].split()[2:]
            num_modes = len(freqs)
            coord_start = line_idx + 5
            for offset in range(1, 15):
                if "Atom" in self.lines[line_idx + offset]:
                    coord_start = line_idx + offset + 1
                    break
            block_vectors = [[] for _ in range(num_modes)]
            for k in range(self.natoms):
                parts = self.lines[coord_start + k].split()
                try:
                    values = [float(v) for v in parts[2:]]
                    for m in range(num_modes):
                        block_vectors[m].append(values[m*3:(m+1)*3])
                except (ValueError, IndexError): continue
            for m in range(num_modes):
                self.modes.append({"frequency": float(freqs[m]), "vector": np.array(block_vectors[m])})

class IntermediateIO:
    @staticmethod
    def save(data, filename):
        with open(filename, 'w') as f:
            f.write("MOLECULE_DATA\n")
            f.write(f"NATOMS {len(data['atoms'])}\n")
            f.write(f"ATOMS {' '.join(data['atoms'])}\n")
            
            f.write("COORDINATES (Angstrom)\n")
            for c in data['coords']:
                f.write(f"{c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n")
            
            f.write("BONDS (AtomIndex1 AtomIndex2) [1-based]\n")
            if not data['bonds']:
                f.write("# No bonds detected. Please add them below (e.g., '1 2' for bond between atom 1 and 2)\n")
            else:
                for b in data['bonds']:
                    # Convert 0-based (internal) to 1-based (user/chemistry standard)
                    f.write(f"{b[0]+1} {b[1]+1}\n")
            
            f.write(f"NUM_MODES {len(data['modes'])}\n")
            for i, mode in enumerate(data['modes']):
                f.write(f"MODE {i+1}\n")
                f.write(f"FREQ {mode['frequency']:.4f}\n")
                f.write("VECTOR\n")
                for v in mode['vector']:
                    f.write(f"{v[0]:.5f} {v[1]:.5f} {v[2]:.5f}\n")
            f.write("END_MOLECULE_DATA\n")
        print(f"Intermediate file saved to: {filename}")

    @staticmethod
    def load(filename):
        data = {'atoms': [], 'coords': [], 'bonds': [], 'modes': []}
        
        with open(filename, 'r') as f:
            content = f.read().splitlines()
        
        idx = 0
        while idx < len(content):
            line = content[idx].strip()
            parts = line.split()
            
            if not parts or line.startswith('#'): 
                idx += 1
                continue
                
            if parts[0] == "ATOMS":
                data['atoms'] = parts[1:]
                idx += 1
                
            elif parts[0] == "COORDINATES":
                idx += 1
                coords = []
                for _ in range(len(data['atoms'])):
                    coords.append([float(x) for x in content[idx].split()])
                    idx += 1
                data['coords'] = np.array(coords)
                
            elif parts[0] == "BONDS":
                idx += 1
                while idx < len(content):
                    line = content[idx].strip()
                    if "NUM_MODES" in line or not line:
                        break
                    if line.startswith('#'):
                        idx += 1
                        continue
                    
                    try:
                        b_parts = line.split()
                        if len(b_parts) >= 2:
                            # Convert 1-based (user) to 0-based (internal)
                            i1 = int(b_parts[0]) - 1
                            i2 = int(b_parts[1]) - 1
                            data['bonds'].append((i1, i2))
                    except ValueError:
                        break
                    idx += 1
            
            elif parts[0] == "NUM_MODES":
                idx += 1
                
            elif parts[0] == "MODE":
                if idx + 3 >= len(content): break
                freq = float(content[idx+1].split()[1])
                idx += 3 
                vecs = []
                for _ in range(len(data['atoms'])):
                    vecs.append([float(x) for x in content[idx].split()])
                    idx += 1
                data['modes'].append({"frequency": freq, "vector": np.array(vecs)})
            else:
                idx += 1
                
        return data