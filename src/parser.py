import re
import os
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

    def parse(self, parse_modes=True):
        self._parse_standard_orientation()
        if parse_modes:
            self._parse_modes()
        self._parse_connectivity()
        
        return {
            "atoms": self.atom_symbols,
            "coords": np.array(self.coordinates),
            "bonds": self.bonds,
            "modes": self.modes
        }

    def _parse_standard_orientation(self):
        # Look for the LAST "Standard orientation" block
        start_indices = [i for i, line in enumerate(self.lines) if "Standard orientation" in line]
        
        if not start_indices:
            start_indices = [i for i, line in enumerate(self.lines) if "Input orientation" in line]
        
        if not start_indices:
            raise ValueError("No orientation block found in log file.")
        
        start_idx = start_indices[-1] + 5 
        self.geom_line = start_idx 
        
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
        """Looks for .com/.gjf file in data/gjf/ matching the log filename."""
        base_name = os.path.splitext(os.path.basename(self.filepath))[0]
        data_dir = os.path.dirname(os.path.dirname(self.filepath))
        gjf_dir = os.path.join(data_dir, "gjf")
        
        input_file = None
        for ext in ['.com', '.gjf']:
            candidate = os.path.join(gjf_dir, base_name + ext)
            if os.path.exists(candidate):
                input_file = candidate
                break
        
        if not input_file: return

        try:
            with open(input_file, 'r') as f:
                lines = [l.strip() for l in f.readlines()]
            
            idx = 0
            found_charge = False
            for i in range(len(lines)):
                parts = lines[i].split()
                if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
                    idx = i + 1
                    found_charge = True
                    break
            
            if not found_charge: return

            while idx < len(lines) and lines[idx]: idx += 1
            idx += 1
            
            while idx < len(lines):
                line = lines[idx]
                if not line: break 
                parts = line.split()
                if len(parts) >= 2 and parts[0].isdigit():
                    center = int(parts[0]) - 1 
                    for k in range(1, len(parts), 2):
                        try:
                            neighbor = int(parts[k]) - 1
                            pair = tuple(sorted((center, neighbor)))
                            if pair not in self.bonds:
                                self.bonds.append(pair)
                        except ValueError: pass
                idx += 1
        except Exception: pass

    def _parse_modes(self):
        # 1. Determine if High Precision (HP) modes are present
        # Gaussian prints "Coord Atom Element:" for HP modes.
        has_hp = any("Coord Atom Element:" in line for line in self.lines)

        # 2. Collect frequency lines
        # Only scan after the geometry we parsed to avoid initial guess freqs
        freq_lines = []
        for i in range(len(self.lines)):
            if "Frequencies --" in self.lines[i]:
                if i > self.geom_line:
                    freq_lines.append(i)
        
        if not freq_lines:
            freq_lines = [i for i, line in enumerate(self.lines) if "Frequencies --" in line]

        # 3. Parse blocks
        for start_idx in freq_lines:
            # Pass the has_hp flag. If True, we ONLY parse HP blocks and skip Standard ones.
            self._parse_block(start_idx, force_hp_only=has_hp)

    def _parse_block(self, line_idx, force_hp_only=False):
        try:
            parts = self.lines[line_idx].split()
            freqs = [float(x) for x in parts[2:]]
        except ValueError: return
            
        num_modes_in_block = len(freqs)
        data_start = -1
        is_hp = False
        
        # Scan ahead to see what kind of data block follows this frequency line
        for i in range(line_idx + 1, min(len(self.lines), line_idx + 30)):
            line = self.lines[i].strip()
            
            # Check for HP header
            if "Coord Atom Element:" in line:
                data_start = i + 1
                is_hp = True
                break
            
            # Check for Standard header
            if "Atom" in line and "AN" in line:
                data_start = i + 1
                is_hp = False
                break
            
            # Fallback HP pattern check (if header missing but data present)
            parts = line.split()
            if len(parts) > 3 and parts[0].isdigit() and parts[1].isdigit() and parts[2].isdigit():
                try:
                    float(parts[3])
                    # It looks like HP data
                    data_start = i
                    is_hp = True
                    break
                except ValueError: pass

        if data_start == -1: return

        # CRITICAL FIX:
        # If HP modes exist in the file (force_hp_only=True), but this specific block
        # is identified as Standard (is_hp=False), we SKIP it.
        # This prevents reading the Standard block when HP is available.
        if force_hp_only and not is_hp:
            return

        current_line = data_start
        if is_hp:
            temp_data = np.zeros((num_modes_in_block, self.natoms, 3))
            count = 0
            while count < 3 * self.natoms and current_line < len(self.lines):
                parts = self.lines[current_line].split()
                if len(parts) < 3: break
                try:
                    c_idx = int(parts[0]) - 1 
                    a_idx = int(parts[1]) - 1
                    if a_idx >= self.natoms: 
                        current_line += 1
                        continue
                    vals = [float(v) for v in parts[3:]]
                    for m in range(min(len(vals), num_modes_in_block)):
                        temp_data[m, a_idx, c_idx] = vals[m]
                    count += 1
                    current_line += 1
                except ValueError: break
            
            for m in range(num_modes_in_block):
                self.modes.append({"frequency": freqs[m], "vector": temp_data[m], "is_emit": False})
        else:
            block_vectors = [[] for _ in range(num_modes_in_block)]
            count = 0
            while count < self.natoms and current_line < len(self.lines):
                parts = self.lines[current_line].split()
                try:
                    raw_vals = [float(v) for v in parts[2:]]
                    for m in range(num_modes_in_block):
                        start_col = m * 3
                        if start_col + 3 <= len(raw_vals):
                            vec = raw_vals[start_col : start_col+3]
                            block_vectors[m].append(vec)
                    count += 1
                    current_line += 1
                except (ValueError, IndexError): break

            for m in range(num_modes_in_block):
                if len(block_vectors[m]) == self.natoms:
                    self.modes.append({"frequency": freqs[m], "vector": np.array(block_vectors[m]), "is_emit": False})

class EMITParser:
    """Parses EMIT modes and Eigenvalues from a text file."""
    def __init__(self, filepath, num_atoms):
        self.filepath = filepath
        self.natoms = num_atoms
        self.modes = []

    def parse(self):
        with open(self.filepath, 'r') as f:
            lines = f.readlines()
        
        start_idx = 0
        for i, line in enumerate(lines):
            if "CART EMIT modes" in line:
                start_idx = i + 1
                break
        
        dim = 3 * self.natoms
        matrix_values = []
        
        i = start_idx
        while i < len(lines):
            line = lines[i].strip()
            if "Eigenvalues" in line: break
            if not line: 
                i += 1
                continue
            parts = line.split()
            for p in parts:
                try: matrix_values.append(float(p))
                except ValueError: pass
            i += 1
            
        eigenvalues = []
        eig_start = -1
        for j in range(len(lines)):
            if "Eigenvalues" in lines[j]:
                eig_start = j
                break
        
        if eig_start != -1:
            for j in range(eig_start, len(lines)):
                clean_line = lines[j].replace("Eigenvalues:", "").strip()
                if not clean_line: continue
                parts = clean_line.split()
                for p in parts:
                    try: eigenvalues.append(float(p))
                    except ValueError: pass

        expected_size = dim * dim
        current_size = len(matrix_values)
        
        if current_size >= expected_size:
            raw_matrix = np.array(matrix_values[:expected_size]).reshape(dim, dim)
            
            for m in range(dim):
                eigenvec = raw_matrix[:, m]
                mode_vec = np.zeros((self.natoms, 3))
                for a in range(self.natoms):
                    mode_vec[a, 0] = eigenvec[3*a]
                    mode_vec[a, 1] = eigenvec[3*a+1]
                    mode_vec[a, 2] = eigenvec[3*a+2]
                
                eig = eigenvalues[m] if m < len(eigenvalues) else 0.0
                
                self.modes.append({
                    "frequency": eig, 
                    "vector": mode_vec,
                    "label": f"EMIT {m+1}",
                    "is_emit": True
                })
        else:
            print(f"Error parsing EMIT file: Expected {expected_size} values, found {current_size}.")
            
        return self.modes

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
                f.write("# No bonds detected. Please add them below.\n")
            else:
                for b in data['bonds']:
                    f.write(f"{b[0]+1} {b[1]+1}\n")
            
            f.write(f"NUM_MODES {len(data['modes'])}\n")
            for i, mode in enumerate(data['modes']):
                label = mode.get('label', f"MODE {i+1}")
                f.write(f"{label}\n")
                
                is_emit = mode.get('is_emit', False)
                header = "EIGEN" if is_emit else "FREQ"
                f.write(f"{header} {mode['frequency']:.6f}\n")
                
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
                    if "NUM_MODES" in line or not line: break
                    if line.startswith('#'):
                        idx += 1
                        continue
                    try:
                        b_parts = line.split()
                        if len(b_parts) >= 2:
                            center = int(b_parts[0]) - 1
                            for k in range(1, len(b_parts), 2):
                                neighbor = int(b_parts[k]) - 1
                                data['bonds'].append((center, neighbor))
                    except ValueError: break
                    idx += 1
            elif parts[0] == "NUM_MODES":
                idx += 1
            elif "MODE" in parts[0] or "EMIT" in parts[0]:
                label = line
                if idx + 3 >= len(content): break
                
                val_line = content[idx+1].strip()
                val_parts = val_line.split()
                try:
                    val = float(val_parts[1])
                except (IndexError, ValueError): val = 0.0
                is_emit = "EIGEN" in val_parts[0]
                
                idx += 3 
                vecs = []
                for _ in range(len(data['atoms'])):
                    vecs.append([float(x) for x in content[idx].split()])
                    idx += 1
                data['modes'].append({"frequency": val, "vector": np.array(vecs), "label": label, "is_emit": is_emit})
            else:
                idx += 1
        return data