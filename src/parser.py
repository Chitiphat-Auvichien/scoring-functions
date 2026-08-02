import re
import os
import warnings
import numpy as np
from .utils import get_symbol, get_atomic_number, find_file

# --- Shared Gaussian frequency-block helpers ---------------------------
# Both "HP" (Coord Atom Element:, 3-dash) and "Standard" (Atom AN, 2-dash)
# block variants share fixed line offsets around the "Frequencies" line:
# -1 irreps, 0 Frequencies, +1 Reduced masses (AMU), +2 Force constants
# (mDyne/Angstrom). Used by GaussianParser._parse_block() and by
# IntermediateIO's normal-mode format (structurally an HP block).

_DASH_RE = re.compile(r'^-+$')


def _values_after_dash_marker(line):
    """Split `line` on whitespace and return the floats appearing after the
    first pure-dash token (e.g. '---' or '--') -- robust to both the 3-dash
    HP block style ('Frequencies ---') and the 2-dash Standard style
    ('Frequencies --', 'Red. masses --', 'Frc consts  --')."""
    parts = line.split()
    for i, tok in enumerate(parts):
        if _DASH_RE.match(tok):
            return [float(x) for x in parts[i + 1:]]
    raise ValueError(f"no dash-marker token ('---'/'--') found in line: {line!r}")


def _parse_freq_block_header(lines, line_idx):
    """Parse the fixed-offset header lines around a Gaussian frequency-block
    line at `line_idx`; returns parallel lists (freqs, irreps, mus, ks), one
    entry per mode. mus/ks/irreps fall back to None per-entry if their line
    is absent, malformed, or has a mismatched token count."""
    freqs = _values_after_dash_marker(lines[line_idx])
    n = len(freqs)

    irreps = lines[line_idx - 1].split() if line_idx - 1 >= 0 else []
    if len(irreps) != n:
        irreps = [None] * n

    def _optional_values(offset):
        try:
            vals = _values_after_dash_marker(lines[line_idx + offset])
        except (ValueError, IndexError):
            return [None] * n
        return vals if len(vals) == n else [None] * n

    mus = _optional_values(1)
    ks = _optional_values(2)
    return freqs, irreps, mus, ks


def _parse_hp_mode_vectors(lines, data_start, natoms, num_modes_in_block):
    """Read an HP-style 'Coord Atom Element:' displacement table starting at
    `data_start`: rows of 'coord_idx atom_idx atomic_num val1..valK',
    3*natoms rows total. Returns (vectors, next_line) where vectors has
    shape (num_modes_in_block, natoms, 3) and next_line lets the caller
    resume scanning for the next block."""
    temp_data = np.zeros((num_modes_in_block, natoms, 3))
    current_line = data_start
    count = 0
    while count < 3 * natoms and current_line < len(lines):
        parts = lines[current_line].split()
        if len(parts) < 3:
            break
        try:
            c_idx = int(parts[0]) - 1
            a_idx = int(parts[1]) - 1
            if a_idx >= natoms:
                current_line += 1
                continue
            vals = [float(v) for v in parts[3:]]
            for m in range(min(len(vals), num_modes_in_block)):
                temp_data[m, a_idx, c_idx] = vals[m]
            count += 1
            current_line += 1
        except ValueError:
            break
    return temp_data, current_line


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
            # Fail loud: vibrational modes must number 3N-6 (nonlinear) or 3N-5 (linear).
            nm = len(self.modes)
            expected = {3 * self.natoms - 6, 3 * self.natoms - 5}
            if nm not in expected:
                raise ValueError(
                    f"Parsed {nm} vibrational modes for {self.natoms} atoms "
                    f"({os.path.basename(self.filepath)}); expected 3N-6="
                    f"{3*self.natoms-6} (nonlinear) or 3N-5={3*self.natoms-5} (linear). "
                    "The frequency block was likely mis-parsed.")
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

        input_file = find_file(gjf_dir, base_name, ('.com', '.gjf'))
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
        except (OSError, UnicodeDecodeError, IndexError) as e:
            warnings.warn(f"Failed to parse connectivity from {input_file}: {e}")

    def _parse_modes(self):
        # Gaussian prints "Coord Atom Element:" for High-Precision (HP) modes.
        has_hp = any("Coord Atom Element:" in line for line in self.lines)

        # Only scan after the geometry we parsed to avoid initial-guess freqs.
        freq_lines = []
        for i in range(len(self.lines)):
            if "Frequencies --" in self.lines[i]:
                if i > self.geom_line:
                    freq_lines.append(i)

        if not freq_lines:
            freq_lines = [i for i, line in enumerate(self.lines) if "Frequencies --" in line]

        for start_idx in freq_lines:
            self._parse_block(start_idx, force_hp_only=has_hp)

    def _parse_block(self, line_idx, force_hp_only=False):
        try:
            freqs, irreps, mus, ks = _parse_freq_block_header(self.lines, line_idx)
        except (ValueError, IndexError):
            return

        num_modes_in_block = len(freqs)
        data_start = -1
        is_hp = False

        # Scan ahead to see what kind of data block follows this frequency line.
        for i in range(line_idx + 1, min(len(self.lines), line_idx + 30)):
            line = self.lines[i].strip()

            if "Coord Atom Element:" in line:
                data_start = i + 1
                is_hp = True
                break

            if "Atom" in line and "AN" in line:
                data_start = i + 1
                is_hp = False
                break

            # Fallback HP pattern check (if header missing but data present).
            parts = line.split()
            if len(parts) > 3 and parts[0].isdigit() and parts[1].isdigit() and parts[2].isdigit():
                try:
                    float(parts[3])
                    data_start = i
                    is_hp = True
                    break
                except ValueError: pass

        if data_start == -1: return

        # CRITICAL: when HP blocks exist in the file, a Standard block for the
        # same frequencies must be skipped, or modes get double-counted / read
        # at the wrong precision.
        if force_hp_only and not is_hp:
            return

        current_line = data_start
        if is_hp:
            temp_data, current_line = _parse_hp_mode_vectors(
                self.lines, data_start, self.natoms, num_modes_in_block)

            for m in range(num_modes_in_block):
                self.modes.append({
                    "frequency": freqs[m], "vector": temp_data[m], "is_emit": False,
                    "reduced_mass": mus[m], "force_constant": ks[m], "irrep": irreps[m],
                })
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
                    self.modes.append({
                        "frequency": freqs[m], "vector": np.array(block_vectors[m]), "is_emit": False,
                        "reduced_mass": mus[m], "force_constant": ks[m], "irrep": irreps[m],
                    })

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
            raise ValueError(
                f"Error parsing EMIT file {os.path.basename(self.filepath)}: expected "
                f"{expected_size} matrix values (3N x 3N for N={self.natoms}), found {current_size}.")

        # Fail loud: EMIT modes must number exactly 3N.
        if len(self.modes) != 3 * self.natoms:
            raise ValueError(
                f"Parsed {len(self.modes)} EMIT modes for {self.natoms} atoms; "
                f"expected exactly 3N={3*self.natoms}.")
        return self.modes


class IntermediateIO:
    """Save/load the editable intermediate checkpoint file between parsing
    and scoring (see CLAUDE.md / main.load_inputs()).

    Two structurally different formats live behind this interface,
    auto-detected on load() by content (not filename):
      - EMIT intermediates (all modes is_emit=True): the original simple
        MOLECULE_DATA/NATOMS/ATOMS/COORDINATES/BONDS/NUM_MODES/MODE format
        (no mu/k/irrep concept). See _save_emit/_load_emit.
      - Normal-mode intermediates (all modes is_emit=False): a Gaussian-direct
        format -- header 'NATOMS LINEAR NAME', a gjf-style connectivity block,
        a Standard-orientation coordinate block, then raw Gaussian HP
        frequency blocks (5 modes/block, carrying irrep/mu/k). See
        _save_normal/_load_normal, which reuse
        _parse_freq_block_header/_parse_hp_mode_vectors since this format is
        structurally identical to a real HP frequency block.
    """

    @staticmethod
    def save(data, filename):
        modes = data['modes']
        is_emit = bool(modes) and all(m.get('is_emit', False) for m in modes)
        if is_emit:
            IntermediateIO._save_emit(data, filename)
        else:
            IntermediateIO._save_normal(data, filename)

    @staticmethod
    def load(filename):
        with open(filename, 'r') as f:
            content = f.read().splitlines()
        if content and content[0].strip() == "MOLECULE_DATA":
            return IntermediateIO._load_emit(content)
        return IntermediateIO._load_normal(content, filename)

    # --- Normal-mode format (Gaussian-direct) ---------------

    @staticmethod
    def _save_normal(data, filename):
        atoms = data['atoms']
        coords = data['coords']
        bonds = data['bonds']
        modes = data['modes']
        natoms = len(atoms)
        nmodes = len(modes)

        # linear = (nmodes == 3N-5); informational round-trip metadata only
        # -- ModeScorer/classifier detect linearity independently from the
        # moment-of-inertia tensor at scoring time.
        linear = 1 if nmodes == 3 * natoms - 5 else 0

        mol_name = os.path.splitext(os.path.basename(filename))[0]
        for suffix in ("_normal_data", "_data"):
            if mol_name.endswith(suffix):
                mol_name = mol_name[: -len(suffix)]
                break

        with open(filename, 'w') as f:
            f.write(f"{natoms} {linear} {mol_name}\n")

            # Bond block: one line per atom, gjf geom=connectivity style
            # 'atom neighbor1 order1 ...' (order defaults to 1.0, real bond
            # order isn't retained). Only higher-index neighbors are listed
            # per atom, so each bond appears once (matches data/gjf/*.com).
            adjacency = {i: [] for i in range(natoms)}
            for (i, j) in bonds:
                lo, hi = (i, j) if i < j else (j, i)
                if hi not in adjacency[lo]:
                    adjacency[lo].append(hi)
            for i in range(natoms):
                neighbors = sorted(adjacency[i])
                parts = [str(i + 1)] + [f"{j + 1} 1.0" for j in neighbors]
                f.write(" " + " ".join(parts) + "\n")

            # Coordinate block: Standard-orientation style.
            for i in range(natoms):
                atomic_num = get_atomic_number(atoms[i]) or 0
                x, y, z = coords[i]
                f.write(f"{i + 1:>7d}{atomic_num:>11d}{0:>12d}"
                        f"{x:>16.6f}{y:>12.6f}{z:>12.6f}\n")

            # Frequency blocks, 5 modes per block (Gaussian's own convention).
            BLOCK = 5
            LABEL_W = 23
            VAL_W = 10
            for start in range(0, nmodes, BLOCK):
                block_modes = modes[start:start + BLOCK]
                n_blk = len(block_modes)

                idx_row = "".join(f"{start + k + 1:>{VAL_W}d}" for k in range(n_blk))
                f.write(" " * LABEL_W + idx_row + "\n")

                irrep_row = "".join(
                    f"{(m.get('irrep') or 'NA'):>{VAL_W}s}" for m in block_modes)
                f.write(" " * LABEL_W + irrep_row + "\n")

                def _num_row(label, vals):
                    line = f"{label:>{LABEL_W}}"
                    for v in vals:
                        line += f"{(0.0 if v is None else v):>{VAL_W}.4f}"
                    return line + "\n"

                f.write(_num_row("Frequencies ---", [m['frequency'] for m in block_modes]))
                f.write(_num_row("Reduced masses ---", [m.get('reduced_mass') for m in block_modes]))
                f.write(_num_row("Force constants ---", [m.get('force_constant') for m in block_modes]))
                f.write(_num_row("IR Intensities ---", [0.0] * n_blk))
                f.write(" Coord Atom Element:\n")

                for a in range(natoms):
                    atomic_num = get_atomic_number(atoms[a]) or 0
                    for c in range(3):
                        line = f"{c + 1:>4d}{a + 1:>6d}{atomic_num:>6d}"
                        for m in block_modes:
                            line += f"{m['vector'][a][c]:>{VAL_W}.5f}"
                        f.write(line + "\n")
        print(f"Intermediate file saved to: {filename}")

    @staticmethod
    def _load_normal(content, filename):
        if not content:
            raise ValueError(f"Empty intermediate file: {filename}")

        header_parts = content[0].split()
        if len(header_parts) < 2:
            raise ValueError(
                f"Malformed intermediate file '{filename}': header line "
                f"{content[0]!r} does not match 'NATOMS LINEAR NAME'.")
        natoms = int(header_parts[0])
        linear_flag = int(header_parts[1])
        idx = 1

        # Bond block: exactly `natoms` lines, one per atom, positionally (no
        # header/footer). Each line is 'atom [neighbor [order] ...]'. A line
        # with exactly one trailing token is a bare 'atom neighbor' pair
        # (order defaults to 1.0) -- this is the manual bond-repair format
        # users hand-edit in (CLAUDE.md: "add '1 2'-style bond lines");
        # anything else is read as the full 'neighbor order neighbor order...'
        # form save() writes.
        bonds = []
        for _ in range(natoms):
            parts = content[idx].split()
            idx += 1
            if not parts:
                continue
            center = int(parts[0]) - 1
            rest = parts[1:]
            if len(rest) == 1:
                neighbor = int(rest[0]) - 1
                pair = tuple(sorted((center, neighbor)))
                if pair not in bonds:
                    bonds.append(pair)
            else:
                for k in range(0, len(rest) - 1, 2):
                    neighbor = int(rest[k]) - 1
                    pair = tuple(sorted((center, neighbor)))
                    if pair not in bonds:
                        bonds.append(pair)

        # Coordinate block: exactly `natoms` lines, Standard-orientation
        # style 'Center# AtomicNum AtomicType X Y Z'.
        atoms = []
        coords = []
        for _ in range(natoms):
            parts = content[idx].split()
            idx += 1
            atomic_num = int(parts[1])
            atoms.append(get_symbol(atomic_num))
            coords.append([float(parts[3]), float(parts[4]), float(parts[5])])
        coords = np.array(coords)

        # Frequency blocks: structurally identical to a real Gaussian HP
        # block -- reuse the same shared helpers GaussianParser uses.
        modes = []
        while idx < len(content):
            line = content[idx]
            if "Frequencies ---" in line or "Frequencies --" in line:
                freqs, irreps, mus, ks = _parse_freq_block_header(content, idx)
                n_blk = len(freqs)

                data_start = None
                for j in range(idx + 1, min(len(content), idx + 10)):
                    if "Coord Atom Element:" in content[j]:
                        data_start = j + 1
                        break
                if data_start is None:
                    raise ValueError(
                        f"Malformed intermediate file '{filename}': no "
                        f"'Coord Atom Element:' header after line {idx}.")

                temp_data, next_line = _parse_hp_mode_vectors(content, data_start, natoms, n_blk)
                for m in range(n_blk):
                    modes.append({
                        "frequency": freqs[m], "vector": temp_data[m],
                        "reduced_mass": mus[m], "force_constant": ks[m],
                        "irrep": irreps[m], "is_emit": False,
                    })
                idx = next_line
            else:
                idx += 1

        nm = len(modes)
        expected = {3 * natoms - 6, 3 * natoms - 5}
        if nm not in expected:
            raise ValueError(
                f"Loaded {nm} vibrational modes for {natoms} atoms from "
                f"'{filename}'; expected 3N-6={3*natoms-6} (nonlinear) or "
                f"3N-5={3*natoms-5} (linear). The intermediate file's "
                "frequency block was likely mis-parsed or hand-edited "
                "incorrectly.")
        is_linear_actual = (nm == 3 * natoms - 5)
        if bool(linear_flag) != is_linear_actual:
            raise ValueError(
                f"Intermediate file '{filename}' header says LINEAR="
                f"{linear_flag} but {nm} modes for {natoms} atoms implies "
                f"{'linear' if is_linear_actual else 'nonlinear'} -- "
                "fail loud rather than trust a stale/hand-edited header flag.")

        return {"atoms": atoms, "coords": coords, "bonds": bonds, "modes": modes}

    # --- EMIT format (original, unchanged -- out of scope) ---------------

    @staticmethod
    def _save_emit(data, filename):
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
    def _load_emit(content):
        data = {'atoms': [], 'coords': [], 'bonds': [], 'modes': []}
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
