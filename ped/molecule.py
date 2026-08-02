"""Resolve a molecule name to the on-disk input files build_veda_fmt.py needs.

A molecule's inputs live under <repo_root>/data/ (the parent scoring-functions
repo's shared data directory, NOT this ped/ folder):
  - logs/<name>.log   -- required. The sole source of the .fmt file's verbatim
                          geometry and frequency-block text, for BOTH Framework
                          1 (reconstruct) and Framework 2 (fchk).
  - gjf/<name>.com or gjf/<name>.gjf -- optional, connectivity only (not
                          currently used by the .fmt-building pipeline, but
                          resolved here for parity with src/parser.py's own
                          lookup and for callers that want it).
  - fchk/<name>.fchk or logs/<name>.fchk -- optional, only needed by
                          Framework 2. See ped/hessian_fchk.py for how to
                          produce one.
"""
import os
from dataclasses import dataclass
from typing import Optional

# Dual-extension search order for connectivity files -- mirrors
# src/parser.py's GaussianParser._parse_connectivity(), which tries '.com'
# before '.gjf'.
_GJF_EXTENSIONS = ('.com', '.gjf')


@dataclass(frozen=True)
class MoleculePaths:
    name: str
    log_path: str
    gjf_path: Optional[str]
    fchk_path: Optional[str]


def resolve_molecule(name, repo_root, fchk_path_override=None):
    """Resolve `name` (e.g. "C6H6") to a MoleculePaths under `repo_root`
    (the scoring-functions repo root, i.e. the directory containing data/,
    ped/, src/).

    Raises FileNotFoundError if data/logs/<name>.log is missing (it is
    required), or if `fchk_path_override` is given but does not exist on
    disk (an explicit user-provided path that's wrong should fail loud, not
    silently fall back to the default search order).
    """
    data_dir = os.path.join(repo_root, 'data')

    log_path = os.path.join(data_dir, 'logs', f'{name}.log')
    if not os.path.isfile(log_path):
        raise FileNotFoundError(
            f"No Gaussian log found for molecule {name!r} at {log_path} -- "
            "this file is required (it is the sole source of the .fmt "
            "file's geometry/frequency text blocks for both frameworks).")

    gjf_dir = os.path.join(data_dir, 'gjf')
    gjf_path = None
    for ext in _GJF_EXTENSIONS:
        candidate = os.path.join(gjf_dir, f'{name}{ext}')
        if os.path.isfile(candidate):
            gjf_path = candidate
            break

    if fchk_path_override is not None:
        if not os.path.isfile(fchk_path_override):
            raise FileNotFoundError(
                f"--fchk-path override {fchk_path_override!r} does not exist.")
        fchk_path = fchk_path_override
    else:
        fchk_path = None
        for candidate in (
            os.path.join(data_dir, 'fchk', f'{name}.fchk'),
            os.path.join(data_dir, 'logs', f'{name}.fchk'),
        ):
            if os.path.isfile(candidate):
                fchk_path = candidate
                break

    return MoleculePaths(name=name, log_path=log_path, gjf_path=gjf_path,
                          fchk_path=fchk_path)
