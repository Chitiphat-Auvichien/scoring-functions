"""Empirical CPU-time comparison: Gaussian's freq-only job step vs. our
classification algorithm, for every molecule with a log in data/logs/.

Replaces the manuscript's theoretical Big-O "Computational cost" argument
with a measured one: for each molecule (N = atom count), compare

  (a) gaussian_freq_cpu_s -- the CPU time Gaussian itself reports for the
      frequency calculation alone (NOT the geometry optimization that
      precedes it), and
  (b) classifier_cpu_s    -- the CPU time build_scorer_and_final() +
      classify_all_modes() (Algorithm 1, Steps 1-4) takes on the SAME
      molecule, once inputs are already parsed.

Log structure (verified across all 69 data/logs/*.log -- see docstring of
_gaussian_freq_cpu_seconds() below): every log's route line is
"opt freq=... <method>/<basis> [geom=connectivity]" (job step 1, geometry
optimization), followed by a Link1 restart with route
"#N Geom=AllCheck ... Freq" (job step 2, a freq-only single point on the
optimized/checkpointed geometry). Each job step ends with its own
"Job cpu time:  X days  Y hours  Z minutes  W seconds." line, so the LAST
such line in the log is the freq-only job's CPU time -- exactly what we want
to compare against (Gaussian's cost to supply the normal modes our algorithm
then classifies, not the cost of finding the stationary point itself).

Usage (from Github/scoring-functions/):
    py scripts/benchmark_cpu_time.py

Writes data/results/cpu_time_benchmark.csv with columns:
    molecule, N, n_basis, gaussian_freq_cpu_s, classifier_cpu_s,
    classifier_cpu_s_stddev, n_iterations, mp2_321g, method_basis

mp2_321g / method_basis are joined in from data/mol_list_method.csv (the
MP2_3-21G / current_method columns there) so downstream figure/manuscript
code can control for the dataset's mixed methods/basis sets (e.g. AsBr3 ran
at MP2/6-311G, not the library's standard MP2/3-21G). All molecules are kept
in the written CSV -- filtering by mp2_321g, if desired, is left to whatever
reads this file.

n_basis (added 2026-08-02, see gaussian_n_basis()) is Gaussian's own
``NBasis=`` count of AO basis functions -- the covariate that actually
drives SCF/MP2 cost. N (atom count) confounds heavy-element/ECP effects
into "size"; a molecule can have small N but a huge basis (e.g. one Xe
atom with diffuse/polarization functions) or vice versa. Included so
Gaussian's CPU-time fit can be checked against n_basis instead of N
without re-running the (slow) classifier timing. NOT meaningful for
classifier_cpu_s -- the classifier never touches AO basis functions, its
cost depends on N (geometry/mode-vector work) only. Plotting classifier
cost against n_basis would be a category error: n_basis is a property of
Gaussian's method/basis choice, not of the molecule's geometry.
"""
import glob
import os
import re
import statistics
import sys
import time
import timeit

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from main import load_inputs, build_scorer_and_final, resolve_dirs  # noqa: E402
from src.classifier import classify_all_modes, Thresholds  # noqa: E402

# Number of independently-calibrated batches to time per molecule, so we get
# a batch-to-batch stddev (real jitter) rather than a per-call one.
# time.process_time() on Windows is backed by GetProcessTimes(), whose
# actual OS-tick granularity (~15.6 ms) is far coarser than the API's
# reported 100 ns resolution -- timing individual sub-millisecond calls
# directly is dominated by that quantization, not real variance. Batching
# (each batch calibrated by timeit.Timer.autorange() to clear its own 0.2 s
# noise floor, per-call time = batch_total / batch_size) averages the
# quantization out within a batch; repeating several such batches gives a
# stddev that reflects genuine run-to-run variability instead.
_N_BATCHES = 10

_JOB_CPU_TIME_RE = re.compile(
    r"Job cpu time:\s*(\d+)\s*days\s+(\d+)\s*hours\s+(\d+)\s*minutes\s+([\d.]+)\s*seconds"
)

_NBASIS_RE = re.compile(r"NBasis=\s*(\d+)")


def _job_cpu_time_lines_to_seconds(log_text):
    """Return a list of CPU seconds, one per 'Job cpu time' line, in file order."""
    out = []
    for m in _JOB_CPU_TIME_RE.finditer(log_text):
        days, hours, minutes, seconds = m.groups()
        total = (int(days) * 86400 + int(hours) * 3600
                  + int(minutes) * 60 + float(seconds))
        out.append(total)
    return out


def gaussian_freq_cpu_seconds(log_path):
    """Extract the freq-only job step's CPU time (seconds) from a Gaussian log.

    Assumes the verified Link1-split pattern: exactly 2 "Job cpu time" lines
    (job 1 = opt, job 2 = freq-only single point via Link1/AllCheck), and
    returns the LAST one. Raises ValueError if a log doesn't have exactly 2
    lines, so a differently-structured log is reported/skipped rather than
    silently mis-timed.
    """
    with open(log_path, "r", errors="replace") as f:
        text = f.read()
    times = _job_cpu_time_lines_to_seconds(text)
    if len(times) != 2:
        raise ValueError(
            f"expected exactly 2 'Job cpu time' lines (opt job + Link1 freq-only "
            f"job), found {len(times)}"
        )
    return times[-1]


def gaussian_n_basis(log_path):
    """Extract the number of AO basis functions (Gaussian's ``NBasis=``) from
    a log, for use as a covariate that isolates basis-set size (driven by
    heavy/ECP-bearing elements and the chosen basis set) from atom count N.

    ``NBasis=`` is printed once per SCF cycle (dozens of times per job step),
    but is constant for a given molecule/method/basis -- both job steps (opt
    and the Link1 freq-only restart) use the same basis, so every occurrence
    in the file is expected to agree. Raises ValueError if the log has no
    ``NBasis=`` line, or if occurrences disagree (would indicate a
    basis-changing multi-step job this benchmark's assumptions don't cover),
    so a mismatch is reported/skipped rather than silently averaged over.
    """
    with open(log_path, "r", errors="replace") as f:
        text = f.read()
    values = {int(m.group(1)) for m in _NBASIS_RE.finditer(text)}
    if not values:
        raise ValueError("no 'NBasis=' line found")
    if len(values) > 1:
        raise ValueError(f"inconsistent NBasis= values across the log: {sorted(values)}")
    return values.pop()


def benchmark_classifier(mol_name, data_dir="data", thresholds=None,
                          n_batches=_N_BATCHES):
    """Time build_scorer_and_final() + classify_all_modes() (algorithm only,
    no log-parsing/CSV I/O) via time.process_time() (CPU time, Windows-safe).

    Warms the intermediate cache with one untimed load_inputs() call, then
    uses timeit.Timer(timer=time.process_time) to calibrate a batch size
    (autorange(): the smallest loop count whose total process_time clears
    0.2 s) and times `n_batches` such batches (repeat()). Per-call time is
    each batch's total / batch size, so individual sub-millisecond calls are
    never timed directly (see _N_BATCHES docstring re: Windows clock
    quantization). Returns (mean_s, stddev_s, n_iterations) where
    n_iterations = batch_size * n_batches.
    """
    # Warm-up: forces the parse + intermediate-cache write once, untimed.
    load_inputs(mol_name, "normal", data_dir, use_cache=True)

    def _one_call():
        raw, _ = load_inputs(mol_name, "normal", data_dir, use_cache=True)
        scorer, final = build_scorer_and_final(raw, "normal")
        classify_all_modes(scorer, final, thresholds)

    timer = timeit.Timer(stmt=_one_call, timer=time.process_time)
    batch_size, _ = timer.autorange()  # smallest n with total process_time >= 0.2s
    batch_totals = timer.repeat(repeat=n_batches, number=batch_size)
    per_call = [t / batch_size for t in batch_totals]

    mean_s = statistics.fmean(per_call)
    stddev_s = statistics.pstdev(per_call) if len(per_call) > 1 else 0.0
    n_iterations = batch_size * n_batches
    return mean_s, stddev_s, n_iterations


def main():
    data_dir = "data"
    dirs = resolve_dirs(data_dir)
    log_paths = sorted(glob.glob(os.path.join(dirs["logs"], "*.log")))
    thresholds = Thresholds.calibrated()

    rows = []
    skipped = []
    pattern_mismatches = []

    for log_path in log_paths:
        mol_name = os.path.splitext(os.path.basename(log_path))[0]

        try:
            gaussian_cpu_s = gaussian_freq_cpu_seconds(log_path)
            n_basis = gaussian_n_basis(log_path)
        except ValueError as e:
            pattern_mismatches.append((mol_name, str(e)))
            print(f"[SKIP] {mol_name}: log pattern mismatch -- {e}")
            continue

        try:
            raw, _ = load_inputs(mol_name, "normal", data_dir, use_cache=True)
            n_atoms = len(raw["atoms"])
        except (FileNotFoundError, ValueError) as e:
            skipped.append((mol_name, f"load_inputs failed: {e}"))
            print(f"[SKIP] {mol_name}: load_inputs failed -- {e}")
            continue

        try:
            mean_s, stddev_s, n_iter = benchmark_classifier(
                mol_name, data_dir, thresholds)
        except ValueError as e:
            # e.g. "No bond connectivity available" from build_scorer_and_final
            skipped.append((mol_name, f"classification failed: {e}"))
            print(f"[SKIP] {mol_name}: classification failed -- {e}")
            continue

        rows.append({
            "molecule": mol_name,
            "N": n_atoms,
            "n_basis": n_basis,
            "gaussian_freq_cpu_s": gaussian_cpu_s,
            "classifier_cpu_s": mean_s,
            "classifier_cpu_s_stddev": stddev_s,
            "n_iterations": n_iter,
        })
        print(f"[OK]   {mol_name}: N={n_atoms:3d}  n_basis={n_basis:4d}  "
              f"gaussian_freq={gaussian_cpu_s:9.4f}s  "
              f"classifier={mean_s*1e3:9.4f}ms +/- {stddev_s*1e3:.4f}ms  (n={n_iter})")

    import pandas as pd
    df = pd.DataFrame(rows).sort_values("N").reset_index(drop=True)

    method_path = os.path.join(data_dir, "mol_list_method.csv")
    method_df = pd.read_csv(method_path)[["molecule", "MP2_3-21G", "current_method"]]
    df = df.merge(method_df, on="molecule", how="left")
    df["mp2_321g"] = df["MP2_3-21G"].astype(bool)
    df["method_basis"] = df.apply(
        lambda r: "MP2/3-21G" if r["mp2_321g"] else r["current_method"], axis=1)
    df = df.drop(columns=["MP2_3-21G", "current_method"])

    out_path = os.path.join(dirs["results"], "cpu_time_benchmark.csv")
    df.to_csv(out_path, index=False)
    print(f"\nWrote {len(df)}-row benchmark -> {out_path}")

    if pattern_mismatches:
        print(f"\n{len(pattern_mismatches)} log(s) did NOT fit the assumed "
              f"Link1-split pattern (excluded):")
        for mol_name, reason in pattern_mismatches:
            print(f"  - {mol_name}: {reason}")

    if skipped:
        print(f"\n{len(skipped)} molecule(s) skipped (algorithm/input failure):")
        for mol_name, reason in skipped:
            print(f"  - {mol_name}: {reason}")

    return df


if __name__ == "__main__":
    main()
