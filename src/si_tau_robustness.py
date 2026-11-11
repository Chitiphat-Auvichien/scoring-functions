"""SI: robustness of tau_SB and the S/B ambiguity range (JCC revision, referee comments R1.9/R1.10).

    py -m src.si_tau_robustness                 # full run, 10^4 replicates per test
    py -m src.si_tau_robustness --n-rep 500     # quick look

How tau_SB is obtained (main text, ``src.figures.plot_tau_sb_derivation``): for a trial cut tau,
count the reference-labelled S/B modes of the single-centre (SC) calibration set that the rule
"S if s[V_S] >= tau, else B" gets wrong. The count E(tau) is a staircase; its tied minimum is the
SC gap, and tau_SB = its midpoint rounded to 0.01 (0.43, ``Thresholds.calibrated().tau_SB``).

This module perturbs the SC set, recalibrates each perturbed copy (a "replicate") with the same
rule, and records where the re-optimised cut tau* lands, whether the adopted tau_SB is still in the
minimum, and the regret E'(tau_SB) - min E' (the cost of keeping tau_SB). The tests:

    D1     molecule bootstrap: draw the 66 SC molecules with replacement (weights w_i)
    L1     every single relabelling of one flip unit (exhaustive)
    L2-L4  random relabelling, each unit flips with p = 0.01, 0.05, 0.20
    N1-N4  score noise s_i -> clip(s_i + N(0, sigma_s^2), 0, 1), sigma_s = 0.01, 0.02, 0.05, 0.10
    C1     D1 and L3 in the same replicate

The S/B ambiguity range is the 95 % percentile-bootstrap confidence interval of tau_SB from C1
replicates with noise of SD sigma_cal = 0.125 (the RMSE of s[V_S] against the PED stretching
fraction for the mixed test modes) on every calibration score. F(s), the share of calibrations
whose tau* <= s, is the share that call a mode of score s "S"; a mode is ambiguous iff
0.025 < F(s) < 0.975.

Data: the program's own per-member s[V_S] in data/results/library_scores.csv (kind == "internal",
ref_label stretch/bend; SC = mol_type ideal/non-ideal, test = mol_type test). Degenerate members
are scored one by one (no member-mean rule). A flip unit is one degenerate set (rows of
``units_csv``, columns molecule, freq, deg; its members share one reference label and flip
together) or a single mode.

Seeds: 20261029 (D, L, C; consumed in the order D1, L2, L3, L4, C1), 20261030 (N1-N4),
20261109 (ambiguity range; sigma_cal scan 0, 0.02, 0.05, 0.10, 0.125 in that order).

Outputs (results_dir):
    si_tau_robustness_summary.csv    one row per test: design (ID, question, replicate, R) + results
    si_tau_robustness_tests.tex      LaTeX tabular of the test design (booktabs)
    si_tau_robustness_ci.csv         confidence interval of tau_SB vs sigma_cal (adopted row flagged)
    si_tau_robustness_modes.csv      share of calibrations calling a mode S, modes near the cut
Figures (out_dir): fig_si_tau_robustness.{pdf,png}, fig_si_tau_ambiguity.{pdf,png}.
"""
from __future__ import annotations

import argparse
import os
import time

import numpy as np
import pandas as pd

from src.classifier import Thresholds
from src.figures import AMBIGUITY_LABEL, AMBIGUITY_RANGE, COLORS, _savefig, _style, plt

GRID = np.round(np.arange(0.0, 1.0005, 0.001), 3)
NG = len(GRID)
SEED_DLC = 20261029
SEED_N = 20261030
SEED_CI = 20261109
SIGMA_CAL = 0.125
COVERAGE = 0.95
C1_FLIP_P = 0.05
FLIP_RATES = (("L2", 0.01), ("L3", 0.05), ("L4", 0.20))
NOISE_SIGMAS = (("N1", 0.01), ("N2", 0.02), ("N3", 0.05), ("N4", 0.10))
CI_SIGMA_SCAN = (0.0, 0.02, 0.05, 0.10, SIGMA_CAL)
N_REP_CI = 100000        # replicates for the ambiguity range (Monte Carlo error of each end ~ +/-0.0004)
CI_BATCH = 10000         # replicates per batch (memory; batch spread -> Monte Carlo standard error)
BAND_TOL = 0.002        # |this run's CI end - AMBIGUITY_RANGE end| up to which the figures shade the constant
                        # (~5 Monte Carlo standard errors at N_REP_CI = 10^5)
DEFAULT_UNITS_CSV = os.path.join("data", "results", "si_degenerate_sets.csv")
GROUP_COLOR = {"Data": COLORS["sens_single_centre"], "Labels": COLORS["mixed_ext"],
               "Scores": "#E69F00", "Combined": "#000000"}     # Okabe-Ito green / pink / orange / black
TAU_LABEL = r"$\tau_\mathrm{SB}$"
SV_LABEL = r"$s[\mathrm{V_S}]$"


# ---------------------------------------------------------------------------- data
def load_sb_modes(library_csv):
    """Reference-labelled S/B internal modes -> (SC frame, test frame), each with y = 1 for S."""
    lib = pd.read_csv(library_csv)
    lib = lib[(lib["kind"] == "internal") & lib["ref_label"].isin(["stretch", "bend"])].copy()
    lib["y"] = (lib["ref_label"] == "stretch").astype(int)
    sc = lib[lib["mol_type"].isin(["ideal", "non-ideal"])].reset_index(drop=True)
    te = lib[lib["mol_type"] == "test"].reset_index(drop=True)
    return sc, te


def flip_units(df, units_csv, tol=0.01):
    """Integer flip-unit id per row of ``df`` (columns molecule, freq, ref_label; RangeIndex).

    Rows of one degenerate set of ``units_csv`` (columns molecule, freq, deg) share an id; every
    other row is its own unit. Returns (ids 0..n_units-1, info dict)."""
    if not os.path.exists(units_csv):
        raise FileNotFoundError(f"degenerate-set file not found: {units_csv}")
    sets = pd.read_csv(units_csv)
    mol = df["molecule"].to_numpy()
    freq = df["freq"].to_numpy(float)
    unit = np.full(len(df), -1, int)
    uid, mixed = 0, []
    for s in sets.itertuples():
        m = np.flatnonzero((mol == s.molecule) & (np.abs(freq - s.freq) < tol))
        if len(m) == 0:
            continue
        if len(m) > s.deg:
            raise ValueError(f"{s.molecule} {s.freq}: {len(m)} rows matched a {s.deg}-fold set")
        if df["ref_label"].iloc[m].nunique() > 1:
            mixed.append(f"{s.molecule} {s.freq:.1f}")
        unit[m] = uid
        uid += 1
    n_sets = uid
    for i in np.flatnonzero(unit < 0):
        unit[i] = uid
        uid += 1
    return unit, dict(n_units=uid, n_sets=n_sets,
                      n_modes_in_sets=int((np.bincount(unit)[unit] > 1).sum()),
                      mixed_label_sets=mixed)


# ---------------------------------------------------------------------------- recalibration
def error_counts(s, y, w):
    """Weighted error counts E'(tau) on GRID for each row of (s, y, w), all (R, n).

    An S mode is wrong at GRID[g] iff s < GRID[g], i.e. g >= k with k = searchsorted(GRID, s,
    'right'); a B mode is wrong iff g < k. Histogram k per row, then cumulative sums."""
    r, n = s.shape
    k = np.searchsorted(GRID, s, side="right")
    flat = (np.arange(r)[:, None] * (NG + 1) + k).ravel()
    h_s = np.bincount(flat, (w * (y == 1)).ravel(), minlength=r * (NG + 1)).reshape(r, NG + 1)
    h_b = np.bincount(flat, (w * (y == 0)).ravel(), minlength=r * (NG + 1)).reshape(r, NG + 1)
    return np.cumsum(h_s, axis=1)[:, :NG] + h_b.sum(1, keepdims=True) - np.cumsum(h_b, axis=1)[:, :NG]


def optimum_rows(E, ref=0.5, tol=1e-9):
    """Tied-minimum run of each row of E: the contiguous run that holds the tied grid point closest
    to ``ref`` (tie-break of src.calibrate). Returns (lo index, hi index, minimum)."""
    emin = E.min(axis=1, keepdims=True)
    ok = E <= emin + tol
    best = np.where(ok, np.abs(GRID - ref) + 1e-9 * GRID, np.inf).argmin(axis=1)
    idx = np.arange(NG)
    prev_bad = np.maximum.accumulate(np.where(~ok, idx, -1), axis=1)
    next_bad = np.minimum.accumulate(np.where(~ok, idx, NG)[:, ::-1], axis=1)[:, ::-1]
    r = np.arange(len(E))
    return prev_bad[r, best] + 1, next_bad[r, best] - 1, emin[:, 0]


def recalibrate(s, y, w, gi_tau, chunk=1000):
    """Re-optimise the cut for every replicate. s, y, w: (n,) or (R, n). Returns a dict of (R,)
    arrays: tau_star (run midpoint), lo, hi (run ends), emin, e_tau (error at the adopted cut),
    wsum (total weight)."""
    arrs = [np.asarray(a, float) for a in (s, y, w)]
    n_rep = max([a.shape[0] for a in arrs if a.ndim == 2] or [1])
    n = arrs[0].shape[-1]
    out = {k: np.empty(n_rep) for k in ("tau_star", "lo", "hi", "emin", "e_tau", "wsum")}
    for i in range(0, n_rep, chunk):
        j = min(i + chunk, n_rep)
        s_, y_, w_ = (np.broadcast_to(a[i:j] if a.ndim == 2 else a, (j - i, n)) for a in arrs)
        E = error_counts(s_, y_, w_)
        lo, hi, emin = optimum_rows(E)
        out["lo"][i:j], out["hi"][i:j] = GRID[lo], GRID[hi]
        out["tau_star"][i:j] = 0.5 * (GRID[lo] + GRID[hi])
        out["emin"][i:j], out["e_tau"][i:j] = emin, E[:, gi_tau]
        out["wsum"][i:j] = w_.sum(axis=1)
    return out


def molecule_weights(codes, n_mol, n_rep, rng):
    """Cluster-bootstrap weights (n_rep, n_modes): draws of the molecule of each mode."""
    pick = rng.integers(0, n_mol, size=(n_rep, n_mol))
    cnt = np.bincount((np.arange(n_rep)[:, None] * n_mol + pick).ravel(),
                      minlength=n_rep * n_mol).reshape(n_rep, n_mol)
    return cnt[:, codes].astype(float)


def flipped_labels(y0, unit, n_unit, p, n_rep, rng):
    flip = rng.random((n_rep, n_unit)) < p
    return np.where(flip[:, unit], 1 - y0, y0).astype(float)


def summarize(res, e_test=None):
    """Distribution of tau*, share of replicates in which the adopted cut is optimal, regret."""
    t = res["tau_star"]
    reg = res["e_tau"] - res["emin"]
    reg_pct = 100 * reg / res["wsum"]
    row = dict(tau_star_mean=t.mean(), tau_star_median=np.median(t),
               tau_star_p2_5=np.percentile(t, 2.5), tau_star_p25=np.percentile(t, 25),
               tau_star_p75=np.percentile(t, 75), tau_star_p97_5=np.percentile(t, 97.5),
               tau_star_min=t.min(), tau_star_max=t.max(),
               pct_tau_optimal=100 * np.mean(reg < 1e-9),
               regret_p95_pct=np.percentile(reg_pct, 95), regret_max_pct=reg_pct.max(),
               regret_p95_modes=np.percentile(reg, 95), regret_max_modes=reg.max())
    if e_test is not None:
        gi = np.clip(np.searchsorted(GRID, t - 1e-9), 0, NG - 1)   # tau* on the 0.0005 grid -> next point
        row.update(test_err_at_tau_star_median_pct=np.median(e_test[gi]),
                   test_err_at_tau_star_p95_pct=np.percentile(e_test[gi], 95),
                   test_err_at_tau_star_max_pct=e_test[gi].max())
    return row


# ---------------------------------------------------------------------------- tables
def _r_text(n):
    k = np.log10(n)
    return f"$10^{{{int(round(k))}}}$" if n >= 10 and np.isclose(k, round(k)) else f"{n:,}".replace(",", r"\,")


def write_tests_tex(design, path, n_rep, n_unit, n_mol):
    """booktabs tabular of the test design (one line per test family, as in the SI table)."""
    r = _r_text(n_rep)
    rows = [
        ("D1", "other molecules?",
         rf"draw {n_mol} molecules with replacement; $w_i$ = draws of the molecule of mode $i$", r),
        ("L1", "one wrong label?", "flip the label of one unit", f"{n_unit} (all)"),
        ("L2--L4", "random wrong labels?",
         "each unit flips with $p=" + "$, $".join(f"{p:.2f}" for _, p in FLIP_RATES) + "$", f"{r} each"),
        ("N1--N4", "scores slightly off?",
         r"$s_i\to s_i+\varepsilon_i$, $\varepsilon_i\sim N(0,\sigma_s^2)$, $\sigma_s="
         + "$, $".join(f"{sg:.2f}" for _, sg in NOISE_SIGMAS) + "$; clipped to $[0,1]$", f"{r} each"),
        ("C1", "molecules and labels together", "D1 and L3 in the same replicate", r),
    ]
    lines = [
        "% Robustness tests of tau_SB (SI). Written by src/si_tau_robustness.py; do not edit.",
        "% Needs booktabs and array. A unit = one degenerate set (members flip together) or one mode.",
        f"% A replicate keeps everything not named. Seeds {SEED_DLC} (D, L, C) and {SEED_N} (N).",
        r"\begin{tabular}{@{}l>{\raggedright\arraybackslash}p{3.3cm}>{\raggedright\arraybackslash}p{7.0cm}l@{}}",
        r"\toprule",
        r"ID & question & one replicate & $R$\\",
        r"\midrule",
    ] + [" & ".join(x) + r"\\" for x in rows] + [r"\bottomrule", r"\end{tabular}", ""]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    return path


# ---------------------------------------------------------------------------- figures
def plot_robustness(summary, draws, gap, tau, band, out_dir, label="fig_si_tau_robustness"):
    """(a) tau* per test: thin line 2.5-97.5 %, thick bar 25-75 %, open circle median, diamond
    mean; dark grey SC gap, light grey ambiguity range, dashed adopted tau_SB. (b) share of
    replicates in which tau_SB is still optimal. (c) 95th percentile of the regret (% of modes)."""
    _style()
    rows = summary.to_dict("records")
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(9.0, 4.6), sharey=True,
                                     gridspec_kw=dict(width_ratios=[2.4, 1, 1], wspace=0.08))
    pos, y_, ghead = {}, 0.0, {}
    for grp in reversed(list(GROUP_COLOR)):
        for r in reversed([r for r in rows if r["group"] == grp]):
            pos[r["id"]] = y_
            y_ += 1
        ghead[grp] = y_ - 0.15
        y_ += 0.9
    xa = (0.25, 0.60)
    edge = 0.008
    clip = lambda x: min(max(x, xa[0] + edge), xa[1] - edge)
    a1.axvspan(*band, color=COLORS["background"], alpha=0.22, lw=0)
    a1.axvspan(*gap, color=COLORS["background"], alpha=0.6, lw=0)
    for a in (a1,):
        a.axvline(tau, color=COLORS["threshold"], ls="--", lw=1.0, zorder=1)
    reg_max = max(r["regret_p95_pct"] for r in rows)
    xr = max(2.0, np.ceil(reg_max * 1.15 * 2) / 2)
    for r in rows:
        p, c, d = pos[r["id"]], GROUP_COLOR[r["group"]], draws[r["id"]]
        lo, q1, md, q3, hi = np.percentile(d, [2.5, 25, 50, 75, 97.5])
        a1.plot([clip(lo), clip(hi)], [p, p], color=c, lw=1.1, solid_capstyle="butt")
        a1.plot([clip(q1), clip(q3)], [p, p], color=c, lw=5, solid_capstyle="butt")
        for x_, e_, mk in ((lo, xa[0] + edge, "<"), (hi, xa[1] - edge, ">")):
            if (x_ < e_) if mk == "<" else (x_ > e_):
                a1.plot(e_, p, mk, color=c, ms=5)
        a1.plot(clip(md), p, "o", mfc="white", mec=c, mew=1.1, ms=5.5, zorder=3)
        a1.plot(clip(d.mean()), p, "D", mfc="k", mec="white", mew=0.5, ms=3.6, zorder=4)
        fo = r["pct_tau_optimal"]
        a2.barh(p, fo, color=c, height=0.62)
        txt = f"{fo:.1f}" if fo < 99.95 else "100"
        if fo >= 30:
            a2.text(fo - 2, p, txt, color="white", ha="right", va="center", fontsize=9)
        else:
            a2.text(fo + 2, p, txt, color=c, ha="left", va="center", fontsize=9)
        rg = r["regret_p95_pct"]
        a3.plot([0, rg], [p, p], color=c, lw=1.0)
        a3.plot(rg, p, "o", color=c, ms=5.5, clip_on=False, zorder=3)
    a1.plot([], [], "o", mfc="white", mec="k", mew=1.1, ms=5.5, label="median")
    a1.plot([], [], "D", mfc="k", mec="white", ms=3.6, label="mean")
    a1.legend(loc="lower right", frameon=True, framealpha=0.95, edgecolor="none", fontsize=9,
              handletextpad=0.3, borderpad=0.3, labelspacing=0.25)
    ids = [r["id"] for r in rows]
    a1.set_yticks([pos[i] for i in ids])
    a1.set_yticklabels([rf"$\mathbf{{{r['id']}}}$  {r['short']}" for r in rows])
    a1.set_xlim(*xa)
    a1.set_xticks([0.3, 0.4, 0.5, 0.6])
    a1.set_xlabel(r"Re-optimised " + TAU_LABEL + r" ($\tau^*$)")
    a2.set_xlim(0, 108)
    a2.set_xticks([0, 50, 100])
    a2.axvline(100, color=COLORS["background"], lw=0.7)
    a2.set_xlabel(f"{tau:.2f} still\noptimal (%)")
    a3.set_xlim(0, xr)
    a3.set_xticks(np.arange(0, xr + 1e-9, 1.0))
    n_sc = summary["n_sc_modes"].iloc[0]
    a3.axvline(100 / n_sc, color=COLORS["threshold"], lw=0.8, ls=":")
    a3.text(100 / n_sc + 0.06, max(ghead.values()), "1 mode", fontsize=9, color=COLORS["threshold"],
            va="center")
    a3.set_xlabel(f"Cost of keeping {tau:.2f},\n95th pct (% of modes)")
    for grp, yh in ghead.items():
        a1.text(xa[0] + 0.004, yh, grp, color=GROUP_COLOR[grp], fontweight="bold", va="center")
    a1.set_ylim(-0.7, max(ghead.values()) + 0.6)
    for a, t in zip((a1, a2, a3), "abc"):
        a.text(0.0 if a is a1 else 0.03, 1.01, f"({t})", transform=a.transAxes, fontweight="bold",
               va="bottom")
        a.tick_params(axis="y", length=0)
    paths = _savefig(fig, out_dir, label)
    plt.close(fig)
    return paths


def plot_ambiguity(tau_draws, ci, tau, marks, out_dir, label="fig_si_tau_ambiguity"):
    """(a) histogram of the re-optimised cuts with the 95 % CI shaded; (b) share of calibrations
    that call a mode of score s "S" (step curve), dotted 2.5/97.5 % lines, marked modes."""
    _style()
    a_ = (1 - COVERAGE) / 2
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.0, 3.6), gridspec_kw=dict(wspace=0.28))
    bins = np.arange(0.2475, 0.6026, 0.005)
    a1.axvspan(*ci, color=COLORS["background"], alpha=0.35, lw=0)
    a1.hist(tau_draws, bins=bins, color="#777777", edgecolor="white", lw=0.3)
    a1.axvline(tau, color=COLORS["threshold"], ls="--", lw=1.0)
    # three decimals: both ends sit on a two-decimal rounding boundary (see figures.AMBIGUITY_LABEL)
    a1.text(0.5 * (ci[0] + ci[1]), 1.01, f"95% CI {ci[0]:.3f}–{ci[1]:.3f}", ha="center",
            va="bottom", fontsize=10, transform=a1.get_xaxis_transform())
    a1.set_xlim(0.25, 0.60)
    a1.set_xlabel(r"Re-optimised " + TAU_LABEL + r" ($\tau^*$)")
    a1.set_ylabel("Calibrations")
    xx = np.linspace(0.20, 0.65, 901)
    srt = np.sort(tau_draws)
    share = np.searchsorted(srt, xx, side="right") / len(srt)
    a2.axvspan(*ci, color=COLORS["background"], alpha=0.35, lw=0)
    a2.step(xx, share, where="post", color="k", lw=1.2)
    for h in (a_, 1 - a_):
        a2.axhline(h, color=COLORS["threshold"], lw=0.8, ls=":")
    a2.axvline(tau, color=COLORS["threshold"], ls="--", lw=1.0)
    a2.text(0.203, 1 - a_ + 0.015, f"{100 * (1 - a_):.1f}%", fontsize=9, va="bottom",
            color=COLORS["threshold"])
    a2.text(0.203, a_ + 0.015, f"{100 * a_:.1f}%", fontsize=9, va="bottom", color=COLORS["threshold"])
    for m in marks:
        col = COLORS["stretching_ref"] if m["literature"] == "S" else COLORS["bending_ref"]
        a2.plot(m["s"], m["share_S"], "o", ms=5.5, mfc="white", mec=col, mew=1.3, zorder=3)
        if m.get("xy_text") is not None:       # data-coordinate label (shared by several marks)
            a2.annotate(m["text"], (m["s"], m["share_S"]), xytext=m["xy_text"], textcoords="data",
                        fontsize=9, ha="left", va="center", color=col,
                        arrowprops=dict(arrowstyle="-", lw=0.5, color=col, shrinkA=1, shrinkB=3))
        else:
            a2.annotate(m["text"], (m["s"], m["share_S"]), xytext=m["offset"], textcoords="offset points",
                        fontsize=9, ha="left" if m["offset"][0] > 0 else "right", va="center", color=col,
                        arrowprops=dict(arrowstyle="-", lw=0.5, color=col, shrinkA=0, shrinkB=3))
    a2.set_xlim(0.20, 0.65)
    a2.set_ylim(-0.03, 1.05)
    a2.set_yticks([0, 0.25, 0.5, 0.75, 1])
    a2.set_xlabel(SV_LABEL + " of a mode")
    a2.set_ylabel("Share of calibrations\nlabelling the mode S")
    for a, t in zip((a1, a2), "ab"):
        a.text(-0.16 if t == "a" else -0.2, 1.01, f"({t})", transform=a.transAxes, fontweight="bold",
               va="bottom")
    paths = _savefig(fig, out_dir, label)
    plt.close(fig)
    return paths


# ---------------------------------------------------------------------------- benzene marks
BENZENE_MODES = ((1237.29, r"$\nu_{10}$"), (1056.39, r"$\nu_{14}$"), (1532.85, r"$\nu_{13}$"))


def _benzene_rows(data_dir, results_dir):
    """Per-member s[V_S] and literature label of benzene nu10, nu14 and nu13 (Wilson numbering)."""
    bz = pd.read_csv(os.path.join(results_dir, "C6H6_normal.csv"))
    ref = pd.read_csv(os.path.join(data_dir, "characterised_modes.csv"), encoding="utf-8-sig")
    ref = ref[ref["molecule"] == "C6H6"]
    out = []
    for f0, name in BENZENE_MODES:
        lab = ref.loc[(ref["freq"] - f0).abs() < 0.05, "type"]
        for k, r in enumerate(bz[(bz["Freq"] - f0).abs() < 0.05].itertuples()):
            out.append(dict(set="benzene", molecule="C6H6", freq=r.Freq, name=name, member=k + 1,
                            literature=lab.iloc[0] if len(lab) else "", s=float(r.V_Stretch)))
    return out


# ---------------------------------------------------------------------------- driver
def run(data_dir="data", out_dir="data/figures", results_dir="data/results", n_rep=10000,
        n_rep_ci=N_REP_CI, units_csv=DEFAULT_UNITS_CSV, make_figures=True, verbose=False):
    """Run every robustness test and the ambiguity range; write CSVs, the LaTeX table and figures.

    Returns a dict of the headline numbers (and the paths written)."""
    t0 = time.time()
    tau = Thresholds.calibrated(path=os.path.join(results_dir, "thresholds.json")).tau_SB
    gi_tau = int(np.argmin(np.abs(GRID - tau)))
    sc, te = load_sb_modes(os.path.join(results_dir, "library_scores.csv"))
    s0, y0 = sc["V_Stretch"].to_numpy(float), sc["y"].to_numpy(float)
    one = np.ones_like(s0)
    n_sc = len(sc)
    mols = sc["molecule"].unique()
    codes = pd.Categorical(sc["molecule"], categories=mols).codes
    unit, uinfo = flip_units(sc, units_csv)
    n_unit = uinfo["n_units"]

    # unperturbed calibration
    base = recalibrate(s0, y0, one, gi_tau)
    run_grid = (float(base["lo"][0]), float(base["hi"][0]))          # tied grid points
    # the gap in score space: (highest score below the run, lowest score at/above its end]
    gap = (float(s0[s0 < run_grid[0]].max()), float(s0[s0 >= run_grid[1]].min()))
    s_te, y_te = te["V_Stretch"].to_numpy(float), te["y"].to_numpy(float)
    e_test = 100 * error_counts(s_te[None], y_te[None], np.ones((1, len(te))))[0] / len(te)
    sc_err_tau = int(round(base["e_tau"][0]))
    te_err_tau = int(round(e_test[gi_tau] * len(te) / 100))

    design = {
        "D1": ("Data", "molecule bootstrap", "other molecules?",
               f"draw {len(mols)} molecules with replacement; w_i = draws of the molecule of mode i"),
        "L1": ("Labels", "every single flip", "one wrong label?", "flip the label of one unit"),
        "L2": ("Labels", "random flips 1%", "random wrong labels?", "each unit flips with p = 0.01"),
        "L3": ("Labels", "random flips 5%", "random wrong labels?", "each unit flips with p = 0.05"),
        "L4": ("Labels", "random flips 20%", "random wrong labels?", "each unit flips with p = 0.20"),
        "N1": ("Scores", r"noise $\sigma_s=0.01$", "scores slightly off?", "s_i -> s_i + N(0, 0.01^2)"),
        "N2": ("Scores", r"noise $\sigma_s=0.02$", "scores slightly off?", "s_i -> s_i + N(0, 0.02^2)"),
        "N3": ("Scores", r"noise $\sigma_s=0.05$", "scores slightly off?", "s_i -> s_i + N(0, 0.05^2)"),
        "N4": ("Scores", r"noise $\sigma_s=0.10$", "scores slightly off?", "s_i -> s_i + N(0, 0.10^2)"),
        "C1": ("Combined", "mol. bootstrap + 5% flips", "molecules and labels together",
               "D1 and L3 in the same replicate"),
    }
    res, seeds = {}, {}
    rng = np.random.default_rng(SEED_DLC)
    res["D1"] = recalibrate(s0, y0, molecule_weights(codes, len(mols), n_rep, rng), gi_tau)
    res["L1"] = recalibrate(s0, np.where(unit[None, :] == np.arange(n_unit)[:, None], 1 - y0, y0), one,
                            gi_tau)
    for tid, p in FLIP_RATES:
        res[tid] = recalibrate(s0, flipped_labels(y0, unit, n_unit, p, n_rep, rng), one, gi_tau)
    y_c1 = flipped_labels(y0, unit, n_unit, C1_FLIP_P, n_rep, rng)
    res["C1"] = recalibrate(s0, y_c1, molecule_weights(codes, len(mols), n_rep, rng), gi_tau)
    rng_n = np.random.default_rng(SEED_N)
    for tid, sg in NOISE_SIGMAS:
        res[tid] = recalibrate(np.clip(s0 + rng_n.normal(0, sg, (n_rep, n_sc)), 0, 1), y0, one, gi_tau)
    for tid in design:
        seeds[tid] = "" if tid == "L1" else (SEED_N if tid.startswith("N") else SEED_DLC)

    rows = []
    for tid, (grp, short, question, recipe) in design.items():
        rows.append(dict(id=tid, group=grp, short=short, question=question, replicate=recipe,
                         n_rep=len(res[tid]["tau_star"]), seed=seeds[tid], tau_SB=tau,
                         n_sc_modes=n_sc, n_sc_molecules=len(mols), n_units=n_unit,
                         sc_gap_lo=gap[0], sc_gap_hi=gap[1], **summarize(res[tid], e_test)))
    summary = pd.DataFrame(rows)
    # how many single flips move the gap / exclude tau_SB (L1 is exhaustive)
    l1 = res["L1"]
    n_l1_gap = int(np.sum(~(np.isclose(l1["lo"], run_grid[0]) & np.isclose(l1["hi"], run_grid[1]))))
    l1_bad = np.flatnonzero(l1["e_tau"] - l1["emin"] > 1e-9)
    l1_bad_modes = [f"{sc['molecule'].iloc[np.flatnonzero(unit == u)[0]]} "
                    f"{sc['freq'].iloc[np.flatnonzero(unit == u)[0]]:.0f}" for u in l1_bad]

    # ambiguity range: C1 + score noise sigma_cal, sigma_cal scan. n_rep_ci replicates, drawn in
    # batches of CI_BATCH (memory); the spread of the batch percentiles gives the Monte Carlo
    # standard error of each end (batch means: SD of the batch values / sqrt(number of batches)).
    rng_ci = np.random.default_rng(SEED_CI)
    ci_rows, ci_draws = [], {}
    lo_q, hi_q = 100 * (1 - COVERAGE) / 2, 100 * (1 + COVERAGE) / 2
    batch = min(CI_BATCH, n_rep_ci)
    n_batch = max(1, n_rep_ci // batch)
    for sg in CI_SIGMA_SCAN:
        parts = []
        for _ in range(n_batch):
            w = molecule_weights(codes, len(mols), batch, rng_ci)
            yy = flipped_labels(y0, unit, n_unit, C1_FLIP_P, batch, rng_ci)
            ss = np.clip(s0 + rng_ci.normal(0, sg, (batch, n_sc)), 0, 1) if sg > 0 else s0
            parts.append(recalibrate(ss, yy, w, gi_tau)["tau_star"])
        t = np.concatenate(parts)
        ci_draws[sg] = t
        lo95, hi95 = np.percentile(t, [lo_q, hi_q])
        if n_batch > 1:
            bq = np.array([np.percentile(p, [lo_q, hi_q]) for p in parts])
            mcse = bq.std(axis=0, ddof=1) / np.sqrt(n_batch)
        else:
            mcse = (np.nan, np.nan)
        ci_rows.append(dict(sigma_cal=sg, adopted=bool(np.isclose(sg, SIGMA_CAL)), coverage=COVERAGE,
                            n_rep=len(t), n_batches=n_batch, seed=SEED_CI, flip_p=C1_FLIP_P,
                            tau_mean=t.mean(), tau_median=np.median(t), tau_sd=t.std(ddof=1),
                            ci90_lo=np.percentile(t, 5), ci90_hi=np.percentile(t, 95),
                            ci95_lo=lo95, ci95_hi=hi95, ci95_lo_mcse=mcse[0], ci95_hi_mcse=mcse[1],
                            n_sc_modes_inside_ci95=int(((s0 >= lo95) & (s0 <= hi95)).sum())))
    ci_tab = pd.DataFrame(ci_rows)
    adopted = ci_tab[ci_tab["adopted"]].iloc[0]
    ci = (float(adopted["ci95_lo"]), float(adopted["ci95_hi"]))
    tau_ci = np.sort(ci_draws[SIGMA_CAL])
    share_s = lambda s: float(np.searchsorted(tau_ci, s, side="right") / len(tau_ci))

    # modes near the cut: share of calibrations that call them S
    a_ = (1 - COVERAGE) / 2
    near = sc[(sc["V_Stretch"] >= ci[0] - 0.03) & (sc["V_Stretch"] <= ci[1] + 0.03)].sort_values("V_Stretch")
    mrows = [dict(set="SC", molecule=r.molecule, freq=r.freq, name="", member=0,
                  literature="S" if r.ref_label == "stretch" else "B", s=float(r.V_Stretch))
             for r in near.itertuples()]
    mrows += _benzene_rows(data_dir, results_dir)
    modes = pd.DataFrame(mrows)
    modes["share_S"] = [share_s(v) for v in modes["s"]]
    modes["label_at_tau"] = np.where(modes["s"] >= tau, "S", "B")
    modes["reading"] = np.where((modes["share_S"] > a_) & (modes["share_S"] < 1 - a_), "ambiguous",
                                "reliable")
    modes = modes.sort_values(["s", "molecule"]).reset_index(drop=True)

    os.makedirs(results_dir, exist_ok=True)
    paths = dict(summary_csv=os.path.join(results_dir, "si_tau_robustness_summary.csv"),
                 ci_csv=os.path.join(results_dir, "si_tau_robustness_ci.csv"),
                 modes_csv=os.path.join(results_dir, "si_tau_robustness_modes.csv"),
                 tests_tex=os.path.join(results_dir, "si_tau_robustness_tests.tex"))
    summary.to_csv(paths["summary_csv"], index=False)
    ci_tab.to_csv(paths["ci_csv"], index=False)
    modes.to_csv(paths["modes_csv"], index=False)
    write_tests_tex(design, paths["tests_tex"], n_rep, n_unit, len(mols))

    # The figures shade the adopted range (src.figures.AMBIGUITY_RANGE) when this run reproduces it
    # within Monte Carlo error (each end +/- ~0.0004 at 10^5 replicates); otherwise this run's CI.
    band_is_constant = bool(max(abs(ci[0] - AMBIGUITY_RANGE[0]), abs(ci[1] - AMBIGUITY_RANGE[1])) <= BAND_TOL)
    band = tuple(AMBIGUITY_RANGE) if band_is_constant else ci
    if make_figures:
        draws = {tid: r["tau_star"] for tid, r in res.items()}
        paths["fig_robustness"] = plot_robustness(summary, draws, gap, tau, band, out_dir)
        sel = []
        offsets = {"SBr4": (-14, 12), "NBr3": (-16, 10)}
        for m in modes[modes["set"] == "SC"].itertuples():
            key = (m.molecule, round(m.freq))
            if key in (("SBr4", 135), ("NBr3", 188)):
                sel.append(dict(s=m.s, share_S=m.share_S, literature=m.literature,
                                text=f"{m.molecule[:-1]}$_{m.molecule[-1]}$ {m.freq:.0f}",
                                xy_text=None, offset=offsets[m.molecule]))
        bz = modes[modes["set"] == "benzene"]
        for m in bz[bz["name"] == r"$\nu_{10}$"].itertuples():
            sel.append(dict(s=m.s, share_S=m.share_S, literature=m.literature, text=m.name,
                            xy_text=None, offset=(14, -12)))
        nu14 = bz[bz["name"] == r"$\nu_{14}$"]
        for k, m in enumerate(nu14.itertuples()):           # e1u pair: one label, a leader to each member
            sel.append(dict(s=m.s, share_S=m.share_S, literature=m.literature,
                            text=r"$\nu_{14}$ pair" if k == 0 else "", xy_text=(0.55, 0.78), offset=None))
        paths["fig_ambiguity"] = plot_ambiguity(ci_draws[SIGMA_CAL], band, tau, sel, out_dir)

    out = dict(
        tau_SB=tau, sc_gap=gap, n_sc_modes=n_sc, n_sc_molecules=len(mols), n_test_modes=len(te),
        sc_errors_at_tau=sc_err_tau, test_errors_at_tau=te_err_tau, n_units=n_unit,
        n_degenerate_sets=uinfo["n_sets"], mixed_label_sets=uinfo["mixed_label_sets"],
        l1_n_gap_changed=n_l1_gap, l1_excluding_tau=l1_bad_modes,
        ci95=ci, ci90=(float(adopted["ci90_lo"]), float(adopted["ci90_hi"])),
        ci95_mcse=(float(adopted["ci95_lo_mcse"]), float(adopted["ci95_hi_mcse"])), n_rep_ci=int(adopted["n_rep"]),
        ci_tau_sd=float(adopted["tau_sd"]), n_sc_inside_ci=int(adopted["n_sc_modes_inside_ci95"]),
        ambiguity_range_constant=AMBIGUITY_RANGE, figure_band=band, figure_band_is_constant=band_is_constant,
        summary=summary, ci_table=ci_tab, modes=modes,
        n_rep=n_rep, units_csv=units_csv, runtime_s=time.time() - t0, paths=paths,
    )
    if verbose:
        pd.set_option("display.width", 250)
        cols = ["id", "n_rep", "tau_star_median", "tau_star_mean", "tau_star_p2_5", "tau_star_p25",
                "tau_star_p75", "tau_star_p97_5", "pct_tau_optimal", "regret_p95_pct", "regret_p95_modes",
                "regret_max_pct"]
        print(summary[cols].round(4).to_string(index=False))
        print(ci_tab.round(4).to_string(index=False))
        print(modes.round(3).to_string(index=False))
        print({k: v for k, v in out.items() if k not in ("summary", "ci_table", "modes")})
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n-rep", type=int, default=10000)
    ap.add_argument("--n-rep-ci", type=int, default=N_REP_CI)
    ap.add_argument("--units-csv", default=DEFAULT_UNITS_CSV)
    ap.add_argument("--no-figures", action="store_true")
    a = ap.parse_args(argv)
    run(n_rep=a.n_rep, n_rep_ci=a.n_rep_ci, units_csv=a.units_csv, make_figures=not a.no_figures, verbose=True)


if __name__ == "__main__":
    main()
