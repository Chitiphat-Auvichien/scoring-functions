# JCC Scoring/Classification Program — Implementation Plan & Progress Checklist

> Living checklist. Tick `[x]` as parts are completed; pick up unchecked items in any later session.
> Companion to the JCC manuscript `JCC/JCC_man_scoring/JCC_temp_LaTeXtemplate.tex` and the
> plan slides `JCC/Scoring_Manuscript_Plan_2026-06-29.pdf`. Last updated: 2026-06-29.

## Context

The JCC manuscript ("A Unified, Reference-Free Framework for Classifying the 3N modes of molecular
motion") is written, but its tables/figures are mostly `TODO-DATA`. The program currently implements
only **Step 1** (per-mode `Tx..Tz, Rx..Rz, V_Stretch` scoring → CSV). Goal: extend it into the
reference implementation that reproduces every manuscript table and figure. Build order: **core engine
first**, validated on data in hand (water, benzene, gramicidin), then scale out using precomputed
library scores from the Excel file.

### Locked decisions
- Hydride-library logs live off-server; **their scores are already in
  `data/vibrational-scoring-functions.xlsx`** → ingest precomputed scores + reference labels, do not re-score.
- Program **generates figures** (matplotlib), reproducing every manuscript figure.
- Sequence: **core engine first**.

### Authoritative spec (from the `.tex`)
- `s[T_Q] = (1/N) Σ unit(d_A)·Q̂`, atoms with `|d_A|>ε_disp`; range `[-1,1]`.
- `s[R_Q] = 1/(N−N_Q) Σ unit(ω_Q^A)·Q̂`, `ω_Q^A = ((r_A−(r_A·Q̂)Q̂)×d_A)/|r_A−(r_A·Q̂)Q̂|²`; `N_Q` = on-axis atoms excluded.
- `s[V_S] = (1/Σ|Δb|²) Σ |Δb_AB|²·|unit(Δb_AB)·b̂_AB|`, `Δb_AB=d_B−d_A`; range `[0,1]`; `s[V_S]=Σ s_AB`.
- Algorithm 1: Step1 score → Step2 Hungarian one-to-one over `n_T+n_R` external slots (degenerate blocks
  as units; `n_T=3`, `n_R=2 if linear else 3`) → Step3 `|score|≥τ_pure` clean else `mixed_external` (flag
  + dominant slot + `s[V_S]`) → Step4 `s[V_S]≥τ_stretch` stretching / `≤τ_bend` bending / else mixed.
- Conventions: `ε_disp=1e-8`; `unit(0):=0`; degeneracy tolerance groups inertia axes & modes.

### Environment (verified)
- Python 3.13 via `py`; numpy/pandas/scipy/matplotlib/openpyxl available.
  `scipy.optimize.linear_sum_assignment` = Step-2 Hungarian.
- Have: `data/logs/{water,benzene,1grm_MM_UFF}.log`, EMIT for water+benzene, gjf connectivity for all three.
- Excel: sheet `data_score` (per-mode scores + candidate formulas incl. eq:vscore column
  `sum_|d₂-d₁|²/sum(|d₂-d₁|²)*cosθ`); sheet `characterised modes` (reference stretch/bend/irrep labels + lit ref).

---

## Phase 1 — Core classification engine  (data: water, benzene)  — **DO NOW**
- [ ] `src/scoring.py`: add `score_bonds()` exposing per-bond `s_AB` (factor existing `Vscore` loop); assert `Σ s_AB == s[V_S]`.
- [ ] `src/scoring.py`: verify `Rscore` matches eq:rscore (ω formulation, `N_Q` exclusion); fix if divergent.
- [ ] `src/classifier.py` (NEW): inertia degeneracy blocks; `n_T/n_R`; `Thresholds` dataclass (`τ_pure=0.9, τ_stretch=0.9, τ_bend=0.2, ε_disp=1e-8` provisional).
- [ ] `src/classifier.py`: Step 2 global assignment via `linear_sum_assignment` on `|score|`, degenerate blocks as units.
- [ ] `src/classifier.py`: Step 3 purity flag (`mixed_external` + dominant slot + `s[V_S]` annotation).
- [ ] `src/classifier.py`: Step 4 internal split (stretching/bending/mixed) + attach `s_AB`.
- [ ] Validate: reproduce `tab:water` (9 modes, ±1.000 externals) to 3 dp.
- [ ] Validate: benzene EMIT 34–36 flagged `mixed_external`; EMIT 2 vs 9 `s[R]` non-monotonicity reproduced.
- [ ] Output `data/results/<mol>_classified.csv` (scores + label + annotations + `s_AB`).

## Phase 2 — Benzene EMIT stress test (projection reference)  (data: benzene)
- [ ] `src/projection.py` (NEW): EMIT→normal-mode projection `Θ̃=QᵀΘ` (eq:emitproj), stated mass-weighting convention.
- [ ] Reproduce numbers: EMIT 34–36 ≈ 77% translational; EMIT 2 (39% Ry, |s[Ry]|=0.143) vs EMIT 9 (14% Ry, |s[Ry]|=0.215).
- [ ] Figure `fig:benzene`: `s[V_S]` vs freq, and score vs projected NM contribution (highlight EMIT 2/9, 34–36).

## Phase 3 — Hydride library (Excel) + τ calibration + clean-category figures
- [ ] `src/excel_ingest.py` (NEW): read `data_score` + `characterised modes`; map eq:vscore col → `s[V_S]`; join ref labels; emit `data/results/library_scores.csv`.
- [ ] `src/calibrate.py` (NEW): derive `τ_stretch`,`τ_bend` from labeled distributions; `τ_pure` sweep + plateau; write `data/results/thresholds.json`; back-fill `Thresholds`.
- [ ] Figure `fig:confusion` (clean-category confusion matrix + precision/recall vs lit labels).
- [ ] Figure `fig:bondscores` (bond score vs relative Δbond length, ideal vs non-ideal).
- [ ] Figure `fig:boxplots` (freq, Δ|b|, `s[V_S]`; stretch vs bend).
- [ ] Figure `fig:modemixing` (mode score vs averaged Δbond; ideal step vs non-ideal gradient).
- [ ] Figure `fig:sensitivity` (label-change fraction & accuracy vs τ; plateau).
- [ ] Parity check vs Excel `box plots` / `CM` sheets.

## Phase 4 — Gramicidin A: scalability + `s_AB`  (data: 1grm)
- [ ] Run full pipeline on `1grm_MM_UFF.log` (1650 modes) with `1grm.com` connectivity; record wall-clock (fills timing `TODO-DATA`).
- [ ] Output per-bond `s_AB` localization + stretching-character visualization over structure.
- [ ] Comparison-with-projection coverage table + op-count/wall-clock.

## Phase 5 — Orchestration, reproducibility, docs
- [ ] `reproduce.py` (NEW): regenerate every `data/results/*.csv` + `data/figures/*` from inputs.
- [ ] `main.py`: add `--classify`, `--emit-projection`, `--library`, `--figures` subcommands.
- [ ] `README.md`: document classification workflow; resolve **JCE-vs-JCC** mismatch (README cites *J. Chem. Educ.* "paper I"; this is the JCC unified-framework paper).

---

## Verification
- [ ] Water table matches `tab:water` to 3 dp; `Σ s_AB == s[V_S]`; externals reach ±1.000.
- [ ] Benzene EMIT 34–36 flagged + EMIT 2/9 `s[R]` inversion reproduced.
- [ ] Confusion-matrix precision/recall vs lit labels; calibrated τ on the plateau; figures match Excel sheets.
- [ ] Gramicidin run completes with timing + `s_AB`; numbers inserted where `.tex` has `TODO-DATA`.
- [ ] `py reproduce.py` regenerates all CSVs + figures with no manual steps.

## Open items to confirm during execution
- [ ] Exact mass-weighting/projection convention for eq:emitproj (check Excel `Eckart` / `Eckart vs score`).
- [ ] Gramicidin bond-connectivity criterion (`TODO-DATA`): use `1grm.com` topology.
- [ ] Confirm precise `data_score` column equal to `s[V_S]` (candidate `sum_|d₂-d₁|²/sum(|d₂-d₁|²)*cosθ`).
