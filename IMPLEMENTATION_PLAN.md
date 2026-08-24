# JCC Scoring/Classification Program — Implementation Plan & Progress Checklist

> Living checklist. Tick `[x]` as parts are completed; pick up unchecked items in any later session.
> Companion to the JCC manuscript `JCC/JCC_man_scoring/JCC_temp_LaTeXtemplate.tex`, the content plan
> `JCC/JCC_manuscript_structure_JCCformat.pdf` (positioning spine — still current), and the plan slides.
> **Manuscript-structure slides: `JCC/Scoring_Manuscript_Plan_2026-07-05.pdf` is the CURRENT
> authoritative Results & Discussion structure/content doc** (supersedes the `..._2026-07-02/-01/06-29.pdf`
> chain, all kept on disk for history). Point any agent working on Results & Discussion at this file.
>
> **This file was compacted 2026-08-10**: ~2500 lines of dated session-by-session narration (the old
> "RESUME HERE" log and "Changelog") were collapsed into the one-line-per-session **Recent history**
> section at the bottom. Nothing was deleted from the repo — the pre-compaction file (with full narration,
> including every reverted/superseded intermediate decision) is recoverable via
> `git log -p -- IMPLEMENTATION_PLAN.md`. Every decision, formula, and numeric target that was still live
> is preserved below in **Locked decisions** / **Authoritative spec** / the Phase checklists.
>
> Last updated: 2026-08-14 (H2O reverted to mp2/3-21g + VEDA PED added for 8 test molecules — see
> Recent history).

## Status snapshot (current stage)
- **Phases 0–3: DONE.** Core engine, unified classifier, EMIT projection, library ingest, τ-calibration,
  all manuscript figures.
- **Phase 4 (Gramicidin scalability): DEFERRED** to a companion paper — intentional, not a gap.
- **Phase 5 (orchestration/reproducibility/docs): ACTIVE — current frontier.** `reproduce.py`
  orchestrator **built 2026-08-14**. Open: SI Cartesian-geometry export, graphical-TOC image,
  Gaussian revision/year `TODO-DATA`.
- **Phase 6 (strengthen for review): partly open.** 1 of 3 recommended items done (result excluded from
  manuscript by author decision, kept as internal diagnostic); 2 open recommended items are the other
  near-term frontier: mixed-SB bucket validation library-wide, out-of-sample/leave-one-molecule-out
  evaluation. **2026-08-17:** the 18-molecule transferability confusion matrix (T/R/B/SB/S,
  `plot_transferability_confusion`) is now built — partial progress on the out-of-sample item (see its
  checklist entry below); still not a full leave-one-molecule-out re-calibration loop.
- Most recent session (2026-08-17): user finished hand-curating `data/characterised_modes.csv`'s `type`
  column (B/S/SB) for the 8 test molecules that previously lacked it (CH3COCH3, C6H4F2, XeF2Cl2, C2H4,
  C10H8, HOCl, HCOOH, C7H8) — all 18 `test`-tier molecules now have 0 blank `type` cells, unblocking the
  transferability confusion matrix noted as blocked in the 2026-08-14 sessions below. `python main.py
  --library` rerun to pick up the new ground truth into `library_scores.csv` (0 NaN `ref_label` remaining
  on the 18 test molecules' 390 internal rows; no `--calibrate` rerun needed since this tier is outside
  calibration scope). Added `test_tier_molecules()` (`src/library_ingest.py`, mirrors
  `multi_centre_molecules()`) and `plot_transferability_confusion()` (`src/figures.py`, PROPOSED
  `fig:transferabilityconfusion`, wired into `regenerate_all()`) — a pooled joint T/R/B/SB/S confusion
  matrix over all 18 test molecules (n=498), with T/R kept as separate categories (unlike
  `fig:benzeneconfusion`'s collapsed "T/R"). Required adding bare `"T"`/`"R"`/`"B"`/`"S"` entries to
  `REF_LABEL_TO_CATEGORY`/`CATEGORY_LABEL`/`CATEGORY_COLOR` (purely additive, self-mapped to their own
  category rather than routed through the existing `"translation"`/`"rotation"` categories, whose
  `CATEGORY_LABEL` text — "clean translation" — would have been wrong for this table). Verified every
  other figure function still regenerates cleanly (individually invoked, not just via `regenerate_all()`)
  — the one pre-existing `KeyError: 'C2_Tx'` in `plot_benzene_stress_test` is unrelated stale-EMIT-
  projection state in `data/results/C6H6_EMIT.csv` (needs `--emit-projection` rerun), not caused by this
  session's changes; left alone as out of scope. Results: T/R exact (54/54, by construction); B recall
  0.744 (169/227), S recall 0.599 (91/152), SB precision only 0.079 (10/126 — most SB predictions are
  really B or S against literature ground truth), directly corroborating the still-open "validate the
  mixed-SB bucket library-wide" Phase 6 item.
- Prior session (2026-08-14, later same day): H2O.gjf/.log reverted mp2/3-21g* → mp2/3-21g
  (undoing the prior session's rerun — an exact byte revert to the pre-G16-promotion G09 input/output),
  and VEDA4 `.ved`+`.vdf` PED output added for all 9 previously-missing `test`-category molecules (H2O
  + the 8 that had none: CH3COCH3, C6H4F2, XeF2Cl2, C2H4, C10H8, HOCl, HCOOH, C7H8). `reproduce.py` (no
  `--molecules` override needed — `discover_molecules()` found the same 26 `<mol>_normal.csv` files as
  last session) rebuilt everything. **Bug found+fixed in the process:** `data/intermediate/H2O_normal_
  data.txt`'s mtime-based cache didn't invalidate on the revert, because the reverted `H2O.log`/`.com`
  carry an OLDER mtime (2026-07-07, from the original file) than the stale 3-21g*-era intermediate cache
  (2026-08-14) — `_cache_is_fresh()`'s `inter_mtime >= source_mtime` check is blind to a source going
  *backward* in time while staying content-different. First `reproduce.py` pass silently reused the
  stale cache (H2O still showed Tx→Rz=0.0414, freqs 1722.4734/3504.6214/3663.3250, and the `ped` stage
  logged spurious "residual exceeds tolerance" warnings against the newly-added VEDA data, which was
  computed against the *reverted* geometry). Fixed by deleting the stale intermediate file and
  rerunning; confirmed H2O now reverts to the pre-3-21g* values byte-for-byte (Tx→Rz=0.0411, freqs
  1722.4570/3501.5073/3660.7973) with no ped-residual warnings. This is a latent cache-invalidation gap
  (mtime can't detect "reverted to older content") worth a follow-up (e.g. hash-based invalidation), not
  fixed here — noted for a future session. `resync_reference_metadata(write=True)` re-corrected H2O's
  `characterised_modes.csv` freq/k back to the reverted values (irrep/ref_label untouched), clearing the
  label-join warning again. τ_TR/τ_S/τ_B unchanged bit-for-bit (H2O still out of calibration scope). PED
  merge went 10/85 → 19/85 molecules (`combined_ped_vs_scores.csv`); `fig_ped_vs_vscore`/
  `fig_ped_vs_bondscore` are the only figures that materially changed (more `mol_type=='test'` points);
  `fig_cputime*`/`fig_gaussian_nbasis*` changed too but only from run-to-run timing-benchmark noise, not
  from either input change; all other figures byte-identical modulo PDF `CreationDate`/`ModDate`/`ID`
  metadata. The 8 new molecules still lack `characterised_modes.csv` ground-truth rows (only freq/k/μ
  get resynced, not irrep/ref_label) — the 18-molecule transferability confusion-matrix figure remained
  blocked on hand-curation, not on PED data, confirmed unchanged at the time (resolved 2026-08-17, see
  most recent session note above). `pytest` 139/141 green (2 pre-existing
  `test_veda_fmt_regression.py` failures only, unrelated to either input change); 2 golden-value tests
  repinned for H2O's reverted numbers (`tests/test_scores.py::test_water_tab_water` — Vib3 Tx/Rz sign
  flipped back to the pre-G16-promotion convention, `-0.1963`/`-0.2959`; `tests/test_library_ingest.py::
  test_resync_reference_metadata_synthetic_fixture` — fixture's "already exact" mode 2/3 values and
  mode 1's corrected value reverted with the log). Roster-count/exclusion-set golden values from the
  prior session (18-molecule `test` set, 85-molecule roster, 19-name exclusion) untouched, as instructed.
- Prior session (2026-08-14): `test`-category transferability roster expanded 9 → 18 molecules
  (H2O moved `non-ideal` → `test`, rerun at mp2/3-21g*; 8 new molecules added: CH3COCH3, C6H4F2,
  XeF2Cl2, C2H4, C10H8, HOCl, HCOOH, C7H8). Full roster 77 → 85. `reproduce.py --molecules <union>` run
  end-to-end; `resync_reference_metadata(write=True)` corrected H2O's stale
  `characterised_modes.csv` freq/k/μ in place (irrep/ref_label untouched, as designed) — this also
  cleared the "H2O internal rows NOT label-joined" warning the basis-set rerun had introduced.
  τ_TR/τ_S/τ_B **unchanged bit-for-bit** (H2O was never `ideal`, and it's now fully out of calibration
  scope). H2O's own worked-example numbers moved by ≤0.0003 (Tx→Rz 0.0411→0.0414, well within the
  golden test's `5e-4` tolerance) — a real but small basis-set shift, not a code regression. Pooled
  stretch/bend confusion-matrix recall shifted slightly (stretch tp 123→121, n_ref 203→201; bend tp
  193→192, n_ref 209→208) purely from H2O's 3 always-correctly-classified non-ideal modes leaving the
  scope-filtered population — not a threshold or logic change. 8 new molecules have no
  `characterised_modes.csv` rows yet (hand-curation deferred, per standing decision) and no VEDA
  PED output yet — both are open follow-ups, not blockers to this session's pipeline run. `pytest`
  139/141 green (2 pre-existing failures only — fchk/log staleness in `test_veda_fmt_regression.py`);
  9 golden-value tests repinned for the new roster counts (see `tests/test_calibrate.py`,
  `tests/test_flag_validation.py`, `tests/test_library_ingest.py`).
- Prior session (2026-08-11): `mol_list_method.csv` format migration (`basename` column dropped,
  9-molecule `test` category added) adapted through `src/library_ingest.py`/`src/calibrate.py`/
  `src/figures.py`. `pytest` 136/138 green (2 pre-existing failures only — fchk/log staleness in
  `test_veda_fmt_regression.py`; the roster/disk-count-check failures this note used to also cover are
  fixed by this change).

## Context

The JCC manuscript ("A Unified, Reference-Free Framework for Classifying the 3N modes of molecular
motion") is written; the program now implements the full reference pipeline (score → classify →
calibrate → figures) reproducing its tables/figures. Build order was **core engine first**, validated on
data in hand (water, benzene, gramicidin), then scaled out via the full 68-molecule roster
(`data/mol_list_method.csv`), each recomputed from its own `.log`/`.gjf` through the real engine.

### Locked decisions
- **Ingest architecture (FINAL 2026-07-07):** every score column is recomputed from `data/logs/`+`data/gjf/`
  via the real engine (`src/library_ingest.py`, formerly `src/excel_ingest.py`); `src/csv_label_ingest.py`'s
  CSVs (`data/characterised_modes.csv`, `data/ref-label_citation.csv`) supply `ref_label`/`ideal`/`ref_key`
  ONLY (ground-truth ties, never scores). The xlsx workbook (`data/vibrational-scoring-functions.xlsx`) and
  the once-intermediate `data/data_score.csv` are **no longer read by any code path** (both retired
  2026-07-07/09).
- Program **generates figures** (matplotlib), reproducing every manuscript figure.
- Sequence: **core engine first.**
- Validation is **two-tier**: *score-level* checks (threshold-independent) pinned early; *label-level*
  checks (threshold-dependent) pinned only **after τ is frozen** (Phase 3).
- Review-survivability extensions **deferred** to optional Phase 6.
- **Algorithm 1 Step 2**: plain one-to-one `linear_sum_assignment`, **no degenerate-axis block-constraint
  mechanism** (Decision 8, retracted 2026-07-01) — axis-degeneracy is a labeling convention, not an
  assignment ambiguity; T/R references are built directly from geometry (Eckart-Sayvetz), never searched for.
- **EMIT 36 blind spot = DECISION X (LOCKED 2026-06-30):** honest limitation, no new machinery. Scores
  measure directional/geometrical character; EMIT 34/35/36 all have `s[T]=1` (perfect axis alignment) —
  only magnitude variation differs, which encodes 34/35's stretching (caught by `s[V_S]`) vs 36's
  out-of-plane bending (invisible to scores; only projection resolves it). A "rigid-body residual" fix (Y)
  was rejected — it reopens the substitutability objection.
- **Gramicidin/Phase 4 scope (Decision 5, 2026-07-01):** entirely deferred to a companion paper built
  around `s_AB` as a standalone analytical tool; the coverage claim doesn't depend on scale.
- **Manuscript exclusion (2026-07-09):** the 36-mode EMIT flag-confusion precision/recall statistic
  (Phase 6) deliberately excluded from the manuscript — its ground truth is a non-canonical threshold on a
  continuous projection fraction, not rigorous enough to state as a bare accuracy number. Kept as internal
  diagnostic only.
- **Label vocabulary (2026-07-02, author-directed):** classifier labels are axis-specific short codes
  (`Tx`, `Ty*`, `S`, `B`, `SB`, …), not generic constants.
- **CLI/CSV consolidation (2026-08-07, most recent):** one CSV per molecule per mode type
  (`<mol>_{normal,EMIT}.csv`); `--classify` removed (classification always-on); PED-merge/EMIT-projection
  enrich the existing CSV in place (reversal of the prior "never modify the classified CSV" rule).
- **G16 promoted to canonical (LOCKED 2026-08-10):** `data/logs`+`data/gjf`+`data/fchk`+
  `data/characterised_modes.csv` now hold the **Gaussian 16, Revision C.01** rerun (previously
  `data/logs_rerun`/etc.); the original **Gaussian 09** calculation is archived, recoverable via `git mv`
  history, under `data/logs_g09`/`data/gjf_g09`/`data/fchk_g09`/`data/characterised_modes_g09.csv`.
  Promotion followed a clean G09-vs-G16 consistency check (`scripts/compare_rerun.py`,
  `data/results/rerun_consistency_report.csv`, 846 rows): 0 point-group mismatches, 3/846
  threshold-boundary label flips (FBr3 modes 1/2 bend↔SB canceling net, ClH3 mode 4 stretch SB→S), 14
  irrep swaps manually reconciled into `characterised_modes.csv` beforehand. Recalibrating τ against the
  full G16 library moved `τ_S` by ~1.5e-7 (`0.9036817451504533` → `0.9036818966195127`, still XeOH4's
  ideal-stretch min, just under slightly different G16 geometry) — negligible; `τ_TR=0.95` and
  `τ_B=0.17326891344050538` are unchanged bit-for-bit, and the sensitivity plateau range
  `[0.34, 0.995]` is unchanged. H2O and CO2 scores are exact-match (0 delta, up to the usual
  arbitrary-sign/degenerate-axis-labeling convention) between G09 and G16; benzene's EMIT-scored and
  PED-enriched normal-mode CSVs came back byte-identical. Resolves the Phase 5 "Gaussian revision/year
  `TODO-DATA`" item below — the canonical revision is now known and pinned.

### Authoritative spec (from the `.tex`; algorithm = PDF §B6.2 `classify_all_modes`)
- `s[T_Q] = (1/N) Σ unit(d_A)·Q̂`, atoms with `|d_A|>ε_disp`; divisor is **N** (zero-motion atoms dilute via
  `unit(0):=0`); range `[-1,1]`.
- `s[R_Q]` — **CANONICAL = consensus form (FINAL 2026-06-30):**
  `s[R_Q] = (1/(N−N_Q)) Σ (unit(r⊥^A)×unit(d^A))·Q̂`, `r⊥^A = r^A−(r^A·Q̂)Q̂`, `N_Q`=on-axis atoms only.
  Normalize `r⊥` and `d` **separately** then cross (retains `sin φ`, so the score measures *how much* of
  the motion is rotation, down-weighting non-tangential/stretch-like motion) — NOT the ω-form (÷`|ω|`,
  which was tried and reverted as wrong). tab:water: Tx→Rz `0.049`; ν_as→Rz `−0.295`. EMIT 9 `|s[Ry]|=0.215`
  > EMIT 2 `|s[Ry]|=0.143` (inversion, intentional).
- `s[V_S] = (1/Σ w_b|Δb|²) Σ w_AB|Δb_AB|²·|unit(Δb_AB)·b̂_AB^{i}|`, `Δb_AB=d_B−d_A`, **`b̂^i`=INITIAL
  (equilibrium) bond direction**; range `[0,1]`. Per-bond
  `s_AB = w_AB|Δb_AB|²·(unit(Δb_AB)·b̂_AB^i)/Σ_bonds w_b|Δb_b|²` — **SIGNED since 2026-07-30**
  (+stretch/−compress); `s[V_S]` itself non-negative (invariant `Σ|s_AB| == s[V_S]`, not `Σs_AB`).
  `fig:bondscores` plots `|s_AB|`.
- **`w_AB` = the bond weight, NEW 2026-08-14, selectable via `--v-weighting {mu,none}`, default `mu`.**
  `none` = 1 (the original definition). `mu` = the reduced mass `μ_AB = m_A m_B/(m_A+m_B)`, rescaled by
  its max, so a bond counts for as much as the kinetic energy its relative stretching motion carries
  (`½μ|ḃ|²`). `w` appears in numerator AND denominator, so `s[V_S]` stays a weighted mean of per-bond
  `|cos|` — `[0,1]` and `Σ|s_AB|==s[V_S]` both survive by construction. **Key property: μ cancels
  identically when all bonds share one μ, so homoleptic AB_n molecules score BIT-FOR-BIT the same under
  both** (asserted with `==` in `tests/test_scores.py`; `tab:water` is therefore untouched, bond scores
  still `0.4983`/`±0.4992`). Only the 9 heteroleptic library molecules move. Consequences:
  `τ_S 0.9036818966195127 → 0.9012868462693647` (set by XeOH4, the one heteroleptic ideal molecule);
  **`τ_B` UNCHANGED** (set by homoleptic IH3); `τ_TR` and its plateau unchanged; library stretch recall
  `0.6009852 → 0.6059113` (the single homoleptic mover, AsCl3 mode 4, crossed the *lowered* `τ_S` without
  its own score changing); bend row unchanged. Benzene: 6 modes move, all into MIXED — **SB recall
  `0.0 → 1.0`** (the 1532.85 cm⁻¹ E₁u pair rises `0.082/0.091 → 0.279/0.314`), stretch `7/10 → 6/10`,
  bend `16/18 → 13/18`, EMIT stretch calls `7 → 1`. **Manuscript updated 2026-08-14** in
  `JCC_temp_LaTeXtemplate.tex` (the live file; `JCC_man_CA.tex` deliberately left untouched): the
  argument that the 1532.85 cm⁻¹ E₁u pair "should be assigned as bending instead" was reversed, since
  the framework now agrees with the reference's mixed label. SI files updated: rigorous_tier_check
  (τ_S), retention_migration, benzene_precision_recall, computational_cost (eq + **19 → 20 ops/bond**,
  one extra multiply — the weighted square is formed once and reused in numerator and denominator).
  `irrep_coupling` and `sensitivity` needed **no** change (AB₃-only and τ_TR-sweep respectively, both
  verified bit-identical).
- **`fig:ped_vs_vscore` / `fig:ped_vs_bondscore` switched from a quadratic to a straight-line OLS fit
  (2026-08-14, author's call).** The two changes reinforce each other: on the *unweighted* data the
  quadratic was earning its keep (linear `R²=0.9479` vs quadratic `0.9650`, i.e. the `s[V_S]`–`%ν`
  relation was visibly curved), whereas under μ-weighting the relation is essentially straight
  (linear `0.9637` vs quadratic `0.9648` — the curvature term buys `+0.001`). The abstract's and
  §Transferability's `R²` is now `0.96` (n=118 frequency points, 9 test molecules; slope 0.00911,
  intercept 0.04411) — **written into the .tex 2026-08-14**.
- Provenance: `thresholds.json` and `library_scores.csv` carry a `v_weighting` stamp; `classify_all_modes`
  refuses a scorer/threshold mismatch (`Thresholds.bootstrap()`'s `"*"` is the escape used by the first
  `--library` pass after a switch). The pre-2026-08-14 outputs are in `data/{results,figures}/archive_unweighted/`;
  `scripts/compare_weighting.py` diffs the two.
- **Algorithm 1**: Step1 score → Step2 global one-to-one Hungarian (`linear_sum_assignment(...,
  maximize=True)` over `n_T+n_R` external slots vs. all modes, `Σ|score|`) → Step3 **two-gate purity**:
  clean iff `|score|≥τ_TR` AND `s[V_S]≤gate2_bar`, else `MIXED_EXTERNAL_WITH_VIBRATION` (flag + dominant slot +
  `s[V_S]`) → Step4 (unassigned modes only) `s[V_S]≥τ_S` stretching / `≤τ_B` bending / else mixed.
  `n_T=3`, `n_R=2 if linear else 3`. Gate 1 (`τ_TR`) and Step 2's assignment are scheme-independent; gate
  2's bar is **scheme-dependent as of 2026-08 (`gate2_bar` in `src/classifier.py`)**: `τ_purity` (fixed,
  0.05 — decoupled from the three-way split, exploratory like `τ_SB`'s own rollout, not swept) under
  `scheme="binary"` (paper-standard default); `τ_B` (calibrated) under `scheme="threeway"`, unchanged from
  before this split existed. For exact normal-mode externals (`s[V_S]=0`) the gate never fires regardless.
  Benzene EMIT 34/35 (`s[V_S]=0.667/0.577`) correctly flagged under both; EMIT 36 is the Decision-X blind
  spot (see Locked decisions).
- Conventions: `ε_disp=1e-8`; `unit(0):=0`. `DEGEN_TOL` survives only as a general numeric constant (Phase-3
  degenerate-mode sanity check), not an assignment mechanism. `τ_TR` comes from **calibration**, not
  hardcoded.

### Environment (verified)
- Python 3.13 via `py`; numpy/pandas/scipy/matplotlib/openpyxl available.
  `scipy.optimize.linear_sum_assignment` = Step-2 Hungarian (minimizes by default — use `maximize=True`).
- `data/logs/{water,benzene,1grm_MM_UFF}.log` + EMIT (water, benzene) + gjf connectivity, all in repo.

---

## Phase 0 — Pre-build decisions & refactor — **DONE (2026-07-01)**
- [x] Headless pipeline refactor: `main.py` exposes `resolve_dirs()`, `load_inputs()`, `score_modes()`,
      `run_pipeline(mol_name, mode_type, data_dir="data", write=True) -> (df, output_path)` (importable,
      no blocking `input()`); interactive `main()` kept with a `--mode {normal,emit}` bypass. Fail-loud
      raises replace the old catch-all `except`.
- [x] eq:emitproj mass-weighting pinned: `sqrt(mass_A)` per atom, `src/projection.py`-local only —
      confirmed via Gram-matrix orthogonality check (unweighted off-diagonals up to 0.80, mass-weighted
      ≤4e-4, the Eckart-Sayvetz signature). Validated to 3.1e-4 max deviation vs. the pre-existing
      `benzene_EMIT_contributions.csv`.
- [x] Excel `data_score` V_Stretch column verified == engine `s[V_S]`: H2S/SF2 re-scored in-engine, match
      to ~5 sig figs (residual ~1e-5, log-rounding not a formula discrepancy).
- [x] Constants centralized in `scoring.py`: `EPS_DISP=1e-8`, `EPS_NORM=1e-9`, `EPS_DENOM=1e-6`,
      `DEGEN_TOL=1e-3`, `RANGE_TOL=1e-6`.
- [x] `principal_axes()`/`axis_blocks()` accessors added; **update (2026-08 cleanup):** `axis_blocks()`/
      `DEGEN_TOL`-as-assignment-mechanism deleted as confirmed-dead code (Decision 8 already retracted
      their purpose); `principal_axes()` retained (used by `classifier.is_linear()`).

## Phase 1 — Score engine + score-level regression — **DONE (2026-06-30/07-01)**
- [x] `Rscore` consensus form confirmed final (see Authoritative spec).
- [x] Linear-molecule guard verified on CO₂: `s[Rx]=0` (n_R=2), `Ry=Rz=1.000`, no divide-by-zero.
- [x] `score_bonds()` exposes per-bond `s_AB`; range-invariant asserts (`s[T],s[R]∈[-1,1]`, `s[V_S]∈[0,1]`).
- [x] Mode-count fail-loud checks: `GaussianParser.parse` requires 3N−6/3N−5 normal modes;
      `EMITParser.parse` requires 3N EMIT modes.
- [x] Golden-reference regression harness `tests/test_scores.py` (pytest): tab:water to 3 dp,
      `Σ|s_AB|==s[V_S]`, score ranges, CO₂ linear guard, benzene-EMIT targets, parser fail-loud.
- [ ] Remaining fail-loud (low priority, deferred): raise on missing bonds instead of silently scoring
      wrong `V`.

## Phase 2 — Unified classifier + projection reference — **DONE (2026-07-01), 1 item open**
- [x] `src/projection.py`: EMIT→normal-mode projection (`Θ̃=QᵀΘ`, eq:emitproj); emits
      `<mol>_EMIT_contributions.csv` + `<mol>_EMIT_projection_full.csv`; `main.run_projection_pipeline()`.
      Reproduces EMIT 34/35/36 ≈76.7% translational; EMIT 2 (38.7% Ry) vs EMIT 9 (14.1% Ry) inversion.
- [x] `src/classifier.py`: `is_linear()`/`external_slots()`, `Thresholds` dataclass. Step 2 global
      one-to-one Hungarian (maximize `Σ|score|`, no block-constraint — Decision 8). Step 3 two-gate purity
      (Authoritative spec). Step 4 internal split; per-bond `s_AB` attached only for
      `STRETCHING`/`MIXED_STRETCH_BEND` (flagged-external stretching modes, e.g. EMIT 34/35, currently
      carry no per-bond detail — self-consistent with spec as written, open author call for later).
- [x] `<mol>_classified.csv` output for all 4 combos (water/benzene × normal/EMIT); columns Mode,
      Freq/Eigenvalue, Tx..Rz, V_Stretch, label, annotation, s_AB.
- [ ] Figure `fig:benzene` (`s[V_S]` vs freq + score vs projected NM contribution, highlight EMIT 2/9,
      34–36) — caption must frame as flag-behavior/non-monotonicity, NOT a T/R-accuracy benchmark
      (A5 spine guard). Confirm this is wired/current before manuscript finalization.

## Phase 3 — Library ingest + τ-calibration + clean-category figures — **DONE (2026-07-02)**
> Sequence: ingest → classify library → calibrate (freeze τ) → label-level validation → figures.
- [x] `src/library_ingest.py`: recomputes every score column from `.log`/`.gjf` per molecule (see Locked
      decisions) → `data/results/library_scores.csv`.
- [x] Classifier run over the library (`attach_geometry_classification()`), gated on whole-molecule freq
      agreement; ideal-T/R external rows appended independent of that gate (feeds `fig:confusion`).
- [x] `src/calibrate.py`: `τ_S`/`τ_B` from ideal-molecule stretch/bend distributions; `τ_TR` swept
      0.05–0.999 (step 0.005, full curve → `tau_sensitivity_sweep.csv`); **frozen**:
      `τ_TR=0.95, τ_S=0.90368 (0.9036817451504533), τ_B=0.17327 (0.17326891344050538)` →
      `data/results/thresholds.json`. Plateau = longest contiguous zero-label-change run at max accuracy.
- [x] Label-level validation (τ frozen): benzene EMIT 34/35 flagged, EMIT 36 blind spot re-verified.
      `confusion_matrix_stats()` reports **both pooled and ideal/non-ideal-tiered** precision/recall+`n`
      (avoids misreading intrinsic non-ideal mixing as classifier error). 68-molecule final roster:
      T/R precision=recall=1.0; stretch precision=1.0, recall=0.6010 (recall_ideal=1.0,
      recall_nonideal=0.5030); bend precision=1.0, recall=0.9206 (recall_ideal=1.0,
      recall_nonideal=0.8963). `floor_met` (0.95 floor) = False for pooled stretch/bend recall — by
      design, not a defect.
- [x] All 5 figures built: `fig:confusion` (2-tier), `fig:bondscores`, `fig:boxplots`, `fig:modemixing`
      (irrep-degeneracy sub-panel intentionally NOT built — spec unconfirmed), `fig:sensitivity`
      (τ_TR only). `src/figures.py`.
- [x] Parity vs Excel `box plots`/`freq vs score` sheets confirmed (exact match; 2-row-per-group gap on
      `freq vs score` traced to SnO2's Excel-side missing frequency, not a code issue); `CM` sheet is an
      unrelated center-of-mass check, not a confusion-matrix reference.

## Phase 4 — DEFERRED to companion paper (Gramicidin A scalability; out of scope for this manuscript)
> **Scope decision (2026-07-01, Decision 5):** Gramicidin scalability/wall-clock/`fig:gramicidin`,
> transition-state characterization, and isotopic-substitution comparison are ALL deferred to a future
> companion paper built around `s_AB` as a standalone analytical tool. **Why:** the coverage claim ("one
> framework classifies all 3N modes in one pass") doesn't depend on scale — water, the hydride library, and
> benzene (normal + EMIT) fully support it. The manuscript's Conclusion states this explicitly so the
> absence of a large system reads as a decision, not a gap.
>
> **Recorded here, not deleted, so the absence is legible as intentional.** No gramicidin work is active:
> not vectorization, not the `1grm_MM_UFF.log` pipeline run, not `s_AB` localization, not
> `fig:gramicidin`/wall-clock table. `data/logs/1grm_MM_UFF.log` + `data/gjf/1grm.com` stay tracked
> (companion-paper input) but nothing in the current build order consumes them. Revisit only when the
> companion paper begins.

## Phase 5 — Orchestration, reproducibility, docs, submission assets — **ACTIVE (current frontier)**
- [x] `src/figures.py`: one function per figure + shared style helper — all 6 manuscript figures
      implemented, each writes vector PDF + ≥300 dpi PNG. **Not yet wired into an orchestrator.**
- [x] `main.py` subcommands: `--emit-projection`, `--library`, `--calibrate`, `--figures`
      (`--classify` removed 2026-08-07 — classification is always-on).
- [x] `README.md` rewritten: JCC is the sole/first paper (JCE being withdrawn before JCC submission,
      not a "Paper I/Paper II" pair) — zero JCE references remain.
- [x] **`reproduce.py`** — **DONE 2026-08-14.** Regenerates every `data/results/*.csv` +
      `data/figures/*` from inputs, headless, calling the existing entry points in-process.
      Stages, in the order the couplings force: `library → calibrate → library2 → molecules →
      emit → projection → ped → benzene_validation → benchmark → figures` (+ opt-in
      `sync_manuscript`). Flags `--v-weighting/--data-dir/--only/--skip/--molecules/--dry-run`;
      writes `data/results/_run_manifest.json` (weighting, commit, stages, outputs). Two
      non-obvious orderings it encodes: `--calibrate` sits BETWEEN two library passes
      (`predicted_label` is threshold-dependent, the τ derivation is not), and the CPU-time
      benchmark must precede the figures (four read it, and it times *our* classifier).
      `src/benzene_validation.py` is a named stage because it has no `main.py` flag and two
      figures read only its CSVs — skipping it leaves them silently stale.
- [x] **SI Cartesian-geometry export, water/benzene** — **DONE 2026-08-21.** `src/si_geometry_export.py`
      reuses `GaussianParser`'s already-parsed Standard-orientation geometry (no new parsing logic) to
      write `longtable` fragments to `JCC/JCC_SI/tables/tab_cartesian_water.tex`
      (`tab:cartwater`) and `tab_cartesian_benzene.tex` (`tab:cartbenzene`). Alongside it,
      `src/si_tables.py` generates the four full score-table fragments (`tab_water_scores.tex`
      `tab:waterscores`, `tab_library_scores.tex` `tab:libraryscores`, `tab_benzene_normal.tex`
      `tab:benzenenormalscores`, `tab_benzene_emit_full.tex` `tab:benzeneemitscores`) for the SI's
      "Full score tables" section. Both wired as the opt-in `si_tables` stage in `reproduce.py`'s
      `OPTIONAL_STAGES` (writes outside `data/`, same reasoning as `sync_manuscript`).
- [ ] **SI Cartesian-geometry export, CO₂/gramicidin/library** — deferred; same `cartesian_table()`
      helper in `src/si_geometry_export.py` extends directly once those logs are in scope
      (gramicidin is companion-paper scope per Decision 5). B7/B14 reproducibility requirement.
- [ ] **Graphical-TOC image** (B4, submission-required): 50×50 mm, per the structure-doc concept.
- [x] Fill the **Gaussian revision/year** `TODO-DATA` and fix the inconsistent citation — resolved by the
      2026-08-10 G16 promotion (see Locked decisions): canonical revision is **Gaussian 16, Revision
      C.01** (confirmed from `data/logs/H2O.log`'s header). The manuscript `.tex` itself still needs this
      value substituted in and its citation reconciled — that substitution is `lead-author`/manuscript
      work, out of scope for this repo-tracking doc, but the data-side ambiguity this item was blocked on
      is now resolved.

## Phase 6 — Strengthen for review (re-tiered 2026-07-01)
> Re-triaged after Decision 5 (Gramicidin deferred) removed the paper's only scale/robustness
> demonstration, raising the weight the remaining validation carries.

### Recommended before submission
- [x] **Flag precision/recall over all 36 benzene EMIT modes** (`src/flag_validation.py`) — DONE
      2026-07-02: ground truth = projection fraction `M_ext` binarized at 0.05/0.95; result TP=5, FP=0,
      FN=12, TN=19 → precision=1.0, recall≈0.294. Reproduces every named anchor (34/35→TP, 36→FN,
      9→TP vs 2→FN). Finding: FN rate is not isolated to EMIT 36 — Step 2's one-to-one assignment caps
      ever-flaggable modes at `n_T+n_R=6` of 36, so 11 modes with genuine 7–39% external character by
      projection can never be flagged (a second, independent limitation alongside Decision X). Library
      externals (146 rows, 25 molecules): FP=0 confirmed directly. **Excluded from the manuscript
      2026-07-09** (see Locked decisions) — kept as internal diagnostic + tests
      (`tests/test_flag_validation.py`) only.
- [ ] **Validate the mixed-SB bucket library-wide** by irrep-degeneracy + CoM arguments (defends against
      the bending=low-stretch circularity concern; direct evidentiary backbone for the manuscript's
      residual-risk discipline). **Partially discharged for benzene only** (`benzene_validation.py`'s
      `benzene_mixed_bond_diagnostic()`: 5/23≈21.7% of lit-labeled bend modes land in
      `MIXED_STRETCH_BEND`, 0 wrong-category; modes 13/14 & 23/24 detected as near-degenerate D6h E-pairs,
      C-C `s_AB` anti-correlated r=−0.969/−0.997). **Still open:** CoM-softening argument, and extending
      the fraction report to the whole library (pooled `mixed_fraction` exists in
      `confusion_matrix_stats()` but not written up as this item's deliverable).
- [ ] **Out-of-sample / leave-one-molecule-out evaluation** of the τ-calibrated classifier — answers
      "reference-free vs. secretly-fit-τ" objection that the sensitivity plateau alone doesn't.
      **Partially discharged 2026-08-17**: `plot_transferability_confusion()` (`src/figures.py`,
      PROPOSED `fig:transferabilityconfusion`) builds the pooled joint T/R/B/SB/S confusion matrix for
      the 18-molecule `test` tier — genuinely out-of-sample since this tier is excluded from τ_S/τ_B
      calibration by construction, against independent hand-curated `characterised_modes.csv` ground
      truth (n=498: T/R 54/54 exact by construction; B recall 0.744 [169/227], S recall 0.599 [91/152],
      SB precision only 0.079 [10/126] — the mixed-SB bucket over-fires against literature ground truth
      here, corroborating the item above). **Still open:** this is one confusion matrix over the held-out
      tier, not a true leave-one-molecule-out re-calibration loop (τ is still fit once on the
      ideal/non-ideal population and applied unchanged) — the CV-style objection isn't fully answered yet.

### Optional (nice-to-have, not fatal if deferred)
- [ ] N-per-cell + confidence intervals on the confusion matrix (thin bend counts invite a significance
      objection).
- [ ] CoM-conservation evidence (central-atom amplitude/neighbor mass vs. `s[V_S]` degradation, TeH₂ vs.
      Br₂O) backing the "explained, not noisy" difficulty gradient.

---

## Verification
- [ ] **Score-level (Phase 1):** water matches tab:water to 3 dp; `Σ|s_AB|==s[V_S]` (1e-6); externals reach
      `|x|≥0.9995`; ranges hold; CO₂ exercises `n_R=2`; benzene EMIT 34–36 `s[V_S]` and EMIT 2/9 `|s[Ry]|`
      to 3 dp. All frozen as goldens in `tests/test_scores.py`.
- [ ] **Label-level (Phase 3, τ frozen):** benzene EMIT 34/35 flagged `mixed_external`; degenerate mode
      sets get consistent labels as an emergent property of plain Hungarian assignment; confusion-matrix
      precision/recall ≥ floor; calibrated τ on the plateau.
- [x] Figures match Excel `box plots`/`CM` sheets (+ spot-checks) — **DONE 2026-07-02** (see Phase 3).
- [x] `py reproduce.py` regenerates all CSVs + figures with no manual steps — **DONE 2026-08-14**.
      `pytest` = 139 passed / 2 failed, the 2 being pre-existing `test_veda_fmt_regression.py`
      failures unrelated to scoring (a stale "no fchk present" assertion now that `data/fchk/C6H6.fchk`
      exists, and a 4.3e-5 Hessian-reconstruction drift vs `ped/reference/C6H6.fmt`). **Fix those two
      before claiming a fully green suite.**
  (Gramicidin verification target removed 2026-07-01 — deferred to companion paper, Phase 4.)

## Open items to confirm during execution
- [x] eq:emitproj mass-weighting convention (Phase 0) — DONE 2026-07-01.
- [x] `data_score` column == `s[V_S]` (Phase 0) — verified 2026-07-01.
- [x] Library optimized geometries — resolved: author supplies `.log`/`.gjf` on request at the Phase-5 SI step.
- [x] Step-2 objective is maximize `Σ|score|` — implemented via `linear_sum_assignment(..., maximize=True)`.
- [ ] V-score uses the initial bond direction `b̂^i` — implemented; needs the checkbox formally ticked.
- [ ] Degeneracy tolerance numeric value (axes by `λ_i`, modes by freq) — fix and document.
- [x] Flag-criterion two-gate purity — EMIT 36 blind spot resolved as Decision X (locked); no rigid-body
      residual added.
- [ ] Verify the per-mode projected translational % for 34/35/36 (36 is NOT ~100% — it has out-of-plane
      bending despite `s[V_S]=0`; don't assume "77%" applies uniformly).

---

## Recent history
> One line per session, newest first. Full narration for any entry predating 2026-08-10 is in git history:
> `git log -p -- IMPLEMENTATION_PLAN.md`.

- **2026-08-20 (later same day, figures follow-up)** — Closed the one remaining open item from the
  entry directly below: canonical `data/figures/` had still not been regenerated after the `tau_SB`
  0.50→0.42 fix, so it was internally inconsistent with the now fully self-consistent canonical
  `data/results/`. Ran `py reproduce.py --only figures` (the `sync_manuscript` stage — which copies
  into `JCC/JCC_man_scoring/images/` — is opt-in and deliberately left untouched; manuscript images
  stay a hand-synced copy, not live). All 26 figure stems (52 files, PDF+PNG) regenerated. Diffed
  old-vs-new canonical figures: raw byte/PDF comparison is unreliable here (matplotlib's PDF font
  subsetting and `bbox_inches="tight"` sizing are not perfectly deterministic between separate
  `savefig` runs, confirmed by rerunning the figures stage twice back-to-back and finding the
  *decompressed drawing content streams* identical while the embedded font-subset bytes and file-level
  MediaBox/CreationDate metadata still varied), so verification instead used decompressed
  content-stream diffs plus rasterized (PNG) pixel-diff magnitude. Result: the 6 tau_SB-sensitive
  figures (`fig_confusion`, `fig_transferability_confusion`, `fig_benzene_confusion`,
  `fig_benzene_normal`, `fig_modemixing`, `fig_ped_vs_vscore`) moved the most (up to ~23% of pixels
  differing for `fig_benzene_confusion`, where reclassified cells shift categories), while the other
  20 moved by a smaller amount (~1-4% of pixels) consistent with `data/figures/` simply having gone
  stale for longer than just the tau_SB switch (e.g. the CPU-time figures reflect a freshly-timed,
  inherently noisy `cpu_time_benchmark.csv`; `fig_benzene` reflects the same C6H6 EMIT pipeline the
  H2O EMIT staleness fix touched) — not a regeneration error. `data/figures_tau_sb_0.42/` promoted
  from its old "6 regenerated + 20 copied-unchanged" set to a full, byte-identical (`filecmp`-verified)
  mirror of canonical `data/figures/`'s 52 top-level files, mirroring what the results-side entry below
  already did for `data/results_tau_sb_0.42/`; `archive_threeway/`/`archive_unweighted/` excluded for
  the same different-axis reason. Both `_tau_sb_0.42/` READMEs updated to match. Bonus consistency fix
  as a side effect: `data/results/transferability_confusion_misclassified.csv` is written only by the
  figures stage (`plot_transferability_confusion`'s misclassified-rows dump), not by `--skip figures`,
  so it had stayed stale at the 0.50 labeling for C10H8 modes 7/13 (`predicted_label` "B") even after
  the CSV-only pass in the entry below — now correctly "S" at 0.42. `pytest` full suite (`tests/` +
  `vscore-webapp/tests/`): 285 passed, 20 failed, 30 skipped — identical to the baseline recorded in
  the entry below; no new test-visible regressions.
- **2026-08-20 (later same day, follow-up fix)** — The same-day `tau_SB` 0.50→0.42 switch below
  turned out to be two bugs, not one. (1) **Incomplete relabel:** the 9-file cheap relabel missed
  `C6H6_EMIT.csv` and every other per-molecule EMIT/normal CSV — `C6H6_EMIT.csv`'s EMIT 11/32
  (V_Stretch 0.4756/0.4774) were still labeled "B" from the stale 0.50 cutoff. (2) **Deeper root
  cause:** the original switch only hand-edited `data/results/thresholds.json`'s `"tau_SB"` key,
  never the actual source of truth — `Thresholds.tau_SB`'s dataclass default in `src/classifier.py`
  was still 0.50, and `calibrate()` (`src/calibrate.py`) never threaded a `tau_SB` kwarg through to
  `Thresholds(...)`, so it always fell back to that 0.50 default and silently overwrote the whole
  `thresholds.json` (including the hand-added key) on the very next real `--calibrate` run — exactly
  what happened during this session's own `reproduce.py` rerun before the fix. Fixed at the source:
  `Thresholds.tau_SB` default is now genuinely 0.42 (with a docstring explaining tau_SB, unlike
  tau_TR/tau_S/tau_B, has no calibration sweep computing a replacement for it — this field IS the
  sole source of truth); `calibrate()` now writes `"tau_SB"` into every `thresholds.json` it
  produces. Also found `data/results/H2O_EMIT*.csv` on disk but never regenerated by any
  `reproduce.py` run at all (`EMIT_MOLECULES` was hardcoded `["C6H6"]`) — stale from *before* the
  2026-08 binary-scheme-default switch (EMIT 3 was labeled "SB", the old threeway vocabulary).
  Added `"H2O"` to `EMIT_MOLECULES`. Archived the retired tau_SB=0.50 canonical snapshot (133 files,
  git commit `7c18326`, the last commit before the first tau_SB session touched anything) into
  `data/results_tau_sb_0.50/` (new sibling folder, mirrors the `_tau_sb_0.42` naming). Reran
  `py reproduce.py --skip figures` (figures out of scope) + `--only emit projection` (to pick up
  H2O) end-to-end against the fixed code; canonical `data/results/` is now genuinely self-consistent
  at 0.42 (18 files changed: the 9 previously-uncovered files below plus 10 previously-untouched
  per-molecule `*_normal.csv` — B3N3H6, C10H16, C10H8, C2H4, C3O3H6, C4H4, C7H8, CH3COCH3, CHCl3,
  OBr4, each with 1-4 modes crossing bend→stretch — plus `C6H6_EMIT.csv`/`H2O_EMIT.csv`/
  `H2O_EMIT_full_cartesian.csv`; remaining diffs across the other ~55 files are float-repr precision
  noise only, no label changes). Promoted `data/results_tau_sb_0.42/` from its 9-file subset to a
  full 73-file mirror of canonical `data/results/` (byte-identical, confirmed via `filecmp`),
  excluding the unrelated `archive/`/`archive_threeway/`/`archive_unweighted/` subdirs (different
  axis). No golden-value tests needed repinning — the 3 tests flagged as fragile by the incomplete
  first pass, plus `test_gate2_scheme_divergence_synthetic_case`'s stub docstring claim, all still
  pass unchanged, because none of them read the 18 files that actually changed label. One test WAS
  load-bearing on the *old* dataclass default and needed a fix:
  `test_classify_all_modes_binary_scheme_never_produces_sb` used the bare `Thresholds()` constructor
  (implicitly tau_SB=0.50) to test binary-vs-threeway scheme divergence on benzene's Vib 13/14 pair;
  repinned to pass `tau_SB=0.5` explicitly so the test stays about scheme-divergence mechanics, not
  about whatever the current canonical tau_SB happens to be. `pytest` 154/156 in `tests/` (same 2
  pre-existing `test_veda_fmt_regression.py` failures, unrelated Hessian-reconstruction precision);
  full suite (`tests/` + `vscore-webapp/tests/`) 285/305 passed, 20 failed, 30 skipped — matches the
  documented pre-existing baseline (2 + ~18) exactly; the 18 `vscore-webapp` failures are a
  pre-existing S-vs-SB scheme mismatch (that suite's own JS-mirroring implementation still defaults
  to `scheme="threeway"` against a `scheme="binary"` reference) spanning molecules far beyond the 18
  relabeled rows, confirming they predate and are unrelated to this fix.
- **2026-08-20** — Canonical `tau_SB` (binary-scheme stretch/bend cutoff) changed 0.50 → 0.42:
  `data/results/thresholds.json` gained a top-level `"tau_SB": 0.42` key (`Thresholds.calibrated()`
  now reads it instead of falling back to 0.50); 9 canonical `data/results/` CSVs (`library_scores.csv`,
  `C6H6_normal.csv`, `combined_ped_vs_scores.csv`, `benzene_normal_reference_detail.csv`,
  `benzene_normal_reference_summary.csv`, `benzene_worked_examples.csv`,
  `transferability_confusion_matrix.csv`, `transferability_confusion_summary.csv`,
  `transferability_confusion_misclassified.csv`) relabeled at the new cutoff (no re-scoring needed —
  `V_Stretch` doesn't depend on `tau_SB`). Figures NOT regenerated this session (CSV-only task). 3
  golden-value tests repinned (`test_benzene_normal_summary_matches_ad_hoc_session_numbers`,
  `test_confusion_matrix_precision_perfect_recall_explained_by_mixed_bucket`,
  `test_confusion_matrix_ideal_nonideal_recall_split`) — mechanism: benzene's near-degenerate bend pair
  at 1056.3901 cm⁻¹ (modes 13/14) both now cross bend→stretch; library-wide, 7 non-ideal
  literature-stretch modes (NBr3 mode 4, OCl4 modes 7/9, OBr4 modes 8/9, SBr4 modes 8/9) are recovered
  while NBr3 mode 3 (literature bend) crosses the wrong way, so stretch recall reaches a perfect 1.0
  while bend absorbs one non-ideal miss. `pytest` unaffected elsewhere (same ~20 pre-existing failures
  in `vscore-webapp/tests/` and `test_veda_fmt_regression.py`).
- **2026-08-18 (Cartesian-overlap comparison pathway, corrected same day)** — `src/projection.py`
  gained `build_reference_basis_cartesian`/`project_emit_cartesian`, an explicitly non-orthonormal
  counterpart to the locked mass-weighted `build_reference_basis`/`project_emit` (same T/R/vibrational
  reference set, but Q and Theta unit-normalized in the plain Cartesian inner product instead of
  `sum_A m_A (u_A . v_A)`). Requested to show what skipping mass weighting would have produced;
  NOT a replacement for the mass-weighted result (module docstring already documented, pre-existing,
  that the unweighted Gram matrix has off-diagonals up to 0.80 on benzene's real modes — i.e. this
  basis is known not to be orthonormal).
  **Correction:** the first cut squared the overlap (`C2cart_*`) to mirror `project_emit`'s
  `Theta_tilde**2`, but that squaring is only a meaningful "fraction of character" when Q is
  orthonormal (Parseval) — not the case here. Caught and fixed same day: `project_emit_cartesian`
  now reports the raw signed overlap `Q.T @ Theta` itself (a cosine similarity between mode shapes),
  columns renamed `C2cart_*` → `Ocart_*` throughout (`main.py`'s idempotent-merge drop-filter keeps
  catching the old `C2cart_` prefix too, so a stale pre-correction CSV self-heals on the next
  `--emit-projection` run). `Ocart_VS`/`VB`/`VMix` are sums of *signed* overlaps and can partially
  cancel — documented as a coarse diagnostic, not a rigorous decomposition. `Ocart_Sum` has no
  expected target value (contrast `project_emit`'s Parseval-motivated sum~1 check); confirmed on
  benzene EMIT 34: mass-weighted `C2_*` still sums to ~1 as always, raw `Ocart_Tx` = 0.876
  (0.876² ≈ the old, now-removed, `C2cart_Tx` = 0.767), `Ocart_Sum` = 1.89.
  `run_projection_pipeline` (`main.py`) merges `Ocart_Tx..Ocart_Sum` into `<mol>_EMIT.csv` and writes
  `<mol>_EMIT_full_cartesian.csv`. Regenerated `C6H6_EMIT.csv`/`_EMIT_full_cartesian.csv` and
  `H2O_EMIT.csv`/`_EMIT_full_cartesian.csv`. 2 regression tests in `tests/test_projection.py`
  (repinned to the corrected values); full suite 141/143 (2 pre-existing unrelated
  `test_veda_fmt_regression.py` failures from a concurrent session's C6H6.log rerun, not this change).

- **2026-08-24 (Ocart_VS/VB normalized, S+B-only denominator)** — `Ocart_VS`/`Ocart_VB` (`project_emit_cartesian`,
  `src/projection.py`) changed from a plain *signed* sum over each group's reference-mode overlaps (could
  partially cancel, and wasn't bounded to any particular range — benzene `Ocart_Sum` already ran up to ~2.67)
  to `sum(|overlap|) / (S_total + B_total)`, i.e. the same magnitude-weighted-total mechanism `Vscore` itself
  uses (`scoring.py`), not Parseval/squaring. Denominator is **S+B only** — `Mix` is deliberately excluded,
  since `classify_all_modes()`'s default `scheme="binary"` (the 2026-08 binary-classification decision) never
  produces a `MIXED_STRETCH_BEND` reference mode, so `totals[Mix]` is 0 under normal use anyway; `Ocart_VMix`
  itself is still reported (divided by the same S+B-only denominator) for the `scheme="threeway"` case, just
  not part of what the S/B split is normalized against. Result: `Ocart_VS + Ocart_VB == 1` exactly and each
  term ∈[0,1] — same range as `V_Stretch` — for every mode with any internal character; both are 0 (EPS_DENOM
  guard) for a mode that's pure external (no S/B activity to normalize against, e.g. H2O EMIT 1/Ry, EMIT 4/Rx).
  `Ocart_Tx..Ocart_Rz` (raw signed overlap against the external T/R reference axes) are unchanged — already the
  same signed, unnormalized form as `Tscore`/`Rscore`, so already directly comparable. Motivation: `Ocart_*`
  is deliberately kept in plain Cartesian (non-mass-weighted) coordinates specifically so it stays comparable
  to the scores themselves (`Tscore`/`Rscore`/`Vscore` all act on raw Cartesian displacements, not mass-weighted
  ones like `C2_*`'s reference basis) — but the un-normalized signed sum wasn't actually comparable to
  `V_Stretch`'s bounded [0,1] range, which this fixes without touching `Q`'s geometry (no orthogonalization/
  rotation of the non-orthonormal Cartesian basis — that was considered and rejected in favor of this simpler,
  `Vscore`-consistent normalization). Regenerated `C6H6_EMIT.csv`/`H2O_EMIT.csv` (`Ocart_VS`/`VB`/`VMix`/`Sum`
  columns only — `_EMIT_full_cartesian.csv` unchanged, since it stores raw per-reference-mode overlaps).
  `tests/test_projection.py::test_cartesian_overlap_pathway_is_raw_and_not_orthonormal` repinned
  (`Ocart_Sum` for benzene EMIT 34: 1.893708 → 1.875991); new
  `test_ocart_vs_vb_sum_to_one` regression test added. Full suite 156/158 (same 2 pre-existing unrelated
  `test_veda_fmt_regression.py` failures noted above, not this change).

- **2026-08-14 (roster expansion)** — `test`-category transferability set 9 → 18: H2O retagged
  `non-ideal` → `test` (rerun at mp2/3-21g*, new `%mem`/`%nprocshared` header) and 8 new molecules added
  (CH3COCH3, C6H4F2, XeF2Cl2, C2H4, C10H8, HOCl, HCOOH, C7H8) with their own `.log`/`.gjf`/`.fchk` under
  `data/{logs,gjf,fchk}/`. `reproduce.py --molecules <union of discovered + new>` regenerated every
  library/threshold/per-molecule/figure output (no code changes needed — the inclusion-filter design
  from 2026-08-11 handled the new category automatically, as intended).
  `resync_reference_metadata(data_dir='data', write=True)` corrected H2O's stale
  `characterised_modes.csv` freq/k/μ (irrep/ref_label untouched), clearing the `H2O internal rows NOT
  label-joined` warning the rerun had triggered. τ_TR/τ_S/τ_B unchanged bit-for-bit (H2O was never
  `ideal`, and moving it out of `non-ideal` just removes it from calibration scope, same as C6H6 and the
  other test molecules already were). Confusion-matrix stretch/bend pooled recall shifted slightly
  (tp 123→121 stretch, 193→192 bend) purely from H2O's 3 modes leaving the scope-filtered population —
  not a threshold change. 9 golden-value tests repinned for the new 85-molecule/19-name-exclusion/
  506-external-row counts (`tests/test_calibrate.py`, `tests/test_flag_validation.py`,
  `tests/test_library_ingest.py`); one incidental finding while repinning: the historical 7-name
  multi-centre exclusion fixture (`_HISTORICAL_SINGLE_CENTRE_EXCLUDE_7`) contains `'C2H4'`, which is now
  coincidentally also the new transferability molecule's real formula, so `old_dropped` in
  `test_single_centre_only_exclude_matches_scope_decision` legitimately grew from `{C6H6}` to
  `{C6H6, C2H4}` — a name collision, not a logic bug. `pytest` 139/141 (2 pre-existing
  `test_veda_fmt_regression.py` failures only). Still open: hand-curated `characterised_modes.csv` rows
  and VEDA4 PED output for the 8 new molecules (both deliberately deferred, not attempted this session).
- **2026-08-11** — `data/mol_list_method.csv` format migration: the redundant `basename` column dropped
  (roster's `molecule` column now doubles as the on-disk basename) and a 9-molecule held-out
  `mol_type=='test'` transferability set added (CH4, C4H4, C10H16, PCl5, C3H6, B3N3H6, CHCl3, CH3CN,
  C3O3H6; 68 → 77 rows). `src/library_ingest.py` repointed off `basename` onto `molecule` everywhere
  (`load_mol_roster`/`resolve_log_basename`/`_basename_to_molecule_map`/ingest+resync loops), plus a new
  `out_of_calibration_scope_molecules()` (inclusion-filter complement, mol_type not in
  {ideal, non-ideal}) now drives `src.calibrate.SINGLE_CENTRE_ONLY_EXCLUDE` (10 names: C6H6 + the 9 test
  molecules) instead of the old multi-centre-only exclusion, so any future `mol_type` category drops out
  of calibration automatically. `fig:ped_vs_vscore`/`fig:ped_vs_bondscore` scoped to `mol_type=='test'`
  only (excludes C6H6, which appeared there only by accident before). Follow-up same day: CO2 was briefly
  re-tagged `ideal` in `mol_list_method.csv` (diverging from the still-`non-ideal` checked-in
  `library_scores.csv`), then re-tagged back to `non-ideal` to match its linear siblings;
  `library_scores.csv` was regenerated to include the 9 test molecules (68 → 77 molecules, verified
  byte-identical for all 68 pre-existing molecules' scores) and `thresholds.json`/
  `tau_sensitivity_sweep.csv`/`fig:sensitivity` were refreshed (`τ_S`/`τ_B`/`τ_TR` unchanged; the T/R
  ground-truth pool used by the sensitivity sweep grew 404 → 458 external rows since it spans every
  geometry-backed roster molecule, not just the ideal/non-ideal calibration scope — accuracy trace
  unaffected, `label_change_fraction` shifted at exactly one grid point by ~0.00025). Three tests
  (`test_filter_single_centre_library_drops_exactly_the_excluded_molecules`,
  `test_library_external_references_never_false_positive`,
  `test_full_population_has_geometry_for_every_row`, plus
  `test_single_centre_only_exclude_matches_scope_decision`'s golden-comparison block) updated from
  68/404-molecule goldens to 77/458. `pytest` 136/138 (2 pre-existing fchk/log-staleness failures only).
- **2026-08-10** — G16 promoted to canonical (see Locked decisions): G09 archived under `_g09` suffixes;
  `library_scores.csv`/`thresholds.json`/all 19 figures regenerated (`τ_S` shifted ~1.5e-7, `τ_TR`/`τ_B`
  unchanged); 3 golden-value test updates (1 sign-convention flip, 2 from a real 3/846 threshold-boundary
  label flip), each documented with old→new value and cause; `pytest` 132/137 (same 5 pre-existing
  failures as before promotion).
- **2026-08-07** — CLI/CSV consolidation: one CSV per molecule/mode-type, `--classify` removed, PED-merge/
  EMIT-projection now enrich that CSV in place. `pytest` 132/137 (5 pre-existing unrelated failures).
- **2026-08-03** — `scoring.py` vectorized; dead code (`normalize()`, `axis_blocks()`,
  `DEGEN_TOL`-as-mechanism) removed.
- **2026-08-01** — `scripts/benchmark_cpu_time.py` added (CPU-scaling analysis; establishes `scripts/`
  convention).
- **2026-07-30** — Per-bond `s_AB` made signed (+stretch/−compress); `s[V_S]` unchanged; invariant now
  `Σ|s_AB|==s[V_S]`.
- **2026-07-24** — `mol_list_method.csv` roster maintenance (CO2 retag).
- **2026-07-09** — Flag-confusion result (Phase 6) excluded from manuscript (non-canonical ground truth,
  kept as internal diagnostic); `data_score.csv` fully retired.
- **2026-07-08** — `characterised_modes.csv` back-filled by author (+24 molecules); OH4/OF4 removed as
  invalid stationary points (non-ideal confusion → precision=1.0, n=331 modes/58 molecules); all data
  files renamed to match `mol_list_method.csv` molecule names (benzene `_scores.csv` deliberately left
  unrenamed).
- **2026-07-07** — Ingest architecture finalized (see Locked decisions): `library_ingest.py` recomputes
  every score from `.log`/`.gjf`; xlsx workbook retired from all code paths; `mol_list_method.csv`
  (72-row roster + `basename` column) is now canonical.
- **2026-07-05** — `fig:benzeneconfusion` added; `confusion_matrix_stats()` fixed for the literal spec.
- **2026-07-03** — `main.py` subcommands added (`--classify`, `--emit-projection`, `--library`,
  `--figures`, `--calibrate`); README rewritten (JCE/JCC "two-paper" framing retracted — JCC is the sole
  paper, JCE is being withdrawn).
- **2026-07-02** — Phase 2/3 built end-to-end: classifier, projection, library ingest, calibration, all 5
  Phase-3 figures, flag-validation, benzene mixed-bond diagnostic. Label vocabulary renamed to
  axis-specific short codes (`Tx`, `Ty*`, `S`, `B`, `SB`, …).
- **2026-07-01** — Phase 0 refactor; Decision 8 (no degenerate-axis block-constraint mechanism); Decision 5
  (Gramicidin/Phase 4 deferred to companion paper); eq:emitproj mass-weighting locked.
- **2026-06-30** — `s[R]` consensus form finalized (ω-form tried and reverted); two-gate purity design;
  Decision X (EMIT 36 blind spot) locked; golden-reference test harness built; CO₂ linear-molecule guard
  verified.
