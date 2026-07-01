# JCC Scoring/Classification Program — Implementation Plan & Progress Checklist

> Living checklist. Tick `[x]` as parts are completed; pick up unchecked items in any later session.
> Companion to the JCC manuscript `JCC/JCC_man_scoring/JCC_temp_LaTeXtemplate.tex`, the content plan
> `JCC/JCC_manuscript_structure_JCCformat.pdf`, and the plan slides `JCC/Scoring_Manuscript_Plan_2026-06-29.pdf`.
> Last updated: 2026-06-30 (revised after the seven-agent plan review — see Changelog).

> ## ▶ RESUME HERE (session pointer — keep current; update + commit after each increment)
> **Done so far (Phase 0/1):** centralized constants; `principal_axes()`/`axis_blocks()`; `score_bonds()`
> with `Σ s_AB==s[V_S]`; range asserts; `Rscore` = consensus form (settled); manuscript Theory section
> (T/R/V) expanded + `s[R]` eq/`sin φ` discussion. Engine results match tab:water/benzene.
> ~~(5) Phase-0 headless `run_pipeline` refactor~~ ✓ **DONE 2026-07-01** — `main.py` now exposes
> `resolve_dirs()/load_inputs()/score_modes()/run_pipeline()` with no `input()` in that path; fail-loud
> raises on missing bonds/bad mode counts replace the old catch-all `except`; interactive `main()`
> preserved as a thin wrapper (`--mode {normal,emit}` skips the prompt). Verifying it exposed a real bug
> (not axis-frame arbitrariness): `Tscore()`'s `EPS_DISP=1e-8` noise floor was two orders of magnitude
> looser than `Rscore`/`Vscore`'s `EPS_DENOM=1e-6`, so ~1e-8–1e-6 numerical noise on symmetry-required-
> zero atoms in degenerate benzene EMIT eigenvectors (e.g. EMIT 3/4, 13–21) was promoted to full-weight
> unit-vector contributions — previously masked by the old pipeline's incidental `IntermediateIO`
> round-trip (5–6 dp text format crushed the noise to exact zero). Fixed: `Tscore()` now gates on
> `EPS_DENOM` like the other two scores; regression tests added (`test_tscore_ignores_subthreshold_noise`,
> golden EMIT-3 Tz pin). **Decision:** the headless path intentionally never round-trips through
> `IntermediateIO` — that format is only ever used interactively as a bonds-editing hand-off, never a
> numerical filter; the `EPS_DENOM` fix is what now deliberately does the noise-filtering job the
> round-trip used to do by accident.
> **Next, in order:** (6) Phase-0 Excel column verification (re-score one library molecule vs
> `data_score` `s[V_S]` column) → then Phase 2 (classifier/projection) — **but see the 2026-07-01
> scope-revision Changelog entry below first**, which supersedes several Phase 2–4 items per the new
> `JCC/JCC_manuscript_structure_scoped.md`.
> **After each step:** `py -m pytest tests/` should stay green.
> **Resilience rule:** work in small increments; after each, tick the checkbox here + below and `git commit`
> so the plan-in-git always reflects true state. Manuscript `.tex` is outside the repo (not committed).

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
- Sequence: **core engine first.**
- Validation is **two-tier**: *score-level* checks (threshold-independent) are pinned early; *label-level*
  checks (threshold-dependent) are pinned only **after τ is frozen** by calibration (Phase 3).
- Review-survivability extensions are **deferred** to an optional Phase 6 (decided after the core works).

### Authoritative spec (from the `.tex`; algorithm = PDF §B6.2 `classify_all_modes`)
- `s[T_Q] = (1/N) Σ unit(d_A)·Q̂`, atoms with `|d_A|>ε_disp`; **divisor is N** (zero-motion atoms stay in
  the count and dilute via `unit(0):=0`); range `[-1,1]`.
- `s[R_Q]` — **CANONICAL = consensus form (FINAL 2026-06-30):**
  `s[R_Q] = (1/(N−N_Q)) Σ (unit(r⊥^A)×unit(d^A))·Q̂ = (1/(N−N_Q)) Σ (r⊥×d)_Q/(|r⊥|·|d|)`, with
  `r⊥^A = r^A−(r^A·Q̂)Q̂`. `N_Q` = on-axis atoms only. **Rationale (author):** normalize `r⊥` and `d`
  SEPARATELY, then cross — the unit-vector cross product has magnitude `sin φ` (φ = angle between r⊥ and d),
  which is 1 only for purely tangential (ideal-rotation) motion and `<1` as `d` tilts toward radial. Retaining
  `sin φ` (i.e. dividing by `|r⊥||d|`, NOT by `|ω|=|r⊥×d|`) makes `s[R]` measure *how much* of the motion is
  rotation about Q, down-weighting non-tangential in-plane (stretching-like) motion. This is the original
  `scoring.py` form and matches tab:water. **History:** I briefly switched the code to the ω-form
  (÷`|ω|`, which discards `sin φ`) — that was WRONG and is reverted; code restored byte-identical to original.
  JCC `eq:rscore` rewritten to the separate-normalization form + a new paragraph on the `sin φ` factor (the
  point JCE left implicit). tab:water unchanged from original (Tx→Rz `0.049`; ν_as→Rz `−0.295`); EMIT 9
  `|s[Ry]|=0.215`.
- `s[V_S] = (1/Σ|Δb|²) Σ |Δb_AB|²·|unit(Δb_AB)·b̂_AB^{i}|`, `Δb_AB=d_B−d_A`; **`b̂^{i}` = the INITIAL
  (equilibrium-geometry) bond direction**, not the perturbed one; range `[0,1]`. Per-bond
  `s_AB = |Δb_AB|²·|unit(Δb_AB)·b̂_AB^{i}| / Σ_bonds|Δb|²` (global denominator) with `s[V_S]=Σ s_AB`.
- **Algorithm 1** (PDF §B6.2): Step1 score → Step2 global one-to-one assignment over `n_T+n_R` external
  slots **MAXIMIZING `Σ|score|`** (`linear_sum_assignment` minimizes → negate the cost / use
  `maximize=True`), degenerate axis-blocks bound to degenerate mode-blocks **collectively** (not per
  element) → Step3 **clean iff `|score|≥τ_pure` AND `s[V_S]≤τ_bend`** (directionally aligned AND
  internally rigid), else `mixed_external` (flag + dominant slot + `s[V_S]`) →
  Step4 `s[V_S]≥τ_stretch` stretching / `≤τ_bend` bending / else mixed. `n_T=3`, `n_R=2 if linear else 3`.
  **Two-gate purity (refined 2026-06-30):** `s[T]/s[R]` are direction-only (read ±1 for amplitude-varying
  impure modes); the second gate uses `s[V_S]` to catch **stretching-type** external impurity. Benzene
  EMIT: 34/35 (E₁ᵤ, Tx/Ty=1, `s[V_S]=0.667/0.577`) correctly flagged. **KNOWN BLIND SPOT:** EMIT 36 (A₂ᵤ,
  Tz=1, `s[V_S]=0`) is out-of-plane **bending**-mixed (amplitude variation ⊥ the in-plane bonds → bending,
  no stretching), yet its score signature {Tz=1, rest 0, V_S=0} is **identical to a pure z-translation** —
  {s[T],s[R],s[V_S]} cannot distinguish them (no bending observable). The flag detects stretching-type
  impurity only; bending-type is invisible (projection resolves it). Reuses τ_pure/τ_bend; for exact
  normal-mode externals (genuinely rigid, `s[V_S]=0`) the gate never fires (completeness intact).
  **DECISION (LOCKED 2026-06-30): X — honest limitation, no new machinery.** Framing: the scores measure
  **directional/geometrical** character; 34/35/36 genuinely have large translational character (every
  displacement vector aligns perfectly with a Cartesian axis, `s[T]=1`) — only their MAGNITUDES differ.
  Equal magnitudes ⇒ pure translation; magnitude variation encodes the internal residual (stretching for
  34/35 → caught by `s[V_S]`; out-of-plane bending for 36 → invisible). The score is "partially true" — it
  correctly reports 36's dominant translational character; the magnitude-encoded bending residual is left
  to projection. 36 = the worked example of the score/projection boundary. (Y, a rigid-body residual, was
  rejected: it is projection onto the T/R subspace and reopens the substitutability objection.)
- Conventions: `ε_disp=1e-8`; `unit(0):=0`; degeneracy tolerance groups inertia axes (by `λ_i`) & modes
  (by frequency/eigenvalue) — **fix the numeric tolerance** (provisional `1e-3` relative; confirm).
  `τ_pure` is taken from **calibration**, NOT hardcoded (PDF example uses 0.95).

### Environment (verified)
- Python 3.13 via `py`; numpy/pandas/scipy/matplotlib/openpyxl available.
  `scipy.optimize.linear_sum_assignment` = Step-2 Hungarian (remember: minimizes by default).
- Have: `data/logs/{water,benzene,1grm_MM_UFF}.log`, EMIT for water+benzene, gjf connectivity for all three.
- Excel: sheet `data_score` (per-mode scores + candidate formulas incl. eq:vscore column
  `sum_|d₂-d₁|²/sum(|d₂-d₁|²)*cosθ`); sheet `characterised modes` (reference stretch/bend/irrep labels + lit ref).

---

## Phase 0 — Pre-build decisions & refactor (gating; DO FIRST)
- [x] **Headless pipeline refactor — DONE 2026-07-01.** `main.py` exposes `resolve_dirs()`,
      `load_inputs()`, `score_modes()`, `run_pipeline(mol_name, mode_type, data_dir="data", write=True)
      -> (df, output_path)` — importable, no `input()`/blocking ENTER in that path. `main()` keeps the
      interactive prompt + bonds-editing fallback, now with a `--mode {normal,emit}` flag to skip it.
      Fail-loud raises (missing bonds, bad mode counts) replace the old catch-all `except`. Surfaced
      and fixed a real bug in the process — see Changelog.
- [ ] **Pin the projection / mass-weighting convention** (eq:emitproj) as a locked decision *before*
      any projection or comparison work: scores use **unweighted** Cartesian displacements; state how the
      C2 normal-mode reference is mass-weighted so score-vs-projection is like-with-like (check Excel
      `Eckart` / `Eckart vs score`). Every `fig:benzene` and coverage number depends on this.
- [ ] **Verify the Excel column identity:** re-score ONE library molecule in-engine and confirm it equals
      the `data_score` `s[V_S]` candidate column (`sum_|d₂-d₁|²/sum(|d₂-d₁|²)*cosθ`) to 3 dp, before
      trusting the whole ingest/calibration chain.
- [x] **Library geometries: RESOLVED (2026-06-30).** Not a blocker — the author will drop the relevant
      `.log`/`.gjf` files into `data/logs/` and `data/gjf/` on request when the SI-geometry export step
      (Phase 5) needs them. Ask for them at that point.
- [x] **Centralize constants** (`ε_disp=1e-8`, degeneracy `tol`) — done as named module-level constants in
      `scoring.py` (`EPS_DISP=1e-8`, `EPS_NORM=1e-9`, `EPS_DENOM=1e-6`, `DEGEN_TOL=1e-3`, `RANGE_TOL=1e-6`);
      removed the scattered `1e-6`/`1e-9`. Tscore cutoff moved 1e-6→1e-8 (water/benzene outputs unchanged).
      Rscore left untouched. (dataclass deferred to the classifier's `Thresholds`.)
- [x] **Expose inertia data:** added `principal_axes()` / `axis_blocks()` accessors to `scoring.py`
      (share `_build_inertia_tensor()` with `MIT`); classifier reads moments + axis-degeneracy blocks.

## Phase 1 — Score engine + score-level regression  (data: water, benzene, CO₂)  — **DO NOW**
- [x] **`Rscore` = consensus form, confirmed/kept (FINAL 2026-06-30).** `s[R_Q]=(1/(N−N_Q)) Σ
      (unit(r⊥)×unit(d))·Q̂ = (r⊥×d)_Q/(|r⊥||d|)` — separate `r⊥`/`d` normalization retains `sin φ`
      (down-weights non-tangential in-plane motion). I briefly switched to the ω-form (÷`|ω|`) and reverted;
      code is byte-identical to the original. Externals ±1; ranges hold; EMIT 2-vs-9 inversion `|Ry| 0.143
      vs 0.215`. JCC `eq:rscore` rewritten to the separate-normalization form + new `sin φ` paragraph.
- [x] **Linear-molecule guard — VERIFIED on CO₂ (2026-06-30).** `co2_mp2_3-21g` runs cleanly: molecular-axis
      rotation `s[Rx]=0` (n_R=2), `Ry=Rz=1.000`, no divide-by-zero; 2 stretches `V_S=1.000`, 2 degenerate
      bends `V_S=0`. Result in `data/results/co2_mp2_3-21g_normal_scores.csv`. Automated regression → harness below.
- [x] `src/scoring.py`: added `score_bonds()` exposing per-bond `s_AB` (factored `Vscore` loop into
      `_bond_contributions()`); asserts `Σ s_AB == s[V_S]` (tol 1e-6; observed err ≤2.2e-16). Vscore value unchanged.
- [x] **Range-invariant asserts:** `s[T],s[R] ∈ [−1,1]`; `s[V_S] ∈ [0,1]` — enforced in `calculate_scores`.
- [x] **Mode-count checks (fail-loud) — DONE 2026-06-30.** `GaussianParser.parse` raises unless normal
      modes = 3N−6 or 3N−5; `EMITParser.parse` raises unless EMIT modes = 3N (and on malformed matrix size).
      Verified non-breaking on water/benzene/CO₂ (normal+EMIT), results unchanged.
- [ ] **Remaining fail-loud (deferred to headless refactor):** raise on missing bonds rather than silently
      scoring wrong `V`; replace the catch-all `except` swallow in `main` with explicit errors.
- [x] **Add a linear molecule (CO₂)** to exercise the `n_R=2` / `N_Q` branch (otherwise unexercised by
      water/benzene). Done — `co2_mp2_3-21g` log+com in repo, runs clean (see Linear-molecule guard above).
- [x] **Golden-reference regression harness — DONE 2026-06-30.** `tests/test_scores.py` (pytest +
      standalone runner; `pytest` in `requirements-dev.txt`). 6 tests, all green: `tab:water` to 3 dp,
      `Σ s_AB == s[V_S]`, score ranges, CO₂ linear (n_R=2), benzene-EMIT targets, parser fail-loud. Run
      `py -m pytest tests/` or `py tests/test_scores.py`.
- [x] **Score-level benzene EMIT checks — DONE (in the harness).** `test_benzene_emit_targets`: EMIT 34–36
      `s[V_S]=0.667/0.577/0`; EMIT 2 `|s[Ry]|=0.143` < EMIT 9 `|s[Ry]|=0.215` (inversion). Labels = Phase 3.
- [ ] Output `data/results/<mol>_scores.csv` (unchanged format) + frozen goldens.

## Phase 2 — Unified classifier + projection reference  (data: water, benzene)
- [ ] `src/projection.py` (NEW): EMIT→normal-mode projection `Θ̃=QᵀΘ` (eq:emitproj) using the Phase-0
      convention; **emit a projection-coefficients data file** (per-mode projected NM contribution), not
      just numbers. Reproduce: EMIT 34–36 ≈ 77% translational; EMIT 2 (39% Ry) vs EMIT 9 (14% Ry).
- [ ] `src/classifier.py` (NEW): inertia degeneracy blocks; `n_T/n_R`; `Thresholds` (provisional
      `τ_pure=0.95, τ_stretch=0.9, τ_bend=0.2` — τ_pure to be replaced by calibration in Phase 3).
- [ ] Step 2 global assignment via `linear_sum_assignment` **maximizing `Σ|score|`** (negate / `maximize=True`),
      with a **concrete, documented block-constrained mechanism** so a degenerate axis-block binds a
      degenerate mode-block collectively (plain 1-to-1 will mis-assign Tₐ/Oₕ tops & degenerate pairs).
- [ ] **Step 3 two-gate purity (see Authoritative spec):** clean external iff `|s_slot|≥τ_pure` AND
      `s[V_S]≤τ_bend`; else `mixed_external` (flag + dominant slot + `s[V_S]`). Reuses existing constants.
      Correctly flags 34/35 (stretching-mixed). **BLIND SPOT — EMIT 36** (A₂ᵤ out-of-plane bending,
      `s[V_S]=0`) has a score signature identical to a pure z-translation, so the scores cannot flag it
      (no bending observable). **Manuscript refinement (B8.4):** "EMIT 34–36 flagged" → "34/35 flagged
      (stretching); 36 is the honest limit — its out-of-plane bending residual is invisible to the
      reference-free scores and resolved only by projection." **DECISION X locked** (honest limitation;
      no rigid-body residual). Discussion must note the scores capture genuine directional/geometrical
      translational character — only magnitudes differ — so the label is correct about dominant character.
- [ ] Step 4 internal split (stretching/bending/mixed) + attach `s_AB`.
- [ ] Output `data/results/<mol>_classified.csv` (scores + label + annotations + `s_AB`).
- [ ] Figure `fig:benzene`: `s[V_S]` vs freq, and score vs projected NM contribution (highlight EMIT 2/9,
      34–36). **Caption/plot as flag-behavior & s[R] non-monotonicity — NOT a T/R-accuracy benchmark**
      (A5 spine guard: do not render it as a parity/accuracy plot).

## Phase 3 — Library ingest + τ-calibration (freeze τ) + clean-category figures
> Sequence within phase: **ingest → classify library → calibrate (freeze τ) → label-level validation → figures.**
- [ ] `src/excel_ingest.py` (NEW): read `data_score` + `characterised modes`; map eq:vscore col → `s[V_S]`;
      **also emit per-bond `s_AB`, frequency, relative/averaged Δbond-length, an ideal/non-ideal tag, and
      T/R reference labels** (not just stretch/bend) → `data/results/library_scores.csv`. Confirm the
      library covers exactly the `tab:ideal`/`tab:nonideal` molecules (else a manuscript table edit is forced).
- [ ] **Run the classifier over the library** to attach *predicted* labels (incl. clean-T/R cells) —
      required for `fig:confusion` (reference labels alone can't form the matrix).
- [ ] `src/calibrate.py` (NEW): derive `τ_stretch`,`τ_bend` from labeled distributions; `τ_pure` sweep +
      plateau; **persist the full τ-sweep curve** (label-change fraction + accuracy per τ) for
      `fig:sensitivity`; write `data/results/thresholds.json`; back-fill `Thresholds`. Define "plateau"
      quantitatively (give a numeric criterion).
- [ ] **Label-level validation (now τ is frozen):** benzene EMIT 34–36 flagged `mixed_external`;
      confusion-matrix precision/recall vs lit labels (set a numeric acceptance floor); re-pin label
      goldens. Add a **degenerate-block consistency** check (e.g. benzene EMIT 1–9 share eigenvalue → same
      slot/label).
- [ ] Figure `fig:confusion` (clean-category confusion matrix + precision/recall vs lit labels).
- [ ] Figure `fig:bondscores` (bond score vs relative Δbond length, ideal vs non-ideal).
- [ ] Figure `fig:boxplots` (freq, Δ|b|, `s[V_S]`; stretch vs bend).
- [ ] Figure `fig:modemixing` (mode score vs averaged Δbond; ideal step vs non-ideal gradient).
- [ ] Figure `fig:sensitivity` (label-change fraction & accuracy vs τ; plateau).
- [ ] Parity check vs Excel `box plots` / `CM` sheets (and spot-check the other figures' source data).

## Phase 4 — Gramicidin A: scalability + `s_AB`  (data: 1grm)
- [ ] **Budget a NumPy vectorization pass** on Step-1 scoring + `MIT` mode rotation (≈550 atoms × 1650
      modes are currently pure-Python loops); measure early so the wall-clock `TODO-DATA` is a real number.
- [ ] Run full pipeline on `1grm_MM_UFF.log` (1650 modes) with `1grm.com` connectivity; record wall-clock.
- [ ] Output per-bond `s_AB` localization (e.g. carbonyl C–O, Mode 1316).
- [ ] **Figure `fig:gramicidin`** (explicit task): (a) mode scores vs ascending frequency; (b) Mode-1316
      carbonyl C–O `s_AB` localization over structure. Label as **demonstration**, not validated accuracy.
- [ ] Comparison-with-projection coverage table + op-count/wall-clock. **Keep Phase 4 strictly
      scale/cost + `s_AB`** — no label-correctness/accuracy claim (MM/UFF has no ground truth).

## Phase 5 — Orchestration, reproducibility, docs, submission assets
- [ ] `reproduce.py` (NEW): regenerate every `data/results/*.csv` + `data/figures/*` from inputs, headless
      (uses the Phase-0 `run_pipeline`). Wire in the figure module.
- [ ] `src/figures.py` (NEW): one function per figure + a **single shared style helper**; create
      `data/figures/`; output **vector PDF** (LaTeX embed) **+ ≥300 dpi PNG** preview per figure.
- [ ] **SI Cartesian-geometry export** step (optimized coords for water/benzene/CO₂/gramicidin from logs,
      + library if retrievable) — B7/B14 reproducibility requirement.
- [ ] **Graphical-TOC image** (B4, submission-REQUIRED): 50×50 mm; assemble per the structure-doc concept.
- [ ] Fill the **Gaussian revision/year** `TODO-DATA` (line 463) and correct the inconsistent citation.
- [ ] `main.py`: add `--classify`, `--emit-projection`, `--library`, `--figures` subcommands.
- [ ] `README.md`: document classification workflow; resolve **JCE-vs-JCC** mismatch (README cites
      *J. Chem. Educ.* "paper I"; this is the JCC unified-framework paper).

## Phase 6 — Strengthen for review (OPTIONAL / DEFERRED; decide after core works)
> From expert-reviewer-jcc. Not required to reproduce the manuscript, but directly defends load-bearing
> claims if a referee pushes. Revisit once Phases 1–5 produce data.
- [ ] **Flag precision/recall over ALL 36 benzene EMIT modes** (and library externals) against the C2
      projection reference — systematic flag-fidelity, vs the current anecdotal EMIT 2/9/34–36 (A4).
- [ ] **Out-of-sample / leave-one-molecule-out** evaluation of the τ-calibrated classifier — answers the
      "reference-free vs trained-τ" objection that the sensitivity plateau alone does not.
- [ ] **CoM-conservation evidence** (correlate central-atom amplitude / neighbor mass vs `s[V_S]`
      degradation, TeH₂ vs Br₂O) backing the B8.3 "explained, not noisy" difficulty gradient.
- [ ] **N-per-cell + confidence intervals** on the confusion matrix (thin bend counts invite a
      significance objection).
- [ ] **Validate the mixed-SB bucket** by irrep-degeneracy + CoM arguments; report the fraction of
      lit-labeled modes landing in "mixed" (defends against the bending = low-stretch circularity concern).

---

## Verification
- [ ] **Score-level (Phase 1):** water matches `tab:water` to 3 dp; `Σ s_AB == s[V_S]` (1e-6); externals
      reach ±1.000 (`|x|≥0.9995`); ranges hold; CO₂ exercises `n_R=2`; benzene EMIT 34–36 `s[V_S]` and
      EMIT 2/9 `|s[Ry]|` to 3 dp. All frozen as goldens.
- [ ] **Label-level (Phase 3, τ frozen):** benzene EMIT 34–36 flagged `mixed_external`; degenerate blocks
      assigned consistently; confusion-matrix precision/recall ≥ floor; calibrated τ on the plateau.
- [ ] Figures match Excel `box plots`/`CM` sheets (+ spot-checks for the rest).
- [ ] Gramicidin run completes with timing + `s_AB`; numbers inserted where `.tex` has `TODO-DATA`.
- [ ] `py reproduce.py` regenerates all CSVs + figures with no manual steps; `pytest` green.

## Open items to confirm during execution
- [ ] (Phase 0) eq:emitproj mass-weighting convention — pinned pre-Phase-2.
- [ ] (Phase 0) `data_score` column == `s[V_S]` — verified by in-engine re-score.
- [x] (Phase 0) library optimized geometries — RESOLVED: author supplies `.log`/`.gjf` on request at the Phase-5 SI step.
- [ ] Step-2 objective is **maximize** `Σ|score|` (not scipy's default minimize).
- [ ] V-score uses the **initial** bond direction `b̂^{i}`.
- [ ] Degeneracy tolerance numeric value (axes by `λ_i`, modes by freq) — fix and document.
- [ ] Gramicidin bond-connectivity criterion (`TODO-DATA`): use `1grm.com` topology.
- [x] Flag-criterion: two-gate purity (`|s|≥τ_pure` AND `s[V_S]≤τ_bend`) flags stretching-type impurity (34/35). EMIT 36 out-of-plane bending blind spot RESOLVED as **Decision X** (honest limitation; scores report genuine directional/geometrical translational character, magnitudes differ; bending residual → projection). No rigid-body residual added.
- [ ] Verify in Phase 2 the per-mode projected translational % for 34/35/36 (36 is NOT ~100% — it has out-of-plane bending despite `s[V_S]=0`; do not assume the "77%" applies uniformly).

---

## Changelog
- **2026-07-01 — Headless refactor landed; EPS_DISP/EPS_DENOM Tscore bug found + fixed.** Completed the
  Phase-0 `run_pipeline` refactor left mid-verification at the end of the prior session. Verifying it
  (normal-mode CSVs unchanged; benzene/water EMIT CSVs shifted by more than rounding in a few rows, e.g.
  EMIT 3 Tz `-0.1667→-0.6667`) traced to a real latent bug, not the degenerate-axis-frame-arbitrariness
  hypothesis the prior session was chasing: `Tscore()` gated on `EPS_DISP=1e-8`, two orders of magnitude
  looser than `EPS_DENOM=1e-6` used by `Rscore()`/`Vscore()` for the same "is this atom moving" test. The
  old always-round-trip-through-`IntermediateIO` pipeline (5–6 dp text format) incidentally crushed
  ~1e-8–1e-6 Gaussian-EMIT-file noise (present on symmetry-required-zero atoms inside degenerate
  eigenvalue blocks, e.g. EMIT 3/4, 13–21) to exact zero; the new full-float64-precision headless path
  let that noise leak through `Tscore` as full-weight unit-vector contributions (confirmed by direct
  reproduction of both code paths — confined entirely to Tx/Ty/Tz columns and to degenerate blocks,
  exactly as observed, because `Rscore`/`Vscore` already used the stricter `EPS_DENOM` floor).
  `axis_blocks()`/`principal_axes()` are confirmed dead code (not called anywhere in the scoring path),
  ruling out an axis-choice explanation. **Fix:** `Tscore()` now gates on `EPS_DENOM` instead of
  `EPS_DISP`; added `test_tscore_ignores_subthreshold_noise` (synthetic reproduction) and a golden pin on
  benzene EMIT 3's Tz. Regenerated `benzene_EMIT_scores.csv`/`water_EMIT_scores.csv`; diff against the
  prior commit is now sub-0.0001 rounding drift only. **Decision (documented in RESUME HERE):** the
  headless path deliberately never round-trips through `IntermediateIO`; that format remains solely an
  interactive bonds-editing hand-off.
- **2026-06-30 (final) — s[R] = consensus form is canonical (ω-form reverted).** Author's settled
  reasoning: normalize `r⊥` and `d` separately and cross them; the unit-vector cross product has magnitude
  `sin φ`, which must be retained (÷`|r⊥||d|`, not ÷`|ω|`) so the score reflects *how tangential* the motion
  is. The ω-form (÷`|ω|`) discards `sin φ` and was wrong. Reverted `Rscore` to the original (byte-identical
  results: tab:water Tx→Rz `0.049`, ν_as→Rz `−0.295`; EMIT 9 `|s[Ry]|=0.215`). Rewrote JCC `eq:rscore` to the
  separate-normalization form and ADDED a `sin φ` discussion paragraph (the point left implicit in JCE).
  Supersedes the "(late)" entry below.
- **2026-06-30 (late) — s[R] = ω-form is canonical; code fixed.** After discussion, the author confirmed
  the JCE ω-form `unit(r⊥×d)·Q̂` is canonical: rotation = circulation ⊥ r⊥; the cross product rightly
  annihilates radial (breathing) motion, which `s[V_S]` owns; breathing+swirl modes are caught by the
  two-gate purity. The prior code (`|r⊥||d|` normalization) was a bug. Rewrote `Rscore` to the ω-form,
  regenerated water/benzene results, updated JCC `eq:rscore` and tab:water (Tx→Rz 0.049→0.333; ν_as→Rz
  −0.295→+0.333) and EMIT 9 |s[Ry]| 0.215→0.233. Externals/`s[V_S]`/EMIT 2-vs-9 inversion preserved.
  Supersedes the earlier (wrong) "consensus form is canonical" entry below.
- **2026-06-30 — s[R] definition resolved.** formula-auditor + author confirmed the CODE is correct (the
  consensus form `unit(d)·unit(Q̂×r)`, = the JCE-manuscript definition, matches tab:water). The JCC `.tex`
  `eq:omega`/`eq:rscore` (ω-form) is the transcription error and will be corrected to match (lead-author);
  **no code change** to the Rscore formula. Retained only the linear-molecule divide-by-zero guard as a code task.
- **2026-06-30 — flag-criterion + geometry decisions.** Resolved the two open questions from the review:
  (1) library geometries are not a blocker (author supplies on request at the Phase-5 SI step);
  (2) Step-3 purity is now **two-gate** (`|s_slot|≥τ_pure` AND `s[V_S]≤τ_bend`), grounded in the real
  benzene EMIT CSV (A₂ᵤ EMIT 36 clean vs E₁ᵤ 34/35 flagged) — reuses existing constants, keeps
  completeness for normal modes. Noted the B8.4 manuscript refinement ("34/35 flagged, 36 clean") and a
  Phase-2 check on the 77%-translational figure.
- **2026-06-30 — seven-agent plan review.** Added Phase 0 (headless refactor, pinned projection
  convention, Excel/geometry verification, centralized constants, inertia accessors). Split validation
  into score-level (early) vs label-level (after τ freeze). Corrected algorithm spec: Step-2 **maximize**
  `Σ|score|`; block-constrained assignment; initial bond direction `b̂^{i}`; `ε_disp` gates `|ω|`;
  `τ_pure` from calibration. Moved `Rscore` fix to Phase-1 first task; added linear-molecule (CO₂) test
  + golden-reference pytest harness + completeness checks. Augmented `excel_ingest` outputs (s_AB, freq,
  Δbond, ideal/non-ideal tag, T/R labels) + library classification for `fig:confusion`; `calibrate`
  persists the τ-sweep. Promoted `fig:gramicidin` to an explicit figure task; added shared `figures.py`
  + PDF/PNG output, graphical-TOC, SI geometry export, Gaussian rev/year. Flagged `fig:benzene` spine
  risk. Added optional Phase 6 (review-strengthening). Marked `tab:ideal/nonideal/analogy`, `fig:flowchart`
  as intentionally codeless (no plan item; correct).
