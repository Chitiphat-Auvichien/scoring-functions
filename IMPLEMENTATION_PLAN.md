# JCC Scoring/Classification Program — Implementation Plan & Progress Checklist

> Living checklist. Tick `[x]` as parts are completed; pick up unchecked items in any later session.
> Companion to the JCC manuscript `JCC/JCC_man_scoring/JCC_temp_LaTeXtemplate.tex`, the content plan
> `JCC/JCC_manuscript_structure_JCCformat.pdf`, and the plan slides `JCC/Scoring_Manuscript_Plan_2026-06-29.pdf`.
> Last updated: 2026-07-01 (scope revision: Gramicidin/companion-paper split, Decision-8 retraction,
> τ renaming, Phase-6 re-tiering — see Changelog).

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
> ~~Phase 2: `src/classifier.py` (Algorithm 1)~~ ✓ **DONE 2026-07-01** — `classify_all_modes()`/
> `classify_to_rows()` implement Steps 2-4 (plain one-to-one `linear_sum_assignment` maximizing
> `Σ|score|`, two-gate purity, `vib_label` internal split + `s_AB`); `main.py` gained
> `build_scorer_and_final()` (factored out of `score_modes()`, shared by both) and
> `run_classify_pipeline()`. Wrote all 4 required outputs (`water`/`benzene` × `normal`/`EMIT`
> `_classified.csv`); formula-auditor PASS on Steps 2-4 (one DIVERGENT finding, fixed same session — see
> Changelog); score-validator PASS on all label/invariant targets. 13/13 tests green
> (`tests/test_classifier.py` new, 6 tests). **Next, in order:** Phase-0 Excel column verification
> (re-score one library molecule vs `data_score` `s[V_S]` column) → `src/projection.py`
> (EMIT→normal-mode projection, the other half of Phase 2) → Phase 3 (library ingest + calibration).
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
- **Algorithm 1** (PDF §B6.2): Step1 score → Step2 global **plain one-to-one** assignment via
  `linear_sum_assignment` over the `n_T+n_R` external slots vs. all modes, **MAXIMIZING `Σ|score|`**
  (scipy minimizes by default → negate the cost matrix or use `maximize=True`). **No block-constraint
  mechanism** (retracted 2026-07-01 per `JCC_manuscript_structure_scoped.md` Decision 8 — degenerate
  inertia-tensor axis choice is a labeling convention, not an assignment ambiguity, and normal-mode T/R
  references are constructed directly from geometry, never searched for) → Step3 **clean iff
  `|score|≥τ_TR` AND `s[V_S]≤τ_B`** (directionally aligned AND internally rigid), else
  `mixed_external` (flag + dominant slot + `s[V_S]`) → Step4 `s[V_S]≥τ_S` stretching / `≤τ_B`
  bending / else mixed. `n_T=3`, `n_R=2 if linear else 3`.
  **Two-gate purity (refined 2026-06-30):** `s[T]/s[R]` are direction-only (read ±1 for amplitude-varying
  impure modes); the second gate uses `s[V_S]` to catch **stretching-type** external impurity. Benzene
  EMIT: 34/35 (E₁ᵤ, Tx/Ty=1, `s[V_S]=0.667/0.577`) correctly flagged. **KNOWN BLIND SPOT:** EMIT 36 (A₂ᵤ,
  Tz=1, `s[V_S]=0`) is out-of-plane **bending**-mixed (amplitude variation ⊥ the in-plane bonds → bending,
  no stretching), yet its score signature {Tz=1, rest 0, V_S=0} is **identical to a pure z-translation** —
  {s[T],s[R],s[V_S]} cannot distinguish them (no bending observable). The flag detects stretching-type
  impurity only; bending-type is invisible (projection resolves it). Reuses τ_TR/τ_B; for exact
  normal-mode externals (genuinely rigid, `s[V_S]=0`) the gate never fires (completeness intact).
  **DECISION (LOCKED 2026-06-30): X — honest limitation, no new machinery.** Framing: the scores measure
  **directional/geometrical** character; 34/35/36 genuinely have large translational character (every
  displacement vector aligns perfectly with a Cartesian axis, `s[T]=1`) — only their MAGNITUDES differ.
  Equal magnitudes ⇒ pure translation; magnitude variation encodes the internal residual (stretching for
  34/35 → caught by `s[V_S]`; out-of-plane bending for 36 → invisible). The score is "partially true" — it
  correctly reports 36's dominant translational character; the magnitude-encoded bending residual is left
  to projection. 36 = the worked example of the score/projection boundary. (Y, a rigid-body residual, was
  rejected: it is projection onto the T/R subspace and reopens the substitutability objection.)
- Conventions: `ε_disp=1e-8`; `unit(0):=0`. Degeneracy-tolerance grouping of inertia axes/modes is
  **no longer part of the assignment mechanism** (Decision 8, retracted 2026-07-01) — `DEGEN_TOL`
  survives only as a general-purpose numeric constant (e.g. for the Phase-3 degenerate-mode-set sanity
  check), not as an axis-block/mode-block binding requirement. `τ_TR` is taken from **calibration**,
  NOT hardcoded (PDF example uses 0.95).

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
      (share `_build_inertia_tensor()` with `MIT`). **Note (2026-07-01):** `axis_blocks()` is retained
      as a diagnostic/introspection accessor only — Decision 8 retracts the requirement that the
      classifier consume axis-degeneracy blocks for assignment; `classifier.py` uses plain
      `linear_sum_assignment` over the full score table with no block-consumption step.

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
- [x] `src/classifier.py` (NEW) — **DONE 2026-07-01.** `n_T/n_R` via `is_linear()`/`external_slots()`;
      `Thresholds` dataclass (`τ_TR=0.95, τ_S=0.9, τ_B=0.2` — `τ_TR` to be replaced by calibration in
      Phase 3). No degenerate-axis-block data structure (block-handling retracted 2026-07-01,
      `JCC_manuscript_structure_scoped.md` Decision 8) — confirmed by formula-auditor: no `axis_blocks()`/
      `DEGEN_TOL` import anywhere in the module.
- [x] Step 2 global assignment — **DONE.** Plain one-to-one `linear_sum_assignment` over the `n_T+n_R`
      external slots vs. all modes, **maximizing `Σ|score|`** (`maximize=True`, with a sign-correct
      `-cost`/minimize fallback for older scipy). **No block-constraint mechanism** — retracted
      2026-07-01 (`JCC_manuscript_structure_scoped.md` Decision 8): translation invariance holds along
      any axis for any molecule; a degenerate inertia tensor's non-unique axes are a labeling convention
      (any deterministic eigensolver's fixed triple works, scoring proceeds normally); and normal-mode
      T/R references are built directly from geometry (Eckart/Sayvetz), so no assignment ambiguity exists
      there. Global assignment does genuine work only for EMIT mode sets (confirmed: water/benzene EMIT
      2, 6, 9 are non-trivially assigned then flagged `MIXED_EXTERNAL_WITH_VIBRATION`), and plain
      Hungarian handles those correctly with no special-casing.
- [x] **Step 3 two-gate purity (see Authoritative spec) — DONE.** Clean external iff `|s_slot|≥τ_TR` AND
      `s[V_S]≤τ_B`; else `MIXED_EXTERNAL_WITH_VIBRATION` (flag + dominant slot + `vib_label(s[V_S])`
      annotation). Reuses existing constants. Correctly flags 34/35 (stretching-mixed;
      `annotation="dominant_external=Tx; vibration=MIXED_STRETCH_BEND"` / `Ty`). **BLIND SPOT — EMIT 36**
      (A₂ᵤ out-of-plane bending, `s[V_S]=0`) has a score signature identical to a pure z-translation, so
      the scores cannot flag it (no bending observable) — reproduced exactly as `CLEAN_TRANSLATION`, per
      Decision X (this is the correct, intentional output, not a bug). **Manuscript refinement (B8.4):**
      "EMIT 34–36 flagged" → "34/35 flagged (stretching); 36 is the honest limit — its out-of-plane
      bending residual is invisible to the reference-free scores and resolved only by projection."
      **DECISION X locked** (honest limitation; no rigid-body residual). Discussion must note the scores
      capture genuine directional/geometrical translational character — only magnitudes differ — so the
      label is correct about dominant character.
- [x] Step 4 internal split (stretching/bending/mixed) + attach `s_AB` — **DONE.** `vib_label()` applies
      only to modes never claimed by any external slot in Step 2 (`classification is None` after Step 3);
      modes flagged `MIXED_EXTERNAL_WITH_VIBRATION` get their vibrational character only inside the
      annotation, never a second top-level label (formula-auditor confirmed no double-classification).
      Per-bond `s_AB` attached only for `STRETCHING`/`MIXED_STRETCH_BEND` (formula-auditor flag: this
      means stretching-flavored `MIXED_EXTERNAL_WITH_VIBRATION` modes, e.g. benzene EMIT 34/35, currently
      carry no per-bond detail at all — self-consistent with the literal spec wording as written, but
      worth a future author call if `fig:bondscores`/localization ever wants bond-level detail for the
      flagged-external cases too; **not changed this session**, flagged here for later).
- [x] Output `data/results/<mol>_classified.csv` (scores + label + annotations + `s_AB`) — **DONE for all
      4 required combos:** `water_normal_classified.csv`, `water_EMIT_classified.csv`,
      `benzene_normal_classified.csv`, `benzene_EMIT_classified.csv` (also verified for
      `co2_mp2_3-21g_normal` as an extra linear-molecule check — see Changelog bug/fix). Columns: Mode,
      Freq/Eigenvalue, Tx,Ty,Tz,Rx,Ry,Rz,V_Stretch, label, annotation, s_AB (semicolon-joined
      `i-j:value`, 1-based atom indices, blank when not attached).
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
- [ ] `src/calibrate.py` (NEW): derive `τ_S`,`τ_B` from labeled distributions; `τ_TR` sweep +
      plateau; **persist the full τ-sweep curve** (label-change fraction + accuracy per τ) for
      `fig:sensitivity`; write `data/results/thresholds.json`; back-fill `Thresholds`. Define "plateau"
      quantitatively (give a numeric criterion).
- [ ] **Label-level validation (now τ is frozen):** benzene EMIT 34–36 flagged `mixed_external`;
      confusion-matrix precision/recall vs lit labels (set a numeric acceptance floor); re-pin label
      goldens. **Sanity check only (not a mechanism):** confirm plain one-to-one Hungarian assignment
      happens to label degenerate mode sets consistently (e.g. benzene EMIT 1–9, sharing an eigenvalue)
      as an emergent property of the scores — no block-special-casing is implemented or needed
      (retracted 2026-07-01).
- [ ] Figure `fig:confusion` (clean-category confusion matrix + precision/recall vs lit labels).
- [ ] Figure `fig:bondscores` (bond score vs relative Δbond length, ideal vs non-ideal).
- [ ] Figure `fig:boxplots` (freq, Δ|b|, `s[V_S]`; stretch vs bend).
- [ ] Figure `fig:modemixing` (mode score vs averaged Δbond; ideal step vs non-ideal gradient).
- [ ] Figure `fig:sensitivity` (label-change fraction & accuracy vs τ; plateau).
- [ ] Parity check vs Excel `box plots` / `CM` sheets (and spot-check the other figures' source data).

## Phase 4 — DEFERRED to companion paper (Gramicidin A scalability; out of scope for this manuscript)
> **Scope decision (2026-07-01, `JCC_manuscript_structure_scoped.md` Decision 5):** Gramicidin A
> scalability/wall-clock/`fig:gramicidin`/per-bond `s_AB`-at-scale, transition-state/bond-breaking
> characterization, and isotopic-substitution mode comparison are ALL deferred to a future companion
> paper built around `s_AB` as a standalone analytical tool. **Why:** this paper's coverage claim
> ("one framework classifies all 3N modes in one pass") does not depend on scale — water, the hydride
> library, and benzene (normal + EMIT) fully support it. Nothing is lost from the core argument; only
> the *empirical, at-scale* demonstration is deferred, and the manuscript's Conclusion states this
> explicitly so the absence of a large system reads as a decision, not a gap.
>
> **This phase is recorded here, not deleted, so the absence is legible as intentional.** No gramicidin
> work is active for this manuscript: not the NumPy vectorization pass, not the `1grm_MM_UFF.log`
> pipeline run, not `s_AB` localization output, not `fig:gramicidin`, not the wall-clock comparison
> table. `data/logs/1grm_MM_UFF.log` and `data/gjf/1grm.com` remain tracked in the repo (companion-paper
> input) but are not consumed by anything in the current build order. Revisit this phase only when the
> companion paper begins.

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

## Phase 6 — Strengthen for review (RE-TIERED 2026-07-01; see Changelog)
> From expert-reviewer-jcc, re-triaged after Decision 5 (Phase 4/Gramicidin deferred to the companion
> paper) removed the paper's only scale/robustness demonstration — raising the weight the remaining
> validation has to carry. Split into a **recommended-before-submission** tier and an **optional** tier.

### Recommended before submission
- [ ] **Flag precision/recall over ALL 36 benzene EMIT modes** (and library externals) against the C2
      projection reference — systematic flag-fidelity, vs the current anecdotal EMIT 2/9/34–36 (A4).
      **Promoted:** with no large-system demonstration left in this paper, benzene EMIT is now the
      paper's only stress test of the flagging mechanism — a 3-mode anecdotal spot-check can no longer
      carry that load; the full 36-mode confusion set is needed.
- [ ] **Validate the mixed-SB bucket** by irrep-degeneracy + CoM arguments; report the fraction of
      lit-labeled modes landing in "mixed" (defends against the bending = low-stretch circularity
      concern). **Promoted:** this is the direct evidentiary backbone for the manuscript's residual-risk
      discipline that mixed-SB/flagged-external buckets be "validated by characterization and
      consistency, never by an accuracy claim" — without this item that discipline is asserted but not
      discharged.
- [ ] **Out-of-sample / leave-one-molecule-out** evaluation of the τ-calibrated classifier — answers the
      "reference-free vs trained-τ" objection that the sensitivity plateau alone does not. **Promoted:**
      with the computational-cost argument now the primary defense against "why not projection/PED," a
      referee who accepts that argument pivots next to "are your thresholds actually reference-free, or
      secretly fit to the test set" — this is the direct answer.

### Optional (nice-to-have, not fatal if deferred)
- [ ] **N-per-cell + confidence intervals** on the confusion matrix (thin bend counts invite a
      significance objection).
- [ ] **CoM-conservation evidence** (correlate central-atom amplitude / neighbor mass vs `s[V_S]`
      degradation, TeH₂ vs Br₂O) backing the B8.3 "explained, not noisy" difficulty gradient.

---

## Verification
- [ ] **Score-level (Phase 1):** water matches `tab:water` to 3 dp; `Σ s_AB == s[V_S]` (1e-6); externals
      reach ±1.000 (`|x|≥0.9995`); ranges hold; CO₂ exercises `n_R=2`; benzene EMIT 34–36 `s[V_S]` and
      EMIT 2/9 `|s[Ry]|` to 3 dp. All frozen as goldens.
- [ ] **Label-level (Phase 3, τ frozen):** benzene EMIT 34–36 flagged `mixed_external`; degenerate mode
      sets receive consistent labels as an emergent property of plain Hungarian assignment (no block
      mechanism); confusion-matrix precision/recall ≥ floor; calibrated τ on the plateau.
- [ ] Figures match Excel `box plots`/`CM` sheets (+ spot-checks for the rest).
- [ ] ~~Gramicidin run completes with timing + `s_AB`; numbers inserted where `.tex` has `TODO-DATA`.~~
      **REMOVED 2026-07-01** — Gramicidin deferred to the companion paper (see Phase 4 above); this
      manuscript carries no gramicidin verification target.
- [ ] `py reproduce.py` regenerates all CSVs + figures with no manual steps; `pytest` green.

## Open items to confirm during execution
- [ ] (Phase 0) eq:emitproj mass-weighting convention — pinned pre-Phase-2.
- [ ] (Phase 0) `data_score` column == `s[V_S]` — verified by in-engine re-score.
- [x] (Phase 0) library optimized geometries — RESOLVED: author supplies `.log`/`.gjf` on request at the Phase-5 SI step.
- [x] Step-2 objective is **maximize** `Σ|score|` (not scipy's default minimize) — implemented in
      `src/classifier.py` via `linear_sum_assignment(cost, maximize=True)`; formula-auditor confirmed
      sign correctness including the `-cost`/minimize fallback path.
- [ ] V-score uses the **initial** bond direction `b̂^{i}`.
- [ ] Degeneracy tolerance numeric value (axes by `λ_i`, modes by freq) — fix and document.
- [ ] ~~Gramicidin bond-connectivity criterion (`TODO-DATA`): use `1grm.com` topology.~~ **DEFERRED
      2026-07-01** — moot for this manuscript (Phase 4 deferred to companion paper).
- [x] Flag-criterion: two-gate purity (`|s|≥τ_TR` AND `s[V_S]≤τ_B`) flags stretching-type impurity (34/35). EMIT 36 out-of-plane bending blind spot RESOLVED as **Decision X** (honest limitation; scores report genuine directional/geometrical translational character, magnitudes differ; bending residual → projection). No rigid-body residual added.
- [ ] Verify in Phase 2 the per-mode projected translational % for 34/35/36 (36 is NOT ~100% — it has out-of-plane bending despite `s[V_S]=0`; do not assume the "77%" applies uniformly).

---

## Changelog
- **2026-07-01 — Phase 2 classifier landed (`src/classifier.py`, Algorithm 1 Steps 2-4) + a linear-
  molecule pool bug found and fixed.** Added `src/classifier.py`: `Thresholds` dataclass
  (`τ_TR=0.95, τ_S=0.9, τ_B=0.2`); `is_linear()`/`external_slots()` (`n_T=3`, `n_R=2 if linear else 3`);
  `classify_all_modes(scorer, final, thresholds)` — Step 1 rescoring (delegates to
  `ModeScorer.calculate_scores`/`score_bonds`, no new formulas), Step 2 plain one-to-one
  `linear_sum_assignment(cost, maximize=True)` over the `n_T+n_R` slots vs. every mode in the pool (no
  block-constraint mechanism, per the Decision-8 retraction), Step 3 two-gate purity
  (`|s_slot|≥τ_TR` AND `s[V_S]≤τ_B` → `CLEAN_TRANSLATION`/`CLEAN_ROTATION`, else
  `MIXED_EXTERNAL_WITH_VIBRATION` + `dominant_external=<slot>; vibration=<vib_label>` annotation), Step 4
  `vib_label()` for modes never claimed by a slot (`STRETCHING`/`BENDING`/`MIXED_STRETCH_BEND`, per-bond
  `s_AB` attached only for the stretching-flavored two). `classify_to_rows()` flattens to the CSV shape.
  Factored `main.build_scorer_and_final()` out of the existing `score_modes()` (byte-identical output,
  confirmed by score-validator) so both Step-1 scoring and the new `main.run_classify_pipeline()` build
  the identical candidate pool from the same raw parse. Generated all 4 required outputs:
  `data/results/{water,benzene}_{normal,EMIT}_classified.csv`.
  **Validation:** water's 6 ideal T/R references → clean (Tx/Ty/Tz `CLEAN_TRANSLATION`, Rx/Ry/Rz
  `CLEAN_ROTATION`); Vib 1 (bend, `V=0.061`) → `BENDING`; Vib 2/3 (stretches, `V=1.000`/`0.998`) →
  `STRETCHING` with `s_AB` summing to `V` over the 2 O-H bonds. Benzene EMIT 34 (`|Tx|=1, V=0.667`) and
  35 (`|Ty|=1, V=0.577`) → `MIXED_EXTERNAL_WITH_VIBRATION` (pass gate 1, fail gate 2); EMIT 36
  (`|Tz|=1, V=0`) → `CLEAN_TRANSLATION`, the documented Decision-X blind spot (out-of-plane bending is
  invisible to the two-gate test since `s[V_S]=0` for it too) — reproduced exactly as specified, not
  "fixed." formula-auditor: PASS on Steps 2-4 mechanics (Hungarian correctness/direction, no block
  machinery, two-gate exact form, no double-classification), with one **DIVERGENT** finding (see below).
  score-validator: full PASS (all label/invariant targets, `score_modes()` byte-identical, 12/12 tests
  green at review time).
  **Bug found + fixed (formula-auditor):** for a linear molecule, `ModeScorer.construct_R()` always
  builds 3 ideal rotation references, but `MIT()` places the smallest-moment (molecular) axis on the new
  X axis, so the ideal "Rx" reference is an all-zero vector for a linear molecule — not a genuine external
  mode (`n_R=2` should exclude it). Left unfiltered, it entered the candidate pool, was never claimed by
  any Step-2 slot (`external_slots()` correctly omits Rx when linear), and fell through Step 4 to be
  mislabeled `BENDING` (`V=0≤τ_B`) — a spurious `3N+1`-mode pool with a meaningless row (confirmed on
  `co2_mp2_3-21g`, previously unexercised since no test ran a linear molecule through the classifier).
  Fixed in `main.build_scorer_and_final()`: when `is_linear(scorer)`, drop the `"Rx"` entry from
  `scorer.construct_R()`'s output before assembling `final`, so the on-axis placeholder never enters
  either the Step-1 scores or the classifier's pool. Verified: `co2_mp2_3-21g` now yields exactly 9
  rows (`3N`) with `Ry`/`Rz` clean and 2 real bends + 2 real stretches, no spurious `Rx` row; water/
  benzene (neither linear) unaffected (row counts unchanged, byte-identical scores). Added
  `test_co2_linear_no_spurious_onaxis_mode` regression test. Not yet resolved (flagged in the Phase-2
  checklist for a future author call, not a bug per the literal current spec): stretching-flavored
  `MIXED_EXTERNAL_WITH_VIBRATION` modes (e.g. benzene EMIT 34/35) currently get no per-bond `s_AB` at all,
  since bond attachment is gated on the top-level classification being `STRETCHING`/`MIXED_STRETCH_BEND`
  only — self-consistent with the spec as written, but a candidate gap if bond-level detail is later
  wanted for the flagged-external cases too (e.g. for `fig:bondscores`/localization).
  Added `tests/test_classifier.py` (6 tests, all green): water externals clean, water Vib1/2/3 split,
  benzene EMIT 34/35 flagged with correct annotation, EMIT 36 clean (blind spot), the CO2 pool-size
  regression, and `classify_to_rows()` column-shape check. Full suite: 13/13 green
  (`py -m pytest tests/`).
- **2026-07-01 — Manuscript scope revision (Gramicidin deferred; Decision-8 retraction; τ renamed;
  Phase 6 re-tiered).** `JCC/JCC_manuscript_structure_scoped.md` (01 July 2026) supersedes prior
  manuscript planning and forces five changes here:
  (1) **Phase 4 (Gramicidin A) removed as an active phase, deferred to a future companion paper** built
  around `s_AB` as a standalone tool — alongside transition-state/bond-breaking characterization and
  isotopic-substitution mode comparison (Decision 5). This paper's coverage claim ("one framework
  classifies all 3N modes in one pass") does not depend on scale — water + the hydride library +
  benzene (normal + EMIT) fully support it — so nothing is lost from the core argument; only the
  empirical at-scale demonstration is deferred. Recorded as an explicit "DEFERRED" phase, not deleted,
  so the absence reads as a decision, not an oversight; `1grm_MM_UFF.log`/`1grm.com` stay tracked
  (companion-paper input) but are out of scope for this manuscript.
  (2) **"Degenerate-block assignment" retracted (Decision 8).** The earlier requirement that Step 2's
  global assignment special-case degenerate inertia-tensor axis-blocks (for symmetric/spherical tops)
  is INCORRECT and removed: translation never needs principal axes at all; a degenerate inertia
  tensor's non-unique axes are a labeling convention fixed once by the eigensolver, not an assignment
  ambiguity; and normal-mode T/R references are constructed directly from geometry (Eckart/Sayvetz),
  never discovered by search, so no ambiguity exists there either. Step 2 is now **plain one-to-one**
  `linear_sum_assignment` over all `n_T+n_R` external slots vs. all modes, maximizing `Σ|score|` — no
  block-constraint data structure or mechanism, anywhere. It does genuine work only for unlabeled mode
  sets (EMIT); for normal modes it remains confirmatory. `axis_blocks()` (Phase 0) is retained purely as
  a diagnostic accessor, not consumed by the classifier. The Phase-3 "degenerate-block consistency
  check" is downgraded from a required mechanism to a **sanity check**: confirming degenerate mode sets
  (e.g. benzene EMIT 1–9) happen to receive consistent labels as an *emergent property* of plain
  Hungarian assignment, not evidence of any block machinery (there is none).
  (3) **Threshold renaming:** `τ_pure→τ_TR`, `τ_stretch→τ_S`, `τ_bend→τ_B` throughout (spec, Phase 2/3
  checklist items, Verification, Open items) — cosmetic, no numeric/logic change.
  (4) **Phase 6 re-tiered** (expert-reviewer-jcc advisory): with Gramicidin gone as the paper's only
  scale/robustness demonstration, three items are promoted from optional to **recommended before
  submission** — full 36-mode benzene-EMIT flag precision/recall (no longer anecdotal EMIT 2/9/34–36
  only), mixed-SB bucket validation (irrep + CoM), and out-of-sample/leave-one-molecule-out threshold
  evaluation. N-per-cell confidence intervals and CoM-conservation evidence remain optional.
  (5) **Computational-cost section (`JCC_manuscript_structure_scoped.md` §8) is now the primary
  defense** against "why not projection/PED" — with no large system in the paper, the three-stage
  (construct/score/classify) op-count comparison carries more rebuttal weight than it would have
  alongside a Gramicidin scale demo. Phase 5/6 work should treat that section's op-count derivation
  (extending Appendix B to `s[T]`/`s[R]`) as a verification target, not just a manuscript exhibit.
  See `JCC_manuscript_structure_scoped.md` Decisions 5/6/7/8 for full rationale. **Not yet done this
  session (flagged for a future writing session):** the `.tex` itself (`JCC_man_scoring/
  JCC_temp_LaTeXtemplate.tex`) still contains the old gramicidin subsection, `τ_pure/τ_stretch/τ_bend`
  notation, and degenerate-block language — a `tex-data-sync` gap-list and a `lead-author` punch list
  for that edit pass exist from this session's review but were deliberately not applied (the `.tex`
  lives outside this git repo, per the note below, and deserves its own session with `check-tex`/
  `check-figures` run afterward). Also flagged: Fig 4 needs a benzene ~400 cm⁻¹ C–C-stretch worked
  example (currently only exists for ethane at 976 cm⁻¹ — a real content gap, not a relabeling), and
  `fig:modemixing` likely needs additional panels for the irrep-degeneracy argument (currently prose-only,
  using ethane not benzene).
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
