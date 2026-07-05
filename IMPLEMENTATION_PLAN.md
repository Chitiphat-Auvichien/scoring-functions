# JCC Scoring/Classification Program — Implementation Plan & Progress Checklist

> Living checklist. Tick `[x]` as parts are completed; pick up unchecked items in any later session.
> Companion to the JCC manuscript `JCC/JCC_man_scoring/JCC_temp_LaTeXtemplate.tex`, the content plan
> `JCC/JCC_manuscript_structure_JCCformat.pdf` (positioning spine — still current), and the plan slides.
> **Manuscript-structure slides: `JCC/Scoring_Manuscript_Plan_2026-07-05.pdf` is now the CURRENT
> authoritative Results & Discussion structure/content doc** (supersedes `..._2026-07-02.pdf`, which
> superseded `..._2026-07-01.pdf`, which itself superseded `..._2026-06-29.pdf` — all four are kept on
> disk for history; always use the newest-dated one unless told otherwise). Any agent working on the
> manuscript's Results & Discussion section (`lead-author` especially, but also
> `figure-builder`/`lead-engineer` when their output feeds a specific section) should be pointed at this
> file, not an older one.
> Last updated: 2026-07-05 (single-centre-only hydride-library scope fix — non-ideal population 58→51
> molecules / 422→277 modes, rigorous tier 237→231; see RESUME HERE below for the full writeup.)

> ## ▶ RESUME HERE (session pointer — keep current; update + commit after each increment)
> **2026-07-05 (author decision + lead-engineer, single-centre-only hydride-library scope fix):** The
> manuscript's "hydride library" validation (tau_S/tau_B derivation + the ideal/non-ideal confusion-matrix
> statistics in the "Stretching/bending classification" section) is explicitly scoped to single-centre
> AB_n topologies only ("one distinguishable atom anchors all bonds symmetrically" — JCC .tex, that
> subsection; `tab:ideal` lists exactly 11 single-centre AB_n molecules). The non-ideal side of the
> library was violating this same scope: 7 of the 58 non-ideal molecules in `library_scores.csv` are NOT
> single-centre — `C2H2`/`C2H4`/`C2H6` (two-centre), `H2O2` (two-centre), `C6H6` (six-centre, benzene),
> `iso-C4H10`/`n-C4H10` (four-centre) — contributing 145 of 422 non-ideal internal modes (~34%) to
> statistics that were supposed to be single-centre-only. Removing them brings the non-ideal population to
> the **51 molecules / 277 modes** an independent hand-count of `tab:nonideal`'s own AB_n grid already
> implied. `H2O` (single-centre) stays, as a familiar illustrative molecule.
> **Mechanism chosen:** NOT an extension of `src/excel_ingest.py`'s `EXCLUDED_MOLECULES` (that constant is
> `{"Gly5"}`, an ingest-time exclusion for a molecule with zero rows in the library at all — a different
> kind of exclusion). Instead, a new analysis-time filter lives in `src/calibrate.py`:
> `SINGLE_CENTRE_ONLY_EXCLUDE = {"C2H2","C2H4","C2H6","H2O2","C6H6","iso-C4H10","n-C4H10"}` +
> `filter_single_centre_library(lib_df)` (drops rows by `molecule.isin(...)`, any `kind`). Deliberately
> does NOT touch `library_scores.csv` itself or `excel_ingest.py` — every one of these 7 molecules' rows
> stays in the CSV (benzene needs its rows there for its own dedicated confusion matrix). Applied
> unconditionally as the first line of `confusion_matrix_stats()` (so every caller — `plot_confusion_
> matrix`'s rigorous/non-ideal tiers, and the direct-call regression tests in `tests/test_calibrate.py` —
> gets the correct scope automatically, without each needing its own copy of the filter) and directly by
> `src/figures.py`'s `plot_bond_scores`/`plot_boxplots`/`plot_mode_mixing` (which build their populations
> straight from `library_scores.csv`, not through `confusion_matrix_stats`). `derive_stretch_bend_
> thresholds()` was confirmed (not assumed) unaffected: none of the 7 excluded molecules has an
> `ideal=='yes'` row, so `tau_S=0.90368`/`tau_B=0.17327`/`tau_TR=0.95` are bit-identical to before.
> **Rigorous tier also shrank, computed not assumed:** benzene (`C6H6`) is the only one of the 7 that is
> geometry-backed (on-disk `.log`+`.gjf`), so besides its 30 non-ideal internal rows, it also contributed 6
> geometry-backed external T/R rows that were being pooled into the "rigorous" tier (`kind=='external'`
> regardless of `ideal` tag) — n=237→**231** (146→140 external + 91 ideal-internal, unchanged). The other 6
> excluded molecules are Excel-only (no on-disk geometry) and contribute only non-ideal internal rows.
> **Recomputed confusion-matrix numbers (pooled, full filtered library):** stretch recall 0.70815→0.66667
> (mixed_fraction 0.29185→0.33333), bend recall 0.97482→0.98378 (mixed_fraction 0.02518→0.01622);
> precision stays 1.0 for all 4 categories, `floor_met` stays False. Ideal-tier recall/n_ref (stretch
> 1.0/41, bend 1.0/50) are exactly unchanged, as required by construction. Non-ideal-tier: bend
> `n_ref_nonideal` 228→135 (`recall_nonideal` 0.96930→0.97778), stretch `n_ref_nonideal` 192→142
> (`recall_nonideal` 0.64583→0.57042).
> **Benzene's own dedicated confusion matrix (`benzene_internal_confusion_matrix.csv`/`_summary.csv`,
> `fig:benzeneconfusion`, commits `69d549e`/`d298a21`) is completely untouched** — it reads `library_
> scores.csv`'s C6H6 rows directly via `src/benzene_validation.py`, which never calls `confusion_matrix_
> stats()` or the new filter; git diff on those two CSVs is empty.
> **All 6 affected figures regenerated** (`py main.py --figures`): `fig:confusion` (n=231/277 in the two
> panel titles), `fig:bondscores` (62 molecules, was 69), `fig:boxplots`/`fig:modemixing` (368 total
> modes = 91 ideal + 277 non-ideal, was 91+422=513). `fig:benzene`/`fig:benzene_normal`/
> `fig:benzeneconfusion`/`fig:sensitivity` were also regenerated as part of the same `--figures` run but
> are pixel-identical (only their PDF metadata/timestamp changed) since none of them read the filtered
> population.
> **Tests:** `tests/test_calibrate.py` — 2 existing tests re-pinned to the new numbers above (docstrings
> explain why, same convention as prior sessions) + 3 new tests (`test_single_centre_only_exclude_matches_
> scope_decision`, `test_filter_single_centre_library_drops_exactly_the_excluded_molecules` — pins the 151
> dropped rows (145 internal + benzene's 6 external) and the 51/277 non-ideal result,
> `test_confusion_matrix_stats_applies_single_centre_filter_internally` — confirms the filter is applied
> regardless of caller pre-filtering). `py -m pytest tests/` **80/80 green** (up from 77). **score-
> validator PASS**: every number above independently re-derived from `library_scores.csv` directly
> (not just re-read from this writeup), benzene's dedicated matrix confirmed byte-identical to HEAD, and
> `git diff --stat -- src/scoring.py src/classifier.py` confirmed empty (core engine untouched).
> **Not done / out of scope:** `src/csv_label_ingest.py` and its 6-molecule xlsx-fallback list happen to
> include the same `C2H2`/`C2H4`/`C2H6`/`H2O2`/`iso-C4H10`/`n-C4H10` set (a coincidence of which molecules
> the new label CSVs don't cover yet, unrelated to this scope decision) — left untouched, since that list
> is about label/citation SOURCE, not statistical scope.
> **2026-07-05 (author decision + lead-engineer, new ref_label/ideal/citation-key source across ALL
> molecules -- literature relabeling of benzene modes 21/22/19/23/24):** The author replaced the
> label/citation content that used to live only inside `data/vibrational-scoring-functions.xlsx`'s
> `data_score`/`characterised modes` sheets with three tracked CSVs: `data/data_score.csv` (63
> molecules, re-export of the `type`/`ideal` columns, **plus a literal new `"SB"` (mixed) literature
> class** for benzene modes 21/22 -- the first genuine, literature-sourced 3rd reference class anywhere
> in this pipeline, not just the classifier's own predicted mixed bucket), `data/characterised_modes.csv`
> (63 molecules, richer per-mode detail incl. a `ref` citation-key column, only ~11/63 molecules
> currently back-filled), `data/ref-label_citation.csv` (citation key -> doi/author/journal/year, ~9
> distinct keys today). Hand-verified against the CURRENT manuscript text
> (`JCC/JCC_man_scoring/JCC_temp_LaTeXtemplate.tex` lines 794-917): of benzene's 30 internal modes, mode
> 19 (1319.27 cm-1) and the degenerate pair 23/24 (1598.90 cm-1) flip literature `bend`->`stretch`; the
> degenerate pair 21/22 (1532.85 cm-1) flips literature `bend`->`SB`. All other benzene modes and all 63
> non-benzene molecules are unchanged in VALUE (only re-exported to CSV).
> **New module `src/csv_label_ingest.py`** builds `{(molecule, mode): {ref_label, ideal, ref_key}}` from
> the 3 new CSVs, extending the accepted literal `type`/ref_label values from `("bend","stretch")` to
> `("bend","stretch","SB")` (anything else still maps to `None`, unchanged fail-quiet-not-fail-wrong
> policy). **Coverage gap, resolved via fallback, not silently dropped:** 6 molecules the old
> xlsx-driven pipeline scored (`C2H2`,`C2H4`,`C2H6`,`H2O2`,`iso-C4H10`,`n-C4H10`) are entirely absent
> from the new CSVs (not yet migrated by the author) -- `build_label_lookup(tables, fallback_ds=...)`
> falls back to the xlsx `data_score` sheet's own `type`/`ideal` ONLY for molecules the new CSVs don't
> cover at all; every molecule the new CSVs DO cover always wins (full switch, not a merge/override of
> individual fields). **`src/excel_ingest.py` wiring:** `ingest_internal_rows()` (source="excel") and
> `attach_excel_labels()` (source="gaussian") now source `ref_label`/`ideal`/`ref_key` from this lookup
> instead of the xlsx `data_score` sheet directly; `SCHEMA_COLUMNS` gained one new APPENDED column,
> `ref_key` (citation key, e.g. `"Shi1972"`; blank where the author hasn't back-filled a citation yet)
> -- purely additive, existing consumers (`calibrate.py`/`figures.py`) read columns by name so this does
> not disturb them. Freq/V_Stretch/delta_b_mean/s_AB/rel_db/has_geometry/predicted_label/Tx..Rz are
> COMPLETELY untouched by this change and still come from whichever source (Excel or Gaussian) the
> existing `source=` parameter already selected -- only label/citation content moved.
> **Regression-verified:** regenerated `library_scores.csv` (both `source="excel"` -- 659 rows/69
> molecules -- and `source="gaussian"` -- 303 rows/25 molecules, both matching their documented prior
> baselines exactly) and diffed column-by-column against the prior committed golden: **every molecule
> EXCEPT C6H6 is byte-for-byte identical in `ref_label`/`ideal`/every other column**; C6H6's `ref_label`
> differs ONLY at mode_index 19/21/22/23/24, exactly as hand-derived above. **Side-effect noted, not a
> regression:** C6H6's `freq` column also shifted by up to 0.015 cm-1 across all 30 modes (e.g. the
> 13/14 and 23/24 near-degenerate pairs are now EXACTLY degenerate to the printed precision, 1056.3901/
> 1056.3901 and 1598.8987/1598.8987, vs. tiny spurious splittings before) -- this traces to
> `data/vibrational-scoring-functions.xlsx` itself having been resaved by the author today (its on-disk
> mtime is newer than the last `library_scores.csv` commit), independent of this session's code change;
> `freq` sourcing was not touched by this session at all. No other molecule's freq shows any drift.
> **`src/benzene_validation.py` extended for the new 3-class ground truth (NOT touching
> `src/scoring.py`/`src/classifier.py`/any threshold):** `benzene_normal_reference_detail()`'s `correct`
> column now uses `_expected_pred_bucket(ref_label)` (`{"SB": "mixed"}.get(ref_label, ref_label)`) instead
> of comparing `predicted_bucket` to `ref_label` directly, since there is no predicted bucket literally
> named `"SB"` (`classification_bucket()` collapses the engine's own `MIXED_STRETCH_BEND` label to
> `"mixed"`) -- a literature `"SB"` mode counts as correctly classified iff the engine calls it `"mixed"`.
> `benzene_normal_reference_summary()`'s loop extended to 5 categories (added `"SB"`). **New function
> `benzene_internal_confusion_matrix()`**: genuine 3x3 (ref bend/stretch/SB x predicted bend/stretch/mixed)
> confusion table + per-category recall for benzene's 30 internal modes, built via `pd.crosstab` (mirrors
> `calibrate.py::confusion_matrix_stats`'s pattern; not a copy, since that function's external-row/
> ideal-tier assumptions don't apply to this benzene-only view). **Exact numbers, hand-derived by the
> author and independently re-verified by score-validator:**
>   - ref bend (n=18): 16 correct, 2 -> mixed (modes 13,14) — recall 0.8889.
>   - ref stretch (n=10): 7 correct, 3 -> mixed (modes 19,23,24) — recall 0.700.
>   - ref SB (n=2, modes 21,22): 0 correct, BOTH -> predicted clean `bend` (V_Stretch 0.09072/0.08119,
>     both < tau_B=0.17327) — recall 0.000, the tau_B two-gate purity test's "bending blind spot" hitting
>     a genuine literature-mixed mode for the first time. Zero pure bend<->stretch crossings (both
>     categories' `n_crossed_opposite==0`) — that part of the "zero crossings" claim still holds exactly.
> **New function `benzene_sb_vs_stretch_bond_diagnostic()`** (Task B', a sibling to the existing
> `benzene_mixed_bond_diagnostic`, NOT an extension of it, since 21/22 are specifically NOT
> predicted-mixed — finding them via `predicted_bucket=="mixed"` would find nothing): per-bond C-C vs.
> C-H breakdown for 21/22 (`case="blind_spot_bend"`) contrasted directly against the near-degenerate pair
> 23/24 (`case="overflagged_mixed"`, literature stretch but predicted mixed) — the manuscript's two
> contrasting miss mechanisms side by side. The contrast pair (23/24, not lone mode 19 which is already
> the dedicated SB worked example elsewhere per Task E) is derived by keeping only literature-stretch/
> predicted-mixed candidates that share an exact frequency with another candidate (a near-degeneracy
> check), not hardcoded. Shares a refactored `_bond_row_stats(r)` helper with `benzene_mixed_bond_
> diagnostic` (behavior-preserving refactor, confirmed by formula-auditor and by the unchanged
> `benzene_mixed_bond_diagnostic` test suite). Finding (supporting evidence, not the primary test): 21/22
> sit at only ~83-85% C-C fraction of V_Stretch (meaningful C-H contribution, 0.0136/0.0138) vs. >95% for
> 23/24/13/14/19 — a plausible mechanistic reason 21/22 read as more bend-like to the framework than the
> literature's own "mixed" call.
> **Side effect on the EXISTING (non-benzene-specific) hydride-library confusion machinery
> (`src/calibrate.py::confusion_matrix_stats`, `fig:confusion`'s numbers) -- untouched code, but its
> INPUT population changed since benzene is part of the pooled ~69-molecule library:** bend's PRECISION
> is no longer exactly 1.0 -- now 271/273 = 0.99267 (was 1.0). Modes 21/22 still predict "bend"
> (predicted_label/thresholds untouched) but their true reference is now the 3rd class "SB", which this
> pooled 4-category (stretch/bend/translation/rotation only) confusion table does not recognize as a ref
> category at all -- so those 2 modes count toward bend's `n_pred` denominator but not its `tp` numerator,
> a real and honest precision cost of introducing genuine 3-class ground truth for 2 modes, not a
> regression. Pooled bend recall moved 0.96466->0.97482 and stretch recall 0.71739->0.70815 (removing 5
> modes from bend's ref population while redistributing 3 to stretch changes both denominators); floor_met
> still False (stretch recall still well under 0.95). `tau_S`/`tau_B`/`tau_TR` themselves are UNCHANGED
> (derived only from `ideal=='yes'` rows; benzene's `ideal` is `'no'` for every row, so this is purely a
> non-ideal-tier population effect) -- `recall_ideal["stretch"]==1.0`/`recall_ideal["bend"]==1.0` still
> hold exactly, as required by construction. **Not yet done, flagged for `lead-author`/`figure-builder`:**
> `fig:confusion`'s already-typeset caption claims "precision 1.0 all 4 categories" -- that specific claim
> is no longer literally true once this CSV switch is applied and the figure/manuscript prose regenerated
> from the new `library_scores.csv`; the figure itself (`src/figures.py::plot_confusion_matrix`) was
> deliberately NOT touched this session (out of scope; it needs a wording decision, not a code fix).
> **Tests:** `tests/test_csv_label_ingest.py` (new, 8 tests) unit-tests `build_label_lookup`'s SB
> acceptance + fallback-vs-supersede semantics in isolation (no xlsx/CSV I/O needed for most cases).
> `tests/test_excel_ingest.py`'s two `attach_excel_labels()` unit tests updated to build a
> fallback-only `label_lookup` (via a new `_fallback_only_label_lookup()` test helper) since that
> function's signature gained a required `label_lookup` parameter. `tests/test_calibrate.py`'s two
> confusion-matrix tests re-pinned to the new precision/recall numbers above (docstrings explain the
> "why" so a future reader isn't confused about which population a number belongs to).
> `tests/test_benzene_validation.py`: the two Task-A tests updated to the new 5-category numbers; two new
> tests added for `benzene_internal_confusion_matrix` (pins the exact 3x3 table above) and
> `benzene_sb_vs_stretch_bond_diagnostic` (finds exactly {21,22,23,24}, not 19). **`py -m pytest tests/`
> 77/77 green** (up from 66; net +11 new tests, 0 removed, ~118s). **formula-auditor PASS** (data-sourcing
> change only; `_expected_pred_bucket` mapping confirmed sound against `classification_bucket()`;
> `_bond_row_stats` refactor confirmed behavior-preserving; `src/scoring.py`/`src/classifier.py` confirmed
> untouched). **score-validator PASS** (ref_label diff scoped exactly as claimed; water/Σs_AB/external-1.000/
> benzene-EMIT invariants unaffected and reproduced exactly; new 3-class confusion numbers independently
> re-verified). Commits: see `git log` for the exact hashes (csv_label_ingest.py + excel_ingest.py wiring;
> benzene_validation.py 3-class extension; test updates; regenerated data/results/*.csv).
> **2026-07-04 (author decision + lead-engineer, dual-source `excel_ingest.py` — REFINEMENT of, NOT a
> reversal of, the 2026-07-03 "recompute from Gaussian" directive below):** The manuscript's
> ALREADY-TYPESET figures (fig:confusion, fig:bondscores, fig:boxplots, fig:modemixing) were built from
> the OLD Excel-scored ~70-molecule `library_scores.csv`, i.e. the commit-`149fc62` population/values,
> BEFORE the 2026-07-03 disk-driven rearchitecture shrank the population to ~25 molecules (only those
> with on-disk Gaussian `.log`+`.gjf` pairs). The author is having Gaussian trouble right now and will
> keep adding real files for the rest of the hydride library incrementally, but wants the manuscript's
> current figures reproducible in the meantime. **Decision: `src/excel_ingest.py`'s
> `build_library_scores()`/`run_ingest_pipeline()` (and `src/calibrate.py`'s
> `run_calibration_pipeline()`) gained an explicit `source="excel"|"gaussian"` parameter** (CLI:
> `main.py --library`/`--calibrate` gained a matching `--source {excel,gaussian}` flag, default
> `"excel"`). Both sources produce the exact same locked `SCHEMA_COLUMNS` output, so
> `src/calibrate.py`/`src/figures.py` never need to know or care which source produced a given
> `library_scores.csv`.
> - `source="excel"` (**new default**, commit-`149fc62` logic restored verbatim, not rewritten from
>   scratch): precomputed scores read directly from `data/vibrational-scoring-functions.xlsx`
>   (`ingest_internal_rows` from `data_score`+`data_mode&bond`) for the full ~70-molecule library;
>   `attach_geometry_classification` overlays real-engine `predicted_label`/`Tx..Rz`/`has_geometry` onto
>   internal rows AND always appends the ideal external T/R rows for the subset (~25 today) that also has
>   on-disk `.log`/`.gjf` geometry, gated by the same whole-molecule frequency-agreement check as before.
>   `V_Stretch`/`s_AB`/`rel_db`/`delta_b_mean`/`ref_label`/`ideal` on internal rows are ALWAYS Excel's own
>   values under this source, even for geometry-backed molecules — only the overlay columns and the
>   always-geometry-only external rows come from the engine.
> - `source="gaussian"` (2026-07-03 disk-driven rearchitecture, kept alive and fully working, unchanged):
>   every score column recomputed from scratch via the real engine for every molecule with a `.log`+`.gjf`
>   pair on disk; Excel supplies only `ref_label`/`ideal`. This remains the intended EVENTUAL default once
>   the full library has on-disk geometry — **no code changes needed when that happens**, same as before.
> **`data/results/library_scores.csv` regenerated with the new default (`source="excel"`) via
> `py main.py --library`: 659 rows / 69 molecules** (`df["molecule"].nunique() == 69`) — confirmed to
> match the documented pre-2026-07-03 baseline exactly (2026-07-02's DONE entry below cites "the 25
> geometry-backed molecules... get real Algorithm-1 predicted labels" and the 2026-07-03 entry itself
> separately cites "the prior Excel-driven baseline's ~69 molecules / 659 rows" when describing the
> shrinkage this session reverses). The same 4 known frequency-mismatch molecules (H2O/OF2/Cl2O/Br2O)
> reproduce the same warnings as before, just gating the geometry OVERLAY now (internal rows keep
> `ref_label` from Excel; only `has_geometry`/`predicted_label`/`Tx..Rz` are left null for those 4).
> **`py main.py --calibrate` (source="excel" default) re-froze `tau_TR=0.95, tau_S=0.90368,
> tau_B=0.17327`** — identical to the pre-2026-07-03 values (`derive_stretch_bend_thresholds` on the
> 69-molecule ideal subset: `ideal_stretch_n=41, ideal_bend_n=50, gap_width=0.73041`, matching the
> "gap width ~0.73" already on record). `confusion_matrix_stats()` reverted to the pre-2026-07-03 numbers
> too: **precision 1.0 all 4 categories; recall 1.0 T/R, 0.96466 bend (mixed 3.53%), 0.71739 stretch
> (mixed 28.26%)** — the SAME numbers fig:confusion was originally built from, not the 2026-07-03 disk-
> driven session's 0.9375/6.25% (that population-specific result is preserved for `source="gaussian"`,
> re-derivable any time by re-running with that source). `py main.py --figures` regenerated all 7
> figure PDF/PNG pairs from the excel-sourced CSV.
> **`source="gaussian"` re-verified standalone, unaffected:** `build_library_scores(source="gaussian")`
> still returns 303 rows / 25 molecules with `has_geometry` unconditionally `True`, identical to its
> 2026-07-03 behavior — confirmed by direct call, not assumed.
> **Tests:** `tests/test_excel_ingest.py`'s three tests that encoded source="gaussian"-specific
> invariants against the checked-in CSV (`test_library_row_count_and_molecule_count_match_disk_roster`,
> `test_all_rows_are_geometry_backed`, `test_frequency_mismatched_molecules_leave_internal_label_null_
> only` → renamed `..._gaussian`) now build a small in-memory `source="gaussian"` DataFrame via a new
> memoized `_gaussian_df()` helper instead of reading `_load()` (the checked-in CSV, which is
> source="excel" by default now) — decouples those tests from whichever source currently produced the
> committed golden. Two new tests added for the excel-sourced golden's own contract
> (`test_excel_sourced_default_has_full_population`, `test_excel_sourced_frequency_mismatched_molecules_
> keep_ref_label_without_geometry_overlay`). `tests/test_calibrate.py`'s two confusion-matrix tests
> (`test_confusion_matrix_precision_perfect_recall_explained_by_mixed_bucket`,
> `test_confusion_matrix_ideal_nonideal_recall_split`) had their pinned numbers restored to the
> source="excel" values above (docstrings explain both directions of the flip so a future reader isn't
> confused about which population a given number belongs to). **`py -m pytest tests/` 66/66 green**
> (up from 64; net +2 new tests, 0 removed, ~176s — the added `_gaussian_df()` rebuild in
> `test_excel_ingest.py` is the main new cost, run once per session via memoization).
> **score-validator PASS** (this session touched ONLY `src/excel_ingest.py`'s data-sourcing + a `source`
> passthrough param on `src/calibrate.py::run_calibration_pipeline` — `src/scoring.py`/`src/classifier.py`
> confirmed untouched, not even appearing in `git diff --stat`): water `tab:water` exact to 3 dp; `Σ s_AB
> == V_Stretch` holds across all 513 internal rows (max |diff| 4.8e-4, pure `.4f`-string-rounding noise
> pre-dating this session — the same formatting `_bond_string()` used in commit `149fc62`, not a
> regression); all 146 external rows reach exactly 1.000 on their own T/R axis; benzene EMIT 34/35 →
> `Tx*`/`Ty*` (vibration=SB), EMIT 36 clean `Tz` blind spot, reproduced exactly; H2S/SF2 eq:vscore spot-
> checks match the plan's own recorded values exactly; `confusion_matrix_stats()` recomputed live matches
> the numbers above to full precision (stretch recall 0.717391304347826, bend recall 0.9646643109540636).
> **For future sessions:** `source="excel"` is the default ONLY as an interim measure while the author's
> Gaussian file coverage is incomplete — do NOT treat this as a permanent reversal of the 2026-07-03
> "recompute from Gaussian" direction. When the author has dropped in enough `.log`/`.gjf` pairs that
> `source="gaussian"`'s population is close to the full library, flip the default back
> (`build_library_scores`'s `source` parameter default, `main.py --source`'s default, and
> `run_calibration_pipeline`'s default) and re-verify manuscript figures still hold under the
> Gaussian-recomputed numbers before doing so — do not flip defaults silently.
> **2026-07-03 (author decision — REVERSES the "ingest, don't recompute" locked decision below):**
> Author will supply real Gaussian `.log`/`.gjf` and EMIT files for the FULL hydride-library
> molecule set (currently only 24 of ~70 have on-disk geometry; see `data/logs/`), to be dropped into
> `data/logs/`, `data/EMIT/`, `data/gjf/` incrementally over time — no longer "off-server." Directive:
> prepare the program so **every score column** (`V_Stretch`, `freq`, `Tx..Rz`, per-bond `s_AB`,
> `delta_b_mean`/`rel_db`, `predicted_label`/`predicted_annotation`) is **recomputed from the raw
> Gaussian/EMIT files via the real engine** (`main.load_inputs` → `build_scorer_and_final` →
> `classify_all_modes`, `ModeScorer.score_bonds`), for every molecule that has geometry on disk —
> `data/vibrational-scoring-functions.xlsx` is to be used **ONLY** as a ground-truth lookup for
> `ref_label` ('type': bend/stretch) and the `ideal` tag, joined onto the recomputed rows by
> (molecule, mode index) with an engine-vs-Excel frequency sanity check gating the LABEL ATTACHMENT
> only (never gating score computation — scores are always the engine's now). This is dispatched to
> `lead-engineer` this session (see its own progress log below when done); **all subagents must treat
> this as the new locked decision** — `library_scores.csv`'s score columns are no longer Excel-sourced,
> full stop. The subset of Excel molecules still lacking on-disk geometry is expected to shrink to zero
> as the author drops in more files; **no further code changes should be needed** when that happens —
> the ingest loop is keyed off what's physically present in `data/logs/`+`data/gjf/`, not off Excel rows.
> **2026-07-03 (lead-engineer, disk-driven `excel_ingest.py` rearchitecture — DONE, formula-auditor +
> score-validator both PASS):** Executed the author directive immediately above. **Architecture chosen:**
> `src/excel_ingest.py` was rewritten around a disk-driven molecule loop, not an Excel-row loop —
> `discover_geometry_molecules()` lists every basename present in BOTH `data/logs/` (`.log`/`.out`) AND
> `data/gjf/` (`.com`/`.gjf`); `score_geometry_molecule()` runs the real engine
> (`main.load_inputs` → `build_scorer_and_final` → `src.classifier.classify_all_modes`,
> `ModeScorer.score_bonds()`) on each one to build every row (external T/R + internal "Vib i") with real
> `V_Stretch`/`Tx..Rz`/`predicted_label`/`s_AB`/`rel_db`/`delta_b_mean`; `attach_excel_labels()` then
> joins `ref_label`/`ideal` from Excel's `data_score` sheet onto the internal rows ONLY, gated by a
> whole-molecule frequency-agreement check (one mismatched mode nulls the WHOLE molecule's label join,
> not just that mode — same all-or-nothing guarantee as the prior architecture); `build_library_scores()`
> orchestrates all three and `run_ingest_pipeline()` is the headless entry point `main.py --library`
> already wired to. Old API (`ingest_internal_rows`, `attach_geometry_classification`,
> `_bond_string`) is gone; `resolve_log_basename`/`EXCEL_TO_LOG`/`EXCLUDED_MOLECULES` are unchanged and
> kept (still needed by `src/calibrate.py` and the new `resolve_excel_molecule_name` reverse-direction
> helper). **This is genuinely "no code changes needed to scale"**: the loop is keyed off
> `os.listdir(data/logs)` ∩ `os.listdir(data/gjf)`, so dropping in more `.log`/`.gjf` pairs grows the
> library automatically next time `--library`/`--calibrate` is re-run.
> `src/scoring.py`'s `ModeScorer._bond_contributions()`/`score_bonds()` gained a third, additive return
> value `rel_db` (signed per-bond relative bond-length change, `(|b_AB+Δd_AB|-|b_AB|)/|b_AB|`, the
> diagnostic fig:bondscores' x-axis needs) alongside the unchanged, already-audited `terms`/`sqdisps` that
> feed `Vscore()`/`s_AB` — confirmed purely additive, zero risk to eq:vscore/eq:bondscore (see audit
> below).
> **Bug found and fixed while validating, not just renaming:** the new `score_geometry_molecule()`
> initially built `s_AB`/`rel_db` bond-label strings from bare 1-based numeric indices (`"1-2:0.0342"`),
> silently losing the atom-symbol convention (`"C1-C2:0.0342"`) that `src/benzene_validation.py`'s C-C-
> vs-C-H bond-type diagnostic parses via `_is_cc_bond()` — this zeroed out `cc_fraction_of_V` for every
> mixed/worked-example benzene mode and broke 4 `test_benzene_validation.py` tests. Fixed by relabeling
> via `f"{scorer.atoms[idx].symbol}{idx+1}"` (verified against `data/gjf/benzene.com`'s own connectivity
> ordering — C1-C2/…/C6-C1 ring bonds, Ci-H(i+6) C-H bonds — matches exactly).
> **Final library size:** 25 molecules, 303 rows (157 internal + 146 external) — down from the prior
> Excel-driven baseline's ~69 molecules / 659 rows, a real and disclosed shrinkage until the author drops
> in more `.log`/`.gjf` pairs (not a bug; `discover_geometry_molecules()` currently returns exactly
> `water, benzene, co2_mp2_3-21g` + 22 hydride-library molecules, matching `data/logs/`+`data/gjf/` as of
> this session). Recalibrated end-to-end against this new population: `python main.py --calibrate` froze
> `tau_TR=0.95, tau_S=0.9036817451504533, tau_B=0.17326891344050538` (barely shifted in the 5th decimal
> from the stale pre-session values; `tau_sensitivity_sweep.csv` came back BYTE-IDENTICAL, a nice
> incidental reproducibility confirmation) — note the ideal-molecule stretch/bend sample sizes did drop
> (33/42, down from 41/50) since only 8 of the original 11 `tab:ideal` molecules (SnO2/TeH2/TeH4 still
> lack on-disk geometry) are in the current 25. `python main.py --figures` regenerated all 7 figure
> PDF/PNG pairs cleanly. **Confusion-matrix numbers moved in a notable, real direction**: because every
> library row now gets the FULL Algorithm 1 (Steps 2-4), not just Step 4's `vib_label` applied in
> isolation to Excel-only rows (the old architecture's necessary limitation for rows with no geometry),
> stretch recall rose from 0.717 to 0.9375 and bend recall is 0.9383 (both just under the 0.95 floor now,
> for the same reason — residual non-ideal external mixing landing in the MIXED bucket, 0 opposite-
> category crossings either way; `floor_met` still honestly `False`). `tests/test_calibrate.py`'s pinned
> confusion-matrix numbers were updated to match.
> **Test suite:** rewrote `tests/test_excel_ingest.py` from scratch for the new API (old file imported
> the retired `attach_geometry_classification` and failed at collection) — 18 tests covering
> `discover_geometry_molecules`/`resolve_excel_molecule_name`/`resolve_log_basename` (pure-logic, no
> slow xlsx I/O), a direct fast `score_geometry_molecule("water", ...)` check, a fully synthetic unit test
> of `attach_excel_labels()`'s frequency-gating contract (matched/mismatched/no-Excel-counterpart cases,
> including the "missing Excel row counts as a mismatch, not silently skipped" guarantee), and the
> checked-in `library_scores.csv` golden invariants (`Σ s_AB == s[V_S]`, all-external-clean, the 4
> known H2O/OF2/Cl2O/Br2O frequency-mismatch cases left `ref_label`-null but fully scored). Updated
> `tests/test_calibrate.py`'s frozen numbers for the new library size/composition (only 2 of its 10 tests
> needed number changes; the other 8 — including the benzene EMIT 34/35/36 and water-external targets —
> were unaffected). **`python -m pytest tests/` 60/60 green** (up from 53 pre-session; net +7 new tests).
> **Validation dispatched and both PASS:** formula-auditor confirmed `rel_db`'s formula matches its own
> stated definition exactly (same `delDisp` object as the existing `s_AB` numerator, so zero sign/
> convention drift risk), confirmed `src/excel_ingest.py`'s docstring formula is byte-identical to the
> code's, and confirmed `Vscore()`/`score_bonds()`'s existing eq:vscore/eq:bondscore math is untouched
> (two minor non-blocking observations noted: a pre-existing unguarded division in the `s_AB` numerator,
> and a scope-only doc-comment note for `EPS_NORM` — flagged for a future pass, not fixed here).
> score-validator confirmed no regression: water `tab:water` exact to 3 dp, `Σ s_AB == s[V_S]` holds
> library-wide (157/157 internal rows within the CSV's `.4f`-rounding bound; full-precision spot-checks
> diff ≤1e-16), all 146 external rows reach |score|=1.000 exactly, benzene EMIT 34/35 → `Tx*`/`Ty*`
> (vibration=SB) and EMIT 36 → clean `Tz` blind spot reproduced exactly, and the EMIT 2-vs-9 `s[R_y]`
> non-monotonicity (39% projected Ry → |s[Ry]|=0.143 vs. 14% projected Ry → |s[Ry]|=0.215) reproduced
> exactly.
> **For future sessions:** the ingest loop scales with zero code changes as more `.log`/`.gjf` pairs
> arrive — just re-run `python main.py --library && python main.py --calibrate && python main.py
> --figures` (or `python -m pytest tests/`, which reads the checked-in goldens and does not itself
> re-ingest). If the checked-in `library_scores.csv`/`thresholds.json` goldens are ever regenerated,
> re-check `tests/test_calibrate.py`'s pinned confusion-matrix numbers first — those are the ones most
> sensitive to population composition, everything else in the test suite is either pure logic or a
> structural invariant that holds regardless of N.
> **2026-07-03 (lead-engineer, CLI subcommands + README JCE-retraction fix):** Phase 5 checklist item
> **DONE** — `main.py` gained 5 new flags, each wiring an EXISTING pipeline function (no logic
> reimplemented): `--classify` (`run_classify_pipeline`, per-molecule, requires `-m` + `--mode`),
> `--emit-projection` (`run_projection_pipeline`, per-molecule, requires `-m`), `--library`
> (`src.excel_ingest.run_ingest_pipeline`, global, ignores `-m`, prints a slow-warning first), `--calibrate`
> (`src.calibrate.run_calibration_pipeline`, global, ignores `-m`), `--figures` (new
> `src.figures.regenerate_all()`, factored out of that module's former `if __name__ == "__main__":` block
> so `main.py` can call it directly instead of shelling out to `python -m src.figures`). Design: the
> default `-m <mol> --mode {normal,emit}` Step-1 path is UNCHANGED (verified no regression); new flags can
> combine freely with each other and with `-m`/`--mode` in one command (e.g.
> `-m benzene --mode emit --classify --emit-projection` runs both in sequence — tested, works); bad
> combinations fail loud with a clear message (`--classify` with no `--mode`; `--classify`/
> `--emit-projection` with no `-m`; `--emit-projection` against a molecule with no EMIT file) rather than
> silently doing nothing. **Every flag actually run and verified this session, not just written:**
> - `py main.py -m water --mode normal` (plain path, no new flags) → unchanged 9-row Tx..Vib3 table,
>   confirms no regression.
> - `py main.py -m water --mode normal --classify` → `data/results/water_normal_classified.csv` (9 rows,
>   fresh timestamp, correct S/B labels).
> - `py main.py -m water --emit-projection` → `water_EMIT_contributions.csv` +
>   `water_EMIT_projection_full.csv` (9 rows each, fresh timestamps).
> - `py main.py -m benzene --mode emit --classify --emit-projection` → all 3 outputs written in one
>   command (36 rows each).
> - `py main.py -m SF2 --emit-projection` (no EMIT file for SF2) → clean one-line error, no traceback.
> - `py main.py -m water --classify` (no `--mode`) / `py main.py --classify` (no `-m`) / `py main.py` (no
>   molecule, no flags) → each gives its own clear one-line error, exit without a stack trace.
> - `py main.py --library` → **659-row** `library_scores.csv` (ran in background, real ~630s cost
>   confirmed live, not just quoted from docs; 4 known frequency-mismatch skip warnings printed, matching
>   the pre-existing documented H2O/OF2/Cl2O/Br2O finding — nothing new or wrong).
> - `py main.py --calibrate` → `thresholds.json` re-derives the SAME frozen values already on record
>   (`tau_TR=0.95, tau_S=0.90368, tau_B=0.17327`) and a 190-row `tau_sensitivity_sweep.csv` — reproducible,
>   not drifted.
> - `py main.py --figures` → all 7 figure PDF/PNG pairs regenerated (fresh timestamps; content confirmed
>   unchanged in substance — only embedded PDF metadata bytes differ, same as any figure regen).
> All of the above CSV/JSON outputs came back **byte-identical** to what was already committed
> (`git status` showed zero diff for any `data/results/*` file after this whole session's runs) — a nice
> incidental reproducibility confirmation, not just a CLI-wiring check. `py -m pytest tests/` **53/53
> green** throughout (no test file touched; not required by the task, so none added — flag behavior was
> validated by actually running it, per the task's own instruction).
> **README.md also rewritten a SECOND time this date** (superseding the earlier-same-day rewrite logged
> just below, which is now WRONG and should not be treated as current): **direct author instruction**
> received mid-session — the *J. Chem. Educ.* (JCE) submission will be WITHDRAWN before JCC submission,
> so this is the FIRST scoring-functions paper, not a sequel to or extension of an earlier one. Removed
> EVERY mention of JCE/*J. Chem. Educ.*/"Paper I"/"Paper II" from `README.md` (grepped for
> `Educ|JCE|Paper I|Paper II`, zero matches confirmed). The repository is now framed as the reference
> implementation for exactly ONE manuscript: *"A Unified, Reference-Free Framework for Classifying the 3N
> Modes of Molecular Motion,"* in preparation for JCC. Citation section now names only JCC (no citation
> yet — "will be added once submitted", same honest placeholder convention as before, just for one paper
> instead of two). New "Classification, EMIT Projection, Library Calibration, and Figures (CLI flags)"
> section documents the exact tested commands above (not Python-snippet examples) — every command shown
> in the README is one that was actually run and verified this session, per the task's own instruction.
> **For future sessions:** JCC is the first/only paper for this repository going forward — do NOT
> reintroduce "Paper I (JCE)"/"Paper II (JCC)" dual-paper framing anywhere (README, docstrings, comments)
> unless the author explicitly reverses this decision again.
> **QUEUED, not started (author decision 2026-07-03, deliberately deferred until SI content stabilizes):**
> merge the two standalone SI documents (`JCC_SI_computational_cost.tex`, `JCC_SI_sensitivity.tex`) into
> ONE combined `JCC_SI.tex` before submission — Wiley/JCC convention is a single Supporting Information
> file, not several (SI is published as-supplied, not typeset by production). Target shape: one doc with
> continuous section/figure numbering (e.g. "S1. Elementary-Operation Derivation of Computational Cost",
> "S2. τ_TR Threshold-Sensitivity Analysis"), still growing (full score tables, Cartesian coordinates, and
> other content already promised in the main text's "Supporting Information" subsection still need to be
> added too) — do the merge once that content is in, not before, to avoid re-merging repeatedly.
> **2026-07-03 (lead-engineer, README.md rewrite):** Phase 5 checklist item done — `README.md` rewritten
> from scratch against a direct read of the current code (not assumption): confirmed `main.py`'s CLI is
> still only `-m`/`--molecule` + `--mode {normal,emit}` (no `--classify`/`--emit-projection`/`--library`/
> `--figures` subcommands exist yet — still Phase-5 outstanding work, checklist line left unchecked for
> that specific item); the classification/projection/library/calibration pipelines are reachable today
> only as importable Python functions (`main.run_classify_pipeline`/`run_projection_pipeline`,
> `src.excel_ingest.run_ingest_pipeline`, `src.calibrate.run_calibration_pipeline`), and figures via
> `python -m src.figures` (confirmed against `src/figures.py`'s `__main__` block, 8 PNG/PDF pairs).
> Overview section reframed around the current scope (Step-1 scoring + Steps 2-4 classification into 6
> categories, EMIT projection, library calibration, figures) and both papers described as an evolving
> two-manuscript codebase (Paper I / JCE still "submitted", Paper II / JCC "in preparation" — neither
> superseding the other, per `CLAUDE.md`'s "What this repository is"). Kept the existing accurate Step
> 1-5 basic-scoring workflow verbatim; added a new "Classification, EMIT Projection, Calibration, and
> Figures" section with verified code snippets. Citation section now lists both papers honestly (JCC has
> no citation yet, said so plainly). Noted `openpyxl` as an undeclared `requirements.txt` dependency of
> `src/excel_ingest.py` in the Installation section rather than silently editing `requirements.txt`
> itself (out of this task's scope). No code changed, no tests affected.
> **2026-07-03 (lead-author, whole-document precision editing pass):** Author reviewed the compiled PDF
> and requested a batch fix across `JCC_temp_LaTeXtemplate.tex`; all 8 items done, both SI docs still
> compile clean, main doc page count unchanged (32 pages before/after):
> 1. **Table 4 (`tab:water`) sizing** — root cause was a missing `\small` before its `\resizebox`'d
>    tabular (every sibling resizebox table has `\small`; water's didn't, so resizebox scaled the
>    normal-size table up instead of down). Added `\small`; now matches sibling tables' footprint.
> 2. **Paragraph indentation** — root cause found: `\captionof{figure}{...}` (used for the 7
>    non-floated pseudo-figures, `\begin{center}...\end{center}` + `\captionof`, not `\begin{figure}`)
>    calls the `caption` package's `\caption@parboxrestore@light`, which sets `\parindent\z@` with NO
>    surrounding group — so the FIRST `\captionof` in the document permanently zeroed `\parindent` for
>    every subsequent paragraph for the rest of the document. Fixed by wrapping each of the 7
>    `\captionof{...}\label{...}` calls in its own `{...}` group so the assignment is properly scoped
>    and restored. This was the actual bug; there were no stray `\noindent`s to remove.
> 3. **Math notation audit** — (a) all `\mathbf{}`/`\bm{}` converted to `\boldsymbol{}` throughout (main
>    tex only; SI left as-is except the `\Q` macro, item d); (b) bare math-italic T/R/V and $\nu_s$/$\nu_{as}$
>    wrapped `\mathrm{}` (labels), $Q$ left italic (variable); (c) single-atom-index subscripts
>    (`d_A`, `r_A`, `m_A`) converted to superscripts (`d^A` style) matching the SI convention,
>    bond-pair `_{AB}` subscripts left untouched; (d) `\hat{Q}` → `\boldsymbol{\hat{Q}}` everywhere in
>    main tex (inline, no new macro); SI's `\newcommand{\Q}{\hat{Q}}` redefined to
>    `\newcommand{\Q}{\boldsymbol{\hat{Q}}}` (single-point fix, ~15 uses); reworded the $\hat Q$
>    introduction to clarify it denotes a generic axis unit vector (in the $\hat\imath/\hat\jmath/\hat k$
>    sense), not literally "the x/y/z axes"; (e) bare `s[\mathrm{T}]`/`s[\mathrm{R}]` (generic/placeholder
>    uses, ~15 instances) → `s[\mathrm{T}_Q]`/`s[\mathrm{R}_Q]`, `Q` italic.
> 4. **Spacing after figure captions** — added `\medskip` after all 7 wrapped `\captionof` blocks
>    (same edit as item 2 also fixes this).
> 5. **Figure 8 placeholder** — shortened the long 3-panel numeric spec to two lines; box no longer
>    overfills.
> 6. **Figure 2 (water) overfill** — `NormalModes_Water.jpg` width reduced `0.95\columnwidth` →
>    `0.78\columnwidth`; the page-13 `Overfull \vbox (40.4pt too high)` is gone (confirmed in compile log
>    diff against pre-edit log). Table 4 shrinking (item 1) also freed vertical room on that page.
> 7. **Prose style pass** — reduced colon/semicolon/em-dash-heavy sentences to plainer ones across
>    Methodology, Results & Discussion (heaviest edit load, incl. the dense EMIT-stress-test and
>    computational-cost paragraphs), and Conclusions; left correctly-used colons (list intro) and the
>    Introduction's semicolon list (`...unambiguous; on the normal modes of benzene...; and, without
>    retuning...`) alone since those are legitimate uses.
> 8. **`\emph{}` pruning** — removed 10 non-load-bearing instances (routine adjectives/duplicates:
>    "unit", "relative", "degree", duplicate "consensus", duplicate "how much", "exact", "classify",
>    "Water"/"Benzene" dataset-name emphasis, "ideal"); kept ~34 that are genuine term-definitions,
>    logical-weight words (not/only/none), or the intentional italicized-topic-sentence device used for
>    the "design points"/"intrinsic limitations"/"three scoping statements" lists.
> Also copied `figure-builder`'s regenerated `fig_confusion.pdf` (Tx/Ty/Tz-style tick labels shortened to
> T/R) into `JCC/JCC_man_scoring/images/` per that agent's handoff. Verified: zero undefined refs/citations;
> the only remaining overfull-box warnings are the two pre-existing, unrelated ones (GTOC placeholder line
> ~102, and a ~1pt sub-visible one inside the Fig. 8 placeholder box) plus a benign `Float too large`
> notice for the Algorithm float, which drifts to the same end-of-document page in the ORIGINAL,
> unedited PDF too (confirmed by diffing text-extracted originals) — not introduced by this pass.
> **2026-07-02 (figure-builder, author follow-up fix):** `fig:confusion` tick labels for
> translation/rotation shortened to "T"/"R" (figure-local `_CONFUSION_LABEL` override in
> `plot_confusion_matrix`/`_confusion_heatmap`; `CATEGORY_LABEL` itself unchanged) and the
> now-unneeded tick-label rotation + several cramped font sizes restored to normal readable
> values; `data/figures/fig_confusion.{pdf,png}` regenerated, `py -m pytest tests/` 53/53 green.
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
> (`tests/test_classifier.py` new, 6 tests).
> ~~Phase 0: projection mass-weighting convention + Phase 2: `src/projection.py`~~ ✓ **DONE
> 2026-07-01** — convention **LOCKED: mass-weight by `sqrt(mass_A)` per atom (applied to all 3
> Cartesian components), applied only inside `src/projection.py`** (never touches the unweighted
> `s[T]/s[R]/s[V_S]` scores). Empirically confirmed via Gram-matrix check on benzene: Gaussian's
> printed normal modes and the raw EMIT eigenvectors are unit-normalized under the PLAIN Cartesian
> dot product but only mutually ORTHOGONAL under the mass-weighted one (plain off-diagonals up to
> 0.80; mass-weighted off-diagonals ≤4e-4, i.e. Eckart-Sayvetz orthogonality, not coincidence).
> `build_reference_basis()` mass-weights+renormalizes the `final_normal` pool (ideal T/R + real vib
> modes) into an orthonormal-to-noise-level `Q`; `project_emit()` computes `Θ̃=Q^T Θ` (eq:emitproj)
> on the similarly mass-weighted+renormalized raw EMIT eigenvectors, squares for fractional
> contributions, and buckets the vibration fraction into VS/VB/VMix by reusing `classifier.vib_label`
> on each real normal mode's own (unweighted) `s[V_S]` — no new stretch/bend boundary. `main.py`
> gained `run_projection_pipeline()`, writing `<mol>_EMIT_contributions.csv` (grouped, matches the
> pre-existing hand-derived ground-truth file's columns) and `<mol>_EMIT_projection_full.csv`
> (per-individual-reference-mode detail). **Validated against the pre-existing
> `data/results/benzene_EMIT_contributions.csv`: max abs deviation 3.1e-4 across all 36 EMIT modes ×
> 9 grouped columns** — reproduces EMIT 34/35/36 ≈76.7% translational and the EMIT 2 (38.7% Ry) vs
> EMIT 9 (14.1% Ry) inversion to spec. formula-auditor / score-validator dispatched — see Changelog
> for verdicts. `tests/test_projection.py` new (4 tests); 17/17 total green.
> ~~Phase-0 Excel column verification~~ ✓ **DONE 2026-07-01** — H2S/SF2 re-scored in-engine,
> `V_Stretch` matches `data_score`'s eq:vscore candidate column to ~5 sig figs (well within 3 dp); see
> Phase-0 checklist / Changelog.
> ~~Phase 3: `src/excel_ingest.py` + library classification + `src/calibrate.py` + confusion matrix~~
> ✓ **DONE 2026-07-02** (commit `149fc62`, previously built uncommitted in a session that hit its
> usage limit — recovered and committed intact this session, nothing lost). `src/excel_ingest.py`
> ingests `data_score`/`data_mode&bond` into `data/results/library_scores.csv`; the 25 geometry-backed
> molecules (water/benzene/CO2 + 22 hydride-library molecules with `.log`/`.gjf` pairs) get real
> Algorithm-1 predicted labels + ideal T/R external rows via `attach_geometry_classification`.
> **Discovered:** H2O/OF2/Cl2O/Br2O's Excel frequencies don't match this repo's own logs (94–386
> cm⁻¹ off — a different calc/basis, not rounding) — correctly left `has_geometry=False` for those
> molecules' internal rows rather than force-merged (external T/R rows unaffected, still attached).
> `src/calibrate.py`: `tau_S`/`tau_B` derived from the ideal-molecule library subset's non-overlapping
> stretch/bend `V_Stretch` populations (gap width ~0.73, zero-overlap by construction) →
> **frozen `tau_TR=0.95, tau_S=0.90368, tau_B=0.17327`** in `data/results/thresholds.json`;
> `tau_TR` swept 0.05–0.999 (step 0.005) against 25 molecules' real T/R ground truth (100% accuracy
> at every grid point) + benzene's 36 EMIT modes → plateau `(lo,hi)` with only ONE label change in the
> whole grid (benzene EMIT 19's Rz crossing at `tau_TR~0.335`); 0.95 falls inside the plateau, frozen
> there per spec. `confusion_matrix_stats()` (fig:confusion's numbers): **precision = 1.0 for all 4
> clean categories** (stretch/bend/translation/rotation never cross-contaminate); recall 1.0 for
> translation/rotation (exact completeness), ≥0.95 for bend, but **0.717 for stretch** — the 0.283
> shortfall lands entirely in the MIXED bucket (28.3%, non-ideal CoM-softening per B8.3), **zero** in
> BEND — reported honestly (`floor_met=False` overall), not forced to pass.
> `src/classifier.py`: `Thresholds.calibrated()` loads `thresholds.json` when present (falls back to
> the provisional 0.95/0.9/0.2 constants otherwise) and is now `classify_all_modes()`'s default;
> `tests/test_classifier.py` pins the explicit provisional `Thresholds()` so those goldens stay fixed
> across recalibration, `tests/test_calibrate.py` independently re-verifies the same targets (water
> externals CLEAN; benzene EMIT 34/35 flagged, EMIT 36 the documented blind spot) under the calibrated
> values. Also includes a degenerate-EMIT-block sanity check (benzene EMIT 1–9, one eigenvalue):
> deterministic run-to-run, exactly `n_T+n_R=6` modes total get an external-slot label with no over/
> under-assignment from the degeneracy (2 of the 6 winners — EMIT 6→Rx, EMIT 9→Ry — sit inside the
> block; the other 7 correctly get their OWN differing `vib_label`, not copies of one label).
> **34/34 tests green** (test count as currently collected in this repo; the "38/38" figure quoted in
> the previous entry above was the count reported that session — no tests were removed since, this is
> just the actual `pytest tests/` collection size, unaffected by the figures work below, which adds no
> new tests since `src/figures.py` is presentation-only).
> ~~Phase 3: 5 remaining figures (`fig:confusion`, `fig:bondscores`, `fig:boxplots`, `fig:modemixing`,
> `fig:sensitivity`) + Excel `box plots`/`CM`-sheet parity spot-check~~ ✓ **DONE 2026-07-02.** All 5
> added to `src/figures.py` (`plot_confusion_matrix`, `plot_bond_scores`, `plot_boxplots`,
> `plot_mode_mixing`, `plot_sensitivity`), each reading only already-computed
> `data/results/library_scores.csv` / `tau_sensitivity_sweep.csv` / `thresholds.json` /
> `confusion_matrix_stats()` — no scores recomputed. **Cross-figure style pass (coordinator directive,
> same session):** centralized a new `IDEAL_STYLE` dict (filled marker/box = ideal, hollow = non-ideal
> — matching the group's earlier undergraduate report H-02-598's own Figures 1-3 convention, confirmed
> by rendering that PDF) alongside the existing `CATEGORY_COLOR`/`CATEGORY_MARKER`/`CATEGORY_LABEL`
> dicts (unchanged, still fig:benzene's originals) plus `REF_LABEL_TO_CATEGORY`/
> `PRED_BUCKET_TO_CATEGORY` maps so the library's `ref_label`/predicted-bucket strings key into the
> SAME shared color/marker mapping fig:benzene already uses (gray=clean T/R, blue circle=bending,
> vermillion square=stretching, teal triangle=mixed stretch/bend, purple plus=mixed
> external+vibration) — one meaning per color/marker across the whole 6-figure set, not per-function.
> `fig:benzene` itself was NOT touched (no conflicting mapping arose). Each new figure's summary dict
> carries a `shared_categories` line naming exactly which shared encodings it reuses.
> **Parity spot-check vs Excel (this session):** the `box plots` sheet matches `library_scores.csv`
> EXACTLY (count/min/max/mean identical to displayed precision, all 4 ideal x stretch/bend groups) —
> this is fig:modemixing's/fig:bondscores' upstream source data, essentially a full pin, not just a
> spot-check. The `freq vs score` sheet (fig:boxplots' closest analog) matches on `V_Stretch`/`delta_b`
> distributions but its ideal-group counts run 2 short (48 vs our 50 bend; 39 vs our 41 stretch) —
> traced to SnO2's 4 internal rows having `freq=NaN` in the raw `data_score` sheet itself (a pre-existing
> Excel data gap, not a rounding issue); that sheet's own pivot table silently drops NaN-freq rows during
> its freq-grouping, while `excel_ingest.py` reads `data_score` directly and correctly retains them
> (SnO2's `V_Stretch`/`delta_b_mean` values are valid; only `freq` is missing) — `fig:boxplots` panel (a)
> naturally drops these 4 NaN-freq rows via per-group `.dropna()`, panels (b)/(c) keep them, matching the
> underlying Excel `data_score` sheet's actual content rather than the derivative pivot's incidental
> omission. Flagged for lead-author (SnO2 missing frequency in the source workbook), not silently
> patched. The `CM` sheet is NOT a confusion-matrix reference (it's an unrelated center-of-mass-
> conservation check by molecular shape, relevant to the OPTIONAL Phase-6 CoM item, not fig:confusion) —
> no Excel confusion-matrix sheet exists to parity-check against; fig:confusion's numbers were instead
> independently re-run against `confusion_matrix_stats()` this session and matched the task's stated
> targets exactly (precision 1.0 all 4 categories; recall 1.0 T/R, 0.9647 bend, 0.7174 stretch, 28.3%→
> mixed/0%→bend for the stretch shortfall).
> ~~Phase 6 item 1: systematic 36-mode benzene-EMIT flag precision/recall + library-external check~~
> ✓ **DONE 2026-07-02.** New `src/flag_validation.py`. Headline result: **precision=1.0, recall=5/17≈0.294**
> over all 36 modes (TP=5/FP=0/FN=12/TN=19) against a projection-derived ground truth (`M_ext=max(C2_Tx..
> C2_Rz)`; MIXED iff `0.05<M_ext<0.95`). Reproduces every named anchor (34/35→TP, 36→FN, 9→TP vs 2→FN).
> **Key finding:** the false-negative problem is WIDER than the single documented EMIT-36 blind spot — a
> second, structural cause is that Step 2's one-to-one Hungarian assignment can only ever flag exactly
> `n_T+n_R=6` of the 36 modes at all, so 11 further modes with genuine 7-39% external character
> (projection-verified) are never even flag-eligible. Library externals (146 rows/25 molecules): FP=0,
> confirmed not assumed. `tests/test_flag_validation.py` (6 tests) → 40/40 green. See Phase 6 checklist
> above for the full write-up; score-validator dispatched to independently confirm the counts.
> ~~Benzene normal-modes-vs-reference validation (Task A) + bond-contribution/degenerate-pair diagnostic
> (Task B) + `confusion_matrix_stats()` ideal/non-ideal recall split (Task C) + `excel_ingest.py`
> mismatch-gate robustness fix (Task D)~~ ✓ **DONE 2026-07-02.** Per
> `JCC/Scoring_Manuscript_Plan_2026-07-01.pdf`'s mandated Results & Discussion order, benzene's real
> NORMAL modes (not EMIT) vs. their literature/group-theory `ref_label` (already in
> `library_scores.csv` for `C6H6`) is now the manuscript's **PRIMARY** classification-vs-reference
> result — a non-circular ground truth (external literature assignment, unlike EMIT's threshold-cut
> "ground truth"). New `src/benzene_validation.py`: `benzene_normal_reference_detail/_summary` (Task A)
> reproduces the session's ad hoc numbers exactly by re-deriving them from `library_scores.csv` (no
> hardcoding) — **6/6 external T/R correct; 7/7 literature-stretch modes recalled (1.0); 18/23
> literature-bend modes recalled (0.7826), the other 5 (mode_index 13/14/19/23/24) landing in
> MIXED_STRETCH_BEND; ZERO crossings into the opposite clean category in either direction** (computed via
> an explicit `crossed_opposite` column, not eyeballed). `benzene_mixed_bond_diagnostic` (Task B) parses
> those 5 modes' per-bond `s_AB` (already in `library_scores.csv`) and computes: C-H bond contributions
> are ~0 (`ch_total` 0.0088-0.0100 total across all 6 C-H bonds per mode, `cc_fraction_of_V` >0.95 for
> all 5); mode 19's C-C contribution is essentially uniform across all 6 ring bonds (coefficient of
> variation <0.02); modes 13/14 and 23/24 are DETECTED (not assumed) as near-degenerate pairs (freq
> splitting 0.0157/0.0292 cm⁻¹, well inside a 1 cm⁻¹ tolerance) with a computed strong NEGATIVE Pearson
> correlation between their 6-bond C-C `s_AB` vectors (-0.969 / -0.997) — the "complementary alternating
> pattern" claim is now a number, not an eyeballed observation; this is exactly the D6h degenerate
> E-type-pair signature. Both functions raise (fail-loud) if the geometry merge is incomplete or if no
> mixed modes exist, rather than silently validating partial/stale data. Outputs:
> `data/results/benzene_normal_reference_{detail,summary}.csv`,
> `data/results/benzene_mixed_{bond_diagnostic,degenerate_pairs}.csv`.
> `tests/test_benzene_validation.py` (9 new tests) pin every number above.
> **Task C:** `confusion_matrix_stats()` (`src/calibrate.py`) now additionally reports
> `per_category[cat]["recall_ideal"]`/`["recall_nonideal"]` (+ `n_ref_ideal`/`n_ref_nonideal`) using the
> SAME tier masks `src/figures.py::plot_confusion_matrix` already applies ad hoc for `fig:confusion`'s
> two-tier layout (commit `424a666`) — purely additive, existing pooled `precision`/`recall` keys
> unchanged, `src/figures.py` NOT touched. Verified computationally (not just asserted):
> `recall_ideal["stretch"] == 1.0` and `recall_ideal["bend"] == 1.0` EXACTLY, as required by
> construction since `tau_S`/`tau_B` are literally the ideal population's own min/max; `recall_nonideal`
> reproduces the figure's numbers (bend 0.9571, stretch 0.6561). New test in `tests/test_calibrate.py`.
> **Task D:** `excel_ingest.py::attach_geometry_classification`'s `if m is None: continue` branch (a row
> whose "Vib i" name resolves to no engine mode at all) now appends a `(mode_index, None, freq)` sentinel
> to `mismatches` instead of silently skipping — so a molecule hitting this case (not currently known to
> occur for any molecule, per the module's own docstring, but not proven impossible) is excluded and
> reported like any other frequency mismatch, matching the documented all-or-nothing merge guarantee. New
> test in `tests/test_excel_ingest.py` (synthetic bogus `mode_index=999` on water/H2O).
> **Note (per author instruction):** the EMIT systematic 36-mode confusion matrix (`src/flag_validation.py`,
> commit `65f4e92`, Phase 6 item above) remains built and tested but is EXCLUDED from the manuscript by
> author decision — its "ground truth" is a threshold cut (`0.05 < M_ext < 0.95`) on continuous,
> genuinely-mixed EMIT projection fractions, which is circular reasoning for an accuracy claim; benzene
> EMIT is now framed purely as an extreme/rare edge-case stress test (EMIT 34-36 flag behavior, EMIT 2-vs-9
> inversion), never as the paper's systematic classification-accuracy evidence. That role now belongs
> entirely to the benzene-normal-modes validation documented above. 51/51 tests green
> (`py -m pytest tests/`).
> **LABEL-VOCABULARY RENAME (2026-07-02, author-directed, IN PROGRESS as of this update):** classifier
> output labels are now short and axis-specific, replacing the old generic constants. Mapping:
> `CLEAN_TRANSLATION`/`CLEAN_ROTATION` (any axis) → the specific Step-2 slot name the mode actually won
> (`"Tx"`,`"Ty"`,`"Tz"`,`"Rx"`,`"Ry"`,`"Rz"`); `MIXED_EXTERNAL_WITH_VIBRATION` → the same slot name with
> a trailing `"*"` (`"Tx*"` etc.); `STRETCHING`/`BENDING`/`MIXED_STRETCH_BEND` → `"S"`/`"B"`/`"SB"` (the
> Python constant NAMES in `src/classifier.py` are unchanged, only their string VALUES). New helper
> predicates in `src/classifier.py`: `is_external_label()`, `external_axis()`, `is_clean_external()`,
> `is_mixed_external()`, `is_translation()`, `is_rotation()`, `classification_bucket()` (the last maps
> ANY classification label to one of 6 semantic buckets — `"translation"/"rotation"/"stretch"/"bend"/
> "mixed"/"mixed_external"` — for anything needing the old coarse-category behavior, e.g.
> `src/figures.py`'s `CATEGORY_COLOR`/`CATEGORY_MARKER`/`CATEGORY_LABEL`, now keyed by these 6 bucket
> names instead of the old label strings). Mixed-external annotations dropped the now-redundant
> `dominant_external=<slot>` prefix (the slot is already IN the classification string) — annotation is
> now just `"vibration=<S|B|SB>"`.
> **Recovery note for a future session finding this mid-flight:** a `lead-engineer` dispatch did the
> core rename (`classifier.py`, `calibrate.py`'s bucket logic, `projection.py`, `benzene_validation.py`,
> `flag_validation.py`) correctly but was cut off by a session-limit hit before finishing propagation —
> left the repo NOT importable (`tests/test_classifier.py`/`test_calibrate.py` failed to import
> `CLEAN_TRANSLATION`, which no longer exists) and nothing committed. The main session finished the
> propagation directly (not re-dispatched, to conserve session budget): fixed `src/figures.py`'s category
> dicts, updated every test file's assertions to the new scheme, and re-ran the full regeneration chain
> (`main.run_classify_pipeline`/`run_projection_pipeline` for water/benzene/CO2 → `src/excel_ingest.
> run_ingest_pipeline()` [genuinely slow: ~630s just to parse the workbook's two sheets via openpyxl,
> not a hang — this workbook's Data-Validation extension parses slowly; NOT the reported "~1 minute" the
> tests' docstrings claimed, that estimate was wrong/optimistic] → `src/calibrate.run_calibration_pipeline()`
> → `src/benzene_validation.run_benzene_normal_validation()`/`run_benzene_bond_diagnostic()` →
> `src/flag_validation.run_flag_validation_pipeline()`). **If resuming and this note is still here
> unresolved:** check `git status`/`git log` in `Github/scoring-functions/` first — the rename may have
> finished and been committed after this note was written (check for a commit message mentioning the
> label rename after `8c0e27c`), or the CSV regeneration chain above may need re-running if it didn't
> finish (`py -m pytest tests/` will fail loudly with stale-label assertion mismatches if so, not a silent
> corruption). Do NOT re-dispatch a fresh `lead-engineer` for this without checking working-tree state
> first — the pattern this session (and the 2026-07-02 Phase-3 session before it) both confirm:
> interrupted-by-limit work is usually recoverable, not lost.
> ~~`figure-builder`: re-render all 6 figures after the label rename~~ ✓ **DONE 2026-07-02** — ran
> `py -m src.figures`, visually inspected all 6 regenerated PNGs (legend/annotation text unchanged:
> `CATEGORY_LABEL` string VALUES like "clean translation"/"bending" never changed, only the dict keys
> did, so the rendered content is pixel-identical to the pre-rename PDFs — confirmed via `git diff
> --stat` showing 0 insertions/deletions, same byte counts, only embedded PDF metadata differs). No
> stale/raw label strings (e.g. `"Tx*"`, `"SB"`) found anywhere a descriptive English label was
> expected. 51/51 tests still green. Committed + pushed.
> ~~`lead-author`: fix literal old-label quotes in the `.tex`~~ ✓ **DONE 2026-07-02** — every
> `\textsc{clean_translation}`-style literal quote (including ones in `fig:flowchart`/`alg:classify`
> that would otherwise have contradicted the Results section) updated to the new short symbols
> (`\texttt{Tz}`, `\texttt{S}`/`\texttt{B}`/`\texttt{SB}`, `\texttt{Tx*}`), consistent with the
> document's existing `\texttt{}` convention for code-like tokens. Compiles clean, 28 pages.
> **GAP FOUND after a direct read of the manuscript (2026-07-02, main-session read, not just an agent
> report) — the real reason the manuscript still read "like the old version":** the "Ideal vs.
> non-ideal" section and the "Benzene normal modes" section immediately after it both lead with
> confusion-matrix-style recall numbers (95.7%/65.6%, then 6/6/7/7/18-23) — structurally repetitive,
> not a content bug. Per `JCC/Scoring_Manuscript_Plan_2026-07-02.pdf` (the current authoritative
> structure doc), benzene-normal-modes is supposed to serve a DIFFERENT narrative role: a descriptive
> worked-example gallery (full T/R/S/B/SB classification across the whole frequency range, with named
> illustrative modes: a ring-breathing mode, an SB mode, a C-H stretch), not a second accuracy report.
> The existing validation numbers are good evidence and should stay, but as supporting detail under
> the descriptive frame, not the section's lead.
> **CONFIRMED PLAN (author consulted 2026-07-02, both open questions resolved):**
>   1. **H2O section:** wire in `JCC/JCC_man_scoring/images/NormalModes_Water.jpg` (already supplied by
>      the author — a labeled 3x3 grid of Tx/Ty/Tz/Rx/Ry/Rz/δ/νs/νas displacement-vector diagrams,
>      matching `tab:water`'s notation exactly) next to `tab:water`. Pure `lead-author` text/figure-wiring
>      task, no new figure-builder computation needed (the image is a hand-made illustrative diagram).
>   2. **Ideal vs. non-ideal section:** CONFIRMED already aligned (two-tier confusion matrix, the
>      non-true-reference caveat, and the one-central-atom/delocalization scope note are all present and
>      correctly worded) — no further action. `fig:confusion` STAYS COMBINED as one 4-panel figure
>      (author confirmed 2026-07-02: splitting into separate `fig:confusion`/`fig:retention` figures, as
>      the deck's bullet list loosely suggests, would separate closely-related panels with no reader
>      benefit — keep the current 2x2 layout).
>   3. **Benzene normal modes — the real work, strict dependency order, do not parallelize:**
>      (a) ~~`lead-engineer` FIRST: identify the ring-breathing/C-H-stretch mode indices~~ ✓ **DONE
>      2026-07-02.** New `src/benzene_validation.py` Task E: `benzene_worked_examples()` /
>      `run_benzene_worked_examples()`. Primary criterion (author-confirmed mid-session, overriding an
>      earlier bond-uniformity-first heuristic): among C6H6's 7 `S` (STRETCHING)-labeled real normal
>      modes there is a clean ~2200 cm⁻¹ frequency gap; the lowest-frequency one with `V_Stretch`≈1.000
>      is the ring-breathing mode, the highest-frequency one is the representative C-H stretch.
>      **Identified: mode 12 (992.5825 cm⁻¹ — matches the literature ~992 cm⁻¹ ring-breathing assignment
>      almost exactly) = ring-breathing; mode 30 (3223.172 cm⁻¹, the highest-frequency `S` mode) =
>      representative C-H stretch.** Per-bond `s_AB` supporting evidence (reusing Task B's parsing
>      helpers, extended with a new `_CH_BONDS` tuple): mode 12 — `V_Stretch=1.0000`, 99.1% of it on the
>      6 C-C ring bonds (`cc_fraction_of_V=0.991`), uniform to CV=0.0024, C-H total negligible (0.0090);
>      mode 30 — `V_Stretch=1.0000`, 99.3% of it on the 6 C-H bonds, C-C total negligible (0.0072), and
>      **not** part of a near-degenerate pair (unlike modes 26/27 and 28/29 among the same `S` cluster,
>      which sit 0.017/0.016 cm⁻¹ apart — mode 30 is isolated, the natural non-degenerate representative
>      pick, computed via a `near_degenerate_partner`/`partner_freq_diff` check, not assumed). Mode 19
>      (SB, 1319.27 cm⁻¹) untouched, as instructed. Output: `data/results/benzene_worked_examples.csv`
>      (2 rows: mode_index, freq, V_Stretch, role, cc_total, ch_total, cc_fraction_of_V, cc_min/max/cv,
>      ch_min/max/cv, near_degenerate_partner, partner_freq_diff). Two new regression tests in
>      `tests/test_benzene_validation.py` pin mode 12/30 as manuscript claims (same protection level as
>      the existing 13/14/19/23/24 pins). **53/53 tests green.**
>      (b) ~~`figure-builder` SECOND: build a NEW score-vs-frequency figure for benzene's 36 normal
>      modes~~ ✓ **DONE 2026-07-02.** New `src/figures.py::plot_benzene_normal_modes` (standalone, no
>      `fig:` label yet — pending (c)); reads `data/results/benzene_normal_classified.csv` directly
>      (single source of truth for Freq/V_Stretch/label — `benzene_worked_examples.csv` was consulted
>      for context but the plotted numbers come live from the classified CSV, not re-derived). Single-
>      column scatter of all 36 modes (6 external T/R + 30 internal), full T/R/S/B/SB classification via
>      the shared `CATEGORY_COLOR`/`CATEGORY_MARKER`/`CATEGORY_LABEL` dicts (same mapping as fig:benzene
>      panel (a) — routes each raw `label` through `classification_bucket()` first), calibrated
>      `tau_S`/`tau_B` threshold lines via `Thresholds.calibrated()`. Modes 12/19/30 get an enlarged
>      black-outlined marker (white halo behind it first, since mode 30 sits in a tight ~40 cm⁻¹-wide
>      C-H-stretch cluster with modes 25-29 and would otherwise show a sliver of its un-highlighted
>      neighbor peeking out) plus a dashed-leader callout box with frequency + `s[V_S]` readout, same
>      idiom as fig:benzene panel (b)'s EMIT 34-36 callout. Legend anchored between the tau_B/tau_S lines
>      (`loc="center left", bbox_to_anchor=(0.0, 0.52)`) rather than the default "upper left", which
>      otherwise put the tau_S dashed line straight through the legend text (this distribution's empty
>      band differs from fig:benzene panel (a)'s, so "center right" doesn't clear it the same way).
>      `plot_benzene_stress_test`/`fig:benzene` itself is untouched (diff to `src/figures.py` is purely
>      additive). `data/figures/fig_benzene_normal.{pdf,png}`. 53/53 tests still green.
>      (c) ~~`lead-author`: rewrite the benzene-normal-modes section's narrative~~ DONE 2026-07-02 --
>      the dispatched agent completed the rewrite, then dispatched its own nested review sub-agent
>      before finalizing and was manually stopped by the user while idle waiting on it. Resuming it via
>      SendMessage failed (once a user stops a task the harness treats it as cancelled). The
>      coordinating session verified the already-written .tex content directly (Read, not just a
>      report), confirmed it was complete and good, and finished the remaining mechanical steps itself:
>      confirmed fig_benzene_normal.pdf was already copied into JCC/JCC_man_scoring/images/, recompiled
>      with latexmk (OneDrive-safe: copied to scratch, compiled there, copied the PDF back), verified
>      zero undefined references/citations, copied the final PDF back. 31 pages (up from 29). Section
>      now leads with the three worked examples (mode 12 ring-breathing 992.6 cm-1, s[V_S]=1.000; mode
>      30 C-H stretch 3223.2 cm-1, s[V_S]=1.000, inverted bond pattern, illustrating frequency-
>      independence; mode 19 mixed S/B 1319.3 cm-1, the honest-mixed-character case), then explicitly
>      pivots ("These three modes anchor the systematic validation behind them") into the pre-existing
>      6/6, 7/7, 18/23 numbers and tab:benzenemixed as supporting evidence -- exactly per plan.
>      **Lesson:** if a dispatched agent's status comes back killed/stopped-by-user rather than
>      completed, verify its file-level work directly before discarding or blindly retrying -- here the
>      substantive work was already good, only the housekeeping tail (compile + copy-back + bookkeeping)
>      was missing.
>   4. **Benzene EMIT modes:** CONFIRMED already aligned (τ_B-specific two-gate framing and the "extreme
>      case that rarely occurs" limitation note are both explicit) — no further action.
> **Benzene-normal-modes 3-step sequence (identify -> figure -> narrative) now fully complete.**
> ~~5. `figure-builder`: simplify `fig:benzene` to single-panel (drop the now-redundant normal-mode
>    panel)~~ ✓ **DONE 2026-07-02.** `plot_benzene_stress_test` (`src/figures.py`) dropped its former
>    panel (a) (`s[V_S]` vs. frequency, benzene's 36 normal modes) -- that content is now strictly
>    subsumed by the standalone `plot_benzene_normal_modes`/`fig_benzene_normal` figure (same data plus
>    named worked-example callouts, built in the sequence above), and the manuscript's "stress test on
>    benzene EMIT modes" prose never referenced panel (a)/"36 normal modes" at all, only the EMIT 2/9 and
>    EMIT 34-36 content. Having real NORMAL-mode data inside a figure captioned around the EMIT stress
>    test was also conceptually confusing independent of the redundancy. `fig:benzene` is now a single
>    `fig, ax = plt.subplots(figsize=(3.8, 3.6))` panel containing only the former panel (b) (EMIT score
>    vs. projected normal-mode contribution, EMIT 2/9 R_y-inversion diamond callout + EMIT 34-36 flagged-
>    external star callout) -- unchanged content, just no longer alongside the dropped panel. Removed the
>    now-unused `normal_csv` parameter (no caller passed it positionally/by keyword outside this module's
>    own `__main__` block, which calls with no args -- nothing else to fix) and the `panel_a_*`/`panel_b_*`
>    summary-dict key prefixes (now unprefixed, e.g. `emit2_Ry_score`; no test asserted on the old keys).
>    `.tex` caption for `fig:benzene` (`JCC_temp_LaTeXtemplate.tex`) checked and needs NO change -- it
>    already only describes EMIT/flagged-external/R_y-inversion content, never "36 normal modes"/panel
>    (a). Regenerated `data/figures/fig_benzene.{pdf,png}` and copied the PDF into
>    `JCC/JCC_man_scoring/images/fig_benzene.pdf`. 53/53 tests still green. `plot_benzene_normal_modes`/
>    `fig_benzene_normal` itself untouched.
> ~~Two queued figure fixes (author-flagged 2026-07-02): fig:confusion footer-text readability +
> S/B/SB short-notation adoption~~ ✓ **DONE 2026-07-02** (same session, `lead-engineer`).
> **Fix 1 (footer readability):** `plot_confusion_matrix`'s whole-figure `fig.text()` footer sentence
> ("Non-ideal tier (n=422): 0% of bend or stretch..." at 6.8pt on a 7.4x6.6in canvas -- ~5.7pt effective
> once LaTeX's `\includegraphics[width=0.95\columnwidth]` rescales it to ~6.2in, below a readable floor)
> is REMOVED from the raster entirely. The exact sentence is now computed in-function and returned as
> `summary["nonideal_footer_text"]` (printed by `__main__`, so it's easy to copy into the LaTeX
> `\captionof{figure}{...}` text at normal caption font size -- lead-author's job, not done here). Exact
> text (n/percentages recomputed live, so copy this from a fresh run if the underlying data ever
> changes): *"Non-ideal tier (n=422): 0% of bend or stretch reference-labeled modes crossed to the
> OPPOSITE clean category (bend->stretch=0.0%, stretch->bend=0.0%); 100% of the non-retained remainder
> lands in the mixed bucket."* (ASCII `->` used instead of a unicode arrow -- the original had `→`, which
> crashes Python's `print()` under Windows' default cp1252 stdout encoding once the text moved into a
> printed summary dict; reads identically fine in a LaTeX caption). **Checked the other 6 figures for
> the same `fig.text()`-footer-annotation pattern** (the task's explicit ask, not just fig:confusion):
> found ONE more instance -- `plot_benzene_stress_test`/`fig:benzene` had a smaller in-axes prose note
> ("non-monotonic by design (flag mechanism, not a parity check)", 6.3pt on a 3.8x3.6in canvas, drawn via
> `ax_b.text(..., transform=ax_b.transAxes)` rather than `fig.text()` but the same architectural problem
> -- explanatory prose baked into the raster at sub-readable size) -- fixed identically, exposed as
> `summary["nonmonotonicity_note"]`. The enlarged-marker/dashed-leader callout boxes in `fig:benzene`
> (EMIT 34-36) and `fig_benzene_normal` (modes 12/19/30) are a different, fine idiom per the task's own
> carve-out (short per-point data readouts, not paragraph-length figure-level prose) -- left untouched.
> No other figure (`fig_bondscores`, `fig_boxplots`, `fig_modemixing`, `fig_sensitivity`) has any
> `fig.text()` call at all (grepped to confirm, not just skimmed).
> **Fix 2 (S/B/SB notation):** `CATEGORY_LABEL`'s `"bend"`/`"stretch"`/`"mixed"` entries changed from
> the spelled-out `"bending"`/`"stretching"`/`"mixed stretch/bend"` to `"B (bending)"`/`"S (stretching)"`/
> `"SB (mixed S/B)"`, matching `src.classifier`'s actual output strings (`STRETCHING="S"`, `BENDING="B"`,
> `MIXED_STRETCH_BEND="SB"`) and the manuscript's `\texttt{S}`/`\texttt{B}`/`\texttt{SB}`/`tab:benzenemixed`
> notation. `translation`/`rotation` entries (`"clean translation"`/`"clean rotation"`) and
> `mixed_external` (`"mixed external+vibration"`) are UNCHANGED -- the manuscript has no single-letter
> T/R bucket symbol to match (its actual short labels are axis-specific, `"Tx".."Rz"`) and no short
> symbol for external+vibration mixing either, so nothing to rename there.
> **Gloss-form decision (author's call per the plan, made explicitly):** chose the self-contained
> in-figure gloss ("S (stretching)", not bare "S" relying on a once-per-caption gloss) because these
> figures are also viewed as standalone PNGs outside the compiled manuscript, and adding a caption gloss
> is `lead-author`'s job (a `.tex` edit), not something this session's `figures.py`-only change could
> guarantee would exist. Applied everywhere `CATEGORY_LABEL` is read (`fig_confusion`'s heatmap ticks +
> precision/recall bars, `fig_modemixing`'s legend, `fig_benzene_normal`'s legend) AND to the two figures
> that spell out the same vocabulary WITHOUT going through the dict: `fig_bondscores`'s hand-built
> `Line2D` legend (`"S (stretching), ideal"` etc., since that legend also encodes ideal/non-ideal, which
> `CATEGORY_LABEL` alone can't) and `fig_boxplots`'s `group_labels` tick-label list. `fig_benzene` itself
> does NOT use `CATEGORY_LABEL` (checked, per the task's own suggestion) -- its background/highlight
> labels ("other EMIT modes", "EMIT 2/9", "EMIT 34-36") are a different, unrelated legend vocabulary, so
> nothing to change there. **Regression found + fixed by rendering, not just eyeballing the diff:**
> `fig_boxplots`'s longer gloss-form tick labels ("B (bending)"/"S (stretching)" vs. the old bare
> "bending"/"stretching") collided horizontally at the existing tick spacing -- fixed by rotating those
> tick labels 30°/`ha="right"` (same idiom `fig:confusion`'s heatmap ticks already use), re-rendered,
> confirmed collision-free.
> **All 8 figures regenerated** (`py -m src.figures`; the module's `__main__` runs all 7 named figures +
> the no-label-yet benzene-normal-modes gallery = 8 PNG/PDF pairs total): `fig_benzene`, `fig_benzene_normal`,
> `fig_confusion`, `fig_bondscores`, `fig_boxplots`, `fig_modemixing`, `fig_sensitivity` (this last one
> untouched by either fix -- no `CATEGORY_LABEL`/`fig.text()` use -- regenerated anyway for a clean,
> consistent `data/figures/` snapshot, confirmed pixel-appropriate by inspection, no unexpected diff).
> Visually re-verified every changed PNG (not just trusted the code): confusion-figure footer sentence
> confirmed GONE from the image; `fig_benzene`'s bottom-right prose note confirmed GONE; S/B/SB gloss
> labels confirmed rendering correctly (readable, non-overlapping) in `fig_confusion`, `fig_bondscores`,
> `fig_boxplots` (post-rotation-fix), `fig_modemixing`, `fig_benzene_normal`.
> `py -m pytest tests/` **53/53 still green** throughout (no figure-specific tests exist; this module is
> presentation-only, per its own docstring -- confirmed nothing else regressed).
> ~~Manuscript-side half of both queued figure fixes (wiring `nonideal_footer_text`/
> `nonmonotonicity_note` into `.tex` captions + refreshing the embedded images)~~ ✓ **DONE
> 2026-07-02** (same day, `lead-author`). `fig:confusion`'s `\captionof{figure}{...}` (panel (d)
> sentence) now states the exact `nonideal_footer_text` content in place of the old "the figure footer
> records..." cross-reference (which pointed at prose no longer baked into the raster): *"Non-ideal
> tier ($n=422$): 0\% of bend or stretch reference-labeled modes crossed to the opposite clean category
> (bend$\to$stretch $=0.0\%$, stretch$\to$bend $=0.0\%$); 100\% of the non-retained remainder lands in
> the mixed bucket."* (real `$\to$` arrow in the `.tex`, not the engineer's console-safe ASCII `->`,
> per the task's own instruction). `fig:benzene`'s caption gained a trailing clause carrying
> `nonmonotonicity_note`'s substance as a full sentence rather than a dropped-in fragment: *"The EMIT
> score-vs.-contribution relationship plotted here is non-monotonic by design (a flag mechanism, not a
> parity check), so the non-monotonic points should be read as the expected behavior of a threshold
> test, not as noise or an error."* Copied the freshly regenerated PDFs (S/B/SB gloss legends, footer/
> note prose removed from the raster) from `data/figures/` into `JCC/JCC_man_scoring/images/`: `fig_benzene`,
> `fig_benzene_normal`, `fig_bondscores`, `fig_boxplots`, `fig_confusion`, `fig_modemixing`,
> `fig_sensitivity` (all 7 that changed bytes in commit `36d24c7`; `fig_sensitivity` copied too for a
> consistent snapshot even though its bytes were identical). Recompiled via the OneDrive-safe
> scratch-dir latexmk workflow: **32 pages before -> 32 pages after** (no page-count shift), zero
> undefined references/citations, no new LaTeX warnings (the sole pre-existing warning, `caption
> Warning: \setcaptiontype ... outside box or environment on input line 616`, is the unrelated
> water-modes figure, present before this edit too). PDF copied back to
> `JCC/JCC_man_scoring/JCC_temp_LaTeXtemplate.pdf`.
> ~~S/B/SB in-figure gloss reverted to bare symbols (author visual review)~~ ✓ **DONE 2026-07-02**
> (same day, `lead-engineer`). The gloss-form decision recorded just above ("S (stretching)" etc.,
> chosen because these PNGs are also viewed standalone) was reviewed by the author once rendered and
> judged TOO LONG inside the figures themselves. Reverted: `CATEGORY_LABEL`'s `"bend"`/`"stretch"`/
> `"mixed"` entries are back to bare `"B"`/`"S"`/`"SB"` (no parenthetical); `fig_bondscores`'s hand-built
> `Line2D` legend back to `"S, ideal"`/`"S, non-ideal"`/`"B, ideal"`/`"B, non-ideal"`; `fig_boxplots`'s
> `group_labels` back to bare `"B"`/`"S"`. The gloss will instead be added ONCE per figure in the LaTeX
> caption text (`lead-author`'s job), not baked into the raster. Reverting also reopened the tick-label-
> collision question `fig_boxplots` had fixed with a 30°/`ha="right"` rotation: bare single-character
> `"B"`/`"S"` have no collision risk at the existing spacing, so the rotation is no longer needed and was
> removed (ticks now horizontal/unrotated) -- reads cleaner, verified by rendering.
> `translation`/`rotation`/`mixed_external` entries untouched (already long-form, not part of this
> gloss). Regenerated the 5 figures that read `CATEGORY_LABEL`/these hand-written legends:
> `fig_confusion`, `fig_modemixing`, `fig_benzene_normal`, `fig_bondscores`, `fig_boxplots`.
> `fig_benzene` confirmed (grep) to never use `CATEGORY_LABEL` -- left untouched, per the prior entry's
> note. Visually re-verified all 5 regenerated PNGs: bare `S`/`B`/`SB` only, no leftover collision, no
> stray gloss text. `py -m pytest tests/` **53/53 still green**.
> **Next:** Phase 5 (`reproduce.py` orchestrator wiring; SI Cartesian-geometry export; graphical TOC)
> and the remaining two Phase 6 recommended items (mixed-SB bucket CoM-argument half; leave-one-
> molecule-out τ evaluation). The `fig:modemixing` irrep-degeneracy sub-panel gap remains BLOCKED
> pending confirmation, not built. Consider an `expert-reviewer-jcc` pass on the manuscript now that
> the label rename, benzene restructuring, water-figure/δ-notation fixes, and both figure-caption
> wirings have all landed -- this was requested earlier and deliberately deferred until things
> stabilized; they have. **Still outstanding:** the LaTeX caption side of the S/B/SB gloss revert --
> a caption sentence spelling out "S = stretching, B = bending, SB = mixed" once per figure caption
> (`lead-author`'s job) has not yet been written; the figures currently show bare symbols with no gloss
> anywhere until that lands. No other figure fixes queued as of this update.
> **After each step:** `py -m pytest tests/` should stay green.
> **Resilience rule:** work in small increments; after each, tick the checkbox here + below and `git commit`
> so the plan-in-git always reflects true state. Manuscript `.tex` is outside the repo (not committed).
> **2026-07-02 (small fix):** `fig_benzene_normal`'s legend showed "clean translation"/"clean rotation" as
> two separate entries despite both sharing the identical gray-X color/marker — merged to one "clean T/R"
> legend entry (local dedup-key fix in `plot_benzene_normal_modes`; `CATEGORY_LABEL`/`CATEGORY_COLOR`/
> `CATEGORY_MARKER` untouched, since other figures' tick labels still need "translation"/"rotation" kept
> distinct). Regenerated `fig_benzene_normal.{pdf,png}` only; 53/53 tests still green.
> **2026-07-02 (Results & Discussion structural pass, `lead-author`, per `Scoring_Manuscript_Plan_2026-07-02.pdf`
> p.4):** Five structural changes to `JCC_temp_LaTeXtemplate.tex`, no code/data touched (pure reorder +
> one ingested table, per "ingest, don't recompute"). (1) **Reordered** "Stretching/bending
> classification": now establishes `τ_S`/`τ_B` from the IDEAL molecules first (`fig:bondscores` +
> `fig:boxplots` + the numeric-derivation paragraph, moved up), THEN applies to non-ideal molecules
> (`fig:modemixing`), THEN the confusion-matrix result (`fig:confusion`) as the resulting accuracy/
> retention validation (where the `SB` third-class framing now explicitly lands), THEN the non-ideal-
> ground-truth caveat (relocated/reworded, "(below)"→"(above)" self-reference fixed), THEN the
> single-center-topology scope note (unchanged, stays last). (2) Refreshed `fig_benzene_normal.pdf` in
> `images/` from the just-fixed legend-merge regeneration in `data/figures/`. (3) Added a
> `\figplaceholder{...}` (matching the GTOC placeholder idiom) for a new benzene depicted-mode figure
> (3-panel small multiple, modes 12/19/30, green arrows analogous to `fig:watermodes`) right after
> `fig:benzenenormal`'s caption — NOT built, explicitly deferred; caption spec is fully numeric
> (frequencies, `s[V_S]`, `s_AB` values already in the manuscript text) so a future session can render it
> without re-deriving anything. (4) **Removed** `fig:benzene` (the all-36-EMIT-mode scatter, in tension
> with this section's "extreme edge case, not systematic validation" framing) and replaced it with
> `tab:emitselected` (5 rows: EMIT 2/9/34/35/36 — score signature, `s[V_S]`, label, projected contribution
> breakdown), sourced from the same numbers already narrated in the surrounding prose (no new
> computation). `fig_benzene.{pdf,png}` and `src/figures.py::plot_benzene_stress_test` (or equivalent)
> left untouched on disk — only the manuscript embedding was removed. (5) **Moved** `fig:sensitivity` to
> a new standalone SI document, `JCC_man_scoring/JCC_SI_sensitivity.tex` (same house style as the
> existing `JCC_SI_computational_cost.tex`: own `\documentclass`, compiles independently), reproducing the
> figure's unchanged caption as "Figure S1"; the main-text plateau paragraph now cites "Figure~S1 in the
> Supporting Information" instead of embedding the image, and the manuscript's own "Supporting
> Information" subsection lists it explicitly. Because this added a table (`tab:emitselected`) before
> `tab:cost` in reading order, `tab:cost` shifted from Table 6 → Table 7 — updated all 8 hardcoded
> "Table~6" cross-references in `JCC_SI_computational_cost.tex` to "Table~7" (that SI doc predates
> auto-numbering and refers to the main-text table by hardcoded number, not `\ref`). Recompiled all three
> documents from a scratch dir (OneDrive-safe latexmk): main **32 pages before → 32 pages after**, zero
> undefined refs/citations after a second pdflatex pass (bibtex is a latexmk heuristic artifact here —
> the bibliography is a manual `thebibliography` environment, no external `.bib` compile actually
> needed); `JCC_SI_computational_cost.pdf` 7 pages, `JCC_SI_sensitivity.pdf` 2 pages, both clean
> independent builds. All three PDFs copied back to `JCC/JCC_man_scoring/`. No `src/` code changed; no
> test suite impact.
> **2026-07-02 (follow-up correction, same day as the legend-dedup fix above):** author's second visual
> review of `fig_benzene_normal` flagged 3 more issues, all fixed in `plot_benzene_normal_modes`: (1)
> marker shape unified to `marker="o"` for every point (was `CATEGORY_MARKER`'s per-category X/o/s/^/P
> shapes) — color-only category encoding now matches `fig:bondscores`/`fig:boxplots`/`fig:modemixing`'s
> established convention, via `_marker_kwargs()`; (2) all markers now rendered hollow
> (`IDEAL_STYLE["no"]` applied uniformly — no ideal/non-ideal axis within one molecule's own normal
> modes, just cross-figure visual consistency); (3) `τ_S`/`τ_B` threshold-line text labels repositioned
> off the dashed lines (`τ_B` nudged below its line; `τ_S` nudged above its line AND moved to the
> left/upper-left corner, since a right-anchored `τ_S` label collided with the high-frequency S/stretch
> data cluster sitting at s[V_S]≈1.0 — checked by rendering, not assumed); (4) removed the enlarged-
> marker + dashed-leader-line + callout-text-box annotation block for modes 12/19/30 entirely — figure is
> now a plain unannotated scatter (`_WORKED_EXAMPLE_MODES` left defined, unused, no other reference in
> the codebase). Regenerated `fig_benzene_normal.{pdf,png}` only; 53/53 tests still green. **Flag for
> next agent (`lead-author`/manuscript-side):** `fig:benzenenormal`'s LaTeX caption still names the 3
> worked-example callouts ("mode 12... mode 19... mode 30...") that no longer exist in the image —
> caption needs rewording to note that pointing/callout annotations for those modes will be added
> manually later, not describe callouts that aren't there.
> **2026-07-02 (stale-caption fix, `lead-author`) -- DONE, closes the flag above.** Copied the freshly
> regenerated `fig_benzene_normal.pdf` into `JCC/JCC_man_scoring/images/` (overwriting the stale copy).
> Reworded `fig:benzenenormal`'s caption: replaced the false "are called out" claim with "are discussed
> in the text as worked examples ... ; [TODO: manually add enlarged-marker/leader-line callouts for
> these three points in the figure]", matching the doc's existing `[TODO: ...]` idiom (acknowledgments,
> ORCID). Verified the body-text paragraphs narrating modes 12/19/30 read correctly standalone (prose
> only, no dependence on figure annotations) -- left unchanged. Recompiled (OneDrive-safe scratch-dir
> latexmk): 32 pages before -> 32 pages after, zero undefined references/citations.

## Context

The JCC manuscript ("A Unified, Reference-Free Framework for Classifying the 3N modes of molecular
motion") is written, but its tables/figures are mostly `TODO-DATA`. The program currently implements
only **Step 1** (per-mode `Tx..Tz, Rx..Rz, V_Stretch` scoring → CSV). Goal: extend it into the
reference implementation that reproduces every manuscript table and figure. Build order: **core engine
first**, validated on data in hand (water, benzene, gramicidin), then scale out using precomputed
library scores from the Excel file.

### Locked decisions
- **SUPERSEDED 2026-07-03 (see RESUME HERE at top):** ~~Hydride-library logs live off-server; their
  scores are already in `data/vibrational-scoring-functions.xlsx` → ingest precomputed scores +
  reference labels, do not re-score.~~ The author is now supplying real `.log`/`.gjf`/EMIT files for
  the full library. **Current rule:** every score column is recomputed from those files via the real
  engine; the Excel workbook supplies `ref_label`/`ideal` ONLY (ground-truth ties, not scores).
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
- [x] **Pin the projection / mass-weighting convention (eq:emitproj) — DONE 2026-07-01.** LOCKED:
      mass-weight by `sqrt(mass_A)` per atom (all 3 Cartesian components scaled equally), applied
      ONLY inside `src/projection.py` — the unweighted `s[T]/s[R]/s[V_S]` scores are untouched.
      Determined empirically (Gram-matrix orthogonality check on benzene's real normal modes + raw
      EMIT eigenvectors: plain-Cartesian off-diagonals up to 0.80, mass-weighted off-diagonals
      ≤4e-4 — the Eckart-Sayvetz signature), not from the Excel `Eckart`/`Eckart vs score` sheets
      (those turned out to hold an unrelated per-mode Eckart-condition residual check — net
      unweighted ΣΔd and Σ(r×Δd)-like columns near zero — not the projection reference itself; the
      Gram-matrix check was more direct and decisive). See `src/projection.py` module docstring for
      full derivation. Validated to 3.1e-4 max abs deviation against
      `data/results/benzene_EMIT_contributions.csv` (all 36 modes × 9 columns).
- [x] **Verify the Excel column identity — DONE 2026-07-01.** Re-scored H2S and SF2 (hydride-library
      logs+gjf already in repo) in-engine via `run_pipeline(mol, "normal")` and compared each vibrational
      mode's `V_Stretch` to the `data_score` sheet's `sum_|d₂-d₁|²/sum(|d₂-d₁|²)*cosθ` column (col index 14;
      sheet has ONLY the 3N-6 internal vibrational rows, no T/R). **MATCH** to well within 3 dp for all
      6 modes checked (H2S: 0.06317/0.99945/0.99954 vs engine 0.063170/0.999451/0.999539; SF2: 0.05399/
      0.91442/0.9446 vs engine 0.053990/0.914422/0.944605 — agreement to ~5 significant figures, residual
      ~1e-5 consistent with log-file coordinate rounding, not a formula discrepancy). Confirms the Excel
      candidate column IS `s[V_S]` as defined in code; ingest chain (Phase 3) can trust it.
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
- [x] `src/projection.py` (NEW) — **DONE 2026-07-01.** EMIT→normal-mode projection `Θ̃=QᵀΘ`
      (eq:emitproj) using the Phase-0 mass-weighting convention; `build_reference_basis()` +
      `project_emit()`; **emits two data files** — `<mol>_EMIT_contributions.csv` (grouped
      Tx..Rz/VS/VB/VMix fractions, matching the pre-existing ground-truth file's columns) and
      `<mol>_EMIT_projection_full.csv` (per-individual-reference-mode `Θ̃²` detail, not just grouped
      numbers). `main.run_projection_pipeline()` is the headless orchestrator. Reproduced: EMIT
      34/35/36 ≈76.7% translational; EMIT 2 (38.7% Ry) vs EMIT 9 (14.1% Ry) inversion. Max abs
      deviation from the pre-existing hand-derived `benzene_EMIT_contributions.csv`: 3.1e-4 across
      all 36 modes × 9 columns — confirms that file's convention (not a stale/wrong artifact).
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
- [x] `src/excel_ingest.py` (NEW) — **DONE 2026-07-02** (commit `149fc62`). Reads `data_score` +
      `data_mode&bond` (the `characterised modes` sheet's `type` column was cross-checked
      byte-identical to `data_score`'s own `type` column, so no separate read was needed); maps
      eq:vscore col → `s[V_S]`; emits per-bond `s_AB`, frequency, mode-averaged + per-bond Δbond-length,
      the `ideal` tag, and `ref_label` (stretch/bend/translation/rotation) → `data/results/library_scores.csv`.
      Library covers all 11 `tab:ideal` + the full `tab:nonideal` roster except 5 bromides absent from
      the workbook entirely (flagged for lead-author/tex-data-sync, not a forced table edit).
- [x] **Run the classifier over the library** — **DONE 2026-07-02.** `attach_geometry_classification()`
      runs real Algorithm 1 on the 25 geometry-backed molecules, attaching `predicted_label`/
      `predicted_annotation`/T,R scores to internal rows (gated on whole-molecule freq agreement) and
      appending `n_T+n_R` ideal-T/R external rows (independent of that gate) — required for
      `fig:confusion`, now available.
- [x] `src/calibrate.py` (NEW) — **DONE 2026-07-02.** `τ_S`/`τ_B` derived from the ideal-molecule
      library's non-overlapping stretch/bend distributions; `τ_TR` swept 0.05–0.999 (step 0.005) with
      the full curve persisted to `data/results/tau_sensitivity_sweep.csv` (`fig:sensitivity` source);
      frozen to `data/results/thresholds.json` (`τ_TR=0.95, τ_S=0.90368, τ_B=0.17327`);
      `Thresholds.calibrated()` back-fills `Thresholds`. **Plateau defined quantitatively:** longest
      contiguous grid run with zero label-change fraction vs. the previous point AND accuracy at its
      run maximum (see `src/calibrate.py` module docstring for the full statement).
- [x] **Label-level validation (now τ is frozen)** — **DONE 2026-07-02.** Benzene EMIT 34/35 flagged
      `MIXED_EXTERNAL_WITH_VIBRATION`, EMIT 36 the documented `CLEAN_TRANSLATION` blind spot — both
      re-verified under the calibrated (not just provisional) thresholds. Confusion-matrix
      precision/recall via `confusion_matrix_stats()` **pooled over the WHOLE library**: precision 1.0
      all 4 categories; recall 1.0 (T/R), ≥0.95 (bend), 0.717 (stretch, 28.3% into MIXED, 0% into BEND)
      — numeric floor 0.95 applied and reported honestly as not fully met (stretch recall only). Label
      goldens re-pinned in `tests/test_calibrate.py` (this pooled number is unchanged/still correct and
      still tested there). **Sanity check done:** benzene EMIT 1–9 (degenerate 9-fold block) —
      deterministic, exactly 6 external-slot winners total, no over/under-assignment from the
      degeneracy — confirms no block-mechanism is needed (Decision 8), not evidence one exists.
      **Superseded for manuscript/figure-facing purposes 2026-07-02** by the two-tier split directly
      below (`ideal=='yes'` internal rows + all externals = rigorous ground truth; `ideal=='no'`
      internal rows = non-ideal characterization) — see that Changelog entry for why pooling the two
      populations into one accuracy number risks misreading intrinsic non-ideal stretch/bend mixing as
      classifier error. **Formalized into `confusion_matrix_stats()` itself, 2026-07-02 (Task C):** the
      function's return value now ALSO carries `per_category[cat]["recall_ideal"]`/`["recall_nonideal"]`
      (+`n_ref_ideal`/`n_ref_nonideal`), using the identical tier masks `plot_confusion_matrix` applies —
      additive only, pooled keys unchanged, `src/figures.py` untouched. Verified computationally:
      `recall_ideal["stretch"]==1.0`/`recall_ideal["bend"]==1.0` exactly (guaranteed by `tau_S`/`tau_B`'s
      own derivation from this same population's min/max); `recall_nonideal` reproduces the figure's
      95.7%/65.6% numbers. See RESUME HERE for the full write-up; pinned in `tests/test_calibrate.py`.
- [x] Figure `fig:confusion` (two-tier: rigorous ground truth confusion matrix + precision/recall,
      vs. non-ideal label retention/migration-to-mixed) — **REBUILT 2026-07-02** (originally a single
      pooled-library matrix, **DONE 2026-07-02** same day; rebuilt same day after the ideal-vs-non-ideal
      ground-truth-strength distinction surfaced). `src/figures.py::plot_confusion_matrix` now calls
      `confusion_matrix_stats()` TWICE on two ad-hoc `ideal`-column filters of `library_scores.csv`
      (no `calibrate.py` change needed — see Changelog for the full split and numbers).
      `data/figures/fig_confusion.{pdf,png}` (same filename, overwritten; 2x2 layout, was 1x2).
- [x] Figure `fig:bondscores` (bond score vs relative Δbond length, ideal vs non-ideal) —
      **DONE 2026-07-02.** `src/figures.py::plot_bond_scores`; parses `library_scores.csv`'s
      semicolon-joined per-bond `s_AB`/`rel_db` strings via new `_explode_bonds()` helper.
      `data/figures/fig_bondscores.{pdf,png}`.
- [x] Figure `fig:boxplots` (freq, Δ|b|, `s[V_S]`; stretch vs bend) — **DONE 2026-07-02.**
      `src/figures.py::plot_boxplots`; 3-panel box plots, 4 groups per panel (bend/stretch x
      ideal/non-ideal). `data/figures/fig_boxplots.{pdf,png}`.
- [x] Figure `fig:modemixing` (mode score vs averaged Δbond; ideal step vs non-ideal gradient) —
      **DONE 2026-07-02** (core 2-panel content only). `src/figures.py::plot_mode_mixing`. The
      irrep-degeneracy sub-panel content (pending gap, molecule/spec unconfirmed) is explicitly NOT
      built — see the function's docstring and the standing pending-gaps note above.
      `data/figures/fig_modemixing.{pdf,png}`.
- [x] Figure `fig:sensitivity` (label-change fraction & accuracy vs τ; plateau) — **DONE 2026-07-02.**
      `src/figures.py::plot_sensitivity`; tau_TR only (the persisted sweep data covers tau_TR, not
      tau_S/tau_B — no sweep fabricated for those). `data/figures/fig_sensitivity.{pdf,png}`.
- [x] Parity check vs Excel `box plots` / `CM` sheets (and spot-check the other figures' source data) —
      **DONE 2026-07-02.** `box plots` sheet matches `library_scores.csv` exactly; `freq vs score`
      sheet matches `V_Stretch`/`delta_b` distributions with a 2-row-per-group count difference traced
      to SnO2's Excel-side missing frequency (not a code issue); `CM` sheet is an unrelated
      center-of-mass check, not a confusion-matrix reference — see Changelog for the full writeup.

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
- [x] `src/figures.py`: one function per figure + a **single shared style helper** — **DONE 2026-07-02**,
      all 6 manuscript figures now implemented (`fig:benzene` from an earlier session; the 5 Phase-3
      figures this session); `data/figures/` created; each writes **vector PDF** (LaTeX embed)
      **+ ≥300 dpi PNG** preview. Still open for Phase 5 proper: wiring these 6 functions into
      `reproduce.py` (not yet built) as a single orchestrated entry point.
- [ ] **SI Cartesian-geometry export** step (optimized coords for water/benzene/CO₂/gramicidin from logs,
      + library if retrievable) — B7/B14 reproducibility requirement.
- [ ] **Graphical-TOC image** (B4, submission-REQUIRED): 50×50 mm; assemble per the structure-doc concept.
- [ ] Fill the **Gaussian revision/year** `TODO-DATA` (line 463) and correct the inconsistent citation.
- [x] `main.py`: add `--classify`, `--emit-projection`, `--library`, `--figures` subcommands. **DONE
      2026-07-03** — also added `--calibrate` (not in the original bullet text but the same class of
      pipeline-wiring flag, author-requested same session). See RESUME HERE for the full writeup.
- [x] `README.md`: document classification workflow; resolve **JCE-vs-JCC** mismatch (README cites
      *J. Chem. Educ.* "paper I"; this is the JCC unified-framework paper). **DONE 2026-07-03** (rewritten
      TWICE this date — see RESUME HERE: the first pass, earlier the same day, still framed this as a
      two-paper "Paper I (JCE)/Paper II (JCC)" codebase; that framing was explicitly retracted by the
      author later the same day — JCE is being withdrawn before JCC submission, so JCC is the first/only
      paper, not a sequel. The current README has zero JCE/*J. Chem. Educ.* references.)

## Phase 6 — Strengthen for review (RE-TIERED 2026-07-01; see Changelog)
> From expert-reviewer-jcc, re-triaged after Decision 5 (Phase 4/Gramicidin deferred to the companion
> paper) removed the paper's only scale/robustness demonstration — raising the weight the remaining
> validation has to carry. Split into a **recommended-before-submission** tier and an **optional** tier.

### Recommended before submission
- [x] **Flag precision/recall over ALL 36 benzene EMIT modes** (and library externals) against the C2
      projection reference — systematic flag-fidelity, vs the current anecdotal EMIT 2/9/34–36 (A4).
      **DONE 2026-07-02.** New `src/flag_validation.py` (`benzene_emit_flag_confusion()`,
      `library_external_flag_confusion()`, `run_flag_validation_pipeline()`). Ground truth (no canonical
      one exists a priori, same honesty as the task's own framing): per EMIT mode, `M_ext = max(C2_Tx..
      C2_Rz)` (projection fractions); **MIXED iff `0.05 < M_ext < 0.95`, else CLEAN** (`GT_EXT_HI=0.95`
      mirrors `tau_TR` itself; `GT_EXT_LO=0.05` is a small "negligible" floor well above the ~1e-4-1e-3
      projection-orthonormality noise floor). Classifier's predicted flag = `classification ==
      MIXED_EXTERNAL_WITH_VIBRATION` (everything else, including Step-4-only internal labels, is
      "predicted negative" — a mode never assigned an external slot in Step 2 has no chance to be
      flagged at all). **Result, all 36 modes: TP=5, FP=0, FN=12, TN=19 → precision=1.0,
      recall=5/17≈0.294.** Reproduces every named anchor: EMIT 34/35 → TP; EMIT 36 → FN (the documented
      Decision-X blind spot, correctly landing as a miss, not a bug); EMIT 9 → TP vs EMIT 2 → FN despite
      EMIT 2 having MORE genuine Ry character by projection (38.7% vs 14.1%) — direct evidence the
      score/projection ranking inversion produces a wrong Step-2 winner. **Finding: the false-negative
      rate is NOT isolated to EMIT 36's amplitude-invariance blind spot — it is substantially more
      widespread, and mechanistically distinct.** Precision is perfect (the flag never fires on a
      genuinely clean mode) but recall is low (0.294) because Step 2's plain one-to-one
      `linear_sum_assignment` structurally caps the number of EVER-flaggable modes at exactly
      `n_T+n_R=6` out of 36 — the other 30 modes (11 of which have genuine 7-39% external character by
      projection: EMIT 1,2,5,7,8,10,11,12,13,14,18) fall straight to Step 4 and can never be flagged
      regardless of their true external content, because they are not the single best-scoring Hungarian
      assignee for any slot. This is a second, independent false-negative mechanism (structural
      assignment-competition loss) alongside the previously-documented amplitude-invariance blind spot
      (EMIT 36 specifically) — both are honest limitations of a reference-free, one-to-one-assignment
      design, not bugs, but the combined effect is more widespread than the single EMIT-36 anecdote
      suggested. **Library externals (the parenthetical):** all 146 geometry-backed real normal-mode T/R
      reference rows (25 molecules) verified — not assumed — to have FP=0 (ground truth trivially CLEAN
      for all, exact Eckart-Sayvetz completeness); confirms the calibration sweep's 100%-accuracy finding
      from a direct flag-confusion angle. `data/results/benzene_EMIT_flag_confusion.csv` (36-row detail
      table) added; `tests/test_flag_validation.py` (6 new tests, pins TP/FP/FN/TN and the EMIT 2/9/34-36
      anchors) — 40/40 tests green. score-validator dispatched to independently re-run and confirm all
      counts — see Changelog for verdict.
- [ ] **Validate the mixed-SB bucket** by irrep-degeneracy + CoM arguments; report the fraction of
      lit-labeled modes landing in "mixed" (defends against the bending = low-stretch circularity
      concern). **Promoted:** this is the direct evidentiary backbone for the manuscript's residual-risk
      discipline that mixed-SB/flagged-external buckets be "validated by characterization and
      consistency, never by an accuracy claim" — without this item that discipline is asserted but not
      discharged. **Partially discharged for benzene, 2026-07-02:** `src/benzene_validation.py`'s
      `benzene_mixed_bond_diagnostic()` reports the fraction (5/23 ≈ 21.7% of benzene's lit-labeled bend
      modes land in MIXED_STRETCH_BEND, zero in the wrong clean category) AND supplies the
      irrep-degeneracy argument as a computed fact, not an eyeballed one: modes 13/14 and 23/24 are
      detected as near-degenerate pairs (freq splitting <0.03 cm⁻¹) whose 6-bond C-C `s_AB` patterns are
      strongly anti-correlated (r=-0.969/-0.997) — the D6h E-type-pair complementary signature. **Still
      open:** the CoM-softening argument (this item's other half; see the OPTIONAL CoM-conservation item
      below, not yet built) and extending the lit-labeled-fraction report beyond benzene to the whole
      library (the pooled/non-ideal `mixed_fraction` already exists in `confusion_matrix_stats()`, Phase
      3, but has not been explicitly written up as this item's deliverable).
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
- [x] Figures match Excel `box plots`/`CM` sheets (+ spot-checks for the rest) — **DONE 2026-07-02**,
      see Phase-3 changelog entry for the full parity writeup (exact match on `box plots`; explained
      2-row-per-group gap vs `freq vs score`; `CM` sheet is unrelated to fig:confusion).
- [ ] ~~Gramicidin run completes with timing + `s_AB`; numbers inserted where `.tex` has `TODO-DATA`.~~
      **REMOVED 2026-07-01** — Gramicidin deferred to the companion paper (see Phase 4 above); this
      manuscript carries no gramicidin verification target.
- [ ] `py reproduce.py` regenerates all CSVs + figures with no manual steps; `pytest` green.

## Open items to confirm during execution
- [x] (Phase 0) eq:emitproj mass-weighting convention — pinned pre-Phase-2, **DONE 2026-07-01**
      (`sqrt(mass_A)` per-atom weighting, `src/projection.py`-local only; see Phase 0/2 checklist).
- [x] (Phase 0) `data_score` column == `s[V_S]` — **verified 2026-07-01** by in-engine re-score (H2S, SF2; see Phase 0 checklist entry above / Changelog).
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
- **2026-07-05 — New `fig:benzeneconfusion` figure + `confusion_matrix_stats()` fixed for the literal
  "SB" reference label.** Follow-up to commit `69d549e` (lead-engineer), which gave benzene modes 21/22
  a genuine literature 3-class ground truth (`ref_label=="SB"`, an E1u degenerate pair at 1532.85 cm⁻¹),
  built from `src/benzene_validation.py`'s new `benzene_internal_confusion_matrix`/`_sb_vs_stretch_bond_
  diagnostic` (Task A'/B', already landed in the same commit — no scoring/classifier code touched here).
  1. **`plot_benzene_internal_confusion` (`fig:benzeneconfusion`, `src/figures.py`), an 8th figure:**
     benzene's own 3x3 (ref bend/stretch/SB x predicted bend/stretch/mixed) internal-mode confusion
     matrix + per-category precision/recall, reusing the shared `_confusion_heatmap`/`CATEGORY_COLOR`/
     `CATEGORY_LABEL` machinery via a new `REF_LABEL_TO_CATEGORY["SB"] = "mixed"` entry (routes the
     literal "SB" reference row to the same color/label already used for the classifier's own
     MIXED_STRETCH_BEND bucket). Wired into `regenerate_all()`. Numbers match the raw CSVs exactly:
     bend 16/0/2, stretch 0/7/3, SB 2/0/0; precision B=0.889/S=1.000/mixed=0.000, recall
     B=0.889/S=0.700/SB=0.000 (both blind-spot modes 21/22 predicted a clean "B", never "mixed" — the
     tau_B two-gate purity blind spot already on record). Per-bond evidence for the 21/22-vs-23/24
     contrast stays in `benzene_sb_vs_stretch_bond_diagnostic.csv` for a companion LaTeX table, not
     re-plotted inside this heatmap.
  2. **`confusion_matrix_stats()` bug fix (`src/calibrate.py`):** the function's 4-category
     (translation/rotation/stretch/bend) precision/recall accounting had no way to hold a foreign 3rd
     reference label — the 2 "SB" rows (predicted bucket "bend" both times) were silently inflating
     "bend"'s `n_pred` (precision's denominator) without ever counting as a true positive, corrupting
     bend precision from 1.000 to 0.99267/0.9910 (whole-pooled-library / non-ideal-tier calls
     respectively) for a reason unrelated to classifier error. Diagnosed both `fig:confusion` tiers
     first, confirming the fix's actual scope: **rigorous tier (n=237) was never affected at all**
     (benzene is `ideal=='no'` throughout, excluded from that tier by construction) — still exactly
     1.000/1.000 for all 4 categories. **Non-ideal tier (n=420 after excluding SB, was 422 raw)**: bend/
     stretch *retention* (recall) was actually never corrupted either (0.9693/0.6458 before and after —
     `ref_mask` only ever matches the literal string "bend"/"stretch", never "SB"); only *precision* was
     wrong, and only inside `confusion_matrix_stats()`'s own return value (fig:confusion's panel (d)
     never plots non-ideal-tier precision, only recall/migration, so the rendered figure was unaffected
     either way). **Fix:** `confusion_matrix_stats()` now drops any row whose `ref_label` is outside the
     4 recognized categories before building the crosstab/per-category stats (one-line guard + comment;
     general, not benzene-specific — any other future unrecognized label would hit the same silent
     miscount). Restores bend precision to exactly 1.000 in both the pooled-library and non-ideal-tier
     calls; rigorous tier, recall, mixed_fraction, floor_met all unchanged. Updated the one pinned test
     that had (correctly, at the time) documented the 0.99267 number as "expected, not a regression"
     (`test_confusion_matrix_precision_perfect_recall_explained_by_mixed_bucket`,
     `tests/test_calibrate.py`) to assert exactly 1.000 for all 4 categories, with a comment explaining
     the reconsidered reasoning (a literal literature "SB" row answers a different question than the
     nominal-stretch/bend migration-under-mass-effects question this statistic is for; see
     `src/benzene_validation.py::benzene_internal_confusion_matrix` for that 3-class question's own
     dedicated table). 77/77 tests green after the fix (was 77/77 green before, including the
     since-corrected pinned number).
- **2026-07-02 — Visual-consistency pass across 5 figures (author feedback: "not high-quality,"
  inconsistent labels, legend/tick-label overlap, redundant marker-shape encoding).** `src/figures.py`
  only; no scoring/data changes. Four fixes, each visually verified by reading the regenerated PNG (not
  just code inspection):
  1. **Redundant shape encoding removed in `plot_bond_scores` (fig:bondscores) and `plot_mode_mixing`
     (fig:modemixing).** Both used to vary marker SHAPE (square=stretching, circle=bending) on top of
     COLOR for the same stretch/bend distinction — over-encoding. `_marker_kwargs()` gained an optional
     `marker=` override (default `None` → falls back to `CATEGORY_MARKER[category]`, so every other
     caller — there are none besides these two — is unaffected); both functions now pass `marker="o"`
     so every point is a circle regardless of stretch/bend, leaving color (stretch/bend) and fill
     (`IDEAL_STYLE`, filled=ideal/hollow=non-ideal) as the only two encodings. Legends rebuilt to match
     (circle-only proxy handles). `CATEGORY_MARKER` itself is untouched — `fig:confusion`/
     `fig_benzene_normal` still legitimately vary shape across >2 mutually-exclusive buckets.
  2. **Label vocabulary standardized in `plot_boxplots` (fig:boxplots).** Was the only figure using
     short forms ("bend"/"stretch"); now uses the long forms ("bending"/"stretching") matching
     `CATEGORY_LABEL` and every other figure's prose.
  3. **Legend/tick-label overlap fixed in two figures.**
     - `plot_confusion_matrix` panel (b): the `loc="lower left"` precision/recall legend sat directly
       above the rotated x-tick labels ("clean translation" etc.), visually crowding them. Moved to
       `loc="upper left"`; also dropped the 4 identical per-bar "1.000" value labels (precision=recall=
       1.000 for every category — already stated once in the panel title "(all = 1.000)", so they were
       pure redundant clutter *and* the exact space the legend needed) — fixes both problems at once.
     - `plot_boxplots`: the 4 two-line tick labels ("bend\n(ideal)" etc.) collided at every panel width
       tried (uniform spacing, paired spacing — the "(non-ideal)" line is simply too wide for a
       3-panels-per-7.4in figure at any reasonable font size/spacing, since JCC's `\includegraphics[
       width=0.95\columnwidth]` rescales the whole PDF to a fixed ~6.2in target regardless of the
       matplotlib figsize chosen, so growing figsize to "solve" the collision would have silently
       shrunk the final-print font below the 8pt-readability floor). Redesigned as a two-level tick
       scheme instead: short primary labels ("bending"/"stretching" only) plus a single shared "ideal"/
       "non-ideal" bracket+label spanning each pair (drawn once per pair via `ax.get_xaxis_transform()`
       at a fixed axes-fraction y-offset, so it doesn't depend on each panel's y-data range) — appears
       once per pair instead of once per box, so it never needs to be as wide as 4 repeated suffixes.
  4. **Benzene highlight colors de-clashed in `plot_benzene_stress_test` (fig:benzene).**
     `COLORS["highlight_r"]` (`#E69F00` orange, too close to `stretching`'s `#D55E00` vermillion) and
     `COLORS["highlight_t"]` (`#0072B2`, an exact duplicate of `bending`'s blue, despite the callout
     being about translation) both clashed with the `CATEGORY_COLOR` vocabulary. Replaced with the two
     remaining unclaimed hues in the extended Okabe-Ito palette — `highlight_r` → `#F0E442` yellow,
     `highlight_t` → `#000000` black — genuinely distinct from all 5 category colors and from each
     other; every other Okabe-Ito hue is already claimed by a category or a fig:sensitivity color. Since
     a bare yellow line has poor contrast on white, added a black `path_effects` halo to the EMIT-2/9
     connecting arrow (keeps the arrow legible without changing its color) and moved the "EMIT 9" label
     off to the upper-left of its marker (its old lower-right offset now visually collided with the
     thicker haloed arrow).
  Regenerated `fig_bondscores`, `fig_modemixing`, `fig_boxplots`, `fig_confusion`, `fig_benzene`
  ({pdf,png}); copied the 5 PDFs into `JCC/JCC_man_scoring/images/`, overwriting the stale copies.
  `fig_benzene_normal`/`fig_sensitivity` intentionally untouched (audit found no issues there) — spot-
  checked their regenerated PNGs anyway to confirm the shared `COLORS`/`_style()` edits didn't leak into
  them; pixel content is unaffected (only the two colors nobody else references and boxplots-local
  layout constants changed). `py -m pytest tests/` stayed 53/53 green throughout.
- **2026-07-02 — Simplified `fig:benzene` to a single panel; dropped the redundant normal-mode panel.**
  `plot_benzene_stress_test` (`src/figures.py`) used to render 2 panels: (a) `s[V_S]` vs. frequency for
  benzene's 36 real normal modes, (b) EMIT score vs. projected normal-mode contribution (EMIT 2/9
  R_y-inversion + EMIT 34-36 flagged-external callouts). Panel (a) became fully redundant once
  `plot_benzene_normal_modes`/`fig_benzene_normal` was built (same data, better: named worked-example
  callouts for modes 12/19/30) — confirmed the manuscript's "stress test on benzene EMIT modes" prose
  never referenced panel (a)/"36 normal modes" at all, and having normal-mode data inside a figure
  captioned around the EMIT stress test was conceptually confusing regardless of redundancy. Removed
  panel (a) entirely (the `normal_csv` param, the `ax_a` subplot, the `normal.iterrows()` loop, the
  tau_S/tau_B reference lines and their panel-(a)-specific text); `fig:benzene` is now one
  `fig, ax_b = plt.subplots(figsize=(3.8, 3.6))` panel with the former panel (b) content unchanged, no
  "(a)"/"(b)" title prefixes. Dropped the `panel_a_*`/`panel_b_*`-prefixed summary-dict keys (now
  unprefixed: `n_background_points`, `emit2_Ry_score`, `emit34_36_scores`, etc.) — no test or other
  caller referenced the old keys or the `normal_csv` parameter, so nothing else needed updating.
  `.tex` caption for `fig:benzene` checked and needs no change (already EMIT-only content). Regenerated
  `data/figures/fig_benzene.{pdf,png}` and copied the PDF to
  `JCC/JCC_man_scoring/images/fig_benzene.pdf`. 53/53 tests green throughout.
  `plot_benzene_normal_modes`/`fig_benzene_normal` untouched.
- **2026-07-02 — Re-rendered all 6 manuscript figures after the label-vocabulary rename
  (presentation-only regen, not a new figure/analysis).** Ran `py -m src.figures`; all 6 regenerated
  without error. Visually inspected every PNG: legend/annotation text is unchanged from before the
  rename (`CATEGORY_LABEL` string VALUES like `"clean translation"`/`"bending"` never changed, only
  the dict keys did — `git diff --stat` on the 6 PDFs shows 0 insertions/deletions and identical byte
  counts, confirming the rendered content is pixel-identical, only embedded PDF metadata differs).
  No stale/raw label strings found. `py -m pytest tests/` stayed 51/51 green throughout.
- **2026-07-02 — Label-vocabulary rename: classifier output labels are now short and axis-specific
  (author-directed design decision, not a bug fix). Recovered from a session-limit interruption
  mid-rename; nothing was lost.** Old scheme: `CLEAN_TRANSLATION`/`CLEAN_ROTATION` (generic, axis-blind)/
  `MIXED_EXTERNAL_WITH_VIBRATION`/`STRETCHING`/`BENDING`/`MIXED_STRETCH_BEND`. New scheme: clean external
  → the specific Step-2 slot the mode won (`"Tx"`,`"Ty"`,`"Tz"`,`"Rx"`,`"Ry"`,`"Rz"`); mixed external →
  same slot + trailing `"*"` (`"Tx*"` etc.); internal vibration → `"S"`/`"B"`/`"SB"` (the Python constant
  NAMES `STRETCHING`/`BENDING`/`MIXED_STRETCH_BEND` in `src/classifier.py` are unchanged, only their
  string VALUES). `src/classifier.py` gained reusable predicates (`is_external_label`, `external_axis`,
  `is_clean_external`, `is_mixed_external`, `is_translation`, `is_rotation`, `classification_bucket`) so
  downstream code checks label PATTERNS instead of exact-matching a small fixed set of constants.
  Mixed-external annotations dropped the now-redundant `dominant_external=<slot>` prefix (the slot is
  already in the classification string itself) — annotation is now just `"vibration=<S|B|SB>"`.
  **What happened this session:** a `lead-engineer` dispatch did the core rename correctly in
  `classifier.py`/`calibrate.py`/`projection.py`/`benzene_validation.py`/`flag_validation.py` (including
  proactively fixing a latent desync risk in `projection.py` by importing the vibration-label constants
  instead of hardcoding their string values) but was cut off by a session-limit hit before finishing
  propagation to `src/figures.py`, 5 test files, and regenerating the committed CSVs — left the repo
  NOT EVEN IMPORTABLE (`ImportError: cannot import name 'CLEAN_TRANSLATION'`) and nothing committed.
  The main session verified this (ran `git diff`/`git status`, confirmed via `py -m pytest` that
  collection itself failed) and finished the work directly rather than re-dispatching, to conserve
  session budget: reworked `src/figures.py`'s `CATEGORY_COLOR`/`CATEGORY_MARKER`/`CATEGORY_LABEL` to key
  on the 6 semantic buckets `classification_bucket()` returns (`translation`/`rotation`/`stretch`/`bend`/
  `mixed`/`mixed_external`) instead of the old label strings; updated every assertion in
  `tests/test_classifier.py`, `test_calibrate.py`, `test_excel_ingest.py`, `test_flag_validation.py`,
  `test_benzene_validation.py` to the new scheme; then re-ran the full regeneration chain in dependency
  order (`main.run_classify_pipeline`/`run_projection_pipeline` for water/benzene/CO2 →
  `src/excel_ingest.run_ingest_pipeline()` → `src/calibrate.run_calibration_pipeline()` →
  `src/benzene_validation.run_benzene_normal_validation()`/`run_benzene_bond_diagnostic()` →
  `src/flag_validation.run_flag_validation_pipeline()`) to refresh every committed CSV that stores these
  labels as data, not just code.
  **A genuine surprise along the way, not a bug:** the Excel ingest step took roughly **630 seconds
  just to parse the workbook's two sheets** via `openpyxl` (its `Data Validation extension is not
  supported` warning correlates with a known slow path) — nearly 15x the "on the order of a minute"
  estimate in `tests/test_excel_ingest.py`'s/`test_calibrate.py`'s own docstrings. This first looked
  like a hung process (near-zero CPU on the `py` launcher for 10+ minutes) and was killed once
  prematurely before the mistake was caught: `py.exe`(the launcher) shows ~0 CPU because the REAL work
  happens in a child `python3.13.exe` process, which was climbing steadily the whole time — always check
  the child process tree, not just the top-level launcher PID, before concluding a job is stuck. Those
  test docstrings' timing estimate should be corrected in a future pass (flagged here, not fixed this
  session — out of scope for a documentation-only nit during an active recovery).
  **Not yet done this session (queued next, in order):** `figure-builder` needs to re-render all 6
  figures (their legends currently bake in the OLD label text from before this rename) and `lead-author`
  needs to update the `.tex`'s literal `\textsc{clean_translation}`-style references + the benzene table
  to the new short symbols. Historical Changelog entries below this one, and Phase 2's original
  checklist prose, still use the OLD label names as an accurate record of what was true when they were
  written — not updated for this rename, matching this file's own established precedent from the
  earlier τ-renaming (which WAS propagated throughout; this rename's historical-prose mentions were
  judged lower-value to rewrite given the volume of scattered occurrences versus the session budget
  already spent recovering from the interruption above).
- **2026-07-02 — Benzene normal-modes-vs-reference validation formalized as the manuscript's PRIMARY
  classification-vs-reference result (new `src/benzene_validation.py`) + `confusion_matrix_stats()`
  ideal/non-ideal recall split (Task C) + `excel_ingest.py` mismatch-gate robustness fix (Task D).**
  **Manuscript positioning (per `JCC/Scoring_Manuscript_Plan_2026-07-01.pdf`, one directory above this
  repo, which mandates the Results & Discussion order):** "Benzene normal modes -> low-frequency
  stretching" is now the paper's PRIMARY classification-vs-reference validation. Benzene EMIT is
  repositioned as an extreme/rare edge-case stress test only (flag behavior on EMIT 34-36, the EMIT
  2-vs-9 `s[R]` inversion) — NOT a systematic accuracy claim. The systematic 36-mode EMIT confusion
  matrix built earlier this session (`src/flag_validation.py`, commit `65f4e92`) is fully built and
  tested (`tests/test_flag_validation.py`, 6 tests) but is EXCLUDED from the manuscript by author
  decision: its "ground truth" is a threshold cut (`0.05 < M_ext < 0.95`) on continuous,
  genuinely-mixed EMIT projection fractions — circular reasoning for an accuracy claim (the very
  continuum being classified is used to manufacture its own ground truth). This decision is recorded
  here, not re-litigated, and `src/flag_validation.py` itself is untouched.
  **Why benzene's real NORMAL modes are non-circular, unlike EMIT:** the literature/group-theory
  vibrational assignment for each of benzene's 30 real normal modes (`ref_label` in
  `data/results/library_scores.csv`, molecule `C6H6`) is an INDEPENDENT external reference (predates
  this code entirely) — comparing the classifier's own `predicted_label` against it is a genuine
  accuracy check.
  **Task A — `benzene_normal_reference_detail`/`_summary`/`run_benzene_normal_validation`:** reproduces
  the session's ad hoc findings exactly, re-derived from `library_scores.csv` with no hardcoded numbers:
  **6/6 external (T/R) modes correct** (`CLEAN_TRANSLATION`/`CLEAN_ROTATION` matching `ref_label`);
  **7/7 literature-stretch modes recalled** (mode_index 12, 25-30, recall=1.0); **18/23 literature-bend
  modes recalled** (recall=0.7826), the other 5 (mode_index 13,14,19,23,24) landing in
  `MIXED_STRETCH_BEND`, never in the wrong clean category; **ZERO crossings into the opposite clean
  category in either direction** (0/7 stretch->bend, 0/23 bend->stretch), computed via an explicit
  `crossed_opposite` boolean column rather than eyeballed. Raises (fail-loud) if any `ref_label` row
  lacks a `predicted_label` (incomplete geometry merge) or if the molecule has no `ref_label` rows at
  all, rather than silently validating partial/stale data.
  **Task B — `benzene_mixed_bond_diagnostic`/`run_benzene_bond_diagnostic`:** for the 5 modes Task A's
  own output identifies as MIXED_STRETCH_BEND (not hardcoded), parses `library_scores.csv`'s
  semicolon-joined per-bond `s_AB` string (already present, no new score computed) and reports, per mode:
  C-C ring-bond total, C-H total, and `cc_fraction_of_V`. Confirms computationally: C-H contributions are
  ~0 (`ch_total` 0.0088-0.0100 across all 6 C-H bonds per mode; `cc_fraction_of_V`>0.95 for all 5 modes —
  essentially all of `V_Stretch` sits on the C-C ring); mode 19 (1319.2678 cm⁻¹) is uniform across all 6
  ring bonds (coefficient of variation <0.02, i.e. genuinely 6-fold-symmetric, not noise). A
  frequency-proximity check (default tolerance 1 cm⁻¹) DETECTS — rather than assumes — that modes 13/14
  (splitting 0.0157 cm⁻¹) and 23/24 (splitting 0.0292 cm⁻¹) are near-degenerate pairs, and computes the
  Pearson correlation of their 6-bond C-C `s_AB` vectors: **-0.969 (13/14) and -0.997 (23/24)** — a
  strong, computed NEGATIVE (complementary/anti-correlated) pattern, the expected signature of a D6h
  doubly-degenerate (E-type) mode pair. This turns the "complementary alternating pattern"/"genuine
  ring-stretch-bend combination mode, not classifier noise" claim from an eyeballed observation into a
  number. Mode 19 correctly pairs with nothing else in the 5 (its nearest same-bucket neighbor is >100
  cm⁻¹ away).
  Outputs: `data/results/benzene_normal_reference_{detail,summary}.csv`,
  `data/results/benzene_mixed_{bond_diagnostic,degenerate_pairs}.csv`.
  `tests/test_benzene_validation.py` (9 new tests): pins every number above, plus the two fail-loud paths
  (incomplete merge; molecule absent) and a defensive check that mutating away all MIXED_STRETCH_BEND
  labels correctly raises rather than silently reporting an empty diagnostic.
  **Task C — `confusion_matrix_stats()` ideal/non-ideal recall split** (`src/calibrate.py`, independently
  recommended earlier this session by both formula-auditor and lead-author): the pooled recall numbers
  (e.g. stretch=0.7174) mix ideal-molecule rigorous ground truth with non-ideal nominal literature labels
  into one statistic with only a prose claim, no computed/asserted split. Added
  `per_category[cat]["recall_ideal"]`/`["recall_nonideal"]` (+`n_ref_ideal`/`n_ref_nonideal`) — purely
  additive (existing `precision`/`recall`/`tp`/`n_ref`/`n_pred`/`mixed_fraction` keys unchanged). Tier
  masks copied verbatim from `src/figures.py::plot_confusion_matrix`'s existing ad hoc split (commit
  `424a666`) for identical semantics: `ideal` tier = every `kind=='external'` row (regardless of its own
  `ideal` tag — Eckart-Sayvetz completeness makes those exact either way) OR any internal row with
  `ideal=='yes'`; `nonideal` tier = internal rows with `ideal=='no'` only. **Verified computationally, not
  just assumed:** `recall_ideal["stretch"]==1.0` and `recall_ideal["bend"]==1.0` EXACTLY — guaranteed by
  construction since `tau_S`/`tau_B` (`derive_stretch_bend_thresholds`) are literally the min/max of this
  same `ideal=='yes'` population, so no ideal-tier row can land on the wrong side of its own defining
  boundary. `recall_nonideal` reproduces the `fig:confusion` changelog's own numbers exactly (bend
  0.9571, stretch 0.6561). `src/figures.py` and `data/figures/*` were explicitly NOT touched (per
  instruction — that figure already has its own correct, independently-computed ad hoc split wired into
  the manuscript; this is an additive formalization of the same split into the stats function's return
  value, not a replacement for the figure's plotting code). New test in `tests/test_calibrate.py`.
  **Task D — `excel_ingest.py` mismatch-gate robustness fix** (formula-auditor finding, low priority):
  `attach_geometry_classification()`'s `if m is None: continue` branch silently skipped a row whose
  "Vib i" name resolves to no engine mode at all, without counting it toward the mismatch gate — meaning
  such a row could in principle cause a partial (half-merged) attachment with no report, contradicting
  the module's documented all-or-nothing guarantee. Not currently known to happen for any of the 25
  geometry-backed molecules (per the module's own docstring), but not proven impossible. Fixed: `m is
  None` now appends a `(mode_index, None, freq)` sentinel to `mismatches`, so it trips the same gate as a
  frequency mismatch and the molecule is excluded-and-reported like any other bad merge. New synthetic
  regression test in `tests/test_excel_ingest.py` (bogus `mode_index=999` on water/H2O, a real
  geometry-backed molecule) confirms the row is reported and NOT half-merged, while the molecule's
  external T/R rows (independent of this gate) still attach normally.
  **Full suite: 51/51 tests green** (`py -m pytest tests/`; was 40/40 before this session — +1
  `test_excel_ingest.py`, +1 `test_calibrate.py`, +9 new `test_benzene_validation.py`).
- **2026-07-02 — `fig:confusion` rebuilt two-tier (rigorous ground truth vs. non-ideal
  characterization) + `fig:benzene` tau_S/tau_B drift fix, per figure-builder consistency audit.**
  **Why:** the library's literature `ref_label` is exact group-theory ground truth only for
  `ideal=='yes'` internal modes (and for every external T/R row, geometry-exact regardless of the
  `ideal` tag) — for `ideal=='no'` molecules the label is a nominal/dominant-character literature
  assignment, since genuine intrinsic stretch/bend mixing in non-ideal molecules is exactly the effect
  this framework is built to detect (B8.3 CoM-softening). The prior single pooled confusion matrix
  (precision 1.0 / recall 1.0,≥0.95,0.717 — see the Phase-3 label-level-validation bullet above) risked
  being read as one rigorous accuracy claim. lead-author had already patched the manuscript prose with
  a two-tier split; this change makes the FIGURE match the prose, replicating the same ad-hoc
  `ideal`-column filter directly in `src/figures.py::plot_confusion_matrix` (calls
  `src.calibrate.confusion_matrix_stats()` TWICE on filtered slices of `library_scores.csv` — no
  `calibrate.py` code change needed, confirmed by formula-auditor).
  - **Rigorous tier** (n=237: 146 external T/R rows [`kind=='external'`, geometry-exact regardless of
    `ideal`] + 91 `ideal=='yes'` internal rows [50 bend, 41 stretch]): precision AND recall = 1.000 in
    all 4 categories (translation/rotation/stretch/bend) — verified by direct re-run, not assumed.
  - **Non-ideal tier** (n=422: `ideal=='no'` internal rows, 233 bend + 189 stretch): bend retains its
    label 95.7% of the time (223/233), stretch 65.6% (124/189); **0% cross to the OPPOSITE clean
    category** in either direction (0/233 bend→stretch, 0/189 stretch→bend) — 100% of the non-retained
    remainder lands in the mixed bucket (10/233 bend, 65/189 stretch). Framed as validation-by-
    characterization (directionality evidence for CoM-softening), NOT an accuracy claim.
  - **New layout** (2x2, `data/figures/fig_confusion.{pdf,png}`, same filename overwritten): (a)
    rigorous confusion matrix [top-left], (b) rigorous per-category precision/recall bars, all bars =
    1.000 [top-right], (c) non-ideal confusion matrix (bend/stretch reference x
    bend/mixed/stretch predicted) [bottom-left], (d) non-ideal stacked retained-vs-migrated-to-mixed
    bars, with a whole-figure footer stating the 0%-opposite-crossing finding explicitly (avoids a
    per-panel annotation that collided with the panel title at these bar heights during layout
    iteration) [bottom-right]. Same shared `CATEGORY_COLOR`/`CATEGORY_LABEL`/colormap conventions as
    the original single-panel version and every other figure.
  - **Consistency audit (Task A) findings:** `COLORS`/`CATEGORY_COLOR`/`CATEGORY_MARKER`/`IDEAL_STYLE`
    dicts are reused unchanged across all 6 figures — no drift found. **One real drift found and
    fixed:** `fig:benzene` panel (a)'s reference dashed lines were hardcoded to the OLD provisional
    `Thresholds()` class defaults (`tau_S=0.9`, `tau_B=0.2`, annotated "≈0.9"/"≈0.2") from before
    Phase-3 calibration existed, while every other figure that shows tau (`fig:boxplots` panel (c),
    `fig:modemixing`, `fig:sensitivity`) already used the frozen calibrated values
    (`tau_S=0.90368`, `tau_B=0.17327`) — and `benzene_normal_classified.csv`'s own category colors in
    that SAME panel were already computed under the calibrated thresholds (`classify_all_modes`
    defaults to `Thresholds.calibrated()`), so the dashed lines no longer matched the marker colors
    they were meant to explain. Fixed: `src/figures.py`'s module-level `TAU_S`/`TAU_B` now read from
    `Thresholds.calibrated()`; panel-(a) annotations show the exact values (`tau_S=0.904`,
    `tau_B=0.173`) instead of the stale rounded approximations. Panel (b) (EMIT 2/9 inversion, EMIT
    34-36 flagged-external highlights, EMIT-36 Decision-X blind-spot framing) was NOT touched.
    `fig:boxplots`/`fig:modemixing` distribution-only captions were checked and need no change (they
    never compute precision/recall, so the ideal-vs-non-ideal ground-truth-strength distinction does
    not apply to their framing).
  - Regenerated all 6 figures via `python -m src.figures`; `pytest tests/` still 40/40 green (no
    `calibrate.py`/`classifier.py` changes were made, so no golden-test churn).
  - **Action item for lead-author** (not done here — out of scope for figure-builder): copy the
    regenerated `fig_confusion.pdf` into `JCC/JCC_man_scoring/images/`; update the `fig:confusion`
    caption to describe the new 2x2 rigorous/non-ideal layout instead of the old pooled 1x2 one; decide
    whether the prose's own two-tier numeric writeup can now be shortened since the split is shown
    directly in the figure.
- **2026-07-02 — Phase 6 item 1: systematic flag precision/recall over ALL 36 benzene EMIT modes +
  library externals.** New `src/flag_validation.py` (`benzene_emit_flag_confusion()`,
  `library_external_flag_confusion()`, `run_flag_validation_pipeline()`), promoting the prior
  anecdotal EMIT 2/9/34-36 spot-check into a full confusion count. **Ground-truth criterion** (no
  canonical one exists a priori — same honesty the task calls for): per EMIT mode, `M_ext = max(C2_Tx,
  C2_Ty, C2_Tz, C2_Rx, C2_Ry, C2_Rz)` (the projection's own fractional-contribution columns,
  `benzene_EMIT_contributions.csv`); **ground truth = MIXED iff `0.05 < M_ext < 0.95`, else CLEAN**.
  `GT_EXT_HI=0.95` reuses `tau_TR` itself (same "dominant" bar); `GT_EXT_LO=0.05` is a small negligible
  floor well clear of the ~1e-4-1e-3 projection-orthonormality noise documented in `projection.py`.
  Verified empirically no benzene EMIT mode's `M_ext` exceeds ~0.767 (the 34/35/36 triad), so the
  CLEAN-via-`M_ext>=0.95` branch never fires for these 36 modes — all 19 ground-truth-CLEAN modes are
  CLEAN via the `M_ext<=0.05` (purely internal) branch, documented rather than hidden. Classifier's
  predicted flag = `classification == MIXED_EXTERNAL_WITH_VIBRATION` (every other label, including
  Step-4-only internal ones, is "not flagged" — a mode never assigned an external slot in Step 2 has no
  chance to be flagged regardless of its true content). **Result over all 36 modes: TP=5, FP=0, FN=12,
  TN=19 → precision=1.0, recall=5/17≈0.294.** Every previously-anecdotal anchor reproduced exactly:
  EMIT 34/35 → TP; EMIT 36 → FN (Decision-X blind spot, correctly landing as a miss); EMIT 9 → TP vs
  EMIT 2 → FN despite EMIT 2 having the LARGER genuine Ry projection fraction (38.7% vs 14.1%) — because
  EMIT 9's raw `|s[Ry]|` SCORE is larger (0.215 vs 0.143), so EMIT 9, not EMIT 2, wins the Ry slot in
  Step 2's Hungarian assignment — direct evidence the documented score/projection ranking inversion
  actively causes a wrong flag outcome, not just a curiosity. **Finding (the honest answer to the
  task's own question): the false-negative rate is wider than the single previously-documented EMIT-36
  blind spot, and the additional cause is structurally distinct.** EMIT 36's blind spot is an
  *amplitude-invariance* problem (its score signature is indistinguishable from pure translation, no
  bending observable — Decision X, unchanged). The newly-quantified SECOND cause is a *Step-2
  assignment-capacity* problem: plain one-to-one `linear_sum_assignment` can only ever flag exactly
  `n_T+n_R=6` of the 36 modes as external candidates at all (by algorithm design, not a bug) — the
  other 30 fall straight to Step 4, and 11 of those (EMIT 1,2,5,7,8,10,11,12,13,14,18) demonstrably
  carry genuine 7-39% external character by projection yet can NEVER be flagged, because they are not
  the single best-scoring Hungarian winner for any slot. Both mechanisms are honest, reference-free-
  design limitations (not defects), but their combined effect (12/17 genuinely-mixed modes missed) is
  more widespread than the original 3-mode anecdote suggested — reported here plainly, per the
  manuscript's own "flag detects, does not quantify" framing. **Library externals** (task's
  parenthetical): all 146 geometry-backed real normal-mode T/R reference rows across 25 molecules
  verified — not assumed — to have FP=0 (ground truth trivially CLEAN for every row, exact
  Eckart-Sayvetz completeness; consistent with, and a more direct restatement of, `src/calibrate.py`'s
  100%-accuracy tau_TR-sweep finding). Output: `data/results/benzene_EMIT_flag_confusion.csv` (36-row
  per-mode detail table: Mode, classifier_label, M_ext, ground_truth, predicted_positive,
  ground_truth_positive, cell). `tests/test_flag_validation.py` (6 new tests, pins TP/FP/FN/TN plus the
  EMIT 2/9/34-36 anchors) → 40/40 tests green. score-validator dispatched to independently re-run and
  confirm every count; verdict: PASS (see below).
- **2026-07-02 — Phase 3's 5 remaining figures built (`fig:confusion`, `fig:bondscores`,
  `fig:boxplots`, `fig:modemixing`, `fig:sensitivity`) + Excel `box plots`/`CM`-sheet parity
  spot-check.** All 5 added to `src/figures.py` (`plot_confusion_matrix`, `plot_bond_scores`,
  `plot_boxplots`, `plot_mode_mixing`, `plot_sensitivity`), matching the `.tex` caption text
  (`JCC_temp_LaTeXtemplate.tex` lines ~727-790) and reading only already-computed
  `data/results/library_scores.csv` / `tau_sensitivity_sweep.csv` / `thresholds.json` /
  `src.calibrate.confusion_matrix_stats()` — no scores recomputed, no numbers invented. Each writes a
  vector PDF + ≥300 dpi PNG to `data/figures/` and returns a summary dict (n points, axis ranges,
  a `shared_categories` line) for sanity-checking without opening the file; all 6 (`fig:benzene` +
  these 5) were rendered and visually inspected this session, none showed NaN-only or empty axes.
  **`fig:confusion`:** independently reran `confusion_matrix_stats(library_scores.csv,
  Thresholds.calibrated())` and got precision 1.0 (all 4 categories), recall 1.0 (translation/
  rotation), 0.9647 (bend), 0.7174 (stretch, 28.3% of stretch → MIXED, 0% → BEND) — matches the task's
  stated targets exactly. Heatmap (reference label x predicted bucket) + grouped precision/recall bars,
  both colored via the shared `CATEGORY_COLOR` mapping (fixed an initial x-tick-label collision in the
  precision/recall panel by rotating labels 20°). **`fig:bondscores`:** parsed `library_scores.csv`'s
  semicolon-joined per-bond `s_AB`/`rel_db` strings via a new `_explode_bonds()` helper (2755 bonds
  across 69 molecules); plotted `s_AB` vs. `|rel_db|` (absolute relative bond-length change), ideal
  (filled) vs. non-ideal (hollow), stretch (vermillion square) vs. bend (blue circle) — reproduces the
  qualitative shape of the group's earlier undergraduate report's Figure 1 (bending bond scores
  ≤~0.15, ideal-stretch bond scores rising with `|Δb|/|b|`); no fitted/forced quadratic guide curve was
  overlaid (the per-bond score is normalized by a per-MODE bond-count-dependent denominator, so a
  literal `y=x²` line would not actually be a correct reference for multi-bond modes — the manuscript's
  qualitative "quadratic" claim is left to the caption text and the data's own visible curvature, not
  fabricated as an overlay). **`fig:boxplots`:** 3 panels (frequency / mode-averaged `|Δb|/|b|` /
  `s[V_S]`) x 4 groups (bend/stretch x ideal/non-ideal), `τ_S`/`τ_B` reference lines on panel (c);
  513 internal library modes with a literature stretch/bend label (50/41/233/189 per group).
  **`fig:modemixing`:** 2 panels, ideal (filled, clean step at the `τ_S`/`τ_B` gap) vs. non-ideal
  (hollow, graded transition) — visually reproduces the step-vs-gradient contrast described in the
  `.tex` prose and the group's earlier report's Figure 3a/b. The irrep-degeneracy sub-panel content
  (pending gap flagged in an earlier session) was deliberately NOT built — molecule/panel-form still
  unconfirmed with lead-author/tex-data-sync; the function's docstring states this explicitly rather
  than inventing a panel. **`fig:sensitivity`:** single panel, twin y-axes (label-change fraction /
  accuracy) vs. `τ_TR`, shaded plateau band `[0.34, 0.995]` from `thresholds.json`, dashed frozen-value
  line at `τ_TR=0.95`; only `τ_TR` is rendered (the persisted sweep data covers `τ_TR` only — `τ_S`/`τ_B`
  are derived analytically from the ideal-molecule gap in `src/calibrate.py`, not grid-swept — so no
  `τ_S`/`τ_B` sensitivity curve was fabricated to fill out the `.tex` prose's broader "sweeping
  `τ_TR`, `τ_S`, `τ_B`" sentence).
  **Cross-figure style/consistency pass (coordinator directive, same session):** read the group's
  earlier undergraduate research report (`Undergr Res Pj I/H-02-598 Report.pdf`, rendered via
  `pdftoppm` since the PDF-page-reader tool lacked poppler) as the style/quality floor; its Figures 1-3
  are close prior versions of `fig:bondscores`/`fig:boxplots`/`fig:modemixing` and confirmed the
  filled-vs-hollow ideal/non-ideal convention independently of the coordinator's instruction. Added a
  new centralized `IDEAL_STYLE` dict (`{"yes": filled, "no": hollow}`) plus `REF_LABEL_TO_CATEGORY`/
  `PRED_BUCKET_TO_CATEGORY` maps to `src/figures.py` so every new figure keys its `ref_label`/
  predicted-bucket strings into the SAME `CATEGORY_COLOR`/`CATEGORY_MARKER`/`CATEGORY_LABEL` dict
  `fig:benzene` already defined (gray=clean T/R, blue circle=bending, vermillion square=stretching,
  teal triangle=mixed stretch/bend, purple plus=mixed external+vibration) -- one meaning per color/
  marker across the whole 6-figure set, never redefined per-function. Added 3 new `COLORS` entries
  (`sens_accuracy`, `sens_change`, `plateau_band`, `confusion_cmap`) for the genuinely new visual
  elements (`fig:sensitivity`'s dual curves/plateau shading, `fig:confusion`'s heatmap colormap) with
  inline comments on what each means. `fig:benzene` itself was NOT touched -- no conflicting mapping
  arose (its provisional `τ_S≈0.9`/`τ_B≈0.2` dashed-line constants are a separate, already-documented
  concern from the calibrated `τ_S=0.90368`/`τ_B=0.17327` used in the new figures, not a duplicate
  definition). Every new figure's summary dict carries a `shared_categories` line naming exactly which
  shared encodings it reuses, so a reviewer/lead-author wiring captions can see the cross-figure
  consistency was deliberate, not incidental.
  **Excel parity spot-check (IMPLEMENTATION_PLAN.md Phase-3 item):** the `box plots` sheet's own
  pivot (mode score / averaged `|Δb|/|b|`, split ideal x type) matches `library_scores.csv` EXACTLY --
  count/min/max/mean identical to displayed precision for all 4 groups (ideal-bend n=50, ideal-stretch
  n=41, non-ideal-bend n=233, non-ideal-stretch n=189) -- effectively a full pin, not just a spot-check,
  for `fig:modemixing`'s and (partially) `fig:bondscores`' upstream numbers. The `freq vs score` sheet
  (closest analog for `fig:boxplots`) matches on `V_Stretch`/`delta_b` distributions but its ideal-group
  counts run 2 short per group (48 vs. our 50 bend; 39 vs. our 41 stretch) -- traced to SnO2's 4 internal
  rows having `freq=NaN` in the raw `data_score` sheet itself (a genuine, pre-existing gap in the Excel
  workbook, not a rounding/ingest issue); that sheet's own PivotTable silently drops NaN-freq rows
  during its freq-based grouping, while `excel_ingest.py` reads `data_score` directly and correctly
  keeps SnO2's valid `V_Stretch`/`delta_b_mean` values -- `fig:boxplots` panel (a) naturally drops these
  4 NaN-freq rows via per-group `.dropna()` (matching Excel's own displayed frequency range), panels
  (b)/(c) keep them (matching Excel's own `V_Stretch`/`delta_b` ranges, which also include SnO2).
  Flagged for lead-author (SnO2 missing frequency in the source workbook), not silently patched. The
  `CM` sheet is NOT a confusion-matrix reference -- direct inspection showed it holds an unrelated
  center-of-mass-conservation check tabulated by molecular shape (relevant to the OPTIONAL Phase-6
  CoM-conservation item, not `fig:confusion`) -- no Excel confusion-matrix sheet exists in the workbook
  to parity-check against; `fig:confusion`'s numbers were instead independently re-verified against
  `confusion_matrix_stats()` directly (see above), which is authoritative per this task's own framing.
- **2026-07-02 — Phase 3 core landed (ingest, library classification, τ-calibration, confusion
  matrix); recovered intact from a session that hit its usage limit mid-build.** The prior session's
  `lead-engineer` build agent finished writing and testing all of Phase 3's core pieces
  (`src/excel_ingest.py`, `src/calibrate.py`, `tests/test_excel_ingest.py`, `tests/test_calibrate.py`,
  plus `src/classifier.py`'s `Thresholds.calibrated()` addition) but the session hit its session limit
  before the work was committed — everything was left sitting correct-and-tested but uncommitted in the
  working tree. This session verified the full state first (`py -m pytest tests/` → 34/34 green before
  touching anything), read both new modules end to end to confirm they matched the plan's Phase-3 spec,
  then committed everything as-is with no rework needed (commit `149fc62`, 38/38 tests green after
  staging). Nothing was lost. Summary of what shipped (see updated Phase-3 checklist above for detail):
  `library_scores.csv` (25 geometry-backed molecules + Excel-only rows for the rest);
  `thresholds.json` (`τ_TR=0.95, τ_S=0.90368, τ_B=0.17327`) + `tau_sensitivity_sweep.csv`;
  confusion-matrix stats (precision 1.0 all 4 categories; stretch recall 0.717, fully explained by the
  MIXED bucket, 0% bend-confusion). **Left for a future session:** the 5 remaining Phase-3 figures
  (`fig:confusion`, `fig:bondscores`, `fig:boxplots`, `fig:modemixing`, `fig:sensitivity`) and the Excel
  `box plots`/`CM` sheet parity spot-check — figure-builder has everything it needs
  (`library_scores.csv`, `tau_sensitivity_sweep.csv`, `confusion_matrix_stats()`) to build them without
  re-deriving any numbers.
- **2026-07-01 — Phase-0 Excel column identity verified.** Re-scored two hydride-library molecules
  (H2S, SF2 — chosen for simple 3-vibrational-mode C2v/bent-triatomic structure, logs+gjf already in
  repo) headlessly via `run_pipeline(mol, "normal")` and compared per-mode `V_Stretch` against the
  `data_score` sheet's candidate eq:vscore column (`sum_|d₂-d₁|²/sum(|d₂-d₁|²)*cosθ`, column index 14 —
  confirmed empirically by header inspection, not assumed from the sheet name, per the sibling
  Eckart-sheet caution). Sheet contains only the `3N-6` internal vibrational rows per molecule (no T/R
  externals), identified by `molecule`+`mode`(1-based sequential)+`freq` columns; molecule/mode ordering
  in the sheet matches the engine's Vib-1/2/3 order by ascending frequency for both test molecules, so no
  reordering was needed.
  **Result: MATCH for all 6 modes checked, well within the 3 dp tolerance** (H2S bend/sym-stretch/
  asym-stretch: Excel 0.06317/0.99945/0.99954 vs engine 0.063170/0.999451/0.999539; SF2: Excel
  0.05399/0.91442/0.9446 vs engine 0.053990/0.914422/0.944605). Agreement is to ~5 significant figures
  (residual ≤1e-5), i.e. an order of magnitude tighter than the 3 dp requirement — consistent with
  minor geometry-precision differences between the Gaussian log's printed coordinates and whatever
  precision produced the Excel workbook, not a formula/convention mismatch. **Conclusion: the Excel
  `data_score` candidate column IS `s[V_S]` as implemented in `scoring.py`'s `Vscore`/`_bond_contributions`
  — Phase 3's `excel_ingest.py` can map it directly with no transformation.** Added
  `data/results/H2S_normal_scores.csv` and `data/results/SF2_normal_scores.csv` (new goldens for these
  two library molecules, first time they've been run through the engine).
- **2026-07-01 — Projection convention pinned + `src/projection.py` landed (Phase 0 + Phase 2's
  other half).** **Convention LOCKED:** the C2 normal-mode reference basis used to compare against
  benzene EMIT modes is mass-weighted by `sqrt(mass_A)` per atom (all 3 Cartesian components of a
  given atom scaled by the same factor), applied ONLY inside `src/projection.py` — the unweighted
  `s[T]/s[R]/s[V_S]` scores in `src/scoring.py` are untouched. Determined empirically, not from the
  Excel `Eckart`/`Eckart vs score` sheets (those turned out to hold an unrelated per-mode
  Eckart-condition residual check, not the projection reference): a direct Gram-matrix computation
  on benzene's 30 real normal modes and 36 raw EMIT eigenvectors showed both sets are unit-normalized
  under the plain (unweighted) Cartesian dot product but are only mutually ORTHOGONAL under the
  mass-weighted inner product `<u,v>=Σ_A m_A u_A·v_A` (plain-Cartesian Gram off-diagonals up to 0.80;
  mass-weighted Gram off-diagonals ≤4e-4, i.e. numerical noise, not real non-orthogonality) — the
  standard Eckart-Sayvetz signature (true harmonic normal modes are orthonormal in mass-weighted
  coordinates; Gaussian, and empirically EMIT too, report them back-transformed to unweighted
  Cartesian and renormalized to unit Euclidean length for display).
  Added `src/projection.py`: `mass_weights_from_scorer()`, `build_reference_basis(scorer, final_normal,
  thresholds)` (mass-weights + renormalizes the `normal`-pool candidates — ideal T/R built the same
  geometric way as `ModeScorer.construct_T/construct_R`, plus the real `3N-6` vibrational normal
  modes — into a basis `Q` that is orthonormal to within the input files' print precision; also
  buckets each real vibrational mode into STRETCHING/BENDING/MIXED_STRETCH_BEND by reusing
  `classifier.vib_label()` on that mode's own unweighted `s[V_S]`, so the stretch/bend boundary is
  defined in exactly one place), and `project_emit(ref, final_emit)` (mass-weights + renormalizes the
  raw EMIT eigenvectors the same way, computes `Θ̃=QᵀΘ` — eq:emitproj — and squares for fractional
  contributions, with an internal Parseval sanity assert). `main.py` gained
  `run_projection_pipeline(mol_name, data_dir="data", thresholds=None, write=True)`, writing
  `<mol>_EMIT_contributions.csv` (grouped `C2_Tx..C2_Rz`/`C2_VS`/`C2_VB`/`C2_VMix`, matching the
  pre-existing hand-derived `data/results/benzene_EMIT_contributions.csv`'s columns/semantics) and
  `<mol>_EMIT_projection_full.csv` (per-individual-reference-mode `Θ̃²` detail — one column per ideal
  T/R slot and per real normal mode — the richer "projection-coefficients data file" the roadmap
  calls for, not just the grouped numbers).
  **Validation:** regenerating benzene's contributions against the pre-existing hand-derived CSV
  (backed up before overwrite) gave a **max absolute deviation of 3.1e-4 across all 36 EMIT modes ×
  9 grouped columns** — confirms that file's convention was already this one (not a stale/wrong
  artifact) rather than contradicting it. Reproduces the manuscript's named targets: EMIT 34/35/36
  ≈76.7% translational (`C2_Tx/Ty/Tz≈0.7674`, headline "~77%"); EMIT 2 (38.7% `C2_Ry`) vs EMIT 9
  (14.1% `C2_Ry`) — the non-monotonicity inversion relative to the `|s[Ry]|` SCORE ordering (EMIT 9
  0.215 > EMIT 2 0.143) is reproduced in both directions. Every EMIT mode's 9 grouped fractions sum
  to 1.000±0.0003 (Parseval, consistent with `Q`'s near-orthonormality). Water also runs cleanly
  (sums to exactly 1.0 per mode; its 3 vib modes split 2 stretch + 1 bend, matching tab:water).
  formula-auditor: **PASS** — confirmed the `Θ̃=QᵀΘ` matrix algebra is transpose-correct, verified
  the T/R/vibration orthogonality claims analytically (not just empirically: `T_i·T_j=0` trivially;
  `T_i·R_j∝ΣmA rA=0` from `COM()`; `R_i·R_j∝` off-diagonal inertia-tensor entries `=0` from `MIT()`
  diagonalizing before `construct_R()` runs; vib-mode orthogonality to T/R is the Eckart-Sayvetz
  theorem), and confirmed no conflict with the `.tex`'s "unweighted scores, stated mass-weighted
  projection reference" framing. Flagged two non-blocking items: (1) the implicit "same-log,
  same-rotation" coupling between the independent `normal`/`emit` parses that makes `Q` and `Θ`
  frame-consistent was relying on parsing determinism silently — **fixed same session**: added a
  runtime `np.allclose` assertion in `run_projection_pipeline()` comparing `scorer_n`'s and
  `scorer_e`'s post-MIT atom coordinates, raising with a diagnostic message if they ever diverge;
  (2) `eq:emitproj` in the `.tex` states only `Θ̃=QᵀΘ` without spelling out that the "contribution
  coefficients" are the SQUARED coefficients (Parseval) — the manuscript's own percentage language
  ("77% translational", "39% Ry") is only consistent with squaring, but the equation's surrounding
  text doesn't say so explicitly; **flagged for lead-author**, not fixed here (out of engineering
  scope; suggested added clause: "with fractional contributions given by the squared coefficients
  `Θ̃ᵢ²` (Parseval, since `Q` is orthonormal)"). score-validator: **PASS** — independently reproduced
  the 3.1e-4 max-deviation figure (found 3.12e-4, tied at EMIT 3/EMIT 15 `C2_VB`), confirmed all spot
  values (EMIT 34/35/36 `C2_T*≈0.767`; EMIT 2 `C2_Ry=0.3867`, EMIT 9 `C2_Ry=0.1411`, ordering
  inverted vs. the `|s[Ry]|` scores 0.1427/0.2151 read from `benzene_EMIT_scores.csv`), and confirmed
  17/17 tests green. Added `tests/test_projection.py` (4 tests, pinning the 77%-translational and
  Ry-inversion targets at a 2e-3 tolerance, well above the observed ≤3.2e-4 deviation).
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
