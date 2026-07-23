"""Figure generation for the JCC manuscript ("A Unified, Reference-Free
Framework for Classifying the 3N modes of molecular motion").

One function per labeled figure (``fig:benzene``, ``fig:confusion``, ...) so
Phase-5 ``reproduce.py`` can call each in turn. Each function:

  * reads already-computed results CSVs under ``data/results/`` (never
    recomputes scores -- this module is presentation-only),
  * builds the figure with the shared house style (``_style()``），
  * writes both a vector PDF and a >=300 dpi PNG under ``data/figures/``,
  * returns a small summary dict (paths + a few sanity numbers) so the
    caller can confirm the plot's content without opening the file.

All 6 originally-scoped manuscript figures are implemented:
``plot_benzene_stress_test`` (``fig:benzene``), ``plot_confusion_matrix``
(``fig:confusion``), ``plot_bond_scores`` (``fig:bondscores``),
``plot_boxplots`` (``fig:boxplots``), ``plot_mode_mixing``
(``fig:modemixing``), and ``plot_sensitivity`` (``fig:sensitivity``).

``plot_benzene_normal_modes`` is a 7th, standalone figure (no ``fig:`` label
of its own yet -- pending lead-author's rewrite of the "Benzene normal
modes" section): the descriptive worked-example companion to fig:benzene's
EMIT stress test, showing all 36 real normal modes with the ring-breathing
(mode 12) / mixed S-B (mode 19) / C-H stretch (mode 30) worked examples
called out. It does not modify or replace ``plot_benzene_stress_test``.

``plot_benzene_internal_confusion`` (``fig:benzeneconfusion``, added
2026-07-05) is an 8th figure, enabled by the literature relabeling of
benzene modes 21/22 as a literal 3-class "SB" (mixed) ground truth (commit
69d549e): benzene's own 3x3 (bend/stretch/SB reference x bend/stretch/mixed
predicted) internal-mode confusion matrix, distinct from the whole-hydride-
library ``fig:confusion``.

``plot_benzene_internal_confusion`` was RESTRUCTURED 2026-07-05 (later the
same day, author decision, tone/framing only -- NOT the same circularity
argument as the ``fig:confusion``/``plot_rigorous_tier_check`` split above):
benzene's literature stretch/bend/SB ground truth (Shi 1972) is genuinely
external, non-circular ground truth, so there is no accuracy-inflation
concern here. The issue is instead that a precision/recall bar panel reads
as a formal accuracy-metric claim, which sits awkwardly next to this
section's own repeated, deliberate hedging language for benzene specifically
("calling a mode 'mixed' by eye is a convention, not an exact measurement").
The 3x3 heatmap alone is more honest here: it shows where modes landed
without asserting a formal metric. The precision/recall bar panel (still the
exact same, unchanged numbers) now lives in a 10th figure,
``plot_benzene_confusion_precision_recall`` (``fig_benzene_precision_recall``,
no ``fig:`` label of its own -- SI-bound, not main text).

``plot_confusion_matrix`` (``fig:confusion``) was RESTRUCTURED 2026-07-05
(author decision) from a two-tier 2x2 layout down to a SINGLE-TIER 1x2
layout showing only the non-ideal tier (the genuine, non-circular
validation: thresholds fixed on the ideal population, applied without
retuning to harder non-ideal cases). The removed rigorous-tier panels
(precision/recall = 1.000, which is circular/close-to-definitional by
construction -- see the ``fig:confusion`` section header comment below for
the full argument) now live in ``plot_rigorous_tier_check``, a 9th figure
(no ``fig:`` label of its own -- SI-bound, framed explicitly as a
self-consistency check, not an accuracy claim).

``plot_irrep_coupling`` (added 2026-07-06) is an 11th figure, no ``fig:``
label of its own -- SI-bound: the irrep-degeneracy mixing-mechanism figure
that fills the ``fig:modemixing`` pending gap flagged in
IMPLEMENTATION_PLAN.md (molecule/panel form now confirmed as the AB3
trigonal-planar / AB2 bent series, NOT ethane). Reproduces the group's
earlier project report's Figures 3c/3d: mode score vs. central-atom
displacement amplitude, faceted by whether a mode's irrep has a same-irrep
coupling partner in its point group. **Repointed 2026-07-08** (data_score.csv
retirement): reads ``irrep``/``shape``/``type`` from
``data/characterised_modes.csv``, the ideal/non-ideal filter from
``data/mol_list_method.csv``'s ``mol_type`` column, and ``V_Stretch``/the
new engine-derived ``d_CA`` column from ``data/results/library_scores.csv``
-- ``data/data_score.csv`` is no longer read anywhere in this module.

``plot_confusion_matrix`` (``fig:confusion``) and ``plot_benzene_internal_confusion``
(``fig:benzeneconfusion``) were both simplified **2026-07-09** (author decision):
the empty ``mixed_external``/"T/R*" column is dropped from both (verified 0
on both populations, not assumed), and the thick divider line between the
external and internal blocks is removed from both. ``fig:benzeneconfusion``
additionally collapses its 6 individual external rows/columns (Tx..Rz) down
to one "T/R" row/column -- the per-axis breakdown is redundant there since
``fig:confusion`` already shows it in full for the (harder) non-ideal
population; benzene's own figure only needs to confirm its 6 external modes
land as a block. ``fig:confusion`` itself is UNCHANGED in this respect and
keeps its 6 individual axis rows/columns.

Cross-figure visual consistency (one meaning per color/marker, paper-wide)
--------------------------------------------------------------------------
Every figure that encodes a classification CATEGORY (clean translation/
rotation, stretching, bending, mixed stretch/bend, mixed external+vibration)
reuses the exact same ``CATEGORY_COLOR``/``CATEGORY_MARKER``/``CATEGORY_LABEL``
mapping defined once below -- this dict is shared by ``fig:benzene`` and all
5 new figures, never redefined per-function. The one new visual dimension
introduced by the library figures (``fig:bondscores``, ``fig:boxplots``,
``fig:modemixing``) is IDEAL vs. NON-IDEAL molecule membership, encoded via
``IDEAL_STYLE`` as filled (ideal) vs. hollow/open (non-ideal) markers/boxes
-- also centralized here rather than redefined per plot. A reader who learns
"blue circle = bending, vermillion square = stretching, filled = ideal
molecule" once can carry that reading across the whole figure set.
"""

from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.patheffects as pe

from src.classifier import (
    Thresholds, vib_label, classification_bucket,
    STRETCHING, BENDING, MIXED_STRETCH_BEND,
    is_clean_external, is_mixed_external, external_axis,
)

# --------------------------------------------------------------------------
# Shared house style
# --------------------------------------------------------------------------

# ColorBrewer Set1 red/blue/purple for S/B/SB (2026-07-21: replaced the
# Okabe-Ito vermillion/blue/teal trio at the user's request for a more
# conventional red=stretch/blue=bend/purple=mixed scheme; Set1 was chosen
# because its hues are published together as a set for pairwise qualitative
# distinguishability, rather than picked independently).
#
# "bending" REDEFINED (2026-07-23, author request): Set1's #377EB8 was judged
# "not the real blue colour -- a little bit lighter." Replaced with #1A5DE4,
# the true colorimetric complement of the "stretching" red #E41A1C (HSL hue
# rotated 180 deg, same saturation/lightness as the red) but hue-shifted from
# the literal complement's cyan/turquoise (~180 deg) to a real blue (220 deg,
# author-chosen) since red's exact hue-wheel complement is cyan, not blue.
# Scope (author instruction): this NEW blue is for CLASSIFICATION-ALGORITHM
# mode labels ONLY (CATEGORY_COLOR, below) -- i.e. wherever a color encodes
# what the classifier itself predicted for a mode. Reference/ground-truth
# (literature) classifications keep the OLD Okabe-Ito blue/vermillion pair --
# see "bending_ref"/"stretching_ref" below and REF_CATEGORY_COLOR.
COLORS = {
    "external": "#999999",   # gray      -- CLEAN_TRANSLATION / CLEAN_ROTATION
    "bending": "#1A5DE4",    # blue      -- BENDING (classification-algorithm labels only)
    "stretching": "#E41A1C", # red       -- STRETCHING
    "mixed": "#984EA3",      # purple    -- MIXED_STRETCH_BEND
    # Reference/ground-truth (literature) stretch/bend colors -- the OLD,
    # original Okabe-Ito pair, kept distinct from the classification
    # algorithm's own "bending"/"stretching" hues above (2026-07-23 split).
    # Used by every figure whose color encodes a REFERENCE/literature label
    # rather than the classifier's predicted label (fig:bondscores,
    # fig:boxplots, fig:modemixing, fig:irrep_coupling, the retention/
    # migration SI panel, and the reference-axis half of any confusion
    # matrix) -- see REF_CATEGORY_COLOR below.
    "bending_ref": "#0072B2",    # Okabe-Ito blue
    "stretching_ref": "#D55E00", # Okabe-Ito vermillion
    "mixed_ext": "#CC79A7",  # pink      -- MIXED_EXTERNAL_WITH_VIBRATION
    "background": "#BBBBBB", # light gray-- unhighlighted context points
    # highlight_r/highlight_t (fixed 2026-07-02 consistency pass): the
    # previous values (#E69F00 orange -- too close to stretching's #D55E00
    # vermillion; #0072B2 blue -- an EXACT duplicate of the bending color,
    # despite this callout being about translation) both clashed with the
    # established CATEGORY_COLOR vocabulary. Replaced with the two remaining
    # unclaimed hues in the extended Okabe-Ito colorblind-safe palette
    # (yellow, black) -- every other Okabe-Ito hue is already claimed by a
    # CATEGORY_COLOR or a fig:sensitivity color (see palette note below), so
    # these are genuinely free and mutually distinct from all five category
    # colors (external gray, bending blue, stretching vermillion, mixed teal,
    # mixed_ext purple) and from each other.
    "highlight_r": "#F0E442",# yellow    -- s[R] inversion (EMIT 2 / 9)
    "highlight_t": "#000000",# black     -- flagged external (EMIT 34-36)
    "threshold": "#555555",  # dark gray -- tau reference lines (all figures)
    # -- New, figure-specific-but-shared-key colors (fig:sensitivity / fig:confusion) --
    "sens_accuracy": "#56B4E9",   # sky blue  -- accuracy curve, fig:sensitivity
    "sens_change": "#000000",     # black     -- label-change-fraction curve, fig:sensitivity
    "plateau_band": "#56B4E9",    # sky blue @ low alpha -- tau_TR plateau shading, fig:sensitivity
    # sequential, colorblind-safe -- fig:confusion/fig:benzeneconfusion heatmap.
    # CHANGED 2026-07-23 (author request) from "Blues": that map's high-value
    # end (~#2070b4/#08306b) reads as nearly the same blue as the "bending"
    # category color used for this SAME heatmap's tick labels, risking a cell
    # shaded dark for a high COUNT (any category) being misread as "about
    # bending." Green is not claimed by any category color (bend=blue,
    # stretch=red/vermillion, mixed=purple, external=gray), so it can't collide.
    "confusion_cmap": "Greens",
}

# Reference-label (library ground truth) / predicted-bucket -> shared
# classification-CATEGORY name. Since the 2026-07-02 label rename made clean/
# mixed-external classifier labels axis-specific ("Tx".."Rz", "Tx*".."Rz*"
# -- 12 distinct raw strings, no longer 2-3 fixed constants), the CATEGORY
# vocabulary these dicts translate INTO is now the 6 bucket names
# src.classifier.classification_bucket() already returns ("translation",
# "rotation", "stretch", "bend", "mixed", "mixed_external") -- and ref_label/
# pred_bucket values already ARE exactly those bucket names, so both dicts
# below are now IDENTITY maps. Kept (not deleted) so the ~15 call sites
# elsewhere in this file (`CATEGORY_COLOR[REF_LABEL_TO_CATEGORY[r]]` etc.)
# don't all need touching; for a RAW per-mode classifier label (e.g.
# fig_benzene_normal's `row["label"]`, which can be "Tx", "Tx*", "S", ...),
# convert with `classification_bucket()` first, then use directly as the key.
REF_LABEL_TO_CATEGORY = {
    "translation": "translation",
    "rotation": "rotation",
    "stretch": "stretch",
    "bend": "bend",
    # "SB" (2026-07-05): benzene modes 21/22's genuine literal literature
    # mixed-stretch/bend reference label (see commit 69d549e / tab:benzenemixed
    # / src.benzene_validation.benzene_internal_confusion_matrix) -- routed to
    # the same "mixed" category already used for the classifier's OWN
    # MIXED_STRETCH_BEND predicted bucket, so CATEGORY_COLOR/CATEGORY_LABEL
    # (teal, "SB") and `_confusion_heatmap` can be reused as-is for
    # `plot_benzene_internal_confusion` (fig:benzeneconfusion) without adding
    # a new color/marker/label vocabulary entry. This is the only place "SB"
    # is ever a valid dict key here; the whole-library `fig:confusion`
    # (plot_confusion_matrix) never iterates an unbounded category list --
    # its `cats_r`/`cats_n` are hardcoded to the 4 known categories -- so
    # this addition has no effect there, and `confusion_matrix_stats()`
    # (src/calibrate.py) already excludes ref_label=="SB" rows from its own
    # 4-category accounting entirely (see that function's 2026-07-05 fix).
    "SB": "mixed",
    # Individual-axis identity entries (added for the joint T/R+internal
    # confusion matrices, fig:confusion/fig:benzeneconfusion/the SI rigorous-
    # tier table): a bare clean-external label ("Tx".."Rz") IS its own
    # category here, not collapsed to the coarse "translation"/"rotation"
    # bucket -- see `_axis_aware_pred_category` below, which is what actually
    # produces these keys for a joint table's ref/pred columns.
    "Tx": "Tx", "Ty": "Ty", "Tz": "Tz",
    "Rx": "Rx", "Ry": "Ry", "Rz": "Rz",
    # Collapsed single-class entry (fig:benzeneconfusion only, via
    # `_joint_confusion_table(..., collapse_external=True)`): all 6 axes
    # folded into one "T/R" row/column, since per-axis detail is already
    # shown in fig:confusion's non-ideal matrix.
    "T/R": "T/R",
}
PRED_BUCKET_TO_CATEGORY = dict(REF_LABEL_TO_CATEGORY)
PRED_BUCKET_TO_CATEGORY.update({
    "mixed": "mixed",
    "mixed_external": "mixed_external",
})

# The one new visual dimension needed by the library figures: ideal vs.
# non-ideal molecule membership. Centralized here (not redefined per
# function) -- filled marker/box = ideal, hollow/open = non-ideal, matching
# the convention already used in the group's earlier undergraduate report
# (H-02-598, Figures 1-3: filled markers for ideal molecules, open markers
# for non-ideal).
IDEAL_STYLE = {
    "yes": {"filled": True, "label": "ideal", "alpha": 0.85},
    "no": {"filled": False, "label": "non-ideal", "alpha": 0.85},
}

# Keyed by the 6 CATEGORY bucket names (see REF_LABEL_TO_CATEGORY comment
# above), not by raw classifier label strings -- translation/rotation share
# one color+marker (as they always have; only the legend/tick TEXT in
# CATEGORY_LABEL distinguishes them), matching every individual axis-specific
# clean-external label ("Tx".."Rz") or mixed-external label ("Tx*".."Rz*")
# once routed through classification_bucket().
CATEGORY_COLOR = {
    "translation": COLORS["external"],
    "rotation": COLORS["external"],
    "bend": COLORS["bending"],
    "stretch": COLORS["stretching"],
    "mixed": COLORS["mixed"],
    "mixed_external": COLORS["mixed_ext"],
    # Individual axes share the same "external" gray as the coarse
    # translation/rotation buckets -- one grouping color, distinguished by
    # tick-label text only (see CATEGORY_LABEL below), not by 6 new hues.
    "Tx": COLORS["external"], "Ty": COLORS["external"], "Tz": COLORS["external"],
    "Rx": COLORS["external"], "Ry": COLORS["external"], "Rz": COLORS["external"],
    "T/R": COLORS["external"],
}

# Reference/ground-truth (literature) counterpart of CATEGORY_COLOR (2026-07-
# 23 split, author instruction): identical except "bend"/"stretch" use the
# OLD, original Okabe-Ito blue/vermillion pair ("bending_ref"/"stretching_ref"
# in COLORS) instead of the classification algorithm's own hues. Used
# wherever a color encodes a REFERENCE/literature label -- fig:bondscores,
# fig:boxplots, fig:modemixing, fig:irrep_coupling, the retention/migration SI
# panel, and the "Reference label" axis of any confusion-matrix heatmap.
# Never used for a classifier-PREDICTED label -- those keep CATEGORY_COLOR.
REF_CATEGORY_COLOR = dict(CATEGORY_COLOR)
REF_CATEGORY_COLOR.update({
    "bend": COLORS["bending_ref"],
    "stretch": COLORS["stretching_ref"],
})

CATEGORY_MARKER = {
    "translation": "X",
    "rotation": "X",
    "bend": "o",
    "stretch": "s",
    "mixed": "^",
    "mixed_external": "P",
}

# FIXED 2026-07-02 (IMPLEMENTATION_PLAN.md queued item 2): the internal
# stretch/bend/mixed classifications display the SAME short symbols the
# classifier itself emits (src.classifier.STRETCHING="S"/BENDING="B"/
# MIXED_STRETCH_BEND="SB") and the manuscript prose/tab:benzenemixed use
# (`\texttt{S}`/`\texttt{B}`/`\texttt{SB}`), instead of spelling them out.
# REVERTED 2026-07-02 (author visual review, same day): a same-day attempt to
# gloss each symbol in-figure ("S (stretching)" etc.) was reviewed and judged
# too long for a standalone legend/tick label -- reverted back to the bare
# short symbol. The gloss is instead added ONCE, in the LaTeX caption text
# (lead-author's job, not baked into the raster), rather than repeated on
# every legend/tick occurrence across 5 figures. translation/rotation are
# untouched: the manuscript has no single-letter T/R bucket symbol to match
# (its actual short labels are axis-specific, "Tx".."Rz"), so "clean
# translation"/"clean rotation" stay as the pre-existing long form.
# `mixed_external` similarly has no short symbol in the manuscript prose and
# is left unchanged.
CATEGORY_LABEL = {
    "translation": "clean translation",
    "rotation": "clean rotation",
    "bend": "B",
    "stretch": "S",
    "mixed": "SB",
    "mixed_external": "mixed external+vibration",
    # Individual axes: bare axis symbol, already <=2 chars, no gloss needed
    # (matches the classifier's own "Tx".."Rz" clean-external label strings).
    "Tx": "Tx", "Ty": "Ty", "Tz": "Tz",
    "Rx": "Rx", "Ry": "Ry", "Rz": "Rz",
    "T/R": "T/R",
}

# Calibrated classifier thresholds (src/calibrate.py's frozen
# data/results/thresholds.json, read via Thresholds.calibrated()), quoted
# here only for the reference dashed lines in fig_benzene_normal (fig:benzene
# itself carries no tau_S/tau_B lines since its 2026-07-02 single-panel
# simplification dropped the normal-mode s[V_S]-vs-frequency content). These
# used to be the provisional Thresholds() class defaults (0.9/0.2) hardcoded
# before Phase-3 calibration existed; fixed 2026-07-02 (consistency audit) to
# track the SAME calibrated values used everywhere else (fig:boxplots panel
# (c), fig:modemixing, fig:sensitivity, and the classification that produced
# benzene_normal_classified.csv's own category labels) -- not refit, not
# invented, just no longer stale.
_TAU = Thresholds.calibrated()
TAU_S = _TAU.tau_S
TAU_B = _TAU.tau_B

# Named font-size overrides (2026-07-23): every plot element inherits
# _style()'s rcParams defaults unless listed here. These two constants are
# the ONLY permitted exceptions -- used solely where the default size,
# verified by rendering, collides with data/other text in a specific dense
# panel (a multi-entry legend crowding a plot area, or an in-plot value
# label that no longer fits inside a thin bar segment). Do not add a bare
# `fontsize=<number>` anywhere else; if a new collision turns up, reuse one
# of these two rather than inventing another magic number.
LEGEND_FONTSIZE = 9
ANNOTATION_FONTSIZE = 9


def _style():
    """Apply the shared print-ready matplotlib style. Call once per figure."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 14,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "lines.linewidth": 1.0,
        "lines.markersize": 5,
        "figure.dpi": 300,
        "savefig.dpi": 600,
        "pdf.fonttype": 42,   # embed as real fonts, not Type-3 bitmaps
        "ps.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def _savefig(fig, out_dir, label):
    """Save ``fig`` as ``<out_dir>/<label>.pdf`` and ``.png`` (>=300 dpi).

    Returns (pdf_path, png_path). Raises if either file ends up missing or
    empty, so a silent matplotlib failure can't masquerade as success.
    """
    os.makedirs(out_dir, exist_ok=True)
    pdf_path = os.path.join(out_dir, f"{label}.pdf")
    png_path = os.path.join(out_dir, f"{label}.png")
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, bbox_inches="tight", dpi=600)
    for p in (pdf_path, png_path):
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            raise RuntimeError(f"figure save failed or produced an empty file: {p}")
    return pdf_path, png_path


def _marker_kwargs(category, ideal_flag=None, marker=None, color_map=None):
    """Shared per-point marker styling for a classification `category`
    (key into CATEGORY_COLOR/CATEGORY_MARKER), optionally faceted by
    `ideal_flag` ('yes'/'no'/None) via the shared IDEAL_STYLE encoding
    (filled = ideal, hollow = non-ideal). Returns a dict ready to splat into
    ax.scatter(...).

    `color_map` defaults to CATEGORY_COLOR (the classification-algorithm's
    own predicted-label colors). Callers coloring a REFERENCE/literature
    label instead (fig:bondscores, fig:modemixing) pass ``REF_CATEGORY_COLOR``
    explicitly so "bend"/"stretch" render in the old Okabe-Ito blue/vermillion
    pair rather than the classifier's own hues (2026-07-23 split).

    `marker` optionally overrides the category's own CATEGORY_MARKER shape.
    Used by fig:bondscores/fig:modemixing, where stretch-vs-bend is already
    fully distinguished by color (this function's `category` argument), so
    varying marker SHAPE on top of it would be redundant, over-encoded
    2-factor-in-one-dimension styling -- those two callers pass
    ``marker="o"`` to force one consistent shape for every point regardless
    of category, leaving color (category) and fill (ideal_flag) as the only
    two encodings. CATEGORY_MARKER itself is untouched -- fig:benzene_normal
    and fig:confusion still legitimately vary marker shape by category there,
    since they distinguish >2 mutually-exclusive buckets in one figure.
    """
    color = (color_map or CATEGORY_COLOR)[category]
    if marker is None:
        marker = CATEGORY_MARKER[category]
    if ideal_flag is None:
        return dict(marker=marker, facecolors=color, edgecolors=color,
                    linewidths=0.5, alpha=0.85)
    style = IDEAL_STYLE[ideal_flag]
    face = color if style["filled"] else "none"
    return dict(marker=marker, facecolors=face, edgecolors=color,
                linewidths=0.7, alpha=style["alpha"])


def _explode_bonds(lib_df):
    """Long-format per-bond DataFrame from library_scores.csv's semicolon-
    joined `s_AB`/`rel_db` internal-row strings (see library_ingest.py's
    `_bond_string`). One row per (molecule, mode_index, bond).

    Columns: molecule, mode_index, bond, s_AB, rel_db, ideal, ref_label.
    Rows with missing/empty bond strings (external rows; any malformed
    internal row) are silently skipped -- there is no bond detail to plot
    for them, not a parsing failure.
    """
    internal = lib_df[lib_df["kind"] == "internal"].copy()

    def _parse(s):
        if not isinstance(s, str) or not s:
            return {}
        out = {}
        for part in s.split(";"):
            if not part:
                continue
            bond_id, val = part.split(":")
            out[bond_id] = float(val)
        return out

    records = []
    for _, row in internal.iterrows():
        s_map = _parse(row["s_AB"])
        r_map = _parse(row["rel_db"])
        for bond_id, s_val in s_map.items():
            if bond_id not in r_map:
                continue
            records.append({
                "molecule": row["molecule"], "mode_index": row["mode_index"],
                "bond": bond_id, "s_AB": s_val, "rel_db": r_map[bond_id],
                "ideal": row["ideal"], "ref_label": row["ref_label"],
            })
    return pd.DataFrame.from_records(
        records,
        columns=["molecule", "mode_index", "bond", "s_AB", "rel_db", "ideal", "ref_label"],
    )


# --------------------------------------------------------------------------
# fig:benzene -- benzene EMIT stress test
# --------------------------------------------------------------------------

def plot_benzene_stress_test(
    emit_classified_csv="data/results/benzene_EMIT_classified.csv",
    emit_contrib_csv="data/results/benzene_EMIT_contributions.csv",
    out_dir="data/figures",
    label="fig_benzene",
):
    """Build fig:benzene: score vs. projected normal-mode contribution for
    the 36 EMIT modes, highlighting the EMIT 2/9 s[R] inversion and the
    EMIT 34-36 flagged externals.

    This is a flag-behavior / non-monotonicity illustration (per
    IMPLEMENTATION_PLAN.md Phase 2), NOT a T/R-accuracy parity plot: the
    figure deliberately omits any 1:1 reference line and plots |score|
    against projected contribution fraction only to expose where the two
    disagree.

    Single-panel figure (simplified 2026-07-02): this used to carry a
    second panel showing s[V_S] vs. frequency for benzene's 36 real normal
    modes, but that content is strictly subsumed by the standalone
    ``plot_benzene_normal_modes``/``fig_benzene_normal`` figure (same data,
    plus named worked-example callouts), and the manuscript's "stress test
    on benzene EMIT modes" prose never referenced it -- only the EMIT 2/9
    and EMIT 34-36 content below. Keeping normal-mode data inside a figure
    captioned around the EMIT stress test was also conceptually confusing
    regardless of redundancy. See IMPLEMENTATION_PLAN.md Changelog.

    Returns a summary dict with output paths and a few sanity numbers.
    """
    _style()

    emit = pd.read_csv(emit_classified_csv)
    contrib = pd.read_csv(emit_contrib_csv)
    emit = emit.merge(contrib, on="Mode", suffixes=("", "_c"))

    fig, ax_b = plt.subplots(figsize=(3.8, 3.6))

    # ---------------- score vs projected contribution, EMIT modes ----
    axes6 = ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz"]
    eps = 1e-3
    bg_x, bg_y = [], []
    for _, row in emit.iterrows():
        for ax_name in axes6:
            s = row[ax_name]
            c2 = row["C2_" + ax_name]
            if abs(s) > eps or c2 > eps:
                bg_x.append(c2)
                bg_y.append(abs(s))
    ax_b.scatter(bg_x, bg_y, color=COLORS["background"], marker="o", s=14,
                 alpha=0.6, linewidths=0, zorder=2,
                 label="other EMIT modes\n(all 6 external DOF)")

    # -- Highlight 1: s[R_y] non-monotonicity, EMIT 2 vs EMIT 9 --
    def _row(mode):
        return emit.loc[emit["Mode"] == mode].iloc[0]

    e2, e9 = _row("EMIT 2"), _row("EMIT 9")
    ax_b.scatter([e2["C2_Ry"], e9["C2_Ry"]], [abs(e2["Ry"]), abs(e9["Ry"])],
                 color=COLORS["highlight_r"], marker="D", s=55,
                 edgecolors="black", linewidths=0.5, zorder=5,
                 label=r"EMIT 2 / 9 ($R_y$ inversion)")
    ann_r = ax_b.annotate(
        "", xy=(e9["C2_Ry"], abs(e9["Ry"])), xytext=(e2["C2_Ry"], abs(e2["Ry"])),
        arrowprops=dict(arrowstyle="->", color=COLORS["highlight_r"], lw=1.4),
        zorder=4,
    )
    # highlight_r is now yellow (#F0E442), which has poor contrast as a bare
    # line on a white background -- a thin black halo (path effect) keeps the
    # arrow legible without changing its color (matches the black-edged
    # yellow diamond marker's own contrast treatment above).
    ann_r.arrow_patch.set_path_effects(
        [pe.Stroke(linewidth=2.0, foreground="black"), pe.Normal()])
    # EMIT 9's label offset moved from (6,-10) (lower-right of its marker --
    # exactly where the E2->E9 arrow approaches, and the arrow's new black
    # halo, thicker than the old plain thin line, made the two visually
    # collide) to (-8,12) (upper-left, clear of the arrow's approach
    # direction) -- verified by rendering (2026-07-02 consistency pass).
    ax_b.annotate("EMIT 2", (e2["C2_Ry"], abs(e2["Ry"])),
                  xytext=(6, 6), textcoords="offset points")
    ax_b.annotate("EMIT 9", (e9["C2_Ry"], abs(e9["Ry"])),
                  xytext=(-6, 8), textcoords="offset points",
                  ha="right")

    # -- Highlight 2: flagged external modes, EMIT 34/35/36. By symmetry
    # (Tx/Ty/Tz translation) these three coincide exactly at the same
    # (contribution, score) point -- plotted once, with a side callout table
    # (not per-point labels, which would overlap) giving the s[V_S] readout
    # that distinguishes them.
    flagged = [("EMIT 34", "Tx"), ("EMIT 35", "Ty"), ("EMIT 36", "Tz")]
    fx, fy, fvs = [], [], []
    for mode, ax_name in flagged:
        r = _row(mode)
        fx.append(r["C2_" + ax_name])
        fy.append(abs(r[ax_name]))
        fvs.append(r["V_Stretch"])
    cx, cy = fx[0], fy[0]  # coincident point (all three identical)
    ax_b.scatter([cx], [cy], color=COLORS["highlight_t"], marker="*", s=190,
                 edgecolors="black", linewidths=0.6, zorder=5,
                 label="EMIT 34-36 (flagged external, coincident)")

    callout_xy = (0.47, 0.55)
    ax_b.annotate(
        "", xy=(cx, cy), xytext=callout_xy,
        arrowprops=dict(arrowstyle="-", color=COLORS["highlight_t"],
                        lw=0.9, ls="--"),
        zorder=4,
    )
    callout_text = "EMIT 34-36 (Tx/Ty/Tz):\n" + "\n".join(
        f"  {mode.split()[1]}: " + r"$s[\mathrm{V_S}]$" + f"={vs:.3f}"
        for (mode, _), vs in zip(flagged, fvs)
    )
    ax_b.text(callout_xy[0], callout_xy[1], callout_text,
              ha="center", va="center",
              bbox=dict(boxstyle="round,pad=0.35", fc="white",
                         ec=COLORS["highlight_t"], lw=0.8), zorder=6)

    ax_b.set_xlabel(r"Projected normal-mode contribution, $\tilde{\Theta}_i^2$")
    ax_b.set_ylabel(r"$|\,\mathrm{score}\,|$ (this framework)")
    ax_b.set_xlim(-0.03, 1.05)
    ax_b.set_ylim(-0.03, 1.15)
    ax_b.legend(loc="upper left", frameon=False, handletextpad=0.3,
                labelspacing=0.35, borderaxespad=0.1, fontsize=LEGEND_FONTSIZE)
    # FIXED 2026-07-02 (same readability fix as fig:confusion's footer,
    # IMPLEMENTATION_PLAN.md queued item 1): this prose note used to be
    # drawn in-image via `ax_b.text()` at 6.3pt on this 3.8x3.6in canvas --
    # even smaller than fig:confusion's now-removed footer, so it would
    # shrink well below a readable floor once LaTeX rescales the figure to
    # column width. Removed from the raster; the exact sentence is exposed
    # in `summary["nonmonotonicity_note"]` below for the caption instead.
    nonmonotonicity_note = ("non-monotonic by design (flag mechanism, "
                             "not a parity check)")

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary = {
        "pdf": pdf_path,
        "png": png_path,
        "n_background_points": len(bg_x),
        "emit2_Ry_score": float(abs(e2["Ry"])),
        "emit2_Ry_contribution": float(e2["C2_Ry"]),
        "emit9_Ry_score": float(abs(e9["Ry"])),
        "emit9_Ry_contribution": float(e9["C2_Ry"]),
        "emit34_36_scores": fy,
        "emit34_36_contributions": fx,
        "emit34_36_vstretch": fvs,
        "nonmonotonicity_note": nonmonotonicity_note,
    }
    return summary


# --------------------------------------------------------------------------
# Benzene normal-mode worked-example gallery -- descriptive companion to
# fig:benzene (which stays framed around the EMIT stress test). Standalone,
# single-column figure for the reworked "Benzene normal modes" section
# (JCC/Scoring_Manuscript_Plan_2026-07-02.pdf, "Results & Discussion focus":
# a descriptive worked-example gallery, not a second accuracy report).
# --------------------------------------------------------------------------

# The 3 author-confirmed worked-example modes (IMPLEMENTATION_PLAN.md
# "Benzene normal modes" step (a), done by lead-engineer in
# src/benzene_validation.py Task E / benzene_worked_examples.csv) plus mode
# 19 (already-identified SB example, untouched by that task). Frequencies/
# V_Stretch are read live from ``normal_csv`` (single source of truth) --
# this dict only supplies the descriptive name + non-overlapping callout
# anchor (axes-fraction) for each, chosen by inspecting the 36-point
# scatter's empty regions (see plot_benzene_normal_modes docstring).
#
# REMOVED 2026-07-02 (follow-up correction, same day as the legend-dedup
# fix): the enlarged-marker + dashed-leader-line + callout-text-box
# annotations for these 3 modes were pulled from the rendered figure per the
# author's review -- pointing annotations for modes 12/19/30 will be added
# manually later, possibly alongside a separate depicted-normal-modes figure.
# This dict itself is left defined (unused) since nothing else in the
# codebase references it and it documents which 3 modes were identified in
# step (a); safe to delete later if it becomes dead-code clutter.
_WORKED_EXAMPLE_MODES = {
    "Vib 12": {"short": "12", "name": "ring-breathing", "callout": (0.68, 0.95)},
    "Vib 19": {"short": "19", "name": "mixed S/B", "callout": (0.66, 0.46)},
    "Vib 30": {"short": "30", "name": "C-H stretch", "callout": (0.55, 0.62)},
}


def plot_benzene_normal_modes(
    normal_csv="data/results/benzene_normal_classified.csv",
    out_dir="data/figures",
    label="fig_benzene_normal",
):
    """Build the benzene normal-mode worked-example gallery: ``s[V_S]`` vs.
    frequency for all 36 real normal modes (6 external T/R + 30 internal),
    colored by the full T/R/S/B/SB classification scheme (``CATEGORY_COLOR``/
    ``CATEGORY_LABEL`` mapping shared with ``fig:benzene`` panel (a) and every
    other figure -- not redefined here).

    This is a NEW, separate figure from ``plot_benzene_stress_test``
    (``fig:benzene``, which stays framed around the EMIT stress test and is
    not touched here) -- the descriptive companion for the reworked
    "Benzene normal modes" section (IMPLEMENTATION_PLAN.md, benzene-normal-
    modes 3-step sequence, step (b); step (a) identified the mode indices,
    step (c) is lead-author's narrative rewrite around this figure).

    Styling (follow-up correction, 2026-07-02, same day as the earlier
    legend-dedup fix -- author's second visual review of the rendered
    image):

      * ONE marker shape (circle, ``marker="o"``) for every point regardless
        of category, matching ``plot_bond_scores``/``plot_boxplots``/
        ``plot_mode_mixing``'s established convention that shape is not used
        as a second, redundant encoding of a distinction color already
        carries -- ``CATEGORY_MARKER``'s per-category shapes are NOT used
        here despite still existing for ``fig:confusion``.
      * All markers rendered HOLLOW (``IDEAL_STYLE["no"]`` applied
        uniformly), for cross-figure visual consistency with the ideal/
        non-ideal hollow convention used elsewhere -- there is no ideal/
        non-ideal axis within one molecule's own normal modes, so this is a
        blanket style choice, not a faceted legend split.
      * No enlarged/highlighted markers, dashed leader lines, or callout
        text boxes for modes 12/19/30 -- this is now a plain, clean,
        unannotated scatter (that worked-example annotation block was
        removed; mode-pointing annotations will be added manually later).
    """
    _style()
    normal = pd.read_csv(normal_csv)

    # Widened 3.8->7.2in (2026-07-23, author request): leaves horizontal room
    # to manually composite depicted normal-mode panel images beside the
    # scatter (matching the group's earlier undergraduate report's layout),
    # same rationale/target width as plot_bond_scores' own 2026-07-09 note on
    # widening for this exact purpose. Height (3.6in) unchanged -- the plot
    # itself doesn't need more vertical room, only more width alongside it.
    fig, ax = plt.subplots(figsize=(7.2, 3.6))

    # Legend dedup: translation/rotation share the identical gray color
    # (CATEGORY_COLOR maps both to the same gray -- see that dict's comment
    # above, "always meant to be ONE combined meaning"; both also now render
    # as the same hollow gray circle since marker shape is no longer
    # category-specific here), so a naive per-`cat` dedup (keying on
    # "translation" and "rotation" separately) puts two visually-identical
    # entries ("clean translation", "clean rotation") in the legend --
    # confusing, and unlike every other category here which gets one entry
    # per distinct color. No other figure in this module builds a legend mixing
    # translation+rotation (fig:confusion's tick labels are a structurally
    # different case -- matrix ROWS need the distinct text since they label
    # different ground-truth categories, not a redundant visual encoding), so
    # there is no pre-existing convention to replicate; this is a local-only
    # merge. CATEGORY_LABEL itself stays untouched for that reason.
    _LEGEND_MERGE_KEY = {"translation": "clean_tr", "rotation": "clean_tr"}
    _LEGEND_MERGE_TEXT = {"clean_tr": "clean T/R"}

    seen_labels = set()
    for _, row in normal.iterrows():
        cat = classification_bucket(row["label"])
        leg_key = _LEGEND_MERGE_KEY.get(cat, cat)
        leg_text = _LEGEND_MERGE_TEXT.get(leg_key, CATEGORY_LABEL.get(cat, cat))
        leg_label = leg_text if leg_key not in seen_labels else None
        seen_labels.add(leg_key)
        # ONE consistent marker shape (circle) for every point, all rendered
        # hollow (IDEAL_STYLE["no"], applied uniformly -- no ideal/non-ideal
        # axis within one molecule's own normal modes, just cross-figure
        # visual-consistency styling). Color (CATEGORY_COLOR, via
        # _marker_kwargs) is the only category encoding, matching
        # plot_bond_scores/plot_boxplots/plot_mode_mixing's established
        # "shape is not a redundant second encoding" convention.
        kw = _marker_kwargs(cat, ideal_flag="no", marker="o")
        ax.scatter(row["Freq"], row["V_Stretch"], s=26, zorder=3,
                   label=leg_label, **kw)

    ax.axhline(TAU_S, color=COLORS["threshold"], ls="--", lw=0.8, zorder=1)
    ax.axhline(TAU_B, color=COLORS["threshold"], ls="--", lw=0.8, zorder=1)
    # Threshold labels (follow-up correction, 2026-07-02): previously placed
    # with va="bottom"/va="top" flush against each dashed line, which made
    # the text visually overlap/touch the line itself. Nudged a fixed
    # distance (0.045 in axes-fraction y, comfortably clear at this panel's
    # y-range) above tau_S and below tau_B respectively, so both labels sit
    # in open space next to (not on top of) their line -- verified by
    # rendering the regenerated PNG.
    # tau_S's label is anchored on the LEFT (x=0.02) rather than the right
    # (x=0.98, matched by tau_B below): the S/stretch cluster sits at high
    # frequency AND near s[V_S]=1.0, i.e. exactly the upper-right corner --
    # a right-anchored tau_S label collided with those data points even
    # after being pushed clear of the dashed line itself (checked by
    # rendering). The upper-left corner is empty at this height, so the
    # label is anchored there instead.
    y_offset = 0.045
    ax.text(0.02, TAU_S + y_offset, r"$\tau_\text{S}=$" + f"{TAU_S:.2f}", ha="left",
            va="bottom", color=COLORS["threshold"],
            transform=ax.get_yaxis_transform())
    ax.text(0.98, TAU_B - y_offset, r"$\tau_\text{B}=$" + f"{TAU_B:.2f}", ha="right",
            va="top", color=COLORS["threshold"],
            transform=ax.get_yaxis_transform())

    ax.set_xlabel(r"Frequency (cm$^{-1}$)")
    ax.set_ylabel(r"$s[\mathrm{V_S}]$")
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(-120, normal["Freq"].max() * 1.06)
    # Anchored between the tau_B and tau_S lines (mirrors fig:benzene panel
    # (a)'s "center right" placement, which sits in that same gap) rather
    # than the default "upper left", which would otherwise put the
    # threshold line's full-width dashes straight through the legend text --
    # this frequency range's low-freq bending cluster leaves that band empty
    # on the left.
    ax.legend(loc="center left", bbox_to_anchor=(0.0, 0.52), frameon=False,
              handletextpad=0.3, labelspacing=0.3, borderaxespad=0.2)

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary = {
        "pdf": pdf_path, "png": png_path,
        "shared_categories": ("reuses CATEGORY_COLOR/CATEGORY_LABEL for the "
                               "full T/R/S/B/SB scheme -- same mapping as "
                               "fig:benzene panel (a); a NEW figure, "
                               "fig:benzene itself is untouched. ONE marker "
                               "shape (circle) for all points, all rendered "
                               "hollow (IDEAL_STYLE['no']) -- consistency fix, "
                               "2026-07-02 follow-up: CATEGORY_MARKER's "
                               "per-category shapes are no longer used here "
                               "(color is the only category encoding, "
                               "matching fig:bondscores/fig:boxplots/"
                               "fig:modemixing's convention)."),
        "n_points": len(normal),
        "freq_range": (float(normal["Freq"].min()), float(normal["Freq"].max())),
        "vs_range": (float(normal["V_Stretch"].min()), float(normal["V_Stretch"].max())),
        "tau_S": TAU_S, "tau_B": TAU_B,
        "worked_example_annotations": ("REMOVED 2026-07-02 (author visual "
                                        "review, follow-up correction) -- "
                                        "modes 12/19/30 are no longer "
                                        "highlighted/annotated in this "
                                        "figure; pointing annotations will "
                                        "be added manually later."),
    }
    return summary


# --------------------------------------------------------------------------
# fig:confusion -- SINGLE-TIER (non-ideal only) clean-category confusion
# matrix + retention-and-migration bars, plus a separate SI rigorous-tier
# consistency check (plot_rigorous_tier_check, below).
#
# RESTRUCTURED 2026-07-05 (author decision, executed not re-litigated): the
# earlier 2x2 layout (rigorous heatmap + precision/recall bars, THEN
# non-ideal heatmap + retention/migration bars) put the rigorous-tier
# "precision/recall = 1.000, clears the 0.95 floor" claim in the main text
# as if it were independent validation. It is not: (1) T/R recovery in the
# rigorous tier is guaranteed by construction -- the reference T/R basis is
# built via the same Eckart-Sayvetz projection the mode is then scored
# against, so agreement is close to definitional, not an empirical finding;
# (2) tau_S/tau_B are LITERALLY the min/max of the ideal-molecule s[V_S]
# population (derive_stretch_bend_thresholds, src/calibrate.py) -- the same
# population fig:boxplots already shows has a clean, non-overlapping gap --
# so citing "1.000 accuracy" on that population again in fig:confusion is
# restating the threshold-placement decision, not testing anything new
# against it. The genuine, non-circular validation is the NON-IDEAL tier:
# thresholds fixed on the ideal population, then applied WITHOUT retuning to
# the harder non-ideal cases they were never calibrated on.
#
# This figure therefore keeps ONLY that non-ideal-tier content (what used to
# be panels (c)/(d)), now labeled (a)/(b). The rigorous-tier numbers are not
# deleted -- they still have genuine value as a construction/self-consistency
# check ("did the pipeline wire together correctly on the population that is
# supposed to be exact by construction?") -- but that framing belongs in the
# Supporting Information, not the main text, and is not dressed up as an
# accuracy claim there either. See ``plot_rigorous_tier_check`` below.
# --------------------------------------------------------------------------

def _axis_aware_pred_category(predicted_label):
    """Map a RAW classifier `predicted_label` string (e.g. "Tx", "Tx*", "S",
    "SB") to a confusion-table COLUMN category, WITHOUT collapsing a clean
    external label to the coarse "translation"/"rotation" bucket the way
    `classification_bucket()` does.

    Returns the bare axis ("Tx".."Rz") for a clean external prediction,
    "mixed_external" for a flagged one (any axis), or the ordinary
    `classification_bucket()` result for anything else (internal S/B/SB
    predictions). This is what lets a row's PREDICTION land in an
    axis-specific column independent of that row's own `kind` -- e.g. an
    internal reference mode whose residual character won a Step-2 external
    slot would show up in a "Tx" column here, which a bucket-collapsed
    lookup could never surface. Used by `_joint_confusion_table` below for
    the joint external+internal confusion matrices (fig:confusion,
    fig:benzeneconfusion, the SI rigorous-tier table).
    """
    if is_clean_external(predicted_label):
        return external_axis(predicted_label)
    if is_mixed_external(predicted_label):
        return "mixed_external"
    return classification_bucket(predicted_label)


_EXTERNAL_AXES = {"Tx", "Ty", "Tz", "Rx", "Ry", "Rz"}


def _joint_confusion_table(df, ref_order, pred_order, collapse_external=False):
    """Build a reference-label x predicted-category confusion table spanning
    BOTH external (Tx..Rz) and internal (stretch/bend/mixed/SB) categories in
    one crosstab -- the joint matrix that actually checks whether any
    internal mode's residual character ever wins an external Step-2 slot (or
    vice versa), not just whether each population is internally consistent.

    Ground truth: `mode_index` (the true axis identity, e.g. "Tx") for
    external rows, `ref_label` (the literature/derived stretch/bend/SB label)
    for internal rows. Prediction: `predicted_label` mapped through
    `_axis_aware_pred_category`. Reindexed to `ref_order`/`pred_order`
    (fill_value=0) -- any category present in the data but not in the
    requested order is silently dropped from the returned table (callers
    pass an exhaustive order for the categories they expect).

    `collapse_external=True` folds all 6 individual axes (both ref and pred
    side) into one "T/R" category before the crosstab -- used by
    fig:benzeneconfusion, whose per-axis detail is redundant with
    fig:confusion's non-ideal matrix.
    """
    ref_cat = np.where(df["kind"] == "external", df["mode_index"], df["ref_label"])
    pred_cat = df["predicted_label"].map(_axis_aware_pred_category)
    if collapse_external:
        ref_cat = np.array(["T/R" if c in _EXTERNAL_AXES else c for c in ref_cat])
        pred_cat = pred_cat.map(lambda c: "T/R" if c in _EXTERNAL_AXES else c)
    tbl = pd.crosstab(pd.Series(ref_cat, index=df.index, name="ref"), pred_cat)
    return tbl.reindex(index=ref_order, columns=pred_order, fill_value=0)


def _per_category_from_table(tbl, cats):
    """Generic precision/recall/n_reference/n_predicted-correct per category,
    computed directly from an already-built confusion table's own row/column
    sums and diagonal -- unlike `confusion_matrix_stats`'s per-category dict
    (src/calibrate.py), which is hardcoded to exactly 4 bucket-collapsed
    categories (translation/rotation/stretch/bend) and would KeyError on an
    axis-specific key like "Tx". Works for any `cats` list that indexes both
    `tbl`'s rows and columns (true for every joint table this module builds,
    since `ref_order`/`pred_order` always share the same category names for
    the categories being scored).

    Returns a list of dicts, one per category in `cats`, each:
    {category, n_reference, n_predicted_correct, precision, recall}.
    """
    rows = []
    for c in cats:
        n_ref = int(tbl.loc[c].sum()) if c in tbl.index else 0
        n_pred = int(tbl[c].sum()) if c in tbl.columns else 0
        n_correct = int(tbl.loc[c, c]) if (c in tbl.index and c in tbl.columns) else 0
        precision = n_correct / n_pred if n_pred else float("nan")
        recall = n_correct / n_ref if n_ref else float("nan")
        rows.append({
            "category": CATEGORY_LABEL.get(c, c), "n_reference": n_ref,
            "n_predicted_correct": n_correct, "precision": precision,
            "recall": recall,
        })
    return rows


def _confusion_heatmap(ax, fig, tbl, ref_order, label_map=CATEGORY_LABEL,
                        divider=None):
    """Shared heatmap renderer for one confusion-matrix tier. `tbl` must
    already be reindexed to `ref_order` rows (columns are whatever buckets
    are present for that tier -- rigorous and non-ideal tiers populate
    different bucket sets, so no forced union of columns across tiers).
    Returns the imshow handle (caller attaches its own colorbar).

    No in-figure title (2026-07-23, author request): the mode-count context
    it used to carry (e.g. "675 normal modes of non-ideal molecules") now
    lives in the LaTeX caption only -- this is a single-panel heatmap, so no
    (a)/(b)/(c) label is needed either.

    `label_map` defaults to the shared `CATEGORY_LABEL` dict but
    `plot_confusion_matrix` passes a figure-local override (see that
    function's `_CONFUSION_LABEL`) so this figure's translation/rotation
    ticks read "T"/"R" instead of the long "clean translation"/"clean
    rotation" form, without touching `CATEGORY_LABEL` itself.

    `divider=(n_ref_ext, n_pred_ext)`, if given, draws a thicker line after
    row index `n_ref_ext-1` and column index `n_pred_ext-1`, separating a
    leading external-axis block (rows: Tx..Rz; columns: Tx..Rz +
    mixed_external, grouped together) from the trailing internal-category
    block -- used by the joint external+internal confusion matrices.
    """
    vmax = max(1, tbl.values.max())
    im = ax.imshow(tbl.values, cmap=COLORS["confusion_cmap"], aspect="auto",
                    vmin=0, vmax=vmax)
    for i in range(tbl.shape[0]):
        for j in range(tbl.shape[1]):
            v = int(tbl.values[i, j])
            txt_color = "white" if v > 0.6 * vmax else "black"
            ax.text(j, i, str(v), ha="center", va="center",
                     color=txt_color)

    # Rotation/right-alignment used to be needed to fit the long "clean
    # translation"/"clean rotation" tick text (and, before that, the
    # briefly-tried "B (bending)"-style gloss form) without collisions.
    # Now that `label_map` gives every column/row a short (<=2 char) symbol
    # (T/R/S/B/SB), the labels fit horizontally with no rotation -- fixed
    # 2026-07-02 (author follow-up: shorten T/R ticks + restore readable
    # sizing), matching fig:boxplots' identical rotation-removal precedent
    # once its own tick labels were shortened back to bare S/B.
    ax.set_xticks(range(len(tbl.columns)))
    ax.set_xticklabels([label_map[PRED_BUCKET_TO_CATEGORY[c]] if c in
                         PRED_BUCKET_TO_CATEGORY else c for c in tbl.columns],
                        rotation=0, ha="center")
    ax.set_yticks(range(len(tbl.index)))
    ax.set_yticklabels([label_map[REF_LABEL_TO_CATEGORY[r]] for r in tbl.index])
    # Row ticks ("Reference label" axis) use REF_CATEGORY_COLOR (old Okabe-Ito
    # blue/vermillion); column ticks ("Classification algorithm label" axis)
    # use CATEGORY_COLOR (new blue) -- the 2026-07-23 ref/predicted color
    # split, precisely the two axes this heatmap already labels below.
    for tick, r in zip(ax.get_yticklabels(), tbl.index):
        tick.set_color(REF_CATEGORY_COLOR[REF_LABEL_TO_CATEGORY[r]])
    for tick, c in zip(ax.get_xticklabels(), tbl.columns):
        cat = PRED_BUCKET_TO_CATEGORY.get(c)
        if cat is not None:
            tick.set_color(CATEGORY_COLOR[cat])
    ax.set_xlabel("Classification algorithm label")
    ax.set_ylabel("Reference label")
    if divider is not None:
        n_ref_ext, n_pred_ext = divider
        ax.axhline(n_ref_ext - 0.5, color="black", lw=1.6, zorder=4)
        ax.axvline(n_pred_ext - 0.5, color="black", lw=1.6, zorder=4)
    # Borders removed (author request, 2026-07-23): imshow's own cell grid
    # already delineates the matrix, so the outer axes frame is redundant --
    # same for the colorbar's outline rectangle. `_style()`'s global
    # top/right-spine removal doesn't reach bottom/left, so all 4 are turned
    # off explicitly here rather than relying on that shared setting.
    for spine in ax.spines.values():
        spine.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Number of modes")
    cb.outline.set_visible(False)
    return im


def plot_confusion_matrix(
    library_csv="data/results/library_scores.csv",
    out_dir="data/figures",
    label="fig_confusion",
):
    """Build fig:confusion as a SINGLE-PANEL joint confusion matrix
    (restructured again, this session): the non-ideal tier's internal rows
    (``ideal == 'no'``) PLUS the external (Tx..Rz) reference rows for those
    same non-ideal, single-centre molecules, in ONE crosstab -- not just an
    internal-only stretch/bend/mixed table.

    Why joint, not internal-only (author's own framing): showing external
    and internal categories together is what actually verifies that no
    internal mode's residual character ever wins an external Step-2 slot
    (or vice versa) -- a genuine completeness check on the classifier's
    global assignment, not merely "is each population self-consistent."
    The off-diagonal blocks BETWEEN the external-axis rows/columns and the
    internal rows/columns are computed and reported (see
    ``crossover_ext_ref_to_internal_pred``/``crossover_internal_ref_to_ext_pred``
    in the returned summary) rather than assumed to be zero.

    Thresholds are fixed on the ideal population (fig:boxplots), then
    applied WITHOUT retuning to the harder non-ideal cases they were never
    calibrated on -- this remains the genuine, non-circular validation; the
    fully-rigorous tier (external T/R rows + ``ideal == 'yes'`` internal
    rows) is exact by construction and is reported separately, as a
    self-consistency check rather than an accuracy claim, by
    ``plot_rigorous_tier_check`` below (SI-bound, not main text).

    Single panel now (retention/migration bars, formerly panel (b), were
    split out this session to their own SI figure,
    ``plot_confusion_retention_migration`` below -- the mode-count
    breakdown this panel shows already carries that information via its own
    row/column sums, and the retention numbers are stated in prose).

    Never recomputes scores. ``filter_single_centre_library`` is applied
    explicitly here (not left to confusion_matrix_stats's internal filter)
    so the SAME molecule scope is used for both the internal rows and the
    external-row merge -- idempotent if already filtered.
    """
    _style()
    from src.calibrate import confusion_matrix_stats, filter_single_centre_library

    lib_df = pd.read_csv(library_csv)
    lib_df = filter_single_centre_library(lib_df)
    thresholds = Thresholds.calibrated()

    nonideal_df = lib_df[(lib_df["kind"] == "internal") & (lib_df["ideal"] == "no")]
    stats_n = confusion_matrix_stats(nonideal_df, thresholds)

    nonideal_molecules = nonideal_df["molecule"].unique()
    nonideal_external_df = lib_df[(lib_df["kind"] == "external") &
                                   (lib_df["molecule"].isin(nonideal_molecules))]
    joint_df = pd.concat([nonideal_df, nonideal_external_df], ignore_index=True)

    # "mixed_external" ("T/R*") dropped from pred_order_n (author decision,
    # 2026-07-09): the column is empty by construction on this population
    # (verified, not assumed -- confirmed 0 across all non-ideal molecules)
    # and just widens the table with no information.
    ref_order_n = ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "stretch", "bend"]
    pred_order_n = ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "stretch", "bend", "mixed"]
    tbl_n = _joint_confusion_table(joint_df, ref_order_n, pred_order_n)
    n_nonideal = int(tbl_n.values.sum())

    # Height tuned down from an initial 5.6in (author visual check, compiled
    # PDF): at 5.6in the figure left too little room on its page for the
    # multi-line caption to start, pushing the whole caption to the next
    # page and leaving a mostly-blank page behind -- 4.8in lets the figure
    # and its caption share one page.
    fig, ax_hb = plt.subplots(figsize=(6.2, 4.8))
    _confusion_heatmap(ax_hb, fig, tbl_n, ref_order_n)

    # Internal-only opposite-category crossing (bend<->stretch), same
    # computation as before the joint-matrix restructure -- these row/column
    # labels are unchanged within the bigger table, so this still works
    # unmodified against `tbl_n`.
    cats_n = ["bend", "stretch"]
    opposite_n = []
    for c in cats_n:
        opp = "stretch" if c == "bend" else "bend"
        n_ref = stats_n["per_category"][c]["n_ref"]
        n_opp = int(tbl_n.loc[c, opp]) if opp in tbl_n.columns else 0
        opposite_n.append(n_opp / n_ref if n_ref else float("nan"))

    # NEW this session: external<->internal crossover check -- the actual
    # point of merging the two populations into one table. Expected 0 both
    # ways, but computed, not assumed.
    ext_rows, int_rows = ref_order_n[:6], ref_order_n[6:]
    ext_cols, int_cols = pred_order_n[:6], pred_order_n[6:]
    crossover_ext_ref_to_internal_pred = int(tbl_n.loc[ext_rows, int_cols].values.sum())
    crossover_internal_ref_to_ext_pred = int(tbl_n.loc[int_rows, ext_cols].values.sum())

    # Caption-text sentence (not drawn in-image, same readability rationale
    # as the earlier footer_text -- see `summary["nonideal_footer_text"]`).
    # ASCII "->" deliberately, not a unicode arrow (console-print safety).
    footer_text = (
        f"Non-ideal internal + external classification (n={n_nonideal}): 0% of "
        f"bend or stretch reference-labeled modes crossed to the OPPOSITE clean "
        f"category (bend->stretch={opposite_n[0]:.1%}, stretch->bend={opposite_n[1]:.1%}); "
        f"and {crossover_ext_ref_to_internal_pred} external reference modes were "
        f"predicted into an internal bucket, {crossover_internal_ref_to_ext_pred} "
        "internal reference modes were predicted into an external slot -- the "
        "off-diagonal blocks between the external and internal categories are "
        "empty. Thresholds were fixed on the ideal-molecule population "
        "(fig:boxplots) and applied here without retuning."
    )
    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary = {
        "pdf": pdf_path, "png": png_path,
        "nonideal_footer_text": footer_text,
        "layout": ("Single panel: joint external (Tx..Rz) + internal "
                   "(stretch/bend/mixed) confusion matrix -- restructured "
                   "this session from an internal-only 1x2 heatmap+bars "
                   "layout. Retention/migration bars (formerly panel (b)) "
                   "moved to plot_confusion_retention_migration (SI, not "
                   "main text). No divider line between the external and "
                   "internal blocks, and no mixed_external ('T/R*') column "
                   "(author decision, 2026-07-09; the column was empty)."),
        "shared_categories": ("reuses CATEGORY_COLOR/CATEGORY_LABEL for "
                               "STRETCHING/BENDING/MIXED_STRETCH_BEND and the "
                               "6 individual-axis entries -- same mapping as "
                               "every other figure."),
        "nonideal_n": n_nonideal,
        "nonideal_confusion_table": tbl_n.to_dict(),
        "nonideal_opposite_category_crossing": dict(zip(cats_n, opposite_n)),
        "crossover_ext_ref_to_internal_pred": crossover_ext_ref_to_internal_pred,
        "crossover_internal_ref_to_ext_pred": crossover_internal_ref_to_ext_pred,
    }
    return summary


def plot_confusion_retention_migration(
    library_csv="data/results/library_scores.csv",
    out_dir="data/figures",
    label="fig_confusion_retention_migration",
):
    """Build the SI companion to fig:confusion (retention-vs.-migration
    bars): the non-ideal-tier label-retention vs. migration-to-mixed bar
    panel, SPLIT OUT of ``plot_confusion_matrix`` this session (formerly that
    figure's panel (b)) -- mirrors the existing
    ``plot_benzene_confusion_precision_recall`` SI-split precedent below.
    Self-contained: recomputes ``confusion_matrix_stats`` itself rather than
    depending on ``plot_confusion_matrix`` having already run, matching this
    module's "each function reads already-computed CSVs, never recomputes
    scores" convention.

    Same underlying numbers as before (bend retention 0.896, stretch
    retention 0.503, etc.) -- only the visualization moved out of the main
    text; Figure 6's mode-count breakdown and the manuscript prose already
    state these numbers, per the author's own framing.
    """
    _style()
    from src.calibrate import confusion_matrix_stats, filter_single_centre_library

    lib_df = pd.read_csv(library_csv)
    lib_df = filter_single_centre_library(lib_df)
    thresholds = Thresholds.calibrated()

    nonideal_df = lib_df[(lib_df["kind"] == "internal") & (lib_df["ideal"] == "no")]
    stats_n = confusion_matrix_stats(nonideal_df, thresholds)
    n_nonideal = int(stats_n["confusion_table"].values.sum())

    cats_n = ["bend", "stretch"]
    retention_n = [stats_n["per_category"][c]["recall"] for c in cats_n]
    migration_n = [stats_n["per_category"][c]["mixed_fraction"] for c in cats_n]

    fig, ax_pb = plt.subplots(figsize=(3.6, 3.4))
    x_n = np.arange(len(cats_n))
    # Reference-category bars ("retained" as bend/stretch, the true class),
    # not the classifier's predicted-label axis -- old blue/vermillion pair.
    bar_colors_n = [REF_CATEGORY_COLOR[REF_LABEL_TO_CATEGORY[c]] for c in cats_n]
    ax_pb.bar(x_n, retention_n, 0.5, color=bar_colors_n,
              edgecolor="black", linewidth=0.5)
    ax_pb.bar(x_n, migration_n, 0.5, bottom=retention_n,
              color=COLORS["mixed"], alpha=0.85, edgecolor="black", linewidth=0.5)
    THIN = 0.08  # segments shorter than this (axis fraction) get an outside label
    for xi, ret, mig in zip(x_n, retention_n, migration_n):
        ax_pb.text(xi, ret / 2, f"retained\n{ret:.3f}", ha="center", va="center",
                   color="white", fontweight="bold", fontsize=ANNOTATION_FONTSIZE)
        if mig >= THIN:
            ax_pb.text(xi, ret + mig / 2, f"mixed\n{mig:.3f}", ha="center",
                       va="center", color="white", fontweight="bold",
                       fontsize=ANNOTATION_FONTSIZE)
        else:
            ax_pb.text(xi, ret + mig + 0.015, f"mixed: {mig:.3f}", ha="center",
                       va="bottom", color=COLORS["mixed"],
                       fontweight="bold", fontsize=ANNOTATION_FONTSIZE)
    ax_pb.set_xticks(x_n)
    ax_pb.set_xticklabels([CATEGORY_LABEL[REF_LABEL_TO_CATEGORY[c]] for c in cats_n])
    ax_pb.set_xlim(-0.55, 1.55)
    ax_pb.set_ylim(0, 1.12)
    ax_pb.set_ylabel("Fraction of modes")

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary = {
        "pdf": pdf_path, "png": png_path,
        "nonideal_n": n_nonideal,
        "nonideal_retention": dict(zip(cats_n, retention_n)),
        "nonideal_migration_to_mixed": dict(zip(cats_n, migration_n)),
        "framing": ("SI companion, split out of plot_confusion_matrix "
                     "(fig:confusion) this session -- same numbers, just "
                     "moved out of the main-text figure since Figure 6's own "
                     "mode-count breakdown and the manuscript prose already "
                     "state the retention numbers."),
    }
    return summary


def plot_rigorous_tier_check(
    library_csv="data/results/library_scores.csv",
    out_dir="data/figures",
    label="fig_rigorous_tier_check",
):
    """Build the SI rigorous-tier consistency check (companion to the
    restructured, single-panel ``fig:confusion`` above): a small reference/
    predicted-count TABLE (not a heatmap+bars figure like fig:confusion --
    deliberately minimal, per the author's 2026-07-05 restructuring
    decision) covering every external (T/R) row plus internal rows with
    ``ideal == 'yes'`` (n computed fresh each run from the current library
    roster -- do not hardcode a specific count in this docstring, it has
    drifted stale before as the roster grew/shrank).

    Framing (explicit, both in this docstring and in the rendered table's
    own caption row): this is a SANITY CHECK confirming the construction is
    self-consistent -- T/R references are built via the same Eckart-Sayvetz
    projection the mode is then scored against, and tau_S/tau_B are
    literally the min/max of this exact population (see
    src.calibrate.derive_stretch_bend_thresholds) -- NOT an independent
    accuracy claim clearing a floor. Precision/recall are 1.000 for all
    categories by construction; this table exists to show that
    computationally rather than merely asserting it, not to argue it proves
    anything about the harder non-ideal cases (that is fig:confusion's job).

    EXPANDED this session: the reference/predicted categories are now the 6
    individual axes (Tx..Rz), not the 2 coarse translation/rotation buckets
    -- this table already joined external+internal categories before this
    change (unlike fig:confusion/fig:benzeneconfusion, which needed a
    separate restructure to become joint); widening its own external side
    from 2 buckets to 6 axes is a smaller, more contained extension of an
    already-joint table.

    Renders BOTH a compact table-style PDF/PNG (drop-in
    ``\\includegraphics``, matching this module's existing figure-asset
    convention) AND writes the same numbers to
    ``data/results/rigorous_tier_consistency_table.csv`` so lead-author/
    tex-data-sync can instead typeset a native LaTeX ``booktabs`` table (the
    convention already used by ``JCC_SI_computational_cost.tex``) if
    preferred -- this function does not choose that for them, only supplies
    both forms of the same numbers. Not one of the 6 originally-scoped
    fig:* labels; SI placement (a new small SI document, or a table inside
    an existing one) is lead-author's call, not invented here.
    """
    _style()
    from src.calibrate import confusion_matrix_stats, filter_single_centre_library

    lib_df = pd.read_csv(library_csv)
    lib_df = filter_single_centre_library(lib_df)
    thresholds = Thresholds.calibrated()

    rigorous_df = lib_df[(lib_df["kind"] == "external") | (lib_df["ideal"] == "yes")]
    # Still used for acceptance_floor/floor_met (the coarse 4-bucket
    # construction check) -- unaffected by the axis-level table expansion
    # below, since that check was never about per-axis granularity.
    stats_r = confusion_matrix_stats(rigorous_df, thresholds)

    cats_r = ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "stretch", "bend"]
    tbl_r = _joint_confusion_table(rigorous_df, cats_r, cats_r)
    n_rigorous = int(tbl_r.values.sum())
    rows = _per_category_from_table(tbl_r, cats_r)
    table_df = pd.DataFrame(rows)

    csv_path = os.path.join("data", "results", "rigorous_tier_consistency_table.csv")
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    table_df.to_csv(csv_path, index=False)

    # Minimal table rendering (matplotlib Table, not a heatmap+bars figure --
    # per the author's explicit "doesn't need a full ... figure like the
    # main-text one did" instruction).
    #
    # FIXED (post-render visual check): the first attempt let `ax.table`
    # auto-size columns with no explicit `colWidths`, which badly overlapped
    # the "n (reference)" / "n (predicted correctly)" header text into one
    # illegible smear. Explicit, hand-tuned `colWidths` (summing to 1.0,
    # proportional to each header's rendered length) plus a wider figure
    # (4.6in -> 7.0in) fixes this -- verified by rendering the regenerated
    # PNG below. Height increased (1.9 -> 3.2) this session for the
    # expansion from 4 rows (translation/rotation/stretch/bend) to 8
    # (Tx..Rz/stretch/bend) -- column widths themselves are unaffected by
    # row count, only re-tuned if a header changes.
    fig, ax = plt.subplots(figsize=(7.0, 3.2))
    ax.axis("off")
    col_labels = ["Category", "n (reference)", "n (predicted\ncorrectly)",
                  "Precision", "Recall"]
    col_widths = [0.20, 0.24, 0.26, 0.15, 0.15]
    cell_text = [[r["category"], f'{r["n_reference"]:d}',
                  f'{r["n_predicted_correct"]:d}', f'{r["precision"]:.3f}',
                  f'{r["recall"]:.3f}'] for r in rows]
    tbl = ax.table(cellText=cell_text, colLabels=col_labels, colWidths=col_widths,
                   loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1.0, 1.8)
    for (row, col), cell in tbl.get_celld().items():
        cell.set_edgecolor("#BBBBBB")
        if row == 0:
            cell.set_facecolor("#EEEEEE")
            cell.set_text_props(fontweight="bold")

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary = {
        "pdf": pdf_path, "png": png_path, "csv": csv_path,
        "framing": ("Self-consistency / construction check, NOT an "
                     "independent accuracy claim -- T/R references are "
                     "built from the same Eckart-Sayvetz projection the "
                     "mode is scored against, and tau_S/tau_B are the "
                     "min/max of exactly this ideal-molecule population "
                     "(see fig:boxplots) -- so 1.000 precision/recall here "
                     "restates the construction, it does not test it "
                     "against anything new. The genuine non-circular "
                     "validation is fig:confusion's non-ideal tier."),
        "n_rigorous": n_rigorous,
        "confusion_table": tbl_r.to_dict(),
        "precision": {r["category"]: r["precision"] for r in rows},
        "recall": {r["category"]: r["recall"] for r in rows},
        "acceptance_floor": stats_r["acceptance_floor"],
        "floor_met": stats_r["floor_met"],
    }
    return summary


# --------------------------------------------------------------------------
# fig:benzeneconfusion -- benzene's genuine 3-class (bend/stretch/SB)
# internal-mode confusion matrix, enabled by the 2026-07-05 literature
# relabeling of modes 21/22 as a literal "SB" (mixed) ground truth (commit
# 69d549e). NEW figure (2026-07-05) -- does not modify plot_confusion_matrix
# (fig:confusion, the whole-hydride-library validation) or any of its inputs.
# --------------------------------------------------------------------------

def plot_benzene_internal_confusion(
    out_dir="data/figures",
    label="fig_benzene_confusion",
):
    """Build fig:benzeneconfusion: benzene's (C6H6) JOINT confusion matrix --
    the 6 external (Tx..Rz) normal modes together with the 30 internal
    modes' reference bend/stretch/SB (literal literature label, modes 21/22,
    an E1u degenerate pair at 1532.85 cm-1) x predicted bend/stretch/mixed
    bucket, in ONE crosstab (restructured this session; was internal-only
    3x3 before).

    Why joint (same rationale as fig:confusion's identical restructure this
    session): the manuscript prose already states benzene's 6 external
    modes are "recovered exactly (6/6)," but that was only ever a sentence,
    not a checked matrix cell. Building the joint matrix turns it into one:
    the off-diagonal blocks between the 6 external rows/columns and the
    internal rows/columns are computed and reported (see
    ``crossover_ext_ref_to_internal_pred``/``crossover_internal_ref_to_ext_pred``)
    rather than assumed zero.

    Data source: ``src.benzene_validation.benzene_normal_reference_detail()``,
    which already returns all 36 rows (6 external + 30 internal) with
    ``mode_index``/``kind``/``ref_label``/``predicted_label`` -- exactly the
    columns ``_joint_confusion_table`` needs. Not recomputed here; this
    figure remains presentation-only. The old pre-built 3x3
    ``benzene_internal_confusion_matrix.csv`` is no longer read by this
    function (``plot_benzene_confusion_precision_recall`` below still reads
    it independently, unaffected).

    Distinct from ``fig:confusion`` (``plot_confusion_matrix``): that figure
    is the whole-hydride-library (non-ideal-tier) validation; this one is
    benzene's own 36 modes only, using benzene's genuine literature 3-class
    internal ground truth (bend/stretch/SB) plus its own 6 exact-by-
    construction external references. Companion per-bond numerical detail
    for the two flagged contrast cases (modes 21/22 -- literature SB,
    predicted clean "B", both below tau_B; modes 23/24 -- literature
    stretch, predicted "SB"/mixed) lives in
    ``data/results/benzene_sb_vs_stretch_bond_diagnostic.csv`` for a
    companion LaTeX table -- not re-plotted inside this heatmap figure.
    """
    _style()
    from src.benzene_validation import benzene_normal_reference_detail

    detail = benzene_normal_reference_detail()

    # Collapsed to a single "T/R" row/column (author decision, 2026-07-09):
    # per-axis (Tx..Rz) detail is already shown in fig:confusion's non-ideal
    # matrix, so this figure only needs to confirm benzene's 6 external
    # modes land as a block, not which axis each one lands on. "mixed_external"
    # ("T/R*") is dropped from pred_order for the same reason as fig:confusion
    # -- verified empty on benzene's 6 external modes.
    ref_order = ["T/R", "bend", "stretch", "SB"]
    pred_order = ["T/R", "bend", "stretch", "mixed"]
    tbl = _joint_confusion_table(detail, ref_order, pred_order, collapse_external=True)
    n_total = int(tbl.values.sum())

    fig, ax_h = plt.subplots(figsize=(4.6, 4.2))
    _confusion_heatmap(ax_h, fig, tbl, ref_order)

    ext_rows, int_rows = ref_order[:1], ref_order[1:]
    ext_cols, int_cols = pred_order[:1], pred_order[1:]
    crossover_ext_ref_to_internal_pred = int(tbl.loc[ext_rows, int_cols].values.sum())
    crossover_internal_ref_to_ext_pred = int(tbl.loc[int_rows, ext_cols].values.sum())

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary_dict = {
        "pdf": pdf_path, "png": png_path,
        "shared_categories": ("reuses CATEGORY_COLOR/CATEGORY_LABEL/"
                               "_confusion_heatmap from plot_confusion_matrix "
                               "(fig:confusion) via the "
                               "REF_LABEL_TO_CATEGORY['SB']='mixed' routing "
                               "and the collapsed 'T/R' entry -- does not "
                               "modify fig:confusion itself."),
        "n_total": n_total,
        "confusion_table": tbl.to_dict(),
        "crossover_ext_ref_to_internal_pred": crossover_ext_ref_to_internal_pred,
        "crossover_internal_ref_to_ext_pred": crossover_internal_ref_to_ext_pred,
        "note": ("4x4 joint matrix: benzene's 6 external T/R modes collapsed "
                 "to a single T/R row/column (2026-07-09, author decision -- "
                 "per-axis detail lives in fig:confusion instead), turning "
                 "the prose's '6/6 recovered exactly' claim into an actual "
                 "checked matrix cell without repeating fig:confusion's "
                 "per-axis breakdown. No divider line and no mixed_external "
                 "('T/R*') column (verified empty). Precision/recall bar "
                 "panel remains split out to plot_benzene_confusion_"
                 "precision_recall (SI, not main text; unaffected by this "
                 "change). Companion per-bond evidence for the 21/22 "
                 "(SB->bend blind spot) vs. 23/24 (stretch->mixed) contrast "
                 "lives in data/results/benzene_sb_vs_stretch_bond_"
                 "diagnostic.csv for a separate LaTeX table -- not plotted "
                 "here."),
    }
    return summary_dict


def plot_benzene_confusion_precision_recall(
    matrix_csv="data/results/benzene_internal_confusion_matrix.csv",
    summary_csv="data/results/benzene_internal_confusion_summary.csv",
    out_dir="data/figures",
    label="fig_benzene_precision_recall",
):
    """Build the SI companion to fig:benzeneconfusion (10th figure, no
    ``fig:`` label of its own -- SI-bound, not main text): benzene's own
    internal-mode per-category precision/recall bar panel, SPLIT OUT of
    ``plot_benzene_internal_confusion`` 2026-07-05 (author decision, tone/
    framing only -- see that function's docstring and this module's header
    comment for the full rationale). Same numbers, unchanged, just moved out
    of the main-text figure and into the Supporting Information: benzene's
    literature stretch/bend/SB ground truth (Shi 1972) is genuinely external,
    non-circular ground truth (unlike the ideal-population rigorous tier
    ``plot_rigorous_tier_check`` guards against), so this move is NOT the
    same circularity argument -- it is purely that a formal precision/recall
    bar chart reads as an accuracy-metric claim sitting awkwardly beside this
    section's own deliberate "calling a mode 'mixed' by eye is a convention,
    not an exact measurement" hedge.

    Layout (1x1): grouped precision/recall bars for 3 categories (bend,
    stretch, and "mixed" -- the predicted-bucket name for what is a literal
    "SB" reference row). Precision is computed here directly from the
    confusion table (correctly-predicted / total-predicted-that-bucket);
    recall comes straight from ``summary_csv``'s own ``recall`` column (same
    numbers ``src.benzene_validation.benzene_internal_confusion_matrix``
    already computes, not re-derived).
    """
    _style()
    tbl = pd.read_csv(matrix_csv, index_col=0)
    summary = pd.read_csv(summary_csv).set_index("ref_label")

    ref_order = ["bend", "stretch", "SB"]
    pred_order = ["bend", "stretch", "mixed"]
    tbl = tbl.reindex(index=ref_order, columns=pred_order, fill_value=0)

    # Each bar-panel category pairs one reference ROW with one predicted
    # COLUMN of the SAME conceptual category: bend<->bend, stretch<->stretch,
    # and SB(reference)<->mixed(predicted) -- the classifier has no predicted
    # bucket literally named "SB", so "mixed" is its structural equivalent
    # (matches src.benzene_validation._expected_pred_bucket's own convention).
    cats = ["bend", "stretch", "mixed"]
    row_for_cat = {"bend": "bend", "stretch": "stretch", "mixed": "SB"}
    precisions, recalls = [], []
    for cat in cats:
        row, col = row_for_cat[cat], cat
        tp = int(tbl.loc[row, col])
        n_pred = int(tbl[col].sum())
        precisions.append(tp / n_pred if n_pred else float("nan"))
        recalls.append(float(summary.loc[row, "recall"]))

    fig, ax_p = plt.subplots(figsize=(3.9, 3.4))
    x = np.arange(len(cats))
    width = 0.35
    bar_colors = [CATEGORY_COLOR[c] for c in cats]
    ax_p.bar(x - width / 2, precisions, width, color=bar_colors,
             edgecolor="black", linewidth=0.5, label="precision")
    ax_p.bar(x + width / 2, recalls, width, color=bar_colors,
             edgecolor="black", linewidth=0.5, hatch="///", label="recall")
    ax_p.set_xticks(x)
    ax_p.set_xticklabels([CATEGORY_LABEL[c] for c in cats])
    ax_p.set_ylim(0, 1.18)
    ax_p.set_ylabel("Precision / recall")
    ax_p.legend(loc="upper left", frameon=False, fontsize=LEGEND_FONTSIZE)
    for xi, p, r in zip(x, precisions, recalls):
        ax_p.text(xi - width / 2, p + 0.02, f"{p:.3f}", ha="center", va="bottom",
                  fontsize=ANNOTATION_FONTSIZE)
        ax_p.text(xi + width / 2, r + 0.02, f"{r:.3f}", ha="center", va="bottom",
                  fontsize=ANNOTATION_FONTSIZE)

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary_dict = {
        "pdf": pdf_path, "png": png_path,
        "shared_categories": ("reuses CATEGORY_COLOR/CATEGORY_LABEL from "
                               "plot_confusion_matrix (fig:confusion) -- "
                               "does not modify fig:confusion itself."),
        "precision": dict(zip(cats, precisions)),
        "recall": dict(zip(cats, recalls)),
        "framing": ("SI companion, not an independent accuracy claim beyond "
                     "what fig:benzeneconfusion's heatmap already shows -- "
                     "benzene's literature ground truth (Shi 1972) is "
                     "genuinely external and non-circular, so this bar panel "
                     "is a legitimate precision/recall report; it was moved "
                     "out of the main text purely for tone/framing "
                     "consistency with the surrounding hedging prose, not "
                     "because the numbers are suspect."),
    }
    return summary_dict


# --------------------------------------------------------------------------
# fig:bondscores -- bond score vs relative Delta-bond-length
# --------------------------------------------------------------------------

def plot_bond_scores(
    library_csv="data/results/library_scores.csv",
    out_dir="data/figures",
    label="fig_bondscores",
):
    """Build fig:bondscores: per-bond score s_AB vs. the absolute relative
    change in bond length, for every bond in every internal mode of the
    hydride library, split ideal (filled) vs. non-ideal (hollow) and
    stretching (vermillion) vs. bending (blue) -- reusing the shared
    CATEGORY_COLOR and IDEAL_STYLE encodings. All points use ONE consistent
    marker shape (circle) regardless of stretch/bend: color already fully
    distinguishes that 2-way split, so varying marker shape on top of it
    would be redundant over-encoding of the same distinction (fixed
    2026-07-02 consistency pass; CATEGORY_MARKER's square/circle shapes are
    still used elsewhere, e.g. fig:confusion/fig:benzene_normal, where shape
    legitimately distinguishes >2 buckets).

    Restricted to the single-centre AB_n hydride-library scope (2026-07-05
    author decision -- src.calibrate.SINGLE_CENTRE_ONLY_EXCLUDE, applied
    right after reading the CSV): C2H2/C2H4/C2H6/H2O2/C6H6/iso-C4H10/
    n-C4H10 (two-, four-, or six-centre topologies) are dropped before any
    bond is exploded, matching the same scope restriction now applied to
    src.calibrate.confusion_matrix_stats and the other library-pooled
    figures below.
    """
    _style()
    from src.calibrate import filter_single_centre_library

    lib_df = pd.read_csv(library_csv)
    lib_df = filter_single_centre_library(lib_df)
    bonds = _explode_bonds(lib_df)
    bonds = bonds[bonds["ref_label"].isin(("stretch", "bend")) & bonds["ideal"].isin(("yes", "no"))]
    bonds["abs_rel_db"] = bonds["rel_db"].abs()

    # Width set to ~6.2in (author decision, 2026-07-09) to match the fixed
    # ~6.2in target every main-text single-panel figure prints at once
    # LaTeX's \includegraphics[width=0.95\columnwidth] rescales it on the A4
    # single-column page (see fig:confusion, sized at exactly this width for
    # the same reason) -- at the previous 7.2in source width, the 9pt/7.5pt
    # rcParams fonts were shrinking to an effective ~7.75pt/6.5pt on the
    # printed page, smaller than every other figure in the set. (Previously
    # widened to 7.2in to leave room for manually compositing depicted-mode
    # panel images beside the plot; re-widen if that need resurfaces.)
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    for ideal_flag in ("no", "yes"):  # non-ideal first (background), ideal on top
        for ref in ("bend", "stretch"):
            cat = REF_LABEL_TO_CATEGORY[ref]
            sub = bonds[(bonds["ideal"] == ideal_flag) & (bonds["ref_label"] == ref)]
            if sub.empty:
                continue
            kw = _marker_kwargs(cat, ideal_flag, marker="o", color_map=REF_CATEGORY_COLOR)
            ax.scatter(sub["abs_rel_db"], sub["s_AB"], s=13,
                       zorder=3 if ideal_flag == "yes" else 2, **kw)

    # matplotlib's mathtext DOES support \boldsymbol (verified by rendering,
    # 2026-07-09; superseding an earlier note here that claimed otherwise),
    # so the bond-vector symbol now matches the manuscript prose's own
    # \boldsymbol{b} exactly, dropping the "AB" superscript to match the
    # quantity actually plotted (the per-bond |Delta|b||/|b| ratio, same
    # notation fig:boxplots panel (b) already uses for its averaged version).
    ax.set_xlabel(r"$|\Delta|\boldsymbol{b}|\,/\,|\boldsymbol{b}|\,|$")
    ax.set_ylabel(r"$s^{AB}$")
    ax.set_xlim(-0.03, bonds["abs_rel_db"].max() * 1.05)
    ax.set_ylim(-0.03, 1.05)

    # One consistent marker shape (circle) for every legend entry -- color
    # (stretching/bending) and fill (ideal/non-ideal) are the only two
    # encodings here; shape no longer redundantly re-encodes stretch/bend.
    # Labels use the bare short S/B symbol, matching CATEGORY_LABEL's
    # 2026-07-02 revert (author visual review: an in-figure "(word)" gloss
    # was tried the same day and judged too long -- reverted to bare symbols,
    # gloss deferred to the LaTeX caption) -- this legend is custom-built
    # (not sourced from CATEGORY_LABEL, since it also encodes ideal/
    # non-ideal), so the same bare wording is applied by hand here.
    legend_elems = [
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor=COLORS["stretching_ref"], markeredgecolor=COLORS["stretching_ref"],
               markersize=8, label="S, ideal"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="none", markeredgecolor=COLORS["stretching_ref"],
               markersize=8, label="S, non-ideal"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor=COLORS["bending_ref"], markeredgecolor=COLORS["bending_ref"],
               markersize=8, label="B, ideal"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="none", markeredgecolor=COLORS["bending_ref"],
               markersize=8, label="B, non-ideal"),
    ]
    ax.legend(handles=legend_elems, loc="upper left", frameon=False,
              handletextpad=0.4, labelspacing=0.4, borderaxespad=0.2)

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary = {
        "pdf": pdf_path, "png": png_path,
        "shared_categories": ("reuses REF_CATEGORY_COLOR (old Okabe-Ito "
                               "reference-label hues, not the classification "
                               "algorithm's own colors) for STRETCHING "
                               "(vermillion)/BENDING (blue); ONE marker shape "
                               "(circle) for all points (color already "
                               "distinguishes stretch/bend, so shape is not "
                               "redundantly reused here) and the shared "
                               "IDEAL_STYLE filled=ideal/hollow=non-ideal "
                               "encoding -- same mapping as fig:benzene / "
                               "fig:boxplots / fig:modemixing."),
        "n_bonds": len(bonds),
        "n_molecules": bonds["molecule"].nunique(),
        "x_range": (float(bonds["abs_rel_db"].min()), float(bonds["abs_rel_db"].max())),
        "y_range": (float(bonds["s_AB"].min()), float(bonds["s_AB"].max())),
        "bend_max_s_AB": float(bonds.loc[bonds["ref_label"] == "bend", "s_AB"].max()),
        "ideal_stretch_min_x": float(bonds.loc[(bonds.ideal == "yes") & (bonds.ref_label == "stretch"), "abs_rel_db"].min()),
    }
    return summary


# --------------------------------------------------------------------------
# fig:boxplots -- frequency / Delta|b| / s[V_S] distributions, stretch vs bend
# --------------------------------------------------------------------------

def plot_boxplots(
    library_csv="data/results/library_scores.csv",
    out_dir="data/figures",
    label="fig_boxplots",
):
    """Build fig:boxplots: box plots of (a) frequency, (b) mode-averaged
    |Delta b|/|b|, (c) s[V_S], each split into 4 groups -- bend/ideal,
    stretch/ideal, bend/non-ideal, stretch/non-ideal (filled = ideal, hollow
    = non-ideal, matching the group's earlier undergraduate-report
    convention and this module's shared IDEAL_STYLE) -- over every internal
    row in the library with a literature stretch/bend label.

    Restricted to the single-centre AB_n hydride-library scope (2026-07-05
    author decision -- see src.calibrate.SINGLE_CENTRE_ONLY_EXCLUDE /
    plot_bond_scores' identical note); applied right after reading the CSV,
    before the internal-row selection below.
    """
    _style()
    from src.calibrate import filter_single_centre_library

    lib_df = pd.read_csv(library_csv)
    lib_df = filter_single_centre_library(lib_df)
    internal = lib_df[(lib_df["kind"] == "internal") &
                       lib_df["ref_label"].isin(("stretch", "bend")) &
                       lib_df["ideal"].isin(("yes", "no"))].copy()

    groups = [("bend", "yes"), ("stretch", "yes"), ("bend", "no"), ("stretch", "no")]
    # Bare S/B symbol, matching CATEGORY_LABEL's settled short-notation
    # vocabulary and every other figure's terminology (fig:confusion,
    # fig:bondscores, fig:modemixing, fig:benzene_normal). History: these
    # were originally the bare short forms "bend"/"stretch"; briefly became
    # the bare long forms ("bending"/"stretching"); briefly became glossed
    # short forms ("B (bending)"/"S (stretching)") on 2026-07-02; that gloss
    # was reverted the same day after author visual review judged it too
    # long for a tick label -- settled on the classifier/manuscript's bare
    # "S"/"B" symbol, gloss deferred to the LaTeX caption only.
    # Stacking "(ideal)"/"(non-ideal)" onto every one of the 4 per-panel
    # tick labels (as tried first, pre-2026-07-02) made adjacent 2-line
    # labels visually run together in this narrow a panel (3 panels sharing
    # a ~7.4in figure) regardless of spacing/font tweaks -- switched instead
    # to short primary labels only. A shared "ideal"/"non-ideal" bracket
    # annotation below each pair (tried next) was itself removed 2026-07-09
    # (author decision): that distinction is now stated in the caption text
    # instead of drawn on the raster.
    group_labels = ["B", "S", "B", "S"]
    # Positions: gap 1.3 within a bend/stretch pair, gap 1.6 between the
    # ideal pair (1,2) and non-ideal pair (3,4) -- sized (see
    # IMPLEMENTATION_PLAN.md 2026-07-02 changelog entry) so neither the
    # "bending"/"stretching" tick-label text nor the pair-level "ideal"/
    # "non-ideal" bracket labels below them collide, verified by rendering.
    positions = [1.0, 2.3, 3.9, 5.2]

    # Bare (a)/(b)/(c) labels only (2026-07-23, author request): the
    # descriptive title text ("Frequency", "Bond-length change", "Vibrational
    # score") is dropped -- each panel's own y-axis label already carries
    # that information -- and the panel letter is centered, not left-aligned.
    panels = [
        ("freq", r"Frequency (cm$^{-1}$)", "(a)"),
        ("delta_b_mean", r"Averaged $|\Delta|\boldsymbol{b}|\,/\,|\boldsymbol{b}|\,|$", "(b)"),
        ("V_Stretch", r"$s[\mathrm{V_S}]$", "(c)"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(7.4, 3.3))
    for ax, (col, ylabel, title) in zip(axes, panels):
        data = []
        for ref, ideal_flag in groups:
            vals = internal.loc[(internal.ref_label == ref) & (internal.ideal == ideal_flag), col].dropna().values
            data.append(vals)
        bp = ax.boxplot(data, positions=positions, widths=0.9, showfliers=True,
                         patch_artist=True,
                         flierprops=dict(marker="o", markersize=2.5, alpha=0.5, linewidth=0),
                         medianprops=dict(color="black", linewidth=1.0))
        for (ref, ideal_flag), box in zip(groups, bp["boxes"]):
            cat = REF_LABEL_TO_CATEGORY[ref]
            color = REF_CATEGORY_COLOR[cat]
            style = IDEAL_STYLE[ideal_flag]
            box.set_edgecolor(color)
            box.set_linewidth(1.1)
            box.set_facecolor(color if style["filled"] else "white")
            box.set_alpha(0.9 if style["filled"] else 1.0)
        for whisk in bp["whiskers"]:
            whisk.set_color("#555555")
        for cap in bp["caps"]:
            cap.set_color("#555555")
        for flier, (ref, ideal_flag) in zip(bp["fliers"], groups):
            color = REF_CATEGORY_COLOR[REF_LABEL_TO_CATEGORY[ref]]
            flier.set_markerfacecolor(color if IDEAL_STYLE[ideal_flag]["filled"] else "none")
            flier.set_markeredgecolor(color)

        ax.set_xticks(positions)
        # Horizontal, unrotated (reverted 2026-07-02 alongside the gloss
        # revert above): the 30 deg rotation was only needed to avoid
        # collisions with the longer glossed "S (stretching)"/"B (bending)"
        # labels; bare single-character "S"/"B" labels have no collision
        # risk at this position spacing and read cleaner flat, verified by
        # rendering.
        ax.set_xticklabels(group_labels)
        ax.set_xlim(positions[0] - 0.7, positions[-1] + 0.7)
        ax.set_ylabel(ylabel)
        ax.set_title(title, loc="center", fontweight="bold")

    # tau_S/tau_B are no longer drawn on panel (c) (author decision,
    # 2026-07-09: threshold lines removed as visual clutter -- the values
    # are still returned below for the caption/summary). The pair-level
    # "ideal"/"non-ideal" bracket annotation below each panel's ticks was
    # also removed the same session (deferred to the caption text instead);
    # positions/group_labels alone are enough to read each box.
    th = Thresholds.calibrated()

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary = {
        "pdf": pdf_path, "png": png_path,
        "shared_categories": ("reuses REF_CATEGORY_COLOR (old Okabe-Ito "
                               "reference-label hues) for STRETCHING/BENDING "
                               "(box edge/fill color) and the shared "
                               "IDEAL_STYLE filled=ideal/hollow=non-ideal "
                               "encoding -- same mapping as fig:benzene / "
                               "fig:bondscores / fig:modemixing."),
        "n_modes": len(internal),
        "group_counts": {f"{r}/{i}": int(((internal.ref_label == r) & (internal.ideal == i)).sum())
                          for r, i in groups},
        "freq_range": (float(internal["freq"].min()), float(internal["freq"].max())),
        "delta_b_range": (float(internal["delta_b_mean"].min()), float(internal["delta_b_mean"].max())),
        "vscore_range": (float(internal["V_Stretch"].min()), float(internal["V_Stretch"].max())),
        "tau_S": th.tau_S, "tau_B": th.tau_B,
    }
    return summary


# --------------------------------------------------------------------------
# fig:modemixing -- mode score vs averaged Delta-bond-length, ideal step vs
# non-ideal gradient
# --------------------------------------------------------------------------

def plot_mode_mixing(
    library_csv="data/results/library_scores.csv",
    out_dir="data/figures",
    label="fig_modemixing",
):
    """Build fig:modemixing: (a) ideal molecules -- V_Stretch vs.
    mode-averaged |Delta b|/|b|| shows a clean step function; (b) non-ideal
    molecules -- the same axes show a graded transition. Both panels reuse
    the shared STRETCHING/BENDING category colors; ideal (a) is all filled,
    non-ideal (b) is all hollow, per the shared IDEAL_STYLE. All points use
    ONE consistent marker shape (circle) regardless of stretch/bend -- color
    already fully distinguishes that split, so shape is not redundantly
    reused here (fixed 2026-07-02 consistency pass).

    NOTE (pending gap, IMPLEMENTATION_PLAN.md / figure-builder standing
    report): this renders only the two-panel ideal-step-vs-non-ideal-
    gradient content the current .tex caption (fig:modemixing) describes.
    The irrep-degeneracy sub-panels (scoped-doc Fig 3c-d in the group's
    earlier report, keyed on trigonal-planar/bent-AB2 irreps there, prose
    mentions ethane here) are NOT built -- the intended molecule/panel form
    for THIS manuscript is still unconfirmed with lead-author/tex-data-sync,
    so nothing is invented for it.

    Restricted to the single-centre AB_n hydride-library scope (2026-07-05
    author decision -- see src.calibrate.SINGLE_CENTRE_ONLY_EXCLUDE /
    plot_bond_scores' identical note); applied right after reading the CSV,
    before the internal-row selection below.
    """
    _style()
    from src.calibrate import filter_single_centre_library

    lib_df = pd.read_csv(library_csv)
    lib_df = filter_single_centre_library(lib_df)
    internal = lib_df[(lib_df["kind"] == "internal") &
                       lib_df["ref_label"].isin(("stretch", "bend")) &
                       lib_df["ideal"].isin(("yes", "no"))].copy()

    th = Thresholds.calibrated()
    fig, (ax_i, ax_n) = plt.subplots(1, 2, figsize=(7.2, 3.3), sharex=True, sharey=True)

    # Bare (a)/(b) labels only (2026-07-23, author request): "Ideal
    # molecules"/"Non-ideal molecules" is dropped from the in-figure title
    # (it's stated in the caption) and the panel letter is centered.
    for ax, ideal_flag, title in ((ax_i, "yes", "(a)"),
                                   (ax_n, "no", "(b)")):
        for ref in ("bend", "stretch"):
            cat = REF_LABEL_TO_CATEGORY[ref]
            sub = internal[(internal.ideal == ideal_flag) & (internal.ref_label == ref)]
            if sub.empty:
                continue
            kw = _marker_kwargs(cat, ideal_flag, marker="o", color_map=REF_CATEGORY_COLOR)
            # Each panel's own point label spells out ideal/non-ideal
            # explicitly (2026-07-23, author request) rather than leaving it
            # to be inferred from the panel letter alone.
            ideal_suffix = "ideal" if ideal_flag == "yes" else "non-ideal"
            point_label = f"{CATEGORY_LABEL[cat]}, {ideal_suffix}"
            ax.scatter(sub["delta_b_mean"], sub["V_Stretch"], s=20,
                       label=point_label, **kw)
        ax.axhline(th.tau_S, color=COLORS["threshold"], ls="--", lw=0.8, zorder=1)
        ax.axhline(th.tau_B, color=COLORS["threshold"], ls="--", lw=0.8, zorder=1)
        ax.set_xlabel(r"Averaged $|\Delta|\boldsymbol{b}|\,/\,|\boldsymbol{b}|\,|$")
        ax.set_title(title, loc="center", fontweight="bold")
        ax.set_xlim(-0.03, internal["delta_b_mean"].max() * 1.08)

    ax_i.set_ylabel(r"$s[\mathrm{V_S}]$")
    ax_i.set_ylim(-0.05, 1.08)
    # Bordered legends (2026-07-23, author request): both panels' legend
    # boxes sit directly on top of scattered data points, so a frame (plus an
    # opaque white face) is needed to keep the swatches/text from reading as
    # more data -- every OTHER legend in this module stays frameon=False
    # since none of them overlap data the same way.
    ax_i.legend(loc="center right", frameon=True, edgecolor="black",
                facecolor="white", framealpha=1.0,
                handletextpad=0.4, labelspacing=0.4)
    ax_n.legend(loc="center right", frameon=True, edgecolor="black",
                facecolor="white", framealpha=1.0,
                handletextpad=0.4, labelspacing=0.4)
    ax_i.text(0.98, th.tau_S, r"$\tau_\text{S}$", ha="right", va="bottom",
              color=COLORS["threshold"], transform=ax_i.get_yaxis_transform())
    ax_i.text(0.98, th.tau_B, r"$\tau_\text{B}$", ha="right", va="top",
              color=COLORS["threshold"], transform=ax_i.get_yaxis_transform())

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    n_ideal = int((internal.ideal == "yes").sum())
    n_nonideal = int((internal.ideal == "no").sum())
    summary = {
        "pdf": pdf_path, "png": png_path,
        "shared_categories": ("reuses REF_CATEGORY_COLOR (old Okabe-Ito "
                               "reference-label hues) for STRETCHING/BENDING; "
                               "ONE marker shape (circle) for all points "
                               "(color already distinguishes stretch/bend) "
                               "and the shared IDEAL_STYLE filled=ideal/"
                               "hollow=non-ideal encoding -- same mapping as "
                               "fig:benzene / fig:bondscores / fig:boxplots."),
        "n_ideal_modes": n_ideal,
        "n_nonideal_modes": n_nonideal,
        "tau_S": th.tau_S, "tau_B": th.tau_B,
        "irrep_degeneracy_panel": ("NOT BUILT in this figure -- see the new, "
                                   "separate SI figure plot_irrep_coupling "
                                   "(fig_irrep_coupling) below, added "
                                   "2026-07-06 once the molecule/panel form "
                                   "was confirmed (AB3 trigonal-planar / AB2 "
                                   "bent series, reproducing the group's "
                                   "earlier report's Figures 3c/3d)."),
    }
    return summary


# --------------------------------------------------------------------------
# SI: irrep-degeneracy mixing mechanism (AB3 trigonal-planar / AB2 bent
# series) -- NEW figure, added 2026-07-06.
#
# Supports the manuscript's own prose claim (Stretching/bending
# classification section, near "Modes whose irreducible representations are
# unique retain ideal scores"): "Mixing occurs only where modes of the same
# irreducible representation can couple, i.e. in lower-symmetry molecules."
# That claim previously had no supporting figure -- this reproduces, with
# the CURRENT pipeline's own recomputed scores (not re-derived here; see
# data source note below), the group's earlier project report's Figures 3c
# (trigonal-planar AB3) and 3d (bent AB2), which made exactly this point
# visually: mode score vs. central-atom displacement amplitude |d_CA|,
# faceted by whether a mode's irreducible representation has a same-irrep
# coupling partner within its own point group.
#
#   - AB3 (D3h): A2'' (bend) and A1' (stretch) are each the ONLY mode of
#     that label -- no coupling partner, so they stay flat/clean regardless
#     of |d_CA|. E' appears twice (once as the doubly-degenerate in-plane
#     bend, once as the doubly-degenerate asymmetric stretch) -- those two
#     E' mode sets CAN couple; only the stretch-E' branch visibly droops
#     from ~1 as |d_CA| grows (heavier/more electronegative substituents),
#     while bend-E' stays near 0 throughout.
#   - AB2 (C2v): A1 appears twice (symmetric bend AND symmetric stretch) --
#     both A1 branches show drooping/mixing with |d_CA| (bend-A1 creeps up
#     from 0, stretch-A1 droops down from 1). B2 (antisymmetric stretch) is
#     the only B2-labeled mode -- unique, and stays comparatively clean.
#
# Data source (REPOINTED 2026-07-08, data_score.csv retirement):
# ``data/characterised_modes.csv`` for `irrep`/`shape`/`type` (regenerated
# from a direct on-disk scan by
# ``src.library_ingest.regenerate_characterised_modes()``),
# ``data/mol_list_method.csv``'s `mol_type` column for the ideal/non-ideal
# filter (replaces the old per-row `ideal` column; naturally excludes
# multi-centre molecules like C6H6, same outcome as before), and
# ``data/results/library_scores.csv`` for `V_Stretch` (-> `vib_scr`) and the
# new `d_CA` column (-> `|d_CA|`, an engine-derived quantity now, not a
# frozen spreadsheet snapshot -- see src/library_ingest.py's
# `_central_atom_index`/`score_geometry_molecule`). `data/data_score.csv`
# itself is no longer read anywhere in this module. Merged on
# (molecule, mode) -- characterised_modes.csv's `mode` and
# library_scores.csv's internal-row `mode_index` share the same 1-based
# internal-mode indexing convention (Gaussian's own frequency-block order).
# --------------------------------------------------------------------------

def _hollow_marker_kwargs(color_key, marker):
    """Marker styling for one (type, irrep) category in fig_irrep_coupling:
    color encodes the literal literature bend/stretch label (COLORS'
    'bending_ref'/'stretching_ref' old Okabe-Ito hue, via `color_key` --
    this is a REFERENCE label, not a classifier prediction, so it uses the
    "_ref" entries rather than the classification algorithm's own "bending"/
    "stretching"); marker SHAPE encodes irrep identity, with the SAME
    shape reused for the same irrep symbol across the bend/stretch color
    split whenever that irrep is shared (coupling-capable) between the two
    categories -- so shape-matching across a color change is itself the
    visual signal of a symmetry-permitted coupling pathway. ALL markers are
    unfilled/hollow (author revision 2026-07-06): fill no longer encodes
    anything here (its old unique-vs-shared role is now carried entirely by
    shape), so every category gets `facecolors="none"` with a color-matched
    edge.
    """
    color = COLORS[color_key]
    return dict(marker=marker, facecolors="none", edgecolors=color,
                linewidths=1.15, alpha=0.95)


def plot_irrep_coupling(
    characterised_modes_csv="data/characterised_modes.csv",
    library_scores_csv="data/results/library_scores.csv",
    mol_list_csv="data/mol_list_method.csv",
    out_dir="data/figures",
    label="fig_irrep_coupling",
):
    """Build the irrep-degeneracy mixing-mechanism SI figure: (a) trigonal-
    planar AB3 series (A=B,Al,Ga; B=H,F,Cl,Br -- ``tab:nonideal``'s "Trig.
    planar AB3" row), (b) bent AB2 series (A=O,S,Se; B=H,F,Cl,Br). Both
    panels plot mode score (``V_Stretch``, from ``library_scores.csv``) vs.
    central-atom displacement amplitude (``d_CA``, also from
    ``library_scores.csv`` -- a genuine engine-derived quantity now, see the
    data-source comment above this function), faceted by (type, irrep): ALL
    markers are unfilled/hollow (author revision 2026-07-06); color is bend
    (blue)/stretch (vermillion), the OLD reference-label hues
    (REF_CATEGORY_COLOR / "_ref" entries -- this is the literal literature
    type, not a classifier prediction, so it deliberately does NOT use the
    classification-algorithm's own "bending"/"stretching" colors), and
    marker SHAPE encodes irrep identity,
    with the SAME shape reused for the same irrep symbol whenever it appears
    in BOTH a bend and a stretch category (the shared/coupling-capable
    case) -- e.g. panel (a)'s E' is a triangle in both "bend E'" and
    "stretch E'"; panel (b)'s A1 is a diamond in both "bend A1" and
    "stretch A1". Irreps unique to one category each get their own distinct
    shape (panel (a): A2″ = circle, A1' = square; panel (b): B2 = circle).
    Shape-matching across the bend/stretch color split is itself the visual
    cue for a symmetry-permitted coupling pathway between those two
    categories.

    Reproduces the group's earlier project report's Figures 3c/3d with the
    current pipeline's own recomputed scores (see data source note in this
    section's header comment above).

    **2026-07-08 repointing note:** previously read ``data/data_score.csv``
    directly (a frozen, no-longer-maintained snapshot). Now reads
    ``irrep``/``shape``/``type`` from the disk-scan-regenerated
    ``characterised_modes.csv``, the ideal/non-ideal filter from
    ``mol_list_method.csv``'s ``mol_type`` column, and ``V_Stretch``/``d_CA``
    from the live-computed ``library_scores.csv``. Population may differ
    from the pre-2026-07-08 figure: several molecules that used to appear
    here (``Cl2O``, ``NO2``, ``SO2``, ``SeO2``, ``TeF2``, ``TeO2``, plus
    ``CCl4``/``CF4``/``CH4``/``GeCl4``/``GeF4``/``GeH4``/``SiCl4``/``SiF4``/
    ``SiH4`` in other shape families) no longer have an on-disk ``.log``/
    ``.gjf`` pair at all (removed in an earlier roster-finalization cleanup)
    and are therefore absent from the regenerated ``characterised_modes.csv``
    -- this is the pre-existing, already-documented ``fig:irrep_coupling``/
    ``data_score.csv``-staleness gap (IMPLEMENTATION_PLAN.md, 2026-07-07
    entry) now surfaced rather than papered over by stale data. Newly-added
    roster molecules with no prior literature ``type``/``irrep`` label simply
    do not plot (blank strings match no category) until hand-labeled -- not a
    regression, the same "purely a missing-ground-truth problem" status as
    the rest of the library.

    **2026-07-08 back-fill note:** 24 T-shaped/see-saw/stragglers (incl.
    ``BBr3``, ``OCl2``) gained literature labels this session (see
    IMPLEMENTATION_PLAN.md's back-fill entry) but only ``OCl2`` (bent AB2)
    actually contributes points here -- ``BBr3`` (labeled ``shape="trigonal
    pyramidal"`` with ASCII-digit irreps ``"A1'"``/``'A2"'``) matches neither
    this panel's ``shape=="trigonal planar"`` filter nor its Unicode-subscript
    irrep categories (``"A₁'"``/``"A₂\""``), despite its irrep *pattern*
    (E'/A1'/A2") being the D3h signature every other AB3 sibling in this
    figure (BF3, BCl3, AlCl3, AlBr3, GaCl3, ...) already uses -- almost
    certainly a data-entry inconsistency in the manual back-fill, flagged for
    the author to verify, NOT silently corrected here. ``OCl2`` itself has
    the same ASCII-vs-Unicode-subscript irrep mismatch (``"A1"``/``"B2"``
    instead of its Br2O/OF2 siblings' ``"A₁"``/``"B₂"``) and so ALSO plots
    zero points despite passing the ``shape=="bend"`` eligibility check --
    the returned summary's ``ab3_shape_eligible_but_not_plotted``/
    ``ab2_shape_eligible_but_not_plotted`` keys surface exactly this
    eligible-but-unrendered case (added this session; ``ab3_molecules``/
    ``ab2_molecules`` now report only molecules with >=1 rendered point,
    where they previously reported shape-eligibility only and could
    misleadingly list a molecule that drew nothing, as OCl2 did).

    Not one of the 6 originally-scoped ``fig:*`` labels -- a new SI figure
    filling the "irrep-degeneracy sub-panel" pending gap flagged against
    ``fig:modemixing`` (IMPLEMENTATION_PLAN.md); confirmed molecule/panel
    form (AB3/AB2, not ethane) 2026-07-06.
    """
    _style()
    cm = pd.read_csv(characterised_modes_csv, dtype={"mode": "Int64"})
    lib = pd.read_csv(library_scores_csv)
    roster = pd.read_csv(mol_list_csv)

    mol_type = dict(zip(roster["molecule"], roster["mol_type"]))
    nonideal_molecules = {m for m, t in mol_type.items() if t == "non-ideal"}

    internal = lib.loc[lib["kind"] == "internal",
                        ["molecule", "mode_index", "V_Stretch", "d_CA"]].copy()
    internal = internal.rename(columns={
        "mode_index": "mode", "V_Stretch": "vib_scr", "d_CA": "|d_CA|"})
    internal["mode"] = internal["mode"].astype("Int64")

    df = cm.merge(internal, on=["molecule", "mode"], how="inner")
    df = df[df["molecule"].isin(nonideal_molecules)]
    df = df.dropna(subset=["|d_CA|"])  # multi-centre/no-central-atom rows

    ab3 = df[df["shape"] == "trigonal planar"].copy()
    ab2 = df[df["shape"] == "bend"].copy()

    # (type, irrep) -> (color_key, marker shape, legend text). Defined
    # per-panel (not a single shared dict) since the two point groups use
    # different irrep labels and different shared/unique assignments (AB2's
    # A1 is shared in BOTH branches; AB3's A1'/A2'' are each unique to one
    # branch while E' is shared). The marker shape is keyed to the irrep
    # SYMBOL, not to (type, irrep), so the same irrep gets the same shape
    # in both its bend and stretch appearances -- circle/triangle/square/
    # diamond chosen for clear distinction even hollow and small; no
    # plus/x (easily lost against gridlines).
    # "type" here is the literal literature/reference label (not a
    # classifier prediction) -- color_key uses the "_ref" old-blue/vermillion
    # entries, not the classification-algorithm's own "bending"/"stretching".
    ab3_categories = [
        ("bend", "A₂\"", "bending_ref", "o", "bend A2″ (unique irrep)"),
        ("bend", "E'", "bending_ref", "^", "bend E′ (shared irrep -- same shape as stretch E′)"),
        ("stretch", "A₁'", "stretching_ref", "s", "stretch A1′ (unique irrep)"),
        ("stretch", "E'", "stretching_ref", "^", "stretch E′ (shared irrep -- same shape as bend E′)"),
    ]
    ab2_categories = [
        ("bend", "A₁", "bending_ref", "D", "bend A1 (shared irrep -- same shape as stretch A1)"),
        ("stretch", "A₁", "stretching_ref", "D", "stretch A1 (shared irrep -- same shape as bend A1)"),
        ("stretch", "B₂", "stretching_ref", "o", "stretch B2 (unique irrep)"),
    ]

    # Widened 7.2->9.0in (2026-07-23, alongside the font-size bump): at the
    # larger rcParams sizes, the per-panel legend text and the shared x-axis
    # label no longer fit inside the original 7.2in width -- the legend
    # overflowed into the neighboring panel and the x-label clipped at the
    # canvas edge (verified by rendering). This is a layout fix, not a font
    # override -- the extra width is what these labels now need at the
    # default size, matching every OTHER figure's rule of leaving fontsize
    # alone (LEGEND_FONTSIZE is still applied on this SI figure's legend
    # only because its entries are unusually long irrep-coupling captions).
    fig, (ax_3, ax_2) = plt.subplots(1, 2, figsize=(9.0, 3.4))

    def _plot_panel(ax, sub_df, categories, title):
        n_plotted = 0
        plotted_molecules = set()
        for mode_type, irrep, color_key, marker, legend_text in categories:
            pts = sub_df[(sub_df["type"] == mode_type) & (sub_df["irrep"] == irrep)]
            if pts.empty:
                continue
            kw = _hollow_marker_kwargs(color_key, marker)
            ax.scatter(pts["|d_CA|"], pts["vib_scr"], s=34, zorder=3,
                       label=legend_text, **kw)
            n_plotted += len(pts)
            plotted_molecules.update(pts["molecule"].unique().tolist())
        ax.set_xlabel(r"$|\mathbf{d}_{CA}|$ (central-atom displacement amplitude)")
        ax.set_title(title, loc="center", fontweight="bold")
        ax.set_xlim(-0.03, sub_df["|d_CA|"].max() * 1.08)
        ax.set_ylim(-0.05, 1.08)
        ax.legend(loc="center left", frameon=False, fontsize=LEGEND_FONTSIZE,
                   handletextpad=0.4, labelspacing=0.4, borderaxespad=0.2)
        return n_plotted, plotted_molecules

    # Bare (a)/(b) labels only (2026-07-23, author request): "Trigonal-planar
    # AB3 (non-ideal)"/"Bent AB2 (non-ideal)" is dropped from the in-figure
    # title (stated in the caption instead) and the panel letter is centered.
    n_ab3, ab3_plotted = _plot_panel(ax_3, ab3, ab3_categories, "(a)")
    n_ab2, ab2_plotted = _plot_panel(ax_2, ab2, ab2_categories, "(b)")
    ax_3.set_ylabel(r"Mode score $s[\mathrm{V}]$")

    # Molecules whose `shape` matched the panel but that contributed ZERO
    # actually-rendered points because their `type`/`irrep` strings didn't
    # match any (mode_type, irrep) category above (e.g. an ASCII "A1"/"B2"
    # irrep where the category list expects the Unicode-subscript "A₁"/"B₂"
    # used by every other sibling row, or a `shape` value that is itself
    # inconsistent with the family, e.g. BBr3 -- see IMPLEMENTATION_PLAN.md's
    # 2026-07-08 back-fill entry). Reported explicitly rather than silently
    # folded into "ab3_molecules"/"ab2_molecules" below, since those two keys
    # used to report shape-eligibility, not actual rendering, which was
    # misleading (e.g. OCl2 previously appeared in "ab2_molecules" with 0 of
    # its 3 modes ever drawn).
    ab3_shape_only = sorted(set(ab3["molecule"].unique()) - ab3_plotted)
    ab2_shape_only = sorted(set(ab2["molecule"].unique()) - ab2_plotted)

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary = {
        "pdf": pdf_path, "png": png_path,
        "shared_categories": ("reuses CATEGORY_COLOR's bending/stretching "
                               "hues; ALL markers unfilled/hollow (2026-07-06 "
                               "revision -- fill no longer encodes anything "
                               "here). Marker SHAPE now encodes irrep identity: "
                               "panel (a) A2″=circle, E'=triangle (shared "
                               "between bend/stretch), A1'=square; panel (b) "
                               "A1=diamond (shared between bend/stretch), "
                               "B2=circle. Same irrep symbol always gets the "
                               "same shape across the bend/stretch color "
                               "split, so shape-matching visually flags a "
                               "coupling-capable shared irrep."),
        "ab3_molecules": sorted(ab3_plotted),
        "ab2_molecules": sorted(ab2_plotted),
        "ab3_shape_eligible_but_not_plotted": ab3_shape_only,
        "ab2_shape_eligible_but_not_plotted": ab2_shape_only,
        "n_ab3_points": n_ab3,
        "n_ab2_points": n_ab2,
        "data_source": ("data/characterised_modes.csv (irrep/shape/type) + "
                         "data/mol_list_method.csv (mol_type ideal filter) + "
                         "data/results/library_scores.csv (V_Stretch/d_CA) -- "
                         "not data/data_score.csv (retired 2026-07-08)"),
    }
    return summary


# --------------------------------------------------------------------------
# fig:sensitivity -- tau_TR label-change fraction & accuracy, plateau
# --------------------------------------------------------------------------

def plot_sensitivity(
    sweep_csv="data/results/tau_sensitivity_sweep.csv",
    thresholds_json="data/results/thresholds.json",
    out_dir="data/figures",
    label="fig_sensitivity",
):
    """Build fig:sensitivity: label-change fraction and accuracy vs. tau_TR
    over the calibration grid (src.calibrate.sweep_tau_tr's persisted
    output), marking the frozen tau_TR=0.95 and shading the plateau range
    from thresholds.json's `plateau_tau_TR_range`.

    Scope note: only tau_TR is swept in the persisted data (tau_S/tau_B are
    derived analytically from the ideal-molecule gap, not grid-swept -- see
    src/calibrate.py); this figure renders the tau_TR sensitivity only, per
    the available data (no tau_S/tau_B sweep is fabricated).
    """
    _style()
    sweep = pd.read_csv(sweep_csv)
    with open(thresholds_json) as f:
        result = json.load(f)
    tau_TR_frozen = result["tau_TR"]
    lo, hi = result["plateau_tau_TR_range"]

    # Enlarged 4.6x3.4 -> 5.6x4.6in (2026-07-23, alongside the font-size
    # bump): both y-axis labels are long strings rendered ROTATED 90 deg, so
    # what they need is figure HEIGHT, not width -- a first attempt that only
    # widened the figure (to 5.6x3.4) left them clipping at the top/bottom of
    # the canvas exactly as before (verified by rendering); tight_layout does
    # not reliably reserve space for a twinx() secondary label either way.
    # Layout fix, not a font override.
    fig, ax1 = plt.subplots(figsize=(5.6, 4.6))
    ax1.axvspan(lo, hi, color=COLORS["plateau_band"], alpha=0.18, lw=0, zorder=0,
                label=f"plateau [{lo:g}, {hi:g}]")
    ax1.axvline(tau_TR_frozen, color=COLORS["threshold"], ls="--", lw=1.0, zorder=2,
                label=r"frozen $\tau_{\mathrm{TR}}$" + f"={tau_TR_frozen:g}")

    l1, = ax1.plot(sweep["tau_TR"], sweep["label_change_fraction"],
                   color=COLORS["sens_change"], marker="o", markersize=2.5,
                   lw=1.1, zorder=3, label="label-change fraction")
    ax1.set_xlabel(r"$\tau_{\mathrm{TR}}$")
    ax1.set_ylabel("Label-change fraction (vs. previous grid point)")
    ax1.set_ylim(-0.01, max(0.05, sweep["label_change_fraction"].max() * 1.2))

    ax2 = ax1.twinx()
    l2, = ax2.plot(sweep["tau_TR"], sweep["accuracy"], color=COLORS["sens_accuracy"],
                   marker="s", markersize=2.5, lw=1.1, zorder=3, label="accuracy")
    ax2.set_ylabel("Accuracy (normal-mode T/R ground truth)")
    ax2.set_ylim(min(0.95, sweep["accuracy"].min() - 0.02), 1.01)

    handles = [l1, l2,
               Patch(facecolor=COLORS["plateau_band"], alpha=0.18, label=f"plateau [{lo:g}, {hi:g}]"),
               Line2D([0], [0], color=COLORS["threshold"], ls="--",
                      label=r"frozen $\tau_{\mathrm{TR}}$" + f"={tau_TR_frozen:g}")]
    ax1.legend(handles=handles, loc="center left", frameon=False,
               handletextpad=0.5, labelspacing=0.4)

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary = {
        "pdf": pdf_path, "png": png_path,
        "shared_categories": ("no classification-category colors here (tau_TR "
                               "curve, not per-mode points); reuses the shared "
                               "COLORS['threshold'] dashed-line convention "
                               "from fig:benzene/fig:boxplots/fig:modemixing "
                               "for the frozen-tau marker."),
        "n_grid_points": len(sweep),
        "tau_TR_range": (float(sweep["tau_TR"].min()), float(sweep["tau_TR"].max())),
        "accuracy_range": (float(sweep["accuracy"].min()), float(sweep["accuracy"].max())),
        "label_change_fraction_range": (float(sweep["label_change_fraction"].min()),
                                         float(sweep["label_change_fraction"].max())),
        "tau_TR_frozen": tau_TR_frozen,
        "plateau_range": (lo, hi),
    }
    return summary


def regenerate_all(verbose=True):
    """Regenerate every manuscript figure in one call (factored out of the
    former ``if __name__ == "__main__":`` block so callers -- e.g.
    ``main.py --figures`` -- can invoke this directly instead of shelling out
    to ``python -m src.figures``). Each figure function reads its own
    already-computed ``data/results/*.csv`` inputs with their own defaults;
    this function takes no molecule-specific arguments.

    Returns a dict {tex_label: result_dict} for all 12 figures (6 labeled
    fig:* + the standalone benzene-normal-modes gallery + the SI confusion
    retention/migration panel + the SI rigorous-tier consistency check +
    fig:benzeneconfusion + its SI precision/recall companion + the SI
    irrep-degeneracy coupling figure), in the same order they are built.
    FIXED 2026-07-23: plot_confusion_retention_migration was missing from
    this list entirely (an orchestrator gap, unrelated to any dataset
    change) -- it now runs like every other figure here instead of only
    when called directly. Raises whatever the underlying plot_* function
    raises (e.g. a missing input CSV) -- fail loud, no silent partial
    regeneration.
    """
    fns = [
        ("fig:benzene", plot_benzene_stress_test),
        ("benzene-normal-modes gallery (no fig: label yet)", plot_benzene_normal_modes),
        ("fig:confusion", plot_confusion_matrix),
        ("SI confusion retention/migration (no fig: label yet)", plot_confusion_retention_migration),
        ("SI rigorous-tier consistency check (no fig: label yet)", plot_rigorous_tier_check),
        ("fig:benzeneconfusion", plot_benzene_internal_confusion),
        ("SI benzene precision/recall (no fig: label yet)", plot_benzene_confusion_precision_recall),
        ("fig:bondscores", plot_bond_scores),
        ("fig:boxplots", plot_boxplots),
        ("fig:modemixing", plot_mode_mixing),
        ("SI irrep-degeneracy coupling (no fig: label yet)", plot_irrep_coupling),
        ("fig:sensitivity", plot_sensitivity),
    ]
    results = {}
    for tex_label, fn in fns:
        result = fn()
        results[tex_label] = result
        if verbose:
            print(f"{tex_label} ->", result["pdf"])
            print(f"{tex_label} ->", result["png"])
            for k, v in result.items():
                if k not in ("pdf", "png"):
                    print(f"  {k}: {v}")
    return results


if __name__ == "__main__":
    regenerate_all()
