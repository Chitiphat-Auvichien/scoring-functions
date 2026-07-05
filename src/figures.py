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
predicted) internal-mode confusion matrix + precision/recall, distinct from
the whole-hydride-library ``fig:confusion``.

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
)

# --------------------------------------------------------------------------
# Shared house style
# --------------------------------------------------------------------------

# Okabe-Ito colorblind-safe palette.
COLORS = {
    "external": "#999999",   # gray      -- CLEAN_TRANSLATION / CLEAN_ROTATION
    "bending": "#0072B2",    # blue      -- BENDING
    "stretching": "#D55E00", # vermillion-- STRETCHING
    "mixed": "#009E73",      # teal      -- MIXED_STRETCH_BEND
    "mixed_ext": "#CC79A7",  # purple    -- MIXED_EXTERNAL_WITH_VIBRATION
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
    "confusion_cmap": "Blues",    # sequential, colorblind-safe -- fig:confusion heatmap
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
}

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


def _style():
    """Apply the shared print-ready matplotlib style. Call once per figure."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7.5,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "lines.linewidth": 1.0,
        "lines.markersize": 5,
        "figure.dpi": 150,
        "savefig.dpi": 300,
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
    fig.savefig(png_path, bbox_inches="tight", dpi=300)
    for p in (pdf_path, png_path):
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            raise RuntimeError(f"figure save failed or produced an empty file: {p}")
    return pdf_path, png_path


def _marker_kwargs(category, ideal_flag=None, marker=None):
    """Shared per-point marker styling for a classification `category`
    (key into CATEGORY_COLOR/CATEGORY_MARKER), optionally faceted by
    `ideal_flag` ('yes'/'no'/None) via the shared IDEAL_STYLE encoding
    (filled = ideal, hollow = non-ideal). Returns a dict ready to splat into
    ax.scatter(...).

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
    color = CATEGORY_COLOR[category]
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
    joined `s_AB`/`rel_db` internal-row strings (see excel_ingest.py's
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
                  xytext=(6, 6), textcoords="offset points", fontsize=7)
    ax_b.annotate("EMIT 9", (e9["C2_Ry"], abs(e9["Ry"])),
                  xytext=(-6, 8), textcoords="offset points", fontsize=7,
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
    ax_b.text(callout_xy[0], callout_xy[1], callout_text, fontsize=6.7,
              ha="center", va="center",
              bbox=dict(boxstyle="round,pad=0.35", fc="white",
                         ec=COLORS["highlight_t"], lw=0.8), zorder=6)

    ax_b.set_xlabel(r"Projected normal-mode contribution, $\tilde{\Theta}_i^2$")
    ax_b.set_ylabel(r"$|\,\mathrm{score}\,|$ (this framework)")
    ax_b.set_xlim(-0.03, 1.05)
    ax_b.set_ylim(-0.03, 1.15)
    ax_b.legend(loc="upper left", frameon=False, handletextpad=0.3,
                labelspacing=0.35, borderaxespad=0.1, fontsize=6.7)
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

    fig, ax = plt.subplots(figsize=(3.8, 3.6))

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
    ax.text(0.02, TAU_S + y_offset, r"$\tau_S=$" + f"{TAU_S:.3f}", ha="left",
            va="bottom", fontsize=7, color=COLORS["threshold"],
            transform=ax.get_yaxis_transform())
    ax.text(0.98, TAU_B - y_offset, r"$\tau_B=$" + f"{TAU_B:.3f}", ha="right",
            va="top", fontsize=7, color=COLORS["threshold"],
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
              handletextpad=0.3, labelspacing=0.3, borderaxespad=0.2, fontsize=7)

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
# fig:confusion -- two-tier clean-category confusion matrix + precision/
# recall (rigorous ground truth) / retention-and-migration (non-ideal
# validation-by-characterization)
# --------------------------------------------------------------------------

def _confusion_heatmap(ax, fig, tbl, ref_order, title, label_map=CATEGORY_LABEL):
    """Shared heatmap renderer for one confusion-matrix tier. `tbl` must
    already be reindexed to `ref_order` rows (columns are whatever buckets
    are present for that tier -- rigorous and non-ideal tiers populate
    different bucket sets, so no forced union of columns across tiers).
    Returns the imshow handle (caller attaches its own colorbar).

    `label_map` defaults to the shared `CATEGORY_LABEL` dict but
    `plot_confusion_matrix` passes a figure-local override (see that
    function's `_CONFUSION_LABEL`) so this figure's translation/rotation
    ticks read "T"/"R" instead of the long "clean translation"/"clean
    rotation" form, without touching `CATEGORY_LABEL` itself.
    """
    vmax = max(1, tbl.values.max())
    im = ax.imshow(tbl.values, cmap=COLORS["confusion_cmap"], aspect="auto",
                    vmin=0, vmax=vmax)
    for i in range(tbl.shape[0]):
        for j in range(tbl.shape[1]):
            v = int(tbl.values[i, j])
            txt_color = "white" if v > 0.6 * vmax else "black"
            ax.text(j, i, str(v), ha="center", va="center", fontsize=8,
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
                        rotation=0, ha="center", fontsize=9)
    ax.set_yticks(range(len(tbl.index)))
    ax.set_yticklabels([label_map[REF_LABEL_TO_CATEGORY[r]] for r in tbl.index],
                        fontsize=9)
    for tick, r in zip(ax.get_yticklabels(), tbl.index):
        tick.set_color(CATEGORY_COLOR[REF_LABEL_TO_CATEGORY[r]])
    for tick, c in zip(ax.get_xticklabels(), tbl.columns):
        cat = PRED_BUCKET_TO_CATEGORY.get(c)
        if cat is not None:
            tick.set_color(CATEGORY_COLOR[cat])
    ax.set_xlabel("Predicted bucket")
    ax.set_ylabel("Reference label")
    ax.set_title(title, loc="left", fontweight="bold", fontsize=9)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("n modes", fontsize=8)
    cb.ax.tick_params(labelsize=7.5)
    return im


def plot_confusion_matrix(
    library_csv="data/results/library_scores.csv",
    out_dir="data/figures",
    label="fig_confusion",
):
    """Build fig:confusion as a TWO-TIER figure (rebuilt 2026-07-02; replaces
    the earlier single pooled-matrix version -- see IMPLEMENTATION_PLAN.md
    Changelog). The library's literature ``ref_label`` is exact group-theory
    ground truth only for ``ideal == 'yes'`` internal modes (and for every
    external T/R row, geometry-exact regardless of the ``ideal`` tag);
    pooling ``ideal == 'no'`` internal modes into the same accuracy number
    risks reading intrinsic stretch/bend mixing (the effect this framework is
    built to detect, B8.3 CoM-softening) as classifier error. So this
    function computes ``src.calibrate.confusion_matrix_stats`` TWICE, once
    per tier, via the same ad-hoc ``ideal``-column filter the manuscript
    prose already uses (no change to calibrate.py):

      * Rigorous tier: every external (T/R) row + internal rows with
        ``ideal == 'yes'`` -- true accuracy claim.
      * Non-ideal tier: internal rows with ``ideal == 'no'`` only --
        reframed as validation-by-characterization (label retention /
        migration-to-mixed), not an accuracy claim.

    Layout (2x2): (a) rigorous confusion matrix [top-left], (b) rigorous
    per-category precision/recall bars, all == 1.000 [top-right],
    (c) non-ideal confusion matrix, bend/stretch reference x
    bend/mixed/stretch predicted bucket [bottom-left], (d) non-ideal
    label-retention vs. migration-to-mixed bars, annotated with the 0%
    opposite-clean-category-crossing finding [bottom-right].

    Never recomputes scores -- only ``confusion_matrix_stats`` (already
    computed from calibrated thresholds) is called, on two filtered slices
    of the same ``library_scores.csv`` this module always reads.
    ``confusion_matrix_stats`` itself additionally restricts both slices to
    the single-centre AB_n hydride-library scope (2026-07-05 author
    decision, ``src.calibrate.SINGLE_CENTRE_ONLY_EXCLUDE`` -- excludes
    C2H2/C2H4/C2H6/H2O2/C6H6/iso-C4H10/n-C4H10) before computing anything,
    so this function does not need its own copy of that filter.
    """
    _style()
    from src.calibrate import confusion_matrix_stats

    # Local-only override, THIS FIGURE ONLY (author flagged 2026-07-02: the
    # spelled-out "clean translation"/"clean rotation" tick text was too
    # long, forcing rotated labels and a generally cramped 2x2 layout).
    # Mirrors `plot_benzene_normal_modes`'s `_LEGEND_MERGE_TEXT` local-dict
    # precedent above -- `CATEGORY_LABEL` itself is untouched (other
    # figures/legends still want the fuller form, or a different merge, for
    # translation/rotation).
    _CONFUSION_LABEL = dict(CATEGORY_LABEL)
    _CONFUSION_LABEL["translation"] = "T"
    _CONFUSION_LABEL["rotation"] = "R"

    lib_df = pd.read_csv(library_csv)
    thresholds = Thresholds.calibrated()

    rigorous_df = lib_df[(lib_df["kind"] == "external") | (lib_df["ideal"] == "yes")]
    nonideal_df = lib_df[(lib_df["kind"] == "internal") & (lib_df["ideal"] == "no")]

    stats_r = confusion_matrix_stats(rigorous_df, thresholds)
    stats_n = confusion_matrix_stats(nonideal_df, thresholds)

    fig, ((ax_ha, ax_pa), (ax_hb, ax_pb)) = plt.subplots(
        2, 2, figsize=(7.4, 6.6), gridspec_kw={"width_ratios": [1.15, 1.0]})

    # ================= Tier 1: rigorous (n=231) =================
    # (was n=237 before the 2026-07-05 single-centre-only scope filter --
    # confusion_matrix_stats() now drops benzene's 6 geometry-backed
    # external T/R rows via src.calibrate.filter_single_centre_library,
    # since C6H6 is a six-centre molecule outside the hydride library's
    # AB_n scope and has its own dedicated confusion matrix elsewhere,
    # fig:benzeneconfusion; 237-6=231, computed not assumed.)
    ref_order_r = ["translation", "rotation", "stretch", "bend"]
    table_r = stats_r["confusion_table"]
    pred_order_r = ["translation", "rotation", "stretch", "bend"]
    pred_cols_r = [c for c in pred_order_r if c in table_r.columns] + \
                  [c for c in table_r.columns if c not in pred_order_r]
    tbl_r = table_r.reindex(index=ref_order_r, columns=pred_cols_r, fill_value=0)
    n_rigorous = int(tbl_r.values.sum())
    _confusion_heatmap(ax_ha, fig, tbl_r, ref_order_r,
                        f"(a) Rigorous ground truth (n={n_rigorous})",
                        label_map=_CONFUSION_LABEL)

    cats_r = ["translation", "rotation", "stretch", "bend"]
    x_r = np.arange(len(cats_r))
    width = 0.35
    precisions_r = [stats_r["per_category"][c]["precision"] for c in cats_r]
    recalls_r = [stats_r["per_category"][c]["recall"] for c in cats_r]
    bar_colors_r = [CATEGORY_COLOR[REF_LABEL_TO_CATEGORY[c]] for c in cats_r]
    bars_p = ax_pa.bar(x_r - width / 2, precisions_r, width, color=bar_colors_r,
                        edgecolor="black", linewidth=0.5, label="precision")
    bars_r = ax_pa.bar(x_r + width / 2, recalls_r, width, color=bar_colors_r,
                        edgecolor="black", linewidth=0.5, hatch="///", label="recall")
    # Per-category "1.000" value labels were dropped (2026-07-02 consistency
    # pass): every category's precision AND recall is exactly 1.000 here --
    # that number is already stated once in the panel title ("all = 1.000"),
    # so 4 identical repeated labels directly above the bars were pure
    # redundant clutter, and their claimed vertical space is exactly what the
    # legend needs (see below) -- removing them fixes both problems at once.
    ax_pa.axhline(stats_r["acceptance_floor"], color=COLORS["threshold"], ls="--", lw=0.8)
    ax_pa.text(len(cats_r) - 0.5, stats_r["acceptance_floor"], "  floor",
               va="bottom", ha="right", fontsize=7, color=COLORS["threshold"])
    ax_pa.set_xticks(x_r)
    # Rotation/right-alignment removed (2026-07-02 follow-up, same fix as
    # `_confusion_heatmap` above): now that `_CONFUSION_LABEL` shortens
    # translation/rotation to "T"/"R" (matching the already-short "S"/"B"),
    # these 4 tick labels fit flat with no collision risk, verified by
    # rendering.
    ax_pa.set_xticklabels([_CONFUSION_LABEL[REF_LABEL_TO_CATEGORY[c]] for c in cats_r],
                          rotation=0, ha="center", fontsize=9)
    ax_pa.set_ylim(0, 1.18)
    ax_pa.set_ylabel("Precision / recall")
    ax_pa.set_title("(b) Rigorous precision/recall (all = 1.000)", loc="left",
                    fontweight="bold", fontsize=9)
    # loc="lower left" (older) sat directly above the previously-rotated
    # x-tick labels ("clean translation" etc.), visually crowding/
    # overlapping them. "upper left" stays clear of the (now flat, short)
    # tick labels and any bar annotation -- verified by rendering.
    ax_pa.legend(loc="upper left", frameon=False, fontsize=8)

    # ================= Tier 2: non-ideal (n=277) =================
    # (was n=422 before the 2026-07-05 single-centre-only scope filter --
    # confusion_matrix_stats() now drops the 145 non-ideal internal modes
    # belonging to the 7 non-single-centre molecules (C2H2, C2H4, C2H6,
    # H2O2, C6H6, iso-C4H10, n-C4H10), which brings the non-ideal molecule
    # count from 58 down to 51 -- matching an independent hand-count of
    # tab:nonideal's own AB_n grid; 422-145=277, computed not assumed.)
    ref_order_n = ["stretch", "bend"]
    table_n = stats_n["confusion_table"]
    pred_order_n = ["stretch", "bend", "mixed"]
    pred_cols_n = [c for c in pred_order_n if c in table_n.columns] + \
                  [c for c in table_n.columns if c not in pred_order_n]
    tbl_n = table_n.reindex(index=ref_order_n, columns=pred_cols_n, fill_value=0)
    n_nonideal = int(tbl_n.values.sum())
    _confusion_heatmap(ax_hb, fig, tbl_n, ref_order_n,
                        f"(c) Non-ideal characterization (n={n_nonideal})",
                        label_map=_CONFUSION_LABEL)

    cats_n = ["bend", "stretch"]
    retention_n = [stats_n["per_category"][c]["recall"] for c in cats_n]
    migration_n = [stats_n["per_category"][c]["mixed_fraction"] for c in cats_n]
    opposite_n = []  # explicit 0% opposite-clean-category crossing, per tbl_n
    for c in cats_n:
        opp = "stretch" if c == "bend" else "bend"
        n_ref = stats_n["per_category"][c]["n_ref"]
        n_opp = int(tbl_n.loc[c, opp]) if opp in tbl_n.columns else 0
        opposite_n.append(n_opp / n_ref if n_ref else float("nan"))

    # Segment labels are numeric-only (no in-panel legend): the "retained"
    # (bend/stretch-colored) vs. "migrated to mixed" (teal) color coding
    # reuses the SAME CATEGORY_COLOR swatches panel (c)'s tick labels just
    # showed two columns over, so a redundant legend here would only add
    # clutter. Thin segments (e.g. bend's 0.043 migrated slice) get their
    # label placed just ABOVE the bar instead of centered inside it, so text
    # never overflows a segment shorter than the label's own height.
    x_n = np.arange(len(cats_n))
    bar_colors_n = [CATEGORY_COLOR[REF_LABEL_TO_CATEGORY[c]] for c in cats_n]
    ax_pb.bar(x_n, retention_n, 0.5, color=bar_colors_n,
              edgecolor="black", linewidth=0.5)
    ax_pb.bar(x_n, migration_n, 0.5, bottom=retention_n,
              color=COLORS["mixed"], alpha=0.85, edgecolor="black", linewidth=0.5)
    THIN = 0.08  # segments shorter than this (axis fraction) get an outside label
    for xi, ret, mig in zip(x_n, retention_n, migration_n):
        ax_pb.text(xi, ret / 2, f"retained\n{ret:.3f}", ha="center", va="center",
                   fontsize=7, color="white", fontweight="bold")
        if mig >= THIN:
            ax_pb.text(xi, ret + mig / 2, f"mixed\n{mig:.3f}", ha="center",
                       va="center", fontsize=7, color="white", fontweight="bold")
        else:
            ax_pb.text(xi, ret + mig + 0.015, f"mixed: {mig:.3f}", ha="center",
                       va="bottom", fontsize=7, color=COLORS["mixed"],
                       fontweight="bold")
    ax_pb.set_xticks(x_n)
    # cats_n is stretch/bend only (already short "S"/"B"; unaffected by the
    # translation/rotation shortening above) -- `_CONFUSION_LABEL` used here
    # too only for consistency with panels (a)/(b)'s tick-label font size.
    ax_pb.set_xticklabels([_CONFUSION_LABEL[REF_LABEL_TO_CATEGORY[c]] for c in cats_n],
                          fontsize=9)
    ax_pb.set_xlim(-0.55, 1.55)
    ax_pb.set_ylim(0, 1.12)
    ax_pb.set_ylabel("Fraction of reference-labeled modes")
    ax_pb.set_title("(d) Non-ideal retention vs. migration", loc="left",
                    fontweight="bold", fontsize=9)

    # Whole-figure footer sentence: the 0%-opposite-crossing finding applies
    # to BOTH non-ideal categories and is the point of tier 2. FIXED
    # 2026-07-02 (readability): this used to be drawn in-image via
    # `fig.text()` at 6.8pt on this 7.4x6.6in canvas; LaTeX's
    # `\includegraphics[width=0.95\columnwidth]` rescales that canvas to
    # ~6.2in, shrinking the effective size to ~5.7pt -- below a readable
    # floor once printed. Removed from the raster entirely; the exact
    # sentence is computed here and returned in `summary["nonideal_footer_text"]`
    # instead, for lead-author to place in the actual LaTeX
    # `\captionof{figure}{...}` text (typeset at normal caption font size,
    # not shrunk with the image).
    # ASCII "->" (not a unicode arrow) deliberately: this string is meant to
    # be easy to print/copy on any console (a literal U+2192 arrow crashes
    # `print()` under Windows' default cp1252 stdout encoding), and reads
    # fine as-is in a LaTeX caption too.
    footer_text = (
        f"Non-ideal tier (n={n_nonideal}): 0% of bend or stretch reference-labeled "
        f"modes crossed to the OPPOSITE clean category "
        f"(bend->stretch={opposite_n[0]:.1%}, stretch->bend={opposite_n[1]:.1%}); "
        "100% of the non-retained remainder lands in the mixed bucket."
    )
    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary = {
        "pdf": pdf_path, "png": png_path,
        "nonideal_footer_text": footer_text,
        "layout": ("2x2: (a) rigorous confusion matrix, (b) rigorous "
                   "precision/recall bars, (c) non-ideal confusion matrix, "
                   "(d) non-ideal retention/migration bars -- REPLACES the "
                   "earlier single pooled-matrix fig:confusion."),
        "shared_categories": ("reuses CATEGORY_COLOR/CATEGORY_LABEL for "
                               "CLEAN_TRANSLATION/CLEAN_ROTATION/STRETCHING/"
                               "BENDING/MIXED_STRETCH_BEND -- same mapping "
                               "as fig:benzene. FIXED 2026-07-02 (author "
                               "follow-up): translation/rotation tick text "
                               "is shortened to 'T'/'R' via a figure-local "
                               "`_CONFUSION_LABEL` override (CATEGORY_LABEL "
                               "itself unchanged); tick-label rotation "
                               "(no longer needed once labels are short) "
                               "and several font sizes were also restored "
                               "to normal, readable values."),
        "rigorous_n": n_rigorous,
        "rigorous_confusion_table": tbl_r.to_dict(),
        "rigorous_precision": dict(zip(cats_r, precisions_r)),
        "rigorous_recall": dict(zip(cats_r, recalls_r)),
        "acceptance_floor": stats_r["acceptance_floor"],
        "rigorous_floor_met": stats_r["floor_met"],
        "nonideal_n": n_nonideal,
        "nonideal_confusion_table": tbl_n.to_dict(),
        "nonideal_retention": dict(zip(cats_n, retention_n)),
        "nonideal_migration_to_mixed": dict(zip(cats_n, migration_n)),
        "nonideal_opposite_category_crossing": dict(zip(cats_n, opposite_n)),
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
    matrix_csv="data/results/benzene_internal_confusion_matrix.csv",
    summary_csv="data/results/benzene_internal_confusion_summary.csv",
    out_dir="data/figures",
    label="fig_benzene_confusion",
):
    """Build fig:benzeneconfusion: benzene's (C6H6) 3x3 internal-mode
    confusion matrix -- reference bend/stretch/SB (literal literature label,
    modes 21/22, an E1u degenerate pair at 1532.85 cm-1) x predicted
    bend/stretch/mixed bucket -- plus a per-category precision/recall bar
    panel. Data comes from
    ``src.benzene_validation.run_benzene_internal_confusion``'s two output
    CSVs (never recomputed here); this figure is presentation-only, exactly
    like every other function in this module.

    Layout (1x2, mirrors ``plot_confusion_matrix``'s top row): (a) the 3x3
    heatmap via the shared ``_confusion_heatmap`` renderer (reusing
    ``CATEGORY_COLOR``/``CATEGORY_LABEL`` -- the literal "SB" reference row
    is routed to the "mixed" category color/label via the
    ``REF_LABEL_TO_CATEGORY["SB"] = "mixed"`` addition, see that dict's
    comment); (b) grouped precision/recall bars for 3 categories (bend,
    stretch, and "mixed" -- the predicted-bucket name for what is a literal
    "SB" reference row). Precision is computed here directly from the
    confusion table (correctly-predicted / total-predicted-that-bucket);
    recall comes straight from ``summary_csv``'s own ``recall`` column (same
    numbers ``src.benzene_validation.benzene_internal_confusion_matrix``
    already computes, not re-derived).

    Distinct from ``fig:confusion`` (``plot_confusion_matrix``): that figure
    is the whole-hydride-library (69-molecule) ideal/non-ideal validation;
    this one is benzene's own internal 30 modes only, using benzene's
    genuine literature 3-class ground truth (bend/stretch/SB) rather than
    the library-wide binary stretch/bend reference. Companion per-bond
    numerical detail for the two flagged contrast cases (modes 21/22 --
    literature SB, predicted clean "B", both below tau_B; modes 23/24 --
    literature stretch, predicted "SB"/mixed) lives in
    ``data/results/benzene_sb_vs_stretch_bond_diagnostic.csv`` for a
    companion LaTeX table -- not re-plotted inside this heatmap figure.
    """
    _style()
    tbl = pd.read_csv(matrix_csv, index_col=0)
    summary = pd.read_csv(summary_csv).set_index("ref_label")

    ref_order = ["bend", "stretch", "SB"]
    pred_order = ["bend", "stretch", "mixed"]
    tbl = tbl.reindex(index=ref_order, columns=pred_order, fill_value=0)
    n_total = int(tbl.values.sum())

    fig, (ax_h, ax_p) = plt.subplots(1, 2, figsize=(7.4, 3.4),
                                      gridspec_kw={"width_ratios": [1.15, 1.0]})

    _confusion_heatmap(ax_h, fig, tbl, ref_order,
                        f"(a) Benzene internal modes (n={n_total})")

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

    x = np.arange(len(cats))
    width = 0.35
    bar_colors = [CATEGORY_COLOR[c] for c in cats]
    ax_p.bar(x - width / 2, precisions, width, color=bar_colors,
             edgecolor="black", linewidth=0.5, label="precision")
    ax_p.bar(x + width / 2, recalls, width, color=bar_colors,
             edgecolor="black", linewidth=0.5, hatch="///", label="recall")
    ax_p.set_xticks(x)
    ax_p.set_xticklabels([CATEGORY_LABEL[c] for c in cats], fontsize=9)
    ax_p.set_ylim(0, 1.18)
    ax_p.set_ylabel("Precision / recall")
    ax_p.set_title("(b) Benzene precision/recall", loc="left",
                    fontweight="bold", fontsize=9)
    ax_p.legend(loc="upper left", frameon=False, fontsize=8)
    for xi, p, r in zip(x, precisions, recalls):
        ax_p.text(xi - width / 2, p + 0.02, f"{p:.3f}", ha="center", va="bottom", fontsize=6.5)
        ax_p.text(xi + width / 2, r + 0.02, f"{r:.3f}", ha="center", va="bottom", fontsize=6.5)

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary_dict = {
        "pdf": pdf_path, "png": png_path,
        "shared_categories": ("reuses CATEGORY_COLOR/CATEGORY_LABEL/"
                               "_confusion_heatmap from plot_confusion_matrix "
                               "(fig:confusion) via the new "
                               "REF_LABEL_TO_CATEGORY['SB']='mixed' routing -- "
                               "does not modify fig:confusion itself."),
        "n_total": n_total,
        "confusion_table": tbl.to_dict(),
        "precision": dict(zip(cats, precisions)),
        "recall": dict(zip(cats, recalls)),
        "note": ("Companion per-bond evidence for the 21/22 (SB->bend blind "
                 "spot) vs. 23/24 (stretch->mixed) contrast lives in "
                 "data/results/benzene_sb_vs_stretch_bond_diagnostic.csv "
                 "for a separate LaTeX table -- not plotted here."),
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

    fig, ax = plt.subplots(figsize=(3.6, 3.4))
    for ideal_flag in ("no", "yes"):  # non-ideal first (background), ideal on top
        for ref in ("bend", "stretch"):
            cat = REF_LABEL_TO_CATEGORY[ref]
            sub = bonds[(bonds["ideal"] == ideal_flag) & (bonds["ref_label"] == ref)]
            if sub.empty:
                continue
            kw = _marker_kwargs(cat, ideal_flag, marker="o")
            ax.scatter(sub["abs_rel_db"], sub["s_AB"], s=14,
                       zorder=3 if ideal_flag == "yes" else 2, **kw)

    ax.set_xlabel(r"$|\Delta|\mathbf{b}|\,/\,|\mathbf{b}|\,|$ (relative bond-length change)")
    ax.set_ylabel(r"Bond score $s_{AB}$")
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
               markerfacecolor=COLORS["stretching"], markeredgecolor=COLORS["stretching"],
               markersize=6, label="S, ideal"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="none", markeredgecolor=COLORS["stretching"],
               markersize=6, label="S, non-ideal"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor=COLORS["bending"], markeredgecolor=COLORS["bending"],
               markersize=6, label="B, ideal"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="none", markeredgecolor=COLORS["bending"],
               markersize=6, label="B, non-ideal"),
    ]
    ax.legend(handles=legend_elems, loc="upper left", frameon=False, fontsize=6.8,
              handletextpad=0.4, labelspacing=0.4, borderaxespad=0.2)

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary = {
        "pdf": pdf_path, "png": png_path,
        "shared_categories": ("reuses CATEGORY_COLOR for STRETCHING "
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
    # to a two-level tick scheme: short primary labels only, comfortably
    # narrow, plus a single shared "ideal"/"non-ideal" group annotation
    # (with an under-bracket) spanning each pair, which only has to appear
    # ONCE per pair rather than once per box.
    group_labels = ["B", "S", "B", "S"]
    # Positions: gap 1.3 within a bend/stretch pair, gap 1.6 between the
    # ideal pair (1,2) and non-ideal pair (3,4) -- sized (see
    # IMPLEMENTATION_PLAN.md 2026-07-02 changelog entry) so neither the
    # "bending"/"stretching" tick-label text nor the pair-level "ideal"/
    # "non-ideal" bracket labels below them collide, verified by rendering.
    positions = [1.0, 2.3, 3.9, 5.2]
    pair_spans = [(positions[0], positions[1], "ideal"),
                  (positions[2], positions[3], "non-ideal")]

    panels = [
        ("freq", r"Frequency (cm$^{-1}$)", "(a) Frequency"),
        ("delta_b_mean", r"Averaged $|\Delta|\mathbf{b}|\,/\,|\mathbf{b}|\,|$", "(b) Bond-length change"),
        ("V_Stretch", r"$s[\mathrm{V_S}]$", "(c) Mode score"),
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
            color = CATEGORY_COLOR[cat]
            style = IDEAL_STYLE[ideal_flag]
            box.set_edgecolor(color)
            box.set_linewidth(1.1)
            box.set_facecolor(color if style["filled"] else "white")
            box.set_alpha(0.9 if style["filled"] else 1.0)
        for whisk in bp["whiskers"]:
            whisk.set_color("#555555")
        for cap in bp["caps"]:
            cap.set_color("#555555")
        for flier, (ref, _) in zip(bp["fliers"], groups):
            flier.set_markerfacecolor(CATEGORY_COLOR[REF_LABEL_TO_CATEGORY[ref]])
            flier.set_markeredgecolor(CATEGORY_COLOR[REF_LABEL_TO_CATEGORY[ref]])

        ax.set_xticks(positions)
        # Horizontal, unrotated (reverted 2026-07-02 alongside the gloss
        # revert above): the 30 deg rotation was only needed to avoid
        # collisions with the longer glossed "S (stretching)"/"B (bending)"
        # labels; bare single-character "S"/"B" labels have no collision
        # risk at this position spacing and read cleaner flat, verified by
        # rendering.
        ax.set_xticklabels(group_labels, fontsize=6.5)
        ax.set_xlim(positions[0] - 0.7, positions[-1] + 0.7)
        ax.set_ylabel(ylabel)
        ax.set_title(title, loc="left", fontweight="bold", fontsize=9)

        # Pair-level "ideal"/"non-ideal" bracket + label, in the axes'
        # x-data/y-axes-fraction mixed transform so it sits at a fixed
        # vertical offset below the primary tick labels regardless of each
        # panel's own y-data range. clip_on=False since this offset is
        # deliberately outside the data area (bbox_inches="tight" on save
        # still captures it).
        trans = ax.get_xaxis_transform()
        for lo, hi, text in pair_spans:
            mid = (lo + hi) / 2
            ax.plot([lo - 0.45, hi + 0.45], [-0.28, -0.28], transform=trans,
                    color="#555555", lw=0.7, clip_on=False)
            ax.plot([lo - 0.45, lo - 0.45], [-0.28, -0.24], transform=trans,
                    color="#555555", lw=0.7, clip_on=False)
            ax.plot([hi + 0.45, hi + 0.45], [-0.28, -0.24], transform=trans,
                    color="#555555", lw=0.7, clip_on=False)
            ax.annotate(text, xy=(mid, -0.34), xycoords=trans, ha="center",
                        va="top", fontsize=7.5, style="italic",
                        color="#333333", annotation_clip=False)

    # tau_S / tau_B reference lines on panel (c) only. Extra right-hand xlim
    # padding (vs. the other two panels) so the tau labels have clear room
    # and don't sit flush against the panel's right edge.
    th = Thresholds.calibrated()
    axes[2].set_xlim(positions[0] - 0.7, positions[-1] + 1.05)
    axes[2].axhline(th.tau_S, color=COLORS["threshold"], ls="--", lw=0.8)
    axes[2].axhline(th.tau_B, color=COLORS["threshold"], ls="--", lw=0.8)
    tau_label_x = positions[-1] + 0.55
    axes[2].text(tau_label_x, th.tau_S, r"$\tau_S$", ha="left", va="center", fontsize=7,
                 color=COLORS["threshold"])
    axes[2].text(tau_label_x, th.tau_B, r"$\tau_B$", ha="left", va="center", fontsize=7,
                 color=COLORS["threshold"])

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary = {
        "pdf": pdf_path, "png": png_path,
        "shared_categories": ("reuses CATEGORY_COLOR for STRETCHING/BENDING "
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

    for ax, ideal_flag, title in ((ax_i, "yes", "(a) Ideal molecules"),
                                   (ax_n, "no", "(b) Non-ideal molecules")):
        for ref in ("bend", "stretch"):
            cat = REF_LABEL_TO_CATEGORY[ref]
            sub = internal[(internal.ideal == ideal_flag) & (internal.ref_label == ref)]
            if sub.empty:
                continue
            kw = _marker_kwargs(cat, ideal_flag, marker="o")
            ax.scatter(sub["delta_b_mean"], sub["V_Stretch"], s=20,
                       label=CATEGORY_LABEL[cat], **kw)
        ax.axhline(th.tau_S, color=COLORS["threshold"], ls="--", lw=0.8, zorder=1)
        ax.axhline(th.tau_B, color=COLORS["threshold"], ls="--", lw=0.8, zorder=1)
        ax.set_xlabel(r"Averaged $|\Delta|\mathbf{b}|\,/\,|\mathbf{b}|\,|$")
        ax.set_title(title, loc="left", fontweight="bold", fontsize=9)
        ax.set_xlim(-0.03, internal["delta_b_mean"].max() * 1.08)

    ax_i.set_ylabel(r"$s[\mathrm{V_S}]$")
    ax_i.set_ylim(-0.05, 1.08)
    ax_i.legend(loc="center right", frameon=False, fontsize=7.5,
                handletextpad=0.4, labelspacing=0.4)
    ax_i.text(0.98, th.tau_S, r"$\tau_S$", ha="right", va="bottom", fontsize=7,
              color=COLORS["threshold"], transform=ax_i.get_yaxis_transform())
    ax_i.text(0.98, th.tau_B, r"$\tau_B$", ha="right", va="top", fontsize=7,
              color=COLORS["threshold"], transform=ax_i.get_yaxis_transform())

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    n_ideal = int((internal.ideal == "yes").sum())
    n_nonideal = int((internal.ideal == "no").sum())
    summary = {
        "pdf": pdf_path, "png": png_path,
        "shared_categories": ("reuses CATEGORY_COLOR for STRETCHING/BENDING; "
                               "ONE marker shape (circle) for all points "
                               "(color already distinguishes stretch/bend) "
                               "and the shared IDEAL_STYLE filled=ideal/"
                               "hollow=non-ideal encoding -- same mapping as "
                               "fig:benzene / fig:bondscores / fig:boxplots."),
        "n_ideal_modes": n_ideal,
        "n_nonideal_modes": n_nonideal,
        "tau_S": th.tau_S, "tau_B": th.tau_B,
        "irrep_degeneracy_panel": "NOT BUILT -- pending lead-author/tex-data-sync confirmation.",
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

    fig, ax1 = plt.subplots(figsize=(4.6, 3.4))
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
    ax1.legend(handles=handles, loc="center left", frameon=False, fontsize=7,
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

    Returns a dict {tex_label: result_dict} for all 7 figures (6 labeled
    fig:* + the standalone benzene-normal-modes gallery), in the same order
    they are built. Raises whatever the underlying plot_* function raises
    (e.g. a missing input CSV) -- fail loud, no silent partial regeneration.
    """
    fns = [
        ("fig:benzene", plot_benzene_stress_test),
        ("benzene-normal-modes gallery (no fig: label yet)", plot_benzene_normal_modes),
        ("fig:confusion", plot_confusion_matrix),
        ("fig:benzeneconfusion", plot_benzene_internal_confusion),
        ("fig:bondscores", plot_bond_scores),
        ("fig:boxplots", plot_boxplots),
        ("fig:modemixing", plot_mode_mixing),
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
