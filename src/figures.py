"""Figure generation for the JCC manuscript ("A Unified, Reference-Free
Framework for Classifying the 3N modes of molecular motion").

One function per figure so ``reproduce.py`` can call each in turn. Each
function reads already-computed results CSVs under ``data/results/`` (never
recomputes scores -- presentation-only), builds the figure with the shared
house style (``_style()``), writes a vector PDF + >=300 dpi PNG under
``data/figures/``, and returns a summary dict (paths + sanity numbers).

CANONICAL SCHEME (2026-08 binary-classification switch): ``binary`` is now
THE scheme for the whole paper (``classifier.classify_all_modes``'s own
default). Consequently the CANONICAL (unsuffixed) figure names below are now
produced by what used to be the "binary sibling" functions;
their three-way counterparts are kept fully functional -- not deleted, not
weakened -- but their outputs moved to ``_threeway``-suffixed names, mirroring
the pre-existing ``_binary`` suffix convention (now retired since binary is
no longer the suffixed variant). Every ``plot_*_binary`` function below now
writes the UNSUFFIXED canonical name; every corresponding non-``_binary``
three-way function now writes a ``_threeway``-suffixed name. See
``src/classifier.py``'s ``vib_label_binary``/
``classify_all_modes(..., scheme="binary")`` for the underlying scheme.

Main-text figures: ``plot_benzene_stress_test`` (fig:benzene),
``plot_benzene_emit_counts`` (fig:benzeneemitcounts, overview bar chart of
how benzene's 36 EMIT modes classify across T/R axes and S/B/SB),
``plot_confusion_matrix_binary`` (fig:confusion, CANONICAL as of the binary
switch -- single-tier non-ideal layout, the genuine non-circular validation;
thresholds are fixed on the ideal population and applied without retuning;
three-way sibling ``plot_confusion_matrix`` now writes
``fig_confusion_threeway``), ``plot_bond_scores``
(fig:bondscores), ``plot_boxplots`` (fig:boxplots), ``plot_mode_mixing``
(fig:modemixing), ``plot_sensitivity`` (fig:sensitivity -- the tau_TR sweep;
NOT renamed/swapped with ``plot_tau_sb_sensitivity``, a genuinely different
sweep over a different threshold, see that function's docstring).

SI/standalone figures (no ``fig:`` label): ``plot_benzene_normal_modes``
(all 36 real normal modes, worked examples called out; companion to
fig:benzene, does not replace it), ``plot_benzene_internal_confusion_binary``
(fig:benzeneconfusion, CANONICAL as of the binary switch -- benzene's own
2x2 bend/stretch internal confusion matrix; three-way sibling
``plot_benzene_internal_confusion`` now writes ``fig_benzene_confusion_threeway``
and keeps the genuine 3x3 bend/stretch/SB breakdown enabled by the
literature "SB" relabeling of modes 21/22),
``plot_benzene_confusion_precision_recall`` (the precision/recall numbers
for the three-way matrix, split into its own figure since a formal-metric
panel sits awkwardly next to this section's deliberate "mixed by eye is a
convention" hedging), ``plot_rigorous_tier_check`` (the removed ideal-tier
confusion panels from fig:confusion -- precision/recall=1.000 there is
circular by construction, kept only as a self-consistency check, not an
accuracy claim), ``plot_transferability_confusion_binary``
(fig:transferabilityconfusion, PROPOSED, CANONICAL as of the binary switch --
the 18-molecule held-out "test" tier's own joint T/R/B/S confusion matrix
[literal literature "SB" reference rows DROPPED, not bucketed], genuinely
non-circular since this tier is excluded from tau_S/tau_B calibration and
its internal ground truth is independently hand-curated in
``characterised_modes.csv``; T and R are kept SEPARATE, unlike
fig:benzeneconfusion's collapsed "T/R"; three-way sibling
``plot_transferability_confusion`` now writes
``fig_transferability_confusion_threeway`` and keeps the 5-way T/R/B/SB/S
breakdown), ``plot_irrep_coupling``
(irrep-degeneracy mixing-mechanism figure: mode
score vs. central-atom displacement amplitude, faceted by same-irrep
coupling partner; reads ``irrep``/``shape``/``type`` from
``characterised_modes.csv``, ideal/non-ideal from ``mol_list_method.csv``'s
``mol_type``, and scores from ``library_scores.csv``), ``plot_cpu_time_benchmark``
(fig:cputime, PROPOSED label not yet wired into the .tex -- empirical CPU-time
figure for the "Computational cost" section as it is rewritten away from a
pure Big-O argument; reads ``data/results/cpu_time_benchmark.csv``),
``plot_tau_sb_sensitivity`` (error_all/error_test/error_single_centre vs.
tau_SB line plot, mirrors ``plot_sensitivity``'s structure but sweeps the
single tau_SB cutoff, not tau_TR -- reads
``data/results/tau_sb_sensitivity_sweep.csv`` from
``src.calibrate.run_tau_sb_error_analysis``).

Both schemes' figures are always registered in ``regenerate_all()`` (see its
``fns`` list) -- ``--figures`` regenerates both the canonical (binary) and
``_threeway``-suffixed outputs every run, so neither scheme is ever silently
stale or disabled.

Cross-figure visual consistency: every figure encoding a classification
category reuses the same ``CATEGORY_COLOR``/``CATEGORY_MARKER``/
``CATEGORY_LABEL`` mapping defined once below, never redefined per-function.
Ideal vs. non-ideal molecule membership is encoded via ``IDEAL_STYLE`` as
filled (ideal) vs. hollow (non-ideal) markers/boxes, also centralized here.
"""

from __future__ import annotations

import dataclasses
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.patheffects as pe
import matplotlib.ticker as mtick
from scipy.optimize import curve_fit

from src.classifier import (
    Thresholds, classification_bucket,
    is_clean_external, is_mixed_external, external_axis,
    vib_label_binary,
)
from src.csv_label_ingest import TYPE_TO_REF_LABEL

# --------------------------------------------------------------------------
# Shared house style
# --------------------------------------------------------------------------

# INVARIANT: "bending"/"stretching" here are for CLASSIFICATION-ALGORITHM
# (predicted) labels only -- CATEGORY_COLOR, wherever a color encodes what
# the classifier itself predicted. Reference/ground-truth (literature)
# labels use the separate Okabe-Ito "bending_ref"/"stretching_ref" pair
# below (REF_CATEGORY_COLOR) -- do not conflate the two.
COLORS = {
    "external": "#999999",   # gray      -- CLEAN_TRANSLATION / CLEAN_ROTATION
    "bending": "#1A5DE4",    # blue      -- BENDING (classification-algorithm labels only)
    "stretching": "#E41A1C", # red       -- STRETCHING
    "mixed": "#984EA3",      # purple    -- MIXED_STRETCH_BEND
    "bending_ref": "#0072B2",    # Okabe-Ito blue -- reference/literature labels only
    "stretching_ref": "#D55E00", # Okabe-Ito vermillion -- reference/literature labels only
    "mixed_ext": "#CC79A7",  # pink      -- MIXED_EXTERNAL_WITH_VIBRATION
    "background": "#BBBBBB", # light gray-- unhighlighted context points
    "highlight_r": "#F0E442",# yellow    -- s[R] inversion (EMIT 2 / 9)
    "highlight_t": "#000000",# black     -- flagged external (EMIT 34-36)
    "threshold": "#555555",  # dark gray -- tau reference lines (all figures)
    "sens_accuracy": "#56B4E9",   # sky blue  -- accuracy curve, fig:sensitivity
    "sens_change": "#000000",     # black     -- label-change-fraction curve, fig:sensitivity
    "sens_single_centre": "#009E73",  # Okabe-Ito bluish green -- single-centre error curve, fig:sensitivity_binary
    "plateau_band": "#56B4E9",    # sky blue @ low alpha -- tau_TR plateau shading, fig:sensitivity
    # Sequential, colorblind-safe -- fig:confusion/fig:benzeneconfusion
    # heatmap. NOTE: YlGnBu's high-value end is blue-ish, a known collision
    # risk with the "bending" tick-label color -- knowingly accepted.
    "confusion_cmap": "YlGnBu",
    # fig:cputime only -- not reused by any classification-category encoding
    # elsewhere in this module (this figure doesn't encode a T/R/S/B
    # category, so the "one color = one category everywhere" rule above
    # doesn't apply to it). Gaussian's own points are colored by a
    # continuous inferno gradient keyed to n_basis (built inline in
    # ``plot_cpu_time_benchmark``, not a fixed COLORS entry -- there is no
    # single "Gaussian color" anymore); its per-N mean trend line is plain
    # black. "cost_classifier" (Okabe-Ito bluish green) is untouched by the
    # gradient, since classifier cost doesn't depend on AO basis-function
    # count.
    "cost_classifier": "#009E73",   # bluish green -- classifier CPU time
}

# Reference-label / predicted-bucket -> shared classification-CATEGORY name.
# These are IDENTITY maps (ref_label/pred_bucket values already ARE the
# CATEGORY bucket names) -- kept, not deleted, so ~15 call sites elsewhere
# in this file (`CATEGORY_COLOR[REF_LABEL_TO_CATEGORY[r]]` etc.) don't need
# touching, and so a future "simplification" doesn't remove them and break
# those call sites. For a RAW per-mode classifier label, convert with
# `classification_bucket()` first, then use the result as the key.
REF_LABEL_TO_CATEGORY = {
    "translation": "translation",
    "rotation": "rotation",
    "stretch": "stretch",
    "bend": "bend",
    # "SB": benzene modes 21/22's literal literature mixed-stretch/bend
    # label, routed to the same "mixed" category as the classifier's own
    # MIXED_STRETCH_BEND bucket so CATEGORY_COLOR/CATEGORY_LABEL and
    # `_confusion_heatmap` can be reused as-is for fig:benzeneconfusion.
    "SB": "mixed",
    # Individual-axis identity entries (joint T/R+internal confusion
    # matrices): a bare clean-external label ("Tx".."Rz") is its own
    # category here, not collapsed to "translation"/"rotation".
    "Tx": "Tx", "Ty": "Ty", "Tz": "Tz",
    "Rx": "Rx", "Ry": "Ry", "Rz": "Rz",
    # Collapsed entry (fig:benzeneconfusion only): all 6 axes folded into
    # one "T/R" row/column, since per-axis detail is already shown elsewhere.
    "T/R": "T/R",
    # Bare T/R (kept SEPARATE, unlike "T/R" above) + bare B/S entries:
    # fig:transferabilityconfusion's 5-way T/R/B/SB/S table
    # (plot_transferability_confusion). "SB" needs no entry of its own --
    # it already routes to "mixed" above, which round-trips through
    # CATEGORY_LABEL back to "SB". "B"/"S" route to the existing
    # "bend"/"stretch" categories (correct short labels + colors already).
    # "T"/"R" are self-mapped to their OWN category (like "Tx".."Rz" are),
    # NOT to "translation"/"rotation" -- that category's CATEGORY_LABEL text
    # is "clean translation"/"clean rotation", wrong for this bare-letter
    # table (needs its own CATEGORY_COLOR/REF_CATEGORY_COLOR entries below).
    "T": "T", "R": "R", "B": "bend", "S": "stretch",
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
    # Bare T/R self-mapped category (fig:transferabilityconfusion) -- same
    # external gray as every other T/R-flavored entry above.
    "T": COLORS["external"], "R": COLORS["external"],
}

# Reference/ground-truth (literature) counterpart of CATEGORY_COLOR:
# identical except "bend"/"stretch" use the Okabe-Ito pair
# ("bending_ref"/"stretching_ref") instead of the algorithm's own hues.
# Used wherever a color encodes a REFERENCE/literature label, never a
# classifier-predicted one (those keep CATEGORY_COLOR).
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

# Internal stretch/bend/mixed use the SAME short symbols the classifier
# itself emits ("S"/"B"/"SB", matching tab:benzenemixed's `\texttt{}`
# labels) -- bare symbols, not glossed in-figure; the gloss belongs once in
# the LaTeX caption, not repeated on every legend/tick across 5 figures.
# translation/rotation have no single-letter equivalent, so they stay long-form.
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
    # Bare T/R/B/S tick text for fig:transferabilityconfusion (see the
    # matching REF_LABEL_TO_CATEGORY entries above).
    "T": "T", "R": "R", "B": "B", "S": "S",
}

# Named font-size overrides: the ONLY permitted exceptions to _style()'s
# rcParams defaults, tuned for dense panels (crowded legends / thin bar
# value labels) where the default size collides with other text. Reuse
# these rather than adding a new bare `fontsize=<number>` elsewhere.
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
    explicitly so "bend"/"stretch" render in the Okabe-Ito blue/vermillion
    pair rather than the classifier's own hues.

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
    emit_csv="data/results/C6H6_EMIT.csv",
    out_dir="data/figures",
    label="fig_benzene",
):
    """Build fig:benzene: score vs. projected normal-mode contribution for
    the 36 EMIT modes, highlighting the EMIT 2/9 s[R] inversion and the
    EMIT 34-36 flagged externals.

    This is a flag-behavior / non-monotonicity illustration, NOT a
    T/R-accuracy parity plot: the figure deliberately omits any 1:1
    reference line and plots |score| against projected contribution
    fraction only, to expose where the two disagree.

    Single-panel: the normal-mode s[V_S]-vs-frequency content lives
    separately in ``plot_benzene_normal_modes``/``fig_benzene_normal``.

    Requires `emit_csv` to already have the C2_* projection columns merged
    in (`python main.py -m benzene --mode emit` then `--emit-projection`).

    Returns a summary dict with output paths and a few sanity numbers.
    """
    _style()

    emit = pd.read_csv(emit_csv)

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
    # highlight_r's yellow has poor contrast as a bare line on white -- a
    # thin black halo keeps the arrow legible without changing its color.
    ann_r.arrow_patch.set_path_effects(
        [pe.Stroke(linewidth=2.0, foreground="black"), pe.Normal()])
    # EMIT 9 label offset upper-left, clear of the arrow's approach angle.
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
    # Caption-text sentence, not drawn in-image (would shrink below a
    # readable floor once LaTeX rescales the figure to column width).
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

def plot_benzene_normal_modes(
    normal_csv="data/results/C6H6_normal.csv",
    out_dir="data/figures",
    label="fig_benzene_normal",
    tau_SB=None,
):
    """Build the benzene normal-mode worked-example gallery: ``s[V_S]`` vs.
    frequency for all 36 real normal modes (6 external T/R + 30 internal),
    colored by the binary T/R/S/B scheme (shared CATEGORY_COLOR/
    CATEGORY_LABEL; C6H6_normal.csv's own label column, never "SB" under the
    default binary classify-scheme). Separate figure from
    ``plot_benzene_stress_test``
    (fig:benzene, untouched), the descriptive companion for the "Benzene
    normal modes" section.

    Styling: one marker shape (circle) for every point, color-only encoding
    (CATEGORY_MARKER's per-category shapes are not used here); all markers
    rendered hollow (IDEAL_STYLE["no"], a blanket style choice -- there is
    no ideal/non-ideal axis within one molecule's own normal modes); plain,
    unannotated scatter (no callout boxes for modes 12/19/30).

    ``tau_SB`` defaults to ``Thresholds.calibrated().tau_SB`` (None ->
    canonical, zero behavior change). When given, the CSV's own baked-in
    ``label`` column is NOT trusted directly for internal rows -- each row's
    category is cheaply re-derived from its own ``V_Stretch`` via
    ``classifier.rescheme_internal_label`` (external T/R rows pass through
    unchanged, since scheme/tau_SB never touches Step 2/3), and the drawn
    threshold line/label use the same override. No re-scoring involved.
    """
    _style()
    from src.classifier import rescheme_internal_label

    normal = pd.read_csv(normal_csv)
    thresholds = Thresholds.calibrated()
    if tau_SB is not None:
        thresholds = dataclasses.replace(thresholds, tau_SB=tau_SB)
    TAU_SB = thresholds.tau_SB

    # Width 7.2in leaves horizontal room to manually composite depicted
    # normal-mode panel images beside the scatter.
    fig, ax = plt.subplots(figsize=(7.2, 3.6))

    # Legend dedup: translation/rotation share the identical gray color and
    # marker here, so a naive per-`cat` dedup would put two visually-
    # identical entries ("clean translation", "clean rotation") in the
    # legend -- merge them into one "clean T/R" entry (local to this figure
    # only; CATEGORY_LABEL itself stays untouched).
    _LEGEND_MERGE_KEY = {"translation": "clean_tr", "rotation": "clean_tr"}
    _LEGEND_MERGE_TEXT = {"clean_tr": "T/R"}

    seen_labels = set()
    for _, row in normal.iterrows():
        row_label = row["label"]
        if tau_SB is not None:
            row_label = rescheme_internal_label(row_label, row["V_Stretch"],
                                                 thresholds, scheme="binary")
        cat = classification_bucket(row_label)
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

    ax.axhline(TAU_SB, color=COLORS["threshold"], ls="--", lw=0.8, zorder=1)
    y_offset = 0.045
    # Right side (ha="right"), clear of the legend which sits on the left.
    ax.text(0.98, TAU_SB + y_offset, r"$\tau_{\mathrm{SB}}=$" + f"{TAU_SB:.2f}",
            ha="right", va="bottom", color=COLORS["threshold"],
            transform=ax.get_yaxis_transform())

    ax.set_xlabel(r"Frequency (cm$^{-1}$)")
    ax.set_ylabel(r"$s[\mathrm{V_S}]$")
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(-120, normal["Freq"].max() * 1.06)
    # Upper-left corner: clear of both the tau_SB line at y=0.50 (a plain
    # "center left" placement puts a legend row right on top of it -- the S
    # row sits at that same height) and of all data (nothing above y~0.7
    # below freq~800).
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.02), frameon=False,
              handletextpad=0.3, labelspacing=0.3, borderaxespad=0.2)

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary = {
        "pdf": pdf_path, "png": png_path,
        "shared_categories": ("reuses CATEGORY_COLOR/CATEGORY_LABEL for the "
                               "binary T/R/S/B scheme -- same mapping as "
                               "fig:benzene panel (a). ONE marker shape "
                               "(circle) for all points, all rendered "
                               "hollow (IDEAL_STYLE['no']); color is the "
                               "only category encoding, matching "
                               "fig:bondscores/fig:boxplots/fig:modemixing's "
                               "convention."),
        "n_points": len(normal),
        "freq_range": (float(normal["Freq"].min()), float(normal["Freq"].max())),
        "vs_range": (float(normal["V_Stretch"].min()), float(normal["V_Stretch"].max())),
        "tau_SB": TAU_SB,
    }
    return summary


# --------------------------------------------------------------------------
# fig:confusion -- single-tier (non-ideal only) clean-category confusion
# matrix + retention-and-migration bars, plus a separate SI rigorous-tier
# consistency check (plot_rigorous_tier_check, below). See
# plot_confusion_matrix's own docstring for the circularity argument that
# motivates restricting the main-text figure to the non-ideal tier only.
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
    are present for that tier). Returns the imshow handle (caller attaches
    its own colorbar). No in-figure title by design -- mode-count context
    lives in the LaTeX caption only. `label_map` lets a caller override
    tick text (e.g. "T"/"R" instead of "clean translation") without
    touching the shared `CATEGORY_LABEL`. `divider=(n_ref_ext, n_pred_ext)`,
    if given, draws a thicker line separating a leading external-axis block
    from the trailing internal-category block.
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

    ax.set_xticks(range(len(tbl.columns)))
    ax.set_xticklabels([label_map[PRED_BUCKET_TO_CATEGORY[c]] if c in
                         PRED_BUCKET_TO_CATEGORY else c for c in tbl.columns],
                        rotation=0, ha="center")
    ax.set_yticks(range(len(tbl.index)))
    ax.set_yticklabels([label_map[REF_LABEL_TO_CATEGORY[r]] for r in tbl.index])
    # Row ticks ("Reference label") use REF_CATEGORY_COLOR; column ticks
    # ("Classification algorithm label") use CATEGORY_COLOR -- see the
    # ref-vs-predicted color-split note on COLORS above.
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
    # imshow's own cell grid already delineates the matrix, so the outer
    # axes frame (all 4 spines) and colorbar outline are redundant.
    for spine in ax.spines.values():
        spine.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Number of modes")
    cb.outline.set_visible(False)
    return im


def plot_confusion_matrix(
    library_csv="data/results/library_scores.csv",
    out_dir="data/figures",
    label="fig_confusion_threeway",
):
    """Three-way sibling of ``plot_confusion_matrix_binary`` (which now owns
    the canonical ``fig_confusion`` name/output, the paper's default scheme
    as of the 2026-08 binary-classification switch) -- kept fully functional
    for comparison/on-demand use, not deleted. Builds a single-panel JOINT
    confusion matrix: the non-ideal tier's internal rows (``ideal == 'no'``)
    plus the external (Tx..Rz) reference rows for those same molecules, in
    one crosstab. Joint (not internal-only) because it verifies that no
    internal mode's residual character ever wins an external Step-2 slot or
    vice versa -- the off-diagonal external<->internal blocks are computed
    and reported
    (``crossover_ext_ref_to_internal_pred``/``crossover_internal_ref_to_ext_pred``),
    not assumed zero.

    CIRCULARITY NOTE: thresholds are fixed on the ideal population
    (fig:boxplots), then applied WITHOUT retuning to these harder non-ideal
    cases -- this is the genuine, non-circular validation. The
    fully-rigorous tier (ideal population) is exact by construction (tau_S/
    tau_B are literally derived from it, so "1.000 accuracy" there is
    circular) and lives separately in ``plot_rigorous_tier_check`` (SI, a
    self-consistency check, not an accuracy claim).

    Retention/migration bars live in the separate SI figure
    ``plot_confusion_retention_migration``. ``filter_single_centre_library``
    is applied explicitly (not left to confusion_matrix_stats's internal
    filter) so internal and external rows share the same molecule scope.

    THREE-WAY RESCHEME: ``library_scores.csv``'s own ``predicted_label`` is
    built under the global default scheme (binary as of the 2026-08 switch),
    which never produces "mixed" -- so this figure re-derives every internal
    row's ``predicted_label`` under ``scheme="threeway"`` from its own
    ``V_Stretch`` (``classifier.rescheme_internal_label``, cheap, no
    re-parse/re-run of Step 1-3) before building the confusion table, so the
    "mixed" column stays genuinely populated regardless of what scheme the
    cached CSV happens to be built under.
    """
    _style()
    from src.calibrate import confusion_matrix_stats, filter_single_centre_library
    from src.classifier import rescheme_internal_label

    lib_df = pd.read_csv(library_csv)
    lib_df = filter_single_centre_library(lib_df)
    thresholds = Thresholds.calibrated()

    nonideal_df = lib_df[(lib_df["kind"] == "internal") & (lib_df["ideal"] == "no")].copy()
    nonideal_df["predicted_label"] = nonideal_df.apply(
        lambda row: rescheme_internal_label(row["predicted_label"], row["V_Stretch"],
                                             thresholds, "threeway"),
        axis=1)
    stats_n = confusion_matrix_stats(nonideal_df, thresholds)

    nonideal_molecules = nonideal_df["molecule"].unique()
    nonideal_external_df = lib_df[(lib_df["kind"] == "external") &
                                   (lib_df["molecule"].isin(nonideal_molecules))]
    joint_df = pd.concat([nonideal_df, nonideal_external_df], ignore_index=True)

    # "mixed_external" ("T/R*") dropped: verified empty on this population.
    ref_order_n = ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "bend", "stretch"]
    pred_order_n = ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "bend", "mixed", "stretch"]
    tbl_n = _joint_confusion_table(joint_df, ref_order_n, pred_order_n)
    n_nonideal = int(tbl_n.values.sum())

    # Height tuned to 4.8in so the figure+caption share one page (5.6in
    # pushed the caption to a mostly-blank next page) -- don't reduce below.
    fig, ax_hb = plt.subplots(figsize=(6.2, 4.8))
    _confusion_heatmap(ax_hb, fig, tbl_n, ref_order_n)

    # Internal-only opposite-category crossing (bend<->stretch).
    cats_n = ["bend", "stretch"]
    opposite_n = []
    for c in cats_n:
        opp = "stretch" if c == "bend" else "bend"
        n_ref = stats_n["per_category"][c]["n_ref"]
        n_opp = int(tbl_n.loc[c, opp]) if opp in tbl_n.columns else 0
        opposite_n.append(n_opp / n_ref if n_ref else float("nan"))

    # External<->internal crossover check -- the point of the joint table.
    # Expected 0 both ways, but computed, not assumed.
    ext_rows, int_rows = ref_order_n[:6], ref_order_n[6:]
    ext_cols, int_cols = pred_order_n[:6], pred_order_n[6:]
    crossover_ext_ref_to_internal_pred = int(tbl_n.loc[ext_rows, int_cols].values.sum())
    crossover_internal_ref_to_ext_pred = int(tbl_n.loc[int_rows, ext_cols].values.sum())

    # Caption-text sentence (not drawn in-image). ASCII "->", not a unicode
    # arrow (console-print safety).
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
                   "(stretch/bend/mixed) confusion matrix. Retention/"
                   "migration bars live separately in "
                   "plot_confusion_retention_migration (SI, not main text). "
                   "No divider line between the external and internal "
                   "blocks, and no mixed_external ('T/R*') column (empty "
                   "for this dataset)."),
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


def plot_confusion_matrix_binary(
    library_csv="data/results/library_scores.csv",
    out_dir="data/figures",
    label="fig_confusion",
    tau_SB=None,
):
    """Binary-scheme sibling of plot_confusion_matrix (fig:confusion,
    threeway version now at ``fig_confusion_threeway``) -- CANONICAL as of
    the 2026-08 binary-classification switch (this function now owns the
    ``fig_confusion`` name/output, the paper's default scheme). Same
    non-ideal single-centre scope (filter_single_centre_library + ideal ==
    'no'), same joint external+internal crosstab machinery
    (_joint_confusion_table/_confusion_heatmap), but every internal row that
    the three-way scheme would classify S/B/SB is RE-LABELED via
    vib_label_binary(V_Stretch, tau_SB) -- so the internal side of the
    matrix is 2-way (bend/stretch only, never "mixed"/"SB").

    An internal REFERENCE row whose three-way predicted_label is actually an
    EXTERNAL slot (a Step-2 crossover -- this mode's residual character won
    a T/R slot) is left untouched, not re-labeled: scheme only changes the
    Step-4 internal vocabulary, never Step 2/3's external assignment (see
    src/classifier.py's classify_all_modes docstring), so overriding it here
    would corrupt the same crossover accounting plot_confusion_matrix reports.

    tau_SB defaults to Thresholds.calibrated().tau_SB. Reference rows with
    ref_label=="SB" are dropped (a literal literature "SB" answers "is this
    genuinely mixed?", a different question than a binary S/B classifier can
    answer) -- none are expected in this tier, but checked rather than
    assumed (n_dropped_sb_ref in the summary).
    """
    _style()
    from src.calibrate import filter_single_centre_library

    lib_df = pd.read_csv(library_csv)
    lib_df = filter_single_centre_library(lib_df)
    tau_SB = tau_SB if tau_SB is not None else Thresholds.calibrated().tau_SB

    nonideal_df = lib_df[(lib_df["kind"] == "internal") & (lib_df["ideal"] == "no")].copy()
    n_dropped_sb_ref = int((nonideal_df["ref_label"] == "SB").sum())
    nonideal_df = nonideal_df[nonideal_df["ref_label"] != "SB"]

    def _binarize(row):
        bucket = classification_bucket(row["predicted_label"])
        if bucket in ("stretch", "bend", "mixed"):
            return vib_label_binary(row["V_Stretch"], tau_SB)
        return row["predicted_label"]  # external crossover -- left untouched

    nonideal_df["predicted_label"] = nonideal_df.apply(_binarize, axis=1)

    nonideal_molecules = nonideal_df["molecule"].unique()
    nonideal_external_df = lib_df[(lib_df["kind"] == "external") &
                                   (lib_df["molecule"].isin(nonideal_molecules))]
    joint_df = pd.concat([nonideal_df, nonideal_external_df], ignore_index=True)

    # 2-way internal: no "mixed" column, unlike plot_confusion_matrix's 3-way.
    ref_order_n = ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "bend", "stretch"]
    pred_order_n = ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "bend", "stretch"]
    tbl_n = _joint_confusion_table(joint_df, ref_order_n, pred_order_n)
    n_nonideal = int(tbl_n.values.sum())

    fig, ax_hb = plt.subplots(figsize=(6.2, 4.8))
    _confusion_heatmap(ax_hb, fig, tbl_n, ref_order_n)

    # Internal-only opposite-category crossing (bend<->stretch) -- n_ref
    # read directly off the table's own row sums (no "mixed" bucket to
    # subtract out, unlike confusion_matrix_stats's threeway-specific dict).
    cats_n = ["bend", "stretch"]
    opposite_n = []
    for c in cats_n:
        opp = "stretch" if c == "bend" else "bend"
        n_ref = int(tbl_n.loc[c].sum()) if c in tbl_n.index else 0
        n_opp = int(tbl_n.loc[c, opp]) if opp in tbl_n.columns else 0
        opposite_n.append(n_opp / n_ref if n_ref else float("nan"))

    ext_rows, int_rows = ref_order_n[:6], ref_order_n[6:]
    ext_cols, int_cols = pred_order_n[:6], pred_order_n[6:]
    crossover_ext_ref_to_internal_pred = int(tbl_n.loc[ext_rows, int_cols].values.sum())
    crossover_internal_ref_to_ext_pred = int(tbl_n.loc[int_rows, ext_cols].values.sum())

    footer_text = (
        f"Binary-scheme (tau_SB={tau_SB:g}) non-ideal internal + external "
        f"classification (n={n_nonideal}): 0% of bend or stretch reference-"
        f"labeled modes crossed to the OPPOSITE category "
        f"(bend->stretch={opposite_n[0]:.1%}, stretch->bend={opposite_n[1]:.1%}); "
        f"and {crossover_ext_ref_to_internal_pred} external reference modes were "
        f"predicted into an internal bucket, {crossover_internal_ref_to_ext_pred} "
        "internal reference modes were predicted into an external slot. No 'SB' "
        "column: the binary scheme forces every internal mode to bend or stretch."
    )
    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary = {
        "pdf": pdf_path, "png": png_path,
        "tau_SB": tau_SB,
        "nonideal_footer_text": footer_text,
        "layout": ("Single panel: joint external (Tx..Rz) + internal "
                   "(stretch/bend only, binary scheme) confusion matrix. "
                   "Binary sibling of fig:confusion -- see that figure for "
                   "the three-way (S/B/SB) version."),
        "nonideal_n": n_nonideal,
        "n_dropped_sb_ref": n_dropped_sb_ref,
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
    """Build the SI companion to fig:confusion: the non-ideal-tier
    label-retention vs. migration-to-mixed bar panel. Self-contained --
    recomputes ``confusion_matrix_stats`` itself rather than depending on
    ``plot_confusion_matrix`` having already run.
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
        "framing": ("SI companion to plot_confusion_matrix (fig:confusion); "
                     "kept out of the main-text figure since Figure 6's own "
                     "mode-count breakdown and the manuscript prose already "
                     "state the retention numbers."),
    }
    return summary


def plot_rigorous_tier_check(
    library_csv="data/results/library_scores.csv",
    out_dir="data/figures",
    label="fig_rigorous_tier_check",
    csv_path="data/results/rigorous_tier_consistency_table.csv",
):
    """Build the SI rigorous-tier consistency check: a small reference/
    predicted-count TABLE (not a heatmap, deliberately minimal) covering
    every external (T/R) row plus internal rows with ``ideal == 'yes'``.

    CIRCULARITY NOTE: this is a SANITY CHECK, not an independent accuracy
    claim -- T/R references are built via the same Eckart-Sayvetz
    projection the mode is scored against, and tau_S/tau_B are literally
    the min/max of this exact population. Precision/recall are 1.000 for
    every category by construction; this table shows that computationally
    rather than merely asserting it. fig:confusion's non-ideal tier is the
    genuine non-circular validation.

    Renders a compact table PDF/PNG and writes the same numbers to
    ``data/results/rigorous_tier_consistency_table.csv`` for an optional
    native LaTeX table instead.
    """
    _style()
    from src.calibrate import confusion_matrix_stats, filter_single_centre_library

    lib_df = pd.read_csv(library_csv)
    lib_df = filter_single_centre_library(lib_df)
    thresholds = Thresholds.calibrated()

    rigorous_df = lib_df[(lib_df["kind"] == "external") | (lib_df["ideal"] == "yes")]
    # For acceptance_floor/floor_met (coarse 4-bucket construction check).
    stats_r = confusion_matrix_stats(rigorous_df, thresholds)

    cats_r = ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "bend", "stretch"]
    tbl_r = _joint_confusion_table(rigorous_df, cats_r, cats_r)
    n_rigorous = int(tbl_r.values.sum())
    rows = _per_category_from_table(tbl_r, cats_r)
    table_df = pd.DataFrame(rows)

    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    table_df.to_csv(csv_path, index=False)

    # Minimal table rendering, not a heatmap+bars figure. colWidths are
    # hand-tuned to fit the header text without overlap -- don't let
    # ax.table auto-size them (verified to badly collide on first attempt).
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
# fig:transferabilityconfusion -- the 18-molecule held-out "test" tier's own
# joint T/R/B/SB/S confusion matrix (Phase 6's out-of-sample-evaluation
# instrument). Genuinely non-circular: this tier is excluded from tau_S/tau_B
# calibration by construction (SINGLE_CENTRE_ONLY_EXCLUDE), and its internal
# ground truth is hand-curated in data/characterised_modes.csv, independent
# of the classifier. Unlike fig:confusion/fig:benzeneconfusion, T and R are
# kept as SEPARATE categories (not collapsed to one "T/R" bucket) per the
# 5-way T/R/B/SB/S spec asked for here -- reuses _confusion_heatmap via the
# bare "T"/"R"/"B"/"S" REF_LABEL_TO_CATEGORY/CATEGORY_LABEL entries added
# above, but builds its own crosstab rather than going through
# _joint_confusion_table (whose collapse_external=True merges T+R together,
# and whose _axis_aware_pred_category discards axis identity for any
# mixed-external prediction -- both wrong for this table).
# --------------------------------------------------------------------------

def plot_transferability_confusion(
    library_csv="data/results/library_scores.csv",
    out_dir="data/figures",
    label="fig_transferability_confusion_threeway",
    matrix_csv_path="data/results/transferability_confusion_threeway_matrix.csv",
    summary_csv_path="data/results/transferability_confusion_threeway_summary.csv",
    misclassified_csv_path="data/results/transferability_confusion_threeway_misclassified.csv",
):
    """Three-way sibling of ``plot_transferability_confusion_binary`` (which
    now owns the canonical ``fig_transferability_confusion`` name/output, the
    paper's default scheme as of the 2026-08 binary-classification switch)
    -- kept fully functional for comparison/on-demand use, not deleted.
    Builds a single 5x5 joint confusion matrix (T, R, B, SB, S) pooling all
    18 ``mol_type=='test'`` molecules' modes -- external (Tx..Rz) rows exact
    by Eckart-Sayvetz construction, internal (B/S/SB) rows from the
    hand-curated ``ref_label`` column ``src/library_ingest.py::attach_labels()``
    joins in from ``data/characterised_modes.csv``.

    THREE-WAY RESCHEME: like ``plot_confusion_matrix``, this figure
    re-derives every internal row's ``predicted_label`` under
    ``scheme="threeway"`` from its own ``V_Stretch``
    (``classifier.rescheme_internal_label``) before use, since
    ``library_scores.csv``'s own column is built under the global default
    scheme (binary) and would otherwise never populate the "SB" prediction.

    T/R ground truth: first character of the external row's ``mode_index``
    ("Tx"->"T", "Rz"->"R"). T/R prediction: ``classifier.external_axis()``
    (strips a mixed-external "*" suffix before taking the first character,
    so an impure "Tx*" still counts as predicted "T" -- confirmed empirically
    absent for this tier's current data, 0/498 rows starred, but handled for
    correctness). B/SB/S on both sides go through the ordinary
    ``ref_label``/``classification_bucket()`` vocabulary. Any internal row
    still missing ``ref_label`` (stale ``library_scores.csv``, i.e. ``python
    main.py --library`` hasn't been rerun since the last
    ``characterised_modes.csv`` edit) is dropped and counted in
    ``n_dropped_missing_ground_truth`` rather than silently ignored.

    Also writes ``misclassified_csv_path``: one row per off-diagonal mode
    (``ref_category != predicted_category``), molecule/mode/freq/irrep plus
    the literature ``description``/``ref`` columns joined in directly from
    ``data/characterised_modes.csv`` (not carried by ``library_scores.csv``
    itself) -- lets a reader open one file and see exactly which modes to
    look at, not just the aggregate counts in ``matrix_csv_path``.
    """
    _style()
    from src.classifier import classification_bucket, external_axis, rescheme_internal_label
    from src.csv_label_ingest import load_label_csvs
    from src.library_ingest import test_tier_molecules

    lib_df = pd.read_csv(library_csv)
    test_mols = test_tier_molecules()
    df = lib_df[lib_df["molecule"].isin(test_mols)].copy()

    thresholds = Thresholds.calibrated()
    internal_mask = df["kind"] == "internal"
    df.loc[internal_mask, "predicted_label"] = df.loc[internal_mask].apply(
        lambda row: rescheme_internal_label(row["predicted_label"], row["V_Stretch"],
                                             thresholds, "threeway"),
        axis=1)

    ext_mask = df["kind"] == "external"
    ref_internal_map = {"bend": "B", "stretch": "S", "SB": "SB"}
    ref_cat = pd.Series(index=df.index, dtype=object)
    ref_cat[ext_mask] = df.loc[ext_mask, "mode_index"].astype(str).str[0]
    ref_cat[~ext_mask] = df.loc[~ext_mask, "ref_label"].map(ref_internal_map)

    pred_internal_map = {"stretch": "S", "bend": "B", "mixed": "SB"}

    def _pred_category(predicted_label):
        axis = external_axis(predicted_label)
        if axis is not None:
            return axis[0]
        return pred_internal_map.get(classification_bucket(predicted_label))

    pred_cat = df["predicted_label"].map(_pred_category)

    keep = ref_cat.notna()
    n_dropped = int((~keep).sum())
    n_folded_mixed_external = int(
        df["predicted_label"].astype(str).str.endswith("*").sum())

    ref_order = ["T", "R", "B", "SB", "S"]
    tbl = pd.crosstab(
        pd.Series(ref_cat[keep], name="ref"),
        pd.Series(pred_cat[keep], name="pred"),
    )
    tbl = tbl.reindex(index=ref_order, columns=ref_order, fill_value=0)
    n_total = int(tbl.values.sum())

    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    _confusion_heatmap(ax, fig, tbl, ref_order)
    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    os.makedirs(os.path.dirname(matrix_csv_path), exist_ok=True)
    tbl.to_csv(matrix_csv_path)
    rows = _per_category_from_table(tbl, ref_order)
    pd.DataFrame(rows).to_csv(summary_csv_path, index=False)

    # Per-mode detail for every off-diagonal (misclassified) mode, joined
    # against characterised_modes.csv's own literature description/ref
    # columns -- library_scores.csv only carries ref_key, not description.
    cm = load_label_csvs()["characterised_modes"][
        ["molecule", "mode", "description", "ref"]].copy()
    cm["mode"] = pd.to_numeric(cm["mode"], errors="coerce")
    mis_mask = keep & (ref_cat != pred_cat)
    mis = df.loc[mis_mask, ["molecule", "mode_index", "kind", "freq",
                             "irrep", "V_Stretch", "predicted_label"]].copy()
    mis["ref_category"] = ref_cat[mis_mask]
    mis["predicted_category"] = pred_cat[mis_mask]
    mis["mode_index_numeric"] = pd.to_numeric(mis["mode_index"], errors="coerce")
    mis = mis.merge(
        cm, how="left", left_on=["molecule", "mode_index_numeric"],
        right_on=["molecule", "mode"],
    ).drop(columns=["mode_index_numeric", "mode"])
    mis = mis.sort_values(
        ["ref_category", "predicted_category", "molecule", "mode_index"]
    ).reset_index(drop=True)
    mis.to_csv(misclassified_csv_path, index=False)
    n_misclassified = len(mis)

    per_molecule_n = df.groupby("molecule").size().to_dict()

    summary = {
        "pdf": pdf_path, "png": png_path,
        "matrix_csv": matrix_csv_path, "summary_csv": summary_csv_path,
        "misclassified_csv": misclassified_csv_path,
        "n_total": n_total,
        "n_misclassified": n_misclassified,
        "n_molecules": len(test_mols),
        "per_molecule_n": per_molecule_n,
        "n_dropped_missing_ground_truth": n_dropped,
        "n_folded_mixed_external": n_folded_mixed_external,
        "confusion_table": tbl.to_dict(),
        "precision": {r["category"]: r["precision"] for r in rows},
        "recall": {r["category"]: r["recall"] for r in rows},
        "framing": ("Genuinely non-circular: the 18-molecule test tier is "
                     "excluded from tau_S/tau_B calibration by construction "
                     "(SINGLE_CENTRE_ONLY_EXCLUDE), so this is the classifier "
                     "applied out-of-sample against independent hand-curated "
                     "ground truth (data/characterised_modes.csv), not a "
                     "restatement of the calibration population like "
                     "plot_rigorous_tier_check."),
    }
    return summary


def plot_transferability_confusion_binary(
    library_csv="data/results/library_scores.csv",
    out_dir="data/figures",
    tau_SB=None,
    label="fig_transferability_confusion",
    matrix_csv_path="data/results/transferability_confusion_matrix.csv",
    summary_csv_path="data/results/transferability_confusion_summary.csv",
    misclassified_csv_path="data/results/transferability_confusion_misclassified.csv",
):
    """Binary-scheme sibling of plot_transferability_confusion
    (fig:transferabilityconfusion, threeway version now at
    ``fig_transferability_confusion_threeway``) -- CANONICAL as of the
    2026-08 binary-classification switch (this function now owns the
    ``fig_transferability_confusion`` name/output, the paper's default
    scheme). The same 18-molecule ``mol_type=='test'``
    tier, but a 4x4 joint confusion matrix (T, R, B, S) -- literal literature
    "SB" reference rows are DROPPED entirely (a different question than a
    binary S/B classifier answers -- "is this mode genuinely mixed?" -- not
    folded into a bucket the way the three-way figure's "SB" column does),
    and every internal prediction goes through
    ``vib_label_binary(V_Stretch, tau_SB)`` instead of the three-way
    ``predicted_label``/``classification_bucket`` vocabulary.

    T/R ground truth/prediction: identical to the three-way figure
    (``external_axis()``, first character). tau_SB defaults to
    ``Thresholds.calibrated().tau_SB``. Same three output CSVs pattern as
    the original (matrix/summary/misclassified), suffixed ``_binary``.
    """
    _style()
    from src.csv_label_ingest import load_label_csvs
    from src.library_ingest import test_tier_molecules

    lib_df = pd.read_csv(library_csv)
    test_mols = test_tier_molecules()
    df = lib_df[lib_df["molecule"].isin(test_mols)].copy()
    tau_SB = tau_SB if tau_SB is not None else Thresholds.calibrated().tau_SB

    ext_mask = df["kind"] == "external"
    # No "SB" entry (unlike the three-way figure's ref_internal_map) -- a
    # literal literature "SB" row maps to NaN here and is dropped by `keep`
    # below, not folded into any bucket.
    ref_internal_map = {"bend": "B", "stretch": "S"}
    ref_cat = pd.Series(index=df.index, dtype=object)
    ref_cat[ext_mask] = df.loc[ext_mask, "mode_index"].astype(str).str[0]
    ref_cat[~ext_mask] = df.loc[~ext_mask, "ref_label"].map(ref_internal_map)

    def _pred_category(row):
        axis = external_axis(row["predicted_label"])
        if axis is not None:
            return axis[0]
        # vib_label_binary already returns the bare "S"/"B" letter.
        return vib_label_binary(row["V_Stretch"], tau_SB)

    pred_cat = df.apply(_pred_category, axis=1)

    keep = ref_cat.notna()
    n_dropped_sb_ref = int((~ext_mask & (df["ref_label"] == "SB")).sum())
    n_dropped_missing_ground_truth = int((~keep).sum()) - n_dropped_sb_ref
    n_folded_mixed_external = int(
        df["predicted_label"].astype(str).str.endswith("*").sum())

    ref_order = ["T", "R", "B", "S"]
    tbl = pd.crosstab(
        pd.Series(ref_cat[keep], name="ref"),
        pd.Series(pred_cat[keep], name="pred"),
    )
    tbl = tbl.reindex(index=ref_order, columns=ref_order, fill_value=0)
    n_total = int(tbl.values.sum())

    fig, ax = plt.subplots(figsize=(4.2, 3.8))
    _confusion_heatmap(ax, fig, tbl, ref_order)
    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    os.makedirs(os.path.dirname(matrix_csv_path), exist_ok=True)
    tbl.to_csv(matrix_csv_path)
    rows = _per_category_from_table(tbl, ref_order)
    pd.DataFrame(rows).to_csv(summary_csv_path, index=False)

    cm = load_label_csvs()["characterised_modes"][
        ["molecule", "mode", "description", "ref"]].copy()
    cm["mode"] = pd.to_numeric(cm["mode"], errors="coerce")
    mis_mask = keep & (ref_cat != pred_cat)
    mis = df.loc[mis_mask, ["molecule", "mode_index", "kind", "freq",
                             "irrep", "V_Stretch", "predicted_label"]].copy()
    mis["ref_category"] = ref_cat[mis_mask]
    mis["predicted_category"] = pred_cat[mis_mask]
    mis["mode_index_numeric"] = pd.to_numeric(mis["mode_index"], errors="coerce")
    mis = mis.merge(
        cm, how="left", left_on=["molecule", "mode_index_numeric"],
        right_on=["molecule", "mode"],
    ).drop(columns=["mode_index_numeric", "mode"])
    mis = mis.sort_values(
        ["ref_category", "predicted_category", "molecule", "mode_index"]
    ).reset_index(drop=True)
    mis.to_csv(misclassified_csv_path, index=False)
    n_misclassified = len(mis)

    per_molecule_n = df.groupby("molecule").size().to_dict()

    summary = {
        "pdf": pdf_path, "png": png_path,
        "tau_SB": tau_SB,
        "matrix_csv": matrix_csv_path, "summary_csv": summary_csv_path,
        "misclassified_csv": misclassified_csv_path,
        "n_total": n_total,
        "n_misclassified": n_misclassified,
        "n_molecules": len(test_mols),
        "per_molecule_n": per_molecule_n,
        "n_dropped_sb_ref": n_dropped_sb_ref,
        "n_dropped_missing_ground_truth": n_dropped_missing_ground_truth,
        "n_folded_mixed_external": n_folded_mixed_external,
        "confusion_table": tbl.to_dict(),
        "precision": {r["category"]: r["precision"] for r in rows},
        "recall": {r["category"]: r["recall"] for r in rows},
        "framing": ("Binary sibling of fig:transferabilityconfusion -- same "
                     "genuinely non-circular 18-molecule test tier, but "
                     "literal literature 'SB' reference rows are DROPPED "
                     "(not bucketed) and internal predictions come from "
                     "vib_label_binary(V_Stretch, tau_SB), never 'SB'."),
    }
    return summary


# --------------------------------------------------------------------------
# fig:benzeneconfusion -- benzene's genuine 3-class (bend/stretch/SB)
# internal-mode confusion matrix, enabled by the 2026-07-05 literature
# relabeling of modes 21/22 as a literal "SB" (mixed) ground truth. Does not
# modify plot_confusion_matrix (fig:confusion) or any of its inputs.
# --------------------------------------------------------------------------

def plot_benzene_internal_confusion(
    out_dir="data/figures",
    label="fig_benzene_confusion_threeway",
):
    """Three-way sibling of ``plot_benzene_internal_confusion_binary`` (which
    now owns the canonical ``fig_benzene_confusion`` name/output, the paper's
    default scheme as of the 2026-08 binary-classification switch) -- kept
    fully functional for comparison/on-demand use, not deleted. Builds
    benzene's JOINT confusion matrix -- the 6 external (Tx..Rz) normal modes
    together with the 30 internal modes' reference bend/stretch/SB (literal
    literature label, modes 21/22) x predicted bend/stretch/mixed bucket, in
    one crosstab. Joint so the manuscript's "recovered exactly (6/6)" claim
    is a checked matrix cell, not just a sentence -- off-diagonal
    external<->internal blocks
    (``crossover_ext_ref_to_internal_pred``/``crossover_internal_ref_to_ext_pred``)
    are computed and reported, not assumed zero.

    Data source: ``src.benzene_validation.benzene_normal_reference_detail(
    scheme="threeway")`` (all 36 rows). The explicit ``scheme="threeway"``
    is required here: that function's own default (``scheme=None``) tracks
    whatever scheme ``library_scores.csv`` was built under (binary, the
    global default as of the 2026-08 switch), which never produces "mixed"
    -- this figure's entire point is a genuine 3-class breakdown, so it
    must request the three-way re-derivation explicitly (see that function's
    docstring for the cheap V_Stretch-based rescheme mechanism). Distinct
    from fig:confusion, which is the whole-hydride-library validation; this
    is benzene's own 36 modes with its genuine literature 3-class internal
    ground truth. Companion per-bond detail for the 21/22 vs. 23/24 contrast
    lives in ``data/results/benzene_sb_vs_stretch_bond_diagnostic.csv``, not
    re-plotted here.
    """
    _style()
    from src.benzene_validation import benzene_normal_reference_detail

    detail = benzene_normal_reference_detail(scheme="threeway")

    # Collapsed to a single "T/R" row/column: per-axis detail is already
    # shown in fig:confusion's non-ideal matrix. "mixed_external" dropped
    # (verified empty on benzene's 6 external modes).
    ref_order = ["T/R", "bend", "SB", "stretch"]
    pred_order = ["T/R", "bend", "mixed", "stretch"]
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
                 "to a single T/R row/column (per-axis detail lives in "
                 "fig:confusion instead), turning the prose's '6/6 recovered "
                 "exactly' claim into an actual checked matrix cell without "
                 "repeating fig:confusion's per-axis breakdown. No divider "
                 "line and no mixed_external "
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


def plot_benzene_internal_confusion_binary(
    library_csv="data/results/library_scores.csv",
    out_dir="data/figures",
    label="fig_benzene_confusion",
    tau_SB=None,
):
    """Binary-scheme sibling of plot_benzene_internal_confusion
    (fig:benzeneconfusion, threeway version now at
    ``fig_benzene_confusion_threeway``) -- CANONICAL as of the 2026-08
    binary-classification switch (this function now owns the
    ``fig_benzene_confusion`` name/output, the paper's default scheme).
    Benzene's 6 external (Tx..Rz, collapsed to one
    T/R row/column) + internal modes, but:

    - internal reference rows literally labeled "SB" (modes 21/22) are
      DROPPED (n_dropped_sb_ref), not folded into a bucket -- a binary
      classifier has no "mixed" answer to compare a genuine literature "is
      this mode mixed?" label against, same convention as
      plot_confusion_matrix_binary/plot_transferability_confusion_binary.
    - every remaining internal row's predicted category comes from
      vib_label_binary(V_Stretch, tau_SB) instead of the three-way
      predicted_label, so the internal side of the matrix is 2-way
      (bend/stretch only, never "mixed").

    tau_SB defaults to Thresholds.calibrated().tau_SB. Reads
    library_scores.csv directly (unlike plot_benzene_internal_confusion,
    which uses benzene_normal_reference_detail() -- that helper doesn't
    carry V_Stretch, which this figure needs) -- self-contained, doesn't
    depend on plot_benzene_internal_confusion having run.
    """
    _style()
    from src.benzene_validation import MOLECULE

    lib_df = pd.read_csv(library_csv)
    b = lib_df[(lib_df["molecule"] == MOLECULE) & lib_df["ref_label"].notna()].copy()
    tau_SB = tau_SB if tau_SB is not None else Thresholds.calibrated().tau_SB

    n_dropped_sb_ref = int((b["ref_label"] == "SB").sum())
    b = b[b["ref_label"] != "SB"]

    def _binarize(row):
        if row["kind"] != "internal":
            return row["predicted_label"]
        bucket = classification_bucket(row["predicted_label"])
        if bucket in ("stretch", "bend", "mixed"):
            return vib_label_binary(row["V_Stretch"], tau_SB)
        return row["predicted_label"]  # external crossover -- left untouched

    b["predicted_label"] = b.apply(_binarize, axis=1)

    ref_order = ["T/R", "bend", "stretch"]
    pred_order = ["T/R", "bend", "stretch"]
    tbl = _joint_confusion_table(b, ref_order, pred_order, collapse_external=True)
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
        "tau_SB": tau_SB,
        "shared_categories": ("reuses CATEGORY_COLOR/CATEGORY_LABEL/"
                               "_confusion_heatmap from plot_confusion_matrix "
                               "(fig:confusion) -- does not modify "
                               "fig:benzeneconfusion itself."),
        "n_total": n_total,
        "n_dropped_sb_ref": n_dropped_sb_ref,
        "confusion_table": tbl.to_dict(),
        "crossover_ext_ref_to_internal_pred": crossover_ext_ref_to_internal_pred,
        "crossover_internal_ref_to_ext_pred": crossover_internal_ref_to_ext_pred,
        "note": (f"3x3 joint matrix (binary scheme, tau_SB={tau_SB:g}): benzene's "
                 "6 external T/R modes collapsed to one T/R row/column (per-axis "
                 "detail lives in fig:confusion instead); the 2 literal literature "
                 "'SB' reference rows (modes 21/22) are dropped rather than folded "
                 "into a bucket. Binary sibling of fig:benzeneconfusion -- see that "
                 "figure for the three-way (S/B/SB) version."),
    }
    return summary_dict


def plot_benzene_confusion_precision_recall(
    matrix_csv="data/results/benzene_internal_confusion_matrix.csv",
    summary_csv="data/results/benzene_internal_confusion_summary.csv",
    out_dir="data/figures",
    label="fig_benzene_precision_recall",
):
    """Build the SI companion to fig:benzeneconfusion: benzene's internal-
    mode per-category precision/recall bar panel, split into the SI since a
    formal accuracy-metric bar chart sits awkwardly next to this section's
    "mixed by eye is a convention" hedging (benzene's literature ground
    truth is genuinely non-circular, so this is a tone choice, not a
    circularity concern -- unlike ``plot_rigorous_tier_check``).

    Layout (1x1): grouped precision/recall bars for bend/stretch/"mixed"
    (the predicted-bucket name for the literal "SB" reference row).
    Precision is computed here from the confusion table; recall comes from
    ``summary_csv``'s own ``recall`` column (not re-derived).
    """
    _style()
    tbl = pd.read_csv(matrix_csv, index_col=0)
    summary = pd.read_csv(summary_csv).set_index("ref_label")

    ref_order = ["bend", "SB", "stretch"]
    pred_order = ["bend", "mixed", "stretch"]
    tbl = tbl.reindex(index=ref_order, columns=pred_order, fill_value=0)

    # Each category pairs one reference row with one predicted column:
    # bend<->bend, stretch<->stretch, SB(reference)<->mixed(predicted).
    cats = ["bend", "mixed", "stretch"]
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
# fig:benzeneemitcounts -- benzene EMIT mode classification overview
# --------------------------------------------------------------------------

def plot_benzene_emit_counts(
    emit_csv="data/results/C6H6_EMIT.csv",
    out_dir="data/figures",
    label="fig_benzene_emit_counts",
):
    """Build fig:benzeneemitcounts: an overview bar chart of how benzene's
    36 EMIT modes classify -- one bar per clean/starred external axis
    (Tx..Rz, Tx*..Rz*) and per internal bucket (S/B/SB).

    Presentation only: reads the already-computed ``label`` column of raw
    classifier-output strings from ``emit_csv`` via ``value_counts()`` --
    never recomputes scores. Categories plotted are exactly whatever is
    present in the data (not a hardcoded set), so the chart stays correct
    if the classification changes as more EMIT diagnostics land.
    """
    _style()
    df = pd.read_csv(emit_csv)
    counts = df["label"].value_counts()

    # Canonical order: T-axis, then R-axis (clean before starred within an
    # axis), then internal B/SB/S in ascending-frequency order (bending is
    # lowest frequency, stretching highest) -- filtered to categories present.
    canonical_order = [
        "Tx", "Tx*", "Ty", "Ty*", "Tz", "Tz*",
        "Rx", "Rx*", "Ry", "Ry*", "Rz", "Rz*",
        "B", "SB", "S",
    ]
    present = [c for c in canonical_order if c in counts.index]
    unexpected = [c for c in counts.index if c not in canonical_order]
    if unexpected:
        raise ValueError(
            f"benzene EMIT label(s) not in the canonical T/R/S/B/SB order: "
            f"{unexpected} -- extend canonical_order, don't silently drop them."
        )

    # T/R (clean or starred) share one gray; S/B/SB use the classification-
    # algorithm's own predicted-label colors (not the *_ref literature pair).
    internal_color = {
        "S": COLORS["stretching"], "B": COLORS["bending"], "SB": COLORS["mixed"],
    }
    bar_colors = [internal_color.get(c, COLORS["external"]) for c in present]
    bar_counts = [int(counts[c]) for c in present]

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    x = np.arange(len(present))
    ax.bar(x, bar_counts, color=bar_colors)
    ax.set_xticks(x)
    ax.set_xticklabels(present)  # matches canonical_order left-to-right
    ax.set_ylabel("Number of EMIT modes")
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.set_ylim(0, max(bar_counts) * 1.18)
    for xi, n in zip(x, bar_counts):
        ax.text(xi, n + max(bar_counts) * 0.02, str(n), va="bottom", ha="center",
                 fontsize=ANNOTATION_FONTSIZE)

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    n_total = int(sum(bar_counts))
    if n_total != len(df):
        raise ValueError(
            f"plotted EMIT mode count {n_total} != len(df) {len(df)} -- "
            "a category/counting bug in this function, not a data change."
        )

    summary_dict = {
        "pdf": pdf_path, "png": png_path,
        "n_total": n_total,
        "counts": {c: int(counts[c]) for c in present},
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
    stretching (vermillion) vs. bending (blue). One consistent marker shape
    (circle) for all points -- color already fully distinguishes stretch/
    bend, so shape would be redundant.

    Restricted to the single-centre AB_n hydride-library scope
    (src.calibrate.SINGLE_CENTRE_ONLY_EXCLUDE), applied right after reading
    the CSV, before any bond is exploded.
    """
    _style()
    from src.calibrate import filter_single_centre_library

    lib_df = pd.read_csv(library_csv)
    lib_df = filter_single_centre_library(lib_df)
    bonds = _explode_bonds(lib_df)
    # s_AB is now signed (positive = stretching, negative = compressing) in
    # the underlying data; this figure plots its magnitude only (unchanged
    # appearance from before s_AB became signed).
    bonds["s_AB"] = bonds["s_AB"].abs()
    bonds = bonds[bonds["ref_label"].isin(("stretch", "bend")) & bonds["ideal"].isin(("yes", "no"))]
    bonds["abs_rel_db"] = bonds["rel_db"].abs()

    # Width fixed at 6.2in to match every other main-text single-panel
    # figure's printed size under LaTeX's \includegraphics rescaling --
    # a wider source width shrinks the rcParams fonts below the rest of
    # the figure set. Don't reduce below this without re-checking that.
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

    # matplotlib's mathtext supports \boldsymbol, matching the manuscript
    # prose's own notation for the per-bond |Delta|b||/|b| ratio.
    ax.set_xlabel(r"$|\Delta|\boldsymbol{b}|\,/\,|\boldsymbol{b}|\,|$")
    ax.set_ylabel(r"$|s^{AB}|$")
    ax.set_xlim(-0.03, bonds["abs_rel_db"].max() * 1.05)
    ax.set_ylim(-0.03, 1.05)

    # One consistent marker shape (circle); color + fill are the only two
    # encodings. Custom-built (not sourced from CATEGORY_LABEL, since it
    # also encodes ideal/non-ideal) using the same bare S/B wording.
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
    stretch/ideal, bend/non-ideal, stretch/non-ideal (filled=ideal,
    hollow=non-ideal, IDEAL_STYLE) -- over every internal row in the
    library with a literature stretch/bend label.

    Restricted to the single-centre AB_n hydride-library scope (see
    plot_bond_scores' identical note), applied right after reading the CSV.
    """
    _style()
    from src.calibrate import filter_single_centre_library

    lib_df = pd.read_csv(library_csv)
    lib_df = filter_single_centre_library(lib_df)
    internal = lib_df[(lib_df["kind"] == "internal") &
                       lib_df["ref_label"].isin(("stretch", "bend")) &
                       lib_df["ideal"].isin(("yes", "no"))].copy()

    groups = [("bend", "yes"), ("stretch", "yes"), ("bend", "no"), ("stretch", "no")]
    # Bare S/B tick labels, matching CATEGORY_LABEL's short-notation
    # vocabulary; the ideal/non-ideal distinction is stated in the caption
    # text instead of drawn on the raster.
    group_labels = ["B", "S", "B", "S"]
    # Gap 1.3 within a bend/stretch pair, gap 1.6 between the ideal pair
    # (1,2) and non-ideal pair (3,4) -- tuned so tick labels don't collide,
    # verified by rendering. Don't shrink without re-checking.
    positions = [1.0, 2.3, 3.9, 5.2]

    # Bare (a)/(b)/(c) labels only -- each panel's y-axis label already
    # carries the descriptive title text.
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
        ax.set_xticklabels(group_labels)
        ax.set_xlim(positions[0] - 0.7, positions[-1] + 0.7)
        ax.set_ylabel(ylabel)
        ax.set_title(title, loc="center", fontweight="bold")

    # tau_SB is not drawn on panel (c) (visual clutter); the value is still
    # returned below for the caption/summary.
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
        "tau_SB": th.tau_SB,
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
    tau_SB=None,
):
    """Build fig:modemixing: (a) ideal molecules -- V_Stretch vs.
    mode-averaged |Delta b|/|b|| shows a clean step function; (b) non-ideal
    molecules -- the same axes show a graded transition. Both panels reuse
    the shared STRETCHING/BENDING category colors; ideal (a) is all filled,
    non-ideal (b) is all hollow (IDEAL_STYLE). One consistent marker shape
    (circle) -- color already fully distinguishes stretch/bend.

    NOTE (pending gap): this renders only the two-panel ideal-step-vs-non-
    ideal-gradient content the current .tex caption (fig:modemixing) describes.
    The irrep-degeneracy sub-panels (scoped-doc Fig 3c-d in the group's
    earlier report, keyed on trigonal-planar/bent-AB2 irreps there, prose
    mentions ethane here) are NOT built -- the intended molecule/panel form
    for THIS manuscript is still unconfirmed with lead-author/tex-data-sync,
    so nothing is invented for it.

    Restricted to the single-centre AB_n hydride-library scope (2026-07-05
    author decision -- see src.calibrate.SINGLE_CENTRE_ONLY_EXCLUDE /
    plot_bond_scores' identical note); applied right after reading the CSV,
    before the internal-row selection below.

    ``tau_SB`` defaults to ``Thresholds.calibrated().tau_SB`` (None ->
    canonical, zero behavior change). Only the drawn threshold line/label use
    the override -- point colors come from ``ref_label`` (literature ground
    truth), never from tau_SB, so they are untouched either way.
    """
    _style()
    from src.calibrate import filter_single_centre_library

    lib_df = pd.read_csv(library_csv)
    lib_df = filter_single_centre_library(lib_df)
    internal = lib_df[(lib_df["kind"] == "internal") &
                       lib_df["ref_label"].isin(("stretch", "bend")) &
                       lib_df["ideal"].isin(("yes", "no"))].copy()

    th = Thresholds.calibrated()
    if tau_SB is not None:
        th = dataclasses.replace(th, tau_SB=tau_SB)
    fig, (ax_i, ax_n) = plt.subplots(1, 2, figsize=(7.2, 3.3), sharex=True, sharey=True)

    # Bare (a)/(b) labels only -- "Ideal/Non-ideal molecules" is stated in
    # the caption, not the in-figure title.
    for ax, ideal_flag, title in ((ax_i, "yes", "(a)"),
                                   (ax_n, "no", "(b)")):
        for ref in ("bend", "stretch"):
            cat = REF_LABEL_TO_CATEGORY[ref]
            sub = internal[(internal.ideal == ideal_flag) & (internal.ref_label == ref)]
            if sub.empty:
                continue
            kw = _marker_kwargs(cat, ideal_flag, marker="o", color_map=REF_CATEGORY_COLOR)
            # Each panel's point label spells out ideal/non-ideal explicitly.
            ideal_suffix = "ideal" if ideal_flag == "yes" else "non-ideal"
            point_label = f"{CATEGORY_LABEL[cat]}, {ideal_suffix}"
            ax.scatter(sub["delta_b_mean"], sub["V_Stretch"], s=20,
                       label=point_label, **kw)
        ax.axhline(th.tau_SB, color=COLORS["threshold"], ls="--", lw=0.8, zorder=1)
        ax.set_xlabel(r"Averaged $|\Delta|\boldsymbol{b}|\,/\,|\boldsymbol{b}|\,|$")
        ax.set_title(title, loc="center", fontweight="bold")
        ax.set_xlim(-0.03, internal["delta_b_mean"].max() * 1.08)

    ax_i.set_ylabel(r"$s[\mathrm{V_S}]$")
    ax_i.set_ylim(-0.05, 1.08)
    # Bordered legends: both panels' legend boxes sit on top of scattered
    # data points, so a frame + opaque face keep swatches from reading as
    # more data (unlike every other legend in this module, frameon=False).
    ax_i.legend(loc="center right", frameon=True, edgecolor="black",
                facecolor="white", framealpha=1.0,
                handletextpad=0.4, labelspacing=0.4)
    ax_n.legend(loc="center right", frameon=True, edgecolor="black",
                facecolor="white", framealpha=1.0,
                handletextpad=0.4, labelspacing=0.4)
    # Left side, clear of the "center right" legend.
    ax_i.text(0.02, th.tau_SB, r"$\tau_{\mathrm{SB}}$", ha="left", va="bottom",
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
        "tau_SB": th.tau_SB,
        "irrep_degeneracy_panel": ("NOT BUILT in this figure -- see the "
                                   "separate SI figure plot_irrep_coupling "
                                   "(fig_irrep_coupling) below (AB3 "
                                   "trigonal-planar / AB2 bent series)."),
    }
    return summary


# --------------------------------------------------------------------------
# SI: irrep-degeneracy mixing mechanism (AB3 trigonal-planar / AB2 bent
# series). Supports the manuscript's prose claim that mixing occurs only
# where modes of the same irrep can couple: mode score vs. central-atom
# displacement amplitude |d_CA|, faceted by whether a mode's irrep has a
# same-irrep coupling partner within its point group.
#
#   - AB3 (D3h): A2'' (bend) and A1' (stretch) are each the only mode of
#     that label -- no partner, stay flat/clean regardless of |d_CA|. E'
#     appears twice (bend and stretch) and CAN couple; only stretch-E'
#     visibly droops from ~1 as |d_CA| grows.
#   - AB2 (C2v): A1 appears twice (bend and stretch) and both branches show
#     drooping/mixing with |d_CA|. B2 (antisymmetric stretch) is unique and
#     stays comparatively clean.
#
# Data source: ``data/characterised_modes.csv`` for irrep/shape/type,
# ``mol_list_method.csv``'s `mol_type` for the ideal/non-ideal filter, and
# ``library_scores.csv`` for V_Stretch and `d_CA`. Merged on (molecule,
# mode) -- both files share the same 1-based internal-mode index.
#
# DATA-QUALITY GOTCHA: some molecules (e.g. OCl2) have an ASCII-vs-Unicode-
# subscript irrep-string mismatch that silently plots zero points for them
# rather than erroring (BBr3 similarly mismatches on its `shape` field) --
# see src/library_ingest.py's regenerate_characterised_modes() docstring for
# the full reasoning/carve-out. Flagged for the author, not auto-corrected.
# --------------------------------------------------------------------------

def _hollow_marker_kwargs(color_key, marker):
    """Marker styling for one (type, irrep) category in fig_irrep_coupling:
    color encodes the literature bend/stretch label (REFERENCE hues, via
    `color_key`, not a classifier prediction); marker SHAPE encodes irrep
    identity, with the same shape reused across the bend/stretch color
    split whenever that irrep is shared (coupling-capable) between the two
    -- shape-matching across a color change is the visual signal of a
    symmetry-permitted coupling pathway. All markers are hollow (fill
    carries no information here).
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
    planar AB3 series, (b) bent AB2 series. Both panels plot mode score
    (``V_Stretch``) vs. central-atom displacement amplitude (``d_CA``),
    faceted by (type, irrep): all markers hollow; color is bend/stretch
    (reference hues, see ``_hollow_marker_kwargs``); marker SHAPE encodes
    irrep identity, with the same shape reused across the bend/stretch
    color split when that irrep is shared between the two (coupling-
    capable) -- shape-matching across a color change is the visual cue for
    a symmetry-permitted coupling pathway.

    Some newly-back-filled molecules (e.g. BBr3, OCl2) are shape-eligible
    but contribute zero rendered points due to the irrep/shape-string data
    mismatch documented in this section's header comment above -- the
    returned summary's ``ab3_shape_eligible_but_not_plotted``/
    ``ab2_shape_eligible_but_not_plotted`` keys surface this explicitly
    (``ab3_molecules``/``ab2_molecules`` report only molecules with >=1
    rendered point, not just shape-eligibility).

    Not one of the 6 originally-scoped ``fig:*`` labels -- an SI figure
    filling the irrep-degeneracy sub-panel gap flagged against
    ``fig:modemixing``.
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
    # cm['type'] is the abbreviated "B"/"S"/"SB" storage form; normalize to
    # the full words the categories tuples below match against.
    df["type"] = df["type"].map(TYPE_TO_REF_LABEL).fillna(df["type"])
    df = df[df["molecule"].isin(nonideal_molecules)]
    df = df.dropna(subset=["|d_CA|"])  # multi-centre/no-central-atom rows

    ab3 = df[df["shape"] == "trigonal planar"].copy()
    ab2 = df[df["shape"] == "bend"].copy()

    # (type, irrep) -> (color_key, marker shape, legend text). Defined
    # per-panel since the two point groups have different shared/unique
    # irrep assignments (AB2's A1 is shared in both branches; AB3's
    # A1'/A2'' are each unique while E' is shared). Marker shape is keyed
    # to the irrep symbol so the same irrep gets the same shape in both its
    # bend and stretch appearances. "type" is the literal literature label,
    # so color_key uses the "_ref" hues, not the algorithm's own colors.
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

    # Width tuned to 9.0in -- narrower overflows the per-panel legend text
    # into the neighboring panel and clips the x-axis label. Don't reduce.
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

    # Bare (a)/(b) labels only -- full titles are stated in the caption.
    n_ab3, ab3_plotted = _plot_panel(ax_3, ab3, ab3_categories, "(a)")
    n_ab2, ab2_plotted = _plot_panel(ax_2, ab2, ab2_categories, "(b)")
    ax_3.set_ylabel(r"Mode score $s[\mathrm{V}]$")

    # Molecules whose `shape` matched the panel but contributed zero
    # rendered points (irrep/shape-string mismatch, see gotcha above).
    # Reported explicitly rather than silently folded into
    # ab3_molecules/ab2_molecules below (shape-eligibility != rendering).
    ab3_shape_only = sorted(set(ab3["molecule"].unique()) - ab3_plotted)
    ab2_shape_only = sorted(set(ab2["molecule"].unique()) - ab2_plotted)

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary = {
        "pdf": pdf_path, "png": png_path,
        "shared_categories": ("reuses CATEGORY_COLOR's bending/stretching "
                               "hues; ALL markers unfilled/hollow (fill does "
                               "not encode anything in this figure). Marker "
                               "SHAPE encodes irrep identity: "
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

    # Height tuned to 4.6in: both y-axis labels are long rotated-90deg
    # strings that clip at the canvas top/bottom at a shorter height
    # (width alone doesn't fix it -- tight_layout doesn't reserve space
    # for a twinx() secondary label). Don't reduce below this.
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


def plot_tau_sb_sensitivity(
    sweep_csv="data/results/tau_sb_sensitivity_sweep.csv",
    thresholds_json="data/results/thresholds.json",
    out_dir="data/figures",
    label="fig_sensitivity_binary",
):
    """Binary-scheme sibling of plot_sensitivity (fig:sensitivity): error_all,
    error_test, and error_single_centre vs. tau_SB over
    src.calibrate.sweep_tau_sb's persisted output
    (data/results/tau_sb_sensitivity_sweep.csv), marking the active tau_SB
    default (Thresholds.calibrated().tau_SB), the fixed tau_SB=0.5 reference
    line (the naive midpoint cutoff), and the three ADVISORY optimal points
    (thresholds.json's "tau_SB_error_sweep" block, written by
    src.calibrate.run_tau_sb_error_analysis -- see that function's docstring
    for why they are advisory-only and never silently adopted as the new
    default).

    Unlike plot_sensitivity's twin accuracy/label-change-fraction axes
    (different natural scales), error_all/error_test/error_single_centre
    share one [0, 1] scale, so this figure uses a single axis.
    """
    _style()
    sweep = pd.read_csv(sweep_csv)
    thresholds = Thresholds.calibrated(path=thresholds_json)
    tau_SB_active = thresholds.tau_SB

    with open(thresholds_json) as f:
        result = json.load(f)
    sb_sweep_info = result.get("tau_SB_error_sweep", {})
    opt_all = sb_sweep_info.get("optimal_tau_all")
    opt_test = sb_sweep_info.get("optimal_tau_test")
    opt_sc = sb_sweep_info.get("optimal_tau_single_centre")

    fig, ax1 = plt.subplots(figsize=(5.6, 4.6))
    ax1.axvline(tau_SB_active, color=COLORS["threshold"], ls="--", lw=1.2, zorder=2,
                label=r"active $\tau_{\mathrm{SB}}$" + f"={tau_SB_active:g}")
    ax1.axvline(0.5, color=COLORS["threshold"], ls="-.", lw=1.0, zorder=2,
                label=r"$\tau_{\mathrm{SB}}$=0.5")
    if opt_all is not None:
        ax1.axvline(opt_all, color=COLORS["sens_change"], ls=":", lw=1.0, zorder=2,
                    label=f"optimal (all)={opt_all:g}")
    if opt_test is not None:
        ax1.axvline(opt_test, color=COLORS["sens_accuracy"], ls=":", lw=1.0, zorder=2,
                    label=f"optimal (test)={opt_test:g}")
    if opt_sc is not None:
        ax1.axvline(opt_sc, color=COLORS["sens_single_centre"], ls=":", lw=1.0, zorder=2,
                    label=f"optimal (single-centre)={opt_sc:g}")

    l1, = ax1.plot(sweep["tau_SB"], sweep["error_all"], color=COLORS["sens_change"],
                   marker="o", markersize=2.5, lw=1.1, zorder=3, label="error (all)")
    l2, = ax1.plot(sweep["tau_SB"], sweep["error_test"], color=COLORS["sens_accuracy"],
                   marker="s", markersize=2.5, lw=1.1, zorder=3, label="error (test)")
    l3, = ax1.plot(sweep["tau_SB"], sweep["error_single_centre"],
                   color=COLORS["sens_single_centre"],
                   marker="^", markersize=2.5, lw=1.1, zorder=3,
                   label="error (single-centre)")
    ax1.set_xlabel(r"$\tau_{\mathrm{SB}}$")
    ax1.set_ylabel("Classification error")
    ax1.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1))
    error_cols = ["error_all", "error_test", "error_single_centre"]
    ax1.set_ylim(-0.01, max(0.05, float(sweep[error_cols].max().max()) * 1.1))
    ax1.legend(loc="best", frameon=False, fontsize=LEGEND_FONTSIZE,
               handletextpad=0.5, labelspacing=0.4)

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary = {
        "pdf": pdf_path, "png": png_path,
        "shared_categories": ("no classification-category colors here (tau_SB "
                               "error curves, not per-mode points); reuses the "
                               "shared COLORS['threshold'] dashed/dash-dot-line "
                               "convention for the active-tau and fixed-0.5 "
                               "markers."),
        "n_grid_points": len(sweep),
        "tau_SB_range": (float(sweep["tau_SB"].min()), float(sweep["tau_SB"].max())),
        "error_all_range": (float(sweep["error_all"].min()), float(sweep["error_all"].max())),
        "error_test_range": (float(sweep["error_test"].min()), float(sweep["error_test"].max())),
        "error_single_centre_range": (float(sweep["error_single_centre"].min()),
                                       float(sweep["error_single_centre"].max())),
        "tau_SB_active": tau_SB_active,
        "tau_SB_fixed_reference": 0.5,
        "optimal_tau_SB_all": opt_all,
        "optimal_tau_SB_test": opt_test,
        "optimal_tau_SB_single_centre": opt_sc,
    }
    return summary


# --------------------------------------------------------------------------
# fig:cputime -- empirical computational-cost figure (replaces/accompanies
# the purely theoretical Big-O argument in "Computational cost", tab:cost).
# PROPOSED LABEL, not yet wired into the .tex (that edit is a separate step
# -- see this module's own docstring: figures.py never touches the .tex).
# --------------------------------------------------------------------------

def _loglog_ols_fit(x, y):
    """Log-log OLS power-law fit y ~ a * x**b, fit as a straight line in
    log-log space (np.polyfit on ln(y) vs. ln(x)); return (slope, intercept,
    r2) with r2 computed on log-space residuals, matching the pattern used
    throughout this module (see plot_gaussian_nbasis_scaling)."""
    log_x, log_y = np.log(x), np.log(y)
    slope, intercept = np.polyfit(log_x, log_y, 1)
    pred = slope * log_x + intercept
    r2 = 1 - np.sum((log_y - pred) ** 2) / np.sum((log_y - log_y.mean()) ** 2)
    return slope, intercept, r2


def _powerlaw_curve_fit(x, y, p0):
    """Nonlinear least squares power-law fit y ~ a * x**b, minimizing raw
    (non-logged) residuals directly (scipy.optimize.curve_fit); return
    (a, b, r2) with r2 computed on linear-space residuals, matching the
    pattern used throughout this module (see plot_gaussian_nbasis_scaling)."""
    def _powerlaw(n, a, b):
        return a * n ** b
    (a, b), _ = curve_fit(_powerlaw, x, y, p0=p0)
    pred = _powerlaw(x, a, b)
    r2 = 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)
    return a, b, r2


def plot_cpu_time_benchmark(
    benchmark_csv="data/results/cpu_time_benchmark.csv",
    out_dir="data/figures",
    label=None,
    scale="linear",
    gaussian_color_by_nbasis=True,
):
    """Build fig:cputime: empirical CPU time vs. atom count N, for the
    classification algorithm (this framework) against the Gaussian
    frequency-calculation step that supplies its input -- one point per
    hydride-library molecule, restricted to the MP2/3-21G subset
    (``mp2_321g == True``, 50 of 68 molecules, N=3..12; benzene/C6H6 at N=12
    is the single largest case in this filtered subset).

    ``scale`` picks the y-axis: "linear" (default, label "fig_cputime") or
    "log" (label "fig_cputime_log" unless overridden). Both are generated
    and shown side by side per standing instruction -- neither alone tells
    the full story. Linear collapses the classifier series to ~0 next to
    Gaussian's, which dramatizes the magnitude gap (classification cost is
    negligible against the frequency calculation that supplies its input)
    but hides the classifier's own N-scaling. Log keeps both series legible
    and closer to showing each series' own N-scaling trend, but undersells
    the magnitude gap. ``label`` overrides the scale-based default if given.

    FITTING: each series (Gaussian, classifier) gets a power-law fit against
    N, but -- unlike plot_gaussian_nbasis_scaling, which draws both a
    log-log OLS fit and a linear-space nonlinear fit together on its linear
    variant -- fig_cputime draws only the fit appropriate to the axes it's
    shown on: the "log" variant draws the log-log OLS fit (``np.polyfit`` on
    ln(t) vs. ln(N)), the conventional way to report a scaling exponent and
    a straight line on log-log axes; the "linear" (default) variant instead
    draws the linear-space nonlinear fit (``scipy.curve_fit``, minimizing
    actual CPU-second residuals), the curve that tracks absolute vertical
    distance on a linear panel. Both fits are always computed regardless of
    ``scale`` and both are returned in the summary dict (``fit_vs_N_loglog_*``
    / ``fit_vs_N_linear_space_*``) even though only one is drawn per call --
    see plot_gaussian_nbasis_scaling for why the two typically disagree.
    Fits use raw per-molecule points, not per-N means -- N=4 (24 molecules)
    dominates each regression, reflecting where the sampling effort
    actually went rather than treating every N as equally well-supported.
    These fit lines replace the previous per-N mean line (a plain
    descriptive average, not a fit); that per-N mean is still computed and
    returned in the summary dict (not drawn) for reference.

    FILTERING: ``cpu_time_benchmark.csv`` carries ``mp2_321g``/
    ``method_basis`` joined from ``data/mol_list_method.csv`` --
    18 of the 68 library molecules were run at a different Gaussian
    method/basis (e.g. MP2/6-311G) for SCF-convergence or symmetry reasons,
    confounding a bare "CPU time vs. N" reading with "CPU time vs. method".
    This function filters to ``mp2_321g == True`` BEFORE plotting or
    computing summary statistics, so both the figure and every number in the
    returned summary describe the controlled, single-method/basis subset
    only. The excluded 18 (including AsBr3/MP2-6-311G, the single largest
    CPU-time point in the full 68-molecule set at 1855 s) are NOT shown and
    NOT folded into any statistic here -- see ``n_excluded_non_mp2_321g``.

    ``gaussian_freq_cpu_s`` is Gaussian's frequency-only CPU time
    (optimization excluded, verified via the Link1 job-step split in each
    .log -- see scripts/benchmark_cpu_time.py); ``classifier_cpu_s`` is this
    framework's classification CPU time (time.process_time(), timeit-
    calibrated against ``n_iterations`` to survive Windows' ~15.6 ms OS-tick
    granularity), with ``classifier_cpu_s_stddev`` plotted as a thin error
    bar (invisible at this linear scale/span -- included anyway, not for
    visual effect).

    Small reproducible x-jitter (fixed seed) separates the many molecules
    sharing the same integer N -- N itself is not perturbed in the
    underlying data, only the plotted x-position.

    COLOR ENCODING: by default (``gaussian_color_by_nbasis=True``), Gaussian
    points are colored by a continuous "inferno" gradient keyed to
    ``n_basis`` (Gaussian's own AO basis-function count, ``NBasis=``), not a
    flat color -- N alone barely explains Gaussian's CPU time (see
    ``plot_gaussian_nbasis_scaling``: log-log R^2~0.04 vs. N, R^2~0.81 vs.
    n_basis), so the gradient gives a reader an at-a-glance reason for the
    vertical scatter within each N without requiring a separate figure,
    with a colorbar carrying the mapping; the legend's Gaussian swatch uses
    a mid-gradient tone with no "colored by..." qualifier in its label,
    since the colorbar alone carries that mapping. Passing
    ``gaussian_color_by_nbasis=False`` abandons this in favor of plain
    black Gaussian markers and drops the colorbar entirely (the figure also
    narrows to match, since there's no colorbar to leave room for) -- a
    simpler alternative when the n_basis story isn't the point of a
    particular use of this figure; the default label gains a ``_bw`` suffix
    (``fig_cputime_bw`` / ``fig_cputime_log_bw``) so both variants can be
    generated side by side without colliding. The classifier series is
    always Okabe-Ito bluish green (``COLORS["cost_classifier"]``), unchanged
    either way since its cost does not depend on n_basis. Marker SHAPE
    distinguishes the two series -- circles for Gaussian, squares for the
    classifier. Markers (and the colorbar, when present) are drawn without
    edge borders (author preference against overusing borders).

    Not one of the 6 named JCC figures in this module's existing scope (see
    module docstring / IMPLEMENTATION_PLAN.md); added because the
    "Computational cost" section (tab:cost, purely theoretical Big-O) is
    being rewritten around this empirical benchmark. Caption text and .tex
    wiring are a separate step -- this function only builds the artifact.
    """
    if scale not in ("linear", "log"):
        raise ValueError(f"scale must be 'linear' or 'log', got {scale!r}")
    if label is None:
        base_label = "fig_cputime" if scale == "linear" else "fig_cputime_log"
        label = base_label if gaussian_color_by_nbasis else f"{base_label}_bw"

    _style()
    df_all = pd.read_csv(benchmark_csv)
    n_total = len(df_all)
    df = df_all[df_all["mp2_321g"] == True].copy()  # noqa: E712 -- explicit bool filter, not truthiness
    n_excluded = n_total - len(df)
    df = df.sort_values(["N", "molecule"]).reset_index(drop=True)
    ratio = df["gaussian_freq_cpu_s"] / df["classifier_cpu_s"]

    N_arr = df["N"].to_numpy(dtype=float)
    gaussian_arr = df["gaussian_freq_cpu_s"].to_numpy(dtype=float)
    classifier_arr = df["classifier_cpu_s"].to_numpy(dtype=float)
    g_slope, g_intercept, g_r2 = _loglog_ols_fit(N_arr, gaussian_arr)
    c_slope, c_intercept, c_r2 = _loglog_ols_fit(N_arr, classifier_arr)
    g_a_lin, g_b_lin, g_r2_lin = _powerlaw_curve_fit(
        N_arr, gaussian_arr, p0=[np.exp(g_intercept), g_slope])
    c_a_lin, c_b_lin, c_r2_lin = _powerlaw_curve_fit(
        N_arr, classifier_arr, p0=[np.exp(c_intercept), c_slope])

    rng = np.random.default_rng(0)
    jitter = rng.uniform(-0.12, 0.12, size=len(df))
    x = df["N"].to_numpy(dtype=float) + jitter

    fig, ax = plt.subplots(figsize=(5.3, 3.9) if gaussian_color_by_nbasis else (4.6, 3.9))

    if gaussian_color_by_nbasis:
        # A dark, fixed point on the same inferno ramp used for the
        # Gaussian points -- ties the trend line to its own series' color
        # family instead of an unrelated black, while staying dark enough
        # (low end of inferno) to read clearly as a line rather than blend
        # into the marker cloud.
        gaussian_trend_color = matplotlib.colormaps["inferno"](0.25)
    else:
        gaussian_trend_color = "black"

    ax.errorbar(x, df["classifier_cpu_s"], yerr=df["classifier_cpu_s_stddev"],
                fmt="none", ecolor=COLORS["cost_classifier"], elinewidth=0.5,
                alpha=0.35, zorder=2, capsize=0)
    if gaussian_color_by_nbasis:
        n_basis_norm = mcolors.Normalize(vmin=df["n_basis"].min(), vmax=df["n_basis"].max())
        gaussian_pts = ax.scatter(x, df["gaussian_freq_cpu_s"], c=df["n_basis"],
                   cmap="inferno", norm=n_basis_norm, marker="o", s=20,
                   edgecolors="none", alpha=1.0, zorder=3)
        cbar = fig.colorbar(gaussian_pts, ax=ax, pad=0.02, fraction=0.06)
        cbar.set_label(r"$N_{\mathrm{basis}}$", fontsize=LEGEND_FONTSIZE)
        cbar.ax.tick_params(labelsize=LEGEND_FONTSIZE)
        cbar.outline.set_visible(False)
    else:
        # Flat black, no colorbar -- carries its own real legend label
        # directly, unlike the gradient version which needs a proxy handle
        # (see gaussian_handle below) since a legend swatch can't show a
        # continuous colormap.
        ax.scatter(x, df["gaussian_freq_cpu_s"], marker="o", s=20,
                   facecolors="black", edgecolors="none", alpha=1.0, zorder=3,
                   label="Freq=hpmodes (MP2/3-21G)")
    ax.scatter(x, df["classifier_cpu_s"], marker="s", s=20,
               facecolors=COLORS["cost_classifier"],
               edgecolors="none", alpha=1.0, zorder=3,
               label="Classification algorithm")

    # Per-N mean (unjittered, true N on the x-axis) -- kept only as a
    # descriptive stat in the summary dict (see per_n_mean below), no longer
    # drawn: a plain average across a wildly uneven per-N sample count
    # (N=4: 24 molecules, N=12: 1 molecule) isn't a fit, and the power-law
    # fit lines below now carry the figure's actual trend claim.
    per_n_mean = df.groupby("N")[["gaussian_freq_cpu_s", "classifier_cpu_s"]].mean()

    xx = np.linspace(N_arr.min() * 0.9, N_arr.max() * 1.1, 100)
    if scale == "log":
        # Log-log OLS -- the conventional way to report a power-law
        # exponent, appears as a straight line on this log-log panel.
        # Labels carry only the equation (no "Gaussian"/"Classifier"
        # prefix, no "log-log fit"/"linear-space fit" descriptor): color
        # alone ties each fit line to its own series' scatter/marker, and
        # keeping the label short is what lets the legend box be placed
        # somewhere it doesn't overlap data or the fit lines themselves.
        ax.plot(xx, np.exp(g_intercept) * xx ** g_slope, color=gaussian_trend_color,
                ls="--", lw=1.2, zorder=4,
                label=f"$t \\propto N^{{{g_slope:.2f}}}$ ($R^2$={g_r2:.2f})")
        ax.plot(xx, np.exp(c_intercept) * xx ** c_slope, color=COLORS["cost_classifier"],
                ls="--", lw=1.2, zorder=4,
                label=f"$t \\propto N^{{{c_slope:.2f}}}$ ($R^2$={c_r2:.2f})")
    else:
        # Linear-space nonlinear fit -- minimizes actual CPU-second
        # residuals, the curve that tracks absolute vertical distance on
        # this linear panel (see plot_gaussian_nbasis_scaling for why this
        # differs from the log-log fit); the log-log fit itself is not
        # drawn here, only computed (see the summary dict) -- fig_cputime
        # is meant to show the one fit judged on the axes it's drawn on.
        ax.plot(xx, g_a_lin * xx ** g_b_lin, color=gaussian_trend_color,
                ls="--", lw=1.2, zorder=4,
                label=f"$t \\propto N^{{{g_b_lin:.2f}}}$ ($R^2$={g_r2_lin:.2f})")
        ax.plot(xx, c_a_lin * xx ** c_b_lin, color=COLORS["cost_classifier"],
                ls="--", lw=1.2, zorder=4,
                label=f"$t \\propto N^{{{c_b_lin:.2f}}}$ ($R^2$={c_r2_lin:.2f})")
    # Proxy legend entry (no real data), only needed when
    # gaussian_color_by_nbasis=True: the gradient series can't show its
    # colormap in a legend swatch, so a mid-inferno-toned circle stands in
    # (the colorbar alone carries the actual n_basis mapping). When False,
    # the flat-black scatter above already carries this same label for
    # real, so no proxy is needed.
    if gaussian_color_by_nbasis:
        gaussian_handle = Line2D([0], [0], marker="o", linestyle="none",
                                  markersize=5, markeredgecolor="none",
                                  markerfacecolor=matplotlib.colormaps["inferno"](0.5),
                                  label="Freq=hpmodes (MP2/3-21G)")

    ax.set_xlabel("Number of atoms, $N$")
    ax.set_ylabel("CPU time (s)")
    n_by_N = df.groupby("N").size()
    ax.set_xticks(sorted(df["N"].unique()))
    ax.set_xlim(df["N"].min() - 0.6, df["N"].max() + 0.6)
    if scale == "linear":
        # Headroom above the tallest point (BrH3, N=4, 55.3s) so the legend
        # box sits in clear space rather than overlapping it. AsBr3/PBr3 (the
        # far larger MP2/6-311G outliers) are excluded from `df` by the
        # mp2_321g filter above.
        # Bottom is a small NEGATIVE offset, not 0: the classifier series
        # (~0.002-0.03s) sits so close to zero that markers centered right
        # on the y=0 axis line get half-swallowed by it and read as
        # invisible -- lowering the axis floor below zero lifts them
        # visually clear of the line (nothing is actually plotted <0).
        y_top = df["gaussian_freq_cpu_s"].max() * 1.28
        ax.set_ylim(-0.06 * y_top, y_top)
    else:
        ax.set_yscale("log")
        # Extra headroom (x20, vs. linear's x1.28) above the tallest point
        # (C10H16, N=26, 460s): Gaussian's low-N cluster (2-55s) already
        # fills most of the panel's upper decade, so the legend needs a
        # dedicated, otherwise-empty decade above ALL data (including the
        # N=26 point) to sit in without overlapping any point -- verified
        # clear by rendering.
        ax.set_ylim(df["classifier_cpu_s"].min() * 0.5,
                    df["gaussian_freq_cpu_s"].max() * 20)
    # "upper left": low-N points/fit-line values stay well under the axis
    # ceiling across both scale variants (linear: tallest low-N point is
    # 55.3s against a >580s ceiling; log: tallest is the same 55.3s against
    # a >800s ceiling), so this corner stays clear even with the C10H16/
    # N=26 point and its steep fit line now in range -- unlike "upper
    # right", which the steep fit line grows into as it approaches N=26.
    # Verified clear by rendering; frameless is fine now that the shortened
    # equation-only labels (see above) keep the box small.
    handles, labels = ax.get_legend_handles_labels()
    if gaussian_color_by_nbasis:
        handles = [gaussian_handle] + handles
        labels = [gaussian_handle.get_label()] + labels
    ax.legend(handles=handles, labels=labels, loc="upper left", frameon=False,
              handletextpad=0.4, labelspacing=0.35, borderaxespad=0.3,
              fontsize=LEGEND_FONTSIZE)

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    gaussian_min_idx = df["gaussian_freq_cpu_s"].idxmin()
    gaussian_max_idx = df["gaussian_freq_cpu_s"].idxmax()
    classifier_min_idx = df["classifier_cpu_s"].idxmin()
    classifier_max_idx = df["classifier_cpu_s"].idxmax()

    summary = {
        "pdf": pdf_path, "png": png_path,
        "scale": scale,
        "gaussian_color_by_nbasis": bool(gaussian_color_by_nbasis),
        "n_molecules": len(df),
        "n_excluded_non_mp2_321g": int(n_excluded),
        "N_range": (int(df["N"].min()), int(df["N"].max())),
        "molecules_per_N": n_by_N.to_dict(),
        "gaussian_cpu_s_range": (float(df["gaussian_freq_cpu_s"].min()),
                                  float(df["gaussian_freq_cpu_s"].max())),
        "gaussian_cpu_s_min_molecule": str(df.loc[gaussian_min_idx, "molecule"]),
        "gaussian_cpu_s_max_molecule": str(df.loc[gaussian_max_idx, "molecule"]),
        "gaussian_cpu_s_median": float(df["gaussian_freq_cpu_s"].median()),
        "classifier_cpu_s_range": (float(df["classifier_cpu_s"].min()),
                                    float(df["classifier_cpu_s"].max())),
        "classifier_cpu_s_min_molecule": str(df.loc[classifier_min_idx, "molecule"]),
        "classifier_cpu_s_max_molecule": str(df.loc[classifier_max_idx, "molecule"]),
        "classifier_cpu_s_median": float(df["classifier_cpu_s"].median()),
        "ratio_gaussian_over_classifier": {
            "min": float(ratio.min()), "min_molecule": str(df.loc[ratio.idxmin(), "molecule"]),
            "median": float(ratio.median()),
            "max": float(ratio.max()), "max_molecule": str(df.loc[ratio.idxmax(), "molecule"]),
        },
        "classifier_mean_cpu_s_by_N": per_n_mean["classifier_cpu_s"].to_dict(),
        "gaussian_mean_cpu_s_by_N": per_n_mean["gaussian_freq_cpu_s"].to_dict(),
        "fit_vs_N_loglog_gaussian": {"exponent": float(g_slope), "r2": float(g_r2)},
        "fit_vs_N_loglog_classifier": {"exponent": float(c_slope), "r2": float(c_r2)},
        "fit_vs_N_linear_space_gaussian": {"exponent": float(g_b_lin), "r2_linear": float(g_r2_lin)},
        "fit_vs_N_linear_space_classifier": {"exponent": float(c_b_lin), "r2_linear": float(c_r2_lin)},
        "framing": ("empirical replacement/companion for the theoretical "
                    "Big-O 'Computational cost' section (tab:cost); PROPOSED "
                    "label fig:cputime, not yet wired into the .tex. Filtered "
                    "to mp2_321g==True, excluding "
                    f"{n_excluded} molecules run at a different Gaussian "
                    "method/basis, to avoid confounding CPU-time-vs-N with "
                    "CPU-time-vs-method. Per-N mean lines replaced with a "
                    "power-law fit per series, fit on raw per-molecule "
                    "points, not per-N means -- log-log OLS drawn on the "
                    "log scale variant, linear-space nonlinear drawn on the "
                    "linear scale variant; both fits are computed and "
                    "returned here regardless of which is drawn."),
    }
    return summary


def plot_gaussian_nbasis_scaling(
    benchmark_csv="data/results/cpu_time_benchmark.csv",
    out_dir="data/figures",
    label=None,
    scale="log",
):
    """SI/diagnostic companion to fig:cputime: Gaussian's freq-only CPU time
    vs. N_basis (AO basis-function count, Gaussian's own ``NBasis=``).

    Restricted to the MP2/3-21G subset (``mp2_321g == True``, 50 of 68
    molecules) -- matching fig:cputime's own filtering, and for the same
    reason: different methods/basis choices (MP2 vs. B3LYP vs. HF,
    correlation treatment, ECPs) have different cost *prefactors and
    exponents* even at matched N_basis, so a single power-law fit across
    mixed methods is confounded by method choice on top of N_basis: the
    excluded "other method/basis" molecules sit systematically ABOVE a fit
    line built from all 68, not scattered around it -- genuinely off-trend
    rather than just noisier. The 18 excluded molecules are dropped
    entirely (not plotted at all, not just excluded from the fit) -- they're
    off a different cost curve, so showing them alongside the MP2/3-21G
    trend doesn't add information, only clutter.

    ``scale`` picks both axes together: "log" (default, label
    "fig_gaussian_nbasis") shows a single fit -- log-log OLS (``np.polyfit``
    on ln(t) vs. ln(N_basis)), which appears as a straight line here, the
    standard way to report a power-law exponent. "linear" (label
    "fig_gaussian_nbasis_linear") shows TWO fits: that same log-log-OLS fit
    (transformed back to linear space) plus a second fit computed by
    nonlinear least squares directly in linear space (``scipy.curve_fit``,
    minimizing actual CPU-second residuals, not log residuals). These
    disagree (exponent ~2.0 vs. ~3.2) because log-log OLS minimizes
    *relative* error uniformly across ~2 orders of magnitude of CPU time, so
    it is not drawn toward matching the largest-N_basis points in absolute
    (second) terms the way a human eye judging a linear plot expects; the
    linear-space fit is the one that actually minimizes visual vertical
    distance on that panel. Both are legitimate, answering different
    questions ("what's the exponent, weighting all points by relative
    error" vs. "what curve best tracks absolute CPU-time on this axis") --
    shown together so the reader sees the disagreement rather than one
    fit presented as if it were the only answer. Both scale variants are
    generated per standing instruction. ``label`` overrides the
    scale-based default if given.

    Motivation: fig:cputime's "CPU time vs. N" comparison is fair to the
    classifier (whose cost genuinely depends on N -- geometry/mode-vector
    work) but not to Gaussian, whose SCF/MP2 cost depends on basis-function
    count, not atom count -- a single heavy/ECP-bearing atom can carry as
    many basis functions as several light atoms. Restricted to MP2/3-21G,
    fitting log(CPU) vs. log(N) for Gaussian gives R^2 ~ 0.04 (pure noise);
    refitting against log(n_basis) gives R^2 ~ 0.81 -- confirming n_basis,
    not N, is Gaussian's real scaling variable within a fixed method/basis.
    NOT extended to classifier_cpu_s: the classifier never touches AO basis
    functions, so n_basis is not a meaningful covariate for it (a category
    error) -- this figure is Gaussian-only by design.
    """
    if scale not in ("linear", "log"):
        raise ValueError(f"scale must be 'linear' or 'log', got {scale!r}")
    if label is None:
        label = "fig_gaussian_nbasis" if scale == "log" else "fig_gaussian_nbasis_linear"

    _style()
    df_all = pd.read_csv(benchmark_csv).copy()
    df = df_all[df_all["mp2_321g"] == True].copy()  # noqa: E712 -- explicit bool filter
    n_excluded = len(df_all) - len(df)

    n_basis_arr = df["n_basis"].to_numpy(float)
    t_arr = df["gaussian_freq_cpu_s"].to_numpy(float)
    log_n, log_t = np.log(n_basis_arr), np.log(t_arr)

    # Fit 1: log-log OLS -- minimizes relative/log-space error, the
    # conventional way to report a power-law exponent.
    slope, intercept = np.polyfit(log_n, log_t, 1)
    pred = slope * log_n + intercept
    r2 = 1 - np.sum((log_t - pred) ** 2) / np.sum((log_t - log_t.mean()) ** 2)

    # Fit 2: nonlinear least squares directly in linear (CPU-second) space --
    # minimizes the actual vertical distance a reader judges on a linear
    # plot. p0 seeded from Fit 1 so curve_fit starts near the right basin.
    def _powerlaw(n, a, b):
        return a * n ** b
    (a_lin, b_lin), _ = curve_fit(_powerlaw, n_basis_arr, t_arr,
                                   p0=[np.exp(intercept), slope])
    pred_lin = _powerlaw(n_basis_arr, a_lin, b_lin)
    r2_lin = 1 - np.sum((t_arr - pred_lin) ** 2) / np.sum((t_arr - t_arr.mean()) ** 2)

    fig, ax = plt.subplots(figsize=(4.6, 3.9))

    ax.scatter(df["n_basis"], df["gaussian_freq_cpu_s"], marker="s", s=24,
               facecolors="black", edgecolors="black",
               alpha=1.0, zorder=3, label="Freq=hpmodes (MP2/3-21G)")

    xx = np.linspace(df["n_basis"].min() * 0.9, df["n_basis"].max() * 1.1, 100)
    ax.plot(xx, np.exp(intercept) * xx ** slope, ls="--", lw=1.2,
            color=COLORS["threshold"], zorder=4,
            label=f"log-log fit: $t \\propto N_{{\\mathrm{{basis}}}}^{{{slope:.2f}}}$ ($R^2$={r2:.2f})")
    if scale == "linear":
        # Distinct dash pattern (dotted, not dashed) + plain black (not a
        # COLORS entry -- this is a one-off diagnostic overlay, not a
        # reusable classification-category encoding) so the two fits stay
        # visually distinguishable.
        ax.plot(xx, _powerlaw(xx, a_lin, b_lin), ls=":", lw=1.6,
                color="black", zorder=4,
                label=f"linear-space fit: $t \\propto N_{{\\mathrm{{basis}}}}^{{{b_lin:.2f}}}$ ($R^2$={r2_lin:.2f})")

    if scale == "log":
        ax.set_xscale("log")
        ax.set_yscale("log")
    else:
        ax.set_ylim(0, df["gaussian_freq_cpu_s"].max() * 1.15)
        ax.set_xlim(0, df["n_basis"].max() * 1.08)
    ax.set_xlabel(r"$N_{\mathrm{basis}}$")
    ax.set_ylabel("CPU time (s)")
    ax.legend(loc="upper left", frameon=False, handletextpad=0.4,
              labelspacing=0.35, borderaxespad=0.3, fontsize=LEGEND_FONTSIZE)

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    log_N = np.log(df["N"].to_numpy(float))
    slope_N, intercept_N = np.polyfit(log_N, log_t, 1)
    pred_N = slope_N * log_N + intercept_N
    r2_N = 1 - np.sum((log_t - pred_N) ** 2) / np.sum((log_t - log_t.mean()) ** 2)

    return {
        "pdf": pdf_path, "png": png_path,
        "scale": scale,
        "n_molecules_fit": len(df),
        "n_excluded_non_mp2_321g": int(n_excluded),
        "n_basis_range_fit": (int(df["n_basis"].min()), int(df["n_basis"].max())),
        "fit_vs_n_basis_loglog": {"exponent": float(slope), "r2": float(r2)},
        "fit_vs_n_basis_linear_space": {"exponent": float(b_lin), "r2_linear": float(r2_lin)},
        "fit_vs_N_for_comparison": {"exponent": float(slope_N), "r2": float(r2_N)},
        "framing": ("SI/diagnostic companion to fig:cputime, PROPOSED label "
                    "fig:gaussian_nbasis, not yet wired into the .tex. Fit "
                    "restricted to mp2_321g==True -- mixing methods "
                    "confounds the N_basis fit the same way mixing methods "
                    "confounded the N fit in fig:cputime. Shows n_basis "
                    "(not N) is Gaussian's real cost-scaling variable "
                    "within a fixed method/basis."),
    }


def plot_cpu_time_scaling_comparison(
    benchmark_csv="data/results/cpu_time_benchmark.csv",
    out_dir="data/figures",
    label="fig_cputime_scaling_comparison",
):
    """Dedicated log-log companion to fig:cputime: Gaussian frequency-calc
    CPU time and classification-algorithm CPU time, both plotted directly
    against atom count N on the same log-log panel, each with its own
    log-log OLS power-law fit -- built to make the two series' N-scaling
    exponents directly comparable at a glance, which fig:cputime itself
    isn't built for (it carries a continuous n_basis colorbar, error bars,
    and x-jitter for its own reasons -- see plot_cpu_time_benchmark -- that
    compete for attention with the fit comparison this figure exists for).

    Restricted to the MP2/3-21G subset (``mp2_321g == True``) -- identical
    filter and rationale to plot_cpu_time_benchmark and
    plot_gaussian_nbasis_scaling: mixing methods/bases confounds a CPU-time-
    vs-N fit with CPU-time-vs-method, so the excluded molecules are dropped
    entirely rather than diluting the trend.

    Log-log-only by design (no linear-scale variant, unlike the other two
    cpu-time figures): this figure's whole purpose is exponent comparison,
    which is a log-log-native question -- a linear-space fit answers "what
    curve tracks absolute CPU time," a question already covered by
    plot_cpu_time_benchmark's linear-scale variant. No jitter, colorbar, or
    error bars either: those exist elsewhere to explain vertical scatter
    within a given N (n_basis for Gaussian) or measurement uncertainty
    (classifier stddev), neither of which this figure is trying to convey --
    it deliberately trades that detail for legibility of the two exponents.

    Fits use raw per-molecule points, not per-N means, for consistency with
    plot_cpu_time_benchmark's and plot_gaussian_nbasis_scaling's own fits.
    """
    _style()
    df_all = pd.read_csv(benchmark_csv).copy()
    df = df_all[df_all["mp2_321g"] == True].copy()  # noqa: E712 -- explicit bool filter
    n_excluded = len(df_all) - len(df)

    N_arr = df["N"].to_numpy(float)
    gaussian_arr = df["gaussian_freq_cpu_s"].to_numpy(float)
    classifier_arr = df["classifier_cpu_s"].to_numpy(float)
    g_slope, g_intercept, g_r2 = _loglog_ols_fit(N_arr, gaussian_arr)
    c_slope, c_intercept, c_r2 = _loglog_ols_fit(N_arr, classifier_arr)

    fig, ax = plt.subplots(figsize=(4.6, 3.9))

    # Same color family as plot_cpu_time_benchmark's own Gaussian fit line
    # and classifier series, so a reader moving between the two figures
    # doesn't have to re-learn which color is which.
    gaussian_color = matplotlib.colormaps["inferno"](0.25)
    ax.scatter(N_arr, gaussian_arr, marker="o", s=24,
               facecolors=gaussian_color, edgecolors="none",
               alpha=1.0, zorder=3, label="Freq=hpmodes (MP2/3-21G)")
    ax.scatter(N_arr, classifier_arr, marker="s", s=24,
               facecolors=COLORS["cost_classifier"], edgecolors="none",
               alpha=1.0, zorder=3, label="Classification algorithm")

    xx = np.linspace(N_arr.min() * 0.9, N_arr.max() * 1.1, 100)
    ax.plot(xx, np.exp(g_intercept) * xx ** g_slope, color=gaussian_color,
            ls="--", lw=1.2, zorder=4,
            label=f"Gaussian log-log fit: $t \\propto N^{{{g_slope:.2f}}}$ ($R^2$={g_r2:.2f})")
    ax.plot(xx, np.exp(c_intercept) * xx ** c_slope, color=COLORS["cost_classifier"],
            ls="--", lw=1.2, zorder=4,
            label=f"Classifier log-log fit: $t \\propto N^{{{c_slope:.2f}}}$ ($R^2$={c_r2:.2f})")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of atoms, $N$")
    ax.set_ylabel("CPU time (s)")
    # Bordered legend (see plot_cpu_time_benchmark's identical choice): the
    # Gaussian scatter is dense enough at low N that a frameless legend
    # text ends up sitting directly over several points.
    ax.legend(loc="upper left", frameon=True, edgecolor="black",
              facecolor="white", framealpha=1.0, handletextpad=0.4,
              labelspacing=0.35, borderaxespad=0.3, fontsize=LEGEND_FONTSIZE)

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    return {
        "pdf": pdf_path, "png": png_path,
        "n_molecules": len(df),
        "n_excluded_non_mp2_321g": int(n_excluded),
        "N_range": (int(df["N"].min()), int(df["N"].max())),
        "fit_vs_N_loglog_gaussian": {"exponent": float(g_slope), "r2": float(g_r2)},
        "fit_vs_N_loglog_classifier": {"exponent": float(c_slope), "r2": float(c_r2)},
        "framing": ("dedicated log-log companion to fig:cputime, isolating "
                    "each series' own N-scaling exponent without the "
                    "colorbar/error-bar/jitter apparatus needed there; "
                    "PROPOSED label fig:cputime_scaling_comparison, not yet "
                    "referenced in the .tex. Filtered to mp2_321g==True for "
                    "the same reason as fig:cputime/fig:gaussian_nbasis."),
    }


def regenerate_all(verbose=True):
    """Regenerate every manuscript figure in one call. Each figure function
    reads its own already-computed ``data/results/*.csv`` inputs with their
    own defaults; this function takes no molecule-specific arguments.

    Returns a dict {tex_label: result_dict} for all 13 figures, in the same
    order they are built. Raises whatever the underlying plot_* function
    raises (e.g. a missing input CSV) -- fail loud, no silent partial
    regeneration.
    """
    fns = [
        ("fig:benzene", plot_benzene_stress_test),
        ("benzene-normal-modes gallery (no fig: label yet)", plot_benzene_normal_modes),
        ("fig:confusion (binary scheme, CANONICAL)", plot_confusion_matrix_binary),
        ("fig:confusion threeway scheme (fig_confusion_threeway)", plot_confusion_matrix),
        ("SI confusion retention/migration (no fig: label yet)", plot_confusion_retention_migration),
        ("SI rigorous-tier consistency check (no fig: label yet)", plot_rigorous_tier_check),
        ("fig:transferabilityconfusion (proposed; binary scheme, CANONICAL)",
         plot_transferability_confusion_binary),
        ("fig:transferabilityconfusion threeway scheme (fig_transferability_confusion_threeway)",
         plot_transferability_confusion),
        ("fig:benzeneconfusion (binary scheme, CANONICAL)", plot_benzene_internal_confusion_binary),
        ("fig:benzeneconfusion threeway scheme (fig_benzene_confusion_threeway)",
         plot_benzene_internal_confusion),
        ("SI benzene precision/recall (no fig: label yet)", plot_benzene_confusion_precision_recall),
        ("fig:benzeneemitcounts", plot_benzene_emit_counts),
        ("fig:bondscores", plot_bond_scores),
        ("fig:boxplots", plot_boxplots),
        ("fig:modemixing", plot_mode_mixing),
        ("SI irrep-degeneracy coupling (no fig: label yet)", plot_irrep_coupling),
        ("fig:sensitivity", plot_sensitivity),
        ("fig:sensitivity binary scheme (no fig: label yet)", plot_tau_sb_sensitivity),
        ("fig:cputime linear (proposed, not yet in .tex)",
         lambda: plot_cpu_time_benchmark(scale="linear")),
        ("fig:cputime log (proposed, not yet in .tex)",
         lambda: plot_cpu_time_benchmark(scale="log")),
        ("fig:cputime linear, black Gaussian markers (alt., not yet in .tex)",
         lambda: plot_cpu_time_benchmark(scale="linear", gaussian_color_by_nbasis=False)),
        ("fig:gaussian_nbasis log (SI/diagnostic, not yet in .tex)",
         lambda: plot_gaussian_nbasis_scaling(scale="log")),
        ("fig:gaussian_nbasis linear (SI/diagnostic, not yet in .tex)",
         lambda: plot_gaussian_nbasis_scaling(scale="linear")),
        ("fig:cputime_scaling_comparison (SI/diagnostic, not yet in .tex)",
         plot_cpu_time_scaling_comparison),
        ("fig:ped_vs_vscore (proposed, not yet in .tex)", plot_ped_vs_vscore),
        ("fig:ped_vs_bondscore (exploratory, not a manuscript figure)",
         plot_ped_vs_bondscore_by_type),
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


# --------------------------------------------------------------------------
# fig:ped_vs_vscore / fig:ped_vs_bondscore -- VEDA4 PED-based %nu vs. this
# framework's own s[V_S] / s^AB scores, from ped/merge_ped_scores.py's
# data/results/combined_ped_vs_scores.csv.
# --------------------------------------------------------------------------

# combined_ped_vs_scores.csv's own `label` column stores classify_all_modes'
# short S/B/SB codes directly (see ped/merge_ped_scores.py) -- map to the
# shared CATEGORY_COLOR/CATEGORY_MARKER/CATEGORY_LABEL keys (classifier-
# PREDICTED colors, not REF_CATEGORY_COLOR, since this label is this
# framework's own output, not a literature/reference label).
_LABEL_CODE_TO_CATEGORY = {"S": "stretch", "B": "bend", "SB": "mixed"}

# Gap between the two members of a genuinely degenerate pair/triple is 0 or
# a ~1e-4 cm^-1 numerical-diagonalization artifact; the next-smallest gap
# between distinct (non-degenerate) modes in this dataset is ~0.37 cm^-1 --
# a clean two-orders-of-magnitude separation, so 0.01 cm^-1 safely clusters
# degenerate partners without merging genuinely different modes.
_DEGENERATE_FREQ_TOL_CM1 = 0.01


def _collapse_degenerate_freqs(df, value_cols, freq_col="Freq", label_col=None,
                                tol=_DEGENERATE_FREQ_TOL_CM1):
    """Average ``value_cols`` over sets of rows sharing one physical
    frequency within a molecule (``freq_col`` values equal to within
    ``tol`` cm^-1). VEDA reports PED for each individual member of a
    degenerate mode independently, but an individual member's PED is not
    physically meaningful on its own -- it depends on the arbitrary linear
    combination diagonalization happened to pick within the degenerate
    subspace; only the set-averaged value is. Returns one row per
    degenerate group (Molecule + clustered Freq), with non-degenerate modes
    (group size 1) passed through unchanged. If ``label_col`` is given, the
    most common label in the group is carried through (ties broken by
    first occurrence) for scatter coloring.
    """
    df = df.sort_values(["Molecule", freq_col]).reset_index(drop=True)
    freqs = df[freq_col].to_numpy(float)
    mols = df["Molecule"].to_numpy()
    gid = np.zeros(len(df), dtype=int)
    g = -1
    for i in range(len(df)):
        if i == 0 or mols[i] != mols[i - 1] or freqs[i] - freqs[i - 1] > tol:
            g += 1
        gid[i] = g
    df = df.assign(_gid=gid)

    agg = {c: "mean" for c in value_cols}
    if label_col is not None:
        agg[label_col] = lambda s: s.mode().iloc[0]
    out = df.groupby(["_gid", "Molecule"], as_index=False, sort=False).agg(agg)
    return out.drop(columns=["_gid"])


def _linear_fit_r2(x, y, degree=1):
    """Ordinary least-squares polynomial fit (default degree 1, a straight
    line y ~ slope*x + intercept); return (coeffs, r2), matching the
    np.polyfit + manual-R^2 pattern used throughout this module (see
    plot_gaussian_nbasis_scaling).

    Returned in np.polyfit/np.polyval coefficient order (highest power
    first), so callers plot it with np.polyval exactly as before.
    """
    coeffs = np.polyfit(x, y, degree)
    pred = np.polyval(coeffs, x)
    r2 = 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)
    return coeffs, r2


def plot_ped_vs_vscore(
    csv_input="data/results/combined_ped_vs_scores.csv",
    mol_list_csv="data/mol_list_method.csv",
    out_dir="data/figures",
    label="fig_ped_vs_vscore",
    tau_SB=None,
    fit_degree=1,
    ylabel=r"$s[\mathrm{V_S}]$",
    label_source_csv=None,
):
    """VEDA4's PED-based %nu (``PED_Stretch_pct``) vs. this framework's own
    molecule-level stretch score s[V_S] (``V_Stretch``), one point per
    physical frequency (degenerate modes averaged together via
    _collapse_degenerate_freqs -- see its docstring) of
    combined_ped_vs_scores.csv, fit with a straight line (ordinary least
    squares, reported with R^2) summarizing the overall trend. Points colored/
    colored by this framework's own S/B/SB classify_all_modes label (the
    csv's own `label` column), reusing CATEGORY_COLOR/CATEGORY_LABEL as-is
    via _LABEL_CODE_TO_CATEGORY. One marker shape (circle) for all points --
    color alone distinguishes the category, so varying marker shape too
    would be a redundant second encoding of the same distinction (matches
    fig:bondscores' marker="o" override).

    Scoped to mol_type=='test' molecules only (the 9-molecule held-out
    transferability set): combined_ped_vs_scores.csv also happens to carry
    C6H6 (a PED cross-check run for a different purpose), which must NOT
    appear here -- joined in via mol_list_csv on
    combined_ped_vs_scores.csv's `Molecule` column (note the capitalization
    mismatch vs. the roster's lowercase `molecule`).

    ``tau_SB`` defaults to ``Thresholds.calibrated().tau_SB`` (None ->
    canonical, zero behavior change). When given, every row's ``label`` is
    re-derived from its own (pre-degenerate-averaging) ``V_Stretch`` via
    ``vib_label_binary`` -- this CSV carries internal (stretch/bend) rows
    only (VEDA's PED comparison has no external T/R analogue), so the plain
    binary label is used directly rather than the external-aware
    ``rescheme_internal_label``. Re-derived before
    ``_collapse_degenerate_freqs`` so the degenerate-group majority-label
    vote sees the overridden per-row labels, matching how the canonical path
    already collapses labels baked in at tau_SB=0.50.

    ``fit_degree`` selects the OLS polynomial order for the summary curve
    (default 1, matching the canonical figure exactly -- zero behavior
    change). 2 fits a quadratic instead; the legend label switches from
    "linear fit" to "quadratic fit" accordingly.

    ``ylabel`` overrides the y-axis label (default the canonical s[V_S]
    math label; e.g. archive_unweighted callers pass a label that says
    "Unweighted" to distinguish this from the canonical mu-weighted figure
    of the same name).

    ``label_source_csv``, if given, points at a DIFFERENT combined_ped_vs_
    scores.csv (matched to this one's rows by Molecule+Mode) whose own
    V_Stretch column is used to derive the S/B coloring via
    ``vib_label_binary(v, tau_SB)`` -- ``tau_SB`` still applies here (falls
    back to ``Thresholds.calibrated().tau_SB`` the same as the no-override
    path). The plotted x/y (``PED_Stretch_pct``/``V_Stretch``) and the
    fitted curve always come from ``csv_input`` -- this only swaps which
    definition's classification colors the points, e.g. coloring an
    unweighted-definition scatter by the canonical mu-weighted S/B split.
    Raises ValueError if any row fails to find a Molecule+Mode match.
    """
    _style()
    df = pd.read_csv(csv_input)
    roster = pd.read_csv(mol_list_csv)
    test_molecules = set(roster.loc[roster["mol_type"] == "test", "molecule"])
    df = df[df["Molecule"].isin(test_molecules)]
    df = df.dropna(subset=["PED_Stretch_pct", "V_Stretch"])
    n_raw_modes = len(df)
    if label_source_csv is not None:
        df = df.copy()
        src_tau_SB = tau_SB if tau_SB is not None else Thresholds.calibrated().tau_SB
        src = pd.read_csv(label_source_csv).dropna(subset=["V_Stretch"])
        src = src.assign(label=src["V_Stretch"].map(lambda v: vib_label_binary(v, src_tau_SB)))
        df = df.drop(columns=["label"]).merge(
            src[["Molecule", "Mode", "label"]], on=["Molecule", "Mode"], how="left")
        if df["label"].isna().any():
            missing = df.loc[df["label"].isna(), ["Molecule", "Mode"]]
            raise ValueError(
                f"label_source_csv={label_source_csv!r} has no Molecule+Mode "
                f"match for {len(missing)} row(s), e.g. {missing.iloc[0].to_dict()}")
    elif tau_SB is not None:
        df = df.copy()
        df["label"] = df["V_Stretch"].map(lambda v: vib_label_binary(v, tau_SB))
    df = _collapse_degenerate_freqs(
        df, value_cols=["PED_Stretch_pct", "V_Stretch"], label_col="label")

    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    for code, cat in _LABEL_CODE_TO_CATEGORY.items():
        sub = df[df["label"] == code]
        if sub.empty:
            continue
        kw = _marker_kwargs(cat, marker="o")
        ax.scatter(sub["PED_Stretch_pct"], sub["V_Stretch"], s=16,
                   zorder=3, label=CATEGORY_LABEL[cat], **kw)

    x = df["PED_Stretch_pct"].to_numpy(float)
    y = df["V_Stretch"].to_numpy(float)
    coeffs, r2 = _linear_fit_r2(x, y, degree=fit_degree)
    fit_name = {1: "linear", 2: "quadratic"}.get(fit_degree, f"degree-{fit_degree}")
    xx = np.linspace(x.min(), x.max(), 200)
    ax.plot(xx, np.polyval(coeffs, xx), ls="--", lw=1.2,
            color=COLORS["threshold"], zorder=4,
            label=f"{fit_name} fit ($R^2$={r2:.2f})")

    ax.set_xlabel(r"$\%\nu$")
    ax.set_ylabel(ylabel)
    ax.set_xlim(-3, 103)
    ax.set_ylim(-0.03, 1.05)
    ax.legend(loc="upper left", frameon=False, handletextpad=0.4,
              labelspacing=0.35, borderaxespad=0.3, fontsize=LEGEND_FONTSIZE)

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    return {
        "pdf": pdf_path, "png": png_path,
        "n_frequency_points": len(df),
        "n_raw_modes_before_degenerate_averaging": n_raw_modes,
        "n_molecules": df["Molecule"].nunique(),
        "fit_degree": fit_degree,
        "fit_coeffs": tuple(float(c) for c in coeffs),
        "r2": float(r2),
        "tau_SB": tau_SB if tau_SB is not None else Thresholds.calibrated().tau_SB,
    }


def plot_ped_vs_bondscore_by_type(
    csv_input="data/results/combined_ped_vs_scores.csv",
    mol_list_csv="data/mol_list_method.csv",
    out_dir="data/figures",
    label="fig_ped_vs_bondscore",
    ylabel=r"$|s^{AB}|$",
    label_source_csv=None,
    tau_SB=None,
    show_fit=True,
):
    """Exploratory per-bond-type grid: VEDA4's ``PED_S_<T>_pct`` (%nu^AB) vs.
    this framework's own raw ``BondScore_<T>`` (s^AB), one panel per bond
    species-pair type T present in combined_ped_vs_scores.csv's own
    BondScore_<T> columns (currently B-H, B-N, C-C, C-Cl, C-H, C-N, C-O,
    Cl-P, H-N -- see ped/merge_ped_scores.py's ``_canonical_bond_type``).
    Each panel gets its own straight-line OLS fit + R^2 (same pattern as
    plot_ped_vs_vscore), for spotting which bond types/molecules diverge
    from the molecule-level trend. Exploratory, not a manuscript figure:
    by default one plain marker color, no S/B/SB faceting (see
    ``label_source_csv`` below to turn that on). One point per physical
    frequency -- degenerate modes are averaged together first (see
    _collapse_degenerate_freqs), same as plot_ped_vs_vscore.

    Scoped to mol_type=='test' molecules only -- see plot_ped_vs_vscore's
    docstring for why (excludes C6H6, joined in via mol_list_csv on
    combined_ped_vs_scores.csv's `Molecule` column).

    ``ylabel`` overrides the shared y-axis label (default the canonical
    |s^AB| math label).

    ``label_source_csv``, if given, points at a DIFFERENT combined_ped_vs_
    scores.csv (matched to this one's rows by Molecule+Mode) whose own
    V_Stretch column is used to derive a per-MODE S/B split via
    ``vib_label_binary(v, tau_SB)`` (``tau_SB`` falls back to
    ``Thresholds.calibrated().tau_SB`` same as plot_ped_vs_vscore) --
    every bond of a given mode is colored by that mode's own S/B, same
    convention as plot_ped_vs_vscore's coloring. Default None keeps the
    original plain-black scheme (zero behavior change).

    ``show_fit`` (default True, zero behavior change) draws each panel's
    OLS fit line + R^2 annotation ("insufficient data" for n<3 panels).
    False skips the fit entirely -- just the raw scatter, every panel,
    regardless of point count.
    """
    _style()
    df = pd.read_csv(csv_input)
    roster = pd.read_csv(mol_list_csv)
    test_molecules = set(roster.loc[roster["mol_type"] == "test", "molecule"])
    df = df[df["Molecule"].isin(test_molecules)]

    color_by_label = label_source_csv is not None
    if color_by_label:
        df = df.copy()
        src_tau_SB = tau_SB if tau_SB is not None else Thresholds.calibrated().tau_SB
        src = pd.read_csv(label_source_csv).dropna(subset=["V_Stretch"])
        src = src.assign(label=src["V_Stretch"].map(lambda v: vib_label_binary(v, src_tau_SB)))
        df = df.drop(columns=["label"], errors="ignore").merge(
            src[["Molecule", "Mode", "label"]], on=["Molecule", "Mode"], how="left")
        if df["label"].isna().any():
            missing = df.loc[df["label"].isna(), ["Molecule", "Mode"]]
            raise ValueError(
                f"label_source_csv={label_source_csv!r} has no Molecule+Mode "
                f"match for {len(missing)} row(s), e.g. {missing.iloc[0].to_dict()}")

    bond_types = [c[len("BondScore_"):] for c in df.columns
                  if c.startswith("BondScore_") and not c.endswith("_pct")]
    value_cols = [f"PED_S_{bt}_pct" for bt in bond_types] + \
                 [f"BondScore_{bt}" for bt in bond_types]
    df = _collapse_degenerate_freqs(
        df, value_cols=value_cols, label_col="label" if color_by_label else None)

    n = len(bond_types)
    ncols = 3
    nrows = -(-n // ncols)  # ceil division
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.8 * nrows),
                              squeeze=False)

    legend_handles = {}
    per_bond_type = {}
    for i, bt in enumerate(bond_types):
        ax = axes.flat[i]
        x_col, y_col = f"PED_S_{bt}_pct", f"BondScore_{bt}"
        cols = [x_col, y_col] + (["label"] if color_by_label else [])
        sub = df[cols].dropna(subset=[x_col, y_col])
        x = sub[x_col].to_numpy(float)
        y = sub[y_col].to_numpy(float)

        if color_by_label:
            for code, cat in _LABEL_CODE_TO_CATEGORY.items():
                mask = sub["label"] == code
                if not mask.any():
                    continue
                kw = _marker_kwargs(cat, marker="o")
                h = ax.scatter(sub.loc[mask, x_col], sub.loc[mask, y_col], s=14,
                               zorder=3, label=CATEGORY_LABEL[cat], **kw)
                legend_handles.setdefault(cat, h)
        else:
            ax.scatter(x, y, s=14, marker="o", facecolors="black",
                       edgecolors="black", alpha=0.75, zorder=3)

        if show_fit:
            # >=3, not >=2: two points define a line exactly, so R^2 would
            # be a meaningless 1.0.
            if len(sub) >= 3:
                coeffs, r2 = _linear_fit_r2(x, y)
                xx = np.linspace(x.min(), x.max(), 100)
                ax.plot(xx, np.polyval(coeffs, xx), ls="--", lw=1.1,
                        color=COLORS["threshold"], zorder=4)
                ax.text(0.05, 0.92, f"$R^2$={r2:.2f}", transform=ax.transAxes,
                        fontsize=ANNOTATION_FONTSIZE, va="top")
                per_bond_type[bt] = {"n": len(sub), "r2": float(r2)}
            else:
                ax.text(0.5, 0.5, "insufficient data", transform=ax.transAxes,
                        ha="center", va="center", fontsize=ANNOTATION_FONTSIZE,
                        color=COLORS["threshold"])
                per_bond_type[bt] = {"n": len(sub), "r2": None}
        else:
            per_bond_type[bt] = {"n": len(sub), "r2": None}

        ax.set_title(bt.replace("-", "–"), fontsize=11)

    for j in range(n, nrows * ncols):
        axes.flat[j].axis("off")

    fig.supxlabel(r"$\%\nu^{AB}$")
    fig.supylabel(ylabel)
    if color_by_label and legend_handles:
        fig.legend(legend_handles.values(), [h.get_label() for h in legend_handles.values()],
                   loc="upper right", frameon=False, fontsize=LEGEND_FONTSIZE)
    fig.tight_layout(rect=(0.02, 0.02, 1, 1))
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    return {
        "pdf": pdf_path, "png": png_path,
        "bond_types": bond_types,
        "per_bond_type": per_bond_type,
    }


if __name__ == "__main__":
    regenerate_all()
