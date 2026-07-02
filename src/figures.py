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

All 6 manuscript figures are implemented: ``plot_benzene_stress_test``
(``fig:benzene``), ``plot_confusion_matrix`` (``fig:confusion``),
``plot_bond_scores`` (``fig:bondscores``), ``plot_boxplots``
(``fig:boxplots``), ``plot_mode_mixing`` (``fig:modemixing``), and
``plot_sensitivity`` (``fig:sensitivity``).

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

from src.classifier import (
    Thresholds, vib_label,
    STRETCHING, BENDING, MIXED_STRETCH_BEND,
    CLEAN_TRANSLATION, CLEAN_ROTATION, MIXED_EXTERNAL_WITH_VIBRATION,
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
    "highlight_r": "#E69F00",# orange    -- s[R] inversion (EMIT 2 / 9)
    "highlight_t": "#0072B2",# blue      -- flagged external (EMIT 34-36)
    "threshold": "#555555",  # dark gray -- tau reference lines (all figures)
    # -- New, figure-specific-but-shared-key colors (fig:sensitivity / fig:confusion) --
    "sens_accuracy": "#56B4E9",   # sky blue  -- accuracy curve, fig:sensitivity
    "sens_change": "#000000",     # black     -- label-change-fraction curve, fig:sensitivity
    "plateau_band": "#56B4E9",    # sky blue @ low alpha -- tau_TR plateau shading, fig:sensitivity
    "confusion_cmap": "Blues",    # sequential, colorblind-safe -- fig:confusion heatmap
}

# Reference-label (library ground truth) -> shared classification-category
# constant, so fig:confusion/fig:bondscores/fig:boxplots/fig:modemixing all
# color/mark "stretch"/"bend"/"translation"/"rotation" identically to
# fig:benzene's STRETCHING/BENDING/CLEAN_TRANSLATION/CLEAN_ROTATION points.
REF_LABEL_TO_CATEGORY = {
    "translation": CLEAN_TRANSLATION,
    "rotation": CLEAN_ROTATION,
    "stretch": STRETCHING,
    "bend": BENDING,
}
# Predicted-bucket (confusion_matrix_stats' "_pred_bucket") -> category,
# extending the above with the two bucket names that only appear as
# predictions, never as reference labels.
PRED_BUCKET_TO_CATEGORY = dict(REF_LABEL_TO_CATEGORY)
PRED_BUCKET_TO_CATEGORY.update({
    "mixed": MIXED_STRETCH_BEND,
    "mixed_external": MIXED_EXTERNAL_WITH_VIBRATION,
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

CATEGORY_COLOR = {
    "CLEAN_TRANSLATION": COLORS["external"],
    "CLEAN_ROTATION": COLORS["external"],
    "BENDING": COLORS["bending"],
    "STRETCHING": COLORS["stretching"],
    "MIXED_STRETCH_BEND": COLORS["mixed"],
    "MIXED_EXTERNAL_WITH_VIBRATION": COLORS["mixed_ext"],
}

CATEGORY_MARKER = {
    "CLEAN_TRANSLATION": "X",
    "CLEAN_ROTATION": "X",
    "BENDING": "o",
    "STRETCHING": "s",
    "MIXED_STRETCH_BEND": "^",
    "MIXED_EXTERNAL_WITH_VIBRATION": "P",
}

CATEGORY_LABEL = {
    "CLEAN_TRANSLATION": "clean translation",
    "CLEAN_ROTATION": "clean rotation",
    "BENDING": "bending",
    "STRETCHING": "stretching",
    "MIXED_STRETCH_BEND": "mixed stretch/bend",
    "MIXED_EXTERNAL_WITH_VIBRATION": "mixed external+vibration",
}

# Provisional classifier thresholds (src/classifier.py Thresholds defaults;
# quoted here only for the reference dashed lines in fig:benzene panel (a) --
# not refit, not invented).
TAU_S = 0.9
TAU_B = 0.2


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


def _marker_kwargs(category, ideal_flag=None):
    """Shared per-point marker styling for a classification `category`
    (key into CATEGORY_COLOR/CATEGORY_MARKER), optionally faceted by
    `ideal_flag` ('yes'/'no'/None) via the shared IDEAL_STYLE encoding
    (filled = ideal, hollow = non-ideal). Returns a dict ready to splat into
    ax.scatter(...).
    """
    color = CATEGORY_COLOR[category]
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
    normal_csv="data/results/benzene_normal_classified.csv",
    emit_classified_csv="data/results/benzene_EMIT_classified.csv",
    emit_contrib_csv="data/results/benzene_EMIT_contributions.csv",
    out_dir="data/figures",
    label="fig_benzene",
):
    """Build fig:benzene: (a) s[V_S] vs. frequency for the 36 normal modes;
    (b) score vs. projected normal-mode contribution for the 36 EMIT modes,
    highlighting the EMIT 2/9 s[R] inversion and the EMIT 34-36 flagged
    externals.

    This is a flag-behavior / non-monotonicity illustration (per
    IMPLEMENTATION_PLAN.md Phase 2), NOT a T/R-accuracy parity plot: panel
    (b) deliberately omits any 1:1 reference line and plots |score| against
    projected contribution fraction only to expose where the two disagree.

    Returns a summary dict with output paths and a few sanity numbers.
    """
    _style()

    normal = pd.read_csv(normal_csv)
    emit = pd.read_csv(emit_classified_csv)
    contrib = pd.read_csv(emit_contrib_csv)
    emit = emit.merge(contrib, on="Mode", suffixes=("", "_c"))

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.0, 3.2))

    # ---------------- Panel (a): s[V_S] vs frequency, 36 normal modes -----
    seen_labels = set()
    for _, row in normal.iterrows():
        cat = row["label"]
        color = CATEGORY_COLOR.get(cat, "black")
        marker = CATEGORY_MARKER.get(cat, "o")
        plot_kwargs = dict(color=color, marker=marker, s=26,
                            edgecolors="black", linewidths=0.3, zorder=3)
        leg_label = CATEGORY_LABEL.get(cat, cat) if cat not in seen_labels else None
        seen_labels.add(cat)
        ax_a.scatter(row["Freq"], row["V_Stretch"], label=leg_label, **plot_kwargs)

    ax_a.axhline(TAU_S, color=COLORS["threshold"], ls="--", lw=0.8, zorder=1)
    ax_a.axhline(TAU_B, color=COLORS["threshold"], ls="--", lw=0.8, zorder=1)
    ax_a.text(0.98, TAU_S, r"$\tau_{\mathrm{stretch}}\approx0.9$", ha="right",
               va="bottom", fontsize=7, color=COLORS["threshold"],
               transform=ax_a.get_yaxis_transform())
    ax_a.text(0.98, TAU_B, r"$\tau_{\mathrm{bend}}\approx0.2$", ha="right",
               va="top", fontsize=7, color=COLORS["threshold"],
               transform=ax_a.get_yaxis_transform())

    ax_a.set_xlabel(r"Frequency (cm$^{-1}$)")
    ax_a.set_ylabel(r"$s[\mathrm{V_S}]$")
    ax_a.set_ylim(-0.05, 1.08)
    ax_a.set_xlim(-120, normal["Freq"].max() * 1.05)
    ax_a.set_title("(a) 36 normal modes", loc="left", fontweight="bold", fontsize=9)
    ax_a.legend(loc="center right", frameon=False, handletextpad=0.3,
                labelspacing=0.3, borderaxespad=0.1)

    # ---------------- Panel (b): score vs projected contribution, EMIT ----
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
    ax_b.annotate(
        "", xy=(e9["C2_Ry"], abs(e9["Ry"])), xytext=(e2["C2_Ry"], abs(e2["Ry"])),
        arrowprops=dict(arrowstyle="->", color=COLORS["highlight_r"], lw=1.1),
        zorder=4,
    )
    ax_b.annotate("EMIT 2", (e2["C2_Ry"], abs(e2["Ry"])),
                  xytext=(6, 6), textcoords="offset points", fontsize=7)
    ax_b.annotate("EMIT 9", (e9["C2_Ry"], abs(e9["Ry"])),
                  xytext=(6, -10), textcoords="offset points", fontsize=7)

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
    ax_b.set_title("(b) EMIT modes: score vs. projection", loc="left",
                    fontweight="bold", fontsize=9)
    ax_b.legend(loc="upper left", frameon=False, handletextpad=0.3,
                labelspacing=0.35, borderaxespad=0.1, fontsize=6.7)
    ax_b.text(0.98, 0.02,
              "non-monotonic by design\n(flag mechanism, not a parity check)",
              transform=ax_b.transAxes, ha="right", va="bottom",
              fontsize=6.3, style="italic", color="#333333")

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary = {
        "pdf": pdf_path,
        "png": png_path,
        "panel_a_n_points": len(normal),
        "panel_a_freq_range": (float(normal["Freq"].min()), float(normal["Freq"].max())),
        "panel_a_vs_range": (float(normal["V_Stretch"].min()), float(normal["V_Stretch"].max())),
        "panel_b_n_background_points": len(bg_x),
        "panel_b_emit2_Ry_score": float(abs(e2["Ry"])),
        "panel_b_emit2_Ry_contribution": float(e2["C2_Ry"]),
        "panel_b_emit9_Ry_score": float(abs(e9["Ry"])),
        "panel_b_emit9_Ry_contribution": float(e9["C2_Ry"]),
        "panel_b_emit34_36_scores": fy,
        "panel_b_emit34_36_contributions": fx,
        "panel_b_emit34_36_vstretch": fvs,
    }
    return summary


# --------------------------------------------------------------------------
# fig:confusion -- clean-category confusion matrix + precision/recall
# --------------------------------------------------------------------------

def plot_confusion_matrix(
    library_csv="data/results/library_scores.csv",
    out_dir="data/figures",
    label="fig_confusion",
):
    """Build fig:confusion: (a) reference-label x predicted-bucket confusion
    matrix over the whole ingested hydride library (both Excel-only and
    geometry-backed rows), (b) per-category precision/recall against the
    acceptance floor. Numbers come straight from
    ``src.calibrate.confusion_matrix_stats`` (calibrated thresholds) -- this
    function only renders them, never recomputes.

    Shared-category reuse: row/column tick labels and panel-(b) bars are
    colored with the SAME CATEGORY_COLOR used by fig:benzene (gray=clean
    T/R, blue=bending, vermillion=stretching, teal=mixed stretch/bend,
    purple=mixed external+vibration).
    """
    _style()
    from src.calibrate import confusion_matrix_stats

    lib_df = pd.read_csv(library_csv)
    thresholds = Thresholds.calibrated()
    stats = confusion_matrix_stats(lib_df, thresholds)
    table = stats["confusion_table"]

    ref_order = ["translation", "rotation", "stretch", "bend"]
    pred_order = ["translation", "rotation", "stretch", "bend", "mixed", "mixed_external"]
    pred_cols = [c for c in pred_order if c in table.columns] + \
                [c for c in table.columns if c not in pred_order]
    tbl = table.reindex(index=ref_order, columns=pred_cols, fill_value=0)

    fig, (ax_h, ax_p) = plt.subplots(1, 2, figsize=(7.2, 3.3),
                                      gridspec_kw={"width_ratios": [1.2, 1.0]})

    vmax = tbl.values.max()
    im = ax_h.imshow(tbl.values, cmap=COLORS["confusion_cmap"], aspect="auto",
                      vmin=0, vmax=vmax)
    for i in range(tbl.shape[0]):
        for j in range(tbl.shape[1]):
            v = int(tbl.values[i, j])
            txt_color = "white" if v > 0.6 * vmax else "black"
            ax_h.text(j, i, str(v), ha="center", va="center", fontsize=7.5,
                      color=txt_color)

    ax_h.set_xticks(range(len(tbl.columns)))
    ax_h.set_xticklabels([CATEGORY_LABEL[PRED_BUCKET_TO_CATEGORY[c]] if c in
                           PRED_BUCKET_TO_CATEGORY else c for c in tbl.columns],
                          rotation=30, ha="right")
    ax_h.set_yticks(range(len(tbl.index)))
    ax_h.set_yticklabels([CATEGORY_LABEL[REF_LABEL_TO_CATEGORY[r]] for r in tbl.index])
    for tick, r in zip(ax_h.get_yticklabels(), tbl.index):
        tick.set_color(CATEGORY_COLOR[REF_LABEL_TO_CATEGORY[r]])
    for tick, c in zip(ax_h.get_xticklabels(), tbl.columns):
        cat = PRED_BUCKET_TO_CATEGORY.get(c)
        if cat is not None:
            tick.set_color(CATEGORY_COLOR[cat])
    ax_h.set_xlabel("Predicted bucket")
    ax_h.set_ylabel("Reference label (symmetry / literature)")
    ax_h.set_title("(a) Confusion matrix", loc="left", fontweight="bold", fontsize=9)
    cb = fig.colorbar(im, ax=ax_h, fraction=0.046, pad=0.04)
    cb.set_label("n modes", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    cats = ["translation", "rotation", "stretch", "bend"]
    x = np.arange(len(cats))
    width = 0.35
    precisions = [stats["per_category"][c]["precision"] for c in cats]
    recalls = [stats["per_category"][c]["recall"] for c in cats]
    bar_colors = [CATEGORY_COLOR[REF_LABEL_TO_CATEGORY[c]] for c in cats]
    ax_p.bar(x - width / 2, precisions, width, color=bar_colors,
             edgecolor="black", linewidth=0.5, label="precision")
    ax_p.bar(x + width / 2, recalls, width, color=bar_colors,
             edgecolor="black", linewidth=0.5, hatch="///", label="recall")
    ax_p.axhline(stats["acceptance_floor"], color=COLORS["threshold"], ls="--", lw=0.8)
    ax_p.text(len(cats) - 0.5, stats["acceptance_floor"], "  floor",
              va="bottom", ha="right", fontsize=7, color=COLORS["threshold"])
    ax_p.set_xticks(x)
    ax_p.set_xticklabels([CATEGORY_LABEL[REF_LABEL_TO_CATEGORY[c]] for c in cats],
                         rotation=20, ha="right")
    ax_p.set_ylim(0, 1.12)
    ax_p.set_ylabel("Precision / recall")
    ax_p.set_title("(b) Per-category precision/recall", loc="left",
                    fontweight="bold", fontsize=9)
    ax_p.legend(loc="lower left", frameon=False, fontsize=7.5)

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary = {
        "pdf": pdf_path, "png": png_path,
        "shared_categories": ("reuses CATEGORY_COLOR/CATEGORY_LABEL for "
                               "CLEAN_TRANSLATION/CLEAN_ROTATION/STRETCHING/"
                               "BENDING/MIXED_STRETCH_BEND/"
                               "MIXED_EXTERNAL_WITH_VIBRATION -- same mapping "
                               "as fig:benzene."),
        "confusion_table": tbl.to_dict(),
        "precision": dict(zip(cats, precisions)),
        "recall": dict(zip(cats, recalls)),
        "acceptance_floor": stats["acceptance_floor"],
        "floor_met": stats["floor_met"],
    }
    return summary


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
    stretching (vermillion square) vs. bending (blue circle) -- reusing the
    shared CATEGORY_COLOR/CATEGORY_MARKER and IDEAL_STYLE encodings.
    """
    _style()
    lib_df = pd.read_csv(library_csv)
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
            kw = _marker_kwargs(cat, ideal_flag)
            ax.scatter(sub["abs_rel_db"], sub["s_AB"], s=14,
                       zorder=3 if ideal_flag == "yes" else 2, **kw)

    ax.set_xlabel(r"$|\Delta|\mathbf{b}|\,/\,|\mathbf{b}|\,|$ (relative bond-length change)")
    ax.set_ylabel(r"Bond score $s_{AB}$")
    ax.set_xlim(-0.03, bonds["abs_rel_db"].max() * 1.05)
    ax.set_ylim(-0.03, 1.05)

    legend_elems = [
        Line2D([0], [0], marker=CATEGORY_MARKER[STRETCHING], color="none",
               markerfacecolor=COLORS["stretching"], markeredgecolor=COLORS["stretching"],
               markersize=6, label="stretching, ideal"),
        Line2D([0], [0], marker=CATEGORY_MARKER[STRETCHING], color="none",
               markerfacecolor="none", markeredgecolor=COLORS["stretching"],
               markersize=6, label="stretching, non-ideal"),
        Line2D([0], [0], marker=CATEGORY_MARKER[BENDING], color="none",
               markerfacecolor=COLORS["bending"], markeredgecolor=COLORS["bending"],
               markersize=6, label="bending, ideal"),
        Line2D([0], [0], marker=CATEGORY_MARKER[BENDING], color="none",
               markerfacecolor="none", markeredgecolor=COLORS["bending"],
               markersize=6, label="bending, non-ideal"),
    ]
    ax.legend(handles=legend_elems, loc="upper left", frameon=False, fontsize=6.8,
              handletextpad=0.4, labelspacing=0.4, borderaxespad=0.2)

    fig.tight_layout()
    pdf_path, png_path = _savefig(fig, out_dir, label)
    plt.close(fig)

    summary = {
        "pdf": pdf_path, "png": png_path,
        "shared_categories": ("reuses CATEGORY_COLOR/CATEGORY_MARKER for "
                               "STRETCHING (vermillion square)/BENDING (blue "
                               "circle) and the shared IDEAL_STYLE filled="
                               "ideal/hollow=non-ideal encoding -- same "
                               "mapping as fig:benzene / fig:boxplots / "
                               "fig:modemixing."),
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
    """
    _style()
    lib_df = pd.read_csv(library_csv)
    internal = lib_df[(lib_df["kind"] == "internal") &
                       lib_df["ref_label"].isin(("stretch", "bend")) &
                       lib_df["ideal"].isin(("yes", "no"))].copy()

    groups = [("bend", "yes"), ("stretch", "yes"), ("bend", "no"), ("stretch", "no")]
    group_labels = ["bend\n(ideal)", "stretch\n(ideal)", "bend\n(non-ideal)", "stretch\n(non-ideal)"]

    panels = [
        ("freq", r"Frequency (cm$^{-1}$)", "(a) Frequency"),
        ("delta_b_mean", r"Averaged $|\Delta|\mathbf{b}|\,/\,|\mathbf{b}|\,|$", "(b) Bond-length change"),
        ("V_Stretch", r"$s[\mathrm{V_S}]$", "(c) Mode score"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(7.4, 3.2))
    for ax, (col, ylabel, title) in zip(axes, panels):
        data = []
        for ref, ideal_flag in groups:
            vals = internal.loc[(internal.ref_label == ref) & (internal.ideal == ideal_flag), col].dropna().values
            data.append(vals)
        bp = ax.boxplot(data, positions=range(1, 5), widths=0.6, showfliers=True,
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

        ax.set_xticks(range(1, 5))
        ax.set_xticklabels(group_labels, fontsize=7)
        ax.set_ylabel(ylabel)
        ax.set_title(title, loc="left", fontweight="bold", fontsize=9)

    # tau_S / tau_B reference lines on panel (c) only.
    th = Thresholds.calibrated()
    axes[2].axhline(th.tau_S, color=COLORS["threshold"], ls="--", lw=0.8)
    axes[2].axhline(th.tau_B, color=COLORS["threshold"], ls="--", lw=0.8)
    axes[2].text(4.55, th.tau_S, r"$\tau_S$", ha="left", va="center", fontsize=7,
                 color=COLORS["threshold"])
    axes[2].text(4.55, th.tau_B, r"$\tau_B$", ha="left", va="center", fontsize=7,
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
    the shared STRETCHING/BENDING category colors+markers; ideal (a) is all
    filled, non-ideal (b) is all hollow, per the shared IDEAL_STYLE.

    NOTE (pending gap, IMPLEMENTATION_PLAN.md / figure-builder standing
    report): this renders only the two-panel ideal-step-vs-non-ideal-
    gradient content the current .tex caption (fig:modemixing) describes.
    The irrep-degeneracy sub-panels (scoped-doc Fig 3c-d in the group's
    earlier report, keyed on trigonal-planar/bent-AB2 irreps there, prose
    mentions ethane here) are NOT built -- the intended molecule/panel form
    for THIS manuscript is still unconfirmed with lead-author/tex-data-sync,
    so nothing is invented for it.
    """
    _style()
    lib_df = pd.read_csv(library_csv)
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
            kw = _marker_kwargs(cat, ideal_flag)
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
        "shared_categories": ("reuses CATEGORY_COLOR/CATEGORY_MARKER for "
                               "STRETCHING/BENDING and the shared IDEAL_STYLE "
                               "filled=ideal/hollow=non-ideal encoding -- "
                               "same mapping as fig:benzene / fig:bondscores "
                               "/ fig:boxplots."),
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


if __name__ == "__main__":
    fns = [
        ("fig:benzene", plot_benzene_stress_test),
        ("fig:confusion", plot_confusion_matrix),
        ("fig:bondscores", plot_bond_scores),
        ("fig:boxplots", plot_boxplots),
        ("fig:modemixing", plot_mode_mixing),
        ("fig:sensitivity", plot_sensitivity),
    ]
    for tex_label, fn in fns:
        result = fn()
        print(f"{tex_label} ->", result["pdf"])
        print(f"{tex_label} ->", result["png"])
        for k, v in result.items():
            if k not in ("pdf", "png"):
                print(f"  {k}: {v}")
