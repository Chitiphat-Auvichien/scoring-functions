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

Only ``plot_benzene_stress_test`` (``fig:benzene``) is implemented so far;
the other figures are BLOCKED on upstream data (see
``IMPLEMENTATION_PLAN.md`` and the figure-builder agent's report) and are
stubbed below with a clear NotImplementedError pointing at the missing input.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

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
    "threshold": "#555555",
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
# Remaining JCC figures -- BLOCKED pending upstream data (see
# IMPLEMENTATION_PLAN.md and the figure-builder agent's standing report).
# Stubs kept here so reproduce.py has a stable import surface to grow into.
# --------------------------------------------------------------------------

def plot_confusion_matrix(*args, **kwargs):
    """fig:confusion -- BLOCKED: needs data/results/library_scores.csv
    (Phase 3 output of excel_ingest.py) with reference stretch/bend labels."""
    raise NotImplementedError(
        "fig:confusion is BLOCKED -- needs library_scores.csv (Phase 3, "
        "excel_ingest.py). Not yet produced."
    )


def plot_bond_scores(*args, **kwargs):
    """fig:bondscores -- BLOCKED: needs data/results/library_scores.csv."""
    raise NotImplementedError(
        "fig:bondscores is BLOCKED -- needs library_scores.csv (Phase 3, "
        "excel_ingest.py). Not yet produced."
    )


def plot_boxplots(*args, **kwargs):
    """fig:boxplots -- BLOCKED: needs data/results/library_scores.csv."""
    raise NotImplementedError(
        "fig:boxplots is BLOCKED -- needs library_scores.csv (Phase 3, "
        "excel_ingest.py). Not yet produced."
    )


def plot_mode_mixing(*args, **kwargs):
    """fig:modemixing -- BLOCKED: needs data/results/library_scores.csv, and
    the irrep-degeneracy panel content/molecule needs confirming with
    lead-author/tex-data-sync before building (see figure-builder agent's
    pending-gaps note)."""
    raise NotImplementedError(
        "fig:modemixing is BLOCKED -- needs library_scores.csv (Phase 3, "
        "excel_ingest.py), and the irrep-degeneracy panel spec is "
        "unconfirmed (ethane vs benzene). Not yet produced."
    )


def plot_sensitivity(*args, **kwargs):
    """fig:sensitivity -- BLOCKED: needs data/results/thresholds.json
    (Phase 3 output of calibrate.py)."""
    raise NotImplementedError(
        "fig:sensitivity is BLOCKED -- needs thresholds.json (Phase 3, "
        "calibrate.py). Not yet produced."
    )


if __name__ == "__main__":
    result = plot_benzene_stress_test()
    print("fig:benzene ->", result["pdf"])
    print("fig:benzene ->", result["png"])
    for k, v in result.items():
        if k not in ("pdf", "png"):
            print(f"  {k}: {v}")
