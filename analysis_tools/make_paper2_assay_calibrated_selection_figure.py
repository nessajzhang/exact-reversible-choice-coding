#!/usr/bin/env python3
"""Build the Paper 2 assay-calibrated selection figure and source data."""

from __future__ import annotations

import os
from pathlib import Path

# Stabilize Matplotlib PDF/SVG metadata for byte-reproducible figure exports.
os.environ.setdefault("SOURCE_DATE_EPOCH", "1784332800")

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import patches

from make_submission_figures import (
    PALETTE,
    add_panel_label,
    save_figure,
    setup_style,
    style_axis,
)


ROOT = Path(__file__).resolve().parents[1]
PAPER_DIR = ROOT / "papers" / "paper2_thermodynamic_risk_coding"
SELECTION_DIR = PAPER_DIR / "bioinformatics_reframe" / "assay_calibrated_selection"
PUBLIC_DIR = PAPER_DIR / "bioinformatics_reframe" / "public_experimental_validation"
FIGURE_BASE = PAPER_DIR / "figures" / "assay_calibrated_selection"
SOURCE_DATA = ROOT / "data" / "paper2_assay_calibrated_selection_figure_source_data.tsv"


MODEL_ORDER = [
    "F0_frozen_low_proxy",
    "P0_composition",
    "P1_variable_structure",
    "P2_assay_context",
    "P3_assay_calibrated",
    "P4_full_structural_context",
    "P5_combined_context",
]
MODEL_LABELS = [
    "Frozen\nlow-W",
    "Composition",
    "Variable\nstructure",
    "Adapter\ncontext",
    "+ variable\nW",
    "+ full-oligo\nstructure",
    "+ full-oligo\nW",
]

SERIES = [
    ("repeated_nested_oof", "GCall", "GCall OOF", PALETTE["blue"], "o"),
    ("repeated_nested_oof", "GCfix", "GCfix OOF", PALETTE["teal"], "o"),
    (
        "source_only_transfer",
        "GCall_to_GCfix",
        "GCall→GCfix",
        PALETTE["gold"],
        "s",
    ),
    (
        "source_only_transfer",
        "GCfix_to_GCall",
        "GCfix→GCall",
        PALETTE["salmon"],
        "s",
    ),
]


def add_flow_box(
    ax,
    x: float,
    y: float,
    width: float,
    height: float,
    text: str,
    facecolor: str,
    edgecolor: str,
) -> None:
    box = patches.FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        linewidth=0.8,
        facecolor=facecolor,
        edgecolor=edgecolor,
        transform=ax.transAxes,
    )
    ax.add_patch(box)
    ax.text(
        x + width / 2,
        y + height / 2,
        text,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=5.7,
        color=PALETTE["neutral_0"],
        linespacing=1.18,
    )


def add_arrow(ax, start: tuple[float, float], end: tuple[float, float], dashed: bool = False) -> None:
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        xycoords=ax.transAxes,
        textcoords=ax.transAxes,
        arrowprops={
            "arrowstyle": "-|>",
            "color": PALETTE["neutral_2"],
            "lw": 0.8,
            "linestyle": (0, (2, 2)) if dashed else "solid",
            "shrinkA": 2,
            "shrinkB": 2,
        },
    )


def schematic_panel(ax) -> None:
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    add_panel_label(ax, "a", x=-0.035, y=1.03)
    ax.text(
        0.0,
        0.93,
        "Retrospective source-only\ncross-pool transfer",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.8,
        fontweight="bold",
        color=PALETTE["neutral_0"],
    )

    ax.text(0.0, 0.73, "direction 1", transform=ax.transAxes, fontsize=5.7, fontweight="bold", color=PALETTE["blue"])
    first = [
        (0.01, "GCall\n11,998"),
        (0.32, "train FullContext\nsource only"),
        (0.67, "test GCfix\n11,994"),
    ]
    for x, label in first:
        add_flow_box(ax, x, 0.52, 0.265, 0.14, label, PALETTE["blue_soft"], PALETTE["blue_mid"])
    for left, right in zip(first[:-1], first[1:]):
        add_arrow(ax, (left[0] + 0.265, 0.59), (right[0], 0.59))

    ax.text(0.0, 0.39, "direction 2", transform=ax.transAxes, fontsize=5.7, fontweight="bold", color=PALETTE["teal"])
    second = [
        (0.01, "GCfix\n11,994"),
        (0.32, "train FullContext\nsource only"),
        (0.67, "test GCall\n11,998"),
    ]
    for x, label in second:
        add_flow_box(ax, x, 0.17, 0.265, 0.14, label, PALETTE["teal_soft"], PALETTE["teal"])
    for left, right in zip(second[:-1], second[1:]):
        add_arrow(ax, (left[0] + 0.265, 0.24), (right[0], 0.24))

def performance_panel(ax, prediction: pd.DataFrame) -> pd.DataFrame:
    add_panel_label(ax, "b", x=-0.08, y=1.10)
    style_axis(ax)
    # F0 is an unfitted proxy, not the first member of the nested P0--P5
    # ridge ladder.  Give it a separate grey marker and a visible gap so the
    # figure does not imply a fitted-model trajectory from F0 to P0.
    x = np.asarray([0.0, 1.35, 2.35, 3.35, 4.35, 5.35, 6.35])
    used: list[pd.DataFrame] = []
    for evaluation, pool, label, color, marker in SERIES:
        subset = prediction.loc[
            prediction["evaluation"].eq(evaluation)
            & prediction["pool_or_direction"].eq(pool)
            & prediction["model"].isin(MODEL_ORDER)
        ].set_index("model")
        values = np.asarray([subset.loc[model, "spearman"] for model in MODEL_ORDER])
        ax.plot(
            x[1:],
            values[1:],
            color=color,
            marker=marker,
            markersize=3.7,
            linewidth=1.2,
            markeredgecolor="white",
            markeredgewidth=0.45,
            label=label,
        )
        ax.plot(
            x[0],
            values[0],
            linestyle="none",
            marker="X",
            markersize=4.0,
            markerfacecolor=PALETTE["neutral_2"],
            markeredgecolor="white",
            markeredgewidth=0.45,
            zorder=3,
        )
        used.append(subset.reset_index())
    ax.axhline(0, color=PALETTE["neutral_3"], linewidth=0.75)
    ax.axvline(0.68, color=PALETTE["neutral_3"], linewidth=0.7, linestyle=(0, (2, 2)))
    ax.text(
        0.0,
        -0.035,
        "unfitted\nbaseline",
        ha="center",
        va="center",
        fontsize=5.2,
        color=PALETTE["neutral_1"],
    )
    ax.set_xticks(x, MODEL_LABELS)
    ax.set_xlim(-0.35, 6.7)
    ax.set_ylabel("Spearman with relative PCR efficiency")
    ax.set_ylim(-0.17, 0.39)
    ax.legend(loc="upper left", ncols=2, handlelength=1.6, columnspacing=0.9)
    ax.set_title("Assay context improves within-pool and transferred ranking", loc="left")
    return pd.concat(used, ignore_index=True)


def retention_panel(ax, selection: pd.DataFrame) -> pd.DataFrame:
    add_panel_label(ax, "c", x=-0.10, y=1.10)
    style_axis(ax)
    used = selection.loc[
        selection["dataset"].eq("Gimpel2025_PCR")
        & selection["model"].eq("P5_combined_context")
    ].copy()
    for evaluation, pool, label, color, marker in SERIES:
        subset = used.loc[
            used["evaluation"].eq(evaluation)
            & used["pool_or_direction"].eq(pool)
        ].sort_values("retention_fraction")
        ax.plot(
            subset["retention_fraction"] * 100,
            subset["standardized_gain_vs_pool"],
            color=color,
            marker=marker,
            markersize=3.8,
            linewidth=1.2,
            markeredgecolor="white",
            markeredgewidth=0.45,
        )
    ax.axhline(0, color=PALETTE["neutral_3"], linewidth=0.75)
    ax.set_xticks([10, 25, 50], ["10\n0.0308", "25\n0.0185", "50\n0.0093"])
    ax.set_xlabel("Selected pool (%)\nfinite-pool penalty (bits/variable-region nt)")
    ax.set_ylabel("Selected-mean gain (outcome s.d.)")
    ax.set_ylim(0, 0.42)
    ax.set_title("Selection utility persists across fixed retention", loc="left")
    ax.text(
        0.98,
        0.04,
        "25%: 7/240 and 9/240\nlow-efficiency labels retained",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=5.5,
        color=PALETTE["neutral_1"],
    )
    return used


def increment_panel(ax, comparisons: pd.DataFrame) -> pd.DataFrame:
    add_panel_label(ax, "d", x=-0.10, y=1.10)
    style_axis(ax)
    wanted = ["full_structural_context_increment", "full_weighted_increment"]
    used = comparisons.loc[
        comparisons["dataset"].eq("Gimpel2025_PCR")
        & comparisons["retention_fraction"].eq(0.25)
        & comparisons["comparison"].isin(wanted)
    ].copy()
    labels = [item[2] for item in SERIES]
    y = np.arange(len(labels))[::-1]
    offsets = {
        "full_structural_context_increment": 0.10,
        "full_weighted_increment": -0.10,
    }
    style = {
        "full_structural_context_increment": (PALETTE["teal"], "o", "+ full-oligo structure"),
        "full_weighted_increment": (PALETTE["gold"], "s", "+ full-oligo weighted score"),
    }
    for comparison in wanted:
        color, marker, label = style[comparison]
        values = []
        lows = []
        highs = []
        for evaluation, pool, *_rest in SERIES:
            row = used.loc[
                used["evaluation"].eq(evaluation)
                & used["pool_or_direction"].eq(pool)
                & used["comparison"].eq(comparison)
            ].iloc[0]
            values.append(float(row["delta_selected_mean"]) * 1e4)
            lows.append(float(row["delta_selected_mean_ci_2p5"]) * 1e4)
            highs.append(float(row["delta_selected_mean_ci_97p5"]) * 1e4)
        values = np.asarray(values)
        errors = np.vstack([values - np.asarray(lows), np.asarray(highs) - values])
        ax.errorbar(
            values,
            y + offsets[comparison],
            xerr=errors,
            fmt=marker,
            color=color,
            ecolor=color,
            markersize=4.0,
            capsize=1.7,
            linewidth=1.0,
            markeredgecolor="white",
            markeredgewidth=0.45,
            label=label,
        )
    ax.axvline(0, color=PALETTE["neutral_3"], linewidth=0.75)
    ax.set_yticks(y, labels)
    ax.set_xlabel(r"$\Delta$ selected mean efficiency ($\times 10^{-4}$)")
    ax.set_title("Full-oligo features add matched utility", loc="left")
    ax.set_ylim(-0.5, 4.2)
    ax.legend(loc="upper left", fontsize=5.4, handlelength=1.2)
    ax.xaxis.grid(True, color=PALETTE["neutral_4"], linewidth=0.46, alpha=0.54)
    ax.yaxis.grid(False)
    return used


def dt_boundary_panel(ax, public_comparisons: pd.DataFrame) -> pd.DataFrame:
    add_panel_label(ax, "e", x=-0.10, y=1.10)
    style_axis(ax)
    used = public_comparisons.loc[
        public_comparisons["analysis"].eq("dt_continuous")
        & public_comparisons["comparison"].eq("D3_minus_D1")
    ].copy()
    order = [
        ("Genscript_GCall", "mean_log2_cpm_plus1", "Genscript mean"),
        ("Twist_GCall", "mean_log2_cpm_plus1", "Twist mean"),
        ("Genscript_GCall", "day7_vs_day0_log2fc_plus1", "Genscript day 7"),
        ("Twist_GCall", "day7_vs_day0_log2fc_plus1", "Twist day 7"),
    ]
    y = np.arange(len(order))[::-1]
    for yi, (pool, endpoint, _label) in zip(y, order):
        row = used.loc[used["pool_or_direction"].eq(pool) & used["endpoint"].eq(endpoint)].iloc[0]
        value = float(row["delta_extended_minus_baseline"])
        low = float(row["bootstrap_2p5"])
        high = float(row["bootstrap_97p5"])
        ax.errorbar(
            value,
            yi,
            xerr=[[value - low], [high - value]],
            fmt="o",
            color=PALETTE["neutral_1"],
            ecolor=PALETTE["neutral_2"],
            markersize=3.8,
            capsize=1.7,
            linewidth=1.0,
            markeredgecolor="white",
            markeredgewidth=0.45,
        )
    ax.axvline(0, color=PALETTE["red"], linewidth=0.8, linestyle=(0, (3, 2)))
    ax.set_yticks(y, [item[2] for item in order])
    ax.set_xlabel("Δ OOF Spearman: weighted + structural\nminus structural rules")
    ax.set_title("DT4DDS adds no weighted-score increment", loc="left")
    ax.xaxis.grid(True, color=PALETTE["neutral_4"], linewidth=0.46, alpha=0.54)
    ax.yaxis.grid(False)
    return used


def source_rows(panel: str, frame: pd.DataFrame) -> pd.DataFrame:
    output = frame.copy()
    output.insert(0, "panel", panel)
    return output


def main() -> None:
    setup_style()
    mpl.rcParams["svg.hashsalt"] = "paper2-20260718"
    prediction = pd.read_csv(SELECTION_DIR / "prediction_metrics.tsv", sep="\t")
    selection = pd.read_csv(SELECTION_DIR / "selection_metrics.tsv", sep="\t")
    comparisons = pd.read_csv(SELECTION_DIR / "selection_model_comparisons.tsv", sep="\t")
    public_comparisons = pd.read_csv(PUBLIC_DIR / "model_comparisons.tsv", sep="\t")

    fig = plt.figure(figsize=(183 / 25.4, 151 / 25.4))
    grid = fig.add_gridspec(
        2,
        3,
        width_ratios=[1.05, 1.12, 1.05],
        height_ratios=[0.95, 1.05],
        left=0.055,
        right=0.985,
        top=0.96,
        bottom=0.105,
        hspace=0.46,
        wspace=0.45,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1:])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])
    ax_e = fig.add_subplot(grid[1, 2])

    schematic_panel(ax_a)
    b_data = performance_panel(ax_b, prediction)
    c_data = retention_panel(ax_c, selection)
    d_data = increment_panel(ax_d, comparisons)
    e_data = dt_boundary_panel(ax_e, public_comparisons)

    FIGURE_BASE.parent.mkdir(parents=True, exist_ok=True)
    save_figure(fig, FIGURE_BASE.with_suffix(".pdf"))
    plt.close(fig)

    source = pd.concat(
        [
            source_rows("b", b_data),
            source_rows("c", c_data),
            source_rows("d", d_data),
            source_rows("e", e_data),
        ],
        ignore_index=True,
        sort=False,
    )
    SOURCE_DATA.parent.mkdir(parents=True, exist_ok=True)
    source.to_csv(SOURCE_DATA, sep="\t", index=False, float_format="%.12g")
    print(f"Wrote {FIGURE_BASE}.svg/.pdf/.png/.tiff")
    print(f"Wrote {SOURCE_DATA}")


if __name__ == "__main__":
    main()
