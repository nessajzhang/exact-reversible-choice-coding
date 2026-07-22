#!/usr/bin/env python3
"""Create the submission figure for the Paper 2 reversible choice codec."""

from __future__ import annotations

import csv
import os
from pathlib import Path
from typing import Any

# Stabilize Matplotlib PDF/SVG metadata for byte-reproducible figure exports.
os.environ.setdefault("SOURCE_DATE_EPOCH", "1784332800")

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parent.parent
PAPER_DIR = ROOT / "papers" / "paper2_thermodynamic_risk_coding"
DATA_DIR = PAPER_DIR / "bioinformatics_reframe" / "reversible_choice_codec"
SOTA_DIR = PAPER_DIR / "bioinformatics_reframe" / "sota_and_external_validation"
FIG_DIR = PAPER_DIR / "figures"
OUT_BASE = FIG_DIR / "reversible_choice_codec"
SOURCE_DATA = FIG_DIR / "reversible_choice_codec_source_data.tsv"
SUPP_OUT_BASE = FIG_DIR / "paired_two_stage_sensitivity"
SUPP_SOURCE_DATA = FIG_DIR / "paired_two_stage_sensitivity_source_data.tsv"

BLUE = "#2C6EAA"
BLUE_LIGHT = "#B8D5EA"
ORANGE = "#D9822B"
ORANGE_LIGHT = "#F2C79D"
PURPLE = "#6F5AA8"
TEAL = "#3C8D8D"
GREY = "#62666A"
LIGHT_GREY = "#E8EAEC"
DARK = "#24272A"
GREEN = "#31824A"
HARD_GREY = "#8A8D91"
Q5_RED = "#B24A4A"

METHOD_COLORS = {
    "P5_combined_context": BLUE,
    "P2_assay_context": TEAL,
    "published_1dcnn": ORANGE,
    "published_sota_hard_filter": HARD_GREY,
}
METHOD_MARKERS = {
    "P5_combined_context": "o",
    "P2_assay_context": "s",
    "published_1dcnn": "D",
    "published_sota_hard_filter": "^",
}
METHOD_LABELS = {
    "P5_combined_context": "FullContext ridge",
    "P2_assay_context": "AssayContext ridge",
    "published_1dcnn": "Released 1D-CNN score",
    "published_sota_hard_filter": "Published hard rule",
}


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
        "svg.fonttype": "none",
        "svg.hashsalt": "paper2-20260718",
        "pdf.fonttype": 42,
        "font.size": 6.5,
        "axes.labelsize": 6.5,
        "axes.titlesize": 7.2,
        "xtick.labelsize": 5.8,
        "ytick.labelsize": 5.8,
        "legend.fontsize": 5.6,
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "lines.linewidth": 1.35,
        "lines.markersize": 4.2,
    }
)


def read(name: str) -> pd.DataFrame:
    return pd.read_csv(DATA_DIR / name, sep="\t")


def read_sota(name: str) -> pd.DataFrame:
    return pd.read_csv(SOTA_DIR / name, sep="\t")


def panel_label(ax: plt.Axes, label: str, x: float = -0.10, y: float = 1.04) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=8.2,
        fontweight="bold",
        ha="left",
        va="bottom",
        color=DARK,
    )


def rounded_box(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    text: str,
    facecolor: str = "white",
    edgecolor: str = GREY,
    linewidth: float = 0.8,
    fontsize: float = 5.7,
    fontweight: str = "normal",
) -> FancyBboxPatch:
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
    )
    ax.add_patch(patch)
    ax.text(
        xy[0] + width / 2,
        xy[1] + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        fontweight=fontweight,
        color=DARK,
        linespacing=1.15,
    )
    return patch


def arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str = GREY,
    connectionstyle: str = "arc3",
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=7,
            linewidth=0.8,
            color=color,
            connectionstyle=connectionstyle,
        )
    )


def draw_method_panel(ax: plt.Axes) -> list[dict[str, Any]]:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    panel_label(ax, "a", x=-0.02, y=1.00)
    ax.set_title("Fixed candidate fibers separate assay selection from payload inversion", pad=2)

    rounded_box(ax, (0.015, 0.65), 0.14, 0.16, "payload m\n(213 − r bits)", BLUE_LIGHT, BLUE)
    rounded_box(ax, (0.205, 0.65), 0.17, 0.16, "logical ranks q\n(m << r) | j", "#F3EEF8", PURPLE)
    rounded_box(ax, (0.425, 0.65), 0.15, 0.16, "R(q) = aq + c mod N\ngcd(a,N) = 1", "#F3EEF8", PURPLE, fontsize=5.35)
    rounded_box(
        ax,
        (0.625, 0.65),
        0.18,
        0.16,
        "exact base codec E\n108 nt; GC 49–59; HP ≤ 3",
        "#E5F2EE",
        TEAL,
    )
    arrow(ax, (0.155, 0.73), (0.205, 0.73))
    arrow(ax, (0.375, 0.73), (0.425, 0.73))
    arrow(ax, (0.575, 0.73), (0.625, 0.73))

    candidates = [
        (0.805, 0.038, r"x$_0$", False),
        (0.852, 0.038, r"x$_1$", True),
        (0.899, 0.030, r"$\cdots$", False),
        (0.938, 0.057, r"x$_{2^r-1}$", False),
    ]
    for x, width, label, selected in candidates:
        rounded_box(
            ax,
            (x, 0.64),
            width,
            0.18,
            label,
            "#DCEEDC" if selected else "white",
            GREEN if selected else GREY,
            linewidth=1.1 if selected else 0.7,
            fontsize=5.0,
            fontweight="bold" if selected else "normal",
        )
    arrow(ax, (0.795, 0.73), (0.805, 0.73))
    ax.text(0.90, 0.86, "2^r exact candidates", ha="center", va="bottom", fontsize=5.6)
    ax.text(
        0.90,
        0.59,
        "choose smallest j maximizing fixed key κ_z(x_j)",
        ha="center",
        va="top",
        fontsize=5.8,
        color=GREEN,
    )

    rounded_box(ax, (0.79, 0.27), 0.13, 0.14, "selected x_j*", "#DCEEDC", GREEN, linewidth=1.0)
    rounded_box(ax, (0.63, 0.27), 0.12, 0.14, "rank D(x_j*)", "#E5F2EE", TEAL)
    rounded_box(ax, (0.47, 0.27), 0.12, 0.14, "R^-1 mod N", "#F3EEF8", PURPLE)
    rounded_box(ax, (0.31, 0.27), 0.12, 0.14, "require q < 2^213", "#F3EEF8", PURPLE, fontsize=5.25)
    rounded_box(ax, (0.15, 0.27), 0.12, 0.14, "right shift r bits", "#F3EEF8", PURPLE, fontsize=5.25)
    rounded_box(ax, (0.015, 0.27), 0.08, 0.14, "m", BLUE_LIGHT, BLUE, fontweight="bold")
    arrow(ax, (0.871, 0.64), (0.81, 0.41), GREEN)
    arrow(ax, (0.79, 0.34), (0.75, 0.34))
    arrow(ax, (0.63, 0.34), (0.59, 0.34))
    arrow(ax, (0.47, 0.34), (0.43, 0.34))
    arrow(ax, (0.31, 0.34), (0.27, 0.34))
    arrow(ax, (0.15, 0.34), (0.095, 0.34))

    rounded_box(
        ax,
        (0.80, 0.055),
        0.17,
        0.10,
        "optional canonical verifier\n(fixed score key + tie rule)",
        "white",
        GREEN,
        fontsize=5.15,
    )
    ax.add_patch(
        FancyArrowPatch(
            (0.855, 0.27),
            (0.875, 0.155),
            arrowstyle="-|>",
            mutation_scale=7,
            linewidth=0.75,
            linestyle="--",
            color=GREEN,
        )
    )

    ax.text(
        0.38,
        0.105,
        r"payload decoder omits the assay score    •    exact cost: r bits/oligo",
        ha="center",
        va="center",
        fontsize=5.9,
        fontweight="bold",
        color=DARK,
    )
    return [
        {
            "panel": "a",
            "series": "formal_audit",
            "category": "generated_candidates",
            "x": "",
            "y": 65024,
            "low": "",
            "high": "",
            "n": 65024,
            "unit": "candidate rank-unrank checks",
            "note": "zero failures; exact selector-bit cost",
        }
    ]


def draw_public_effect_panel(ax: plt.Axes, public: pd.DataFrame) -> list[dict[str, Any]]:
    data = public.loc[
        public["selector_model"].eq("P5_combined_context")
        & public["selector_bits"].isin([2, 4])
    ].copy()
    order = [
        ("GCall_to_GCfix", 2),
        ("GCall_to_GCfix", 4),
        ("GCfix_to_GCall", 2),
        ("GCfix_to_GCall", 4),
    ]
    labels = ["GCall→GCfix, r=2", "GCall→GCfix, r=4", "GCfix→GCall, r=2", "GCfix→GCall, r=4"]
    y = np.arange(len(order))[::-1]
    source_rows: list[dict[str, Any]] = []
    for yi, key, label in zip(y, order, labels):
        row = data.loc[
            data["direction"].eq(key[0]) & data["selector_bits"].eq(key[1])
        ].iloc[0]
        estimate = row["paired_mean_outcome_gain"] * 1000
        low = row["paired_gain_ci_2p5"] * 1000
        high = row["paired_gain_ci_97p5"] * 1000
        color = BLUE if key[0].startswith("GCall") else ORANGE
        ax.plot([low, high], [yi, yi], color=color, lw=1.8, solid_capstyle="round")
        ax.plot(estimate, yi, "o", color=color, mec="white", mew=0.5, ms=5.0, zorder=3)
        ax.text(high + 0.08, yi, f"n={int(row['payloads'])}", va="center", ha="left", fontsize=5.1, color=GREY)
        source_rows.append(
            {
                "panel": "b",
                "series": "FullContext measured gain",
                "category": label,
                "x": estimate,
                "y": yi,
                "low": low,
                "high": high,
                "n": int(row["payloads"]),
                "unit": "relative PCR efficiency x 1e-3",
                "note": "95% paired fiber-bootstrap interval; target outcomes held out from fitting and library construction",
            }
        )
    ax.axvline(0, color="#9A9A9A", lw=0.8, ls="--", zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Measured relative-PCR gain (x 10^-3)")
    ax.set_title("Held-out public codebooks select higher measured PCR outcomes")
    ax.set_xlim(-0.1, 4.25)
    ax.set_ylim(-0.6, 3.6)
    ax.grid(axis="x", color=LIGHT_GREY, lw=0.6)
    panel_label(ax, "b", x=-0.05, y=1.10)
    return source_rows


def draw_identical_fiber_benchmark(
    ax: plt.Axes,
    benchmark: pd.DataFrame,
    *,
    dataset: str,
    panel: str,
    title: str,
    methods: list[str],
    show_legend: bool,
) -> list[dict[str, Any]]:
    data = benchmark.loc[
        benchmark["dataset"].eq(dataset)
        & benchmark["selector_bits"].isin([2, 4])
        & benchmark["selector_model"].isin(methods)
    ].copy()
    directions = list(dict.fromkeys(data["direction"].tolist()))
    order = [(direction, selector_bits) for direction in directions for selector_bits in (2, 4)]
    row_labels = []
    for direction, selector_bits in order:
        source = direction.split("_to_")[0]
        target = direction.split("_to_")[1]
        target_label = "KAPA" if target == "external_Taq" else target
        row_labels.append(f"{source}→{target_label}, r={selector_bits}")

    if len(methods) == 4:
        offsets = [0.24, 0.08, -0.08, -0.24]
    else:
        offsets = [0.18, 0.0, -0.18]
    base_y = np.arange(len(order))[::-1]
    source_rows: list[dict[str, Any]] = []
    labels = dict(METHOD_LABELS)
    for method, offset in zip(methods, offsets):
        color = METHOD_COLORS[method]
        marker = METHOD_MARKERS[method]
        for yi, (direction, selector_bits), row_label in zip(base_y, order, row_labels):
            row = data.loc[
                data["direction"].eq(direction)
                & data["selector_bits"].eq(selector_bits)
                & data["selector_model"].eq(method)
            ].iloc[0]
            estimate = float(row["paired_mean_outcome_gain"]) * 1000
            low_value = row.get("paired_gain_bca_ci_2p5", np.nan)
            high_value = row.get("paired_gain_bca_ci_97p5", np.nan)
            if pd.isna(low_value) or pd.isna(high_value):
                low_value = row["paired_gain_ci_2p5"]
                high_value = row["paired_gain_ci_97p5"]
            low = float(low_value) * 1000
            high = float(high_value) * 1000
            yy = yi + offset
            ax.plot([low, high], [yy, yy], color=color, lw=1.25, solid_capstyle="round")
            ax.plot(
                estimate,
                yy,
                marker=marker,
                ms=4.0,
                color=color,
                mec="white" if method != "published_sota_hard_filter" else color,
                mew=0.45,
                linestyle="none",
                zorder=3,
            )
            source_rows.append(
                {
                    "panel": panel,
                    "series": labels[method],
                    "category": row_label,
                    "x": estimate,
                    "y": yy,
                    "low": low,
                    "high": high,
                    "n": int(row["payloads"]),
                    "unit": "relative PCR efficiency x 1e-3",
                    "note": "mean selected-minus-fiber-mean; 95% BCa target-fiber bootstrap interval",
                }
            )

    ax.axvline(0, color="#9A9A9A", lw=0.8, ls="--", zorder=0)
    ax.set_yticks(base_y)
    ax.set_yticklabels(row_labels)
    ax.set_xlabel(r"Measured selection gain ($\times 10^{-3}$)")
    ax.set_title(title, pad=4)
    ax.set_ylim(-0.55, len(order) - 0.45)
    ax.grid(axis="x", color=LIGHT_GREY, lw=0.6)
    ax.margins(x=0.08)
    panel_label(ax, panel, x=-0.08, y=1.08)
    if show_legend:
        handles = [
            Line2D(
                [0],
                [0],
                color=METHOD_COLORS[method],
                marker=METHOD_MARKERS[method],
                lw=1.2,
                ms=3.8,
                label=labels[method],
            )
            for method in methods
        ]
        ax.legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.27),
            ncol=2,
            handlelength=1.4,
            columnspacing=1.0,
            fontsize=5.1,
        )
    return source_rows


def draw_two_stage_panel(
    ax: plt.Axes,
    benchmark: pd.DataFrame,
    two_stage: pd.DataFrame,
) -> list[dict[str, Any]]:
    ts = two_stage.loc[
        two_stage["estimand"].eq("P5_combined_context_mean_gain")
        & two_stage["selector_bits"].isin([2, 4])
    ].copy()
    order = [("GCall", 2), ("GCall", 4), ("GCfix", 2), ("GCfix", 4)]
    base_y = np.arange(len(order))[::-1]
    source_rows: list[dict[str, Any]] = []
    datasets = [
        (
            "Gimpel2025_cross_pool_public_codebook",
            0.18,
            "o",
            True,
            BLUE,
            "Cross-pool",
            "matched assay",
        ),
        (
            "Gimpel2025_external_laboratory_Taq",
            0.0,
            "D",
            False,
            BLUE,
            "External KAPA",
            "matched assay; same source publication",
        ),
        (
            "Gimpel2025_external_laboratory_Q5_sensitivity",
            -0.18,
            "v",
            True,
            Q5_RED,
            "Q5 altered protocol",
            "secondary exploratory altered-protocol sensitivity",
        ),
    ]
    for dataset, offset, marker, filled, color, label, evidence_note in datasets:
        for yi, (source_pool, selector_bits) in zip(base_y, order):
            row = ts.loc[
                ts["target_dataset"].eq(dataset)
                & ts["source_pool"].eq(source_pool)
                & ts["selector_bits"].eq(selector_bits)
            ].iloc[0]
            fixed = benchmark.loc[
                benchmark["dataset"].eq(dataset)
                & benchmark["direction"].str.startswith(f"{source_pool}_to_")
                & benchmark["selector_bits"].eq(selector_bits)
                & benchmark["selector_model"].eq("P5_combined_context")
            ].iloc[0]
            estimate = float(fixed["paired_mean_outcome_gain"]) * 1000
            low = float(row["two_stage_ci_2p5"]) * 1000
            high = float(row["two_stage_ci_97p5"]) * 1000
            yy = yi + offset
            ax.plot([low, high], [yy, yy], color=color, lw=1.5, solid_capstyle="round")
            ax.plot(
                estimate,
                yy,
                marker=marker,
                ms=4.2,
                mfc=color if filled else "white",
                mec=color,
                mew=0.8,
                linestyle="none",
                zorder=3,
            )
            source_rows.append(
                {
                    "panel": "d",
                    "series": label,
                    "category": f"{source_pool} source, r={selector_bits}",
                    "x": estimate,
                    "y": yy,
                    "low": low,
                    "high": high,
                    "n": int(row["replicates"]),
                    "unit": "relative PCR efficiency x 1e-3",
                    "note": (
                        "point is frozen full-data FullContext estimate; interval uses "
                        f"2,000 grouped-source plus target-fiber bootstrap replicates; {evidence_note}"
                    ),
                }
            )
    ax.axvline(0, color="#9A9A9A", lw=0.8, ls="--", zorder=0)
    ax.set_yticks(base_y)
    ax.set_yticklabels([f"{source}, r={selector_bits}" for source, selector_bits in order])
    ax.set_xlabel(r"FullContext measured selection gain ($\times 10^{-3}$)")
    ax.set_title("Matched assays are positive; Q5 altered protocol reverses", pad=4)
    ax.set_xlim(-30.5, 6.0)
    ax.set_ylim(-0.5, 3.5)
    ax.grid(axis="x", color=LIGHT_GREY, lw=0.6)
    ax.margins(x=0.08)
    ax.legend(
        handles=[
            Line2D([0], [0], color=BLUE, marker="o", mfc=BLUE, label="Cross-pool public"),
            Line2D(
                [0],
                [0],
                color=BLUE,
                marker="D",
                mfc="white",
                label="External KAPA",
            ),
            Line2D(
                [0],
                [0],
                color=Q5_RED,
                marker="v",
                mfc=Q5_RED,
                label="Q5 altered protocol",
            ),
        ],
        loc="center left",
        bbox_to_anchor=(0.02, 0.49),
        fontsize=4.6,
        handlelength=1.2,
        labelspacing=0.35,
    )
    panel_label(ax, "d", x=-0.08, y=1.08)
    return source_rows


def draw_seed_panel(ax: plt.Axes, sensitivity: pd.DataFrame) -> list[dict[str, Any]]:
    data = sensitivity.loc[
        sensitivity["selector_bits"].isin([2, 4])
        & sensitivity["selector_model"].isin(
            ["P5_combined_context", "P5_minus_P2_selected_outcome"]
        )
    ].copy()
    groups = [
        ("GCall_to_GCfix", 2, "GCall→GCfix, r=2"),
        ("GCall_to_GCfix", 4, "GCall→GCfix, r=4"),
        ("GCfix_to_GCall", 2, "GCfix→GCall, r=2"),
        ("GCfix_to_GCall", 4, "GCfix→GCall, r=4"),
    ]
    y = np.arange(len(groups))[::-1]
    source_rows: list[dict[str, Any]] = []
    for yi, (direction, selector_bits, label) in zip(y, groups):
        for model, offset, color, marker, series_label in (
            ("P5_combined_context", 0.12, PURPLE, "o", "FullContext gain"),
            ("P5_minus_P2_selected_outcome", -0.12, TEAL, "s", "FullContext − AssayContext"),
        ):
            row = data.loc[
                data["direction"].eq(direction)
                & data["selector_bits"].eq(selector_bits)
                & data["selector_model"].eq(model)
            ].iloc[0]
            lo = row["minimum_standardized_gain"]
            med = row["median_standardized_gain"]
            hi = row["maximum_standardized_gain"]
            yy = yi + offset
            ax.plot([lo, hi], [yy, yy], color=color, lw=1.5, solid_capstyle="round")
            ax.plot(med, yy, marker=marker, color=color, ms=3.7, mec="white", mew=0.35)
            source_rows.append(
                {
                    "panel": "c",
                    "series": series_label,
                    "category": label,
                    "x": med,
                    "y": yy,
                    "low": lo,
                    "high": hi,
                    "n": int(row["seeds"]),
                    "unit": "target outcome SD",
                    "note": "minimum-median-maximum over outcome-blind mappings; not a confidence interval",
                }
            )
    ax.axvline(0, color="#9A9A9A", lw=0.8, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([group[2] for group in groups])
    ax.set_xlabel("Standardized measured gain")
    ax.set_title("Positive gain persists across 32 mappings")
    ax.set_xlim(-0.025, 0.53)
    ax.set_ylim(-0.55, 3.55)
    ax.grid(axis="x", color=LIGHT_GREY, lw=0.6)
    ax.legend(
        handles=[
            Line2D([0], [0], color=PURPLE, marker="o", label="FullContext gain"),
            Line2D([0], [0], color=TEAL, marker="s", label="FullContext − AssayContext"),
        ],
        loc="lower right",
        handlelength=1.6,
    )
    panel_label(ax, "c")
    return source_rows


def draw_rate_utility_panel(ax: plt.Axes, generated: pd.DataFrame) -> list[dict[str, Any]]:
    data = generated.loc[generated["selector_model"].eq("P5_combined_context")].copy()
    source_rows: list[dict[str, Any]] = []
    for source_pool, color in (("GCall", BLUE), ("GCfix", ORANGE)):
        part = data.loc[data["selector_source_pool"].eq(source_pool)].sort_values("selector_bits")
        rate = part["payload_bits_per_variable_nt"].to_numpy(float)
        own = part["own_model_mean_gain_vs_fiber"].to_numpy(float) * 1000
        cross = part["cross_model_mean_gain_vs_fiber"].to_numpy(float) * 1000
        ax.plot(rate, own, "-o", color=color, label=f"{source_pool} selector: own model")
        ax.plot(rate, cross, "--D", color=color, alpha=0.78, ms=3.3, label=f"{source_pool} selector: other-pool model")
        for _, row in part.iterrows():
            if source_pool == "GCall" and int(row["selector_bits"]) in {0, 2, 4, 6}:
                ax.text(
                    row["payload_bits_per_variable_nt"],
                    row["own_model_mean_gain_vs_fiber"] * 1000 + 0.08,
                    f"r={int(row['selector_bits'])}",
                    color=color,
                    fontsize=4.8,
                    ha="center",
                    va="bottom",
                )
            for series, value in (
                ("own model", row["own_model_mean_gain_vs_fiber"] * 1000),
                ("other-pool model", row["cross_model_mean_gain_vs_fiber"] * 1000),
            ):
                source_rows.append(
                    {
                        "panel": "e",
                        "series": f"{source_pool} selector: {series}",
                        "category": f"r={int(row['selector_bits'])}",
                        "x": row["payload_bits_per_variable_nt"],
                        "y": value,
                        "low": "",
                        "high": "",
                        "n": int(row["sampled_payloads"]),
                        "unit": "x=exact payload bits/variable-region nt; y=model-score gain x 1e-3",
                        "note": "generated sequences; computational prediction only",
                    }
                )
    ax.set_xlabel("Exact retained payload rate\n(bits/variable-region nt)")
    ax.set_ylabel("Predicted score gain (x 10^-3)")
    ax.set_title("Generated predicted score gain has an exact rate cost")
    ax.set_xlim(1.907, 1.978)
    ax.set_ylim(-0.1, 5.0)
    ax.grid(color=LIGHT_GREY, lw=0.6)
    ax.legend(loc="upper right", handlelength=2.0, fontsize=4.7)
    panel_label(ax, "e", x=-0.08, y=1.08)
    return source_rows


def draw_failure_panel(ax: plt.Axes, rejection: pd.DataFrame) -> list[dict[str, Any]]:
    data = rejection.loc[rejection["selector_model"].eq("P5_combined_context")].copy()
    source_rows: list[dict[str, Any]] = []
    for source_pool, color in (("GCall", BLUE), ("GCfix", ORANGE)):
        part = data.loc[data["selector_source_pool"].eq(source_pool)].sort_values("selector_bits")
        x = part["selector_bits"].to_numpy(int)
        y = part["failure_fraction"].to_numpy(float)
        ax.plot(x, y, "-o", color=color, label=f"threshold, {source_pool}")
        for _, row in part.iterrows():
            source_rows.append(
                {
                    "panel": "e",
                    "series": f"threshold, {source_pool}",
                    "category": f"r={int(row['selector_bits'])}",
                    "x": int(row["selector_bits"]),
                    "y": row["failure_fraction"],
                    "low": "",
                    "high": "",
                    "n": int(row["sampled_payloads"]),
                    "unit": "payload failure fraction",
                    "note": "public top-quartile score threshold applied to generated candidates",
                }
            )
    x = np.arange(0, 7)
    ax.plot(x, np.zeros_like(x), color=DARK, ls=":", lw=1.3, label="fixed-choice codec")
    ax.text(3.1, 0.035, "choice codec: 0 failures", color=DARK, fontsize=5.1, ha="center")
    ax.set_xticks(x)
    ax.set_xlabel("Selector bits r")
    ax.set_ylabel("Payload failure fraction")
    ax.set_title("Fixed choice is total; threshold rejection is not")
    ax.set_ylim(-0.03, 0.83)
    ax.grid(axis="y", color=LIGHT_GREY, lw=0.6)
    ax.legend(loc="upper right", fontsize=4.8)
    top = ax.secondary_xaxis("top")
    top.set_xticks(x)
    top.set_xticklabels([str(2**int(value)) for value in x])
    top.set_xlabel("Candidate cap", labelpad=2)
    top.tick_params(length=2, pad=1)
    panel_label(ax, "e")
    for selector_bits in x:
        source_rows.append(
            {
                "panel": "e",
                "series": "fixed-choice codec",
                "category": f"r={selector_bits}",
                "x": selector_bits,
                "y": 0,
                "low": "",
                "high": "",
                "n": 512,
                "unit": "payload failure fraction",
                "note": "sampled generated payloads; exact construction is total on its declared domain",
            }
        )
    return source_rows


def write_source_data(
    rows: list[dict[str, Any]], path: Path = SOURCE_DATA
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["panel", "series", "category", "x", "y", "low", "high", "n", "unit", "note"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def draw_paired_two_stage_supplement(two_stage: pd.DataFrame) -> None:
    fig, (ax_a, ax_b) = plt.subplots(
        1,
        2,
        figsize=(7.01, 3.85),
        gridspec_kw={"width_ratios": [1.24, 1.0]},
    )
    fig.subplots_adjust(left=0.15, right=0.985, top=0.90, bottom=0.21, wspace=0.60)
    source_rows: list[dict[str, Any]] = []

    primary_order = [
        ("Gimpel2025_cross_pool_public_codebook", "GCall", 2, "GCall→GCfix, r=2"),
        ("Gimpel2025_cross_pool_public_codebook", "GCall", 4, "GCall→GCfix, r=4"),
        ("Gimpel2025_cross_pool_public_codebook", "GCfix", 2, "GCfix→GCall, r=2"),
        ("Gimpel2025_cross_pool_public_codebook", "GCfix", 4, "GCfix→GCall, r=4"),
        ("Gimpel2025_external_laboratory_Taq", "GCall", 2, "GCall→KAPA, r=2"),
        ("Gimpel2025_external_laboratory_Taq", "GCall", 4, "GCall→KAPA, r=4"),
        ("Gimpel2025_external_laboratory_Taq", "GCfix", 2, "GCfix→KAPA, r=2"),
        ("Gimpel2025_external_laboratory_Taq", "GCfix", 4, "GCfix→KAPA, r=4"),
    ]
    primary_styles = [
        ("P5_minus_P2_selected_outcome", 0.11, TEAL, "s", "FullContext − AssayContext"),
        (
            "P5_minus_released_1dcnn_selected_outcome",
            -0.11,
            ORANGE,
            "D",
            "FullContext − 1D-CNN score",
        ),
    ]
    y_a = np.arange(len(primary_order))[::-1]
    for yi, (dataset, source, selector_bits, category) in zip(y_a, primary_order):
        for estimand, offset, color, marker, label in primary_styles:
            row = two_stage.loc[
                two_stage["target_dataset"].eq(dataset)
                & two_stage["source_pool"].eq(source)
                & two_stage["selector_bits"].eq(selector_bits)
                & two_stage["estimand"].eq(estimand)
            ].iloc[0]
            estimate = float(row["two_stage_mean_estimate"]) * 1000
            low = float(row["two_stage_ci_2p5"]) * 1000
            high = float(row["two_stage_ci_97p5"]) * 1000
            yy = yi + offset
            ax_a.plot([low, high], [yy, yy], color=color, lw=1.35, solid_capstyle="round")
            ax_a.plot(estimate, yy, marker=marker, color=color, mec="white", mew=0.45, ms=4.0)
            source_rows.append(
                {
                    "panel": "S1a",
                    "series": label,
                    "category": category,
                    "x": estimate,
                    "y": yy,
                    "low": low,
                    "high": high,
                    "n": int(row["replicates"]),
                    "unit": "relative PCR efficiency x 1e-3",
                    "note": "two-stage source-sequence plus target-fiber percentile interval",
                }
            )
    ax_a.axvline(0, color="#9A9A9A", lw=0.8, ls="--")
    ax_a.set_yticks(y_a)
    ax_a.set_yticklabels([item[3] for item in primary_order])
    ax_a.set_xlabel(r"Paired selected-outcome difference ($\times 10^{-3}$)")
    ax_a.set_title("Continuous-efficiency contrasts\nremain positive", pad=5)
    ax_a.set_ylim(-0.55, len(primary_order) - 0.45)
    ax_a.set_xlim(-0.25, 3.80)
    ax_a.grid(axis="x", color=LIGHT_GREY, lw=0.6)
    ax_a.legend(
        handles=[
            Line2D([0], [0], color=color, marker=marker, label=label)
            for _estimand, _offset, color, marker, label in primary_styles
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=2,
        handlelength=1.4,
        columnspacing=0.9,
        fontsize=5.1,
    )
    panel_label(ax_a, "a", x=-0.20, y=1.05)

    q5_order = [
        ("GCall", 2, "GCall, r=2"),
        ("GCall", 4, "GCall, r=4"),
        ("GCfix", 2, "GCfix, r=2"),
        ("GCfix", 4, "GCfix, r=4"),
    ]
    q5_styles = [
        ("P5_combined_context_mean_gain", 0.18, PURPLE, "o", "FullContext gain"),
        ("P5_minus_P2_selected_outcome", 0.0, TEAL, "s", "FullContext − AssayContext"),
        (
            "P5_minus_released_1dcnn_selected_outcome",
            -0.18,
            ORANGE,
            "D",
            "FullContext − 1D-CNN score",
        ),
    ]
    y_b = np.arange(len(q5_order))[::-1]
    q5_dataset = "Gimpel2025_external_laboratory_Q5_sensitivity"
    for yi, (source, selector_bits, category) in zip(y_b, q5_order):
        for estimand, offset, color, marker, label in q5_styles:
            row = two_stage.loc[
                two_stage["target_dataset"].eq(q5_dataset)
                & two_stage["source_pool"].eq(source)
                & two_stage["selector_bits"].eq(selector_bits)
                & two_stage["estimand"].eq(estimand)
            ].iloc[0]
            estimate = float(row["two_stage_mean_estimate"]) * 1000
            low = float(row["two_stage_ci_2p5"]) * 1000
            high = float(row["two_stage_ci_97p5"]) * 1000
            yy = yi + offset
            ax_b.plot([low, high], [yy, yy], color=color, lw=1.35, solid_capstyle="round")
            ax_b.plot(estimate, yy, marker=marker, color=color, mec="white", mew=0.45, ms=4.0)
            source_rows.append(
                {
                    "panel": "S1b",
                    "series": label,
                    "category": category,
                    "x": estimate,
                    "y": yy,
                    "low": low,
                    "high": high,
                    "n": int(row["replicates"]),
                    "unit": "relative PCR efficiency x 1e-3",
                    "note": "secondary exploratory Q5 sensitivity; two-stage percentile interval",
                }
            )
    ax_b.axvline(0, color="#9A9A9A", lw=0.8, ls="--")
    ax_b.set_yticks(y_b)
    ax_b.set_yticklabels([item[2] for item in q5_order])
    ax_b.set_xlabel(r"Q5 gain or paired difference ($\times 10^{-3}$)")
    ax_b.set_title("Q5 workflow shift reverses\nthe KAPA-context result", pad=5)
    ax_b.set_ylim(-0.55, len(q5_order) - 0.45)
    ax_b.set_xlim(-31.5, 1.0)
    ax_b.grid(axis="x", color=LIGHT_GREY, lw=0.6)
    ax_b.legend(
        handles=[
            Line2D([0], [0], color=color, marker=marker, label=label)
            for _estimand, _offset, color, marker, label in q5_styles
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=1,
        handlelength=1.4,
        fontsize=5.0,
    )
    panel_label(ax_b, "b", x=-0.22, y=1.05)

    write_source_data(source_rows, SUPP_SOURCE_DATA)
    fig.savefig(SUPP_OUT_BASE.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(SUPP_OUT_BASE.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(SUPP_OUT_BASE.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(SUPP_OUT_BASE.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {SUPP_OUT_BASE}.{{svg,pdf,tiff,png}}")
    print(f"source_data={SUPP_SOURCE_DATA} rows={len(source_rows)}")


def main() -> None:
    generated = read("generated_fiber_results.tsv")
    benchmark = read_sota("fiber_benchmark_summary.tsv")
    two_stage = read_sota("two_stage_bootstrap_summary.tsv")

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(7.01, 6.10), constrained_layout=False)
    grid = fig.add_gridspec(
        3,
        2,
        height_ratios=[0.88, 1.18, 1.14],
        left=0.085,
        right=0.985,
        bottom=0.115,
        top=0.965,
        wspace=0.45,
        hspace=0.70,
    )
    ax_a = fig.add_subplot(grid[0, :])
    ax_b = fig.add_subplot(grid[1, 0])
    ax_c = fig.add_subplot(grid[1, 1])
    ax_d = fig.add_subplot(grid[2, 0])
    ax_e = fig.add_subplot(grid[2, 1])

    source_rows: list[dict[str, Any]] = []
    source_rows.extend(draw_method_panel(ax_a))
    source_rows.extend(
        draw_identical_fiber_benchmark(
            ax_b,
            benchmark,
            dataset="Gimpel2025_cross_pool_public_codebook",
            panel="b",
            title="Cross-pool benchmark on identical fibers\n(released CNN regression prediction)",
            methods=[
                "P5_combined_context",
                "P2_assay_context",
                "published_1dcnn",
                "published_sota_hard_filter",
            ],
            show_legend=True,
        )
    )
    source_rows.extend(
        draw_identical_fiber_benchmark(
            ax_c,
            benchmark,
            dataset="Gimpel2025_external_laboratory_Taq",
            panel="c",
            title=(
                "Analysis-plan-locked external KAPA benchmark\n"
                "(same publication; released CNN low-efficiency probability)"
            ),
            methods=["P5_combined_context", "P2_assay_context", "published_1dcnn"],
            show_legend=False,
        )
    )
    source_rows.extend(draw_two_stage_panel(ax_d, benchmark, two_stage))
    source_rows.extend(draw_rate_utility_panel(ax_e, generated))

    fig.text(
        0.5,
        0.022,
        "Panels b–d: retrospective measured selection; panel d contrasts matched assays with exploratory Q5 transfer. "
        "Panel e: generated model scores only; no emitted-codeword PCR experiment.",
        ha="center",
        va="bottom",
        fontsize=5.2,
        color=GREY,
    )
    write_source_data(source_rows)
    fig.savefig(OUT_BASE.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(OUT_BASE.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(OUT_BASE.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(OUT_BASE.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    draw_paired_two_stage_supplement(two_stage)
    print(f"wrote {OUT_BASE}.{{svg,pdf,tiff,png}}")
    print(f"source_data={SOURCE_DATA} rows={len(source_rows)}")


if __name__ == "__main__":
    main()
