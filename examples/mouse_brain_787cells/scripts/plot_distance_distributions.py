#!/usr/bin/env python
"""Plot selected RNA-pair distance distributions for the SCRIN example."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


EXAMPLE_DIR = Path(__file__).resolve().parents[1]
DEFAULT_DATA = EXAMPLE_DIR / "data" / "Mouse_brain_CosMx_787cells_selected_pair_distances.parquet"
DEFAULT_OUTPUT = EXAMPLE_DIR / "figures" / "distance-distribution-kde.png"
DEFAULT_PAIRS = ["Apoe_Clu", "Gabra2_Gabrb1"]
COLORS = ["#d55e00", "#0072b2", "#009e73", "#cc79a7", "#f0e442"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot KDE curves for selected undirected RNA pairs."
    )
    parser.add_argument("--distance-path", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--pairs", nargs="+", default=DEFAULT_PAIRS)
    parser.add_argument("--r-dist", type=float, default=1.0)
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 12,
            "axes.titlesize": 12,
            "axes.labelsize": 12,
            "legend.fontsize": 12,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def pair_label(pair: str) -> str:
    return pair.replace("_", "–", 1)


def main() -> None:
    args = parse_args()
    if args.r_dist <= 0:
        raise ValueError("--r-dist must be positive.")

    distances = pd.read_parquet(args.distance_path)
    required = {"pair", "distance"}
    missing_columns = required.difference(distances.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    configure_matplotlib()
    x_grid = np.linspace(0, args.r_dist, 600)
    fig, ax = plt.subplots(figsize=(6.4, 5.0), constrained_layout=True)

    for index, pair in enumerate(args.pairs):
        values = distances.loc[distances["pair"] == pair, "distance"].to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        values = values[(values >= 0) & (values <= args.r_dist)]
        if len(values) < 2 or np.unique(values).size < 2:
            raise ValueError(f"Pair {pair!r} does not have enough variable distances for KDE.")

        density = gaussian_kde(values)(x_grid)
        area = np.trapz(density, x_grid)
        if area > 0:
            density = density / area
        ax.plot(
            x_grid,
            density,
            color=COLORS[index % len(COLORS)],
            linewidth=2.0,
            label=f"{pair_label(pair)} (n={len(values):,})",
        )

    ax.plot(
        x_grid,
        2 * x_grid / (args.r_dist**2),
        color="#4d4d4d",
        linewidth=1.8,
        linestyle="--",
        label="Random expectation",
    )
    ax.set_xlim(0, args.r_dist)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Pairwise molecule distance (µm)")
    ax.set_ylabel("Density")
    ax.grid(False)
    ax.legend(frameon=False, loc="best")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, facecolor="white")
    plt.close(fig)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
