#!/usr/bin/env python
"""Plot cell-level colocalization-event counts for the SCRIN example."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


EXAMPLE_DIR = Path(__file__).resolve().parents[1]
DEFAULT_DATA = EXAMPLE_DIR / "data" / "Mouse_brain_CosMx_787cells_Gabra2_Gabrb1_cell_events.csv"
DEFAULT_OUTPUT = EXAMPLE_DIR / "figures" / "gabra2-gabrb1-tissue-projection.png"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot per-cell colocalization-event counts on tissue coordinates."
    )
    parser.add_argument("--event-path", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument(
        "--no-flip-x",
        action="store_false",
        dest="flip_x",
        help="Retain the input x-axis orientation instead of flipping it.",
    )
    parser.set_defaults(flip_x=True)
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


def main() -> None:
    args = parse_args()
    events = pd.read_csv(args.event_path)
    required = {"cell_x", "cell_y", "gene_A", "gene_B", "colocalization_count"}
    missing_columns = required.difference(events.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    pairs = events[["gene_A", "gene_B"]].drop_duplicates()
    if len(pairs) != 1:
        raise ValueError("The event file must contain exactly one undirected gene pair.")
    gene_a, gene_b = pairs.iloc[0]

    x_coordinates = -events["cell_x"] if args.flip_x else events["cell_x"]
    events = events.assign(
        x_relative=x_coordinates - x_coordinates.min(),
        y_relative=events["cell_y"] - events["cell_y"].min(),
    )
    positive = events[events["colocalization_count"] > 0]
    if positive.empty:
        raise ValueError("The event file does not contain any positive cells.")

    configure_matplotlib()
    fig, ax = plt.subplots(figsize=(6.2, 5.5), constrained_layout=True)
    ax.scatter(
        events["x_relative"],
        events["y_relative"],
        s=13,
        c="#dedede",
        edgecolors="none",
    )
    points = ax.scatter(
        positive["x_relative"],
        positive["y_relative"],
        s=21,
        c=positive["colocalization_count"],
        cmap="viridis",
        vmin=1,
        vmax=positive["colocalization_count"].max(),
        edgecolors="none",
    )
    colorbar = fig.colorbar(points, ax=ax, pad=0.025, fraction=0.045)
    colorbar.set_label("Colocalization events per cell")
    colorbar.set_ticks(range(1, int(positive["colocalization_count"].max()) + 1))

    ax.set_xlabel("Relative x coordinate (µm)")
    ax.set_ylabel("Relative y coordinate (µm)")
    ax.set_title(f"Cell-level {gene_a}–{gene_b} colocalization events", pad=10)
    ax.set_aspect("equal", adjustable="box")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(False)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        args.output,
        dpi=args.dpi,
        facecolor="white",
        bbox_inches="tight",
        pad_inches=0.12,
    )
    plt.close(fig)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
