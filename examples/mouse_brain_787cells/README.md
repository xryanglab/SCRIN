# Mouse brain CosMx example (787 cells)

This example uses a spatially contiguous, approximately 500 × 500 µm mouse-brain region containing 787 cells. Transcript coordinates and all distance parameters are in micrometers. The full transcript CSV is available from [Zenodo](https://zenodo.org/records/21486759).

## Contents

| Path | Description |
| --- | --- |
| `data/Mouse_brain_CosMx_787cells_hyper_result_cb_dedup_1e-05_post_proc.csv` | Final SCRIN colocalization result after directional deduplication and q-value filtering. |
| `data/Mouse_brain_CosMx_787cells_selected_pair_distances.parquet` | Distance observations for the undirected Apoe-Clu and Gabra2-Gabrb1 pairs, retained within `r_dist = 1.0 µm`. |
| `data/Mouse_brain_CosMx_787cells_Gabra2_Gabrb1_cell_events.csv` | Per-cell Gabra2-Gabrb1 colocalization-event counts and median cell coordinates, calculated with `r_check = 0.5 µm`. |
| `scripts/plot_distance_distributions.py` | Recreates the selected-pair KDE curves and the two-dimensional random reference. |
| `scripts/plot_cell_colocalization_map.py` | Recreates the cell-level tissue projection. |
| `figures/*.png` | Figures displayed below. |

The full transcript CSV and SCRIN intermediate directories are intentionally omitted. They are not required to run the plotting scripts.

## Recreate the figures

Install `numpy`, `pandas`, `scipy`, `matplotlib`, and a Parquet engine such as `pyarrow`, then run:

```bash
python scripts/plot_distance_distributions.py
python scripts/plot_cell_colocalization_map.py
```

Both scripts resolve their default inputs relative to this example directory, so they can be launched from any working directory. Use `--help` to provide another input or output path. The tissue-projection script flips the x axis by default to match the displayed tissue orientation; pass `--no-flip-x` to retain the input orientation.

## Distance distributions

The observed KDE curves show different distance preferences for Apoe-Clu and Gabra2-Gabrb1. The dashed line is the theoretical radial density expected under uniform two-dimensional placement within `R = 1.0 µm`, `f(r) = 2r/R²`.

![Selected-pair distance distributions](figures/distance-distribution-kde.png)

## Cell-level tissue projection

Cells are positioned using the median x/y coordinates of all transcripts assigned to each cell. Cells with positive Gabra2-Gabrb1 colocalization-event counts are colored by event count; cells without detected events are shown in light gray.

![Gabra2-Gabrb1 cell-level tissue projection](figures/gabra2-gabrb1-tissue-projection.png)
