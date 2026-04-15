# CMML3_ICA1_Morris_water_maze
# Reproducibility (ICA1 Water Maze Models)

This folder contains three original model scripts and a set of wrapper scripts used to reproduce the reported results and figures.

## Prerequisites

- R (tested with R 4.5.3)
- R package: `ggplot2` (only needed for figure generation)

## What not to edit

The original model scripts are kept unchanged:

- `place-cell model.R`
- `distance-cell model.R`
- `place distance cells combined model.R`

All experiments are run via wrappers that read these scripts, override parameters (e.g., `Nruns`, plotting flags), execute them, and then save standardized CSV outputs.

## Baseline results (6 conditions, Nruns = 50)

Run the six wrappers below (they overwrite the corresponding `*_batch.csv` files if they already exist):

```bash
Rscript place_fixed.R
Rscript place_variable.R
Rscript distance_fixed.R
Rscript distance_variable.R
Rscript combined_fixed_w02.R
Rscript combined_variable_w02.R
```

Expected outputs:

- `place_fixed_batch.csv`
- `place_variable_batch.csv`
- `distance_fixed_batch.csv`
- `distance_variable_batch.csv`
- `combined_fixed_w02_batch.csv`
- `combined_variable_w02_batch.csv`

Each CSV has 8 rows (days 1–8) and columns:

`day, latency, dist, target_quadrant, opposite_quadrant, wall_zone`

## Combined ratio extension (Nruns = 50)

These additional wrappers vary the combined model balance:

```bash
Rscript combined_fixed_w05.R
Rscript combined_variable_w05.R
Rscript combined_fixed_w08.R
Rscript combined_variable_w08.R
```

Expected outputs:

- `combined_fixed_w05_batch.csv`
- `combined_variable_w05_batch.csv`
- `combined_fixed_w08_batch.csv`
- `combined_variable_w08_batch.csv`

## Figures (mean curves)

After the CSVs exist, generate figures with:

```bash
Rscript make_main_figures.R
```

Expected outputs:

- `fig1_latency_vs_day.pdf`
- `fig2_wall_zone_vs_day.pdf`
- `fig3_combined_ratio_tradeoff.pdf` (generated only if all ratio CSVs exist)

## Run-level export and 95% CI figures

To quantify between-run variability, export a run-level long table and re-draw figures with 95% confidence intervals.

1) Export run-level plot data:

```bash
Rscript export_plot_data.R --nruns=50 --out=plot_data_runs.csv
```

This produces a long-format table with columns:

`task, model, weight_wall, day, run, latency, wall_zone, target_quadrant`

2) Plot figures with 95% CI (uses `plot_data_runs.csv` if present):

```bash
Rscript plot_figures_with_variability.R
```

Expected outputs:

- `fig1_latency_vs_day_ci.pdf`
- `figS1_wall_zone_vs_day_ci.pdf`
- `fig2_combined_ratio_tradeoff_ci.pdf`

## Example trajectories (qualitative)

Representative fixed-task trajectories can be generated with:

```bash
Rscript make_example_trajectories.R
```

Expected output:

- `figS2_example_trajectories_fixed.pdf`

## Notes

- All wrappers call `set.seed(1)` for reproducibility.
- Plotting inside the original model scripts is disabled in batch wrappers to avoid slowdowns.

