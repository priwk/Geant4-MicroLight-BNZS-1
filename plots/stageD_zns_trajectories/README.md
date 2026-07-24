# Stage D alpha/Li trajectory figures

This folder contains a standalone plotting workflow for Stage B alpha/Li
particle trajectories. It is intended for Stage D-style trajectory
visualization and cross-ratio summary figures.

The script reads:

- `Output/stageB/<ratio>/*_alpha_li_steps.csv`
- `Output/stageB/<ratio>/*_zns_track_steps.csv`

For each ratio/thickness pair, full `*_alpha_li_steps.csv` files are preferred
because they contain all recorded alpha/Li steps and all phases. If a full file
is not available, the script falls back to `*_zns_track_steps.csv`, which is
ZnS-only.

It is written to scan all available ratio folders. A local checkout may contain
only `Output/stageB/1-2`; the same code will process all ratios on a server with
complete Stage B outputs.

## Run

```bash
cd ~/g4work/B2
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r plots/stageD_zns_trajectories/requirements.txt
python3 plots/stageD_zns_trajectories/generate_stageD_zns_trajectories.py
```

For a quick smoke test on large data:

```bash
python3 plots/stageD_zns_trajectories/generate_stageD_zns_trajectories.py \
  --max-files 2 \
  --summary-stride 20 \
  --max-showcase-steps 20000
```

## Generated outputs

The script writes figures and a derived summary table to:

- `plots/stageD_zns_trajectories/output/`

Main outputs:

- `stageD_zns_trajectory_showcase_3d.png`
- `stageD_origin_aligned_trajectory_showcase_3d.png`
- `stageD_macro_trajectory_showcase_3d.png`
- `stageD_zns_trajectory_projection_panels.png`
- `stageD_zns_trajectory_statistics_vs_thickness.png`
- `stageD_1000um_mean_total_track_length_by_ratio.png`
- `stageD_zns_trajectory_summary.csv`

The showcase figures select representative particle trajectories from one track
file. Full files are plotted with segment color set by material phase:

- BN phase
- ZnS phase
- Binder/porosity phase
- Outside domain

Alpha and Li7 are separated by line width. By default, the showcase uses the
largest processed thickness for the first processed ratio. You can select a
specific case:

```bash
python3 plots/stageD_zns_trajectories/generate_stageD_zns_trajectories.py \
  --showcase-ratio 1-2 \
  --showcase-thickness 1000 \
  --sample-trajectories 120 \
  --max-showcase-steps 800000
```

New Stage B files provide continuous `unwrapped_*` coordinates, so
`stageD_zns_trajectory_showcase_3d.png` no longer folds trajectories back into
one RVE after a periodic crossing. The original `x_*` columns remain folded
local coordinates for geometry diagnostics. If the selected input contains
explicit `screen_*` columns, the script uses them directly for
`stageD_macro_trajectory_showcase_3d.png`. That figure maps each local step back
to the macroscopic capture location. For legacy files it falls back to:

```text
macro_x = capture_x_um + (local_x - local_capture_x_um)
macro_y = capture_y_um + (local_y - local_capture_y_um)
macro_depth = depth_um - (local_z - local_capture_z_um)
```

`stageD_origin_aligned_trajectory_showcase_3d.png` is a display-only transform:
each selected particle trajectory is translated so its first recorded point is
at `(0, 0, 0)`. This does not modify the source CSV. Alpha and Li7 trajectories
are shown as separate 3D panels, with phase colors preserved and equal axes so
the approximate travel distance is visible at a glance.

For the cloud server case where all six ratios have `1000_alpha_li_steps.csv`,
run only the 1000 um files:

```bash
python3 plots/stageD_zns_trajectories/generate_stageD_zns_trajectories.py \
  --thicknesses 1000 \
  --showcase-ratio 1-2 \
  --showcase-thickness 1000 \
  --sample-trajectories 160 \
  --reservoir-size 2000 \
  --max-showcase-steps 1000000
```

To regenerate only `stageD_origin_aligned_trajectory_showcase_3d.png`, skip
the summary and other figures:

```bash
python3 plots/stageD_zns_trajectories/generate_stageD_zns_trajectories.py \
  --only-origin-aligned \
  --ratios 1-2 \
  --thicknesses 1000 \
  --showcase-ratio 1-2 \
  --showcase-thickness 1000 \
  --sample-trajectories 160 \
  --reservoir-size 2000 \
  --max-showcase-steps 1000000
```

The `stageD_1000um_mean_total_track_length_by_ratio.png` figure uses only
`t = 1000 um` records and plots the mean total recorded trajectory length for
alpha and Li7 separately. With full input this is the full recorded alpha/Li
path length. With slim fallback data it is labeled as ZnS-only.

Useful options:

- `--ratios`: process selected ratio folders while keeping the code generic.
- `--thicknesses`: process selected thicknesses only; use `--thicknesses 1000`
  for the six-ratio 1000 um comparison figure.
- `--summary-stride`: read every Nth row for faster approximate diagnostics.
  Use the default `1` for exact local summaries.
- `--sample-trajectories`: number of source trajectories shown in the 3D and
  projection figures.
- `--reservoir-size`: candidate pool for the showcase trajectory sampler.
- `--max-showcase-steps`: maximum row count scanned for the showcase sampler.
- `--only-origin-aligned`: regenerate only the origin-aligned 3D showcase image.
