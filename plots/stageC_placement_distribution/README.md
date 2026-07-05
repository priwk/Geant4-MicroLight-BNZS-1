# Stage C placement distribution figures

This folder contains a standalone plotting workflow for representative Stage C
placement-model distributions.

The script reads placement CSV files from:

- `Input/placements/<ratio>/**/*.csv`

It plots sphere centers only, because the placement files contain many
particles. BN and ZnS(Ag) centers are separated by color and marker.

## Run

```bash
cd ~/g4work/B2
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r plots/stageC_placement_distribution/requirements.txt
python3 plots/stageC_placement_distribution/generate_stageC_placement_distribution.py
```

By default, the script scans all ratio folders under `Input/placements` and uses
the first sorted placement CSV for each ratio. The built-in display order is:

- `2-1`
- `1-1`
- `1-1.5`
- `1-2`
- `1-2.5`
- `1-3`
- `1-3.5`
- `1-4`

Additional valid ratio folders are appended after these in numeric order.

Use another representative placement by passing a zero-based index:

```bash
python3 plots/stageC_placement_distribution/generate_stageC_placement_distribution.py \
  --placement-index 4
```

The 2D slice defaults to `z = 0 +/- 1 um`. To change it:

```bash
python3 plots/stageC_placement_distribution/generate_stageC_placement_distribution.py \
  --slice-axis z \
  --slice-center 0 \
  --slice-half-width 2
```

## Generated outputs

The script writes figures and a summary table to:

- `plots/stageC_placement_distribution/output/`

Main outputs:

- `stageC_placement_centers_3d_panels.png`
- `stageC_placement_center_slice_panels.png`
- `stageC_placement_distribution_summary.csv`

Useful options:

- `--max-points-per-phase`: cap plotted centers per phase in each panel. The
  default is `2500`; use `0` to plot every center.
- `--ratios`: select a subset or a custom order of ratio folders.
- `--placement-root`: read placement CSVs from a non-default directory.
