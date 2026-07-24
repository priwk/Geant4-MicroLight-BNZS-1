# Stage B scintillation-source figures

This folder contains a standalone plotting workflow for Stage B ZnS(Ag)
effective energy-deposition metrics.

The script reads:

- `Output/stageB/*/*_capture_anchors.csv`
- `Output/stageB/*/*_zns_track_steps.csv`
- `Input/stageA/*/neutron_capture_absorption/*.csv`
- `Input/stageA/*/neutron_capture_positions/*.csv`

The Stage B slim files are interpreted as replay-weighted trajectory samples.
`capture_anchors.csv` defines the full denominator, including trajectories with
zero ZnS(Ag) deposition. `zns_track_steps.csv` contributes ZnS(Ag) deposition
from alpha, Li-7, and their transported secondaries. This is important for
Fig. 3(b): zero-deposition capture outcomes remain in the average.

## Run

```bash
cd ~/g4work/B2
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r plots/stageB_scintillation_strength/requirements.txt
python3 plots/stageB_scintillation_strength/generate_stageB_scintillation_strength.py
```

If Birks quenching should be applied to the step energy before averaging:

```bash
python3 plots/stageB_scintillation_strength/generate_stageB_scintillation_strength.py \
  --birks-kb 0.01
```

`--birks-kb` is in `um/keV` because `dE/dx` is computed from `keV/um`.
The default is `0`, so effective energy equals the raw ZnS(Ag) deposited
energy.

## Generated outputs

The script writes derived CSV tables and figures to:

- `plots/stageB_scintillation_strength/output/`

Main figures:

- `fig3_stageB_scintillation_metrics.png`
- `fig3a_effective_deposition_fraction_vs_thickness.png`
- `fig3a_effective_deposition_fraction_heatmap.png`
- `fig3b_mean_effective_edep_vs_thickness.png`
- `fig3b_mean_effective_edep_by_ratio_with_placement_dispersion.png`
- `fig3c_relative_initial_scintillation_strength_vs_thickness.png`
- `fig3c_relative_initial_scintillation_strength_heatmap.png`

Main tables:

- `stageB_effective_deposition_metrics_by_ratio_and_thickness.csv`
- `stageB_effective_deposition_metrics_by_placement.csv`
