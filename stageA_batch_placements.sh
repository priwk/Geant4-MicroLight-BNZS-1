#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="$PROJECT_ROOT/build"
PLACEMENT_ROOT="$PROJECT_ROOT/Input/output_pbc"
if [ ! -d "$PLACEMENT_ROOT" ]; then
  PLACEMENT_ROOT="$PROJECT_ROOT/Input/placements"
fi
OUTPUT_ROOT="$PROJECT_ROOT/Output/stageA"
LOG_DIR="$PROJECT_ROOT/logs/stageA"
EVENTS="${STAGEA_EVENTS:-100000}"

project_relpath() {
  python3 - "$PROJECT_ROOT" "$1" <<'PY'
from pathlib import Path
import sys
root = Path(sys.argv[1]).resolve()
target = Path(sys.argv[2]).resolve()
print(target.relative_to(root).as_posix())
PY
}

placement_tag() {
  python3 - "$1" "$2" <<'PY'
from pathlib import Path
import sys
root = Path(sys.argv[1]).resolve()
target = Path(sys.argv[2]).resolve()
rel = target.relative_to(root)
parts = list(rel.parts)
parts[-1] = rel.stem
print("__".join(parts))
PY
}

mkdir -p "$BUILD_DIR" "$LOG_DIR"

echo "=== Build ==="
cd "$BUILD_DIR"
cmake ..
make -j"$(nproc)"

if [ "$#" -gt 0 ]; then
  ratios=("$@")
else
  shopt -s nullglob
  ratios=()
  for dir in "$PLACEMENT_ROOT"/*; do
    [ -d "$dir" ] && ratios+=("$(basename "$dir")")
  done
fi

if [ "${#ratios[@]}" -eq 0 ]; then
  echo "No placement ratio folders found in: $PLACEMENT_ROOT"
  exit 1
fi

for ratio in "${ratios[@]}"; do
  ratio_dir="$PLACEMENT_ROOT/$ratio"
  if [ ! -d "$ratio_dir" ]; then
    echo "Missing placement ratio folder: $ratio_dir"
    exit 1
  fi

  bn_wt="${ratio%%-*}"
  zns_wt="${ratio#*-}"

  mapfile -t placements < <(
    find "$ratio_dir" -type f \( -name '*.csv' -o -name '*.txt' \) \
      ! -name '*_pbc_images.csv' ! -name '*_radius_stats.csv' | sort
  )
  if [ "${#placements[@]}" -eq 0 ]; then
    echo "No placement files found in: $ratio_dir"
    continue
  fi

  summary="$OUTPUT_ROOT/$ratio/neutron_transport_summary.csv"
  mkdir -p "$(dirname "$summary")"
  rm -f "$summary"

  echo
  echo "=== Stage A ratio $ratio: ${#placements[@]} placement files, $EVENTS events each ==="

  for placement in "${placements[@]}"; do
    base="$(basename "$placement")"
    tag="$(placement_tag "$ratio_dir" "$placement")"
    macro="$BUILD_DIR/stageA_batch_current.mac"
    log="$LOG_DIR/${ratio}_${tag}.log"
    placement_rel="$(project_relpath "$placement")"

    cat > "$macro" <<EOF
/control/verbose 2
/run/verbose 1
/event/verbose 0
/tracking/verbose 0

/cfg/setRunMode StageA_NeutronPatch
/cfg/setWeightRatio $bn_wt $zns_wt
/cfg/setPlacementFile $placement_rel

/run/initialize
/run/beamOn $EVENTS
EOF

    echo ">>> $ratio / $base"
    BNZS_RUN_MODE=StageA_NeutronPatch \
      BNZS_PLACEMENT_FILE="$placement_rel" \
      BNZS_USE_RANDOM_PLACEMENT=0 \
      ./Geant4-MicroLight-BNZS "$macro" > "$log" 2>&1
  done

  echo "Summary: $summary"
done

echo
echo "=== Stage A placement batch done ==="
