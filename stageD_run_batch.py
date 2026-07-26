#!/usr/bin/env python3
import argparse
import random
import subprocess
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate and optionally run StageD_OpticalHomogenization macros "
            "for placement CSV files under one BN/ZnS ratio."
        )
    )
    parser.add_argument(
        "--ratio-tag",
        required=True,
        help="Ratio folder under Input/output_pbc, e.g. 1-2 or 2-1.",
    )
    parser.add_argument(
        "--project-root",
        default=None,
        help="Project root. Defaults to the directory containing this script.",
    )
    parser.add_argument(
        "--build-dir",
        default=None,
        help="Build directory containing Geant4-MicroLight-BNZS. Defaults to <project-root>/build.",
    )
    parser.add_argument(
        "--placement-dir",
        default=None,
        help="Placement directory. Defaults to Input/output_pbc/<ratio-tag>.",
    )
    parser.add_argument(
        "--placements",
        nargs="*",
        default=None,
        help=(
            "Specific placement basenames or stems to run, e.g. "
            "placement_f_0.64_0004.csv placement_f_0.64_0008. "
            "If omitted, all placements in placement-dir are candidates."
        ),
    )
    parser.add_argument(
        "--max-placements",
        type=int,
        default=0,
        help="Run only the first N selected placements. Use 0 for all.",
    )
    parser.add_argument(
        "--start-index",
        type=int,
        default=0,
        help="Skip the first N sorted selected placements.",
    )
    parser.add_argument(
        "--shuffle",
        action="store_true",
        help="Shuffle placement order before applying start-index/max-placements.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260507,
        help="Shuffle seed when --shuffle is used.",
    )
    parser.add_argument(
        "--beam-on",
        type=int,
        default=10000,
        help="Geant4 /run/beamOn value for every placement.",
    )
    parser.add_argument(
        "--wavelength-nm",
        type=float,
        default=450.0,
        help="StageD photon wavelength in nm.",
    )
    parser.add_argument(
        "--optical-params",
        default=None,
        help=(
            "Six values: matrix_n matrix_abs_um bn_n bn_abs_um zns_n zns_abs_um. "
            "If omitted, StageD uses the code-side defaults/estimates."
        ),
    )
    parser.add_argument(
        "--source-mode",
        default="uniform_all_phase",
        choices=["uniform_ZnS", "uniform_all_phase"],
        help="StageD source mode.",
    )
    parser.add_argument(
        "--boundary-mode",
        default="periodic_wrap",
        choices=["periodic_wrap", "escape", "same_phase_reentry"],
        help="StageD boundary mode.",
    )
    parser.add_argument(
        "--reentry-mode",
        default="state_matched",
        choices=["state_matched", "same_phase_rho_over_R", "same_phase_random"],
        help="StageD same-phase re-entry mode.",
    )
    parser.add_argument(
        "--particle-reentry-mode",
        default="sphere_q_mu",
        choices=["sphere_q_mu", "same_phase_rho_over_R", "same_phase_random"],
        help="StageD particle re-entry mode.",
    )
    parser.add_argument(
        "--matrix-reentry-mode",
        default="clearance_binned_portal",
        choices=["clearance_binned_portal", "random_matrix_debug", "distance_matched_matrix", "random_matrix"],
        help="StageD matrix re-entry mode.",
    )
    parser.add_argument(
        "--scatter-metric",
        default="particle_encounter_no_threshold",
        choices=[
            "angle_threshold",
            "particle_encounter_angle_threshold",
            "particle_encounter_no_threshold",
        ],
        help="StageD primary scatter metric. Formal homogenization uses particle_encounter_no_threshold; angle_threshold is diagnostic.",
    )
    parser.add_argument(
        "--target-primary-scatter",
        type=int,
        default=0,
        help="Stop a photon once it accumulates this many primary scatter events. Use 0 to disable.",
    )
    parser.add_argument(
        "--theta-threshold-deg",
        type=float,
        default=0.0,
        help="Diagnostic effective-scatter angle threshold; raw outputs are always unthresholded.",
    )
    parser.add_argument(
        "--max-reentry",
        type=int,
        default=10000,
        help="StageD maximum re-entry count per photon.",
    )
    parser.add_argument(
        "--max-steps",
        type=int,
        default=100000,
        help="StageD maximum Geant4 step count per photon.",
    )
    parser.add_argument(
        "--max-path-length-um",
        type=float,
        default=5000.0,
        help="StageD maximum physical path length per photon in um.",
    )
    parser.add_argument(
        "--write-event-csv",
        action="store_true",
        help="Write sampled per-photon Stage D event rows. Disabled by default.",
    )
    parser.add_argument(
        "--event-sampling-rate",
        type=float,
        default=1.0,
        help="Deterministic event CSV sampling rate in [0,1].",
    )
    parser.add_argument(
        "--max-event-rows",
        type=int,
        default=100000,
        help="Maximum event CSV rows when --write-event-csv is enabled.",
    )
    parser.add_argument(
        "--macro-dir",
        default=None,
        help="Directory for generated macros. Defaults to Output/stageD_macros/<ratio-tag>.",
    )
    parser.add_argument(
        "--log-dir",
        default=None,
        help="Directory for run logs. Defaults to logs/stageD/<ratio-tag>.",
    )
    parser.add_argument(
        "--run",
        action="store_true",
        help="Actually run Geant4-MicroLight-BNZS for each generated macro.",
    )
    parser.add_argument(
        "--executable-name",
        default="Geant4-MicroLight-BNZS",
        help="Executable name inside build-dir. Default: Geant4-MicroLight-BNZS.",
    )
    return parser.parse_args()


def as_posix(path: Path) -> str:
    return path.as_posix()


def rel_to_project(path: Path, project_root: Path) -> str:
    return path.resolve().relative_to(project_root.resolve()).as_posix()


def ratio_parts(ratio_tag: str):
    if "-" not in ratio_tag:
        raise SystemExit(f"Invalid ratio tag: {ratio_tag}")
    lhs, rhs = ratio_tag.split("-", 1)
    try:
        float(lhs)
        float(rhs)
    except ValueError as exc:
        raise SystemExit(f"Invalid ratio tag: {ratio_tag}") from exc
    return lhs, rhs


def normalize_requested_placement(raw: str) -> str:
    return raw if raw.endswith(".csv") else raw + ".csv"


def is_main_placement_csv(path: Path) -> bool:
    name = path.name
    return (
        path.suffix.lower() == ".csv"
        and not name.endswith("_pbc_images.csv")
        and not name.endswith("_radius_stats.csv")
    )


def placement_display_tag(path: Path, placement_dir: Path) -> str:
    rel = path.resolve().relative_to(placement_dir.resolve())
    suffix = rel.suffix
    parts = list(rel.parts)
    if suffix:
        parts[-1] = rel.stem
    return "__".join(parts)


def resolve_placements(placement_dir: Path, requested):
    if not requested:
        return sorted(
            path for path in placement_dir.rglob("*.csv")
            if is_main_placement_csv(path)
        )

    resolved = []
    for item in requested:
        name = normalize_requested_placement(item).replace("\\", "/")
        candidate = placement_dir / name
        if candidate.is_file() and is_main_placement_csv(candidate):
            resolved.append(candidate)
            continue

        matches = sorted(
            path for path in placement_dir.rglob("*.csv")
            if is_main_placement_csv(path) and path.name == Path(name).name
        )
        if not matches:
            raise SystemExit(f"Placement not found under {placement_dir}: {item}")
        if len(matches) > 1:
            raise SystemExit(
                f"Placement name is ambiguous under {placement_dir}: {item}. "
                "Pass a relative subpath from placement-dir."
            )
        resolved.append(matches[0])
    return resolved


def macro_text(
    placement_rel_to_build: str,
    bn_wt: str,
    zns_wt: str,
    beam_on: int,
    optical_params: str | None,
    wavelength_nm: float,
    source_mode: str,
    boundary_mode: str,
    reentry_mode: str,
    particle_reentry_mode: str,
    matrix_reentry_mode: str,
    scatter_metric: str,
    target_primary_scatter: int,
    theta_threshold_deg: float,
    max_reentry: int,
    max_steps: int,
    max_path_length_um: float,
    write_event_csv: bool,
    event_sampling_rate: float,
    max_event_rows: int,
    random_seed: int,
):
    return "\n".join(
        [
            "/run/verbose 0",
            "/event/verbose 0",
            "/tracking/verbose 0",
            f"/random/setSeeds {random_seed} {random_seed + 1}",
            "",
            "/cfg/setRunMode StageD_OpticalHomogenization",
            f"/cfg/setWeightRatio {bn_wt} {zns_wt}",
            f"/cfg/setPlacementFile {placement_rel_to_build}",
            *([f"/cfg/setOpticalParams {optical_params}"] if optical_params else []),
            f"/cfg/stageD/setWavelengthNm {wavelength_nm}",
            f"/cfg/stageD/setSourceMode {source_mode}",
            f"/cfg/stageD/setBoundaryMode {boundary_mode}",
            f"/cfg/stageD/setReentryMode {reentry_mode}",
            f"/cfg/stageD/setParticleReentryMode {particle_reentry_mode}",
            f"/cfg/stageD/setMatrixReentryMode {matrix_reentry_mode}",
            f"/cfg/stageD/setScatterMetric {scatter_metric}",
            f"/cfg/stageD/setTargetPrimaryScatter {target_primary_scatter}",
            f"/cfg/stageD/setThetaThresholdDeg {theta_threshold_deg}",
            f"/cfg/stageD/setMaxReentry {max_reentry}",
            f"/cfg/stageD/setMaxSteps {max_steps}",
            f"/cfg/stageD/setMaxPathLengthUm {max_path_length_um}",
            f"/cfg/stageD/setWriteEventCsv {'true' if write_event_csv else 'false'}",
            f"/cfg/stageD/setEventSamplingRate {event_sampling_rate}",
            f"/cfg/stageD/setMaxEventRows {max_event_rows}",
            "",
            "/run/initialize",
            f"/run/beamOn {beam_on}",
            "",
        ]
    )


def main():
    args = parse_args()
    if args.max_placements < 0 or args.start_index < 0:
        raise SystemExit("--max-placements and --start-index must be >= 0")
    if args.seed <= 0:
        raise SystemExit("--seed must be > 0")
    if args.beam_on <= 0:
        raise SystemExit("--beam-on must be > 0")
    if args.wavelength_nm <= 0.0:
        raise SystemExit("--wavelength-nm must be > 0")
    if args.theta_threshold_deg < 0.0:
        raise SystemExit("--theta-threshold-deg must be >= 0")
    if args.target_primary_scatter < 0:
        raise SystemExit("--target-primary-scatter must be >= 0")
    if args.max_reentry <= 0 or args.max_steps <= 0 or args.max_path_length_um <= 0.0:
        raise SystemExit("--max-reentry, --max-steps, and --max-path-length-um must be > 0")
    if not 0.0 <= args.event_sampling_rate <= 1.0:
        raise SystemExit("--event-sampling-rate must be in [0,1]")
    if args.max_event_rows < 0:
        raise SystemExit("--max-event-rows must be >= 0")

    wavelength_tag = (
        f"lambda_{int(round(args.wavelength_nm))}nm"
        if abs(args.wavelength_nm - round(args.wavelength_nm)) < 1.0e-9
        else f"lambda_{args.wavelength_nm:.3f}nm"
    )
    project_root = (
        Path(args.project_root).resolve()
        if args.project_root
        else Path(__file__).resolve().parent
    )
    build_dir = (
        Path(args.build_dir).resolve()
        if args.build_dir
        else project_root / "build"
    )
    placement_dir = (
        Path(args.placement_dir).resolve()
        if args.placement_dir
        else project_root / "Input" / "output_pbc" / args.ratio_tag
    )
    macro_dir = (
        Path(args.macro_dir).resolve()
        if args.macro_dir
        else project_root / "Output" / "stageD_macros" / args.ratio_tag / wavelength_tag
    )
    log_dir = (
        Path(args.log_dir).resolve()
        if args.log_dir
        else project_root / "logs" / "stageD" / args.ratio_tag / wavelength_tag
    )
    executable = build_dir / args.executable_name

    if not build_dir.is_dir():
        raise SystemExit(f"Build directory not found: {build_dir}")
    if not placement_dir.is_dir():
        raise SystemExit(f"Placement directory not found: {placement_dir}")
    if args.run and not executable.is_file():
        raise SystemExit(f"Executable not found: {executable}")

    bn_wt, zns_wt = ratio_parts(args.ratio_tag)
    placements = resolve_placements(placement_dir, args.placements)

    if args.shuffle:
        rng = random.Random(args.seed)
        placements = list(placements)
        rng.shuffle(placements)
    else:
        placements = sorted(placements)

    placements = placements[args.start_index :]
    if args.max_placements > 0:
        placements = placements[: args.max_placements]
    if not placements:
        raise SystemExit("No placements selected.")

    macro_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    print(f"Project root: {project_root}")
    print(f"Build dir:    {build_dir}")
    print(f"Ratio tag:    {args.ratio_tag}")
    print(f"Placement dir:{placement_dir}")
    print(f"Placements:   {len(placements)}")
    print(f"Optical:      {args.optical_params if args.optical_params else '<code defaults/estimates>'}")
    print(f"Scatter:      {args.scatter_metric}")
    print(f"Target scat:  {args.target_primary_scatter}")
    print(f"Run enabled:  {'yes' if args.run else 'no'}")

    for idx, placement in enumerate(placements, start=1):
        placement_rel_to_build = "../" + rel_to_project(placement, project_root)
        placement_tag = placement_display_tag(placement, placement_dir)
        macro_path = macro_dir / f"{placement_tag}_StageD_OpticalHomogenization.mac"
        macro_rel_to_build = "../" + rel_to_project(macro_path, project_root)
        log_path = log_dir / f"p{idx:04d}_{placement_tag}.log"

        macro_path.write_text(
            macro_text(
                placement_rel_to_build=placement_rel_to_build,
                bn_wt=bn_wt,
                zns_wt=zns_wt,
                beam_on=args.beam_on,
                optical_params=args.optical_params,
                wavelength_nm=args.wavelength_nm,
                source_mode=args.source_mode,
                boundary_mode=args.boundary_mode,
                reentry_mode=args.reentry_mode,
                particle_reentry_mode=args.particle_reentry_mode,
                matrix_reentry_mode=args.matrix_reentry_mode,
                scatter_metric=args.scatter_metric,
                target_primary_scatter=args.target_primary_scatter,
                theta_threshold_deg=args.theta_threshold_deg,
                max_reentry=args.max_reentry,
                max_steps=args.max_steps,
                max_path_length_um=args.max_path_length_um,
                write_event_csv=args.write_event_csv,
                event_sampling_rate=args.event_sampling_rate,
                max_event_rows=args.max_event_rows,
                random_seed=args.seed + args.start_index + idx,
            ),
            encoding="utf-8",
        )

        print(f"[{idx}/{len(placements)}] {placement.relative_to(placement_dir)}")
        print(f"  macro: {macro_path}")
        print(f"  log:   {log_path}")
        print(f"  cmd:   cd {build_dir} && ./{args.executable_name} {macro_rel_to_build}")

        if args.run:
            with log_path.open("w", encoding="utf-8") as log_file:
                subprocess.run(
                    [str(executable), macro_rel_to_build],
                    cwd=str(build_dir),
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    check=True,
                )

    if not args.run:
        print("Macros generated. Re-run with --run to execute them.")


if __name__ == "__main__":
    main()
