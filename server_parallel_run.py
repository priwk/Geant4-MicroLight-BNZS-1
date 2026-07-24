#!/usr/bin/env python3
import argparse
import os
import random
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import stageB_balanced_cycle as stageb
import stageD_run_batch as staged


@dataclass
class StageDTask:
    index: int
    total: int
    placement: Path
    placement_tag: str
    macro_path: Path
    macro_rel_to_build: str
    log_path: Path


@dataclass
class StageBTask:
    ratio: str
    placement_pos: int
    placement_total: int
    placement_tag: str
    event_count: int
    chunk_count: int
    macro_path: Path
    chunk_dir: Path
    log_path: Path
    env: dict


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Server-oriented parallel runner for Stage B and Stage D. "
            "This script parallelizes by launching multiple independent "
            "Geant4 processes, because the current executable uses "
            "single-threaded G4RunManager."
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    stage_d = subparsers.add_parser(
        "staged",
        help="Parallel Stage D runs across placements.",
    )
    stage_d.add_argument("--ratio-tag", required=True)
    stage_d.add_argument("--project-root", default=None)
    stage_d.add_argument("--build-dir", default=None)
    stage_d.add_argument("--placement-dir", default=None)
    stage_d.add_argument("--placements", nargs="*", default=None)
    stage_d.add_argument("--max-placements", type=int, default=0)
    stage_d.add_argument("--start-index", type=int, default=0)
    stage_d.add_argument("--shuffle", action="store_true")
    stage_d.add_argument("--seed", type=int, default=20260507)
    stage_d.add_argument("--beam-on", type=int, default=10000)
    stage_d.add_argument("--wavelength-nm", type=float, default=450.0)
    stage_d.add_argument("--optical-params", default=None)
    stage_d.add_argument(
        "--source-mode",
        default="uniform_all_phase",
        choices=["uniform_ZnS", "uniform_all_phase"],
    )
    stage_d.add_argument(
        "--boundary-mode",
        default="periodic_wrap",
        choices=["periodic_wrap", "escape", "same_phase_reentry"],
    )
    stage_d.add_argument(
        "--reentry-mode",
        default="state_matched",
        choices=["state_matched", "same_phase_rho_over_R", "same_phase_random"],
    )
    stage_d.add_argument(
        "--particle-reentry-mode",
        default="sphere_q_mu",
        choices=["sphere_q_mu", "same_phase_rho_over_R", "same_phase_random"],
    )
    stage_d.add_argument(
        "--matrix-reentry-mode",
        default="clearance_binned_portal",
        choices=[
            "clearance_binned_portal",
            "random_matrix_debug",
            "distance_matched_matrix",
            "random_matrix",
        ],
    )
    stage_d.add_argument(
        "--scatter-metric",
        default="angle_threshold",
        choices=[
            "angle_threshold",
            "particle_encounter_angle_threshold",
            "particle_encounter_no_threshold",
        ],
    )
    stage_d.add_argument("--target-primary-scatter", type=int, default=0)
    stage_d.add_argument("--theta-threshold-deg", type=float, default=0.5)
    stage_d.add_argument("--max-reentry", type=int, default=10000)
    stage_d.add_argument("--max-steps", type=int, default=100000)
    stage_d.add_argument("--max-path-length-um", type=float, default=5000.0)
    stage_d.add_argument("--macro-dir", default=None)
    stage_d.add_argument("--log-dir", default=None)
    stage_d.add_argument("--executable-name", default="Geant4-MicroLight-BNZS")
    stage_d.add_argument("--jobs", type=int, default=0)
    stage_d.add_argument("--merge-after", action="store_true")
    stage_d.add_argument("--python-exe", default=sys.executable)
    stage_d.add_argument("--dry-run", action="store_true")
    stage_d.add_argument(
        "--no-pin-libs",
        action="store_true",
        help="Do not force OMP/OpenBLAS/MKL-style helper libraries to one thread.",
    )

    stage_b = subparsers.add_parser(
        "stageb",
        help="Parallel Stage B grouped-by-placement runs.",
    )
    stage_b.add_argument("ratios", nargs="*")
    stage_b.add_argument("--replay-multiplier", type=int, default=1)
    stage_b.add_argument("--seed", type=int, default=20260427)
    stage_b.add_argument("--min-thickness-um", type=float, default=30.0)
    stage_b.add_argument("--thicknesses", default=None)
    stage_b.add_argument("--project-root", default=None)
    stage_b.add_argument("--build-dir", default=None)
    stage_b.add_argument("--executable", default=None)
    stage_b.add_argument("--macro", default=None)
    stage_b.add_argument("--keep-chunks", action="store_true")
    stage_b.add_argument("--keep-part-outputs", action="store_true")
    stage_b.add_argument("--merge-only", action="store_true")
    stage_b.add_argument("--dry-run", action="store_true")
    stage_b.add_argument("--jobs", type=int, default=0)
    stage_b.add_argument(
        "--no-pin-libs",
        action="store_true",
        help="Do not force OMP/OpenBLAS/MKL-style helper libraries to one thread.",
    )

    return parser.parse_args()


def project_root_from_args(raw_value):
    return Path(raw_value).resolve() if raw_value else Path(__file__).resolve().parent


def normalize_jobs(requested_jobs, task_count):
    if task_count <= 0:
        return 1
    if requested_jobs is None or requested_jobs <= 0:
        requested_jobs = os.cpu_count() or 1
    return max(1, min(requested_jobs, task_count))


def pin_helper_libraries_to_one_thread(env):
    for key in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "BLIS_NUM_THREADS",
    ):
        env[key] = "1"


def server_env(pin_libs):
    env = os.environ.copy()
    if pin_libs:
        pin_helper_libraries_to_one_thread(env)
    return env


def execute_process(executable, build_dir, macro_arg, log_path, env, dry_run):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [str(executable), str(macro_arg)]

    if dry_run:
        print("DRY RUN:", " ".join(cmd))
        print("  cwd=", build_dir)
        print("  log=", log_path)
        return 0

    with log_path.open("w", encoding="utf-8") as log_file:
        completed = subprocess.run(
            cmd,
            cwd=str(build_dir),
            env=env,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            check=False,
        )
    return completed.returncode


def validate_stage_d_args(args):
    if args.max_placements < 0 or args.start_index < 0:
        raise SystemExit("--max-placements and --start-index must be >= 0")
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


def resolve_stage_d_paths(args):
    project_root = project_root_from_args(args.project_root)
    build_dir = Path(args.build_dir).resolve() if args.build_dir else project_root / "build"
    placement_dir = (
        Path(args.placement_dir).resolve()
        if args.placement_dir
        else project_root / "Input" / "output_pbc" / args.ratio_tag
    )
    wavelength_tag = (
        f"lambda_{int(round(args.wavelength_nm))}nm"
        if abs(args.wavelength_nm - round(args.wavelength_nm)) < 1.0e-9
        else f"lambda_{args.wavelength_nm:.3f}nm"
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
    return project_root, build_dir, placement_dir, macro_dir, log_dir, executable


def select_stage_d_placements(args, placement_dir):
    placements = staged.resolve_placements(placement_dir, args.placements)
    if args.shuffle:
        placements = list(placements)
        random.Random(args.seed).shuffle(placements)
    else:
        placements = sorted(placements)

    placements = placements[args.start_index :]
    if args.max_placements > 0:
        placements = placements[: args.max_placements]
    if not placements:
        raise SystemExit("No placements selected.")
    return placements


def build_stage_d_tasks(args, project_root, build_dir, placement_dir, macro_dir, log_dir):
    bn_wt, zns_wt = staged.ratio_parts(args.ratio_tag)
    placements = select_stage_d_placements(args, placement_dir)

    macro_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    tasks = []
    for idx, placement in enumerate(placements, start=1):
        placement_rel_to_build = "../" + staged.rel_to_project(placement, project_root)
        placement_tag = staged.placement_display_tag(placement, placement_dir)
        macro_path = macro_dir / f"{placement_tag}_StageD_OpticalHomogenization.mac"
        macro_rel_to_build = "../" + staged.rel_to_project(macro_path, project_root)
        log_path = log_dir / f"p{idx:04d}_{placement_tag}.log"

        macro_path.write_text(
            staged.macro_text(
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
            ),
            encoding="utf-8",
        )

        tasks.append(
            StageDTask(
                index=idx,
                total=len(placements),
                placement=placement,
                placement_tag=placement_tag,
                macro_path=macro_path,
                macro_rel_to_build=macro_rel_to_build,
                log_path=log_path,
            )
        )

    return tasks


def run_stage_d(args):
    validate_stage_d_args(args)
    project_root, build_dir, placement_dir, macro_dir, log_dir, executable = resolve_stage_d_paths(args)

    if not build_dir.is_dir():
        raise SystemExit(f"Build directory not found: {build_dir}")
    if not placement_dir.is_dir():
        raise SystemExit(f"Placement directory not found: {placement_dir}")
    if not args.dry_run and not executable.is_file():
        raise SystemExit(f"Executable not found: {executable}")

    tasks = build_stage_d_tasks(args, project_root, build_dir, placement_dir, macro_dir, log_dir)
    jobs = normalize_jobs(args.jobs, len(tasks))
    env = server_env(pin_libs=not args.no_pin_libs)

    print(f"Stage D placements: {len(tasks)}")
    print(f"Stage D jobs:       {jobs}")
    print(f"Build dir:          {build_dir}")
    print(f"Logs:               {log_dir}")
    print(f"Macros:             {macro_dir}")

    failures = []
    if args.dry_run:
        for task in tasks:
            print(f"[{task.index}/{task.total}] {task.placement.relative_to(placement_dir)}")
            execute_process(executable, build_dir, task.macro_rel_to_build, task.log_path, env, True)
    else:
        with ThreadPoolExecutor(max_workers=jobs) as executor:
            future_to_task = {
                executor.submit(
                    execute_process,
                    executable,
                    build_dir,
                    task.macro_rel_to_build,
                    task.log_path,
                    env,
                    False,
                ): task
                for task in tasks
            }

            done_count = 0
            for future in as_completed(future_to_task):
                task = future_to_task[future]
                code = future.result()
                done_count += 1
                if code == 0:
                    print(
                        f"[done {done_count}/{len(tasks)}] "
                        f"{task.placement_tag} -> {task.log_path}"
                    )
                else:
                    failures.append((task, code))
                    print(
                        f"[failed {done_count}/{len(tasks)}] "
                        f"{task.placement_tag} code={code} log={task.log_path}"
                    )

    if failures:
        print(f"Stage D failed runs: {len(failures)}")
        for task, code in failures:
            print(f"  code={code} placement={task.placement_tag} log={task.log_path}")
        raise SystemExit(1)

    if args.merge_after and not args.dry_run:
        merge_script = project_root / "scripts" / "merge_stageD_optical_params.py"
        subprocess.run(
            [
                args.python_exe,
                str(merge_script),
                "--ratio",
                args.ratio_tag,
                "--project-root",
                str(project_root),
            ],
            check=True,
        )
        print(f"Merged Stage D optical params for ratio {args.ratio_tag}")


def resolve_stage_b_paths(args):
    project_root = project_root_from_args(args.project_root)
    build_dir = Path(args.build_dir).resolve() if args.build_dir else project_root / "build"
    executable = (
        Path(args.executable).resolve()
        if args.executable
        else build_dir / "Geant4-MicroLight-BNZS"
    )
    macro = Path(args.macro).resolve() if args.macro else project_root / "run.mac"
    input_root, stagea_layout = stageb.resolve_capture_input_root(project_root)
    placement_root = stageb.resolve_placement_root(project_root)
    chunk_root = build_dir / "stageB_balanced_chunks_parallel"
    log_root = project_root / "logs" / "stageB" / "balanced_parallel"
    return project_root, build_dir, executable, macro, input_root, stagea_layout, placement_root, chunk_root, log_root


def resolve_stage_b_ratios(args, input_root, stagea_layout):
    ratios = list(args.ratios)
    if ratios:
        return ratios

    return sorted(
        p.name
        for p in input_root.iterdir()
        if p.is_dir() and "-" in p.name and stageb.ratio_capture_dir(input_root, p.name, stagea_layout).is_dir()
    )


def prepare_stage_b_tasks_for_ratio(
    args,
    project_root,
    build_dir,
    macro,
    chunk_root,
    log_root,
    ratio,
    bn_wt,
    zns_wt,
    placement_files,
    assignments,
    base_env,
):
    tasks = []
    ratio_placement_dir = stageb.resolve_placement_root(project_root) / ratio
    placement_index_lookup = {path.resolve(): idx for idx, path in enumerate(placement_files)}

    for placement_file in placement_files:
        placement_key = placement_file.resolve()
        placement_tag = stageb.placement_tag_for_ratio(placement_file, ratio_placement_dir)
        placement_pos = placement_index_lookup[placement_key] + 1
        placement_chunk_dir = chunk_root / ratio / f"placement_{placement_pos:04d}_{placement_tag}"
        stageb.remove_tree_if_exists(placement_chunk_dir)

        event_count = 0
        chunk_count = 0

        for assignment in assignments:
            if not assignment["rows"]:
                continue

            source_tag = assignment["source_tag"]
            for replay in assignment["replays"]:
                replay_idx = replay["replay_idx"]
                placement_chunk_rows = None

                for replay_placement_file, chunk_rows in zip(replay["placements"], replay["chunks"]):
                    if replay_placement_file.resolve() == placement_key:
                        placement_chunk_rows = chunk_rows
                        break

                if not placement_chunk_rows:
                    continue

                chunk_name = (
                    f"{source_tag}_m{replay_idx + 1:02d}_"
                    f"p{placement_pos:04d}_{placement_tag}_"
                    "neutron_capture_positions.csv"
                )
                chunk_path = placement_chunk_dir / chunk_name
                stageb.write_chunk(
                    chunk_path,
                    assignment["header"],
                    placement_chunk_rows,
                    placement_replay_index=replay_idx,
                )
                event_count += len(placement_chunk_rows)
                chunk_count += 1

        if event_count <= 0:
            stageb.remove_tree_if_exists(placement_chunk_dir)
            continue

        env = stageb.build_stageb_env(
            base_env,
            build_dir,
            project_root,
            bn_wt,
            zns_wt,
            ratio,
            placement_file,
        )
        env["BNZS_INPUT_DIR"] = stageb.relative_to_dir(placement_chunk_dir, build_dir)

        replay_per_capture = stageb.replay_per_capture_from_env(env)
        beam_on_count = event_count * replay_per_capture
        macro_to_run = stageb.make_beamon_macro(build_dir, macro, beam_on_count)
        log_file = log_root / ratio / "grouped" / f"p{placement_pos:04d}_{placement_tag}.log"

        tasks.append(
            StageBTask(
                ratio=ratio,
                placement_pos=placement_pos,
                placement_total=len(placement_files),
                placement_tag=placement_tag,
                event_count=event_count,
                chunk_count=chunk_count,
                macro_path=macro_to_run,
                chunk_dir=placement_chunk_dir,
                log_path=log_file,
                env=env,
            )
        )

    return tasks


def cleanup_stage_b_task(task, keep_chunks):
    task.macro_path.unlink(missing_ok=True)
    if not keep_chunks:
        stageb.remove_tree_if_exists(task.chunk_dir)


def run_stage_b(args):
    if args.replay_multiplier <= 0:
        raise SystemExit("--replay-multiplier must be > 0")

    (
        project_root,
        build_dir,
        executable,
        macro,
        input_root,
        stagea_layout,
        placement_root,
        chunk_root,
        log_root,
    ) = resolve_stage_b_paths(args)

    ratios = resolve_stage_b_ratios(args, input_root, stagea_layout)
    if not ratios:
        raise SystemExit(f"No ratio directories found under {input_root}")

    if args.merge_only:
        stageb.merge_existing_outputs(
            project_root,
            ratios,
            remove_parts=not args.keep_part_outputs,
            dry_run=args.dry_run,
        )
        return

    thickness_whitelist = stageb.parse_thickness_whitelist(
        args.thicknesses,
        os.environ.get("STAGEB_THICKNESSES"),
    )

    if not args.dry_run and not executable.exists():
        raise SystemExit(f"Executable not found: {executable}")

    base_env = server_env(pin_libs=not args.no_pin_libs)

    if not args.keep_chunks and chunk_root.exists():
        stageb.remove_tree_if_exists(chunk_root)

    total_runs = 0
    total_failures = 0

    for ratio in ratios:
        bn_wt, zns_wt = stageb.ratio_parts(ratio)
        ratio_input_dir = stageb.ratio_capture_dir(input_root, ratio, stagea_layout)
        ratio_placement_dir = placement_root / ratio

        capture_files = sorted(
            (
                path
                for path in ratio_input_dir.glob("*_neutron_capture_positions.csv")
                if thickness_whitelist is None
                or stageb.normalize_numeric_tag(stageb.natural_float_tag(path)) in thickness_whitelist
            ),
            key=lambda p: (stageb.natural_float_tag(p), p.name),
        )
        placement_files = stageb.discover_placements(ratio_placement_dir)

        if not capture_files:
            print(f">>> Skip {ratio}: no capture CSV files in {ratio_input_dir}")
            continue
        if not placement_files:
            print(f">>> Skip {ratio}: no placement CSV files in {ratio_placement_dir}")
            continue

        print()
        print(f"=== Ratio {ratio}: {len(capture_files)} capture files, {len(placement_files)} placements ===")
        assignments = stageb.collect_ratio_capture_assignments(
            ratio,
            capture_files,
            placement_files,
            args.min_thickness_um,
            args.replay_multiplier,
            args.seed,
        )

        tasks = prepare_stage_b_tasks_for_ratio(
            args,
            project_root,
            build_dir,
            macro,
            chunk_root,
            log_root,
            ratio,
            bn_wt,
            zns_wt,
            placement_files,
            assignments,
            base_env,
        )
        if not tasks:
            print(f">>> Skip {ratio}: no grouped placement jobs after filtering")
            continue

        jobs = normalize_jobs(args.jobs, len(tasks))
        print(f"Stage B jobs for {ratio}: {jobs}")
        print(f"Logs:   {log_root / ratio / 'grouped'}")

        failures = []
        if args.dry_run:
            for task in tasks:
                print(
                    f">>> {ratio} grouped placement {task.placement_pos}/{task.placement_total} "
                    f"events={task.event_count} chunks={task.chunk_count} placement={task.placement_tag}"
                )
                execute_process(
                    executable,
                    build_dir,
                    task.macro_path,
                    task.log_path,
                    task.env,
                    True,
                )
        else:
            with ThreadPoolExecutor(max_workers=jobs) as executor:
                future_to_task = {
                    executor.submit(
                        execute_process,
                        executable,
                        build_dir,
                        task.macro_path,
                        task.log_path,
                        task.env,
                        False,
                    ): task
                    for task in tasks
                }

                done_count = 0
                for future in as_completed(future_to_task):
                    task = future_to_task[future]
                    code = future.result()
                    done_count += 1
                    total_runs += 1
                    if code == 0:
                        print(
                            f"[done {done_count}/{len(tasks)}] "
                            f"{ratio} {task.placement_tag} events={task.event_count}"
                        )
                        cleanup_stage_b_task(task, args.keep_chunks)
                    else:
                        total_failures += 1
                        failures.append((task, code))
                        print(
                            f"[failed {done_count}/{len(tasks)}] "
                            f"{ratio} {task.placement_tag} code={code} log={task.log_path}"
                        )

        if failures:
            print(f">>> Skip merge for {ratio}: {len(failures)} placement jobs failed")
            for task, code in failures:
                print(f"  code={code} placement={task.placement_tag} log={task.log_path}")
            continue

        if not args.dry_run:
            source_tags = [assignment["source_tag"] for assignment in assignments]
            stageb.merge_ratio_outputs(project_root, ratio, source_tags, args.keep_part_outputs)
            if not args.keep_chunks:
                stageb.remove_tree_if_exists(chunk_root / ratio)

    print()
    print(f"=== Stage B parallel grouped run complete: runs={total_runs}, failed={total_failures} ===")
    print(f"Logs:   {log_root}")
    print(f"Output: {project_root / 'Output' / 'stageB'}")

    if total_failures > 0:
        raise SystemExit(1)


def main():
    args = parse_args()
    if args.command == "staged":
        run_stage_d(args)
        return
    if args.command == "stageb":
        run_stage_b(args)
        return
    raise SystemExit(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    main()
