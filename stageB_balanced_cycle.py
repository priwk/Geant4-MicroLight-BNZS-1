#!/usr/bin/env python3
import argparse
import csv
import hashlib
import math
import os
import re
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


DEFAULT_STAGEB_THICKNESSES = (
    "30,40,50,70,100,125,150,175,200,225,250,275,300,325,350,375,400,500,700,1000"
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Balanced Stage B replay across placement files."
    )
    parser.add_argument("ratios", nargs="*", help="Ratio folders such as 1-2 or 1-3.")
    parser.add_argument("--replay-multiplier", type=int, default=1)
    parser.add_argument("--seed", type=int, default=20260427)
    parser.add_argument("--min-thickness-um", type=float, default=30.0)
    parser.add_argument(
        "--thicknesses",
        default=None,
        help=(
            "Comma-separated thickness whitelist in um. When set, only matching "
            "*_neutron_capture_positions.csv files are processed."
        ),
    )
    parser.add_argument("--project-root", default=None)
    parser.add_argument("--build-dir", default=None)
    parser.add_argument("--executable", default=None)
    parser.add_argument("--macro", default=None)
    parser.add_argument("--keep-chunks", action="store_true")
    parser.add_argument("--keep-part-outputs", action="store_true")
    parser.add_argument("--merge-only", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--group-by-placement",
        action="store_true",
        help=(
            "Batch all selected thickness chunks for a placement into one Geant4 run "
            "instead of launching one run per thickness x placement chunk."
        ),
    )
    return parser.parse_args()


def natural_float_tag(path):
    name = Path(path).name
    key = "_neutron_capture_positions.csv"
    if not name.endswith(key):
        return math.inf
    try:
        return float(name[: -len(key)])
    except ValueError:
        return math.inf


def normalize_numeric_tag(value):
    return f"{float(value):.12g}"


def parse_thickness_whitelist(args_value, env_value):
    raw = args_value if args_value is not None else env_value
    if raw is None:
        return None

    values = set()
    for item in str(raw).split(","):
        token = item.strip()
        if not token:
            continue
        try:
            values.add(normalize_numeric_tag(token))
        except ValueError as exc:
            raise SystemExit(f"Invalid thickness value in whitelist: {token!r}") from exc

    return values if values else None


def stable_seed(base_seed, *parts):
    h = hashlib.sha256()
    h.update(str(base_seed).encode("utf-8"))
    for part in parts:
        h.update(b"\0")
        h.update(str(part).encode("utf-8"))
    return int.from_bytes(h.digest()[:8], "big")


def resolve_placement_root(project_root):
    pbc_root = project_root / "Input" / "output_pbc"
    if pbc_root.is_dir():
        return pbc_root
    return project_root / "Input" / "placements"


def is_main_placement_csv(path):
    name = Path(path).name
    return (
        Path(path).suffix.lower() == ".csv"
        and not name.endswith("_pbc_images.csv")
        and not name.endswith("_radius_stats.csv")
    )


def discover_placements(ratio_placement_dir):
    return sorted(
        path for path in ratio_placement_dir.rglob("*.csv")
        if is_main_placement_csv(path)
    )


def ensure_record_index_header(header):
    if "record_index" in header:
        return header, header.index("record_index")
    return header + ["record_index"], len(header)


def ensure_header_column(header, name):
    if name in header:
        return header, header.index(name)
    return header + [name], len(header)


def capture_file_uid(csv_path):
    hash_value = 14695981039346656037
    with open(csv_path, "rb") as source:
        while True:
            block = source.read(8192)
            if not block:
                break
            for byte in block:
                hash_value ^= byte
                hash_value = (hash_value * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return f"f{hash_value:016x}"


def read_capture_rows(csv_path, min_thickness_um):
    rows = []
    with open(csv_path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if header is None:
            return [], []
        header, record_index_col = ensure_record_index_header(list(header))
        header, input_uid_col = ensure_header_column(header, "input_file_uid")
        header, placement_replay_col = ensure_header_column(
            header, "placement_replay_index"
        )
        source_uid = capture_file_uid(csv_path)

        for idx, row in enumerate(reader):
            if not row:
                continue
            try:
                thickness = float(row[1])
            except (ValueError, IndexError):
                continue
            if thickness + 1.0e-12 < min_thickness_um:
                continue
            required_size = max(record_index_col, input_uid_col, placement_replay_col) + 1
            if len(row) < required_size:
                row = list(row) + [""] * (required_size - len(row))
            if not row[record_index_col].strip():
                row[record_index_col] = str(idx)
            if not row[input_uid_col].strip():
                row[input_uid_col] = source_uid
            if not row[placement_replay_col].strip():
                row[placement_replay_col] = "0"
            rows.append(row)
    return header, rows


def write_chunk(path, header, rows, placement_replay_index=0):
    path.parent.mkdir(parents=True, exist_ok=True)
    header, replay_col = ensure_header_column(list(header), "placement_replay_index")
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for source_row in rows:
            row = list(source_row)
            if len(row) <= replay_col:
                row.extend([""] * (replay_col + 1 - len(row)))
            row[replay_col] = str(placement_replay_index)
            writer.writerow(row)


def remove_empty_parent_dirs(path, stop_at):
    path = Path(path)
    stop_at = Path(stop_at).resolve()
    current = path
    while True:
        try:
            current_resolved = current.resolve()
        except FileNotFoundError:
            current_resolved = current

        if current_resolved == stop_at or current == current.parent:
            break

        try:
            current.rmdir()
        except OSError:
            break

        current = current.parent


def remove_path_and_empty_parents(path, stop_at):
    path = Path(path)
    if path.is_file():
        path.unlink(missing_ok=True)
        remove_empty_parent_dirs(path.parent, stop_at)


def remove_tree_if_exists(path):
    path = Path(path)
    try:
        shutil.rmtree(path)
    except FileNotFoundError:
        return


def cleanup_chunk_file(chunk_path, chunk_root, keep_chunks):
    if keep_chunks:
        return
    remove_path_and_empty_parents(chunk_path, chunk_root)


def cleanup_chunk_directory(chunk_dir, chunk_root, keep_chunks):
    if keep_chunks:
        return
    remove_tree_if_exists(chunk_dir)
    remove_empty_parent_dirs(Path(chunk_dir).parent, chunk_root)


def build_stageb_env(base_env, build_dir, project_root, bn_wt, zns_wt, ratio, placement_file):
    env = base_env.copy()
    env["BNZS_RUN_MODE"] = "StageB_ReplayAlphaLi"
    env.pop("BNZS_INPUT_CSV", None)
    env.pop("BNZS_INPUT_DIR", None)
    env["BNZS_BN_WT"] = bn_wt
    env["BNZS_ZNS_WT"] = zns_wt
    env["BNZS_OUTPUT_RATIO"] = ratio
    env["BNZS_PLACEMENT_FILE"] = relative_to_dir(placement_file, project_root)
    if "BNZS_STAGEB_OUTPUT_MODE" in base_env:
        env["BNZS_STAGEB_OUTPUT_MODE"] = base_env["BNZS_STAGEB_OUTPUT_MODE"]
    if "BNZS_ALPHALI_REPLAY_PER_CAPTURE" in base_env:
        env["BNZS_ALPHALI_REPLAY_PER_CAPTURE"] = base_env["BNZS_ALPHALI_REPLAY_PER_CAPTURE"]
    return env


def replay_per_capture_from_env(env):
    return int(env.get("BNZS_ALPHALI_REPLAY_PER_CAPTURE", "1"))


def collect_ratio_capture_assignments(
    ratio,
    capture_files,
    placement_files,
    min_thickness_um,
    replay_multiplier,
    seed,
):
    assignments = []

    for capture_file in capture_files:
        header, rows = read_capture_rows(capture_file, min_thickness_um)
        if not rows:
            assignments.append(
                {
                    "capture_file": capture_file,
                    "source_tag": capture_file.name.replace(
                        "_neutron_capture_positions.csv", ""
                    ),
                    "header": header,
                    "rows": [],
                    "replays": [],
                }
            )
            continue

        source_tag = capture_file.name.replace("_neutron_capture_positions.csv", "")
        shuffled = list(rows)
        rng = random.Random(stable_seed(seed, ratio, capture_file.name))
        rng.shuffle(shuffled)

        placements = list(placement_files)
        rng.shuffle(placements)

        replay_assignments = []
        for replay_idx in range(replay_multiplier):
            chunks = [[] for _ in placements]
            for idx, row in enumerate(shuffled):
                # A replay round shifts each record to a different placement
                # while preserving near-equal chunk sizes.
                chunks[(idx + replay_idx) % len(placements)].append(row)

            replay_assignments.append(
                {
                    "replay_idx": replay_idx,
                    "placements": placements,
                    "chunks": chunks,
                }
            )

        assignments.append(
            {
                "capture_file": capture_file,
                "source_tag": source_tag,
                "header": header,
                "rows": rows,
                "replays": replay_assignments,
            }
        )

    return assignments


def merge_ratio_outputs(project_root, ratio, source_tags, keep_part_outputs):
    output_ratio_dir = project_root / "Output" / "stageB" / ratio

    for source_tag in source_tags:
        merged_results = [
            merge_step_outputs(
                output_ratio_dir,
                source_tag,
                remove_parts=not keep_part_outputs,
            ),
            merge_capture_anchor_outputs(
                output_ratio_dir,
                source_tag,
                remove_parts=not keep_part_outputs,
            ),
            merge_zns_track_outputs(
                output_ratio_dir,
                source_tag,
                remove_parts=not keep_part_outputs,
            ),
            merge_unexpected_boundary_exit_outputs(
                output_ratio_dir,
                source_tag,
                remove_parts=not keep_part_outputs,
            ),
            merge_boundary_summary_outputs(
                output_ratio_dir,
                source_tag,
                remove_parts=not keep_part_outputs,
            ),
        ]
        for rows_merged, merged_path in merged_results:
            if merged_path is None:
                continue
            print(
                f">>> Merged {ratio} {source_tag}: "
                f"rows={rows_merged} -> {merged_path.name}"
            )


def ratio_parts(ratio):
    if "-" not in ratio:
        raise ValueError(f"Invalid ratio folder name: {ratio}")
    lhs, rhs = ratio.split("-", 1)
    float(lhs)
    float(rhs)
    return lhs, rhs


def relative_to_dir(path, start_dir):
    path = Path(path)
    start_dir = Path(start_dir).resolve()
    candidate = path.resolve() if path.is_absolute() else (start_dir / path).resolve()
    try:
        return candidate.relative_to(start_dir).as_posix()
    except ValueError:
        return os.path.relpath(candidate, start_dir)


def resolve_capture_input_root(project_root):
    preferred = project_root / "Input" / "stageA"
    legacy = project_root / "Input" / "neutron_capture_positions"
    if preferred.is_dir():
        return preferred, True
    return legacy, False


def ratio_capture_dir(input_root, ratio, stagea_layout):
    if stagea_layout:
        return input_root / ratio / "neutron_capture_positions"
    return input_root / ratio


def placement_tag_for_ratio(path, ratio_placement_dir):
    rel = Path(path).resolve().relative_to(Path(ratio_placement_dir).resolve())
    parts = list(rel.parts)
    parts[-1] = rel.stem
    return "__".join(parts)


def run_one(build_dir, executable, macro, env, log_file, dry_run):
    log_file.parent.mkdir(parents=True, exist_ok=True)
    cmd = [str(executable), str(macro)]

    if dry_run:
        print("DRY RUN:", " ".join(cmd))
        print("  BNZS_INPUT_CSV=", env.get("BNZS_INPUT_CSV"))
        print("  BNZS_INPUT_DIR=", env.get("BNZS_INPUT_DIR"))
        print("  BNZS_PLACEMENT_FILE=", env.get("BNZS_PLACEMENT_FILE"))
        print("  BNZS_STAGEB_OUTPUT_MODE=", env.get("BNZS_STAGEB_OUTPUT_MODE"))
        print("  BNZS_ALPHALI_REPLAY_PER_CAPTURE=", env.get("BNZS_ALPHALI_REPLAY_PER_CAPTURE"))
        print("  log=", log_file)
        return 0

    with open(log_file, "w") as log:
        return subprocess.call(cmd, cwd=build_dir, env=env, stdout=log, stderr=subprocess.STDOUT)


def make_beamon_macro(build_dir, base_macro, beam_on_count):
    with open(base_macro, encoding="utf-8") as src:
        lines = src.readlines()

    filtered = [
        line for line in lines
        if not line.lstrip().startswith("/run/beamOn")
    ]
    filtered.append(f"/run/beamOn {beam_on_count}\n")

    fd, path = tempfile.mkstemp(
        prefix="stageb_beamon_",
        suffix=".mac",
        dir=build_dir,
        text=True,
    )
    os.close(fd)
    with open(path, "w", encoding="utf-8") as out:
        out.writelines(filtered)
    return Path(path)


def merge_step_outputs(output_ratio_dir, source_tag, remove_parts):
    part_files = sorted(output_ratio_dir.glob(f"{source_tag}_m*_p*_alpha_li_steps.csv"))
    if not part_files:
        return 0, None

    merged_path = output_ratio_dir / f"{source_tag}_alpha_li_steps.csv"
    tmp_path = output_ratio_dir / f".{source_tag}_alpha_li_steps.tmp"

    total_rows = 0
    header_written = False

    with open(tmp_path, "w", newline="") as out:
        writer = None

        for part_file in part_files:
            with open(part_file, newline="") as f:
                reader = csv.reader(f)
                header = next(reader, None)
                if header is None:
                    continue

                if not header_written:
                    writer = csv.writer(out)
                    writer.writerow(header)
                    header_written = True

                for row in reader:
                    if not row:
                        continue
                    writer.writerow(row)
                    total_rows += 1

    if not header_written:
        tmp_path.unlink(missing_ok=True)
        return 0, None

    tmp_path.replace(merged_path)

    if remove_parts:
        for part_file in part_files:
            part_file.unlink(missing_ok=True)

    return total_rows, merged_path


def merge_capture_anchor_outputs(output_ratio_dir, source_tag, remove_parts):
    part_files = sorted(output_ratio_dir.glob(f"{source_tag}_m*_p*_capture_anchors.csv"))
    if not part_files:
        return 0, None

    merged_path = output_ratio_dir / f"{source_tag}_capture_anchors.csv"
    tmp_path = output_ratio_dir / f".{source_tag}_capture_anchors.tmp"

    total_rows = 0
    header_written = False
    with open(tmp_path, "w", newline="") as out:
        writer = None
        for part_file in part_files:
            with open(part_file, newline="") as f:
                reader = csv.reader(f)
                header = next(reader, None)
                if header is None:
                    continue
                if not header_written:
                    writer = csv.writer(out)
                    writer.writerow(header)
                    header_written = True
                for row in reader:
                    if not row:
                        continue
                    writer.writerow(row)
                    total_rows += 1

    if not header_written:
        tmp_path.unlink(missing_ok=True)
        return 0, None

    tmp_path.replace(merged_path)
    if remove_parts:
        for part_file in part_files:
            part_file.unlink(missing_ok=True)
    return total_rows, merged_path


def merge_zns_track_outputs(output_ratio_dir, source_tag, remove_parts):
    part_files = sorted(output_ratio_dir.glob(f"{source_tag}_m*_p*_zns_track_steps.csv"))
    if not part_files:
        return 0, None

    merged_path = output_ratio_dir / f"{source_tag}_zns_track_steps.csv"
    tmp_path = output_ratio_dir / f".{source_tag}_zns_track_steps.tmp"

    total_rows = 0
    header_written = False
    with open(tmp_path, "w", newline="") as out:
        writer = None
        for part_file in part_files:
            with open(part_file, newline="") as f:
                reader = csv.reader(f)
                header = next(reader, None)
                if header is None:
                    continue
                if not header_written:
                    writer = csv.writer(out)
                    writer.writerow(header)
                    header_written = True
                for row in reader:
                    if not row:
                        continue
                    writer.writerow(row)
                    total_rows += 1

    if not header_written:
        tmp_path.unlink(missing_ok=True)
        return 0, None

    tmp_path.replace(merged_path)
    if remove_parts:
        for part_file in part_files:
            part_file.unlink(missing_ok=True)
    return total_rows, merged_path


def merge_unexpected_boundary_exit_outputs(output_ratio_dir, source_tag, remove_parts):
    part_files = sorted(
        output_ratio_dir.glob(f"{source_tag}_m*_p*_unexpected_boundary_exits.csv")
    )
    if not part_files:
        return 0, None

    merged_path = output_ratio_dir / f"{source_tag}_unexpected_boundary_exits.csv"
    tmp_path = output_ratio_dir / f".{source_tag}_unexpected_boundary_exits.tmp"

    total_rows = 0
    header_written = False
    with open(tmp_path, "w", newline="") as out:
        writer = None
        for part_file in part_files:
            with open(part_file, newline="") as f:
                reader = csv.reader(f)
                header = next(reader, None)
                if header is None:
                    continue
                if not header_written:
                    writer = csv.writer(out)
                    writer.writerow(header)
                    header_written = True
                for row in reader:
                    if not row:
                        continue
                    writer.writerow(row)
                    total_rows += 1

    if not header_written:
        tmp_path.unlink(missing_ok=True)
        return 0, None

    tmp_path.replace(merged_path)
    if remove_parts:
        for part_file in part_files:
            part_file.unlink(missing_ok=True)
    return total_rows, merged_path


def merge_boundary_summary_outputs(output_ratio_dir, source_tag, remove_parts):
    part_files = sorted(output_ratio_dir.glob(f"{source_tag}_m*_p*_boundary_stop_summary.csv"))
    if not part_files:
        return 0, None

    merged_path = output_ratio_dir / f"{source_tag}_boundary_stop_summary.csv"
    fieldnames = [
        "thickness_um",
        "placement_file",
        "n_physical_surface_exit",
        "sum_physical_surface_exit_ekin_post_keV",
        "n_unexpected_artificial_exit",
        "sum_unexpected_artificial_exit_ekin_post_keV",
        "max_unexpected_artificial_exit_ekin_post_keV",
        "n_unexpected_artificial_exit_alpha",
        "sum_unexpected_artificial_exit_alpha_ekin_post_keV",
        "n_unexpected_artificial_exit_Li7",
        "sum_unexpected_artificial_exit_Li7_ekin_post_keV",
        "n_unexpected_bulk_exit",
        "sum_unexpected_bulk_exit_ekin_post_keV",
    ]
    int_fields = {
        "n_physical_surface_exit",
        "n_unexpected_artificial_exit",
        "n_unexpected_artificial_exit_alpha",
        "n_unexpected_artificial_exit_Li7",
        "n_unexpected_bulk_exit",
    }
    max_field = "max_unexpected_artificial_exit_ekin_post_keV"

    merged = {}
    for part_file in part_files:
        with open(part_file, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                key = (
                    normalize_numeric_tag(row.get("thickness_um", "0")),
                    row.get("placement_file", "").strip(),
                )
                acc = merged.setdefault(
                    key,
                    {
                        "thickness_um": key[0],
                        "placement_file": key[1],
                        **{name: 0 for name in int_fields},
                        **{
                            name: 0.0
                            for name in fieldnames
                            if name not in {"thickness_um", "placement_file"} | int_fields
                        },
                    },
                )
                for name in fieldnames:
                    if name in {"thickness_um", "placement_file"}:
                        continue
                    value = row.get(name, "").strip() or "0"
                    if name in int_fields:
                        acc[name] += int(float(value))
                    elif name == max_field:
                        acc[name] = max(acc[name], float(value))
                    else:
                        acc[name] += float(value)

    with open(merged_path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()
        for key in sorted(merged, key=lambda item: (float(item[0]), item[1])):
            record = merged[key]
            writer.writerow(record)

    if remove_parts:
        for part_file in part_files:
            part_file.unlink(missing_ok=True)

    return len(merged), merged_path


def merge_existing_outputs(project_root, ratios, remove_parts, dry_run=False):
    output_root = project_root / "Output" / "stageB"
    if not ratios:
        ratios = sorted(p.name for p in output_root.iterdir() if p.is_dir()) if output_root.exists() else []

    pattern = re.compile(
        r"^(?P<tag>.+)_m\d+_p\d+_(?:.+_)?"
        r"(?P<kind>alpha_li_steps|capture_anchors|zns_track_steps|boundary_stop_summary)\.csv$"
    )
    merged_any = False

    for ratio in ratios:
        ratio_dir = output_root / ratio
        if not ratio_dir.is_dir():
            print(f">>> Skip {ratio}: output directory not found")
            continue

        source_tags = sorted(
            {
                match.group("tag")
                for path in ratio_dir.glob("*.csv")
                for match in [pattern.match(path.name)]
                if match is not None
            }
        )

        for source_tag in source_tags:
            if dry_run:
                part_count = len(list(ratio_dir.glob(f"{source_tag}_m*_p*_*.csv")))
                print(
                    f"DRY RUN: would merge {part_count} files for "
                    f"{ratio}/{source_tag}"
                )
                continue

            merged_results = [
                merge_step_outputs(ratio_dir, source_tag, remove_parts),
                merge_capture_anchor_outputs(ratio_dir, source_tag, remove_parts),
                merge_zns_track_outputs(ratio_dir, source_tag, remove_parts),
                merge_unexpected_boundary_exit_outputs(ratio_dir, source_tag, remove_parts),
                merge_boundary_summary_outputs(ratio_dir, source_tag, remove_parts),
            ]
            for rows_merged, merged_path in merged_results:
                if merged_path is None:
                    continue
                merged_any = True
                print(
                    f">>> Merged {ratio} {source_tag}: "
                    f"rows={rows_merged} -> {merged_path.name}"
                )

    if not merged_any:
        print("No Stage B part outputs found to merge.")


def run_grouped_by_chunk(
    args,
    project_root,
    build_dir,
    executable,
    macro,
    chunk_root,
    log_root,
    ratio,
    bn_wt,
    zns_wt,
    assignments,
):
    total_runs = 0
    failed_runs = 0
    source_tags = []

    for assignment in assignments:
        source_tag = assignment["source_tag"]
        source_tags.append(source_tag)
        if not assignment["rows"]:
            print(f">>> Skip {ratio}/{assignment['capture_file'].name}: no compatible records")
            continue

        for replay in assignment["replays"]:
            replay_idx = replay["replay_idx"]
            placements = replay["placements"]
            chunks = replay["chunks"]

            for placement_idx, (placement_file, chunk_rows) in enumerate(zip(placements, chunks)):
                if not chunk_rows:
                    continue

                placement_tag = placement_tag_for_ratio(placement_file, resolve_placement_root(project_root) / ratio)
                chunk_name = (
                    f"{source_tag}_m{replay_idx + 1:02d}_"
                    f"p{placement_idx + 1:04d}_{placement_tag}_"
                    "neutron_capture_positions.csv"
                )
                chunk_path = (
                    chunk_root / ratio / source_tag / f"m{replay_idx + 1:02d}" / chunk_name
                )
                write_chunk(
                    chunk_path,
                    assignment["header"],
                    chunk_rows,
                    placement_replay_index=replay_idx,
                )

                log_file = (
                    log_root
                    / ratio
                    / source_tag
                    / f"m{replay_idx + 1:02d}"
                    / f"p{placement_idx + 1:04d}_{placement_tag}.log"
                )

                env = build_stageb_env(
                    os.environ,
                    build_dir,
                    project_root,
                    bn_wt,
                    zns_wt,
                    ratio,
                    placement_file,
                )
                env["BNZS_INPUT_CSV"] = relative_to_dir(chunk_path, build_dir)

                replay_per_capture = replay_per_capture_from_env(env)
                beam_on_count = len(chunk_rows) * replay_per_capture
                macro_to_run = make_beamon_macro(build_dir, macro, beam_on_count)

                total_runs += 1
                print(
                    f">>> {ratio} {source_tag} m{replay_idx + 1:02d} "
                    f"placement {placement_idx + 1}/{len(placements)} "
                    f"events={len(chunk_rows)} beamOn={beam_on_count} placement={placement_tag}"
                )

                try:
                    code = run_one(build_dir, executable, macro_to_run, env, log_file, args.dry_run)
                finally:
                    macro_to_run.unlink(missing_ok=True)
                if code != 0:
                    failed_runs += 1
                    print(f"!!! Failed code={code}: {log_file}")
                    raise SystemExit(code)
                cleanup_chunk_file(chunk_path, chunk_root, args.keep_chunks)

    merge_ratio_outputs(project_root, ratio, source_tags, args.keep_part_outputs)
    return total_runs, failed_runs


def run_grouped_by_placement(
    args,
    project_root,
    build_dir,
    executable,
    macro,
    chunk_root,
    log_root,
    ratio,
    bn_wt,
    zns_wt,
    placement_files,
    assignments,
):
    total_runs = 0
    failed_runs = 0
    source_tags = [assignment["source_tag"] for assignment in assignments]
    placement_index_lookup = {path.resolve(): idx for idx, path in enumerate(placement_files)}
    ratio_placement_dir = resolve_placement_root(project_root) / ratio

    for assignment in assignments:
        if not assignment["rows"]:
            print(f">>> Skip {ratio}/{assignment['capture_file'].name}: no compatible records")

    for placement_file in placement_files:
        placement_key = placement_file.resolve()
        placement_tag = placement_tag_for_ratio(placement_file, ratio_placement_dir)
        placement_pos = placement_index_lookup[placement_key] + 1
        placement_chunk_dir = chunk_root / ratio / f"placement_{placement_pos:04d}_{placement_tag}"
        remove_tree_if_exists(placement_chunk_dir)

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
                write_chunk(
                    chunk_path,
                    assignment["header"],
                    placement_chunk_rows,
                    placement_replay_index=replay_idx,
                )
                event_count += len(placement_chunk_rows)
                chunk_count += 1

        if event_count <= 0:
            cleanup_chunk_directory(placement_chunk_dir, chunk_root, keep_chunks=False)
            continue

        env = build_stageb_env(
            os.environ,
            build_dir,
            project_root,
            bn_wt,
            zns_wt,
            ratio,
            placement_file,
        )
        env["BNZS_INPUT_DIR"] = relative_to_dir(
            placement_chunk_dir,
            build_dir,
        )

        replay_per_capture = replay_per_capture_from_env(env)
        beam_on_count = event_count * replay_per_capture
        macro_to_run = make_beamon_macro(build_dir, macro, beam_on_count)
        log_file = log_root / ratio / "grouped" / f"p{placement_pos:04d}_{placement_tag}.log"

        total_runs += 1
        print(
            f">>> {ratio} grouped placement {placement_pos}/{len(placement_files)} "
            f"events={event_count} beamOn={beam_on_count} placement={placement_tag} "
            f"chunks={chunk_count}"
        )

        try:
            code = run_one(build_dir, executable, macro_to_run, env, log_file, args.dry_run)
        finally:
            macro_to_run.unlink(missing_ok=True)
        if code != 0:
            failed_runs += 1
            print(f"!!! Failed code={code}: {log_file}")
            raise SystemExit(code)
        cleanup_chunk_directory(placement_chunk_dir, chunk_root, args.keep_chunks)

    merge_ratio_outputs(project_root, ratio, source_tags, args.keep_part_outputs)
    return total_runs, failed_runs


def main():
    args = parse_args()
    if args.replay_multiplier <= 0:
        raise SystemExit("--replay-multiplier must be > 0")
    thickness_whitelist = parse_thickness_whitelist(
        args.thicknesses,
        os.environ.get("STAGEB_THICKNESSES"),
    )

    project_root = Path(args.project_root).resolve() if args.project_root else Path(__file__).resolve().parent
    build_dir = Path(args.build_dir).resolve() if args.build_dir else project_root / "build"
    executable = (
        Path(args.executable).resolve()
        if args.executable
        else build_dir / "Geant4-MicroLight-BNZS"
    )
    macro = Path(args.macro).resolve() if args.macro else project_root / "run.mac"

    input_root, stagea_layout = resolve_capture_input_root(project_root)
    placement_root = resolve_placement_root(project_root)
    chunk_root = build_dir / "stageB_balanced_chunks"
    log_root = project_root / "logs" / "stageB" / "balanced"

    ratios = args.ratios
    if not ratios:
        ratios = sorted(
            p.name
            for p in input_root.iterdir()
            if p.is_dir() and "-" in p.name and ratio_capture_dir(input_root, p.name, stagea_layout).is_dir()
        )

    if not ratios:
        raise SystemExit(f"No ratio directories found under {input_root}")

    if args.merge_only:
        merge_existing_outputs(
            project_root,
            ratios,
            remove_parts=not args.keep_part_outputs,
            dry_run=args.dry_run,
        )
        return

    if not args.dry_run and not executable.exists():
        raise SystemExit(f"Executable not found: {executable}")

    if not args.keep_chunks and chunk_root.exists():
        remove_tree_if_exists(chunk_root)

    total_runs = 0
    failed_runs = 0

    for ratio in ratios:
        bn_wt, zns_wt = ratio_parts(ratio)
        ratio_input_dir = ratio_capture_dir(input_root, ratio, stagea_layout)
        ratio_placement_dir = placement_root / ratio

        capture_files = sorted(
            (
                path
                for path in ratio_input_dir.glob("*_neutron_capture_positions.csv")
                if thickness_whitelist is None
                or normalize_numeric_tag(natural_float_tag(path)) in thickness_whitelist
            ),
            key=lambda p: (natural_float_tag(p), p.name),
        )
        placement_files = discover_placements(ratio_placement_dir)

        if not capture_files:
            print(f">>> Skip {ratio}: no capture CSV files in {ratio_input_dir}")
            continue
        if not placement_files:
            print(f">>> Skip {ratio}: no placement CSV files in {ratio_placement_dir}")
            continue

        print()
        print(f"=== Ratio {ratio}: {len(capture_files)} capture files, {len(placement_files)} placements ===")
        assignments = collect_ratio_capture_assignments(
            ratio,
            capture_files,
            placement_files,
            args.min_thickness_um,
            args.replay_multiplier,
            args.seed,
        )

        if args.group_by_placement:
            runs, failures = run_grouped_by_placement(
                args,
                project_root,
                build_dir,
                executable,
                macro,
                chunk_root,
                log_root,
                ratio,
                bn_wt,
                zns_wt,
                placement_files,
                assignments,
            )
        else:
            runs, failures = run_grouped_by_chunk(
                args,
                project_root,
                build_dir,
                executable,
                macro,
                chunk_root,
                log_root,
                ratio,
                bn_wt,
                zns_wt,
                assignments,
            )

        total_runs += runs
        failed_runs += failures

    print()
    print(f"=== Stage B balanced cycling complete: runs={total_runs}, failed={failed_runs} ===")
    print(f"Logs:   {log_root}")
    print(f"Output: {project_root / 'Output' / 'stageB'}")


if __name__ == "__main__":
    main()
