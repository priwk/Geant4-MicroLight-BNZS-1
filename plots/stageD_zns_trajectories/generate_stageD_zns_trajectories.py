#!/usr/bin/env python3
"""Generate Stage D-style alpha/Li trajectory figures from Stage B tracks."""

from __future__ import annotations

import argparse
import csv
import math
import os
import random
import re
from collections import OrderedDict, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

DEPENDENCY_ERROR = """
Failed to import the plotting dependencies.

Use a clean virtual environment from the project root:

    cd ~/g4work/B2
    python3 -m venv .venv
    . .venv/bin/activate
    python3 -m pip install --upgrade pip
    python3 -m pip install -r plots/stageD_zns_trajectories/requirements.txt
    python3 plots/stageD_zns_trajectories/generate_stageD_zns_trajectories.py
""".strip()

try:
    import matplotlib

    matplotlib.use("Agg")

    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import font_manager
    from matplotlib.collections import LineCollection
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
except ImportError as exc:
    raise SystemExit(f"{DEPENDENCY_ERROR}\n\nOriginal error: {exc}") from None
except AttributeError as exc:
    if "_ARRAY_API" in str(exc):
        raise SystemExit(f"{DEPENDENCY_ERROR}\n\nOriginal error: {exc}") from None
    raise


SCRIPT_DIR = Path(__file__).resolve().parent
FONT_DIR = SCRIPT_DIR / "fonts"

DESIRED_RATIO_ORDER = ["2-1", "1-1", "1-1.5", "1-2", "1-2.5", "1-3"]
RATIO_RE = re.compile(
    r"^(?P<bn>[+]?(?:\d+(?:\.\d*)?|\.\d+))-(?P<zns>[+]?(?:\d+(?:\.\d*)?|\.\d+))$"
)
SLIM_TRACK_RE = re.compile(
    r"^(?P<thickness>[+]?(?:\d+(?:\.\d*)?|\.\d+))_zns_track_steps\.csv$"
)
FULL_TRACK_RE = re.compile(
    r"^(?P<thickness>[+]?(?:\d+(?:\.\d*)?|\.\d+))_alpha_li_steps\.csv$"
)

PARTICLE_COLORS = {
    "alpha": "#c0392b",
    "Li7": "#0b3c79",
}
PARTICLE_LINEWIDTH = {
    "alpha": 1.45,
    "Li7": 1.05,
}
PHASE_COLORS = {
    "BN": "#0b3c79",
    "ZnS": "#c0392b",
    "binder_void": "#6b7280",
    "outside": "#111827",
    "world": "#111827",
    "unknown": "#8e44ad",
}
PHASE_ORDER = ["BN", "ZnS", "binder_void", "outside", "world", "unknown"]
PHASE_LABELS = {
    "BN": "BN phase",
    "ZnS": "ZnS phase",
    "binder_void": "Binder/porosity phase",
    "outside": "Outside domain",
    "world": "Outside domain",
    "unknown": "Unclassified phase",
}


@dataclass(frozen=True)
class RatioKey:
    tag: str
    bn_wt: float
    zns_wt: float

    @property
    def display_tag(self) -> str:
        return self.tag.replace("-", ":")


@dataclass
class TrackStep:
    particle: str
    source_event_uid: str
    physical_event_uid: str
    thickness_um: float
    phase_pre: str
    phase_post: str
    start: tuple[float, float, float]
    end: tuple[float, float, float]
    macro_start: tuple[float, float, float] | None
    macro_end: tuple[float, float, float] | None
    step_len_um: float
    edep_kev: float
    ekin_pre_kev: float
    ekin_post_kev: float
    trajectory_weight: float


@dataclass
class Trajectory:
    source_event_uid: str
    physical_event_uid: str
    particle: str
    ratio_tag: str
    thickness_um: float
    steps: list[TrackStep] = field(default_factory=list)
    total_edep_kev: float = 0.0
    total_len_um: float = 0.0
    trajectory_weight: float = 1.0
    input_mode: str = "unknown"

    @property
    def n_steps(self) -> int:
        return len(self.steps)

    @property
    def end_to_end_um(self) -> float:
        if not self.steps:
            return 0.0
        start = np.asarray(self.steps[0].start, dtype=float)
        end = np.asarray(self.steps[-1].end, dtype=float)
        return float(np.linalg.norm(end - start))


def configure_fonts() -> None:
    font_files = [
        FONT_DIR / "Arial.ttf",
        FONT_DIR / "Arial-Bold.ttf",
        FONT_DIR / "Arial-Italic.ttf",
        FONT_DIR / "Arial-BoldItalic.ttf",
    ]
    for path in font_files:
        if path.is_file():
            font_manager.fontManager.addfont(str(path))

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "legend.fontsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "axes.unicode_minus": False,
            "figure.dpi": 140,
        }
    )


configure_fonts()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot alpha/Li particle trajectories from Stage B CSVs. Full "
            "alpha_li_steps files are preferred; slim zns_track_steps files are "
            "used as a ZnS-only fallback."
        )
    )
    parser.add_argument(
        "--project-root",
        default=None,
        help="Project root. Defaults to the repository containing this script.",
    )
    parser.add_argument(
        "--stageb-root",
        default=None,
        help="Stage B output root. Defaults to <project-root>/Output/stageB.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory. Defaults to this plotting folder's output/.",
    )
    parser.add_argument(
        "--ratios",
        nargs="+",
        default=None,
        help="Optional ratio folders to process, such as 1-2 1-3. Defaults to all found ratios.",
    )
    parser.add_argument(
        "--thicknesses",
        nargs="+",
        type=float,
        default=None,
        help="Optional thicknesses in um to process, such as 1000. Defaults to all found thicknesses.",
    )
    parser.add_argument(
        "--showcase-ratio",
        default=None,
        help="Ratio used for trajectory showcase figures. Defaults to the first processed ratio.",
    )
    parser.add_argument(
        "--showcase-thickness",
        type=float,
        default=None,
        help="Thickness used for trajectory showcase. Defaults to the largest available thickness.",
    )
    parser.add_argument(
        "--sample-trajectories",
        type=int,
        default=80,
        help="Number of particle trajectories to show in the 3D/projection figures. Default: 80.",
    )
    parser.add_argument(
        "--reservoir-size",
        type=int,
        default=1200,
        help="Reservoir size for candidate trajectories per showcase file. Default: 1200.",
    )
    parser.add_argument(
        "--max-showcase-steps",
        type=int,
        default=500000,
        help="Maximum raw step rows scanned into candidate trajectories for showcase. Default: 500000.",
    )
    parser.add_argument(
        "--only-origin-aligned",
        action="store_true",
        help=(
            "Only generate stageD_origin_aligned_trajectory_showcase_3d.png. "
            "This skips summary statistics and the other Stage D figures."
        ),
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=0,
        help="Process only the first N track files after sorting. Use 0 for all. Default: 0.",
    )
    parser.add_argument(
        "--summary-stride",
        type=int,
        default=1,
        help="Read every Nth row for summary statistics. Use 1 for exact summary. Default: 1.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260617,
        help="Random seed for trajectory sampling. Default: 20260617.",
    )
    return parser.parse_args()


def parse_ratio_tag(tag: str) -> RatioKey:
    match = RATIO_RE.fullmatch(tag)
    if match is None:
        raise ValueError(f"Unsupported ratio folder name: {tag}")
    return RatioKey(
        tag=tag,
        bn_wt=float(match.group("bn")),
        zns_wt=float(match.group("zns")),
    )


def ratio_display_sort_key(ratio: RatioKey) -> tuple[int, float, float]:
    try:
        return (DESIRED_RATIO_ORDER.index(ratio.tag), ratio.bn_wt, ratio.zns_wt)
    except ValueError:
        return (len(DESIRED_RATIO_ORDER), ratio.bn_wt, ratio.zns_wt)


def parse_track_thickness(path: Path) -> float:
    match = FULL_TRACK_RE.fullmatch(path.name) or SLIM_TRACK_RE.fullmatch(path.name)
    if match is None:
        raise ValueError(f"Unsupported Stage B track CSV name: {path.name}")
    return float(match.group("thickness"))


def input_mode_for_path(path: Path) -> str:
    if FULL_TRACK_RE.fullmatch(path.name):
        return "full"
    if SLIM_TRACK_RE.fullmatch(path.name):
        return "slim_zns_only"
    raise ValueError(f"Unsupported Stage B track CSV name: {path.name}")


def thickness_sort_key(value: float) -> tuple[int, float]:
    rounded = round(value)
    is_int_like = math.isclose(value, rounded, rel_tol=0.0, abs_tol=1.0e-9)
    return (0 if is_int_like else 1, value)


def mu_m_text() -> str:
    return "μm"


def safe_float(value: str | None, *, default: float = 0.0) -> float:
    if value is None or value == "":
        return default
    try:
        number = float(value)
    except ValueError:
        return default
    if not math.isfinite(number):
        return default
    return number


def fmt(value: object) -> object:
    if isinstance(value, float):
        return f"{value:.12g}"
    return value


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key, "")) for key in fieldnames})


def has_macro_anchor(row: dict[str, str]) -> bool:
    required = {
        "capture_x_um",
        "capture_y_um",
        "depth_um",
        "local_capture_x_um",
        "local_capture_y_um",
        "local_capture_z_um",
    }
    return required.issubset(row.keys())


def macro_points_from_row(
    row: dict[str, str],
    start: tuple[float, float, float],
    end: tuple[float, float, float],
) -> tuple[tuple[float, float, float] | None, tuple[float, float, float] | None]:
    explicit_fields = {
        "screen_x_pre_um",
        "screen_y_pre_um",
        "screen_depth_pre_um",
        "screen_x_post_um",
        "screen_y_post_um",
        "screen_depth_post_um",
    }
    if explicit_fields.issubset(row.keys()):
        return (
            safe_float(row.get("screen_x_pre_um")),
            safe_float(row.get("screen_y_pre_um")),
            safe_float(row.get("screen_depth_pre_um")),
        ), (
            safe_float(row.get("screen_x_post_um")),
            safe_float(row.get("screen_y_post_um")),
            safe_float(row.get("screen_depth_post_um")),
        )

    if not has_macro_anchor(row):
        return None, None

    capture = (
        safe_float(row.get("capture_x_um"), default=math.nan),
        safe_float(row.get("capture_y_um"), default=math.nan),
        safe_float(row.get("depth_um"), default=math.nan),
    )
    local_capture = (
        safe_float(row.get("local_capture_x_um"), default=math.nan),
        safe_float(row.get("local_capture_y_um"), default=math.nan),
        safe_float(row.get("local_capture_z_um"), default=math.nan),
    )
    if not all(math.isfinite(value) for value in capture + local_capture):
        return None, None

    def to_macro(point: tuple[float, float, float]) -> tuple[float, float, float]:
        return (
            capture[0] + point[0] - local_capture[0],
            capture[1] + point[1] - local_capture[1],
            capture[2] - (point[2] - local_capture[2]),
        )

    return to_macro(start), to_macro(end)


def discover_track_files(
    stageb_root: Path,
    requested_ratios: list[str] | None,
    requested_thicknesses: list[float] | None,
    max_files: int,
) -> tuple[list[RatioKey], list[tuple[RatioKey, float, Path, str]]]:
    if not stageb_root.is_dir():
        raise SystemExit(f"Stage B root not found: {stageb_root}")

    requested_set = set(requested_ratios or [])
    requested_thickness_set = set(requested_thicknesses or [])
    ratios: dict[str, RatioKey] = {}
    entries: list[tuple[RatioKey, float, Path, str]] = []

    for ratio_dir in sorted(stageb_root.iterdir()):
        if not ratio_dir.is_dir():
            continue
        if requested_set and ratio_dir.name not in requested_set:
            continue
        try:
            ratio = parse_ratio_tag(ratio_dir.name)
        except ValueError:
            continue
        full_by_thickness: dict[float, Path] = {}
        slim_by_thickness: dict[float, Path] = {}
        for path in sorted(ratio_dir.glob("*_alpha_li_steps.csv")):
            full_by_thickness[parse_track_thickness(path)] = path
        for path in sorted(ratio_dir.glob("*_zns_track_steps.csv")):
            slim_by_thickness[parse_track_thickness(path)] = path

        thicknesses = sorted(
            set(full_by_thickness).union(slim_by_thickness),
            key=thickness_sort_key,
        )
        if requested_thickness_set:
            thicknesses = [
                thickness
                for thickness in thicknesses
                if any(
                    math.isclose(thickness, requested, rel_tol=0.0, abs_tol=1.0e-9)
                    for requested in requested_thickness_set
                )
            ]
        if max_files > 0:
            thicknesses = thicknesses[:max_files]
        for thickness in thicknesses:
            if thickness in full_by_thickness:
                path = full_by_thickness[thickness]
                input_mode = "full"
            else:
                path = slim_by_thickness[thickness]
                input_mode = "slim_zns_only"
            entries.append((ratio, thickness, path, input_mode))
            ratios[ratio.tag] = ratio

    if requested_set:
        missing = sorted(requested_set.difference(ratios))
        if missing:
            print(f"Warning: requested ratios not found under {stageb_root}: {', '.join(missing)}")

    if not entries:
        raise SystemExit(
            f"No Stage B *_alpha_li_steps.csv or *_zns_track_steps.csv files found in {stageb_root}"
        )

    entries.sort(key=lambda item: (ratio_display_sort_key(item[0]), thickness_sort_key(item[1])))
    ratio_list = sorted(ratios.values(), key=ratio_display_sort_key)
    return ratio_list, entries


def read_step(row: dict[str, str], ratio: RatioKey, thickness_um: float) -> TrackStep:
    has_unwrapped = all(
        name in row
        for name in {
            "unwrapped_x_pre_um",
            "unwrapped_y_pre_um",
            "unwrapped_z_pre_um",
            "unwrapped_x_post_um",
            "unwrapped_y_post_um",
            "unwrapped_z_post_um",
        }
    )
    prefix = "unwrapped_" if has_unwrapped else ""
    start = (
        safe_float(row.get(f"{prefix}x_pre_um")),
        safe_float(row.get(f"{prefix}y_pre_um")),
        safe_float(row.get(f"{prefix}z_pre_um")),
    )
    end = (
        safe_float(row.get(f"{prefix}x_post_um")),
        safe_float(row.get(f"{prefix}y_post_um")),
        safe_float(row.get(f"{prefix}z_post_um")),
    )
    macro_start, macro_end = macro_points_from_row(row, start, end)
    return TrackStep(
        particle=row["particle"].strip(),
        source_event_uid=row["source_event_uid"].strip(),
        physical_event_uid=row["physical_event_uid"].strip(),
        thickness_um=thickness_um,
        phase_pre=row.get("phase_pre", "unknown").strip() or "unknown",
        phase_post=row.get("phase_post", "unknown").strip() or "unknown",
        start=start,
        end=end,
        macro_start=macro_start,
        macro_end=macro_end,
        step_len_um=safe_float(row.get("step_len_um")),
        edep_kev=safe_float(row.get("edep_keV")),
        ekin_pre_kev=safe_float(row.get("ekin_pre_keV")),
        ekin_post_kev=safe_float(row.get("ekin_post_keV")),
        trajectory_weight=safe_float(row.get("trajectory_weight"), default=1.0),
    )


def update_trajectory(traj: Trajectory, step: TrackStep) -> None:
    traj.steps.append(step)
    traj.total_edep_kev += step.edep_kev
    traj.total_len_um += step.step_len_um
    traj.trajectory_weight = step.trajectory_weight


def collect_summary(
    entries: list[tuple[RatioKey, float, Path, str]],
    summary_stride: int,
) -> list[dict[str, object]]:
    if summary_stride <= 0:
        raise SystemExit("--summary-stride must be >= 1")

    records: list[dict[str, object]] = []
    required = {
        "source_event_uid",
        "particle",
        "step_len_um",
        "edep_keV",
        "ekin_pre_keV",
        "trajectory_weight",
    }

    for ratio, thickness_um, path, input_mode in entries:
        particle_rows: defaultdict[str, int] = defaultdict(int)
        particle_step_len: defaultdict[str, float] = defaultdict(float)
        particle_edep: defaultdict[str, float] = defaultdict(float)
        particle_weighted_edep: defaultdict[str, float] = defaultdict(float)
        particle_weight_sum: defaultdict[str, float] = defaultdict(float)
        particle_phase_len: defaultdict[str, defaultdict[str, float]] = defaultdict(lambda: defaultdict(float))
        particle_source_ids: defaultdict[str, set[str]] = defaultdict(set)
        particle_source_edep: defaultdict[str, defaultdict[str, float]] = defaultdict(lambda: defaultdict(float))
        particle_source_len: defaultdict[str, defaultdict[str, float]] = defaultdict(lambda: defaultdict(float))
        particle_source_phase_len: defaultdict[str, defaultdict[str, defaultdict[str, float]]] = defaultdict(
            lambda: defaultdict(lambda: defaultdict(float))
        )
        particle_source_weight: defaultdict[str, dict[str, float]] = defaultdict(dict)
        source_edep: defaultdict[str, float] = defaultdict(float)
        source_len: defaultdict[str, float] = defaultdict(float)
        source_phase_len: defaultdict[str, defaultdict[str, float]] = defaultdict(lambda: defaultdict(float))
        source_weight: dict[str, float] = {}
        rows_scanned = 0
        rows_used = 0

        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            missing = required.difference(reader.fieldnames or [])
            if missing:
                raise SystemExit(f"{path} is missing columns: {', '.join(sorted(missing))}")

            for rows_scanned, row in enumerate(reader, start=1):
                if (rows_scanned - 1) % summary_stride != 0:
                    continue
                particle = row["particle"].strip()
                if particle not in PARTICLE_COLORS:
                    continue
                source_uid = row["source_event_uid"].strip()
                step_len = safe_float(row.get("step_len_um"))
                edep = safe_float(row.get("edep_keV"))
                weight = safe_float(row.get("trajectory_weight"), default=1.0)
                phase = row.get("phase_pre", "unknown").strip() or "unknown"

                rows_used += 1
                particle_rows[particle] += 1
                particle_step_len[particle] += step_len
                particle_phase_len[particle][phase] += step_len
                particle_edep[particle] += edep
                particle_weighted_edep[particle] += edep * weight
                particle_weight_sum[particle] += weight
                particle_source_ids[particle].add(source_uid)
                particle_source_edep[particle][source_uid] += edep
                particle_source_len[particle][source_uid] += step_len
                particle_source_phase_len[particle][source_uid][phase] += step_len
                particle_source_weight[particle][source_uid] = weight
                source_edep[source_uid] += edep
                source_len[source_uid] += step_len
                source_phase_len[source_uid][phase] += step_len
                source_weight[source_uid] = weight

        all_source_count = len(source_edep)
        all_edep_values = np.asarray(list(source_edep.values()), dtype=float)
        all_len_values = np.asarray(list(source_len.values()), dtype=float)
        all_zns_len_values = np.asarray(
            [source_phase_len[uid].get("ZnS", 0.0) for uid in source_edep],
            dtype=float,
        )
        all_bn_len_values = np.asarray(
            [source_phase_len[uid].get("BN", 0.0) for uid in source_edep],
            dtype=float,
        )
        all_binder_len_values = np.asarray(
            [source_phase_len[uid].get("binder_void", 0.0) for uid in source_edep],
            dtype=float,
        )
        all_weighted_edep = sum(source_edep[uid] * source_weight.get(uid, 1.0) for uid in source_edep)
        all_weight_sum = sum(source_weight.get(uid, 1.0) for uid in source_edep)
        records.append(
            {
                "ratio_tag": ratio.tag,
                "bn_wt": ratio.bn_wt,
                "zns_wt": ratio.zns_wt,
                "thickness_um": thickness_um,
                "input_mode": input_mode,
                "particle": "all",
                "track_file": path.as_posix(),
                "rows_scanned": rows_scanned,
                "rows_used": rows_used,
                "n_source_trajectories": all_source_count,
                "total_step_length_um": float(np.sum(all_len_values)) if all_len_values.size else 0.0,
                "mean_step_length_um": (
                    sum(particle_step_len.values()) / rows_used if rows_used > 0 else math.nan
                ),
                "total_edep_keV": float(np.sum(all_edep_values)) if all_edep_values.size else 0.0,
                "weighted_total_edep_keV": all_weighted_edep,
                "weighted_mean_edep_per_source_keV": (
                    all_weighted_edep / all_weight_sum if all_weight_sum > 0 else math.nan
                ),
                "mean_edep_per_source_keV": float(np.mean(all_edep_values)) if all_edep_values.size else math.nan,
                "median_edep_per_source_keV": (
                    float(np.median(all_edep_values)) if all_edep_values.size else math.nan
                ),
                "mean_path_length_per_source_um": (
                    float(np.mean(all_len_values)) if all_len_values.size else math.nan
                ),
                "mean_zns_path_length_per_source_um": (
                    float(np.mean(all_zns_len_values)) if all_zns_len_values.size else math.nan
                ),
                "mean_bn_path_length_per_source_um": (
                    float(np.mean(all_bn_len_values)) if all_bn_len_values.size else math.nan
                ),
                "mean_binder_path_length_per_source_um": (
                    float(np.mean(all_binder_len_values)) if all_binder_len_values.size else math.nan
                ),
            }
        )

        for particle in sorted(particle_rows):
            source_ids = particle_source_ids[particle]
            edep_values = np.asarray(
                [particle_source_edep[particle][uid] for uid in source_ids],
                dtype=float,
            )
            len_values = np.asarray(
                [particle_source_len[particle][uid] for uid in source_ids],
                dtype=float,
            )
            zns_len_values = np.asarray(
                [particle_source_phase_len[particle][uid].get("ZnS", 0.0) for uid in source_ids],
                dtype=float,
            )
            bn_len_values = np.asarray(
                [particle_source_phase_len[particle][uid].get("BN", 0.0) for uid in source_ids],
                dtype=float,
            )
            binder_len_values = np.asarray(
                [particle_source_phase_len[particle][uid].get("binder_void", 0.0) for uid in source_ids],
                dtype=float,
            )
            weight_sum = particle_weight_sum[particle]
            source_weight_sum = sum(particle_source_weight[particle].get(uid, 1.0) for uid in source_ids)
            source_weighted_edep = sum(
                particle_source_edep[particle][uid] * particle_source_weight[particle].get(uid, 1.0)
                for uid in source_ids
            )
            records.append(
                {
                    "ratio_tag": ratio.tag,
                    "bn_wt": ratio.bn_wt,
                    "zns_wt": ratio.zns_wt,
                    "thickness_um": thickness_um,
                    "input_mode": input_mode,
                    "particle": particle,
                    "track_file": path.as_posix(),
                    "rows_scanned": rows_scanned,
                    "rows_used": particle_rows[particle],
                    "n_source_trajectories": len(source_ids),
                    "total_step_length_um": particle_step_len[particle],
                    "mean_step_length_um": particle_step_len[particle] / particle_rows[particle],
                    "total_edep_keV": particle_edep[particle],
                    "weighted_total_edep_keV": particle_weighted_edep[particle],
                    "weighted_mean_edep_per_source_keV": (
                        source_weighted_edep / source_weight_sum if source_weight_sum > 0 else math.nan
                    ),
                    "mean_edep_per_source_keV": (
                        float(np.mean(edep_values)) if edep_values.size else math.nan
                    ),
                    "median_edep_per_source_keV": (
                        float(np.median(edep_values)) if edep_values.size else math.nan
                    ),
                    "mean_path_length_per_source_um": (
                        float(np.mean(len_values)) if len_values.size else math.nan
                    ),
                    "mean_zns_path_length_per_source_um": (
                        float(np.mean(zns_len_values)) if zns_len_values.size else math.nan
                    ),
                    "mean_bn_path_length_per_source_um": (
                        float(np.mean(bn_len_values)) if bn_len_values.size else math.nan
                    ),
                    "mean_binder_path_length_per_source_um": (
                        float(np.mean(binder_len_values)) if binder_len_values.size else math.nan
                    ),
                }
            )

        print(
            f"Summarized {ratio.tag} t={thickness_um:g} um ({input_mode}): "
            f"{rows_used:,} used rows from {rows_scanned:,} scanned rows, "
            f"{all_source_count:,} source trajectories"
        )

    return records


def choose_showcase_entry(
    entries: list[tuple[RatioKey, float, Path, str]],
    showcase_ratio: str | None,
    showcase_thickness: float | None,
) -> tuple[RatioKey, float, Path, str]:
    ratio_entries = entries
    if showcase_ratio is not None:
        ratio_entries = [item for item in ratio_entries if item[0].tag == showcase_ratio]
        if not ratio_entries:
            raise SystemExit(f"--showcase-ratio not found in processed Stage B data: {showcase_ratio}")

    if showcase_thickness is not None:
        exact = [
            item
            for item in ratio_entries
            if math.isclose(item[1], showcase_thickness, rel_tol=0.0, abs_tol=1.0e-9)
        ]
        if not exact:
            available = ", ".join(f"{item[0].tag}:{item[1]:g}" for item in ratio_entries)
            raise SystemExit(
                f"--showcase-thickness {showcase_thickness:g} not found. "
                f"Available processed entries: {available}"
            )
        return exact[-1]

    return sorted(ratio_entries, key=lambda item: thickness_sort_key(item[1]))[-1]


def collect_showcase_trajectories(
    ratio: RatioKey,
    thickness_um: float,
    path: Path,
    input_mode: str,
    reservoir_size: int,
    sample_count: int,
    max_steps: int,
    seed: int,
) -> list[Trajectory]:
    if reservoir_size <= 0:
        raise SystemExit("--reservoir-size must be > 0")
    if sample_count <= 0:
        raise SystemExit("--sample-trajectories must be > 0")
    if max_steps <= 0:
        raise SystemExit("--max-showcase-steps must be > 0")

    required = {
        "physical_event_uid",
        "source_event_uid",
        "particle",
        "x_pre_um",
        "y_pre_um",
        "z_pre_um",
        "x_post_um",
        "y_post_um",
        "z_post_um",
        "step_len_um",
        "edep_keV",
        "ekin_pre_keV",
        "ekin_post_keV",
        "trajectory_weight",
    }
    rng = random.Random(seed)
    active: OrderedDict[str, Trajectory] = OrderedDict()
    seen_sources = 0
    rows_scanned = 0

    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"{path} is missing columns: {', '.join(sorted(missing))}")

        for rows_scanned, row in enumerate(reader, start=1):
            if rows_scanned > max_steps:
                break
            source_uid = row["source_event_uid"].strip()
            particle = row["particle"].strip()
            if particle not in PARTICLE_COLORS:
                continue
            track_id = row.get("trackID", "").strip()
            trajectory_key = f"{source_uid}|{particle}|{track_id}"
            traj = active.get(trajectory_key)
            if traj is None:
                seen_sources += 1
                keep = len(active) < reservoir_size
                if not keep:
                    replacement_index = rng.randrange(seen_sources)
                    keep = replacement_index < reservoir_size
                if not keep:
                    continue
                if len(active) >= reservoir_size:
                    remove_key = rng.choice(list(active.keys()))
                    active.pop(remove_key, None)
                traj = Trajectory(
                    source_event_uid=source_uid,
                    physical_event_uid=row["physical_event_uid"].strip(),
                    particle=particle,
                    ratio_tag=ratio.tag,
                    thickness_um=thickness_um,
                    input_mode=input_mode,
                )
                active[trajectory_key] = traj

            step = read_step(row, ratio, thickness_um)
            update_trajectory(traj, step)

    candidates = [traj for traj in active.values() if traj.n_steps > 0]
    if not candidates:
        raise SystemExit(f"No showcase trajectories collected from {path}")

    candidates.sort(key=lambda traj: (traj.total_len_um, traj.total_edep_kev, traj.n_steps), reverse=True)
    selected: list[Trajectory] = []
    per_particle_target = max(1, sample_count // 2)
    for particle in ["alpha", "Li7"]:
        selected.extend([traj for traj in candidates if traj.particle == particle][:per_particle_target])
    if len(selected) < sample_count:
        selected_ids = {id(traj) for traj in selected}
        selected.extend([traj for traj in candidates if id(traj) not in selected_ids][: sample_count - len(selected)])
    selected = selected[:sample_count]
    selected.sort(key=lambda traj: (traj.particle, -traj.total_len_um))

    print(
        f"Showcase source: {ratio.tag} t={thickness_um:g} um, {path.name}; "
        f"selected {len(selected)} trajectories from {len(candidates)} candidates "
        f"after scanning {rows_scanned:,} rows"
    )
    return selected


def trajectory_segments_3d(trajectories: list[Trajectory]) -> tuple[np.ndarray, list[str], list[str], np.ndarray]:
    segments = []
    particles = []
    phases = []
    edeps = []
    for traj in trajectories:
        for step in traj.steps:
            segments.append([step.start, step.end])
            particles.append(traj.particle)
            phases.append(step.phase_pre)
            edeps.append(step.edep_kev)
    return np.asarray(segments, dtype=float), particles, phases, np.asarray(edeps, dtype=float)


def trajectory_bounds(trajectories: list[Trajectory]) -> tuple[np.ndarray, np.ndarray]:
    points = []
    for traj in trajectories:
        for step in traj.steps:
            points.append(step.start)
            points.append(step.end)
    arr = np.asarray(points, dtype=float)
    if arr.size == 0:
        return np.zeros(3), np.ones(3)
    mins = np.min(arr, axis=0)
    maxs = np.max(arr, axis=0)
    pad = np.maximum((maxs - mins) * 0.08, 0.5)
    return mins - pad, maxs + pad


def origin_aligned_segments(
    trajectories: list[Trajectory],
) -> tuple[np.ndarray, list[str], list[str], np.ndarray]:
    segments = []
    particles = []
    phases = []
    edeps = []
    for traj in trajectories:
        if not traj.steps:
            continue
        origin = np.asarray(traj.steps[0].start, dtype=float)
        for step in traj.steps:
            start = tuple(np.asarray(step.start, dtype=float) - origin)
            end = tuple(np.asarray(step.end, dtype=float) - origin)
            segments.append([start, end])
            particles.append(traj.particle)
            phases.append(step.phase_pre)
            edeps.append(step.edep_kev)
    return np.asarray(segments, dtype=float), particles, phases, np.asarray(edeps, dtype=float)


def origin_aligned_bounds(segments: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if segments.size == 0:
        return np.asarray([-1.0, -1.0, -1.0]), np.asarray([1.0, 1.0, 1.0])
    points = segments.reshape((-1, 3))
    max_abs = float(np.max(np.abs(points)))
    half = max(max_abs * 1.12, 1.0)
    return np.asarray([-half, -half, -half]), np.asarray([half, half, half])


def macro_trajectory_bounds(trajectories: list[Trajectory]) -> tuple[np.ndarray, np.ndarray] | None:
    points = []
    for traj in trajectories:
        for step in traj.steps:
            if step.macro_start is None or step.macro_end is None:
                continue
            points.append(step.macro_start)
            points.append(step.macro_end)
    if not points:
        return None
    arr = np.asarray(points, dtype=float)
    mins = np.min(arr, axis=0)
    maxs = np.max(arr, axis=0)
    pad = np.maximum((maxs - mins) * 0.08, np.asarray([50.0, 50.0, 2.0]))
    return mins - pad, maxs + pad


def set_3d_equalish(ax, mins: np.ndarray, maxs: np.ndarray) -> None:
    centers = 0.5 * (mins + maxs)
    ranges = np.maximum(maxs - mins, 1.0)
    half = 0.5 * float(np.max(ranges))
    ax.set_xlim(centers[0] - half, centers[0] + half)
    ax.set_ylim(centers[1] - half, centers[1] + half)
    ax.set_zlim(centers[2] - half, centers[2] + half)
    try:
        ax.set_box_aspect((1, 1, 1))
    except AttributeError:
        pass


def style_axes(ax) -> None:
    ax.tick_params(direction="in", which="both", top=True, right=True, width=1.2)
    for spine in ax.spines.values():
        spine.set_linewidth(1.4)


def plot_trajectory_showcase_3d(
    trajectories: list[Trajectory],
    output_path: Path,
    ratio: RatioKey,
    thickness_um: float,
) -> None:
    segments, particles, phases, edeps = trajectory_segments_3d(trajectories)
    mins, maxs = trajectory_bounds(trajectories)
    edep_norm = Normalize(vmin=0.0, vmax=max(float(np.max(edeps)) if edeps.size else 1.0, 1.0))

    fig = plt.figure(figsize=(10.2, 8.0))
    ax = fig.add_subplot(111, projection="3d")
    for phase in PHASE_ORDER:
        mask = np.asarray([value == phase for value in phases], dtype=bool)
        if not np.any(mask):
            continue
        widths = np.asarray(
            [
                PARTICLE_LINEWIDTH.get(particle, 1.0) + 1.2 * edep_norm(edep)
                for particle, edep, keep in zip(particles, edeps, mask)
                if keep
            ],
            dtype=float,
        )
        collection = Line3DCollection(
            segments[mask],
            colors=PHASE_COLORS.get(phase, PHASE_COLORS["unknown"]),
            linewidths=widths,
            alpha=0.64 if phase == "binder_void" else 0.82,
        )
        ax.add_collection3d(collection)

    starts = np.asarray([traj.steps[0].start for traj in trajectories if traj.steps], dtype=float)
    if starts.size:
        ax.scatter(starts[:, 0], starts[:, 1], starts[:, 2], c="#111827", s=18, marker="o", label="start")

    set_3d_equalish(ax, mins, maxs)
    ax.set_xlabel(f"x ({mu_m_text()})")
    ax.set_ylabel(f"y ({mu_m_text()})")
    ax.set_zlabel(f"z ({mu_m_text()})")
    ax.set_title(
        f"Alpha and Li-7 Particle Trajectories\n"
        f"BN:ZnS = {ratio.display_tag}, t = {thickness_um:g} {mu_m_text()}"
    )
    ax.view_init(elev=22, azim=-52)
    ax.grid(True, alpha=0.16)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    phase_handles = [
        Line2D([0], [0], color=PHASE_COLORS[phase], lw=3.0, label=PHASE_LABELS.get(phase, phase))
        for phase in PHASE_ORDER
        if phase in set(phases)
    ]
    particle_handles = [
        Line2D([0], [0], color="#2d3436", lw=3.2, label="alpha width"),
        Line2D([0], [0], color="#2d3436", lw=2.0, label="Li7 width"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#111827", markersize=5, label="track start"),
    ]
    ax.legend(handles=phase_handles + particle_handles, frameon=False, loc="upper left", fontsize=8.5)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_origin_aligned_trajectory_showcase_3d(
    trajectories: list[Trajectory],
    output_path: Path,
    ratio: RatioKey,
    thickness_um: float,
) -> None:
    segments, particles, phases, edeps = origin_aligned_segments(trajectories)
    mins, maxs = origin_aligned_bounds(segments)
    half_range = float(maxs[0])
    edep_norm = Normalize(vmin=0.0, vmax=max(float(np.max(edeps)) if edeps.size else 1.0, 1.0))
    ticks = np.linspace(-half_range, half_range, 5)

    fig = plt.figure(figsize=(13.2, 6.7))
    axes = [
        fig.add_subplot(1, 2, 1, projection="3d"),
        fig.add_subplot(1, 2, 2, projection="3d"),
    ]
    particle_titles = [("alpha", "Alpha Particle Trajectories"), ("Li7", "Li-7 Particle Trajectories")]

    for ax, (particle_name, panel_title) in zip(axes, particle_titles):
        particle_mask = np.asarray([value == particle_name for value in particles], dtype=bool)
        for phase in PHASE_ORDER:
            mask = particle_mask & np.asarray([value == phase for value in phases], dtype=bool)
            if not np.any(mask):
                continue
            widths = np.asarray(
                [0.72 + 0.72 * edep_norm(edep) for edep, keep in zip(edeps, mask) if keep],
                dtype=float,
            )
            collection = Line3DCollection(
                segments[mask],
                colors=PHASE_COLORS.get(phase, PHASE_COLORS["unknown"]),
                linewidths=widths,
                alpha=0.48 if phase == "binder_void" else 0.72,
            )
            ax.add_collection3d(collection)

        ax.set_xlim(mins[0], maxs[0])
        ax.set_ylim(mins[1], maxs[1])
        ax.set_zlim(mins[2], maxs[2])
        try:
            ax.set_box_aspect((1, 1, 1))
        except AttributeError:
            pass
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_zticks(ticks)
        ax.set_xlabel(f"relative x ({mu_m_text()})", labelpad=9, fontsize=11, weight="bold")
        ax.set_ylabel(f"relative y ({mu_m_text()})", labelpad=9, fontsize=11, weight="bold")
        ax.set_zlabel(f"relative z ({mu_m_text()})", labelpad=9, fontsize=11, weight="bold")
        ax.set_title(panel_title)
        ax.view_init(elev=24, azim=-45)
        ax.grid(True, alpha=0.28, linewidth=0.8)
        ax.tick_params(labelsize=8, pad=1)
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False

    fig.suptitle(
        f"Origin-Aligned Alpha and Li-7 Particle Trajectories\n"
        f"BN:ZnS = {ratio.display_tag}, t = {thickness_um:g} {mu_m_text()}",
        y=0.98,
        fontsize=14,
    )

    phases_present = set(phases)
    phase_handles = [
        Line2D([0], [0], color=PHASE_COLORS[phase], lw=2.4, label=PHASE_LABELS.get(phase, phase))
        for phase in PHASE_ORDER
        if phase in phases_present
    ]
    fig.legend(
        handles=phase_handles,
        loc="lower center",
        ncol=min(len(phase_handles), 6),
        frameon=False,
        bbox_to_anchor=(0.5, 0.01),
        fontsize=9,
    )
    fig.subplots_adjust(left=0.03, right=0.99, top=0.84, bottom=0.12, wspace=0.08)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_macro_trajectory_showcase_3d(
    trajectories: list[Trajectory],
    output_path: Path,
    ratio: RatioKey,
    thickness_um: float,
) -> bool:
    bounds = macro_trajectory_bounds(trajectories)
    if bounds is None:
        print("Skipping macro trajectory showcase: selected trajectories do not contain macro anchor columns")
        return False
    mins, maxs = bounds

    segments_by_phase: defaultdict[str, list[list[tuple[float, float, float]]]] = defaultdict(list)
    widths_by_phase: defaultdict[str, list[float]] = defaultdict(list)
    starts = []
    for traj in trajectories:
        if not traj.steps or traj.steps[0].macro_start is None:
            continue
        starts.append(traj.steps[0].macro_start)
        for step in traj.steps:
            if step.macro_start is None or step.macro_end is None:
                continue
            segments_by_phase[step.phase_pre].append([step.macro_start, step.macro_end])
            widths_by_phase[step.phase_pre].append(PARTICLE_LINEWIDTH.get(traj.particle, 1.0) + 0.6)

    if not any(segments_by_phase.values()):
        print("Skipping macro trajectory showcase: no macro segments after conversion")
        return False

    fig = plt.figure(figsize=(10.6, 8.0))
    ax = fig.add_subplot(111, projection="3d")
    for phase in PHASE_ORDER:
        segments = segments_by_phase.get(phase, [])
        if not segments:
            continue
        collection = Line3DCollection(
            np.asarray(segments, dtype=float),
            colors=PHASE_COLORS.get(phase, PHASE_COLORS["unknown"]),
            linewidths=widths_by_phase[phase],
            alpha=0.58 if phase == "binder_void" else 0.78,
        )
        ax.add_collection3d(collection)

    if starts:
        arr = np.asarray(starts, dtype=float)
        ax.scatter(arr[:, 0], arr[:, 1], arr[:, 2], c="#111827", s=14, marker="o", label="track start")

    ax.set_xlim(mins[0], maxs[0])
    ax.set_ylim(mins[1], maxs[1])
    ax.set_zlim(mins[2], maxs[2])
    try:
        ax.set_box_aspect(
            (
                max(float(maxs[0] - mins[0]), 1.0),
                max(float(maxs[1] - mins[1]), 1.0),
                max(float(maxs[2] - mins[2]), 1.0),
            )
        )
    except AttributeError:
        pass
    ax.invert_zaxis()
    ax.set_xlabel(f"macro x ({mu_m_text()})")
    ax.set_ylabel(f"macro y ({mu_m_text()})")
    ax.set_zlabel(f"capture depth ({mu_m_text()})")
    ax.set_title(
        f"Macroscopic Alpha and Li-7 Particle Trajectories\n"
        f"BN:ZnS = {ratio.display_tag}, t = {thickness_um:g} {mu_m_text()}"
    )
    ax.view_init(elev=24, azim=-54)
    ax.grid(True, alpha=0.16)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    phases_present = {phase for phase, segments in segments_by_phase.items() if segments}
    handles = [
        Line2D([0], [0], color=PHASE_COLORS[phase], lw=3.0, label=PHASE_LABELS.get(phase, phase))
        for phase in PHASE_ORDER
        if phase in phases_present
    ]
    handles.extend(
        [
            Line2D([0], [0], color="#2d3436", lw=3.2, label="alpha width"),
            Line2D([0], [0], color="#2d3436", lw=2.0, label="Li7 width"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor="#111827", markersize=5, label="track start"),
        ]
    )
    ax.legend(handles=handles, frameon=False, loc="upper left", fontsize=8.5)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True


def plot_trajectory_projections(
    trajectories: list[Trajectory],
    output_path: Path,
    ratio: RatioKey,
    thickness_um: float,
) -> None:
    plane_defs = [
        ("x", "y", 0, 1),
        ("x", "z", 0, 2),
        ("y", "z", 1, 2),
    ]
    mins, maxs = trajectory_bounds(trajectories)
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.2), constrained_layout=False)

    for ax, (x_label, y_label, x_idx, y_idx) in zip(axes, plane_defs):
        for traj in trajectories:
            phase_segments: defaultdict[str, list[list[tuple[float, float]]]] = defaultdict(list)
            phase_widths: defaultdict[str, list[float]] = defaultdict(list)
            for step in traj.steps:
                phase_segments[step.phase_pre].append(
                    [
                        (step.start[x_idx], step.start[y_idx]),
                        (step.end[x_idx], step.end[y_idx]),
                    ]
                )
                phase_widths[step.phase_pre].append(
                    PARTICLE_LINEWIDTH.get(traj.particle, 1.0) + 0.9 * min(step.edep_kev / 500.0, 1.0)
                )
            for phase, segments in phase_segments.items():
                ax.add_collection(
                    LineCollection(
                        segments,
                        colors=PHASE_COLORS.get(phase, PHASE_COLORS["unknown"]),
                        linewidths=phase_widths[phase],
                        alpha=0.62 if phase == "binder_void" else 0.8,
                    )
                )
            start = traj.steps[0].start
            ax.scatter(start[x_idx], start[y_idx], c="#111827", s=10, marker="o", zorder=4)

        ax.set_xlim(mins[x_idx], maxs[x_idx])
        ax.set_ylim(mins[y_idx], maxs[y_idx])
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(f"{x_label} ({mu_m_text()})")
        ax.set_ylabel(f"{y_label} ({mu_m_text()})")
        ax.set_title(f"{x_label}-{y_label} projection")
        style_axes(ax)

    phases_present = {step.phase_pre for traj in trajectories for step in traj.steps}
    handles = [
        Line2D([0], [0], color=PHASE_COLORS[phase], lw=2.8, label=PHASE_LABELS.get(phase, phase))
        for phase in PHASE_ORDER
        if phase in phases_present
    ]
    handles.extend(
        [
            Line2D([0], [0], color="#2d3436", lw=3.2, label="alpha width"),
            Line2D([0], [0], color="#2d3436", lw=2.0, label="Li7 width"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor="#111827", markersize=5, label="start"),
        ]
    )
    fig.legend(handles=handles, loc="lower center", ncol=min(len(handles), 7), frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle(
        f"Alpha and Li-7 Trajectory Projections, BN:ZnS = {ratio.display_tag}, "
        f"t = {thickness_um:g} {mu_m_text()}",
        y=0.98,
        fontsize=14,
    )
    fig.subplots_adjust(left=0.06, right=0.99, top=0.82, bottom=0.18, wspace=0.24)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def records_to_series(
    records: list[dict[str, object]],
    particle: str,
    field: str,
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    grouped: defaultdict[str, list[tuple[float, float]]] = defaultdict(list)
    for row in records:
        if row["particle"] != particle:
            continue
        value = safe_float(str(row.get(field, "")), default=math.nan)
        if not math.isfinite(value):
            continue
        grouped[str(row["ratio_tag"])].append((float(row["thickness_um"]), value))
    result: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for ratio_tag, pairs in grouped.items():
        pairs.sort(key=lambda item: thickness_sort_key(item[0]))
        result[ratio_tag] = (
            np.asarray([item[0] for item in pairs], dtype=float),
            np.asarray([item[1] for item in pairs], dtype=float),
        )
    return result


def plot_summary_lines(
    ratios: list[RatioKey],
    records: list[dict[str, object]],
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.9), constrained_layout=False)
    fields = [
        ("mean_edep_per_source_keV", "Mean edep per source trajectory (keV)"),
        ("mean_path_length_per_source_um", f"Mean total recorded path length ({mu_m_text()})"),
    ]
    line_styles = {"alpha": "-", "Li7": "--", "all": ":"}
    particle_labels = {"alpha": "alpha", "Li7": "Li7", "all": "all"}
    ratio_colors = {
        ratio.tag: color
        for ratio, color in zip(
            ratios,
            ["#0b3c79", "#c0392b", "#188977", "#8e44ad", "#d68910", "#2d3436"],
        )
    }

    for ax, (field_name, ylabel) in zip(axes, fields):
        for particle in ["alpha", "Li7", "all"]:
            series = records_to_series(records, particle, field_name)
            for ratio in ratios:
                if ratio.tag not in series:
                    continue
                x, y = series[ratio.tag]
                label = f"{ratio.display_tag} {particle_labels[particle]}"
                ax.plot(
                    x,
                    y,
                    linestyle=line_styles[particle],
                    color=ratio_colors[ratio.tag],
                    linewidth=1.9 if particle != "all" else 1.4,
                    alpha=0.88 if particle != "all" else 0.62,
                    label=label,
                )
        ax.set_xscale("log")
        ax.set_xlabel(f"Thickness t ({mu_m_text()})")
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both", alpha=0.18, linewidth=0.7)
        style_axes(ax)

    axes[0].set_title("Energy deposition")
    axes[1].set_title("Track length")
    handles, labels = axes[0].get_legend_handles_labels()
    if len(handles) <= 12:
        axes[1].legend(handles, labels, frameon=False, fontsize=8.2, loc="best")
    else:
        phase_handles = [
            Line2D([0], [0], color="#2d3436", linestyle=line_styles["alpha"], lw=2, label="alpha"),
            Line2D([0], [0], color="#2d3436", linestyle=line_styles["Li7"], lw=2, label="Li7"),
            Line2D([0], [0], color="#2d3436", linestyle=line_styles["all"], lw=2, label="all"),
        ]
        axes[1].legend(handles=phase_handles, frameon=False, loc="best")

    fig.suptitle("Alpha and Li-7 Trajectory Statistics", y=0.98, fontsize=14)
    fig.subplots_adjust(left=0.07, right=0.99, top=0.82, bottom=0.14, wspace=0.24)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_1000um_mean_track_length(
    ratios: list[RatioKey],
    records: list[dict[str, object]],
    output_path: Path,
) -> None:
    target_thickness = 1000.0
    ratio_tags = [ratio.tag for ratio in ratios]
    ratio_labels = [ratio.display_tag for ratio in ratios]
    alpha_values = []
    li_values = []
    mode_labels = []

    for ratio in ratios:
        alpha = math.nan
        li7 = math.nan
        mode = ""
        for row in records:
            if str(row["ratio_tag"]) != ratio.tag:
                continue
            if not math.isclose(float(row["thickness_um"]), target_thickness, rel_tol=0.0, abs_tol=1.0e-9):
                continue
            if row["particle"] == "alpha":
                alpha = safe_float(str(row.get("mean_path_length_per_source_um", "")), default=math.nan)
                mode = str(row.get("input_mode", "")) or mode
            elif row["particle"] == "Li7":
                li7 = safe_float(str(row.get("mean_path_length_per_source_um", "")), default=math.nan)
                mode = str(row.get("input_mode", "")) or mode
        alpha_values.append(alpha)
        li_values.append(li7)
        mode_labels.append(mode)

    if not any(math.isfinite(value) for value in alpha_values + li_values):
        print("Skipping 1000 um mean track-length plot: no matching records")
        return

    x = np.arange(len(ratio_tags), dtype=float)
    width = 0.34
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    ax.bar(
        x - width / 2,
        alpha_values,
        width=width,
        color=PARTICLE_COLORS["alpha"],
        label="alpha",
        alpha=0.9,
    )
    ax.bar(
        x + width / 2,
        li_values,
        width=width,
        color=PARTICLE_COLORS["Li7"],
        label="Li7",
        alpha=0.9,
    )

    for idx, mode in enumerate(mode_labels):
        if not mode:
            continue
        label = "full" if mode == "full" else "ZnS-only"
        ymax = max(
            alpha_values[idx] if math.isfinite(alpha_values[idx]) else 0.0,
            li_values[idx] if math.isfinite(li_values[idx]) else 0.0,
        )
        ax.text(
            x[idx],
            ymax * 1.015 if ymax > 0 else 0.02,
            label,
            ha="center",
            va="bottom",
            fontsize=8,
            color="#374151",
        )

    ax.set_xticks(x)
    ax.set_xticklabels(ratio_labels)
    ax.set_xlabel("BN:ZnS weight ratio")
    ax.set_ylabel(f"Mean total particle track length ({mu_m_text()})")
    ax.set_title(f"Mean alpha/Li trajectory length at t = 1000 {mu_m_text()}")
    ax.grid(True, axis="y", alpha=0.18, linewidth=0.7)
    ax.legend(frameon=False, loc="best")
    style_axes(ax)
    finite = [value for value in alpha_values + li_values if math.isfinite(value)]
    if finite:
        ax.set_ylim(0.0, max(finite) * 1.18)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    if args.sample_trajectories <= 0:
        raise SystemExit("--sample-trajectories must be > 0")
    if args.max_files < 0:
        raise SystemExit("--max-files must be >= 0")

    project_root = Path(args.project_root).resolve() if args.project_root else SCRIPT_DIR.parents[1]
    stageb_root = Path(args.stageb_root).resolve() if args.stageb_root else project_root / "Output" / "stageB"
    output_dir = Path(args.output_dir).resolve() if args.output_dir else SCRIPT_DIR / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    ratios, entries = discover_track_files(stageb_root, args.ratios, args.thicknesses, args.max_files)
    if args.only_origin_aligned:
        showcase_ratio, showcase_thickness, showcase_path, showcase_input_mode = choose_showcase_entry(
            entries,
            args.showcase_ratio,
            args.showcase_thickness,
        )
        trajectories = collect_showcase_trajectories(
            showcase_ratio,
            showcase_thickness,
            showcase_path,
            showcase_input_mode,
            args.reservoir_size,
            args.sample_trajectories,
            args.max_showcase_steps,
            args.seed,
        )
        origin_path = output_dir / "stageD_origin_aligned_trajectory_showcase_3d.png"
        plot_origin_aligned_trajectory_showcase_3d(
            trajectories,
            origin_path,
            showcase_ratio,
            showcase_thickness,
        )
        print(f"Wrote {origin_path}")
        return 0

    records = collect_summary(entries, args.summary_stride)

    summary_fields = [
        "ratio_tag",
        "bn_wt",
        "zns_wt",
        "thickness_um",
        "input_mode",
        "particle",
        "track_file",
        "rows_scanned",
        "rows_used",
        "n_source_trajectories",
        "total_step_length_um",
        "mean_step_length_um",
        "total_edep_keV",
        "weighted_total_edep_keV",
        "weighted_mean_edep_per_source_keV",
        "mean_edep_per_source_keV",
        "median_edep_per_source_keV",
        "mean_path_length_per_source_um",
        "mean_zns_path_length_per_source_um",
        "mean_bn_path_length_per_source_um",
        "mean_binder_path_length_per_source_um",
    ]
    summary_csv = output_dir / "stageD_zns_trajectory_summary.csv"
    write_csv(summary_csv, records, summary_fields)

    showcase_ratio, showcase_thickness, showcase_path, showcase_input_mode = choose_showcase_entry(
        entries,
        args.showcase_ratio,
        args.showcase_thickness,
    )
    trajectories = collect_showcase_trajectories(
        showcase_ratio,
        showcase_thickness,
        showcase_path,
        showcase_input_mode,
        args.reservoir_size,
        args.sample_trajectories,
        args.max_showcase_steps,
        args.seed,
    )

    plot_trajectory_showcase_3d(
        trajectories,
        output_dir / "stageD_zns_trajectory_showcase_3d.png",
        showcase_ratio,
        showcase_thickness,
    )
    plot_origin_aligned_trajectory_showcase_3d(
        trajectories,
        output_dir / "stageD_origin_aligned_trajectory_showcase_3d.png",
        showcase_ratio,
        showcase_thickness,
    )
    wrote_macro_showcase = plot_macro_trajectory_showcase_3d(
        trajectories,
        output_dir / "stageD_macro_trajectory_showcase_3d.png",
        showcase_ratio,
        showcase_thickness,
    )
    plot_trajectory_projections(
        trajectories,
        output_dir / "stageD_zns_trajectory_projection_panels.png",
        showcase_ratio,
        showcase_thickness,
    )
    plot_summary_lines(
        ratios,
        records,
        output_dir / "stageD_zns_trajectory_statistics_vs_thickness.png",
    )
    plot_1000um_mean_track_length(
        ratios,
        records,
        output_dir / "stageD_1000um_mean_total_track_length_by_ratio.png",
    )

    print(f"Wrote {summary_csv}")
    print(f"Wrote {output_dir / 'stageD_zns_trajectory_showcase_3d.png'}")
    print(f"Wrote {output_dir / 'stageD_origin_aligned_trajectory_showcase_3d.png'}")
    if wrote_macro_showcase:
        print(f"Wrote {output_dir / 'stageD_macro_trajectory_showcase_3d.png'}")
    print(f"Wrote {output_dir / 'stageD_zns_trajectory_projection_panels.png'}")
    print(f"Wrote {output_dir / 'stageD_zns_trajectory_statistics_vs_thickness.png'}")
    print(f"Wrote {output_dir / 'stageD_1000um_mean_total_track_length_by_ratio.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
