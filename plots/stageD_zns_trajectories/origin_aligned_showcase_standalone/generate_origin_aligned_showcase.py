#!/usr/bin/env python3
"""Build the standalone origin-aligned Stage D trajectory showcase."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
STAGED_DIR = SCRIPT_DIR.parent
os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR / ".matplotlib"))

DEPENDENCY_ERROR = """
Failed to import plotting dependencies.

Install the required packages from this folder:

    python3 -m pip install -r requirements.txt
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


FONT_DIRS = [SCRIPT_DIR / "fonts", STAGED_DIR / "fonts"]
DEFAULT_PROJECT_ROOT = SCRIPT_DIR.parents[2]
DESIRED_RATIO_ORDER = ["2-1", "1-1", "1-1.5", "1-2", "1-2.5", "1-3"]
RATIO_RE = re.compile(
    r"^(?P<bn>[+]?(?:\d+(?:\.\d*)?|\.\d+))-(?P<zns>[+]?(?:\d+(?:\.\d*)?|\.\d+))$"
)
FULL_TRACK_RE = re.compile(r"^(?P<thickness>[+]?(?:\d+(?:\.\d*)?|\.\d+))_alpha_li_steps\.csv$")
SLIM_TRACK_RE = re.compile(r"^(?P<thickness>[+]?(?:\d+(?:\.\d*)?|\.\d+))_zns_track_steps\.csv$")

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

ALIGNED_FIELDS = [
    "trajectory_id",
    "particle",
    "phase_pre",
    "phase_post",
    "rel_x_pre_um",
    "rel_y_pre_um",
    "rel_z_pre_um",
    "rel_x_post_um",
    "rel_y_post_um",
    "rel_z_post_um",
    "step_len_um",
    "edep_keV",
]


@dataclass(frozen=True)
class RatioKey:
    tag: str
    bn_wt: float
    zns_wt: float

    @property
    def display_tag(self) -> str:
        return self.tag.replace("-", ":")


@dataclass
class AlignedStep:
    trajectory_id: str
    particle: str
    phase_pre: str
    phase_post: str
    start: tuple[float, float, float]
    end: tuple[float, float, float]
    step_len_um: float
    edep_kev: float


@dataclass
class AlignedTrajectory:
    trajectory_id: str
    particle: str
    steps: list[AlignedStep] = field(default_factory=list)


def configure_fonts() -> None:
    for font_dir in FONT_DIRS:
        for path in [
            font_dir / "Arial.ttf",
            font_dir / "Arial-Bold.ttf",
            font_dir / "Arial-Italic.ttf",
            font_dir / "Arial-BoldItalic.ttf",
        ]:
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot stageD_origin_aligned_trajectory_showcase_3d.png from the standalone "
            "origin-aligned CSV in this folder. Use --extract-from-stageb on the server "
            "to rebuild that CSV from Stage B raw tracks."
        )
    )
    parser.add_argument(
        "--extract-from-stageb",
        action="store_true",
        help="Server-side mode: rebuild the standalone aligned CSV from Stage B raw track files before plotting.",
    )
    parser.add_argument(
        "--project-root",
        default=None,
        help=(
            "Server-side extraction root. Defaults to the repository root inferred from this script path. "
            "Use . only when running from the project root and overriding this behavior intentionally."
        ),
    )
    parser.add_argument(
        "--stageb-root",
        default="Output/stageB",
        help="Stage B output root relative to project root unless absolute. Default: Output/stageB.",
    )
    parser.add_argument(
        "--ratios",
        nargs="+",
        default=["1-2"],
        help="Ratio folders to process. Default: 1-2.",
    )
    parser.add_argument(
        "--thicknesses",
        nargs="+",
        type=float,
        default=[1000.0],
        help="Optional thicknesses in um to process. Default: 1000.",
    )
    parser.add_argument(
        "--showcase-ratio",
        default="1-2",
        help="Ratio used for the origin-aligned showcase. Default: 1-2.",
    )
    parser.add_argument(
        "--showcase-thickness",
        type=float,
        default=1000.0,
        help="Thickness used for the origin-aligned showcase. Default: 1000.",
    )
    parser.add_argument(
        "--ratio-label",
        default=None,
        help="Optional display label for BN:ZnS ratio. Defaults to the selected ratio folder label.",
    )
    parser.add_argument(
        "--row-limit",
        type=int,
        default=10000,
        help="Number of aligned step rows written to the standalone CSV. Default: 10000.",
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=0,
        help="Process only the first N track files after sorting. Use 0 for all. Default: 0.",
    )
    parser.add_argument(
        "--aligned-csv",
        default="origin_aligned_particle_trajectories_first10000.csv",
        help=(
            "Standalone aligned CSV path, relative to this script folder unless absolute. "
            "Default: origin_aligned_particle_trajectories_first10000.csv."
        ),
    )
    parser.add_argument(
        "--figure",
        default="stageD_origin_aligned_trajectory_showcase_3d.png",
        help=(
            "Output PNG path, relative to this script folder unless absolute. "
            "Default: stageD_origin_aligned_trajectory_showcase_3d.png."
        ),
    )
    return parser.parse_args([arg for arg in sys.argv[1:] if arg.strip()])


def resolve_relative_path(path_text: str, base_dir: Path) -> Path:
    path = Path(path_text)
    if path.is_absolute():
        return path
    return base_dir / path


def display_path(path: Path, base_dir: Path) -> str:
    try:
        return path.resolve().relative_to(base_dir.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def parse_ratio_tag(tag: str) -> RatioKey:
    match = RATIO_RE.fullmatch(tag)
    if match is None:
        raise ValueError(f"Unsupported ratio folder name: {tag}")
    return RatioKey(tag=tag, bn_wt=float(match.group("bn")), zns_wt=float(match.group("zns")))


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

        thicknesses = sorted(set(full_by_thickness).union(slim_by_thickness), key=thickness_sort_key)
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


def source_point(row: dict[str, str]) -> tuple[float, float, float]:
    prefix = "unwrapped_" if "unwrapped_x_pre_um" in row else ""
    return (
        safe_float(row.get(f"{prefix}x_pre_um")),
        safe_float(row.get(f"{prefix}y_pre_um")),
        safe_float(row.get(f"{prefix}z_pre_um")),
    )


def end_point(row: dict[str, str]) -> tuple[float, float, float]:
    prefix = "unwrapped_" if "unwrapped_x_post_um" in row else ""
    return (
        safe_float(row.get(f"{prefix}x_post_um")),
        safe_float(row.get(f"{prefix}y_post_um")),
        safe_float(row.get(f"{prefix}z_post_um")),
    )


def trajectory_key(row: dict[str, str]) -> str:
    source_uid = row.get("source_event_uid", "").strip()
    particle = row.get("particle", "").strip()
    track_id = row.get("trackID", "").strip()
    return f"{source_uid}|{particle}|{track_id}"


def required_source_columns() -> set[str]:
    return {
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
    }


def write_origin_aligned_csv(source_path: Path, output_path: Path, row_limit: int) -> int:
    if row_limit <= 0:
        raise SystemExit("--row-limit must be > 0")

    origins: dict[str, tuple[float, float, float]] = {}
    trajectory_ids: dict[str, str] = {}
    rows_written = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with source_path.open(newline="", encoding="utf-8-sig") as source_handle:
        reader = csv.DictReader(source_handle)
        missing = required_source_columns().difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"{source_path} is missing columns: {', '.join(sorted(missing))}")

        with output_path.open("w", newline="", encoding="utf-8") as output_handle:
            writer = csv.DictWriter(output_handle, fieldnames=ALIGNED_FIELDS)
            writer.writeheader()

            for row in reader:
                if row.get("particle", "").strip() not in {"alpha", "Li7"}:
                    continue
                key = trajectory_key(row)
                if key not in origins:
                    origins[key] = source_point(row)
                    trajectory_ids[key] = f"T{len(trajectory_ids) + 1:06d}"
                origin = np.asarray(origins[key], dtype=float)
                start = np.asarray(source_point(row), dtype=float) - origin
                end = np.asarray(end_point(row), dtype=float) - origin
                writer.writerow(
                    {
                        "trajectory_id": trajectory_ids[key],
                        "particle": row.get("particle", "").strip(),
                        "phase_pre": row.get("phase_pre", "unknown").strip() or "unknown",
                        "phase_post": row.get("phase_post", "unknown").strip() or "unknown",
                        "rel_x_pre_um": f"{start[0]:.8g}",
                        "rel_y_pre_um": f"{start[1]:.8g}",
                        "rel_z_pre_um": f"{start[2]:.8g}",
                        "rel_x_post_um": f"{end[0]:.8g}",
                        "rel_y_post_um": f"{end[1]:.8g}",
                        "rel_z_post_um": f"{end[2]:.8g}",
                        "step_len_um": f"{safe_float(row.get('step_len_um')):.8g}",
                        "edep_keV": f"{safe_float(row.get('edep_keV')):.8g}",
                    }
                )
                rows_written += 1
                if rows_written >= row_limit:
                    break

    if rows_written == 0:
        raise SystemExit(f"No trajectory rows written from {source_path}")
    return rows_written


def read_aligned_csv(path: Path) -> list[AlignedTrajectory]:
    trajectories: dict[str, AlignedTrajectory] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        missing = set(ALIGNED_FIELDS).difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"{path} is missing columns: {', '.join(sorted(missing))}")

        for row in reader:
            trajectory_id = row["trajectory_id"].strip()
            particle = row["particle"].strip()
            traj = trajectories.get(trajectory_id)
            if traj is None:
                traj = AlignedTrajectory(trajectory_id=trajectory_id, particle=particle)
                trajectories[trajectory_id] = traj

            traj.steps.append(
                AlignedStep(
                    trajectory_id=trajectory_id,
                    particle=particle,
                    phase_pre=row.get("phase_pre", "unknown").strip() or "unknown",
                    phase_post=row.get("phase_post", "unknown").strip() or "unknown",
                    start=(
                        safe_float(row.get("rel_x_pre_um")),
                        safe_float(row.get("rel_y_pre_um")),
                        safe_float(row.get("rel_z_pre_um")),
                    ),
                    end=(
                        safe_float(row.get("rel_x_post_um")),
                        safe_float(row.get("rel_y_post_um")),
                        safe_float(row.get("rel_z_post_um")),
                    ),
                    step_len_um=safe_float(row.get("step_len_um")),
                    edep_kev=safe_float(row.get("edep_keV")),
                )
            )

    result = [traj for traj in trajectories.values() if traj.steps]
    if not result:
        raise SystemExit(f"No aligned trajectories found in {path}")
    return result


def aligned_segments(
    trajectories: list[AlignedTrajectory],
) -> tuple[np.ndarray, list[str], list[str], np.ndarray]:
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


def aligned_bounds(segments: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if segments.size == 0:
        return np.asarray([-1.0, -1.0, -1.0]), np.asarray([1.0, 1.0, 1.0])
    points = segments.reshape((-1, 3))
    max_abs = float(np.max(np.abs(points)))
    half = max(max_abs * 1.12, 1.0)
    return np.asarray([-half, -half, -half]), np.asarray([half, half, half])


def plot_origin_aligned_showcase(
    trajectories: list[AlignedTrajectory],
    output_path: Path,
    ratio_label: str,
    thickness_um: float,
) -> None:
    segments, particles, phases, edeps = aligned_segments(trajectories)
    mins, maxs = aligned_bounds(segments)
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
        f"BN:ZnS = {ratio_label}, t = {thickness_um:g} {mu_m_text()}",
        y=0.98,
        fontsize=14,
    )

    phases_present = set(phases)
    phase_handles = [
        Line2D([0], [0], color=PHASE_COLORS[phase], lw=2.4, label=PHASE_LABELS.get(phase, phase))
        for phase in PHASE_ORDER
        if phase in phases_present
    ]
    if phase_handles:
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


def main() -> int:
    configure_fonts()
    args = parse_args()
    if args.row_limit <= 0:
        raise SystemExit("--row-limit must be > 0")
    if args.max_files < 0:
        raise SystemExit("--max-files must be >= 0")

    project_root = Path(args.project_root).resolve() if args.project_root else DEFAULT_PROJECT_ROOT
    aligned_csv = resolve_relative_path(args.aligned_csv, SCRIPT_DIR)
    figure = resolve_relative_path(args.figure, SCRIPT_DIR)
    ratio_label = args.ratio_label or args.showcase_ratio.replace("-", ":")
    showcase_thickness = args.showcase_thickness

    if args.extract_from_stageb:
        stageb_root = resolve_relative_path(args.stageb_root, project_root)
        showcase_ratio_arg = args.showcase_ratio
        if showcase_ratio_arg is None and args.ratios and len(args.ratios) == 1:
            showcase_ratio_arg = args.ratios[0]

        _ratios, entries = discover_track_files(stageb_root, args.ratios, args.thicknesses, args.max_files)
        showcase_ratio, showcase_thickness, source_path, input_mode = choose_showcase_entry(
            entries,
            showcase_ratio_arg,
            args.showcase_thickness,
        )
        ratio_label = args.ratio_label or showcase_ratio.display_tag
        rows_written = write_origin_aligned_csv(source_path, aligned_csv, args.row_limit)
        print(f"Source: {display_path(source_path, project_root)} ({input_mode})")
        print(f"Selected: BN:ZnS={showcase_ratio.display_tag}, t={showcase_thickness:g} um")
        print(f"Wrote {rows_written:,} aligned step rows to {display_path(aligned_csv, SCRIPT_DIR)}")
    elif not aligned_csv.is_file():
        raise SystemExit(
            f"Standalone aligned CSV not found: {display_path(aligned_csv, SCRIPT_DIR)}\n"
            "Run this script in the folder containing the CSV, or run on the server with --extract-from-stageb first."
        )

    trajectories = read_aligned_csv(aligned_csv)
    plot_origin_aligned_showcase(trajectories, figure, ratio_label, showcase_thickness)

    particle_counts: defaultdict[str, int] = defaultdict(int)
    for traj in trajectories:
        particle_counts[traj.particle] += len(traj.steps)

    print(
        "Aligned rows by particle: "
        + ", ".join(f"{particle}={count:,}" for particle, count in sorted(particle_counts.items()))
    )
    print(f"Read {display_path(aligned_csv, SCRIPT_DIR)}")
    print(f"Wrote {display_path(figure, SCRIPT_DIR)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
