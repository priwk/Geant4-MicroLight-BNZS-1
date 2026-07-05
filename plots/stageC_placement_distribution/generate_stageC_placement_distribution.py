#!/usr/bin/env python3
"""Generate Stage C placement-model distribution figures."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

DEPENDENCY_ERROR = """
Failed to import the plotting dependencies.

Use a clean virtual environment from the project root:

    cd ~/g4work/B2
    python3 -m venv .venv
    . .venv/bin/activate
    python3 -m pip install --upgrade pip
    python3 -m pip install -r plots/stageC_placement_distribution/requirements.txt
    python3 plots/stageC_placement_distribution/generate_stageC_placement_distribution.py
""".strip()

try:
    import matplotlib

    matplotlib.use("Agg")

    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import font_manager
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

DESIRED_RATIO_ORDER = ["2-1", "1-1", "1-1.5", "1-2", "1-2.5", "1-3", "1-3.5", "1-4"]
RATIO_RE = re.compile(
    r"^(?P<bn>[+]?(?:\d+(?:\.\d*)?|\.\d+))-(?P<zns>[+]?(?:\d+(?:\.\d*)?|\.\d+))$"
)
META_RE = re.compile(r"^#\s*(?P<key>[A-Za-z0-9_]+)\s*=\s*(?P<value>.+?)\s*$")

PHASE_COLORS = {
    "BN": "#0b3c79",
    "ZnS": "#c0392b",
}
PHASE_MARKERS = {
    "BN": "o",
    "ZnS": "^",
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
class PlacementData:
    ratio: RatioKey
    path: Path
    metadata: dict[str, str]
    bn: np.ndarray
    zns: np.ndarray

    @property
    def all_centers(self) -> np.ndarray:
        if self.bn.size == 0:
            return self.zns
        if self.zns.size == 0:
            return self.bn
        return np.vstack([self.bn, self.zns])

    @property
    def patch_xy_um(self) -> float:
        return safe_float(self.metadata.get("patchXY_um"), default=extent_from_points(self.all_centers, 0))

    @property
    def local_thickness_um(self) -> float:
        return safe_float(
            self.metadata.get("localThickness_um"),
            default=extent_from_points(self.all_centers, 2),
        )

    @property
    def bn_radius_um(self) -> float:
        return safe_float(self.metadata.get("bnRadius_um"), default=math.nan)

    @property
    def zns_radius_um(self) -> float:
        return safe_float(self.metadata.get("znsRadius_um"), default=math.nan)

    @property
    def phi_achieved(self) -> float:
        return safe_float(self.metadata.get("phiAchieved"), default=math.nan)

    @property
    def placement_label(self) -> str:
        parent = self.path.parent.name
        return f"{parent}/{self.path.stem}"


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
            "axes.labelsize": 10,
            "axes.titlesize": 11,
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
        description="Plot representative Stage C BN/ZnS placement center distributions."
    )
    parser.add_argument(
        "--project-root",
        default=None,
        help="Project root. Defaults to the repository containing this script.",
    )
    parser.add_argument(
        "--placement-root",
        default=None,
        help="Placement root. Defaults to <project-root>/Input/placements.",
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
        help="Ratio folders to plot. Defaults to all ratio folders found under the placement root.",
    )
    parser.add_argument(
        "--placement-index",
        type=int,
        default=0,
        help="Sorted placement index to use per ratio. Default: 0.",
    )
    parser.add_argument(
        "--max-points-per-phase",
        type=int,
        default=2500,
        help="Maximum plotted centers per phase in each panel. Use 0 for all. Default: 2500.",
    )
    parser.add_argument(
        "--slice-axis",
        choices=["x", "y", "z"],
        default="z",
        help="Axis normal for the slice panels. Default: z.",
    )
    parser.add_argument(
        "--slice-center",
        type=float,
        default=0.0,
        help="Slice center coordinate in um. Default: 0.",
    )
    parser.add_argument(
        "--slice-half-width",
        type=float,
        default=1.0,
        help="Half-width of slice slab in um. Default: 1.",
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


def discover_ratio_tags(placement_root: Path) -> list[str]:
    ratio_tags: list[str] = []
    for ratio_dir in sorted(placement_root.iterdir()):
        if not ratio_dir.is_dir():
            continue
        try:
            parse_ratio_tag(ratio_dir.name)
        except ValueError:
            continue
        if any(ratio_dir.rglob("*.csv")):
            ratio_tags.append(ratio_dir.name)
    if not ratio_tags:
        raise SystemExit(f"No placement ratio folders with CSV files found under {placement_root}")
    return ratio_tags


def safe_float(value: str | None, *, default: float) -> float:
    if value is None or value == "":
        return default
    try:
        number = float(value)
    except ValueError:
        return default
    if not math.isfinite(number):
        return default
    return number


def extent_from_points(points: np.ndarray, axis: int) -> float:
    if points.size == 0:
        return math.nan
    return float(np.max(points[:, axis]) - np.min(points[:, axis]))


def ratio_to_dir_name(path: Path, placement_root: Path) -> str:
    rel = path.resolve().relative_to(placement_root.resolve())
    return rel.parts[0]


def find_placement_for_ratio(placement_root: Path, ratio_tag: str, placement_index: int) -> Path:
    ratio_dir = placement_root / ratio_tag
    if not ratio_dir.is_dir():
        raise SystemExit(f"Placement ratio directory not found: {ratio_dir}")
    candidates = sorted(ratio_dir.rglob("*.csv"))
    if not candidates:
        raise SystemExit(f"No placement CSV files found under {ratio_dir}")
    if placement_index >= len(candidates):
        raise SystemExit(
            f"--placement-index {placement_index} is out of range for {ratio_tag}; "
            f"only {len(candidates)} placement files are available"
        )
    return candidates[placement_index]


def read_placement_csv(path: Path, placement_root: Path) -> PlacementData:
    ratio = parse_ratio_tag(ratio_to_dir_name(path, placement_root))
    metadata: dict[str, str] = {}
    bn_points: list[tuple[float, float, float]] = []
    zns_points: list[tuple[float, float, float]] = []

    with path.open(newline="", encoding="utf-8-sig") as handle:
        content_lines = []
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                match = META_RE.match(line)
                if match is not None:
                    metadata[match.group("key")] = match.group("value")
                continue
            content_lines.append(line)

    if not content_lines:
        raise SystemExit(f"Placement CSV has no center rows: {path}")

    reader = csv.DictReader(content_lines)
    required = {"type", "x_um", "y_um", "z_um"}
    missing = required.difference(reader.fieldnames or [])
    if missing:
        raise SystemExit(f"{path} is missing columns: {', '.join(sorted(missing))}")

    for row in reader:
        phase = row["type"].strip()
        point = (float(row["x_um"]), float(row["y_um"]), float(row["z_um"]))
        if phase == "BN":
            bn_points.append(point)
        elif phase == "ZnS":
            zns_points.append(point)
        else:
            raise SystemExit(f"{path} contains unsupported phase type: {phase!r}")

    bn = np.asarray(bn_points, dtype=float).reshape((-1, 3))
    zns = np.asarray(zns_points, dtype=float).reshape((-1, 3))
    return PlacementData(ratio=ratio, path=path, metadata=metadata, bn=bn, zns=zns)


def downsample(points: np.ndarray, max_points: int, seed: int) -> np.ndarray:
    if max_points <= 0 or points.shape[0] <= max_points:
        return points
    rng = np.random.default_rng(seed)
    indices = rng.choice(points.shape[0], size=max_points, replace=False)
    indices.sort()
    return points[indices, :]


def style_2d_axes(ax) -> None:
    ax.tick_params(direction="in", which="both", top=True, right=True, width=1.1)
    for spine in ax.spines.values():
        spine.set_linewidth(1.3)


def set_equal_3d_box(ax, data: PlacementData) -> None:
    half_xy = 0.5 * data.patch_xy_um
    half_z = 0.5 * data.local_thickness_um
    ax.set_xlim(-half_xy, half_xy)
    ax.set_ylim(-half_xy, half_xy)
    ax.set_zlim(-half_z, half_z)
    try:
        ax.set_box_aspect((data.patch_xy_um, data.patch_xy_um, data.local_thickness_um))
    except AttributeError:
        pass


def box_edges(data: PlacementData) -> list[list[tuple[float, float, float]]]:
    half_xy = 0.5 * data.patch_xy_um
    half_z = 0.5 * data.local_thickness_um
    corners = [
        (x, y, z)
        for x in (-half_xy, half_xy)
        for y in (-half_xy, half_xy)
        for z in (-half_z, half_z)
    ]
    corner_set = set(corners)
    edges = []
    for c0 in corners:
        for c1 in corner_set:
            diff = sum(1 for a, b in zip(c0, c1) if not math.isclose(a, b))
            if diff == 1 and c0 < c1:
                edges.append([c0, c1])
    return edges


def panel_grid(n_panels: int) -> tuple[int, int]:
    if n_panels <= 0:
        return 1, 1
    cols = min(4, max(1, math.ceil(math.sqrt(n_panels))))
    rows = math.ceil(n_panels / cols)
    return rows, cols


def plot_3d_panels(placements: list[PlacementData], output_path: Path, max_points: int) -> None:
    rows, cols = panel_grid(len(placements))
    fig = plt.figure(figsize=(4.2 * cols, 3.9 * rows + 0.45))
    axes = [fig.add_subplot(rows, cols, idx + 1, projection="3d") for idx in range(len(placements))]

    for idx, (ax, data) in enumerate(zip(axes, placements)):
        bn = downsample(data.bn, max_points, seed=1000 + idx)
        zns = downsample(data.zns, max_points, seed=2000 + idx)

        ax.scatter(
            bn[:, 0],
            bn[:, 1],
            bn[:, 2],
            s=1.4,
            c=PHASE_COLORS["BN"],
            marker=PHASE_MARKERS["BN"],
            alpha=0.36,
            linewidths=0,
            depthshade=False,
        )
        ax.scatter(
            zns[:, 0],
            zns[:, 1],
            zns[:, 2],
            s=3.2,
            c=PHASE_COLORS["ZnS"],
            marker=PHASE_MARKERS["ZnS"],
            alpha=0.62,
            linewidths=0,
            depthshade=False,
        )

        ax.add_collection3d(Line3DCollection(box_edges(data), colors="#6b7280", linewidths=0.55, alpha=0.45))
        ax.view_init(elev=20, azim=-48)
        set_equal_3d_box(ax, data)
        ax.set_xlabel(f"x ({mu_m_text()})", labelpad=2)
        ax.set_ylabel(f"y ({mu_m_text()})", labelpad=2)
        ax.set_zlabel(f"z ({mu_m_text()})", labelpad=2)
        ax.set_title(
            f"BN:ZnS = {data.ratio.display_tag}\n"
            f"BN {data.bn.shape[0]:,}, ZnS {data.zns.shape[0]:,}",
            pad=0,
        )
        ax.grid(False)
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        ax.tick_params(pad=0, labelsize=8)

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker=PHASE_MARKERS["BN"],
            color="none",
            label="BN centers",
            markerfacecolor=PHASE_COLORS["BN"],
            markeredgewidth=0,
            markersize=6,
            alpha=0.8,
        ),
        Line2D(
            [0],
            [0],
            marker=PHASE_MARKERS["ZnS"],
            color="none",
            label="ZnS(Ag) centers",
            markerfacecolor=PHASE_COLORS["ZnS"],
            markeredgewidth=0,
            markersize=7,
            alpha=0.85,
        ),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 0.012),
    )
    fig.suptitle("Stage C representative placement center distributions", y=0.985, fontsize=14)
    fig.subplots_adjust(left=0.02, right=0.99, top=0.91, bottom=0.07, wspace=0.03, hspace=0.16)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def slice_mask(points: np.ndarray, axis_idx: int, center: float, half_width: float) -> np.ndarray:
    if points.size == 0:
        return np.zeros(0, dtype=bool)
    return np.abs(points[:, axis_idx] - center) <= half_width


def plot_slice_panels(
    placements: list[PlacementData],
    output_path: Path,
    max_points: int,
    axis: str,
    center: float,
    half_width: float,
) -> None:
    axis_idx = {"x": 0, "y": 1, "z": 2}[axis]
    plane_axes = [idx for idx in range(3) if idx != axis_idx]
    labels = ["x", "y", "z"]

    rows, cols = panel_grid(len(placements))
    fig, axes = plt.subplots(rows, cols, figsize=(4.15 * cols, 3.55 * rows + 0.45), constrained_layout=False)
    axes_arr = np.asarray(axes, dtype=object).reshape(-1)

    for idx, (ax, data) in enumerate(zip(axes_arr, placements)):
        bn_slice = data.bn[slice_mask(data.bn, axis_idx, center, half_width), :]
        zns_slice = data.zns[slice_mask(data.zns, axis_idx, center, half_width), :]
        bn_plot = downsample(bn_slice, max_points, seed=3000 + idx)
        zns_plot = downsample(zns_slice, max_points, seed=4000 + idx)

        if bn_plot.size:
            ax.scatter(
                bn_plot[:, plane_axes[0]],
                bn_plot[:, plane_axes[1]],
                s=2.0,
                c=PHASE_COLORS["BN"],
                marker=PHASE_MARKERS["BN"],
                alpha=0.42,
                linewidths=0,
            )
        if zns_plot.size:
            ax.scatter(
                zns_plot[:, plane_axes[0]],
                zns_plot[:, plane_axes[1]],
                s=5.0,
                c=PHASE_COLORS["ZnS"],
                marker=PHASE_MARKERS["ZnS"],
                alpha=0.68,
                linewidths=0,
            )

        half_xy = 0.5 * data.patch_xy_um
        half_z = 0.5 * data.local_thickness_um
        extents = [(-half_xy, half_xy), (-half_xy, half_xy), (-half_z, half_z)]
        ax.set_xlim(*extents[plane_axes[0]])
        ax.set_ylim(*extents[plane_axes[1]])
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(f"{labels[plane_axes[0]]} ({mu_m_text()})")
        ax.set_ylabel(f"{labels[plane_axes[1]]} ({mu_m_text()})")
        ax.set_title(
            f"BN:ZnS = {data.ratio.display_tag}\n"
            f"slice BN {bn_slice.shape[0]:,}, ZnS {zns_slice.shape[0]:,}"
        )
        style_2d_axes(ax)

    for ax in axes_arr[len(placements) :]:
        ax.set_visible(False)

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker=PHASE_MARKERS["BN"],
            color="none",
            label="BN centers",
            markerfacecolor=PHASE_COLORS["BN"],
            markeredgewidth=0,
            markersize=6,
            alpha=0.8,
        ),
        Line2D(
            [0],
            [0],
            marker=PHASE_MARKERS["ZnS"],
            color="none",
            label="ZnS(Ag) centers",
            markerfacecolor=PHASE_COLORS["ZnS"],
            markeredgewidth=0,
            markersize=7,
            alpha=0.85,
        ),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 0.01),
    )
    fig.suptitle(
        f"Stage C placement center slices: {axis} = {center:g} +/- {half_width:g} {mu_m_text()}",
        y=0.985,
        fontsize=14,
    )
    fig.subplots_adjust(left=0.06, right=0.99, top=0.90, bottom=0.09, wspace=0.24, hspace=0.34)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def mu_m_text() -> str:
    return "μm"


def fmt(value: object) -> object:
    if isinstance(value, float):
        return f"{value:.12g}"
    return value


def write_summary_csv(path: Path, placements: list[PlacementData]) -> None:
    fields = [
        "ratio_tag",
        "bn_wt",
        "zns_wt",
        "placement_file",
        "bn_centers",
        "zns_centers",
        "patchXY_um",
        "localThickness_um",
        "bnRadius_um",
        "znsRadius_um",
        "phiAchieved",
        "actualBNToZnSMass",
        "actualZnSToBNMass",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for data in placements:
            row = {
                "ratio_tag": data.ratio.tag,
                "bn_wt": data.ratio.bn_wt,
                "zns_wt": data.ratio.zns_wt,
                "placement_file": data.path.as_posix(),
                "bn_centers": data.bn.shape[0],
                "zns_centers": data.zns.shape[0],
                "patchXY_um": data.patch_xy_um,
                "localThickness_um": data.local_thickness_um,
                "bnRadius_um": data.bn_radius_um,
                "znsRadius_um": data.zns_radius_um,
                "phiAchieved": data.phi_achieved,
                "actualBNToZnSMass": safe_float(data.metadata.get("actualBNToZnSMass"), default=math.nan),
                "actualZnSToBNMass": safe_float(data.metadata.get("actualZnSToBNMass"), default=math.nan),
            }
            writer.writerow({field: fmt(row[field]) for field in fields})


def main() -> None:
    args = parse_args()
    if args.placement_index < 0:
        raise SystemExit("--placement-index must be >= 0")
    if args.max_points_per_phase < 0:
        raise SystemExit("--max-points-per-phase must be >= 0")
    if args.slice_half_width <= 0.0:
        raise SystemExit("--slice-half-width must be > 0")

    project_root = Path(args.project_root).resolve() if args.project_root else SCRIPT_DIR.parents[1]
    placement_root = (
        Path(args.placement_root).resolve()
        if args.placement_root
        else project_root / "Input" / "placements"
    )
    output_dir = Path(args.output_dir).resolve() if args.output_dir else SCRIPT_DIR / "output"

    if not placement_root.is_dir():
        raise SystemExit(f"Placement root not found: {placement_root}")

    ratio_tags = args.ratios if args.ratios is not None else discover_ratio_tags(placement_root)
    requested_ratios = [parse_ratio_tag(tag) for tag in ratio_tags]
    requested_ratios = sorted(requested_ratios, key=ratio_display_sort_key)

    placements: list[PlacementData] = []
    for ratio in requested_ratios:
        placement_path = find_placement_for_ratio(placement_root, ratio.tag, args.placement_index)
        data = read_placement_csv(placement_path, placement_root)
        placements.append(data)
        print(
            f"Loaded {ratio.tag}: {placement_path.relative_to(project_root)} "
            f"({data.bn.shape[0]:,} BN, {data.zns.shape[0]:,} ZnS centers)"
        )

    if not placements:
        raise SystemExit("No placements loaded.")

    summary_csv = output_dir / "stageC_placement_distribution_summary.csv"
    figure_3d = output_dir / "stageC_placement_centers_3d_panels.png"
    figure_slice = output_dir / "stageC_placement_center_slice_panels.png"

    write_summary_csv(summary_csv, placements)
    plot_3d_panels(placements, figure_3d, args.max_points_per_phase)
    plot_slice_panels(
        placements,
        figure_slice,
        args.max_points_per_phase,
        args.slice_axis,
        args.slice_center,
        args.slice_half_width,
    )

    print(f"Wrote {summary_csv}")
    print(f"Wrote {figure_3d}")
    print(f"Wrote {figure_slice}")


if __name__ == "__main__":
    main()
