#!/usr/bin/env python3

import argparse
import csv
import json
import math
from pathlib import Path


def mean(values):
    return sum(values) / len(values) if values else 0.0


def stddev(values):
    if len(values) < 2:
        return 0.0
    mu = mean(values)
    return math.sqrt(sum((v - mu) ** 2 for v in values) / (len(values) - 1))


def load_rows(base_dir: Path, ratio: str):
    root = base_dir / "Output" / "stageD_optical_homogenization" / ratio
    rows = []
    summary_paths = sorted(root.rglob("stageD_summary.csv"))
    if not summary_paths:
        summary_paths = sorted(root.rglob("optical_homogenization_summary.csv"))
    for path in summary_paths:
        with path.open(newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                row["_path"] = str(path)
                rows.append(row)
    return rows


def load_phase_function_rows(summary_rows):
    grouped = {}
    for row in summary_rows:
        summary_path = Path(row["_path"])
        phase_name = row.get("phase_function_file", "phase_function.csv")
        phase_path = summary_path.parent / phase_name
        if not phase_path.is_file():
            continue
        with phase_path.open(newline="") as f:
            reader = csv.DictReader(f)
            for pf_row in reader:
                key = int(pf_row["bin_id"])
                grouped.setdefault(key, {
                    "lambda_nm": float(pf_row["lambda_nm"]),
                    "cos_theta_min": float(pf_row["cos_theta_min"]),
                    "cos_theta_max": float(pf_row["cos_theta_max"]),
                    "counts": [],
                    "probabilities": [],
                    "densities": [],
                })
                grouped[key]["counts"].append(float(pf_row["count"]))
                grouped[key]["probabilities"].append(float(pf_row["probability"]))
                grouped[key]["densities"].append(float(pf_row["probability_density"]))
    return grouped


def numeric_series(rows, key):
    return [float(r[key]) for r in rows]


def primary_key(rows, preferred, fallback):
    if preferred in rows[0]:
        return preferred
    return fallback


def optional_numeric_series(rows, preferred_key, fallback_key):
    values = []
    for row in rows:
        key = preferred_key if preferred_key in row else fallback_key
        values.append(float(row[key]))
    return values


def row_config_signature(row):
    return (
        row.get("source_mode", ""),
        row.get("boundary_mode", ""),
        row.get("reentry_mode", ""),
        row.get("particle_reentry_mode", ""),
        row.get("matrix_reentry_mode", ""),
        row.get("scatter_metric", ""),
        row.get("target_primary_scatter", ""),
        row.get("theta_threshold_deg", ""),
        row.get("matrix_n", ""),
        row.get("matrix_abs_um", ""),
        row.get("bn_n", ""),
        row.get("bn_abs_um", ""),
        row.get("zns_n", ""),
        row.get("zns_abs_um", ""),
        row.get("wavelength_nm", ""),
    )


def write_scalar_summary(out_path: Path, args, wavelengths, photons,
                         mu_a_count, mu_a_expected, mu_s_primary, g_primary, g2_primary,
                         mu_sp_primary, mu_s_total, g_total, mu_sp_total, mu_s_particle,
                         g_particle, mu_sp_particle, mu_s_boundary_primary,
                         g_boundary_primary, mu_sp_boundary_primary, mu_s_boundary,
                         g_boundary, mu_sp_boundary, absorbed, path_length, reentry):
    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "ratio",
                "wavelength_nm",
                "n_placements",
                "n_photons_total",
                "mu_a_mean_per_um",
                "mu_a_std_per_um",
                "mu_a_expected_mean_per_um",
                "mu_a_expected_std_per_um",
                "mu_s_total_mean_per_um",
                "mu_s_total_std_per_um",
                "g_total_mean",
                "g_total_std",
                "mu_s_prime_total_mean_per_um",
                "mu_s_prime_total_std_per_um",
                "mu_s_mean_per_um",
                "mu_s_std_per_um",
                "g_mean",
                "g_std",
                "g2_mean",
                "g2_std",
                "mu_s_prime_mean_per_um",
                "mu_s_prime_std_per_um",
                "mu_s_particle_mean_per_um",
                "mu_s_particle_std_per_um",
                "g_particle_mean",
                "g_particle_std",
                "mu_s_prime_particle_mean_per_um",
                "mu_s_prime_particle_std_per_um",
                "mu_s_boundary_primary_mean_per_um",
                "mu_s_boundary_primary_std_per_um",
                "g_boundary_primary_mean",
                "g_boundary_primary_std",
                "mu_s_prime_boundary_primary_mean_per_um",
                "mu_s_prime_boundary_primary_std_per_um",
                "mu_s_boundary_mean_per_um",
                "mu_s_boundary_std_per_um",
                "g_boundary_mean",
                "g_boundary_std",
                "mu_s_prime_boundary_mean_per_um",
                "mu_s_prime_boundary_std_per_um",
                "mean_absorbed_fraction",
                "mean_path_length_um",
                "mean_reentry_count",
                "notes",
            ]
        )
        writer.writerow(
            [
                args.ratio,
                mean(wavelengths),
                len(wavelengths),
                sum(photons),
                mean(mu_a_count),
                stddev(mu_a_count),
                mean(mu_a_expected),
                stddev(mu_a_expected),
                mean(mu_s_total),
                stddev(mu_s_total),
                mean(g_total),
                stddev(g_total),
                mean(mu_sp_total),
                stddev(mu_sp_total),
                mean(mu_s_primary),
                stddev(mu_s_primary),
                mean(g_primary),
                stddev(g_primary),
                mean(g2_primary),
                stddev(g2_primary),
                mean(mu_sp_primary),
                stddev(mu_sp_primary),
                mean(mu_s_particle),
                stddev(mu_s_particle),
                mean(g_particle),
                stddev(g_particle),
                mean(mu_sp_particle),
                stddev(mu_sp_particle),
                mean(mu_s_boundary_primary),
                stddev(mu_s_boundary_primary),
                mean(g_boundary_primary),
                stddev(g_boundary_primary),
                mean(mu_sp_boundary_primary),
                stddev(mu_sp_boundary_primary),
                mean(mu_s_boundary),
                stddev(mu_s_boundary),
                mean(g_boundary),
                stddev(g_boundary),
                mean(mu_sp_boundary),
                stddev(mu_sp_boundary),
                mean(absorbed),
                mean(path_length),
                mean(reentry),
                (
                    "Ensemble average over StageD placements; one Geant4 run equals one placement. "
                    "Primary mu_s/g/g2/mu_s_prime columns follow the StageD primary scatter metric; "
                    "raw columns remain particle encounter without threshold for comparison. "
                    "Use phase_function_mean_by_ratio.csv for tabulated scattering instead of HG when possible."
                ),
            ]
        )


def write_phase_function_summary(out_path: Path, grouped):
    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "lambda_nm",
                "bin_id",
                "cos_theta_min",
                "cos_theta_max",
                "count_mean",
                "count_std",
                "probability_mean",
                "probability_std",
                "probability_density_mean",
                "probability_density_std",
            ]
        )
        for bin_id in sorted(grouped):
            entry = grouped[bin_id]
            writer.writerow(
                [
                    entry["lambda_nm"],
                    bin_id,
                    entry["cos_theta_min"],
                    entry["cos_theta_max"],
                    mean(entry["counts"]),
                    stddev(entry["counts"]),
                    mean(entry["probabilities"]),
                    stddev(entry["probabilities"]),
                    mean(entry["densities"]),
                    stddev(entry["densities"]),
                ]
            )


def write_mc_input_csv(out_path: Path, ratio: str, wavelength_nm: float,
                       mu_a_expected_mean: float, mu_s_mean: float, g_mean: float,
                       g2_mean: float, mu_sp_mean: float,
                       phase_function_file: str):
    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "ratio",
                "wavelength_nm",
                "param_model",
                "recommended_mu_a_per_um",
                "recommended_mu_s_per_um",
                "recommended_g1",
                "recommended_g2",
                "recommended_mu_s_prime_per_um",
                "phase_function_file",
                "recommended_scattering_model",
                "notes",
            ]
        )
        writer.writerow(
            [
                ratio,
                wavelength_nm,
                "GO_RVE",
                mu_a_expected_mean,
                mu_s_mean,
                g_mean,
                g2_mean,
                mu_sp_mean,
                phase_function_file,
                "tabulated_phase_function",
                "Use the tabulated phase function when the downstream MC supports it; otherwise HG with g1 is the fallback approximation.",
            ]
        )


def write_mc_input_json(out_path: Path, ratio: str, wavelength_nm: float,
                        mu_a_expected_mean: float, mu_s_mean: float, g_mean: float,
                        g2_mean: float, mu_sp_mean: float, phase_function_file: str):
    payload = {
        "ratio": ratio,
        "wavelength_nm": wavelength_nm,
        "param_model": "GO_RVE",
        "recommended_inputs": {
            "mu_a_per_um": mu_a_expected_mean,
            "mu_s_per_um": mu_s_mean,
            "g1": g_mean,
            "g2": g2_mean,
            "mu_s_prime_per_um": mu_sp_mean,
            "phase_function_file": phase_function_file,
            "scattering_model": "tabulated_phase_function",
        },
        "fallback_inputs": {
            "scattering_model": "henyey_greenstein",
            "g1": g_mean,
        },
        "notes": [
            "Primary statistics follow the StageD primary scatter metric.",
            "Raw encounter columns remain available without threshold for comparison.",
            "mu_a_per_um is the expected absorption estimate from phase path lengths and ABSLENGTH.",
            "Prefer the tabulated phase function over HG when the downstream Monte Carlo supports it.",
        ],
    }
    out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(
        description="Merge StageD raw optical parameters over placements for one ratio."
    )
    parser.add_argument("--ratio", required=True, help="Ratio tag such as 1-2 or 2-1.")
    parser.add_argument(
        "--project-root",
        default=".",
        help="Project root containing Output/.",
    )
    parser.add_argument(
        "--scatter-metric",
        default=None,
        help="Only merge rows with this scatter_metric.",
    )
    parser.add_argument(
        "--target-primary-scatter",
        default=None,
        help="Only merge rows with this target_primary_scatter value.",
    )
    parser.add_argument(
        "--allow-mixed-config",
        action="store_true",
        help="Allow merging rows from multiple StageD run configurations.",
    )
    args = parser.parse_args()

    base_dir = Path(args.project_root).resolve()
    rows = load_rows(base_dir, args.ratio)
    if not rows:
        raise SystemExit(f"No StageD summary CSV files found for ratio {args.ratio}")

    if args.scatter_metric is not None:
        rows = [row for row in rows if row.get("scatter_metric", "") == args.scatter_metric]
    if args.target_primary_scatter is not None:
        rows = [
            row for row in rows
            if row.get("target_primary_scatter", "") == str(args.target_primary_scatter)
        ]
    if not rows:
        raise SystemExit(
            f"No StageD summary CSV rows match the requested filters for ratio {args.ratio}"
        )

    configs = {}
    for row in rows:
        signature = row_config_signature(row)
        configs.setdefault(signature, 0)
        configs[signature] += 1

    if len(configs) > 1 and not args.allow_mixed_config:
        lines = [
            "Mixed StageD configurations detected. Re-run with filters or a clean output directory."
        ]
        for signature, count in sorted(configs.items(), key=lambda item: (-item[1], item[0])):
            (
                source_mode,
                boundary_mode,
                reentry_mode,
                particle_reentry_mode,
                matrix_reentry_mode,
                scatter_metric,
                target_primary_scatter,
                theta_threshold_deg,
                matrix_n,
                matrix_abs_um,
                bn_n,
                bn_abs_um,
                zns_n,
                zns_abs_um,
                wavelength_nm,
            ) = signature
            lines.append(
                "  "
                f"count={count} "
                f"source={source_mode} boundary={boundary_mode} "
                f"reentry={reentry_mode}/{particle_reentry_mode}/{matrix_reentry_mode} "
                f"scatter_metric={scatter_metric} target_primary_scatter={target_primary_scatter} "
                f"theta={theta_threshold_deg} wavelength={wavelength_nm} "
                f"optical=({matrix_n},{matrix_abs_um}; {bn_n},{bn_abs_um}; {zns_n},{zns_abs_um})"
            )
        raise SystemExit("\n".join(lines))

    required_any = [
        ("mu_a_count_per_um", "mu_a_raw_per_um"),
        ("mu_a_expected_per_um", "mu_a_raw_per_um"),
        ("mu_s_per_um", "mu_s_raw_per_um"),
        ("g1", "g_raw"),
        ("g2",),
        ("mu_s_prime_direct_per_um", "mu_s_prime_raw_per_um"),
        ("mu_s_bulk_raw_per_um",),
        ("g_bulk_raw",),
        ("mu_s_prime_bulk_raw_per_um",),
        ("mu_s_boundary_raw_per_um",),
        ("g_boundary_raw",),
        ("mu_s_prime_boundary_raw_per_um",),
        ("mu_s_particle_raw_per_um", "mu_s_raw_per_um"),
        ("g_particle_raw", "g_raw"),
        ("mu_s_prime_particle_raw_per_um", "mu_s_prime_raw_per_um"),
        ("mu_s_boundary_primary_raw_per_um", "mu_s_boundary_raw_per_um"),
        ("g_boundary_primary_raw", "g_boundary_raw"),
        ("mu_s_prime_boundary_primary_raw_per_um", "mu_s_prime_boundary_raw_per_um"),
        ("absorbed_fraction",),
        ("mean_path_length_um",),
        ("mean_num_reentry",),
        ("n_photons",),
        ("wavelength_nm",),
    ]
    missing_groups = [group for group in required_any if not any(key in rows[0] for key in group)]
    if missing_groups:
        missing_text = ", ".join("/".join(group) for group in missing_groups)
        raise SystemExit(
            "StageD summary CSV is missing required columns. "
            "Re-run StageD with the current code. Missing: "
            f"{missing_text}"
        )

    mu_a_count_key = primary_key(rows, "mu_a_count_per_um", "mu_a_raw_per_um")
    mu_a_expected_key = primary_key(rows, "mu_a_expected_per_um", "mu_a_raw_per_um")
    mu_s_primary_key = primary_key(rows, "mu_s_per_um", "mu_s_raw_per_um")
    g_primary_key = primary_key(rows, "g1", "g_raw")
    g2_primary_key = primary_key(rows, "g2", "g_raw")
    mu_sp_primary_key = primary_key(rows, "mu_s_prime_direct_per_um", "mu_s_prime_raw_per_um")
    mu_s_total_key = primary_key(rows, "mu_s_step_total_raw_per_um", "mu_s_raw_per_um")
    g_total_key = primary_key(rows, "g_step_total_raw", "g_raw")
    mu_sp_total_key = primary_key(rows, "mu_s_prime_step_total_raw_per_um", "mu_s_prime_raw_per_um")

    mu_a_count = numeric_series(rows, mu_a_count_key)
    mu_a_expected = numeric_series(rows, mu_a_expected_key)
    mu_s_primary = numeric_series(rows, mu_s_primary_key)
    g_primary = numeric_series(rows, g_primary_key)
    g2_primary = numeric_series(rows, g2_primary_key)
    mu_sp_primary = numeric_series(rows, mu_sp_primary_key)
    mu_s_total = numeric_series(rows, mu_s_total_key)
    g_total = numeric_series(rows, g_total_key)
    mu_sp_total = numeric_series(rows, mu_sp_total_key)
    mu_s_bulk = numeric_series(rows, "mu_s_bulk_raw_per_um")
    g_bulk = numeric_series(rows, "g_bulk_raw")
    mu_sp_bulk = numeric_series(rows, "mu_s_prime_bulk_raw_per_um")
    mu_s_particle = optional_numeric_series(rows, "mu_s_particle_raw_per_um", "mu_s_raw_per_um")
    g_particle = optional_numeric_series(rows, "g_particle_raw", "g_raw")
    mu_sp_particle = optional_numeric_series(rows, "mu_s_prime_particle_raw_per_um", "mu_s_prime_raw_per_um")
    mu_s_boundary_primary = optional_numeric_series(rows, "mu_s_boundary_primary_raw_per_um", "mu_s_boundary_raw_per_um")
    g_boundary_primary = optional_numeric_series(rows, "g_boundary_primary_raw", "g_boundary_raw")
    mu_sp_boundary_primary = optional_numeric_series(rows, "mu_s_prime_boundary_primary_raw_per_um", "mu_s_prime_boundary_raw_per_um")
    mu_s_boundary = numeric_series(rows, "mu_s_boundary_raw_per_um")
    g_boundary = numeric_series(rows, "g_boundary_raw")
    mu_sp_boundary = numeric_series(rows, "mu_s_prime_boundary_raw_per_um")
    absorbed = numeric_series(rows, "absorbed_fraction")
    path_length = numeric_series(rows, "mean_path_length_um")
    reentry = numeric_series(rows, "mean_num_reentry")
    photons = [int(float(r["n_photons"])) for r in rows]
    wavelengths = numeric_series(rows, "wavelength_nm")

    out_dir = base_dir / "Output" / "optical_params" / args.ratio
    out_dir.mkdir(parents=True, exist_ok=True)

    scalar_path = out_dir / "rve_raw_optical_params_by_ratio.csv"
    phase_path = out_dir / "phase_function_mean_by_ratio.csv"
    mc_csv_path = out_dir / "monte_carlo_recommended_inputs.csv"
    mc_json_path = out_dir / "monte_carlo_recommended_inputs.json"

    write_scalar_summary(
        scalar_path,
        args,
        wavelengths,
        photons,
        mu_a_count,
        mu_a_expected,
        mu_s_primary,
        g_primary,
        g2_primary,
        mu_sp_primary,
        mu_s_total,
        g_total,
        mu_sp_total,
        mu_s_particle,
        g_particle,
        mu_sp_particle,
        mu_s_boundary_primary,
        g_boundary_primary,
        mu_sp_boundary_primary,
        mu_s_boundary,
        g_boundary,
        mu_sp_boundary,
        absorbed,
        path_length,
        reentry,
    )

    grouped_phase = load_phase_function_rows(rows)
    write_phase_function_summary(phase_path, grouped_phase)

    write_mc_input_csv(
        mc_csv_path,
        args.ratio,
        mean(wavelengths),
        mean(mu_a_expected),
        mean(mu_s_primary),
        mean(g_primary),
        mean(g2_primary),
        mean(mu_sp_primary),
        phase_path.name,
    )
    write_mc_input_json(
        mc_json_path,
        args.ratio,
        mean(wavelengths),
        mean(mu_a_expected),
        mean(mu_s_primary),
        mean(g_primary),
        mean(g2_primary),
        mean(mu_sp_primary),
        phase_path.name,
    )

    print(scalar_path)
    print(phase_path)
    print(mc_csv_path)
    print(mc_json_path)


if __name__ == "__main__":
    main()
