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
    all_raw = True
    for row in summary_rows:
        summary_path = Path(row["_path"])
        phase_name = row.get("phase_function_raw_file")
        if not phase_name:
            all_raw = False
            phase_name = row.get("phase_function_file") or "phase_function.csv"
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
    return grouped, all_raw


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


def numeric_value(row, *keys, default=0.0):
    for key in keys:
        value = row.get(key, "")
        if value != "":
            return float(value)
    return default


def numeric_text_equal(lhs, rhs):
    try:
        return math.isclose(float(lhs), float(rhs), rel_tol=0.0, abs_tol=1.0e-12)
    except (TypeError, ValueError):
        return False


def normalize_scatter_metric(value):
    if value in ("angle_threshold", "particle_encounter_angle_threshold"):
        return "particle_encounter_angle_threshold"
    return value


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
                         g_boundary, mu_sp_boundary, absorbed, path_length, reentry,
                         global_raw):
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
                "global_total_medium_path_um",
                "global_num_encounter_raw",
                "global_mu_s_encounter_raw_per_um",
                "global_g_encounter_raw",
                "global_g2_encounter_raw",
                "global_mu_s_prime_direct_raw_per_um",
                "global_mu_s_prime_from_g_raw_per_um",
                "global_mu_s_prime_consistency_relative_error",
                "global_mu_a_expected_per_um",
                "global_post_first_encounter_mu_s_raw_per_um",
                "global_post_first_encounter_g_raw",
                "global_post_first_encounter_mu_s_prime_raw_per_um",
                "global_post_first_mu_a_count_per_um",
                "global_post_first_mu_a_expected_per_um",
                "global_post_first_medium_path_um",
                "global_post_first_absorption_count",
                *[
                    column
                    for phase in ("Matrix", "BN", "ZnS", "Invalid")
                    for column in (
                        f"global_source_{phase}_post_first_medium_path_um",
                        f"global_source_{phase}_post_first_encounter_count",
                        f"global_post_first_mu_a_source_{phase}_per_um",
                        f"global_post_first_mu_s_source_{phase}_per_um",
                        f"global_post_first_g_source_{phase}",
                        f"global_post_first_mu_s_prime_source_{phase}_per_um",
                    )
                ],
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
                global_raw["total_medium_path_um"],
                global_raw["num_encounter_raw"],
                global_raw["mu_s_raw_per_um"],
                global_raw["g_raw"],
                global_raw["g2_raw"],
                global_raw["mu_s_prime_raw_per_um"],
                global_raw["mu_s_prime_from_g_raw_per_um"],
                global_raw["mu_s_prime_consistency_relative_error"],
                global_raw["mu_a_expected_per_um"],
                global_raw["post_first_mu_s_raw_per_um"],
                global_raw["post_first_g_raw"],
                global_raw["post_first_mu_s_prime_raw_per_um"],
                global_raw["post_first_mu_a_count_per_um"],
                global_raw["post_first_mu_a_expected_per_um"],
                global_raw["post_first_medium_path_um"],
                global_raw["post_first_absorption_count"],
                *[
                    value
                    for phase in ("Matrix", "BN", "ZnS", "Invalid")
                    for value in (
                        global_raw["source_post_first"][phase]["path_um"],
                        global_raw["source_post_first"][phase]["encounter_count"],
                        global_raw["source_post_first"][phase]["mu_a_per_um"],
                        global_raw["source_post_first"][phase]["mu_s_per_um"],
                        global_raw["source_post_first"][phase]["g"],
                        global_raw["source_post_first"][phase]["mu_s_prime_per_um"],
                    )
                ],
                (
                    "Placement mean/std columns are descriptive. Global raw columns use ratio-of-sums "
                    "over all placements and are the formal StageD estimates. Use the aggregated raw "
                    "phase function for native transport, or the explicitly exported g=0 equivalent."
                ),
            ]
        )


def write_phase_function_summary(out_path: Path, grouped):
    total_count = sum(sum(entry["counts"]) for entry in grouped.values())
    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "lambda_nm",
                "bin_id",
                "cos_theta_min",
                "cos_theta_max",
                "count_total",
                "probability",
                "probability_density",
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
            bin_count_total = sum(entry["counts"])
            probability_global = bin_count_total / total_count if total_count > 0.0 else 0.0
            bin_width = entry["cos_theta_max"] - entry["cos_theta_min"]
            writer.writerow(
                [
                    entry["lambda_nm"],
                    bin_id,
                    entry["cos_theta_min"],
                    entry["cos_theta_max"],
                    bin_count_total,
                    probability_global,
                    probability_global / bin_width if bin_width > 0.0 else 0.0,
                    mean(entry["counts"]),
                    stddev(entry["counts"]),
                    mean(entry["probabilities"]),
                    stddev(entry["probabilities"]),
                    mean(entry["densities"]),
                    stddev(entry["densities"]),
                ]
            )


def write_mc_input_csv(out_path: Path, ratio: str, wavelength_nm: float,
                       global_raw, phase_function_file: str,
                       phase_function_is_raw: bool):
    phase_model = (
        "tabulated_phase_function"
        if phase_function_is_raw
        else "legacy_primary_phase_function"
    )
    phase_note = (
        "Native StageD encounter model. Prefer the aggregated raw phase function; HG with g1 is only a fallback approximation."
        if phase_function_is_raw
        else "Legacy runs did not store a separate raw histogram. The referenced phase function follows the old primary metric; re-run StageD for a formal raw phase function."
    )
    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "ratio",
                "wavelength_nm",
                "param_model",
                "recommended_model",
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
                "GO_RVE_NATIVE",
                "explicit_rve_native_phase_function",
                global_raw["mu_a_expected_per_um"],
                global_raw["mu_s_raw_per_um"],
                global_raw["g_raw"],
                global_raw["g2_raw"],
                global_raw["mu_s_prime_raw_per_um"],
                phase_function_file,
                phase_model,
                phase_note,
            ]
        )
        writer.writerow(
            [
                ratio,
                wavelength_nm,
                "GO_RVE_G0_EQUIVALENT",
                "isotropic_g0_full_path",
                global_raw["mu_a_expected_per_um"],
                global_raw["mu_s_prime_raw_per_um"],
                0.0,
                0.0,
                global_raw["mu_s_prime_raw_per_um"],
                "",
                "isotropic",
                "Full-path g=0 equivalent: mu_s equals the StageD direct raw full-path mu_s_prime.",
            ]
        )
        writer.writerow(
            [
                ratio,
                wavelength_nm,
                "GO_RVE_G0_POST_FIRST",
                "isotropic_g0_post_first_moving_photons",
                global_raw["post_first_mu_a_expected_per_um"],
                global_raw["post_first_mu_s_prime_raw_per_um"],
                0.0,
                0.0,
                global_raw["post_first_mu_s_prime_raw_per_um"],
                "",
                "isotropic",
                "Conditional moving-photon model after one complete particle encounter; combine with a separate initial extraction/trapping model when sources originate in ZnS.",
            ]
        )


def write_mc_input_json(out_path: Path, ratio: str, wavelength_nm: float,
                        global_raw, phase_function_file: str,
                        phase_function_is_raw: bool):
    payload = {
        "ratio": ratio,
        "wavelength_nm": wavelength_nm,
        "recommended_models": {
            "explicit_rve_native_phase_function": {
                "param_model": "GO_RVE_NATIVE",
                "mu_a_per_um": global_raw["mu_a_expected_per_um"],
                "mu_s_per_um": global_raw["mu_s_raw_per_um"],
                "g1": global_raw["g_raw"],
                "g2": global_raw["g2_raw"],
                "mu_s_prime_per_um": global_raw["mu_s_prime_raw_per_um"],
                "phase_function_file": phase_function_file,
                "scattering_model": (
                    "tabulated_phase_function"
                    if phase_function_is_raw
                    else "legacy_primary_phase_function"
                ),
                "phase_function_is_raw": phase_function_is_raw,
            },
            "isotropic_g0_full_path": {
                "param_model": "GO_RVE_G0_FULL_PATH",
                "mu_a_per_um": global_raw["mu_a_expected_per_um"],
                "mu_s_per_um": global_raw["mu_s_prime_raw_per_um"],
                "g1": 0.0,
                "g2": 0.0,
                "mu_s_prime_per_um": global_raw["mu_s_prime_raw_per_um"],
                "phase_function_file": "",
                "scattering_model": "isotropic",
            },
            "isotropic_g0_post_first_moving_photons": {
                "param_model": "GO_RVE_G0_POST_FIRST",
                "mu_a_per_um": global_raw["post_first_mu_a_expected_per_um"],
                "mu_s_per_um": global_raw["post_first_mu_s_prime_raw_per_um"],
                "g1": 0.0,
                "g2": 0.0,
                "mu_s_prime_per_um": global_raw["post_first_mu_s_prime_raw_per_um"],
                "phase_function_file": "",
                "scattering_model": "isotropic",
                "conditional_population": "photons_after_one_complete_particle_encounter",
            },
        },
        "fallback_inputs": {
            "scattering_model": "henyey_greenstein",
            "g1": global_raw["g_raw"],
        },
        "notes": [
            "Recommended values use global ratio-of-sums over placements.",
            "StageD native transport is not forced to g=0.",
            "The full-path g=0 equivalent uses mu_s = direct raw full-path mu_s_prime.",
            "The post-first g=0 model is conditional on entering the moving-photon population and does not include initial ZnS extraction/trapping loss.",
            "mu_a_per_um is the expected absorption estimate from phase path lengths and ABSLENGTH.",
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
    parser.add_argument("--source-mode", default=None)
    parser.add_argument("--boundary-mode", default=None)
    parser.add_argument("--reentry-mode", default=None)
    parser.add_argument("--particle-reentry-mode", default=None)
    parser.add_argument("--matrix-reentry-mode", default=None)
    parser.add_argument("--theta-threshold-deg", default=None)
    parser.add_argument("--wavelength-nm", default=None)
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
        requested_metric = normalize_scatter_metric(args.scatter_metric)
        rows = [
            row for row in rows
            if normalize_scatter_metric(row.get("scatter_metric", "")) == requested_metric
        ]
    if args.target_primary_scatter is not None:
        rows = [
            row for row in rows
            if numeric_text_equal(row.get("target_primary_scatter", ""),
                                  args.target_primary_scatter)
        ]
    exact_filters = {
        "source_mode": args.source_mode,
        "boundary_mode": args.boundary_mode,
        "reentry_mode": args.reentry_mode,
        "particle_reentry_mode": args.particle_reentry_mode,
        "matrix_reentry_mode": args.matrix_reentry_mode,
    }
    for key, value in exact_filters.items():
        if value is not None:
            rows = [row for row in rows if row.get(key, "") == str(value)]
    numeric_filters = {
        "theta_threshold_deg": args.theta_threshold_deg,
        "wavelength_nm": args.wavelength_nm,
    }
    for key, value in numeric_filters.items():
        if value is not None:
            rows = [
                row for row in rows
                if numeric_text_equal(row.get(key, ""), value)
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

    total_medium_path = sum(
        numeric_value(row, "total_medium_path_length_um") for row in rows
    )
    total_encounter_raw = sum(
        numeric_value(row, "num_encounter_total", "total_complete_encounter")
        for row in rows
    )
    total_sum_cos_raw = sum(
        numeric_value(
            row,
            "stageD_sum_cos_theta_encounter_raw",
            default=(
                numeric_value(row, "stageD_g_encounter_raw", "g_raw")
                * numeric_value(row, "num_encounter_total", "total_complete_encounter")
            ),
        )
        for row in rows
    )
    total_sum_one_minus_cos_raw = sum(
        numeric_value(
            row,
            "stageD_sum_one_minus_cos_theta_encounter_raw",
            default=(
                numeric_value(
                    row,
                    "stageD_mu_s_prime_direct_raw_per_um",
                    "mu_s_prime_raw_per_um",
                )
                * numeric_value(row, "total_medium_path_length_um")
            ),
        )
        for row in rows
    )
    total_sum_cos2_raw = sum(
        ((2.0 * numeric_value(row, "stageD_g2_encounter_raw", "g2") + 1.0) / 3.0)
        * numeric_value(row, "num_encounter_total", "total_complete_encounter")
        for row in rows
    )
    global_mu_s_raw = (
        total_encounter_raw / total_medium_path if total_medium_path > 0.0 else 0.0
    )
    global_g_raw = (
        total_sum_cos_raw / total_encounter_raw if total_encounter_raw > 0.0 else 0.0
    )
    global_g2_raw = (
        0.5 * (3.0 * total_sum_cos2_raw / total_encounter_raw - 1.0)
        if total_encounter_raw > 0.0
        else 0.0
    )
    global_mu_sp_raw = (
        total_sum_one_minus_cos_raw / total_medium_path
        if total_medium_path > 0.0
        else 0.0
    )
    global_mu_sp_from_g = global_mu_s_raw * (1.0 - global_g_raw)
    global_consistency_scale = max(
        1.0e-15, abs(global_mu_sp_raw), abs(global_mu_sp_from_g)
    )

    expected_absorption_numerator = sum(
        numeric_value(row, "mu_a_expected_per_um", "mu_a_raw_per_um")
        * numeric_value(row, "total_medium_path_length_um")
        for row in rows
    )
    global_mu_a_expected = (
        expected_absorption_numerator / total_medium_path
        if total_medium_path > 0.0
        else mean(mu_a_expected)
    )

    post_first_path = sum(
        numeric_value(row, "post_first_encounter_medium_path_um") for row in rows
    )
    post_first_count = sum(
        numeric_value(row, "post_first_encounter_count_raw") for row in rows
    )
    post_first_sum_cos = sum(
        numeric_value(row, "post_first_encounter_g_raw")
        * numeric_value(row, "post_first_encounter_count_raw")
        for row in rows
    )
    post_first_sum_one_minus_cos = sum(
        numeric_value(row, "post_first_encounter_mu_s_prime_raw_per_um")
        * numeric_value(row, "post_first_encounter_medium_path_um")
        for row in rows
    )
    post_first_absorption_count = sum(
        numeric_value(row, "post_first_absorption_count") for row in rows
    )
    post_first_expected_absorption_numerator = sum(
        numeric_value(row, "post_first_mu_a_expected_per_um")
        * numeric_value(row, "post_first_encounter_medium_path_um")
        for row in rows
    )
    source_post_first = {}
    for phase in ("Matrix", "BN", "ZnS", "Invalid"):
        path = sum(
            numeric_value(row, f"source_{phase}_post_first_medium_path_um")
            for row in rows
        )
        encounter_count = sum(
            numeric_value(row, f"source_{phase}_post_first_encounter_count")
            for row in rows
        )
        absorption_count = sum(
            numeric_value(row, f"source_{phase}_post_first_absorption_count")
            for row in rows
        )
        sum_cos = sum(
            numeric_value(row, f"source_{phase}_post_first_sum_cos_theta")
            for row in rows
        )
        sum_one_minus_cos = sum(
            numeric_value(row, f"source_{phase}_post_first_sum_one_minus_cos_theta")
            for row in rows
        )
        source_post_first[phase] = {
            "path_um": path,
            "encounter_count": encounter_count,
            "mu_a_per_um": absorption_count / path if path > 0.0 else 0.0,
            "mu_s_per_um": encounter_count / path if path > 0.0 else 0.0,
            "g": sum_cos / encounter_count if encounter_count > 0.0 else 0.0,
            "mu_s_prime_per_um": (
                sum_one_minus_cos / path if path > 0.0 else 0.0
            ),
        }
    global_raw = {
        "total_medium_path_um": total_medium_path,
        "num_encounter_raw": total_encounter_raw,
        "mu_s_raw_per_um": global_mu_s_raw,
        "g_raw": global_g_raw,
        "g2_raw": global_g2_raw,
        "mu_s_prime_raw_per_um": global_mu_sp_raw,
        "mu_s_prime_from_g_raw_per_um": global_mu_sp_from_g,
        "mu_s_prime_consistency_relative_error": (
            abs(global_mu_sp_raw - global_mu_sp_from_g) / global_consistency_scale
        ),
        "mu_a_expected_per_um": global_mu_a_expected,
        "post_first_mu_s_raw_per_um": (
            post_first_count / post_first_path if post_first_path > 0.0 else 0.0
        ),
        "post_first_g_raw": (
            post_first_sum_cos / post_first_count if post_first_count > 0.0 else 0.0
        ),
        "post_first_mu_s_prime_raw_per_um": (
            post_first_sum_one_minus_cos / post_first_path
            if post_first_path > 0.0
            else 0.0
        ),
        "post_first_mu_a_count_per_um": (
            post_first_absorption_count / post_first_path
            if post_first_path > 0.0
            else 0.0
        ),
        "post_first_mu_a_expected_per_um": (
            post_first_expected_absorption_numerator / post_first_path
            if post_first_path > 0.0
            else 0.0
        ),
        "post_first_medium_path_um": post_first_path,
        "post_first_absorption_count": post_first_absorption_count,
        "source_post_first": source_post_first,
    }

    out_dir = base_dir / "Output" / "optical_params" / args.ratio
    out_dir.mkdir(parents=True, exist_ok=True)

    scalar_path = out_dir / "rve_raw_optical_params_by_ratio.csv"
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
        global_raw,
    )

    grouped_phase, all_phase_functions_raw = load_phase_function_rows(rows)
    phase_path = out_dir / (
        "phase_function_raw_by_ratio.csv"
        if all_phase_functions_raw
        else "phase_function_legacy_primary_by_ratio.csv"
    )
    legacy_phase_path = out_dir / "phase_function_mean_by_ratio.csv"
    write_phase_function_summary(phase_path, grouped_phase)
    legacy_phase_path.write_bytes(phase_path.read_bytes())

    write_mc_input_csv(
        mc_csv_path,
        args.ratio,
        mean(wavelengths),
        global_raw,
        phase_path.name,
        all_phase_functions_raw,
    )
    write_mc_input_json(
        mc_json_path,
        args.ratio,
        mean(wavelengths),
        global_raw,
        phase_path.name,
        all_phase_functions_raw,
    )

    print(scalar_path)
    print(phase_path)
    print(mc_csv_path)
    print(mc_json_path)


if __name__ == "__main__":
    main()
