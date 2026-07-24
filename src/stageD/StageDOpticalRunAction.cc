#include "StageDOpticalRunAction.hh"

#include "AnalysisConfig.hh"
#include "DetectorConstruction.hh"

#include "G4Exception.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <sstream>

namespace
{
  std::string WeightPartToTagString(G4double value)
  {
    const G4double rounded = std::round(value);
    std::ostringstream oss;
    if (std::abs(value - rounded) < 1.0e-9)
      oss << static_cast<long long>(rounded);
    else
      oss << value;
    return oss.str();
  }

  std::string CsvQuote(const std::string &value)
  {
    if (value.find_first_of(",\"") == std::string::npos)
      return value;
    std::string out = "\"";
    for (char ch : value)
    {
      if (ch == '"')
        out += "\"\"";
      else
        out += ch;
    }
    out += "\"";
    return out;
  }

  G4bool UseThresholdedEncounterMetric(const AnalysisConfig *config)
  {
    return config != nullptr &&
           config->stageD_scatter_metric == "particle_encounter_angle_threshold";
  }

  std::string PrimaryScatterMetricLabel(const AnalysisConfig *config)
  {
    if (!UseThresholdedEncounterMetric(config))
      return "particle_encounter_no_threshold";
    return config->stageD_scatter_metric.empty()
               ? "particle_encounter_angle_threshold"
               : config->stageD_scatter_metric;
  }

  G4bool KeepReentryDiagnostic(const StageDReentryDiagnosticRecord &record,
                               G4double samplingRate)
  {
    if (samplingRate <= 0.0)
      return false;
    if (samplingRate >= 1.0)
      return true;

    std::uint64_t value =
        (static_cast<std::uint64_t>(static_cast<std::uint32_t>(record.event_id)) << 32) ^
        static_cast<std::uint32_t>(record.reentry_index);
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    value ^= value >> 31;
    const G4double unitValue =
        static_cast<G4double>(value >> 11) * (1.0 / 9007199254740992.0);
    return unitValue < samplingRate;
  }

  std::string SafePathComponent(const std::string &value)
  {
    std::string result = value.empty() ? "unset" : value;
    std::transform(result.begin(), result.end(), result.begin(), [](unsigned char ch)
                   { return (std::isalnum(ch) || ch == '-' || ch == '_')
                                ? static_cast<char>(ch)
                                : '_'; });
    return result;
  }

  std::uint64_t StableHash(const std::string &value)
  {
    std::uint64_t hash = 1469598103934665603ULL;
    for (const unsigned char ch : value)
    {
      hash ^= ch;
      hash *= 1099511628211ULL;
    }
    return hash;
  }
}

StageDOpticalRunAction::StageDOpticalRunAction(AnalysisConfig *config)
    : G4UserRunAction(),
      fConfig(config),
      fEventsCsv(),
      fReentryDiagnosticsCsv(),
      fEventsCsvPath(""),
      fSummaryCsvPath(""),
      fPhaseFunctionCsvPath(""),
      fReentryDiagnosticsCsvPath(""),
      fOutputDir(""),
      fRatioTag(""),
      fPlacementFile(""),
      fPlacementStem(""),
      fConfigurationTag(""),
      fRunId(0),
      fRandomSeed0(0),
      fRandomSeed1(0),
      fAccumulator(),
      fReentryPortalSummary(),
      fReentryDiagnosticRowsWritten(0)
{
}

StageDOpticalRunAction::~StageDOpticalRunAction()
{
  if (fEventsCsv.is_open())
    fEventsCsv.close();
  if (fReentryDiagnosticsCsv.is_open())
    fReentryDiagnosticsCsv.close();
}

std::string StageDOpticalRunAction::MakeRatioTag() const
{
  if (fConfig == nullptr)
    return "unknown";
  return WeightPartToTagString(fConfig->bnWt) + "-" +
         WeightPartToTagString(fConfig->znsWt);
}

std::string StageDOpticalRunAction::MakePlacementStem() const
{
  if (fPlacementFile.empty())
    return "unknown_placement";

  const std::filesystem::path placementPath(fPlacementFile);
  const std::filesystem::path normalized = placementPath.lexically_normal();
  std::vector<std::string> suffixParts;
  bool underPlacementsRoot = false;

  std::vector<std::string> parts;
  for (const auto &part : normalized)
    parts.push_back(part.string());

  for (std::size_t i = 0; i + 1 < parts.size(); ++i)
  {
    if (parts[i] == "placements")
    {
      underPlacementsRoot = true;
      for (std::size_t j = i + 2; j < parts.size(); ++j)
        suffixParts.push_back(parts[j]);
      break;
    }
  }

  if (!underPlacementsRoot || suffixParts.empty())
    suffixParts.push_back(normalized.filename().string());

  if (!suffixParts.empty())
  {
    const std::filesystem::path tailPath = suffixParts.back();
    suffixParts.back() = tailPath.stem().string();
  }

  std::ostringstream oss;
  for (std::size_t i = 0; i < suffixParts.size(); ++i)
  {
    if (i > 0)
      oss << "__";
    oss << suffixParts[i];
  }
  return oss.str();
}

std::string StageDOpticalRunAction::MakeWavelengthTag() const
{
  const G4double wavelengthNm = fConfig ? fConfig->stageD_wavelength_nm : 0.0;
  const G4double rounded = std::round(wavelengthNm);
  std::ostringstream oss;
  oss << "lambda_";
  if (std::abs(wavelengthNm - rounded) < 1.0e-9)
    oss << static_cast<long long>(rounded);
  else
    oss << std::fixed << std::setprecision(3) << wavelengthNm;
  oss << "nm";
  return oss.str();
}

std::string StageDOpticalRunAction::MakeConfigurationTag() const
{
  if (fConfig == nullptr)
    return "cfg_unknown";

  std::ostringstream signature;
  signature << std::setprecision(15)
            << fConfig->stageD_source_mode << '|'
            << fConfig->stageD_boundary_mode << '|'
            << fConfig->stageD_reentry_mode << '|'
            << fConfig->stageD_particle_reentry_mode << '|'
            << fConfig->stageD_matrix_reentry_mode << '|'
            << fConfig->stageD_scatter_metric << '|'
            << fConfig->stageD_target_primary_scatter << '|'
            << fConfig->stageD_theta_threshold_deg << '|'
            << fConfig->stageD_max_reentry << '|'
            << fConfig->stageD_max_steps << '|'
            << fConfig->stageD_max_path_length_um << '|'
            << fConfig->stageD_portal_nu << '|'
            << fConfig->stageD_portal_nv << '|'
            << fConfig->stageD_portal_margin_um << '|'
            << fConfig->stageD_clearance_bin0_um << '|'
            << fConfig->stageD_clearance_bin1_um << '|'
            << fConfig->stageD_clearance_bin2_um << '|'
            << fConfig->stageD_max_particle_reentry_trials << '|'
            << fConfig->stageD_max_portal_fallback_level << '|'
            << fConfig->opticalMatrixRIndex << '|'
            << fConfig->opticalMatrixAbsLengthUm << '|'
            << fConfig->opticalBnRIndex << '|'
            << fConfig->opticalBnAbsLengthUm << '|'
            << fConfig->opticalZnsRIndex << '|'
            << fConfig->opticalZnsAbsLengthUm;

  std::ostringstream tag;
  tag << "source_" << SafePathComponent(fConfig->stageD_source_mode)
      << "__boundary_" << SafePathComponent(fConfig->stageD_boundary_mode)
      << "__cfg_" << std::hex << std::setw(12) << std::setfill('0')
      << (StableHash(signature.str()) & 0xffffffffffffULL);
  return tag.str();
}

std::string StageDOpticalRunAction::ResolveOutputDirectory() const
{
  if (fConfig != nullptr && !fConfig->stageD_output_dir.empty())
    return fConfig->stageD_output_dir;

  return (AnalysisConfig::ProjectRootPath() /
          "Output" /
          "stageD_optical_homogenization" /
          fRatioTag /
          fPlacementStem /
          MakeWavelengthTag() /
          fConfigurationTag /
          ("run_" + std::to_string(fRunId) +
           "_seeds_" + std::to_string(fRandomSeed0) +
           "_" + std::to_string(fRandomSeed1)))
      .string();
}

void StageDOpticalRunAction::WriteEventHeader()
{
  fEventsCsv
      << "photonID,"
      << "ratio,"
      << "placement_file,"
      << "source_mode,"
      << "boundary_mode,"
      << "reentry_mode,"
      << "particle_reentry_mode,"
      << "matrix_reentry_mode,"
      << "wavelength_nm,"
      << "source_phase,"
      << "source_x_um,"
      << "source_y_um,"
      << "source_z_um,"
      << "final_status,"
      << "absorbed,"
      << "total_path_length_um,"
      << "path_length_bn_um,"
      << "path_length_zns_um,"
      << "path_length_matrix_um,"
      << "path_length_world_um,"
      << "num_steps,"
      << "num_absorbed_total,"
      << "num_absorbed_BN,"
      << "num_absorbed_ZnS,"
      << "num_absorbed_Matrix,"
      << "num_absorbed_World,"
      << "num_encounter_total,"
      << "num_encounter_BN,"
      << "num_encounter_ZnS,"
      << "num_encounter_effective_total,"
      << "num_encounter_effective_BN,"
      << "num_encounter_effective_ZnS,"
      << "sum_cos_theta_encounter,"
      << "sum_cos_theta_encounter_BN,"
      << "sum_cos_theta_encounter_ZnS,"
      << "sum_cos_theta_encounter_effective,"
      << "sum_cos_theta_encounter_effective_BN,"
      << "sum_cos_theta_encounter_effective_ZnS,"
      << "sum_one_minus_cos_theta_encounter,"
      << "sum_one_minus_cos_theta_encounter_BN,"
      << "sum_one_minus_cos_theta_encounter_ZnS,"
      << "sum_one_minus_cos_theta_encounter_effective,"
      << "sum_one_minus_cos_theta_encounter_effective_BN,"
      << "sum_one_minus_cos_theta_encounter_effective_ZnS,"
      << "sum_cos2_theta_encounter,"
      << "sum_cos2_theta_encounter_BN,"
      << "sum_cos2_theta_encounter_ZnS,"
      << "sum_cos2_theta_encounter_effective,"
      << "sum_cos2_theta_encounter_effective_BN,"
      << "sum_cos2_theta_encounter_effective_ZnS,"
      << "g1_encounter_for_this_photon,"
      << "g2_encounter_for_this_photon,"
      << "mu_s_prime_direct_encounter_per_um_for_this_photon,"
      << "g1_encounter_raw_for_this_photon,"
      << "g2_encounter_raw_for_this_photon,"
      << "mu_s_prime_direct_encounter_raw_per_um_for_this_photon,"
      << "num_particle_scatter_legacy,"
      << "num_particle_scatter_BN_legacy,"
      << "num_particle_scatter_ZnS_legacy,"
      << "num_real_scatter_debug,"
      << "num_bulk_scatter_debug,"
      << "num_boundary_scatter_debug,"
      << "num_boundary_scatter_BN_debug,"
      << "num_boundary_scatter_ZnS_debug,"
      << "num_material_boundary,"
      << "num_complete_encounter_total,"
      << "num_complete_encounter_BN,"
      << "num_complete_encounter_ZnS,"
      << "num_surface_reflection_encounter,"
      << "num_incomplete_initial_particle_exit,"
      << "num_censored_particle_encounter,"
      << "num_inconsistent_encounter_state,"
      << "num_particle_to_particle_boundary,"
      << "num_unknown_particle_reflection,"
      << "num_outer_boundary_hits,"
      << "num_outer_boundary_reentry_success,"
      << "num_outer_boundary_reentry_failed,"
      << "num_outer_boundary_escape,"
      << "num_outer_boundary_max_reentry_stop,"
      << "num_outer_boundary_reentry_algorithm_failed,"
      << "num_outer_boundary_periodic_phase_mismatch,"
      << "num_outer_boundary_fresnel_reflection,"
      << "num_outer_boundary_total_internal_reflection,"
      << "num_outer_boundary_refraction,"
      << "num_outer_boundary_transmission,"
      << "num_outer_boundary_other_status,"
      << "num_reentry,"
      << "num_reentry_BN,"
      << "num_reentry_ZnS,"
      << "num_reentry_matrix,"
      << "num_reentry_particle_q_mu,"
      << "num_reentry_particle_q_only,"
      << "num_reentry_particle_q_only_fallback,"
      << "num_reentry_particle_volume_random,"
      << "num_reentry_matrix_clearance_portal,"
      << "num_reentry_fallback_same_bin,"
      << "num_reentry_fallback_adjacent_bin,"
      << "num_reentry_fallback_any_bin,"
      << "num_reentry_fallback_any_phase_same_bin,"
      << "num_reentry_fallback_any_portal,"
      << "num_reentry_random_matrix_debug,"
      << "num_reentry_failed,"
      << "sum_cos_theta_particle_legacy,"
      << "sum_cos_theta_particle_BN_legacy,"
      << "sum_cos_theta_particle_ZnS_legacy,"
      << "sum_cos_theta_debug,"
      << "sum_cos_theta_bulk_debug,"
      << "sum_cos_theta_boundary_debug,"
      << "sum_cos_theta_boundary_BN_debug,"
      << "sum_cos_theta_boundary_ZnS_debug,"
      << "mean_cos_theta_particle_for_this_photon_legacy,"
      << "mean_cos_theta_for_this_photon_debug,"
      << "mean_cos_theta_bulk_for_this_photon_debug,"
      << "mean_cos_theta_boundary_for_this_photon_debug,"
      << "weight"
      << "\n";
}

void StageDOpticalRunAction::WriteReentryDiagnosticsHeader()
{
  fReentryDiagnosticsCsv
      << "event_id,"
      << "reentry_index,"
      << "strategy,"
      << "fallback_level,"
      << "exit_phase,"
      << "entry_phase,"
      << "old_dir_x,"
      << "old_dir_y,"
      << "old_dir_z,"
      << "exit_x_um,"
      << "exit_y_um,"
      << "exit_z_um,"
      << "entry_x_um,"
      << "entry_y_um,"
      << "entry_z_um,"
      << "particle_q_exit,"
      << "particle_q_entry,"
      << "particle_mu_exit,"
      << "particle_mu_entry,"
      << "matrix_clearance_exit_um,"
      << "matrix_clearance_entry_um,"
      << "matrix_nearest_phase_exit,"
      << "matrix_nearest_phase_entry,"
      << "matrix_clearance_bin_exit,"
      << "matrix_clearance_bin_entry,"
      << "trials"
      << "\n";
}

void StageDOpticalRunAction::OpenOutputs()
{
  std::error_code ec;
  std::filesystem::create_directories(fOutputDir, ec);
  if (ec)
  {
    G4Exception("StageDOpticalRunAction::OpenOutputs",
                "BNZS_D_RUN_001", FatalException,
                ("Failed to create Stage D output directory: " + fOutputDir).c_str());
    return;
  }

  fEventsCsvPath = (std::filesystem::path(fOutputDir) /
                    "stageD_events.csv")
                       .string();
  fReentryDiagnosticsCsvPath = (std::filesystem::path(fOutputDir) /
                                "stageD_reentry_diagnostics.csv")
                                   .string();
  fSummaryCsvPath = (std::filesystem::path(fOutputDir) /
                     "stageD_summary.csv")
                        .string();
  fPhaseFunctionCsvPath = (std::filesystem::path(fOutputDir) /
                           "phase_function.csv")
                              .string();

  fEventsCsv.open(fEventsCsvPath.c_str(), std::ios::out);
  if (!fEventsCsv)
  {
    G4Exception("StageDOpticalRunAction::OpenOutputs",
                "BNZS_D_RUN_002", FatalException,
                ("Failed to open Stage D events CSV: " + fEventsCsvPath).c_str());
    return;
  }
  fEventsCsv << std::setprecision(15);

  if (fConfig != nullptr && fConfig->stageD_write_reentry_diagnostics &&
      fConfig->stageD_max_diagnostic_rows > 0 &&
      fConfig->stageD_reentry_diagnostics_sampling_rate > 0.0)
  {
    fReentryDiagnosticsCsv.open(fReentryDiagnosticsCsvPath.c_str(), std::ios::out);
    if (!fReentryDiagnosticsCsv)
    {
      G4Exception("StageDOpticalRunAction::OpenOutputs",
                  "BNZS_D_RUN_005", FatalException,
                  ("Failed to open Stage D re-entry diagnostics CSV: " + fReentryDiagnosticsCsvPath).c_str());
      return;
    }
    fReentryDiagnosticsCsv << std::setprecision(15);
    WriteReentryDiagnosticsHeader();
  }

  WriteEventHeader();
}

void StageDOpticalRunAction::WriteGeometryMetadataFile(
    const DetectorConstruction *detector) const
{
  if (detector == nullptr)
    return;
  const std::filesystem::path path =
      std::filesystem::path(fOutputDir) / "rve_geometry_metadata.csv";
  std::ofstream output(path);
  if (!output)
  {
    G4Exception("StageDOpticalRunAction::WriteGeometryMetadataFile",
                "BNZS_D_RUN_006", JustWarning,
                ("Failed to write geometry metadata: " + path.string()).c_str());
    return;
  }
  output << "format_version,seed,placement_file,periodic_images_file,geometry_mode,"
            "box_x_um,box_y_um,box_z_um,unique_particles,geometry_copies,radius_classes,"
            "phi_achieved,zns_to_bn_mass_ratio\n"
         << std::setprecision(15)
         << detector->GetPlacementFormatVersion() << ','
         << detector->GetPlacementSeed() << ','
         << CsvQuote(detector->GetLoadedPlacementFileForRecord()) << ','
         << CsvQuote(detector->GetLoadedPeriodicImagesFileForRecord()) << ','
         << CsvQuote(fConfig != nullptr ? fConfig->placementGeometryMode : "unknown") << ','
         << detector->GetBoxXUm() << ','
         << detector->GetBoxYUm() << ','
         << detector->GetBoxZUm() << ','
         << detector->GetUniqueParticleCount() << ','
         << detector->GetGeometryCopyCount() << ','
         << detector->GetRadiusClassCount() << ','
         << detector->GetPlacementPhiAchieved() << ','
         << detector->GetPlacementZnSToBNMassRatio() << '\n';
}

void StageDOpticalRunAction::BeginOfRunAction(const G4Run *run)
{
  fAccumulator = StageDRunAccumulator{};
  fReentryPortalSummary = StageDReentryPortalSummary{};
  fReentryDiagnosticRowsWritten = 0;
  if (fEventsCsv.is_open())
    fEventsCsv.close();
  if (fReentryDiagnosticsCsv.is_open())
    fReentryDiagnosticsCsv.close();

  const auto *detector = dynamic_cast<const DetectorConstruction *>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  fRatioTag = MakeRatioTag();
  fPlacementFile = detector ? detector->GetLoadedPlacementFileForRecord() : "unknown";
  fPlacementStem = MakePlacementStem();
  fRunId = run != nullptr ? run->GetRunID() : 0;
  const G4long *randomSeeds = G4Random::getTheSeeds();
  fRandomSeed0 = randomSeeds != nullptr ? randomSeeds[0] : G4Random::getTheSeed();
  fRandomSeed1 = randomSeeds != nullptr ? randomSeeds[1] : 0;
  fConfigurationTag = MakeConfigurationTag();
  fOutputDir = ResolveOutputDirectory();

  OpenOutputs();
  WriteGeometryMetadataFile(detector);
}

void StageDOpticalRunAction::RecordPhotonEvent(const StageDPhotonEventRecord &event)
{
  AccumulatePhotonEvent(event);

  if (!fEventsCsv.is_open())
    return;

  fEventsCsv
      << event.photonID << ","
      << CsvQuote(event.ratio) << ","
      << CsvQuote(event.placement_file) << ","
      << CsvQuote(event.source_mode) << ","
      << CsvQuote(event.boundary_mode) << ","
      << CsvQuote(event.reentry_mode) << ","
      << CsvQuote(event.particle_reentry_mode) << ","
      << CsvQuote(event.matrix_reentry_mode) << ","
      << event.wavelength_nm << ","
      << CsvQuote(event.source_phase) << ","
      << event.source_x_um << ","
      << event.source_y_um << ","
      << event.source_z_um << ","
      << CsvQuote(event.final_status) << ","
      << (event.absorbed ? 1 : 0) << ","
      << event.total_path_length_um << ","
      << event.path_length_bn_um << ","
      << event.path_length_zns_um << ","
      << event.path_length_matrix_um << ","
      << event.path_length_world_um << ","
      << event.num_steps << ","
      << event.num_absorbed_total << ","
      << event.num_absorbed_BN << ","
      << event.num_absorbed_ZnS << ","
      << event.num_absorbed_Matrix << ","
      << event.num_absorbed_World << ","
      << event.num_encounter_total << ","
      << event.num_encounter_BN << ","
      << event.num_encounter_ZnS << ","
      << event.num_encounter_effective_total << ","
      << event.num_encounter_effective_BN << ","
      << event.num_encounter_effective_ZnS << ","
      << event.sum_cos_theta_encounter << ","
      << event.sum_cos_theta_encounter_BN << ","
      << event.sum_cos_theta_encounter_ZnS << ","
      << event.sum_cos_theta_encounter_effective << ","
      << event.sum_cos_theta_encounter_effective_BN << ","
      << event.sum_cos_theta_encounter_effective_ZnS << ","
      << event.sum_one_minus_cos_theta_encounter << ","
      << event.sum_one_minus_cos_theta_encounter_BN << ","
      << event.sum_one_minus_cos_theta_encounter_ZnS << ","
      << event.sum_one_minus_cos_theta_encounter_effective << ","
      << event.sum_one_minus_cos_theta_encounter_effective_BN << ","
      << event.sum_one_minus_cos_theta_encounter_effective_ZnS << ","
      << event.sum_cos2_theta_encounter << ","
      << event.sum_cos2_theta_encounter_BN << ","
      << event.sum_cos2_theta_encounter_ZnS << ","
      << event.sum_cos2_theta_encounter_effective << ","
      << event.sum_cos2_theta_encounter_effective_BN << ","
      << event.sum_cos2_theta_encounter_effective_ZnS << ","
      << event.g1_encounter_for_this_photon << ","
      << event.g2_encounter_for_this_photon << ","
      << event.mu_s_prime_direct_encounter_per_um_for_this_photon << ","
      << event.g1_encounter_raw_for_this_photon << ","
      << event.g2_encounter_raw_for_this_photon << ","
      << event.mu_s_prime_direct_encounter_raw_per_um_for_this_photon << ","
      << event.num_particle_scatter << ","
      << event.num_particle_scatter_BN << ","
      << event.num_particle_scatter_ZnS << ","
      << event.num_real_scatter << ","
      << event.num_bulk_scatter << ","
      << event.num_boundary_scatter << ","
      << event.num_boundary_scatter_BN << ","
      << event.num_boundary_scatter_ZnS << ","
      << event.num_material_boundary << ","
      << event.num_complete_encounter_total << ","
      << event.num_complete_encounter_BN << ","
      << event.num_complete_encounter_ZnS << ","
      << event.num_surface_reflection_encounter << ","
      << event.num_incomplete_initial_particle_exit << ","
      << event.num_censored_particle_encounter << ","
      << event.num_inconsistent_encounter_state << ","
      << event.num_particle_to_particle_boundary << ","
      << event.num_unknown_particle_reflection << ","
      << event.num_outer_boundary_hits << ","
      << event.num_outer_boundary_reentry_success << ","
      << event.num_outer_boundary_reentry_failed << ","
      << event.num_outer_boundary_escape << ","
      << event.num_outer_boundary_max_reentry_stop << ","
      << event.num_outer_boundary_reentry_algorithm_failed << ","
      << event.num_outer_boundary_periodic_phase_mismatch << ","
      << event.num_outer_boundary_fresnel_reflection << ","
      << event.num_outer_boundary_total_internal_reflection << ","
      << event.num_outer_boundary_refraction << ","
      << event.num_outer_boundary_transmission << ","
      << event.num_outer_boundary_other_status << ","
      << event.num_reentry << ","
      << event.num_reentry_BN << ","
      << event.num_reentry_ZnS << ","
      << event.num_reentry_matrix << ","
      << event.num_reentry_particle_q_mu << ","
      << event.num_reentry_particle_q_only << ","
      << event.num_reentry_particle_q_only_fallback << ","
      << event.num_reentry_particle_volume_random << ","
      << event.num_reentry_matrix_clearance_portal << ","
      << event.num_reentry_fallback_same_bin << ","
      << event.num_reentry_fallback_adjacent_bin << ","
      << event.num_reentry_fallback_any_bin << ","
      << event.num_reentry_fallback_any_phase_same_bin << ","
      << event.num_reentry_fallback_any_portal << ","
      << event.num_reentry_random_matrix_debug << ","
      << event.num_reentry_failed << ","
      << event.sum_cos_theta_particle << ","
      << event.sum_cos_theta_particle_BN << ","
      << event.sum_cos_theta_particle_ZnS << ","
      << event.sum_cos_theta << ","
      << event.sum_cos_theta_bulk << ","
      << event.sum_cos_theta_boundary << ","
      << event.sum_cos_theta_boundary_BN << ","
      << event.sum_cos_theta_boundary_ZnS << ","
      << event.mean_cos_theta_particle_for_this_photon << ","
      << event.mean_cos_theta_for_this_photon << ","
      << event.mean_cos_theta_bulk_for_this_photon << ","
      << event.mean_cos_theta_boundary_for_this_photon << ","
      << event.weight
      << "\n";
}

void StageDOpticalRunAction::AccumulatePhotonEvent(
    const StageDPhotonEventRecord &event)
{
  auto &acc = fAccumulator;
  ++acc.nPhotons;
  if (event.absorbed)
    ++acc.nAbsorbed;
  acc.nAbsorbedBN += event.num_absorbed_BN;
  acc.nAbsorbedZnS += event.num_absorbed_ZnS;
  acc.nAbsorbedMatrix += event.num_absorbed_Matrix;

  if (event.final_status == "lost")
    ++acc.nLost;
  else if (event.final_status == "reentry_failed")
    ++acc.nReentryFailed;
  else if (event.final_status == "escaped_debug")
    ++acc.nEscapedDebug;
  else if (event.final_status == "periodic_phase_mismatch")
  {
    ++acc.nPeriodicPhaseMismatch;
    ++acc.nReentryFailed;
  }
  else if (event.final_status == "periodic_geometry_required")
  {
    ++acc.nPeriodicGeometryRequired;
    ++acc.nReentryFailed;
  }
  else if (event.final_status == "detected")
    ++acc.nDetected;
  else if (event.final_status == "max_reentry")
    ++acc.nMaxReentry;
  else if (event.final_status == "max_steps")
    ++acc.nMaxSteps;
  else if (event.final_status == "max_path_length")
    ++acc.nMaxPathLength;
  else if (event.final_status == "target_primary_scatter")
    ++acc.nTargetPrimaryScatter;

  acc.totalPathLengthUm += event.total_path_length_um;
  acc.totalPathLengthBNUm += event.path_length_bn_um;
  acc.totalPathLengthZnSUm += event.path_length_zns_um;
  acc.totalPathLengthMatrixUm += event.path_length_matrix_um;
  acc.totalPathLengthWorldUm += event.path_length_world_um;
  acc.totalEncounter += event.num_encounter_total;
  acc.totalEncounterBN += event.num_encounter_BN;
  acc.totalEncounterZnS += event.num_encounter_ZnS;
  acc.totalEncounterEffective += event.num_encounter_effective_total;
  acc.totalEncounterEffectiveBN += event.num_encounter_effective_BN;
  acc.totalEncounterEffectiveZnS += event.num_encounter_effective_ZnS;
  acc.sumCosThetaEncounter += event.sum_cos_theta_encounter;
  acc.sumCosThetaEncounterBN += event.sum_cos_theta_encounter_BN;
  acc.sumCosThetaEncounterZnS += event.sum_cos_theta_encounter_ZnS;
  acc.sumCosThetaEncounterEffective += event.sum_cos_theta_encounter_effective;
  acc.sumCosThetaEncounterEffectiveBN += event.sum_cos_theta_encounter_effective_BN;
  acc.sumCosThetaEncounterEffectiveZnS += event.sum_cos_theta_encounter_effective_ZnS;
  acc.sumOneMinusCosThetaEncounter += event.sum_one_minus_cos_theta_encounter;
  acc.sumOneMinusCosThetaEncounterBN += event.sum_one_minus_cos_theta_encounter_BN;
  acc.sumOneMinusCosThetaEncounterZnS += event.sum_one_minus_cos_theta_encounter_ZnS;
  acc.sumOneMinusCosThetaEncounterEffective += event.sum_one_minus_cos_theta_encounter_effective;
  acc.sumOneMinusCosThetaEncounterEffectiveBN += event.sum_one_minus_cos_theta_encounter_effective_BN;
  acc.sumOneMinusCosThetaEncounterEffectiveZnS += event.sum_one_minus_cos_theta_encounter_effective_ZnS;
  acc.sumCos2ThetaEncounter += event.sum_cos2_theta_encounter;
  acc.sumCos2ThetaEncounterBN += event.sum_cos2_theta_encounter_BN;
  acc.sumCos2ThetaEncounterZnS += event.sum_cos2_theta_encounter_ZnS;
  acc.sumCos2ThetaEncounterEffective += event.sum_cos2_theta_encounter_effective;
  acc.sumCos2ThetaEncounterEffectiveBN += event.sum_cos2_theta_encounter_effective_BN;
  acc.sumCos2ThetaEncounterEffectiveZnS += event.sum_cos2_theta_encounter_effective_ZnS;
  for (std::size_t i = 0; i < acc.phaseFunctionCounts.size(); ++i)
    acc.phaseFunctionCounts[i] += event.phase_function_histogram[i];

  acc.totalParticleScatter += event.num_particle_scatter;
  acc.totalParticleScatterBN += event.num_particle_scatter_BN;
  acc.totalParticleScatterZnS += event.num_particle_scatter_ZnS;
  acc.totalRealScatter += event.num_real_scatter;
  acc.totalBulkScatter += event.num_bulk_scatter;
  acc.totalBoundaryScatter += event.num_boundary_scatter;
  acc.totalBoundaryScatterBN += event.num_boundary_scatter_BN;
  acc.totalBoundaryScatterZnS += event.num_boundary_scatter_ZnS;
  acc.totalMaterialBoundary += event.num_material_boundary;
  acc.totalCompleteEncounter += event.num_complete_encounter_total;
  acc.totalCompleteEncounterBN += event.num_complete_encounter_BN;
  acc.totalCompleteEncounterZnS += event.num_complete_encounter_ZnS;
  acc.totalSurfaceReflectionEncounter += event.num_surface_reflection_encounter;
  acc.totalIncompleteInitialParticleExit += event.num_incomplete_initial_particle_exit;
  acc.totalCensoredParticleEncounter += event.num_censored_particle_encounter;
  acc.totalInconsistentEncounterState += event.num_inconsistent_encounter_state;
  acc.totalParticleToParticleBoundary += event.num_particle_to_particle_boundary;
  acc.totalUnknownParticleReflection += event.num_unknown_particle_reflection;
  acc.totalOuterBoundaryHits += event.num_outer_boundary_hits;
  acc.totalOuterBoundaryReentrySuccess += event.num_outer_boundary_reentry_success;
  acc.totalOuterBoundaryReentryFailed += event.num_outer_boundary_reentry_failed;
  acc.totalOuterBoundaryEscape += event.num_outer_boundary_escape;
  acc.totalOuterBoundaryMaxReentryStop += event.num_outer_boundary_max_reentry_stop;
  acc.totalOuterBoundaryReentryAlgorithmFailed +=
      event.num_outer_boundary_reentry_algorithm_failed;
  acc.totalOuterBoundaryPeriodicPhaseMismatch +=
      event.num_outer_boundary_periodic_phase_mismatch;
  acc.totalOuterBoundaryFresnelReflection += event.num_outer_boundary_fresnel_reflection;
  acc.totalOuterBoundaryTotalInternalReflection +=
      event.num_outer_boundary_total_internal_reflection;
  acc.totalOuterBoundaryRefraction += event.num_outer_boundary_refraction;
  acc.totalOuterBoundaryTransmission += event.num_outer_boundary_transmission;
  acc.totalOuterBoundaryOtherStatus += event.num_outer_boundary_other_status;
  acc.totalReentry += event.num_reentry;
  acc.totalReentryBN += event.num_reentry_BN;
  acc.totalReentryZnS += event.num_reentry_ZnS;
  acc.totalReentryMatrix += event.num_reentry_matrix;
  acc.totalReentryParticleQMu += event.num_reentry_particle_q_mu;
  acc.totalReentryParticleQOnly += event.num_reentry_particle_q_only;
  acc.totalReentryParticleQOnlyFallback +=
      event.num_reentry_particle_q_only_fallback;
  acc.totalReentryParticleVolumeRandom +=
      event.num_reentry_particle_volume_random;
  acc.totalReentryMatrixClearancePortal += event.num_reentry_matrix_clearance_portal;
  acc.totalReentryFallbackSameBin += event.num_reentry_fallback_same_bin;
  acc.totalReentryFallbackAdjacentBin += event.num_reentry_fallback_adjacent_bin;
  acc.totalReentryFallbackAnyBin += event.num_reentry_fallback_any_bin;
  acc.totalReentryFallbackAnyPhaseSameBin +=
      event.num_reentry_fallback_any_phase_same_bin;
  acc.totalReentryFallbackAnyPortal += event.num_reentry_fallback_any_portal;
  acc.totalReentryRandomMatrixDebug += event.num_reentry_random_matrix_debug;
  acc.totalReentryFailed += event.num_reentry_failed;
  acc.sumCosThetaParticleScatter += event.sum_cos_theta_particle;
  acc.sumCosThetaParticleScatterBN += event.sum_cos_theta_particle_BN;
  acc.sumCosThetaParticleScatterZnS += event.sum_cos_theta_particle_ZnS;
  acc.sumCosThetaAllScatter += event.sum_cos_theta;
  acc.sumCosThetaBulkScatter += event.sum_cos_theta_bulk;
  acc.sumCosThetaBoundaryScatter += event.sum_cos_theta_boundary;
  acc.sumCosThetaBoundaryScatterBN += event.sum_cos_theta_boundary_BN;
  acc.sumCosThetaBoundaryScatterZnS += event.sum_cos_theta_boundary_ZnS;
}

void StageDOpticalRunAction::RecordReentryDiagnostic(
    const StageDReentryDiagnosticRecord &record)
{
  if (!fReentryDiagnosticsCsv.is_open() || fConfig == nullptr ||
      !fConfig->stageD_write_reentry_diagnostics ||
      fConfig->stageD_max_diagnostic_rows <= 0 ||
      fReentryDiagnosticRowsWritten >= fConfig->stageD_max_diagnostic_rows ||
      !KeepReentryDiagnostic(
          record, fConfig->stageD_reentry_diagnostics_sampling_rate))
    return;

  fReentryDiagnosticsCsv
      << record.event_id << ","
      << record.reentry_index << ","
      << CsvQuote(record.strategy) << ","
      << CsvQuote(record.fallback_level) << ","
      << CsvQuote(record.exit_phase) << ","
      << CsvQuote(record.entry_phase) << ","
      << record.old_dir.x() << ","
      << record.old_dir.y() << ","
      << record.old_dir.z() << ","
      << record.exit_point.x() / um << ","
      << record.exit_point.y() / um << ","
      << record.exit_point.z() / um << ","
      << record.entry_point.x() / um << ","
      << record.entry_point.y() / um << ","
      << record.entry_point.z() / um << ","
      << record.particle_q_exit << ","
      << record.particle_q_entry << ","
      << record.particle_mu_exit << ","
      << record.particle_mu_entry << ","
      << record.matrix_clearance_exit_um << ","
      << record.matrix_clearance_entry_um << ","
      << CsvQuote(record.matrix_nearest_phase_exit) << ","
      << CsvQuote(record.matrix_nearest_phase_entry) << ","
      << record.matrix_clearance_bin_exit << ","
      << record.matrix_clearance_bin_entry << ","
      << record.trials
      << "\n";
  ++fReentryDiagnosticRowsWritten;
}

void StageDOpticalRunAction::SetReentryPortalSummary(
    const StageDReentryPortalSummary &summary)
{
  fReentryPortalSummary = summary;
}

void StageDOpticalRunAction::WriteSummaryFile() const
{
  std::ofstream fout(fSummaryCsvPath.c_str(), std::ios::out);
  if (!fout)
  {
    G4Exception("StageDOpticalRunAction::WriteSummaryFile",
                "BNZS_D_RUN_003", FatalException,
                ("Failed to open Stage D summary CSV: " + fSummaryCsvPath).c_str());
    return;
  }
  fout << std::setprecision(15);

  const auto &acc = fAccumulator;
  const G4long nPhotons = acc.nPhotons;
  const G4long nAbsorbed = acc.nAbsorbed;
  const G4long nAbsorbedBN = acc.nAbsorbedBN;
  const G4long nAbsorbedZnS = acc.nAbsorbedZnS;
  const G4long nAbsorbedMatrix = acc.nAbsorbedMatrix;
  const G4long nLost = acc.nLost;
  const G4long nReentryFailed = acc.nReentryFailed;
  const G4long nEscapedDebug = acc.nEscapedDebug;
  const G4long nPeriodicPhaseMismatch = acc.nPeriodicPhaseMismatch;
  const G4long nPeriodicGeometryRequired = acc.nPeriodicGeometryRequired;
  const G4long nDetected = acc.nDetected;
  const G4long nMaxReentry = acc.nMaxReentry;
  const G4long nMaxSteps = acc.nMaxSteps;
  const G4long nMaxPathLength = acc.nMaxPathLength;
  const G4long nTargetPrimaryScatter = acc.nTargetPrimaryScatter;
  const G4double totalPathLengthUm = acc.totalPathLengthUm;
  const G4double totalPathLengthBNUm = acc.totalPathLengthBNUm;
  const G4double totalPathLengthZnSUm = acc.totalPathLengthZnSUm;
  const G4double totalPathLengthMatrixUm = acc.totalPathLengthMatrixUm;
  const G4double totalPathLengthWorldUm = acc.totalPathLengthWorldUm;
  const G4long totalEncounter = acc.totalEncounter;
  const G4long totalEncounterBN = acc.totalEncounterBN;
  const G4long totalEncounterZnS = acc.totalEncounterZnS;
  const G4long totalEncounterEffective = acc.totalEncounterEffective;
  const G4long totalEncounterEffectiveBN = acc.totalEncounterEffectiveBN;
  const G4long totalEncounterEffectiveZnS = acc.totalEncounterEffectiveZnS;
  const G4double sumCosThetaEncounter = acc.sumCosThetaEncounter;
  const G4double sumCosThetaEncounterBN = acc.sumCosThetaEncounterBN;
  const G4double sumCosThetaEncounterZnS = acc.sumCosThetaEncounterZnS;
  const G4double sumCosThetaEncounterEffective = acc.sumCosThetaEncounterEffective;
  const G4double sumCosThetaEncounterEffectiveBN = acc.sumCosThetaEncounterEffectiveBN;
  const G4double sumCosThetaEncounterEffectiveZnS = acc.sumCosThetaEncounterEffectiveZnS;
  const G4double sumOneMinusCosThetaEncounter = acc.sumOneMinusCosThetaEncounter;
  const G4double sumOneMinusCosThetaEncounterBN = acc.sumOneMinusCosThetaEncounterBN;
  const G4double sumOneMinusCosThetaEncounterZnS = acc.sumOneMinusCosThetaEncounterZnS;
  const G4double sumOneMinusCosThetaEncounterEffective =
      acc.sumOneMinusCosThetaEncounterEffective;
  const G4double sumOneMinusCosThetaEncounterEffectiveBN =
      acc.sumOneMinusCosThetaEncounterEffectiveBN;
  const G4double sumOneMinusCosThetaEncounterEffectiveZnS =
      acc.sumOneMinusCosThetaEncounterEffectiveZnS;
  const G4double sumCos2ThetaEncounter = acc.sumCos2ThetaEncounter;
  const G4double sumCos2ThetaEncounterBN = acc.sumCos2ThetaEncounterBN;
  const G4double sumCos2ThetaEncounterZnS = acc.sumCos2ThetaEncounterZnS;
  const G4double sumCos2ThetaEncounterEffective = acc.sumCos2ThetaEncounterEffective;
  const G4double sumCos2ThetaEncounterEffectiveBN = acc.sumCos2ThetaEncounterEffectiveBN;
  const G4double sumCos2ThetaEncounterEffectiveZnS = acc.sumCos2ThetaEncounterEffectiveZnS;
  const G4long totalParticleScatter = acc.totalParticleScatter;
  const G4long totalParticleScatterBN = acc.totalParticleScatterBN;
  const G4long totalParticleScatterZnS = acc.totalParticleScatterZnS;
  const G4long totalRealScatter = acc.totalRealScatter;
  const G4long totalBulkScatter = acc.totalBulkScatter;
  const G4long totalBoundaryScatter = acc.totalBoundaryScatter;
  const G4long totalBoundaryScatterBN = acc.totalBoundaryScatterBN;
  const G4long totalBoundaryScatterZnS = acc.totalBoundaryScatterZnS;
  const G4long totalMaterialBoundary = acc.totalMaterialBoundary;
  const G4long totalCompleteEncounter = acc.totalCompleteEncounter;
  const G4long totalCompleteEncounterBN = acc.totalCompleteEncounterBN;
  const G4long totalCompleteEncounterZnS = acc.totalCompleteEncounterZnS;
  const G4long totalSurfaceReflectionEncounter = acc.totalSurfaceReflectionEncounter;
  const G4long totalIncompleteInitialParticleExit =
      acc.totalIncompleteInitialParticleExit;
  const G4long totalCensoredParticleEncounter = acc.totalCensoredParticleEncounter;
  const G4long totalInconsistentEncounterState = acc.totalInconsistentEncounterState;
  const G4long totalParticleToParticleBoundary = acc.totalParticleToParticleBoundary;
  const G4long totalUnknownParticleReflection = acc.totalUnknownParticleReflection;
  const G4long totalOuterBoundaryHits = acc.totalOuterBoundaryHits;
  const G4long totalOuterBoundaryReentrySuccess = acc.totalOuterBoundaryReentrySuccess;
  const G4long totalOuterBoundaryReentryFailed = acc.totalOuterBoundaryReentryFailed;
  const G4long totalOuterBoundaryEscape = acc.totalOuterBoundaryEscape;
  const G4long totalOuterBoundaryMaxReentryStop = acc.totalOuterBoundaryMaxReentryStop;
  const G4long totalOuterBoundaryReentryAlgorithmFailed =
      acc.totalOuterBoundaryReentryAlgorithmFailed;
  const G4long totalOuterBoundaryPeriodicPhaseMismatch =
      acc.totalOuterBoundaryPeriodicPhaseMismatch;
  const G4long totalOuterBoundaryFresnelReflection =
      acc.totalOuterBoundaryFresnelReflection;
  const G4long totalOuterBoundaryTotalInternalReflection =
      acc.totalOuterBoundaryTotalInternalReflection;
  const G4long totalOuterBoundaryRefraction = acc.totalOuterBoundaryRefraction;
  const G4long totalOuterBoundaryTransmission = acc.totalOuterBoundaryTransmission;
  const G4long totalOuterBoundaryOtherStatus = acc.totalOuterBoundaryOtherStatus;
  const G4long totalReentry = acc.totalReentry;
  const G4long totalReentryBN = acc.totalReentryBN;
  const G4long totalReentryZnS = acc.totalReentryZnS;
  const G4long totalReentryMatrix = acc.totalReentryMatrix;
  const G4long totalReentryParticleQMu = acc.totalReentryParticleQMu;
  const G4long totalReentryParticleQOnly = acc.totalReentryParticleQOnly;
  const G4long totalReentryParticleQOnlyFallback =
      acc.totalReentryParticleQOnlyFallback;
  const G4long totalReentryParticleVolumeRandom =
      acc.totalReentryParticleVolumeRandom;
  const G4long totalReentryMatrixClearancePortal =
      acc.totalReentryMatrixClearancePortal;
  const G4long totalReentryFallbackSameBin = acc.totalReentryFallbackSameBin;
  const G4long totalReentryFallbackAdjacentBin = acc.totalReentryFallbackAdjacentBin;
  const G4long totalReentryFallbackAnyBin = acc.totalReentryFallbackAnyBin;
  const G4long totalReentryFallbackAnyPhaseSameBin =
      acc.totalReentryFallbackAnyPhaseSameBin;
  const G4long totalReentryFallbackAnyPortal = acc.totalReentryFallbackAnyPortal;
  const G4long totalReentryRandomMatrixDebug = acc.totalReentryRandomMatrixDebug;
  const G4long totalReentryFailed = acc.totalReentryFailed;
  const G4double sumCosThetaParticleScatter = acc.sumCosThetaParticleScatter;
  const G4double sumCosThetaParticleScatterBN = acc.sumCosThetaParticleScatterBN;
  const G4double sumCosThetaParticleScatterZnS = acc.sumCosThetaParticleScatterZnS;
  const G4double sumCosThetaAllScatter = acc.sumCosThetaAllScatter;
  const G4double sumCosThetaBulkScatter = acc.sumCosThetaBulkScatter;
  const G4double sumCosThetaBoundaryScatter = acc.sumCosThetaBoundaryScatter;
  const G4double sumCosThetaBoundaryScatterBN = acc.sumCosThetaBoundaryScatterBN;
  const G4double sumCosThetaBoundaryScatterZnS = acc.sumCosThetaBoundaryScatterZnS;
  const G4double nPhotonsD = (nPhotons > 0) ? static_cast<G4double>(nPhotons) : 1.0;
  const G4double totalMediumPathLengthUm =
      totalPathLengthBNUm + totalPathLengthZnSUm + totalPathLengthMatrixUm;
  const G4double muACount = (totalMediumPathLengthUm > 0.0)
                                ? static_cast<G4double>(nAbsorbed) / totalMediumPathLengthUm
                                : 0.0;
  const G4double muAExpected =
      (totalMediumPathLengthUm > 0.0 && fConfig != nullptr)
          ? ((totalPathLengthBNUm / std::max(1.0e-12, fConfig->opticalBnAbsLengthUm)) +
             (totalPathLengthZnSUm / std::max(1.0e-12, fConfig->opticalZnsAbsLengthUm)) +
             (totalPathLengthMatrixUm / std::max(1.0e-12, fConfig->opticalMatrixAbsLengthUm))) /
                totalMediumPathLengthUm
          : 0.0;
  const G4double muABNCount = (totalPathLengthBNUm > 0.0)
                                  ? static_cast<G4double>(nAbsorbedBN) / totalPathLengthBNUm
                                  : 0.0;
  const G4double muAZnSCount = (totalPathLengthZnSUm > 0.0)
                                   ? static_cast<G4double>(nAbsorbedZnS) / totalPathLengthZnSUm
                                   : 0.0;
  const G4double muAMatrixCount = (totalPathLengthMatrixUm > 0.0)
                                      ? static_cast<G4double>(nAbsorbedMatrix) / totalPathLengthMatrixUm
                                      : 0.0;

  const G4bool useThresholdedEncounterMetric =
      UseThresholdedEncounterMetric(fConfig);
  const G4long primaryEncounterCount =
      useThresholdedEncounterMetric ? totalEncounterEffective : totalEncounter;
  const G4long primaryEncounterCountBN =
      useThresholdedEncounterMetric ? totalEncounterEffectiveBN : totalEncounterBN;
  const G4long primaryEncounterCountZnS =
      useThresholdedEncounterMetric ? totalEncounterEffectiveZnS : totalEncounterZnS;
  const G4double nan = std::numeric_limits<G4double>::quiet_NaN();
  const G4double primarySumCosThetaEncounter =
      useThresholdedEncounterMetric ? sumCosThetaEncounterEffective : sumCosThetaEncounter;
  const G4double primarySumCosThetaEncounterBN =
      useThresholdedEncounterMetric ? sumCosThetaEncounterEffectiveBN : sumCosThetaEncounterBN;
  const G4double primarySumCosThetaEncounterZnS =
      useThresholdedEncounterMetric ? sumCosThetaEncounterEffectiveZnS : sumCosThetaEncounterZnS;
  const G4double primarySumOneMinusCosThetaEncounter =
      useThresholdedEncounterMetric ? sumOneMinusCosThetaEncounterEffective : sumOneMinusCosThetaEncounter;
  const G4double primarySumCos2ThetaEncounter =
      useThresholdedEncounterMetric ? sumCos2ThetaEncounterEffective : sumCos2ThetaEncounter;

  const G4double muSEncounter = (totalMediumPathLengthUm > 0.0)
                                    ? static_cast<G4double>(primaryEncounterCount) / totalMediumPathLengthUm
                                    : 0.0;
  const G4double muSEncounterBN = (totalMediumPathLengthUm > 0.0)
                                      ? static_cast<G4double>(primaryEncounterCountBN) / totalMediumPathLengthUm
                                      : 0.0;
  const G4double muSEncounterZnS = (totalMediumPathLengthUm > 0.0)
                                       ? static_cast<G4double>(primaryEncounterCountZnS) / totalMediumPathLengthUm
                                       : 0.0;
  const G4double muSEncounterRaw = (totalMediumPathLengthUm > 0.0)
                                    ? static_cast<G4double>(totalEncounter) / totalMediumPathLengthUm
                                    : 0.0;
  const G4double muSEncounterRawBN = (totalMediumPathLengthUm > 0.0)
                                      ? static_cast<G4double>(totalEncounterBN) / totalMediumPathLengthUm
                                      : 0.0;
  const G4double muSEncounterRawZnS = (totalMediumPathLengthUm > 0.0)
                                       ? static_cast<G4double>(totalEncounterZnS) / totalMediumPathLengthUm
                                       : 0.0;
  const G4double g1Encounter = (primaryEncounterCount > 0)
                                   ? primarySumCosThetaEncounter / static_cast<G4double>(primaryEncounterCount)
                                   : nan;
  const G4double g1EncounterBN = (primaryEncounterCountBN > 0)
                                     ? primarySumCosThetaEncounterBN / static_cast<G4double>(primaryEncounterCountBN)
                                     : nan;
  const G4double g1EncounterZnS = (primaryEncounterCountZnS > 0)
                                      ? primarySumCosThetaEncounterZnS / static_cast<G4double>(primaryEncounterCountZnS)
                                      : nan;
  const G4double g1EncounterRaw = (totalEncounter > 0)
                                      ? sumCosThetaEncounter / static_cast<G4double>(totalEncounter)
                                      : nan;
  const G4double g1EncounterRawBN = (totalEncounterBN > 0)
                                        ? sumCosThetaEncounterBN / static_cast<G4double>(totalEncounterBN)
                                        : nan;
  const G4double g1EncounterRawZnS = (totalEncounterZnS > 0)
                                         ? sumCosThetaEncounterZnS / static_cast<G4double>(totalEncounterZnS)
                                         : nan;
  const G4double meanCos2Encounter = (primaryEncounterCount > 0)
                                         ? primarySumCos2ThetaEncounter / static_cast<G4double>(primaryEncounterCount)
                                         : nan;
  const G4double meanCos2EncounterRaw = (totalEncounter > 0)
                                            ? sumCos2ThetaEncounter / static_cast<G4double>(totalEncounter)
                                            : nan;
  const G4double g2Encounter = (primaryEncounterCount > 0)
                                   ? 0.5 * (3.0 * meanCos2Encounter - 1.0)
                                   : nan;
  const G4double g2EncounterRaw = (totalEncounter > 0)
                                      ? 0.5 * (3.0 * meanCos2EncounterRaw - 1.0)
                                      : nan;
  const G4double muSPrimeDirectEncounter = (totalMediumPathLengthUm > 0.0)
                                               ? primarySumOneMinusCosThetaEncounter / totalMediumPathLengthUm
                                               : 0.0;
  const G4double muSPrimeDirectEncounterRaw = (totalMediumPathLengthUm > 0.0)
                                                  ? sumOneMinusCosThetaEncounter / totalMediumPathLengthUm
                                                  : 0.0;
  const G4double muSPrimeEncounterFromG = muSEncounter * (1.0 - g1Encounter);

  const G4double muSParticle = (totalPathLengthUm > 0.0)
                                   ? static_cast<G4double>(totalParticleScatter) / totalPathLengthUm
                                   : 0.0;
  const G4double muSBN = (totalPathLengthUm > 0.0)
                             ? static_cast<G4double>(totalParticleScatterBN) / totalPathLengthUm
                             : 0.0;
  const G4double muSZnS = (totalPathLengthUm > 0.0)
                              ? static_cast<G4double>(totalParticleScatterZnS) / totalPathLengthUm
                              : 0.0;
  const G4double muSBoundaryPrimary = (totalPathLengthUm > 0.0)
                                          ? static_cast<G4double>(totalBoundaryScatter) / totalPathLengthUm
                                          : 0.0;
  const G4double muSBoundaryBN = (totalPathLengthUm > 0.0)
                                     ? static_cast<G4double>(totalBoundaryScatterBN) / totalPathLengthUm
                                     : 0.0;
  const G4double muSBoundaryZnS = (totalPathLengthUm > 0.0)
                                      ? static_cast<G4double>(totalBoundaryScatterZnS) / totalPathLengthUm
                                      : 0.0;
  const G4double muSStepTotal = (totalPathLengthUm > 0.0)
                                    ? static_cast<G4double>(totalRealScatter) / totalPathLengthUm
                                    : 0.0;
  const G4double muSBulk = (totalPathLengthUm > 0.0)
                               ? static_cast<G4double>(totalBulkScatter) / totalPathLengthUm
                               : 0.0;
  const G4double muSBoundary = (totalPathLengthUm > 0.0)
                                   ? static_cast<G4double>(totalBoundaryScatter) / totalPathLengthUm
                                   : 0.0;
  const G4double gParticle = (totalParticleScatter > 0)
                                 ? sumCosThetaParticleScatter / static_cast<G4double>(totalParticleScatter)
                                 : 0.0;
  const G4double gParticleBN = (totalParticleScatterBN > 0)
                                   ? sumCosThetaParticleScatterBN / static_cast<G4double>(totalParticleScatterBN)
                                   : 0.0;
  const G4double gParticleZnS = (totalParticleScatterZnS > 0)
                                    ? sumCosThetaParticleScatterZnS / static_cast<G4double>(totalParticleScatterZnS)
                                    : 0.0;
  const G4double gBoundaryPrimary = (totalBoundaryScatter > 0)
                                        ? sumCosThetaBoundaryScatter / static_cast<G4double>(totalBoundaryScatter)
                                        : 0.0;
  const G4double gBoundaryBN = (totalBoundaryScatterBN > 0)
                                   ? sumCosThetaBoundaryScatterBN / static_cast<G4double>(totalBoundaryScatterBN)
                                   : 0.0;
  const G4double gBoundaryZnS = (totalBoundaryScatterZnS > 0)
                                    ? sumCosThetaBoundaryScatterZnS / static_cast<G4double>(totalBoundaryScatterZnS)
                                    : 0.0;
  const G4double gStepRaw = (totalRealScatter > 0)
                                ? sumCosThetaAllScatter / static_cast<G4double>(totalRealScatter)
                                : 0.0;
  const G4double gBulk = (totalBulkScatter > 0)
                             ? sumCosThetaBulkScatter / static_cast<G4double>(totalBulkScatter)
                             : 0.0;
  const G4double gBoundary = (totalBoundaryScatter > 0)
                                 ? sumCosThetaBoundaryScatter / static_cast<G4double>(totalBoundaryScatter)
                                 : 0.0;
  const G4double muSPrimeParticle = muSParticle * (1.0 - gParticle);
  const G4double muSPrimeBoundaryPrimary = muSBoundaryPrimary * (1.0 - gBoundaryPrimary);
  const G4double muSPrimeStepTotal = muSStepTotal * (1.0 - gStepRaw);
  const G4double muSPrimeBulk = muSBulk * (1.0 - gBulk);
  const G4double muSPrimeBoundary = muSBoundary * (1.0 - gBoundary);
  const G4double particleToParticleBoundaryFraction =
      (totalMaterialBoundary > 0)
          ? static_cast<G4double>(totalParticleToParticleBoundary) /
                static_cast<G4double>(totalMaterialBoundary)
          : 0.0;

  if (particleToParticleBoundaryFraction > 1.0e-6)
  {
    G4cout << "[StageDOpticalRunAction] Warning: direct particle-to-particle boundaries detected."
           << " total_particle_to_particle_boundary=" << totalParticleToParticleBoundary
           << " fraction_of_material_boundaries=" << particleToParticleBoundaryFraction
           << ". Matrix-entry/matrix-exit encounter statistics may be incomplete for this geometry."
           << G4endl;
  }

  if (totalOuterBoundaryHits !=
      totalOuterBoundaryReentrySuccess +
          totalOuterBoundaryEscape +
          totalOuterBoundaryMaxReentryStop +
          totalOuterBoundaryReentryAlgorithmFailed +
          totalOuterBoundaryPeriodicPhaseMismatch)
  {
    G4cout << "[StageDOpticalRunAction] Warning: outer boundary diagnostic mismatch."
           << " hits=" << totalOuterBoundaryHits
           << " success=" << totalOuterBoundaryReentrySuccess
           << " escape=" << totalOuterBoundaryEscape
           << " max_reentry=" << totalOuterBoundaryMaxReentryStop
           << " algorithm_failed=" << totalOuterBoundaryReentryAlgorithmFailed
           << " phase_mismatch=" << totalOuterBoundaryPeriodicPhaseMismatch
           << G4endl;
  }

  if (fConfig != nullptr &&
      (fConfig->stageD_boundary_mode == "same_phase_reentry" ||
       fConfig->stageD_boundary_mode == "periodic_wrap") &&
      (totalOuterBoundaryReentryAlgorithmFailed > 0 ||
       totalOuterBoundaryPeriodicPhaseMismatch > 0))
  {
    G4cout << "[StageDOpticalRunAction] Warning: outer boundary re-entry failures detected."
           << " algorithm_failed=" << totalOuterBoundaryReentryAlgorithmFailed
           << " phase_mismatch=" << totalOuterBoundaryPeriodicPhaseMismatch
           << ". Production periodic runs should have zero failures."
           << G4endl;
  }

  if (fConfig != nullptr &&
      (fConfig->stageD_boundary_mode == "same_phase_reentry" ||
       fConfig->stageD_boundary_mode == "periodic_wrap") &&
      std::abs(totalPathLengthWorldUm) > 1.0e-6)
  {
    G4cout << "[StageDOpticalRunAction] Warning: nonzero World path length in periodic boundary mode."
           << " path_length_World_um=" << totalPathLengthWorldUm
           << G4endl;
  }

  fout
      << "ratio,"
      << "placement_file,"
      << "run_id,"
      << "random_seed_0,"
      << "random_seed_1,"
      << "configuration_tag,"
      << "param_model,"
      << "primary_scatter_metric,"
      << "source_mode,"
      << "boundary_mode,"
      << "reentry_mode,"
      << "particle_reentry_mode,"
      << "matrix_reentry_mode,"
      << "wavelength_nm,"
      << "n_photons,"
      << "num_absorbed_total,"
      << "num_absorbed_BN,"
      << "num_absorbed_ZnS,"
      << "num_absorbed_Matrix,"
      << "n_absorbed,"
      << "absorbed_fraction,"
      << "n_lost,"
      << "lost_fraction,"
      << "n_reentry_failed,"
      << "n_escaped_debug,"
      << "n_periodic_phase_mismatch,"
      << "n_periodic_geometry_required,"
      << "n_detected,"
      << "n_max_reentry,"
      << "n_max_steps,"
      << "n_max_path_length,"
      << "n_target_primary_scatter,"
      << "total_path_length_um,"
      << "path_length_BN_um,"
      << "path_length_ZnS_um,"
      << "path_length_Matrix_um,"
      << "path_length_World_um,"
      << "mean_path_length_um,"
      << "total_medium_path_length_um,"
      << "path_length_closure_error_um,"
      << "num_encounter_total,"
      << "num_encounter_BN,"
      << "num_encounter_ZnS,"
      << "num_encounter_effective_total,"
      << "num_encounter_effective_BN,"
      << "num_encounter_effective_ZnS,"
      << "mu_a_count_per_um,"
      << "mu_a_expected_per_um,"
      << "mu_a_BN_count_per_um,"
      << "mu_a_ZnS_count_per_um,"
      << "mu_a_Matrix_count_per_um,"
      << "mu_s_per_um,"
      << "mu_s_BN_per_um,"
      << "mu_s_ZnS_per_um,"
      << "mu_s_prime_direct_per_um,"
      << "mu_s_prime_from_g_per_um,"
      << "g1,"
      << "g2,"
      << "phase_function_file,"
      << "total_particle_scatter_legacy,"
      << "mean_num_particle_scatter_legacy,"
      << "total_particle_scatter_BN_legacy,"
      << "total_particle_scatter_ZnS_legacy,"
      << "total_real_scatter_debug,"
      << "mean_num_real_scatter_debug,"
      << "total_bulk_scatter_debug,"
      << "mean_num_bulk_scatter_debug,"
      << "total_boundary_scatter_debug,"
      << "mean_num_boundary_scatter_debug,"
      << "total_complete_encounter,"
      << "total_complete_encounter_BN,"
      << "total_complete_encounter_ZnS,"
      << "total_surface_reflection_encounter,"
      << "total_incomplete_initial_particle_exit,"
      << "total_censored_particle_encounter,"
      << "total_inconsistent_encounter_state,"
      << "total_particle_to_particle_boundary,"
      << "particle_to_particle_boundary_fraction,"
      << "total_unknown_particle_reflection,"
      << "total_outer_boundary_hits,"
      << "total_outer_boundary_reentry_success,"
      << "total_outer_boundary_reentry_failed,"
      << "total_outer_boundary_escape,"
      << "total_outer_boundary_max_reentry_stop,"
      << "total_outer_boundary_reentry_algorithm_failed,"
      << "total_outer_boundary_periodic_phase_mismatch,"
      << "total_outer_boundary_fresnel_reflection,"
      << "total_outer_boundary_total_internal_reflection,"
      << "total_outer_boundary_refraction,"
      << "total_outer_boundary_transmission,"
      << "total_outer_boundary_other_status,"
      << "total_reentry,"
      << "mean_num_reentry,"
      << "total_reentry_BN,"
      << "total_reentry_ZnS,"
      << "total_reentry_matrix,"
      << "total_reentry_particle_q_mu,"
      << "total_reentry_particle_q_only,"
      << "total_reentry_particle_q_only_fallback,"
      << "particle_q_only_fallback_fraction,"
      << "total_reentry_particle_volume_random,"
      << "total_reentry_matrix_clearance_portal,"
      << "total_reentry_fallback_same_bin,"
      << "total_reentry_fallback_adjacent_bin,"
      << "total_reentry_fallback_any_bin,"
      << "total_reentry_fallback_any_phase_same_bin,"
      << "total_reentry_fallback_any_portal,"
      << "total_reentry_random_matrix_debug,"
      << "total_reentry_failed,"
      << "matrix_portal_count_total,"
      << "matrix_portal_count_posx,"
      << "matrix_portal_count_negx,"
      << "matrix_portal_count_posy,"
      << "matrix_portal_count_negy,"
      << "matrix_portal_count_posz,"
      << "matrix_portal_count_negz,"
      << "matrix_portal_count_bin0,"
      << "matrix_portal_count_bin1,"
      << "matrix_portal_count_bin2,"
      << "matrix_portal_count_bin3,"
      << "mu_a_raw_per_um,"
      << "mu_s_raw_per_um,"
      << "mu_s_particle_raw_per_um,"
      << "mu_s_boundary_primary_raw_per_um,"
      << "mu_s_raw_BN_per_um,"
      << "mu_s_raw_ZnS_per_um,"
      << "mu_s_step_total_raw_per_um,"
      << "mu_s_bulk_raw_per_um,"
      << "mu_s_boundary_raw_per_um,"
      << "g_raw,"
      << "g_particle_raw,"
      << "g_boundary_primary_raw,"
      << "g_raw_BN,"
      << "g_raw_ZnS,"
      << "g_step_total_raw,"
      << "g_bulk_raw,"
      << "g_boundary_raw,"
      << "mu_s_prime_raw_per_um,"
      << "mu_s_prime_particle_raw_per_um,"
      << "mu_s_prime_boundary_primary_raw_per_um,"
      << "mu_s_prime_raw_BN_per_um,"
      << "mu_s_prime_raw_ZnS_per_um,"
      << "mu_s_prime_step_total_raw_per_um,"
      << "mu_s_prime_bulk_raw_per_um,"
      << "mu_s_prime_boundary_raw_per_um,"
      << "optical_params_provided,"
      << "scatter_metric,"
      << "target_primary_scatter,"
      << "matrix_n,"
      << "matrix_abs_um,"
      << "bn_n,"
      << "bn_abs_um,"
      << "zns_n,"
      << "zns_abs_um,"
      << "theta_threshold_deg,"
      << "max_reentry,"
      << "max_steps,"
      << "max_path_length_um"
      << "\n";

  fout
      << CsvQuote(fRatioTag) << ","
      << CsvQuote(fPlacementFile) << ","
      << fRunId << ","
      << fRandomSeed0 << ","
      << fRandomSeed1 << ","
      << CsvQuote(fConfigurationTag) << ","
      << "GO_RVE,"
      << CsvQuote(PrimaryScatterMetricLabel(fConfig)) << ","
      << CsvQuote(fConfig ? fConfig->stageD_source_mode : "") << ","
      << CsvQuote(fConfig ? fConfig->stageD_boundary_mode : "") << ","
      << CsvQuote(fConfig ? fConfig->stageD_reentry_mode : "") << ","
      << CsvQuote(fConfig ? fConfig->stageD_particle_reentry_mode : "") << ","
      << CsvQuote(fConfig ? fConfig->stageD_matrix_reentry_mode : "") << ","
      << (fConfig ? fConfig->stageD_wavelength_nm : 0.0) << ","
      << nPhotons << ","
      << nAbsorbed << ","
      << nAbsorbedBN << ","
      << nAbsorbedZnS << ","
      << nAbsorbedMatrix << ","
      << nAbsorbed << ","
      << static_cast<G4double>(nAbsorbed) / nPhotonsD << ","
      << nLost << ","
      << static_cast<G4double>(nLost) / nPhotonsD << ","
      << nReentryFailed << ","
      << nEscapedDebug << ","
      << nPeriodicPhaseMismatch << ","
      << nPeriodicGeometryRequired << ","
      << nDetected << ","
      << nMaxReentry << ","
      << nMaxSteps << ","
      << nMaxPathLength << ","
      << nTargetPrimaryScatter << ","
      << totalPathLengthUm << ","
      << totalPathLengthBNUm << ","
      << totalPathLengthZnSUm << ","
      << totalPathLengthMatrixUm << ","
      << totalPathLengthWorldUm << ","
      << totalPathLengthUm / nPhotonsD << ","
      << totalMediumPathLengthUm << ","
      << (totalPathLengthUm - totalMediumPathLengthUm - totalPathLengthWorldUm) << ","
      << totalEncounter << ","
      << totalEncounterBN << ","
      << totalEncounterZnS << ","
      << totalEncounterEffective << ","
      << totalEncounterEffectiveBN << ","
      << totalEncounterEffectiveZnS << ","
      << muACount << ","
      << muAExpected << ","
      << muABNCount << ","
      << muAZnSCount << ","
      << muAMatrixCount << ","
      << muSEncounter << ","
      << muSEncounterBN << ","
      << muSEncounterZnS << ","
      << muSPrimeDirectEncounter << ","
      << muSPrimeEncounterFromG << ","
      << g1Encounter << ","
      << g2Encounter << ","
      << CsvQuote(std::filesystem::path(fPhaseFunctionCsvPath).filename().string()) << ","
      << totalParticleScatter << ","
      << static_cast<G4double>(totalParticleScatter) / nPhotonsD << ","
      << totalParticleScatterBN << ","
      << totalParticleScatterZnS << ","
      << totalRealScatter << ","
      << static_cast<G4double>(totalRealScatter) / nPhotonsD << ","
      << totalBulkScatter << ","
      << static_cast<G4double>(totalBulkScatter) / nPhotonsD << ","
      << totalBoundaryScatter << ","
      << static_cast<G4double>(totalBoundaryScatter) / nPhotonsD << ","
      << totalCompleteEncounter << ","
      << totalCompleteEncounterBN << ","
      << totalCompleteEncounterZnS << ","
      << totalSurfaceReflectionEncounter << ","
      << totalIncompleteInitialParticleExit << ","
      << totalCensoredParticleEncounter << ","
      << totalInconsistentEncounterState << ","
      << totalParticleToParticleBoundary << ","
      << particleToParticleBoundaryFraction << ","
      << totalUnknownParticleReflection << ","
      << totalOuterBoundaryHits << ","
      << totalOuterBoundaryReentrySuccess << ","
      << totalOuterBoundaryReentryFailed << ","
      << totalOuterBoundaryEscape << ","
      << totalOuterBoundaryMaxReentryStop << ","
      << totalOuterBoundaryReentryAlgorithmFailed << ","
      << totalOuterBoundaryPeriodicPhaseMismatch << ","
      << totalOuterBoundaryFresnelReflection << ","
      << totalOuterBoundaryTotalInternalReflection << ","
      << totalOuterBoundaryRefraction << ","
      << totalOuterBoundaryTransmission << ","
      << totalOuterBoundaryOtherStatus << ","
      << totalReentry << ","
      << static_cast<G4double>(totalReentry) / nPhotonsD << ","
      << totalReentryBN << ","
      << totalReentryZnS << ","
      << totalReentryMatrix << ","
      << totalReentryParticleQMu << ","
      << totalReentryParticleQOnly << ","
      << totalReentryParticleQOnlyFallback << ","
      << (totalReentryParticleQMu + totalReentryParticleQOnlyFallback > 0
              ? static_cast<G4double>(totalReentryParticleQOnlyFallback) /
                    static_cast<G4double>(totalReentryParticleQMu +
                                          totalReentryParticleQOnlyFallback)
              : 0.0)
      << ","
      << totalReentryParticleVolumeRandom << ","
      << totalReentryMatrixClearancePortal << ","
      << totalReentryFallbackSameBin << ","
      << totalReentryFallbackAdjacentBin << ","
      << totalReentryFallbackAnyBin << ","
      << totalReentryFallbackAnyPhaseSameBin << ","
      << totalReentryFallbackAnyPortal << ","
      << totalReentryRandomMatrixDebug << ","
      << totalReentryFailed << ","
      << fReentryPortalSummary.total_portal_count << ","
      << fReentryPortalSummary.portal_count_by_face[0] << ","
      << fReentryPortalSummary.portal_count_by_face[1] << ","
      << fReentryPortalSummary.portal_count_by_face[2] << ","
      << fReentryPortalSummary.portal_count_by_face[3] << ","
      << fReentryPortalSummary.portal_count_by_face[4] << ","
      << fReentryPortalSummary.portal_count_by_face[5] << ","
      << fReentryPortalSummary.portal_count_by_bin[0] << ","
      << fReentryPortalSummary.portal_count_by_bin[1] << ","
      << fReentryPortalSummary.portal_count_by_bin[2] << ","
      << fReentryPortalSummary.portal_count_by_bin[3] << ","
      << muACount << ","
      << muSEncounterRaw << ","
      << muSParticle << ","
      << muSBoundaryPrimary << ","
      << muSEncounterRawBN << ","
      << muSEncounterRawZnS << ","
      << muSStepTotal << ","
      << muSBulk << ","
      << muSBoundary << ","
      << g1EncounterRaw << ","
      << gParticle << ","
      << gBoundaryPrimary << ","
      << g1EncounterRawBN << ","
      << g1EncounterRawZnS << ","
      << gStepRaw << ","
      << gBulk << ","
      << gBoundary << ","
      << muSPrimeDirectEncounterRaw << ","
      << muSPrimeParticle << ","
      << muSPrimeBoundaryPrimary << ","
      << (muSEncounterRawBN * (1.0 - g1EncounterRawBN)) << ","
      << (muSEncounterRawZnS * (1.0 - g1EncounterRawZnS)) << ","
      << muSPrimeStepTotal << ","
      << muSPrimeBulk << ","
      << muSPrimeBoundary << ","
      << ((fConfig && fConfig->opticalParamsProvided) ? 1 : 0) << ","
      << CsvQuote(fConfig ? fConfig->stageD_scatter_metric : "") << ","
      << (fConfig ? fConfig->stageD_target_primary_scatter : 0) << ","
      << (fConfig ? fConfig->opticalMatrixRIndex : 0.0) << ","
      << (fConfig ? fConfig->opticalMatrixAbsLengthUm : 0.0) << ","
      << (fConfig ? fConfig->opticalBnRIndex : 0.0) << ","
      << (fConfig ? fConfig->opticalBnAbsLengthUm : 0.0) << ","
      << (fConfig ? fConfig->opticalZnsRIndex : 0.0) << ","
      << (fConfig ? fConfig->opticalZnsAbsLengthUm : 0.0) << ","
      << (fConfig ? fConfig->stageD_theta_threshold_deg : 0.0) << ","
      << (fConfig ? fConfig->stageD_max_reentry : 0) << ","
      << (fConfig ? fConfig->stageD_max_steps : 0) << ","
      << (fConfig ? fConfig->stageD_max_path_length_um : 0.0)
      << "\n";
  fout.close();

  std::error_code copyEc;
  std::filesystem::copy_file(
      fSummaryCsvPath,
      (std::filesystem::path(fOutputDir) / "optical_homogenization_summary.csv").string(),
      std::filesystem::copy_options::overwrite_existing,
      copyEc);
}

void StageDOpticalRunAction::WritePhaseFunctionFile() const
{
  std::ofstream fout(fPhaseFunctionCsvPath.c_str(), std::ios::out);
  if (!fout)
  {
    G4Exception("StageDOpticalRunAction::WritePhaseFunctionFile",
                "BNZS_D_RUN_004", FatalException,
                ("Failed to open Stage D phase function CSV: " + fPhaseFunctionCsvPath).c_str());
    return;
  }
  fout << std::setprecision(15);

  const auto &counts = fAccumulator.phaseFunctionCounts;
  G4long totalCount = 0;
  for (const auto count : counts)
    totalCount += count;

  fout << "lambda_nm,bin_id,cos_theta_min,cos_theta_max,count,probability,probability_density\n";
  const G4double binWidth = 2.0 / static_cast<G4double>(counts.size());
  for (std::size_t i = 0; i < counts.size(); ++i)
  {
    const G4double cosMin = -1.0 + static_cast<G4double>(i) * binWidth;
    const G4double cosMax = cosMin + binWidth;
    const G4double probability =
        (totalCount > 0) ? static_cast<G4double>(counts[i]) / static_cast<G4double>(totalCount) : 0.0;
    const G4double density = (binWidth > 0.0) ? (probability / binWidth) : 0.0;
    fout << (fConfig ? fConfig->stageD_wavelength_nm : 0.0) << ","
         << i << ","
         << cosMin << ","
         << cosMax << ","
         << counts[i] << ","
         << probability << ","
         << density << "\n";
  }
  fout.close();

  std::error_code copyEc;
  std::filesystem::copy_file(
      fPhaseFunctionCsvPath,
      (std::filesystem::path(fOutputDir) / "phase_function_lambda.csv").string(),
      std::filesystem::copy_options::overwrite_existing,
      copyEc);
}

void StageDOpticalRunAction::EndOfRunAction(const G4Run *run)
{
  (void)run;

  if (fEventsCsv.is_open())
    fEventsCsv.close();
  if (fReentryDiagnosticsCsv.is_open())
    fReentryDiagnosticsCsv.close();

  WriteSummaryFile();
  WritePhaseFunctionFile();

  G4cout << "[StageDOpticalRunAction] End run"
         << "\n  events csv  = " << fEventsCsvPath
         << "\n  summary csv = " << fSummaryCsvPath
         << "\n  phase csv   = " << fPhaseFunctionCsvPath
         << "\n  n photons   = " << fAccumulator.nPhotons
         << G4endl;
}
