#include "AnalysisConfig.hh"

#include <cstdlib>
#include <cctype>
#include <string>
#include <algorithm>
#include <filesystem>
#include <vector>

namespace
{
  std::filesystem::path DetectProjectRoot()
  {
    namespace fs = std::filesystem;

    std::error_code ec;
    const fs::path cwd = fs::current_path(ec);
    if (ec)
      return fs::path(".");

    std::vector<fs::path> candidates;
    candidates.push_back(cwd);
    if (cwd.filename() == "build" && cwd.has_parent_path())
      candidates.push_back(cwd.parent_path());
    if (cwd.has_parent_path())
      candidates.push_back(cwd.parent_path());

    for (const auto &candidate : candidates)
    {
      if (candidate.empty())
        continue;
      if (fs::exists(candidate / "CMakeLists.txt", ec) && !ec)
        return candidate;
      ec.clear();
      if (fs::exists(candidate / "README.md", ec) && !ec)
        return candidate;
      ec.clear();
    }

    return (cwd.filename() == "build" && cwd.has_parent_path()) ? cwd.parent_path() : cwd;
  }

  std::string ToLowerCopy(std::string s)
  {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char ch)
                   { return static_cast<char>(std::tolower(ch)); });
    return s;
  }

  bool IsTruthyEnv(const char *value)
  {
    if (value == nullptr)
      return false;

    const std::string v = ToLowerCopy(value);
    return v == "1" || v == "true" || v == "yes" || v == "on";
  }

  std::string NormalizeStageDScatterMetric(std::string raw)
  {
    raw.erase(raw.begin(), std::find_if(raw.begin(), raw.end(), [](unsigned char ch)
                                        { return !std::isspace(ch); }));
    raw.erase(std::find_if(raw.rbegin(), raw.rend(), [](unsigned char ch)
                           { return !std::isspace(ch); })
                 .base(),
             raw.end());
    const std::string value = ToLowerCopy(raw);
    if (value == "particle_encounter_no_threshold" ||
        value == "particleencounternothreshold" ||
        value == "encounter" ||
        value == "particle_encounter")
      return "particle_encounter_no_threshold";
    if (value == "angle_threshold" ||
        value == "particle_encounter_angle_threshold" ||
        value == "particleencounteranglethreshold" ||
        value == "step_angle_threshold" ||
        value == "stepanglethreshold")
      return "particle_encounter_angle_threshold";
    if (value == "boundary_deflection" || value == "boundarydeflection")
      return "boundary_deflection";
    if (value == "particle_exit_deflection" || value == "particleexitdeflection")
      return "particle_exit_deflection";
    return raw;
  }

  bool TryParseRatioFolderName(const std::string &name, double &bnWt, double &znsWt)
  {
    const std::size_t dashPos = name.find('-');
    if (dashPos == std::string::npos)
      return false;

    try
    {
      bnWt = std::stod(name.substr(0, dashPos));
      znsWt = std::stod(name.substr(dashPos + 1));
    }
    catch (...)
    {
      return false;
    }

    return bnWt > 0.0 && znsWt > 0.0;
  }

}

std::filesystem::path AnalysisConfig::ProjectRootPath()
{
  return DetectProjectRoot();
}

std::string AnalysisConfig::PathForRecord(const std::filesystem::path &path)
{
  namespace fs = std::filesystem;

  if (path.empty())
    return "";

  const fs::path projectRoot = ProjectRootPath();
  std::error_code ec;

  fs::path normalized = path;
  if (path.is_absolute())
    normalized = fs::weakly_canonical(path, ec);
  else
    normalized = (projectRoot / path).lexically_normal();

  if (ec)
  {
    ec.clear();
    normalized = path.lexically_normal();
  }

  const fs::path normalizedRoot = fs::weakly_canonical(projectRoot, ec);
  if (ec)
  {
    ec.clear();
    return normalized.generic_string();
  }

  const fs::path relative = normalized.lexically_relative(normalizedRoot);
  if (!relative.empty() && relative.native().find("..") != 0)
    return relative.generic_string();

  return normalized.generic_string();
}

std::string AnalysisConfig::PathForRecord(const std::string &path)
{
  if (path.empty())
    return "";
  return PathForRecord(std::filesystem::path(path));
}

AnalysisConfig::AnalysisConfig()
    : runMode(RunMode::StageB_ReplayAlphaLi),
      patchXY_um(50.0),
      microThickness_um(30.0),
      bnWt(1.0),
      znsWt(2.0),
      useRandomPlacement(true),
      placementFilePath(""),
      placementGeometryMode("pbc_clipped"),
      periodicImagesFilePath(""),
      captureCsvPath(""),
      captureInputDir(""),
      opticalParamsProvided(false),
      opticalMatrixRIndex(1.0),
      opticalMatrixAbsLengthUm(1.0e5),
      opticalBnRIndex(1.92),
      opticalBnAbsLengthUm(3.2e3),
      opticalZnsRIndex(2.47),
      opticalZnsAbsLengthUm(4.1e4),
      stageD_wavelength_nm(450.0),
      stageD_source_mode("uniform_all_phase"),
      stageD_boundary_mode("periodic_wrap"),
      stageD_reentry_mode("state_matched"),
      stageD_particle_reentry_mode("sphere_q_mu"),
      stageD_matrix_reentry_mode("clearance_binned_portal"),
      stageD_scatter_metric("particle_encounter_angle_threshold"),
      stageD_target_primary_scatter(0),
      stageD_theta_threshold_deg(0),
      stageD_max_reentry(50000),
      stageD_max_steps(100000),
      stageD_max_path_length_um(10000.0),
      stageD_output_dir(""),
      stageD_portal_nu(128),
      stageD_portal_nv(128),
      stageD_portal_margin_um(0.0),
      stageD_clearance_bin0_um(0.02),
      stageD_clearance_bin1_um(0.10),
      stageD_clearance_bin2_um(0.50),
      stageD_max_particle_reentry_trials(64),
      stageD_max_portal_fallback_level(4),
      allowThicknessEqualLocalPatch(true)
{
  const char *runModeEnv = std::getenv("BNZS_RUN_MODE");
  if (runModeEnv != nullptr)
  {
    const std::string mode = ToLowerCopy(runModeEnv);
    if (mode == "stagea" || mode == "a" || mode == "stagea_neutronpatch")
    {
      runMode = RunMode::StageA_NeutronPatch;
    }
    else if (mode == "stageb" || mode == "b" || mode == "stageb_replayalphali")
    {
      runMode = RunMode::StageB_ReplayAlphaLi;
    }
    else if (mode == "staged_opticalhomogenization" || mode == "staged" || mode == "d")
    {
      runMode = RunMode::StageD_OpticalHomogenization;
    }
  }

  const char *placementEnv = std::getenv("BNZS_PLACEMENT_FILE");
  if (placementEnv != nullptr && std::string(placementEnv).size() > 0)
  {
    placementFilePath = placementEnv;
    useRandomPlacement = false;
  }

  const char *randomPlacementEnv = std::getenv("BNZS_USE_RANDOM_PLACEMENT");
  if (randomPlacementEnv != nullptr)
  {
    useRandomPlacement = IsTruthyEnv(randomPlacementEnv);
  }

  const char *geometryModeEnv = std::getenv("BNZS_PLACEMENT_GEOMETRY_MODE");
  if (geometryModeEnv != nullptr && std::string(geometryModeEnv).size() > 0)
    placementGeometryMode = geometryModeEnv;

  const char *periodicImagesEnv = std::getenv("BNZS_PBC_IMAGES_CSV");
  if (periodicImagesEnv != nullptr && std::string(periodicImagesEnv).size() > 0)
    periodicImagesFilePath = periodicImagesEnv;

  const char *captureCsvEnv = std::getenv("BNZS_INPUT_CSV");
  if (captureCsvEnv != nullptr && std::string(captureCsvEnv).size() > 0)
  {
    captureCsvPath = captureCsvEnv;

    double dirBnWt = 0.0;
    double dirZnsWt = 0.0;
    const std::string dirName = std::filesystem::path(captureCsvPath).parent_path().filename().string();
    if (TryParseRatioFolderName(dirName, dirBnWt, dirZnsWt))
    {
      bnWt = dirBnWt;
      znsWt = dirZnsWt;
    }
  }

  const char *captureDirEnv = std::getenv("BNZS_INPUT_DIR");
  if (captureDirEnv != nullptr && std::string(captureDirEnv).size() > 0)
  {
    captureInputDir = captureDirEnv;

    double dirBnWt = 0.0;
    double dirZnsWt = 0.0;
    const std::string dirName = std::filesystem::path(captureInputDir).filename().string();
    if (TryParseRatioFolderName(dirName, dirBnWt, dirZnsWt))
    {
      bnWt = dirBnWt;
      znsWt = dirZnsWt;
    }
  }

  auto readOpticalDoubleEnv = [&](const char *name, double &target)
  {
    const char *value = std::getenv(name);
    if (value == nullptr || std::string(value).empty())
      return;
    try
    {
      const double parsed = std::stod(value);
      if (parsed > 0.0)
      {
        target = parsed;
        opticalParamsProvided = true;
      }
    }
    catch (...)
    {
    }
  };

  readOpticalDoubleEnv("BNZS_OPTICAL_MATRIX_RINDEX", opticalMatrixRIndex);
  readOpticalDoubleEnv("BNZS_OPTICAL_MATRIX_ABSLENGTH_UM", opticalMatrixAbsLengthUm);
  readOpticalDoubleEnv("BNZS_OPTICAL_BN_RINDEX", opticalBnRIndex);
  readOpticalDoubleEnv("BNZS_OPTICAL_BN_ABSLENGTH_UM", opticalBnAbsLengthUm);
  readOpticalDoubleEnv("BNZS_OPTICAL_ZNS_RINDEX", opticalZnsRIndex);
  readOpticalDoubleEnv("BNZS_OPTICAL_ZNS_ABSLENGTH_UM", opticalZnsAbsLengthUm);

  readOpticalDoubleEnv("BNZS_STAGED_WAVELENGTH_NM", stageD_wavelength_nm);
  readOpticalDoubleEnv("BNZS_STAGED_THETA_THRESHOLD_DEG", stageD_theta_threshold_deg);
  readOpticalDoubleEnv("BNZS_STAGED_MAX_PATH_LENGTH_UM", stageD_max_path_length_um);

  auto readPositiveIntEnv = [](const char *name, int &target)
  {
    const char *value = std::getenv(name);
    if (value == nullptr || std::string(value).empty())
      return;
    try
    {
      const int parsed = std::stoi(value);
      if (parsed > 0)
        target = parsed;
    }
    catch (...)
    {
    }
  };

  readPositiveIntEnv("BNZS_STAGED_MAX_REENTRY", stageD_max_reentry);
  readPositiveIntEnv("BNZS_STAGED_MAX_STEPS", stageD_max_steps);
  readPositiveIntEnv("BNZS_STAGED_TARGET_PRIMARY_SCATTER", stageD_target_primary_scatter);
  readPositiveIntEnv("BNZS_STAGED_PORTAL_NU", stageD_portal_nu);
  readPositiveIntEnv("BNZS_STAGED_PORTAL_NV", stageD_portal_nv);
  readPositiveIntEnv("BNZS_STAGED_MAX_PARTICLE_REENTRY_TRIALS", stageD_max_particle_reentry_trials);
  readPositiveIntEnv("BNZS_STAGED_MAX_PORTAL_FALLBACK_LEVEL", stageD_max_portal_fallback_level);

  const char *stageDSourceModeEnv = std::getenv("BNZS_STAGED_SOURCE_MODE");
  if (stageDSourceModeEnv != nullptr && std::string(stageDSourceModeEnv).size() > 0)
    stageD_source_mode = stageDSourceModeEnv;

  const char *stageDBoundaryModeEnv = std::getenv("BNZS_STAGED_BOUNDARY_MODE");
  if (stageDBoundaryModeEnv != nullptr && std::string(stageDBoundaryModeEnv).size() > 0)
    stageD_boundary_mode = stageDBoundaryModeEnv;

  const char *stageDReentryModeEnv = std::getenv("BNZS_STAGED_REENTRY_MODE");
  if (stageDReentryModeEnv != nullptr && std::string(stageDReentryModeEnv).size() > 0)
    stageD_reentry_mode = stageDReentryModeEnv;

  const char *stageDParticleReentryModeEnv = std::getenv("BNZS_STAGED_PARTICLE_REENTRY_MODE");
  if (stageDParticleReentryModeEnv != nullptr && std::string(stageDParticleReentryModeEnv).size() > 0)
    stageD_particle_reentry_mode = stageDParticleReentryModeEnv;

  const char *stageDMatrixReentryModeEnv = std::getenv("BNZS_STAGED_MATRIX_REENTRY_MODE");
  if (stageDMatrixReentryModeEnv != nullptr && std::string(stageDMatrixReentryModeEnv).size() > 0)
    stageD_matrix_reentry_mode = stageDMatrixReentryModeEnv;

  const char *stageDScatterMetricEnv = std::getenv("BNZS_STAGED_SCATTER_METRIC");
  if (stageDScatterMetricEnv != nullptr && std::string(stageDScatterMetricEnv).size() > 0)
    stageD_scatter_metric = NormalizeStageDScatterMetric(stageDScatterMetricEnv);

  const char *stageDOutputDirEnv = std::getenv("BNZS_STAGED_OUTPUT_DIR");
  if (stageDOutputDirEnv != nullptr && std::string(stageDOutputDirEnv).size() > 0)
    stageD_output_dir = stageDOutputDirEnv;

  readOpticalDoubleEnv("BNZS_STAGED_PORTAL_MARGIN_UM", stageD_portal_margin_um);
  readOpticalDoubleEnv("BNZS_STAGED_CLEARANCE_BIN0_UM", stageD_clearance_bin0_um);
  readOpticalDoubleEnv("BNZS_STAGED_CLEARANCE_BIN1_UM", stageD_clearance_bin1_um);
  readOpticalDoubleEnv("BNZS_STAGED_CLEARANCE_BIN2_UM", stageD_clearance_bin2_um);

  const char *bnWtEnv = std::getenv("BNZS_BN_WT");
  const char *znsWtEnv = std::getenv("BNZS_ZNS_WT");
  if (bnWtEnv != nullptr && znsWtEnv != nullptr)
  {
    try
    {
      const double envBnWt = std::stod(bnWtEnv);
      const double envZnsWt = std::stod(znsWtEnv);
      if (envBnWt > 0.0 && envZnsWt > 0.0)
      {
        bnWt = envBnWt;
        znsWt = envZnsWt;
      }
    }
    catch (...)
    {
    }
  }
}

AnalysisConfig::~AnalysisConfig() = default;

const char *AnalysisConfig::RunModeName(RunMode mode)
{
  switch (mode)
  {
  case RunMode::StageA_NeutronPatch:
    return "StageA_NeutronPatch";
  case RunMode::StageB_ReplayAlphaLi:
    return "StageB_ReplayAlphaLi";
  case RunMode::StageD_OpticalHomogenization:
    return "StageD_OpticalHomogenization";
  default:
    return "UnknownRunMode";
  }
}
