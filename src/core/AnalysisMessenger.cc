#include "AnalysisMessenger.hh"

#include "AnalysisConfig.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4ios.hh"
#include "G4Exception.hh"

#include <algorithm>
#include <cctype>
#include <string>
#include <sstream>
#include <filesystem>

namespace
{
  std::string Trim(const std::string &s)
  {
    const auto first = std::find_if_not(s.begin(), s.end(),
                                        [](unsigned char ch)
                                        { return std::isspace(ch); });
    const auto last = std::find_if_not(s.rbegin(), s.rend(),
                                       [](unsigned char ch)
                                       { return std::isspace(ch); })
                          .base();
    if (first >= last)
      return "";
    return std::string(first, last);
  }

  std::string ToLowerCopy(std::string s)
  {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char ch)
                   { return static_cast<char>(std::tolower(ch)); });
    return s;
  }

  std::string NormalizeStageDSourceMode(const std::string &raw)
  {
    const std::string value = ToLowerCopy(Trim(raw));
    if (value == "uniform_zns" || value == "uniformzns")
      return "uniform_ZnS";
    if (value == "uniform_all_phase" || value == "uniformallphase")
      return "uniform_all_phase";
    return "";
  }

  std::string NormalizeStageDBoundaryMode(const std::string &raw)
  {
    const std::string value = ToLowerCopy(Trim(raw));
    if (value == "escape")
      return "escape";
    if (value == "periodic_wrap" || value == "periodicwrap")
      return "periodic_wrap";
    if (value == "same_phase_reentry" || value == "samephasereentry")
      return "same_phase_reentry";
    return "";
  }

  std::string NormalizeStageDReentryMode(const std::string &raw)
  {
    const std::string value = ToLowerCopy(Trim(raw));
    if (value == "same_phase_rho_over_r" || value == "samephaserhooverr")
      return "same_phase_rho_over_R";
    if (value == "same_phase_random" || value == "samephaserandom")
      return "same_phase_random";
    if (value == "state_matched" || value == "statematched")
      return "state_matched";
    return "";
  }

  std::string NormalizeStageDParticleReentryMode(const std::string &raw)
  {
    const std::string value = ToLowerCopy(Trim(raw));
    if (value == "sphere_q_mu" || value == "sphereqmu")
      return "sphere_q_mu";
    if (value == "same_phase_rho_over_r" || value == "samephaserhooverr")
      return "same_phase_rho_over_R";
    if (value == "same_phase_random" || value == "samephaserandom")
      return "same_phase_random";
    return "";
  }

  std::string NormalizeStageDMatrixReentryMode(const std::string &raw)
  {
    const std::string value = ToLowerCopy(Trim(raw));
    if (value == "random_matrix" || value == "randommatrix")
      return "random_matrix";
    if (value == "random_matrix_debug" || value == "randommatrixdebug")
      return "random_matrix_debug";
    if (value == "distance_matched_matrix" || value == "distancematchedmatrix")
      return "distance_matched_matrix";
    if (value == "clearance_binned_portal" || value == "clearancebinnedportal")
      return "clearance_binned_portal";
    return "";
  }

  std::string NormalizeStageDScatterMetric(const std::string &raw)
  {
    const std::string value = ToLowerCopy(Trim(raw));
    if (value == "particle_encounter_no_threshold" ||
        value == "particleencounternothreshold" ||
        value == "encounter" ||
        value == "particle_encounter")
      return "particle_encounter_no_threshold";
    if (value == "angle_threshold" ||
        value == "particle_encounter_angle_threshold" ||
        value == "particleencounteranglethreshold")
      return "particle_encounter_angle_threshold";
    return "";
  }

  RunMode ParseRunModeOrThrow(const std::string &raw)
  {
    const std::string v = ToLowerCopy(Trim(raw));

    if (v == "stagea_neutronpatch" || v == "stagea" || v == "a")
      return RunMode::StageA_NeutronPatch;

    if (v == "stageb_replayalphali" || v == "stageb" || v == "b")
      return RunMode::StageB_ReplayAlphaLi;

    if (v == "staged_opticalhomogenization" || v == "staged" || v == "d")
      return RunMode::StageD_OpticalHomogenization;

    G4Exception("AnalysisMessenger::ParseRunModeOrThrow",
                "BNZS_CFG_001", FatalException,
                ("Unknown run mode: " + raw +
                 ". Supported values: StageA_NeutronPatch, StageB_ReplayAlphaLi, StageD_OpticalHomogenization")
                    .c_str());

    return RunMode::StageB_ReplayAlphaLi;
  }

  bool TryParseRatioFolderName(const std::string &name, G4double &bnWt, G4double &znsWt)
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

} // namespace

AnalysisMessenger::AnalysisMessenger(AnalysisConfig *config)
    : G4UImessenger(),
      fConfig(config),
      fCfgDir(nullptr),
      fRunModeCmd(nullptr),
      fCaptureCsvCmd(nullptr),
      fCaptureDirCmd(nullptr),
      fPlacementFileCmd(nullptr),
      fPlacementGeometryModeCmd(nullptr),
      fPeriodicImagesFileCmd(nullptr),
      fUseRandomPlacementCmd(nullptr),
      fAllowThicknessEqualCmd(nullptr),
      fStageDDir(nullptr),
      fStageDWavelengthNmCmd(nullptr),
      fStageDSourceModeCmd(nullptr),
      fStageDBoundaryModeCmd(nullptr),
      fStageDReentryModeCmd(nullptr),
      fStageDMatrixReentryModeCmd(nullptr),
      fStageDThetaThresholdDegCmd(nullptr),
      fStageDMaxReentryCmd(nullptr),
      fStageDMaxStepsCmd(nullptr),
      fStageDMaxPathLengthUmCmd(nullptr),
      fStageDOutputDirCmd(nullptr),
      fStageDWriteReentryDiagnosticsCmd(nullptr),
      fStageDReentryDiagnosticsSamplingRateCmd(nullptr),
      fStageDMaxDiagnosticRowsCmd(nullptr),
      fOpticalParamsCmd(nullptr),
      fWeightRatioCmd(nullptr)
{
  fCfgDir = new G4UIdirectory("/cfg/");
  fCfgDir->SetGuidance("Analysis configuration control.");

  fRunModeCmd = new G4UIcmdWithAString("/cfg/setRunMode", this);
  fRunModeCmd->SetGuidance("Set analysis run mode.");
  fRunModeCmd->SetGuidance("Supported: StageA_NeutronPatch | StageB_ReplayAlphaLi | StageD_OpticalHomogenization");
  fRunModeCmd->SetParameterName("runMode", false);
  fRunModeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fCaptureCsvCmd = new G4UIcmdWithAString("/cfg/setCaptureCsv", this);
  fCaptureCsvCmd->SetGuidance("Set explicit capture CSV path for Stage B replay.");
  fCaptureCsvCmd->SetParameterName("captureCsvPath", false);
  fCaptureCsvCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fCaptureDirCmd = new G4UIcmdWithAString("/cfg/setCaptureDir", this);
  fCaptureDirCmd->SetGuidance("Set explicit capture CSV directory for Stage B replay.");
  fCaptureDirCmd->SetParameterName("captureInputDir", false);
  fCaptureDirCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fPlacementFileCmd = new G4UIcmdWithAString("/cfg/setPlacementFile", this);
  fPlacementFileCmd->SetGuidance("Set explicit placement CSV path for geometry construction.");
  fPlacementFileCmd->SetParameterName("placementFilePath", false);
  fPlacementFileCmd->AvailableForStates(G4State_PreInit);

  fPlacementGeometryModeCmd = new G4UIcmdWithAString("/cfg/setPlacementGeometryMode", this);
  fPlacementGeometryModeCmd->SetGuidance("Set placement geometry mode: pbc_clipped | primary_only.");
  fPlacementGeometryModeCmd->SetParameterName("placementGeometryMode", false);
  fPlacementGeometryModeCmd->AvailableForStates(G4State_PreInit);

  fPeriodicImagesFileCmd = new G4UIcmdWithAString("/cfg/setPeriodicImagesCsv", this);
  fPeriodicImagesFileCmd->SetGuidance("Override the periodic images CSV for a format_version=3 placement.");
  fPeriodicImagesFileCmd->SetParameterName("periodicImagesFilePath", false);
  fPeriodicImagesFileCmd->AvailableForStates(G4State_PreInit);

  fUseRandomPlacementCmd = new G4UIcmdWithABool("/cfg/setUseRandomPlacement", this);
  fUseRandomPlacementCmd->SetGuidance("Enable random placement selection from the current BN:ZnS ratio folder.");
  fUseRandomPlacementCmd->SetParameterName("useRandomPlacement", false);
  fUseRandomPlacementCmd->AvailableForStates(G4State_PreInit);

  fAllowThicknessEqualCmd = new G4UIcmdWithABool("/cfg/setAllowThicknessEqualLocalPatch", this);
  fAllowThicknessEqualCmd->SetGuidance("Allow input thickness == local patch thickness in Stage B.");
  fAllowThicknessEqualCmd->SetParameterName("allowEqual", false);
  fAllowThicknessEqualCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDDir = new G4UIdirectory("/cfg/stageD/");
  fStageDDir->SetGuidance("Stage D optical homogenization configuration.");

  fStageDWavelengthNmCmd = new G4UIcmdWithADouble("/cfg/stageD/setWavelengthNm", this);
  fStageDWavelengthNmCmd->SetGuidance("Set Stage D optical photon wavelength in nm.");
  fStageDWavelengthNmCmd->SetParameterName("wavelengthNm", false);
  fStageDWavelengthNmCmd->SetRange("wavelengthNm>0.");
  fStageDWavelengthNmCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDSourceModeCmd = new G4UIcmdWithAString("/cfg/stageD/setSourceMode", this);
  fStageDSourceModeCmd->SetGuidance("Set Stage D source mode: uniform_ZnS | uniform_all_phase.");
  fStageDSourceModeCmd->SetParameterName("sourceMode", false);
  fStageDSourceModeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDBoundaryModeCmd = new G4UIcmdWithAString("/cfg/stageD/setBoundaryMode", this);
  fStageDBoundaryModeCmd->SetGuidance("Set Stage D boundary mode: periodic_wrap | escape | same_phase_reentry.");
  fStageDBoundaryModeCmd->SetParameterName("boundaryMode", false);
  fStageDBoundaryModeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDReentryModeCmd = new G4UIcmdWithAString("/cfg/stageD/setReentryMode", this);
  fStageDReentryModeCmd->SetGuidance("Set Stage D re-entry family mode: state_matched | same_phase_rho_over_R | same_phase_random.");
  fStageDReentryModeCmd->SetParameterName("reentryMode", false);
  fStageDReentryModeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDParticleReentryModeCmd = new G4UIcmdWithAString("/cfg/stageD/setParticleReentryMode", this);
  fStageDParticleReentryModeCmd->SetGuidance("Set Stage D particle re-entry mode: sphere_q_mu | same_phase_rho_over_R | same_phase_random.");
  fStageDParticleReentryModeCmd->SetParameterName("particleReentryMode", false);
  fStageDParticleReentryModeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDMatrixReentryModeCmd = new G4UIcmdWithAString("/cfg/stageD/setMatrixReentryMode", this);
  fStageDMatrixReentryModeCmd->SetGuidance("Set Stage D matrix re-entry mode: clearance_binned_portal | random_matrix_debug | distance_matched_matrix.");
  fStageDMatrixReentryModeCmd->SetParameterName("matrixReentryMode", false);
  fStageDMatrixReentryModeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDScatterMetricCmd = new G4UIcmdWithAString("/cfg/stageD/setScatterMetric", this);
  fStageDScatterMetricCmd->SetGuidance(
      "Set Stage D primary scatter metric: particle_encounter_angle_threshold | angle_threshold | "
      "particle_encounter_no_threshold.");
  fStageDScatterMetricCmd->SetParameterName("scatterMetric", false);
  fStageDScatterMetricCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDTargetPrimaryScatterCmd = new G4UIcmdWithAnInteger("/cfg/stageD/setTargetPrimaryScatter", this);
  fStageDTargetPrimaryScatterCmd->SetGuidance(
      "Set Stage D target number of primary scatter events per photon before early stop. Use 0 to disable.");
  fStageDTargetPrimaryScatterCmd->SetParameterName("targetPrimaryScatter", false);
  fStageDTargetPrimaryScatterCmd->SetRange("targetPrimaryScatter>=0");
  fStageDTargetPrimaryScatterCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDThetaThresholdDegCmd = new G4UIcmdWithADouble("/cfg/stageD/setThetaThresholdDeg", this);
  fStageDThetaThresholdDegCmd->SetGuidance("Set Stage D minimum direction change angle counted as effective scatter.");
  fStageDThetaThresholdDegCmd->SetParameterName("thetaThresholdDeg", false);
  fStageDThetaThresholdDegCmd->SetRange("thetaThresholdDeg>=0.");
  fStageDThetaThresholdDegCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDMaxReentryCmd = new G4UIcmdWithAnInteger("/cfg/stageD/setMaxReentry", this);
  fStageDMaxReentryCmd->SetGuidance("Set Stage D maximum allowed statistical re-entries per photon.");
  fStageDMaxReentryCmd->SetParameterName("maxReentry", false);
  fStageDMaxReentryCmd->SetRange("maxReentry>0");
  fStageDMaxReentryCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDMaxStepsCmd = new G4UIcmdWithAnInteger("/cfg/stageD/setMaxSteps", this);
  fStageDMaxStepsCmd->SetGuidance("Set Stage D maximum allowed Geant4 steps per photon.");
  fStageDMaxStepsCmd->SetParameterName("maxSteps", false);
  fStageDMaxStepsCmd->SetRange("maxSteps>0");
  fStageDMaxStepsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDMaxPathLengthUmCmd = new G4UIcmdWithADouble("/cfg/stageD/setMaxPathLengthUm", this);
  fStageDMaxPathLengthUmCmd->SetGuidance("Set Stage D maximum accumulated physical path length per photon in um.");
  fStageDMaxPathLengthUmCmd->SetParameterName("maxPathLengthUm", false);
  fStageDMaxPathLengthUmCmd->SetRange("maxPathLengthUm>0.");
  fStageDMaxPathLengthUmCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDOutputDirCmd = new G4UIcmdWithAString("/cfg/stageD/setOutputDir", this);
  fStageDOutputDirCmd->SetGuidance("Override Stage D output directory.");
  fStageDOutputDirCmd->SetParameterName("outputDir", false);
  fStageDOutputDirCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDPortalNuCmd = new G4UIcmdWithAnInteger("/cfg/stageD/setPortalNu", this);
  fStageDPortalNuCmd->SetGuidance("Set Stage D virtual portal sampling count along face U.");
  fStageDPortalNuCmd->SetParameterName("portalNu", false);
  fStageDPortalNuCmd->SetRange("portalNu>0");
  fStageDPortalNuCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDPortalNvCmd = new G4UIcmdWithAnInteger("/cfg/stageD/setPortalNv", this);
  fStageDPortalNvCmd->SetGuidance("Set Stage D virtual portal sampling count along face V.");
  fStageDPortalNvCmd->SetParameterName("portalNv", false);
  fStageDPortalNvCmd->SetRange("portalNv>0");
  fStageDPortalNvCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDPortalMarginUmCmd = new G4UIcmdWithADouble("/cfg/stageD/setPortalMarginUm", this);
  fStageDPortalMarginUmCmd->SetGuidance("Set Stage D portal virtual box margin in um. Use <=0 to auto-pick max particle radius.");
  fStageDPortalMarginUmCmd->SetParameterName("portalMarginUm", false);
  fStageDPortalMarginUmCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDClearanceBinEdgesCmd = new G4UIcommand("/cfg/stageD/setClearanceBinEdgesUm", this);
  fStageDClearanceBinEdgesCmd->SetGuidance("Set Stage D matrix clearance bin edges in um: edge0 edge1 edge2.");
  auto *clearance0 = new G4UIparameter("edge0Um", 'd', false);
  auto *clearance1 = new G4UIparameter("edge1Um", 'd', false);
  auto *clearance2 = new G4UIparameter("edge2Um", 'd', false);
  fStageDClearanceBinEdgesCmd->SetParameter(clearance0);
  fStageDClearanceBinEdgesCmd->SetParameter(clearance1);
  fStageDClearanceBinEdgesCmd->SetParameter(clearance2);
  fStageDClearanceBinEdgesCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDMaxParticleReentryTrialsCmd = new G4UIcmdWithAnInteger("/cfg/stageD/setMaxParticleReentryTrials", this);
  fStageDMaxParticleReentryTrialsCmd->SetGuidance("Set Stage D maximum trials for particle q/mu same-phase re-entry.");
  fStageDMaxParticleReentryTrialsCmd->SetParameterName("maxTrials", false);
  fStageDMaxParticleReentryTrialsCmd->SetRange("maxTrials>0");
  fStageDMaxParticleReentryTrialsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDMaxPortalFallbackLevelCmd = new G4UIcmdWithAnInteger("/cfg/stageD/setMaxPortalFallbackLevel", this);
  fStageDMaxPortalFallbackLevelCmd->SetGuidance("Set Stage D matrix portal fallback depth: 0..4.");
  fStageDMaxPortalFallbackLevelCmd->SetParameterName("maxPortalFallbackLevel", false);
  fStageDMaxPortalFallbackLevelCmd->SetRange("maxPortalFallbackLevel>=0");
  fStageDMaxPortalFallbackLevelCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDWriteReentryDiagnosticsCmd =
      new G4UIcmdWithABool("/cfg/stageD/setWriteReentryDiagnostics", this);
  fStageDWriteReentryDiagnosticsCmd->SetGuidance(
      "Enable Stage D per-reentry diagnostic CSV output.");
  fStageDWriteReentryDiagnosticsCmd->SetParameterName("enabled", false);
  fStageDWriteReentryDiagnosticsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDReentryDiagnosticsSamplingRateCmd =
      new G4UIcmdWithADouble("/cfg/stageD/setReentryDiagnosticsSamplingRate", this);
  fStageDReentryDiagnosticsSamplingRateCmd->SetGuidance(
      "Set deterministic Stage D re-entry diagnostic sampling rate in [0,1].");
  fStageDReentryDiagnosticsSamplingRateCmd->SetParameterName("samplingRate", false);
  fStageDReentryDiagnosticsSamplingRateCmd->SetRange(
      "samplingRate>=0. && samplingRate<=1.");
  fStageDReentryDiagnosticsSamplingRateCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fStageDMaxDiagnosticRowsCmd =
      new G4UIcmdWithAnInteger("/cfg/stageD/setMaxDiagnosticRows", this);
  fStageDMaxDiagnosticRowsCmd->SetGuidance(
      "Set maximum Stage D re-entry diagnostic rows; 0 disables row output.");
  fStageDMaxDiagnosticRowsCmd->SetParameterName("maxRows", false);
  fStageDMaxDiagnosticRowsCmd->SetRange("maxRows>=0");
  fStageDMaxDiagnosticRowsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fOpticalParamsCmd = new G4UIcommand("/cfg/setOpticalParams", this);
  fOpticalParamsCmd->SetGuidance("Set optical material parameters used by Stage D.");
  fOpticalParamsCmd->SetGuidance("Usage: /cfg/setOpticalParams <matrix_n> <matrix_abs_um> <bn_n> <bn_abs_um> <zns_n> <zns_abs_um>");

  auto *matrixNParam = new G4UIparameter("matrixRIndex", 'd', false);
  fOpticalParamsCmd->SetParameter(matrixNParam);
  auto *matrixAbsParam = new G4UIparameter("matrixAbsLengthUm", 'd', false);
  fOpticalParamsCmd->SetParameter(matrixAbsParam);
  auto *bnNParam = new G4UIparameter("bnRIndex", 'd', false);
  fOpticalParamsCmd->SetParameter(bnNParam);
  auto *bnAbsParam = new G4UIparameter("bnAbsLengthUm", 'd', false);
  fOpticalParamsCmd->SetParameter(bnAbsParam);
  auto *znsNParam = new G4UIparameter("znsRIndex", 'd', false);
  fOpticalParamsCmd->SetParameter(znsNParam);
  auto *znsAbsParam = new G4UIparameter("znsAbsLengthUm", 'd', false);
  fOpticalParamsCmd->SetParameter(znsAbsParam);
  fOpticalParamsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fWeightRatioCmd = new G4UIcommand("/cfg/setWeightRatio", this);
  fWeightRatioCmd->SetGuidance("Set BN:ZnS weight ratio used for geometry construction.");
  fWeightRatioCmd->SetGuidance("Usage: /cfg/setWeightRatio <bnWt> <znsWt>");

  auto *bnParam = new G4UIparameter("bnWt", 'd', false);
  bnParam->SetGuidance("BN weight part, e.g. 1");
  fWeightRatioCmd->SetParameter(bnParam);

  auto *znsParam = new G4UIparameter("znsWt", 'd', false);
  znsParam->SetGuidance("ZnS weight part, e.g. 2");
  fWeightRatioCmd->SetParameter(znsParam);

  fWeightRatioCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

AnalysisMessenger::~AnalysisMessenger()
{
  delete fWeightRatioCmd;
  delete fOpticalParamsCmd;
  delete fStageDMaxDiagnosticRowsCmd;
  delete fStageDReentryDiagnosticsSamplingRateCmd;
  delete fStageDWriteReentryDiagnosticsCmd;
  delete fStageDOutputDirCmd;
  delete fStageDMaxPortalFallbackLevelCmd;
  delete fStageDMaxParticleReentryTrialsCmd;
  delete fStageDClearanceBinEdgesCmd;
  delete fStageDPortalMarginUmCmd;
  delete fStageDPortalNvCmd;
  delete fStageDPortalNuCmd;
  delete fStageDMaxPathLengthUmCmd;
  delete fStageDMaxStepsCmd;
  delete fStageDMaxReentryCmd;
  delete fStageDThetaThresholdDegCmd;
  delete fStageDMatrixReentryModeCmd;
  delete fStageDParticleReentryModeCmd;
  delete fStageDReentryModeCmd;
  delete fStageDBoundaryModeCmd;
  delete fStageDSourceModeCmd;
  delete fStageDWavelengthNmCmd;
  delete fStageDDir;
  delete fAllowThicknessEqualCmd;
  delete fUseRandomPlacementCmd;
  delete fPeriodicImagesFileCmd;
  delete fPlacementGeometryModeCmd;
  delete fPlacementFileCmd;
  delete fCaptureDirCmd;
  delete fCaptureCsvCmd;
  delete fRunModeCmd;
  delete fCfgDir;
}

void AnalysisMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
  if (fConfig == nullptr)
  {
    G4Exception("AnalysisMessenger::SetNewValue",
                "BNZS_CFG_002", FatalException,
                "AnalysisConfig pointer is null.");
    return;
  }

  if (command == fRunModeCmd)
  {
    fConfig->runMode = ParseRunModeOrThrow(newValue);

    G4cout << "[AnalysisMessenger] runMode set to "
           << AnalysisConfig::RunModeName(fConfig->runMode)
           << G4endl;
    return;
  }

  if (command == fCaptureCsvCmd)
  {
    fConfig->captureCsvPath = Trim(newValue);

    G4double dirBnWt = 0.0;
    G4double dirZnsWt = 0.0;
    const std::string dirName = std::filesystem::path(fConfig->captureCsvPath).parent_path().filename().string();
    if (TryParseRatioFolderName(dirName, dirBnWt, dirZnsWt))
    {
      fConfig->bnWt = dirBnWt;
      fConfig->znsWt = dirZnsWt;
      G4cout << "[AnalysisMessenger] inferred weight ratio from captureCsvPath: BN:ZnS = "
             << fConfig->bnWt << ":" << fConfig->znsWt
             << G4endl;
    }

    G4cout << "[AnalysisMessenger] captureCsvPath set to "
           << fConfig->captureCsvPath
           << G4endl;
    return;
  }

  if (command == fCaptureDirCmd)
  {
    fConfig->captureInputDir = Trim(newValue);

    G4double dirBnWt = 0.0;
    G4double dirZnsWt = 0.0;
    const std::string dirName = std::filesystem::path(fConfig->captureInputDir).filename().string();
    if (TryParseRatioFolderName(dirName, dirBnWt, dirZnsWt))
    {
      fConfig->bnWt = dirBnWt;
      fConfig->znsWt = dirZnsWt;
      G4cout << "[AnalysisMessenger] inferred weight ratio from captureInputDir: BN:ZnS = "
             << fConfig->bnWt << ":" << fConfig->znsWt
             << G4endl;
    }

    G4cout << "[AnalysisMessenger] captureInputDir set to "
           << fConfig->captureInputDir
           << G4endl;
    return;
  }

  if (command == fPlacementFileCmd)
  {
    const std::string requestedPlacement = Trim(newValue);
    fConfig->placementFilePath = requestedPlacement;
    fConfig->useRandomPlacement = false;

    G4cout << "[AnalysisMessenger] placementFilePath set to "
           << fConfig->placementFilePath
           << " ; useRandomPlacement=false"
           << G4endl;
    return;
  }

  if (command == fUseRandomPlacementCmd)
  {
    fConfig->useRandomPlacement =
        fUseRandomPlacementCmd->GetNewBoolValue(newValue);

    G4cout << "[AnalysisMessenger] useRandomPlacement set to "
           << (fConfig->useRandomPlacement ? "true" : "false")
           << G4endl;
    return;
  }

  if (command == fPlacementGeometryModeCmd)
  {
    const std::string mode = ToLowerCopy(Trim(newValue));
    if (mode != "pbc_clipped" && mode != "primary_only")
    {
      G4Exception("AnalysisMessenger::SetNewValue", "BNZS_CFG_020", FatalException,
                  "placementGeometryMode must be pbc_clipped or primary_only.");
      return;
    }
    fConfig->placementGeometryMode = mode;
    G4cout << "[AnalysisMessenger] placementGeometryMode set to " << mode << G4endl;
    return;
  }

  if (command == fPeriodicImagesFileCmd)
  {
    fConfig->periodicImagesFilePath = Trim(newValue);
    G4cout << "[AnalysisMessenger] periodicImagesFilePath set to "
           << fConfig->periodicImagesFilePath << G4endl;
    return;
  }

  if (command == fAllowThicknessEqualCmd)
  {
    fConfig->allowThicknessEqualLocalPatch =
        fAllowThicknessEqualCmd->GetNewBoolValue(newValue);

    G4cout << "[AnalysisMessenger] allowThicknessEqualLocalPatch set to "
           << (fConfig->allowThicknessEqualLocalPatch ? "true" : "false")
           << G4endl;
    return;
  }
  if (command == fStageDWavelengthNmCmd)
  {
    fConfig->stageD_wavelength_nm =
        fStageDWavelengthNmCmd->GetNewDoubleValue(newValue);
    G4cout << "[AnalysisMessenger] stageD_wavelength_nm set to "
           << fConfig->stageD_wavelength_nm
           << G4endl;
    return;
  }
  if (command == fStageDSourceModeCmd)
  {
    const std::string normalized = NormalizeStageDSourceMode(newValue);
    if (normalized.empty())
    {
      G4Exception("AnalysisMessenger::SetNewValue",
                  "BNZS_CFG_009", FatalException,
                  "Stage D sourceMode must be uniform_ZnS or uniform_all_phase.");
      return;
    }
    fConfig->stageD_source_mode = normalized;
    G4cout << "[AnalysisMessenger] stageD_source_mode set to "
           << fConfig->stageD_source_mode
           << G4endl;
    return;
  }
  if (command == fStageDBoundaryModeCmd)
  {
    const std::string normalized = NormalizeStageDBoundaryMode(newValue);
    if (normalized.empty())
    {
      G4Exception("AnalysisMessenger::SetNewValue",
                  "BNZS_CFG_010", FatalException,
                  "Stage D boundaryMode must be periodic_wrap, escape, or same_phase_reentry.");
      return;
    }
    fConfig->stageD_boundary_mode = normalized;
    G4cout << "[AnalysisMessenger] stageD_boundary_mode set to "
           << fConfig->stageD_boundary_mode
           << G4endl;
    return;
  }
  if (command == fStageDReentryModeCmd)
  {
    const std::string normalized = NormalizeStageDReentryMode(newValue);
    if (normalized.empty())
    {
      G4Exception("AnalysisMessenger::SetNewValue",
                  "BNZS_CFG_011", FatalException,
                  "Stage D reentryMode must be state_matched, same_phase_rho_over_R, or same_phase_random.");
      return;
    }
    fConfig->stageD_reentry_mode = normalized;
    G4cout << "[AnalysisMessenger] stageD_reentry_mode set to "
           << fConfig->stageD_reentry_mode
           << G4endl;
    return;
  }
  if (command == fStageDParticleReentryModeCmd)
  {
    const std::string normalized = NormalizeStageDParticleReentryMode(newValue);
    if (normalized.empty())
    {
      G4Exception("AnalysisMessenger::SetNewValue",
                  "BNZS_CFG_011A", FatalException,
                  "Stage D particleReentryMode must be sphere_q_mu, same_phase_rho_over_R, or same_phase_random.");
      return;
    }
    fConfig->stageD_particle_reentry_mode = normalized;
    G4cout << "[AnalysisMessenger] stageD_particle_reentry_mode set to "
           << fConfig->stageD_particle_reentry_mode
           << G4endl;
    return;
  }
  if (command == fStageDMatrixReentryModeCmd)
  {
    const std::string normalized = NormalizeStageDMatrixReentryMode(newValue);
    if (normalized.empty())
    {
      G4Exception("AnalysisMessenger::SetNewValue",
                  "BNZS_CFG_012", FatalException,
                  "Stage D matrixReentryMode must be clearance_binned_portal, random_matrix_debug, random_matrix, or distance_matched_matrix.");
      return;
    }
    fConfig->stageD_matrix_reentry_mode = normalized;
    G4cout << "[AnalysisMessenger] stageD_matrix_reentry_mode set to "
           << fConfig->stageD_matrix_reentry_mode
           << G4endl;
    return;
  }
  if (command == fStageDScatterMetricCmd)
  {
    const std::string normalized = NormalizeStageDScatterMetric(newValue);
    if (normalized.empty())
    {
      G4Exception("AnalysisMessenger::SetNewValue",
                  "BNZS_CFG_013", FatalException,
                  "Stage D scatterMetric must be particle_encounter_angle_threshold (or angle_threshold) or particle_encounter_no_threshold.");
      return;
    }
    fConfig->stageD_scatter_metric = normalized;
    G4cout << "[AnalysisMessenger] stageD_scatter_metric set to "
           << fConfig->stageD_scatter_metric
           << G4endl;
    return;
  }
  if (command == fStageDTargetPrimaryScatterCmd)
  {
    fConfig->stageD_target_primary_scatter =
        fStageDTargetPrimaryScatterCmd->GetNewIntValue(newValue);
    G4cout << "[AnalysisMessenger] stageD_target_primary_scatter set to "
           << fConfig->stageD_target_primary_scatter
           << G4endl;
    return;
  }
  if (command == fStageDThetaThresholdDegCmd)
  {
    fConfig->stageD_theta_threshold_deg =
        fStageDThetaThresholdDegCmd->GetNewDoubleValue(newValue);
    G4cout << "[AnalysisMessenger] stageD_theta_threshold_deg set to "
           << fConfig->stageD_theta_threshold_deg
           << G4endl;
    return;
  }
  if (command == fStageDMaxReentryCmd)
  {
    fConfig->stageD_max_reentry =
        fStageDMaxReentryCmd->GetNewIntValue(newValue);
    G4cout << "[AnalysisMessenger] stageD_max_reentry set to "
           << fConfig->stageD_max_reentry
           << G4endl;
    return;
  }
  if (command == fStageDMaxStepsCmd)
  {
    fConfig->stageD_max_steps =
        fStageDMaxStepsCmd->GetNewIntValue(newValue);
    G4cout << "[AnalysisMessenger] stageD_max_steps set to "
           << fConfig->stageD_max_steps
           << G4endl;
    return;
  }
  if (command == fStageDMaxPathLengthUmCmd)
  {
    fConfig->stageD_max_path_length_um =
        fStageDMaxPathLengthUmCmd->GetNewDoubleValue(newValue);
    G4cout << "[AnalysisMessenger] stageD_max_path_length_um set to "
           << fConfig->stageD_max_path_length_um
           << G4endl;
    return;
  }
  if (command == fStageDOutputDirCmd)
  {
    fConfig->stageD_output_dir = Trim(newValue);
    G4cout << "[AnalysisMessenger] stageD_output_dir set to "
           << fConfig->stageD_output_dir
           << G4endl;
    return;
  }
  if (command == fStageDPortalNuCmd)
  {
    fConfig->stageD_portal_nu = fStageDPortalNuCmd->GetNewIntValue(newValue);
    G4cout << "[AnalysisMessenger] stageD_portal_nu set to "
           << fConfig->stageD_portal_nu
           << G4endl;
    return;
  }
  if (command == fStageDPortalNvCmd)
  {
    fConfig->stageD_portal_nv = fStageDPortalNvCmd->GetNewIntValue(newValue);
    G4cout << "[AnalysisMessenger] stageD_portal_nv set to "
           << fConfig->stageD_portal_nv
           << G4endl;
    return;
  }
  if (command == fStageDPortalMarginUmCmd)
  {
    fConfig->stageD_portal_margin_um =
        fStageDPortalMarginUmCmd->GetNewDoubleValue(newValue);
    G4cout << "[AnalysisMessenger] stageD_portal_margin_um set to "
           << fConfig->stageD_portal_margin_um
           << G4endl;
    return;
  }
  if (command == fStageDClearanceBinEdgesCmd)
  {
    std::istringstream iss(Trim(newValue));
    G4double edge0 = 0.0;
    G4double edge1 = 0.0;
    G4double edge2 = 0.0;
    if (!(iss >> edge0 >> edge1 >> edge2) ||
        !(edge0 > 0.0 && edge0 < edge1 && edge1 < edge2))
    {
      G4Exception("AnalysisMessenger::SetNewValue",
                  "BNZS_CFG_012A", FatalException,
                  "Stage D clearance bin edges must satisfy 0 < edge0 < edge1 < edge2.");
      return;
    }
    fConfig->stageD_clearance_bin0_um = edge0;
    fConfig->stageD_clearance_bin1_um = edge1;
    fConfig->stageD_clearance_bin2_um = edge2;
    G4cout << "[AnalysisMessenger] stageD clearance bin edges set to "
           << edge0 << ", " << edge1 << ", " << edge2
           << " um" << G4endl;
    return;
  }
  if (command == fStageDMaxParticleReentryTrialsCmd)
  {
    fConfig->stageD_max_particle_reentry_trials =
        fStageDMaxParticleReentryTrialsCmd->GetNewIntValue(newValue);
    G4cout << "[AnalysisMessenger] stageD_max_particle_reentry_trials set to "
           << fConfig->stageD_max_particle_reentry_trials
           << G4endl;
    return;
  }
  if (command == fStageDMaxPortalFallbackLevelCmd)
  {
    fConfig->stageD_max_portal_fallback_level =
        fStageDMaxPortalFallbackLevelCmd->GetNewIntValue(newValue);
    G4cout << "[AnalysisMessenger] stageD_max_portal_fallback_level set to "
           << fConfig->stageD_max_portal_fallback_level
           << G4endl;
    return;
  }
  if (command == fStageDWriteReentryDiagnosticsCmd)
  {
    fConfig->stageD_write_reentry_diagnostics =
        fStageDWriteReentryDiagnosticsCmd->GetNewBoolValue(newValue);
    G4cout << "[AnalysisMessenger] stageD_write_reentry_diagnostics set to "
           << (fConfig->stageD_write_reentry_diagnostics ? "true" : "false")
           << G4endl;
    return;
  }
  if (command == fStageDReentryDiagnosticsSamplingRateCmd)
  {
    fConfig->stageD_reentry_diagnostics_sampling_rate =
        fStageDReentryDiagnosticsSamplingRateCmd->GetNewDoubleValue(newValue);
    G4cout << "[AnalysisMessenger] stageD_reentry_diagnostics_sampling_rate set to "
           << fConfig->stageD_reentry_diagnostics_sampling_rate
           << G4endl;
    return;
  }
  if (command == fStageDMaxDiagnosticRowsCmd)
  {
    fConfig->stageD_max_diagnostic_rows =
        fStageDMaxDiagnosticRowsCmd->GetNewIntValue(newValue);
    G4cout << "[AnalysisMessenger] stageD_max_diagnostic_rows set to "
           << fConfig->stageD_max_diagnostic_rows
           << G4endl;
    return;
  }
  if (command == fOpticalParamsCmd)
  {
    std::istringstream iss(Trim(newValue));
    G4double matrixN = 0.0;
    G4double matrixAbsUm = 0.0;
    G4double bnN = 0.0;
    G4double bnAbsUm = 0.0;
    G4double znsN = 0.0;
    G4double znsAbsUm = 0.0;

    if (!(iss >> matrixN >> matrixAbsUm >> bnN >> bnAbsUm >> znsN >> znsAbsUm) ||
        matrixN <= 0.0 || matrixAbsUm <= 0.0 ||
        bnN <= 0.0 || bnAbsUm <= 0.0 ||
        znsN <= 0.0 || znsAbsUm <= 0.0)
    {
      G4Exception("AnalysisMessenger::SetNewValue",
                  "BNZS_CFG_008", FatalException,
                  "Optical params must be six positive numbers: matrix_n matrix_abs_um bn_n bn_abs_um zns_n zns_abs_um.");
      return;
    }

    fConfig->opticalMatrixRIndex = matrixN;
    fConfig->opticalMatrixAbsLengthUm = matrixAbsUm;
    fConfig->opticalBnRIndex = bnN;
    fConfig->opticalBnAbsLengthUm = bnAbsUm;
    fConfig->opticalZnsRIndex = znsN;
    fConfig->opticalZnsAbsLengthUm = znsAbsUm;
    fConfig->opticalParamsProvided = true;

    G4cout << "[AnalysisMessenger] optical params set:"
           << " matrix n/abs_um=" << matrixN << "/" << matrixAbsUm
           << " BN n/abs_um=" << bnN << "/" << bnAbsUm
           << " ZnS n/abs_um=" << znsN << "/" << znsAbsUm
           << G4endl;
    return;
  }
  if (command == fWeightRatioCmd)
  {
    std::istringstream iss(Trim(newValue));
    G4double bnWt = 0.0;
    G4double znsWt = 0.0;

    if (!(iss >> bnWt >> znsWt))
    {
      G4Exception("AnalysisMessenger::SetNewValue",
                  "BNZS_CFG_003", FatalException,
                  ("Failed to parse /cfg/setWeightRatio arguments: " + std::string(newValue)).c_str());
      return;
    }

    if (bnWt <= 0.0 || znsWt <= 0.0)
    {
      G4Exception("AnalysisMessenger::SetNewValue",
                  "BNZS_CFG_004", FatalException,
                  "BN and ZnS weight parts must both be > 0.");
      return;
    }

    fConfig->bnWt = bnWt;
    fConfig->znsWt = znsWt;

    G4cout << "[AnalysisMessenger] weight ratio set to BN:ZnS = "
           << fConfig->bnWt << ":" << fConfig->znsWt
           << G4endl;
    return;
  }
}
