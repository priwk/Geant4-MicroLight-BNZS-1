#ifndef AnalysisMessenger_h
#define AnalysisMessenger_h 1

#include "G4UImessenger.hh"

class AnalysisConfig;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;
class G4UIcommand;

class AnalysisMessenger : public G4UImessenger
{
public:
  explicit AnalysisMessenger(AnalysisConfig *config);
  ~AnalysisMessenger() override;

  void SetNewValue(G4UIcommand *command, G4String newValue) override;

private:
  AnalysisConfig *fConfig;

  G4UIdirectory *fCfgDir;
  G4UIcmdWithAString *fRunModeCmd;
  G4UIcmdWithAString *fCaptureCsvCmd;
  G4UIcmdWithAString *fCaptureDirCmd;
  G4UIcmdWithAString *fPlacementFileCmd;
  G4UIcmdWithAString *fPlacementGeometryModeCmd;
  G4UIcmdWithAString *fPeriodicImagesFileCmd;
  G4UIcmdWithABool *fUseRandomPlacementCmd;
  G4UIcmdWithABool *fAllowThicknessEqualCmd;
  G4UIdirectory *fStageDDir;
  G4UIcmdWithADouble *fStageDWavelengthNmCmd;
  G4UIcmdWithAString *fStageDSourceModeCmd;
  G4UIcmdWithAString *fStageDBoundaryModeCmd;
  G4UIcmdWithAString *fStageDReentryModeCmd;
  G4UIcmdWithAString *fStageDParticleReentryModeCmd;
  G4UIcmdWithAString *fStageDMatrixReentryModeCmd;
  G4UIcmdWithAString *fStageDScatterMetricCmd;
  G4UIcmdWithAnInteger *fStageDTargetPrimaryScatterCmd;
  G4UIcmdWithADouble *fStageDThetaThresholdDegCmd;
  G4UIcmdWithAnInteger *fStageDMaxReentryCmd;
  G4UIcmdWithAnInteger *fStageDMaxStepsCmd;
  G4UIcmdWithADouble *fStageDMaxPathLengthUmCmd;
  G4UIcmdWithAString *fStageDOutputDirCmd;
  G4UIcmdWithAnInteger *fStageDPortalNuCmd;
  G4UIcmdWithAnInteger *fStageDPortalNvCmd;
  G4UIcmdWithADouble *fStageDPortalMarginUmCmd;
  G4UIcommand *fStageDClearanceBinEdgesCmd;
  G4UIcmdWithAnInteger *fStageDMaxParticleReentryTrialsCmd;
  G4UIcmdWithAnInteger *fStageDMaxPortalFallbackLevelCmd;
  G4UIcommand *fOpticalParamsCmd;
  G4UIcommand *fWeightRatioCmd;
};

#endif
