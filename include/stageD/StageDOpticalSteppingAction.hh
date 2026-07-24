#ifndef StageDOpticalSteppingAction_h
#define StageDOpticalSteppingAction_h 1

#include "DetectorConstruction.hh"
#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4Step;
class G4Track;
class G4OpBoundaryProcess;
class AnalysisConfig;
class StageDOpticalEventAction;
class StageDOpticalRunAction;
class StageDReentrySampler;

class StageDOpticalSteppingAction : public G4UserSteppingAction
{
public:
  StageDOpticalSteppingAction(StageDOpticalRunAction *runAction,
                              StageDOpticalEventAction *eventAction,
                              AnalysisConfig *config);
  ~StageDOpticalSteppingAction() override;

  void UserSteppingAction(const G4Step *step) override;
  void PrepareForNewRun();

private:
  const DetectorConstruction *ResolveDetector() const;
  G4OpBoundaryProcess *ResolveBoundaryProcess() const;
  G4bool HandleBoundaryReentry(const G4Step *step,
                               G4Track *track,
                               const DetectorConstruction *detector,
                               DetectorConstruction::Phase prePhase);
  G4bool HandleHardPathLimits(G4Track *track);
  G4bool HandleLimitKills(const G4Step *step, G4Track *track);

private:
  AnalysisConfig *fConfig;
  StageDOpticalRunAction *fRunAction;
  StageDOpticalEventAction *fEventAction;
  mutable const DetectorConstruction *fDetector;
  mutable G4OpBoundaryProcess *fBoundaryProcess;
  StageDReentrySampler *fReentrySampler;
};

#endif
