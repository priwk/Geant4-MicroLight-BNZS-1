#include "ActionInitialization.hh"

#include "AnalysisConfig.hh"

#include "ModePrimaryGeneratorAction.hh"
#include "ModeRunAction.hh"
#include "ModeEventAction.hh"
#include "ModeSteppingAction.hh"
#include "StageAStackingAction.hh"

#include "G4Exception.hh"
#include "G4ios.hh"

ActionInitialization::ActionInitialization(AnalysisConfig *config)
    : G4VUserActionInitialization(),
      fConfig(config)
{
}

ActionInitialization::~ActionInitialization() = default;

void ActionInitialization::BuildForMaster() const
{
  if (fConfig == nullptr)
  {
    G4Exception("ActionInitialization::BuildForMaster",
                "BNZS_ACT_000", FatalException,
                "AnalysisConfig pointer is null.");
    return;
  }

  auto *runAction = new ModeRunAction(fConfig);
  SetUserAction(runAction);
}

void ActionInitialization::Build() const
{
  if (fConfig == nullptr)
  {
    G4Exception("ActionInitialization::Build",
                "BNZS_ACT_001", FatalException,
                "AnalysisConfig pointer is null.");
    return;
  }

  G4cout << "[ActionInitialization] Build dispatcher actions."
         << " current runMode = "
         << AnalysisConfig::RunModeName(fConfig->runMode)
         << G4endl;

  auto *primaryAction = new ModePrimaryGeneratorAction(fConfig);
  SetUserAction(primaryAction);

  auto *runAction = new ModeRunAction(fConfig);
  runAction->SetStageBPrimaryAction(primaryAction->GetStageBPrimaryAction());
  SetUserAction(runAction);

  auto *eventAction = new ModeEventAction(runAction, primaryAction, fConfig);
  SetUserAction(eventAction);

  auto *steppingAction = new ModeSteppingAction(
      runAction,
      fConfig,
      primaryAction,
      eventAction->GetStageBEventAction(),
      primaryAction->GetStageBPrimaryAction(),
      eventAction->GetStageDEventAction());
  SetUserAction(steppingAction);

  auto *stackingAction = new StageAStackingAction(fConfig);
  SetUserAction(stackingAction);
}
