#include "ModeSteppingAction.hh"

#include "AnalysisConfig.hh"
#include "ModeEventAction.hh"
#include "ModeRunAction.hh"
#include "ModePrimaryGeneratorAction.hh"

#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "SteppingAction.hh"
#include "StageARunAction.hh"
#include "StageASteppingAction.hh"
#include "StageDOpticalEventAction.hh"
#include "StageDOpticalRunAction.hh"
#include "StageDOpticalSteppingAction.hh"

#include "G4Step.hh"
#include "G4EventManager.hh"
#include "G4Exception.hh"
#include "G4ios.hh"

ModeSteppingAction::ModeSteppingAction(ModeRunAction *modeRunAction,
                                       AnalysisConfig *config,
                                       ModePrimaryGeneratorAction *modePrimaryAction,
                                       EventAction *stageBEventAction,
                                       PrimaryGeneratorAction *stageBPrimaryAction,
                                       StageDOpticalEventAction *stageDEventAction)
    : G4UserSteppingAction(),
      fConfig(config),
      fModePrimaryAction(modePrimaryAction),
      fStageDRunAction(nullptr),
      fStageBSteppingAction(nullptr),
      fStageASteppingAction(nullptr),
      fStageDSteppingAction(nullptr)
{
    if (fConfig == nullptr)
    {
        G4Exception("ModeSteppingAction::ModeSteppingAction",
                    "BNZS_MODE_STEP_001", FatalException,
                    "AnalysisConfig pointer is null.");
        return;
    }

    if (modeRunAction == nullptr)
    {
        G4Exception("ModeSteppingAction::ModeSteppingAction",
                    "BNZS_MODE_STEP_002", FatalException,
                    "ModeRunAction pointer is null.");
        return;
    }

    StageARunAction *stageARunAction = modeRunAction->GetStageARunAction();
    if (stageARunAction == nullptr)
    {
        G4Exception("ModeSteppingAction::ModeSteppingAction",
                    "BNZS_MODE_STEP_003", FatalException,
                    "Stage A run action from ModeRunAction is null.");
        return;
    }

    // Stage A stepping implementation
    fStageASteppingAction = new StageASteppingAction(stageARunAction, fConfig);

    // Stage B stepping implementation
    if (stageBEventAction != nullptr && stageBPrimaryAction != nullptr)
    {
        fStageBSteppingAction = new SteppingAction(stageBEventAction, stageBPrimaryAction);
    }

    fStageDRunAction = modeRunAction->GetStageDRunAction();
    if (fStageDRunAction != nullptr && stageDEventAction != nullptr)
    {
        fStageDSteppingAction = new StageDOpticalSteppingAction(
            fStageDRunAction,
            stageDEventAction,
            fConfig);
    }

    G4cout << "[ModeSteppingAction] Dispatcher initialized."
           << " current runMode = "
           << AnalysisConfig::RunModeName(fConfig->runMode)
           << G4endl;
}

ModeSteppingAction::~ModeSteppingAction()
{
    delete fStageASteppingAction;
    delete fStageBSteppingAction;
    delete fStageDSteppingAction;
}

StageDOpticalSteppingAction *ModeSteppingAction::EnsureStageDSteppingAction()
{
    if (fStageDSteppingAction != nullptr)
        return fStageDSteppingAction;

    if (fStageDRunAction == nullptr && fModePrimaryAction != nullptr)
    {
        // fStageDRunAction is owned by ModeRunAction and should already exist
        // whenever dispatcher actions were built; keep this branch defensive.
    }

    auto *eventDispatcher = dynamic_cast<ModeEventAction *>(
        G4EventManager::GetEventManager()->GetUserEventAction());
    StageDOpticalEventAction *stageDEventAction =
        eventDispatcher ? eventDispatcher->GetStageDEventAction() : nullptr;

    if (stageDEventAction == nullptr && eventDispatcher != nullptr)
        stageDEventAction = eventDispatcher->GetStageDEventAction();

    if (fStageDRunAction != nullptr && stageDEventAction != nullptr)
    {
        fStageDSteppingAction = new StageDOpticalSteppingAction(
            fStageDRunAction,
            stageDEventAction,
            fConfig);
    }

    return fStageDSteppingAction;
}

void ModeSteppingAction::UserSteppingAction(const G4Step *step)
{
    if (fConfig == nullptr)
    {
        G4Exception("ModeSteppingAction::UserSteppingAction",
                    "BNZS_MODE_STEP_004", FatalException,
                    "AnalysisConfig pointer is null.");
        return;
    }

    switch (fConfig->runMode)
    {
    case RunMode::StageA_NeutronPatch:
        if (fStageASteppingAction == nullptr)
        {
            G4Exception("ModeSteppingAction::UserSteppingAction",
                        "BNZS_MODE_STEP_005", FatalException,
                        "Stage A stepping action is null.");
            return;
        }
        fStageASteppingAction->UserSteppingAction(step);
        return;

    case RunMode::StageB_ReplayAlphaLi:
        if (fStageBSteppingAction == nullptr)
        {
            G4Exception("ModeSteppingAction::UserSteppingAction",
                        "BNZS_MODE_STEP_006", FatalException,
                        "Stage B stepping action is null.");
            return;
        }
        fStageBSteppingAction->UserSteppingAction(step);
        return;

    case RunMode::StageD_OpticalHomogenization:
        if (EnsureStageDSteppingAction() == nullptr)
        {
            G4Exception("ModeSteppingAction::UserSteppingAction",
                        "BNZS_MODE_STEP_010", FatalException,
                        "Stage D optical stepping action is null.");
            return;
        }
        fStageDSteppingAction->UserSteppingAction(step);
        return;

    default:
        G4Exception("ModeSteppingAction::UserSteppingAction",
                    "BNZS_MODE_STEP_008", FatalException,
                    "Unknown run mode.");
        return;
    }
}
