#include "ModeEventAction.hh"

#include "AnalysisConfig.hh"
#include "ModeRunAction.hh"
#include "ModePrimaryGeneratorAction.hh"

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "StageDOpticalEventAction.hh"
#include "StageDOpticalPrimaryGeneratorAction.hh"
#include "StageDOpticalRunAction.hh"

#include "G4Event.hh"
#include "G4Exception.hh"
#include "G4ios.hh"

ModeEventAction::ModeEventAction(ModeRunAction *modeRunAction,
                                 ModePrimaryGeneratorAction *modePrimaryAction,
                                 AnalysisConfig *config)
    : G4UserEventAction(),
      fConfig(config),
      fModeRunAction(modeRunAction),
      fModePrimaryAction(modePrimaryAction),
      fStageBEventAction(nullptr),
      fStageDEventAction(nullptr)
{
  if (fConfig == nullptr)
  {
    G4Exception("ModeEventAction::ModeEventAction",
                "BNZS_MODE_EVT_001", FatalException,
                "AnalysisConfig pointer is null.");
    return;
  }

  if (fModeRunAction == nullptr)
  {
    G4Exception("ModeEventAction::ModeEventAction",
                "BNZS_MODE_EVT_002", FatalException,
                "ModeRunAction pointer is null.");
    return;
  }

  if (fModePrimaryAction == nullptr)
  {
    G4Exception("ModeEventAction::ModeEventAction",
                "BNZS_MODE_EVT_003", FatalException,
                "ModePrimaryGeneratorAction pointer is null.");
    return;
  }

  // Reuse the existing Stage B event implementation
  RunAction *stageBRunAction = fModeRunAction->GetStageBRunAction();
  PrimaryGeneratorAction *stageBPrimaryAction = fModePrimaryAction->GetStageBPrimaryAction();

  if (stageBRunAction != nullptr && stageBPrimaryAction != nullptr)
  {
    fStageBEventAction = new EventAction(stageBRunAction, stageBPrimaryAction);
  }

  StageDOpticalRunAction *stageDRunAction = fModeRunAction->GetStageDRunAction();
  StageDOpticalPrimaryGeneratorAction *stageDPrimaryAction =
      fModePrimaryAction->GetStageDPrimaryAction();
  if (stageDRunAction != nullptr && stageDPrimaryAction != nullptr)
  {
    fStageDEventAction = new StageDOpticalEventAction(
        stageDRunAction,
        stageDPrimaryAction,
        fConfig);
  }

  G4cout << "[ModeEventAction] Dispatcher initialized."
         << " current runMode = "
         << AnalysisConfig::RunModeName(fConfig->runMode)
         << G4endl;
}

ModeEventAction::~ModeEventAction()
{
  delete fStageBEventAction;
  delete fStageDEventAction;
}

EventAction *ModeEventAction::EnsureStageBEventAction()
{
  if (fStageBEventAction != nullptr)
    return fStageBEventAction;
  if (fModeRunAction == nullptr || fModePrimaryAction == nullptr)
    return nullptr;

  RunAction *stageBRunAction = fModeRunAction->GetStageBRunAction();
  PrimaryGeneratorAction *stageBPrimaryAction = fModePrimaryAction->GetStageBPrimaryAction();
  if (stageBRunAction != nullptr && stageBPrimaryAction != nullptr)
  {
    fStageBEventAction = new EventAction(stageBRunAction, stageBPrimaryAction);
  }
  return fStageBEventAction;
}

StageDOpticalEventAction *ModeEventAction::EnsureStageDEventAction()
{
  if (fStageDEventAction != nullptr)
    return fStageDEventAction;
  if (fModeRunAction == nullptr || fModePrimaryAction == nullptr)
    return nullptr;

  StageDOpticalRunAction *stageDRunAction = fModeRunAction->GetStageDRunAction();
  StageDOpticalPrimaryGeneratorAction *stageDPrimaryAction =
      fModePrimaryAction->GetStageDPrimaryAction();
  if (stageDRunAction != nullptr && stageDPrimaryAction != nullptr)
  {
    fStageDEventAction = new StageDOpticalEventAction(
        stageDRunAction,
        stageDPrimaryAction,
        fConfig);
  }
  return fStageDEventAction;
}

void ModeEventAction::BeginOfEventAction(const G4Event *event)
{
  if (fConfig == nullptr)
  {
    G4Exception("ModeEventAction::BeginOfEventAction",
                "BNZS_MODE_EVT_004", FatalException,
                "AnalysisConfig pointer is null.");
    return;
  }

  switch (fConfig->runMode)
  {
  case RunMode::StageA_NeutronPatch:
    // Stage A currently does not need a dedicated event-layer implementation.
    return;

  case RunMode::StageB_ReplayAlphaLi:
    if (EnsureStageBEventAction() == nullptr)
    {
      G4Exception("ModeEventAction::BeginOfEventAction",
                  "BNZS_MODE_EVT_005", FatalException,
                  "Stage B event action is null.");
      return;
    }
    fStageBEventAction->BeginOfEventAction(event);
    return;

  case RunMode::StageD_OpticalHomogenization:
    if (EnsureStageDEventAction() == nullptr)
    {
      G4Exception("ModeEventAction::BeginOfEventAction",
                  "BNZS_MODE_EVT_012", FatalException,
                  "Stage D event action is null.");
      return;
    }
    fStageDEventAction->BeginOfEventAction(event);
    return;

  default:
    G4Exception("ModeEventAction::BeginOfEventAction",
                "BNZS_MODE_EVT_007", FatalException,
                "Unknown run mode.");
    return;
  }
}

void ModeEventAction::EndOfEventAction(const G4Event *event)
{
  if (fConfig == nullptr)
  {
    G4Exception("ModeEventAction::EndOfEventAction",
                "BNZS_MODE_EVT_008", FatalException,
                "AnalysisConfig pointer is null.");
    return;
  }

  switch (fConfig->runMode)
  {
  case RunMode::StageA_NeutronPatch:
    // Stage A currently does not need a dedicated event-layer implementation.
    return;

  case RunMode::StageB_ReplayAlphaLi:
    if (EnsureStageBEventAction() == nullptr)
    {
      G4Exception("ModeEventAction::EndOfEventAction",
                  "BNZS_MODE_EVT_009", FatalException,
                  "Stage B event action is null.");
      return;
    }
    fStageBEventAction->EndOfEventAction(event);
    return;

  case RunMode::StageD_OpticalHomogenization:
    if (EnsureStageDEventAction() == nullptr)
    {
      G4Exception("ModeEventAction::EndOfEventAction",
                  "BNZS_MODE_EVT_013", FatalException,
                  "Stage D event action is null.");
      return;
    }
    fStageDEventAction->EndOfEventAction(event);
    return;

  default:
    G4Exception("ModeEventAction::EndOfEventAction",
                "BNZS_MODE_EVT_011", FatalException,
                "Unknown run mode.");
    return;
  }
}

EventAction *ModeEventAction::GetStageBEventAction() const
{
  return fStageBEventAction;
}

StageDOpticalEventAction *ModeEventAction::GetStageDEventAction() const
{
  return fStageDEventAction;
}
