#include "EventAction.hh"

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "DetectorConstruction.hh"

// --------------------------------------------------------------------

EventAction::EventAction(RunAction *runAction,
                         const PrimaryGeneratorAction *primaryAction)
    : G4UserEventAction(),
      fRunAction(runAction),
      fPrimaryAction(primaryAction),
      fEdep(0.0),
      fCurrentRecord(),
      fCurrentLocalCapturePosition(),
      fCurrentSelectedBNCenter(),
      fCurrentSurfaceMode(""),
      fCurrentTargetLocalZ(0.0),
      fCurrentUsedLocalZ(0.0),
      fCurrentAlphaLiReplayIndex(0),
      fCurrentAlphaLiReplayCount(1)
{
}

// --------------------------------------------------------------------

EventAction::~EventAction() = default;

// --------------------------------------------------------------------

void EventAction::BeginOfEventAction(const G4Event *event)
{
  (void)event;

  fEdep = 0.0;

  if (fPrimaryAction && !fPrimaryAction->HasValidCurrentReplay())
  {
    return;
  }

  if (fPrimaryAction)
  {
    fCurrentRecord = fPrimaryAction->GetCurrentRecord();
    fCurrentLocalCapturePosition = fPrimaryAction->GetCurrentLocalCapturePosition();
    fCurrentSelectedBNCenter = fPrimaryAction->GetCurrentSelectedBNCenter();
    fCurrentSurfaceMode = fPrimaryAction->GetCurrentSurfaceMode();
    fCurrentTargetLocalZ = fPrimaryAction->GetCurrentTargetLocalZ();
    fCurrentUsedLocalZ = fPrimaryAction->GetCurrentUsedLocalZ();
    fCurrentAlphaLiReplayIndex = fPrimaryAction->GetCurrentAlphaLiReplayIndex();
    fCurrentAlphaLiReplayCount = fPrimaryAction->GetCurrentAlphaLiReplayCount();
  }
  else
  {
    fCurrentRecord = CaptureRecord{};
    fCurrentLocalCapturePosition = G4ThreeVector();
    fCurrentSelectedBNCenter = G4ThreeVector();
    fCurrentSurfaceMode.clear();
    fCurrentTargetLocalZ = 0.0;
    fCurrentUsedLocalZ = 0.0;
    fCurrentAlphaLiReplayIndex = 0;
    fCurrentAlphaLiReplayCount = 1;
  }

  if (fRunAction && fPrimaryAction)
  {
    fRunAction->SwitchOutputCsvForInputPath(fPrimaryAction->GetCurrentRecordInputFile());
    fRunAction->AppendCaptureAnchor(MakeCurrentCaptureAnchorRow());
  }
}

// --------------------------------------------------------------------

RunAction::CaptureAnchorRow EventAction::MakeCurrentCaptureAnchorRow() const
{
  RunAction::CaptureAnchorRow row;
  row.physical_event_uid = fPrimaryAction
                               ? fPrimaryAction->MakeCurrentPhysicalEventUid()
                               : "";
  row.eventID = fCurrentRecord.eventID;
  row.record_index = fCurrentRecord.record_index;
  row.thickness_um = fCurrentRecord.thickness_um;
  row.bn_wt = fCurrentRecord.bn_wt;
  row.zns_wt = fCurrentRecord.zns_wt;
  row.capture_x_um = fCurrentRecord.capture_x_um;
  row.capture_y_um = fCurrentRecord.capture_y_um;
  row.source_x_um = fCurrentRecord.source_x_um;
  row.source_y_um = fCurrentRecord.source_y_um;
  row.depth_um = fCurrentRecord.depth_um;
  row.local_capture_x_um = fCurrentLocalCapturePosition.x() / um;
  row.local_capture_y_um = fCurrentLocalCapturePosition.y() / um;
  row.local_capture_z_um = fCurrentLocalCapturePosition.z() / um;
  row.surface_mode = fCurrentSurfaceMode;
  row.target_local_z_um = fCurrentTargetLocalZ / um;
  row.used_local_z_um = fCurrentUsedLocalZ / um;
  row.bn_center_x_um = fCurrentSelectedBNCenter.x() / um;
  row.bn_center_y_um = fCurrentSelectedBNCenter.y() / um;
  row.bn_center_z_um = fCurrentSelectedBNCenter.z() / um;
  row.alphali_replay_index = fCurrentAlphaLiReplayIndex;
  row.alphali_replay_count = fCurrentAlphaLiReplayCount;
  row.trajectory_weight = (fCurrentAlphaLiReplayCount > 0)
                              ? (1.0 / static_cast<G4double>(fCurrentAlphaLiReplayCount))
                              : 1.0;

  std::string placementFile = "unknown";
  const auto *det = dynamic_cast<const DetectorConstruction *>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  if (det)
  {
    placementFile = det->GetLoadedPlacementFileForRecord();
  }
  row.placement_file = placementFile;
  row.source_event_uid = fPrimaryAction
                             ? fPrimaryAction->MakeCurrentSourceEventUid()
                             : "";
  return row;
}

// --------------------------------------------------------------------

void EventAction::EndOfEventAction(const G4Event *event)
{
  if (fPrimaryAction && !fPrimaryAction->HasValidCurrentReplay())
  {
    return;
  }

  // keep quiet for production; useful summary hook left here
  if (event->GetEventID() < 3)
  {
    G4cout
        << "[EventAction] End event " << event->GetEventID()
        << "  input eventID=" << fCurrentRecord.eventID
        << "  record_index=" << fCurrentRecord.record_index
        << "  replay=" << fCurrentAlphaLiReplayIndex << "/" << fCurrentAlphaLiReplayCount
        << "  mode=" << fCurrentSurfaceMode
        << "  total ZnS edep (all tracks)=" << fEdep / keV << " keV"
        << G4endl;
  }
}
