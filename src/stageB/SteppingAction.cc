#include "SteppingAction.hh"

#include "EventAction.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4DynamicParticle.hh"
#include "G4EventManager.hh"
#include "G4StackManager.hh"
#include "G4VUserTrackInformation.hh"

#include <cmath>
#include <fstream>
#include <string>

// --------------------------------------------------------------------
// helpers
namespace
{
    class PeriodicTrackInfo final : public G4VUserTrackInformation
    {
    public:
        PeriodicTrackInfo(G4int originalTrackId, G4int completedSteps)
            : fOriginalTrackId(originalTrackId),
              fCompletedSteps(completedSteps)
        {
        }

        void Print() const override {}
        G4int OriginalTrackId() const { return fOriginalTrackId; }
        G4int CompletedSteps() const { return fCompletedSteps; }

    private:
        G4int fOriginalTrackId;
        G4int fCompletedSteps;
    };

    const PeriodicTrackInfo *GetPeriodicTrackInfo(const G4Track *track)
    {
        return track != nullptr
                   ? dynamic_cast<const PeriodicTrackInfo *>(track->GetUserInformation())
                   : nullptr;
    }

    G4int OutputTrackId(const G4Track *track)
    {
        const auto *info = GetPeriodicTrackInfo(track);
        return info != nullptr ? info->OriginalTrackId() : track->GetTrackID();
    }

    G4int OutputStepId(const G4Track *track)
    {
        const auto *info = GetPeriodicTrackInfo(track);
        return (info != nullptr ? info->CompletedSteps() : 0) +
               track->GetCurrentStepNumber();
    }

    G4bool IsTrackedHeavyParticle(const G4Track *track)
    {
        const auto *def = track->GetDefinition();
        if (!def)
            return false;

        if (def->GetParticleName() == "alpha")
        {
            return true;
        }

        // Li7 ion
        if (def->GetParticleType() == "nucleus" &&
            def->GetAtomicNumber() == 3 &&
            def->GetAtomicMass() == 7)
        {
            return true;
        }

        return false;
    }

    std::string ParticleLabel(const G4Track *track)
    {
        const auto *def = track->GetDefinition();
        if (!def)
            return "unknown";

        if (def->GetParticleName() == "alpha")
        {
            return "alpha";
        }

        if (def->GetParticleType() == "nucleus" &&
            def->GetAtomicNumber() == 3 &&
            def->GetAtomicMass() == 7)
        {
            return "Li7";
        }

        return def->GetParticleName();
    }

    std::string PhaseLabel(const G4VPhysicalVolume *pv, const DetectorConstruction *detector)
    {
        if (!pv)
            return "outside";

        const auto *lv = pv->GetLogicalVolume();
        if (!lv)
            return "outside";

        const auto phase = detector != nullptr
                               ? detector->GetPhaseFromLogicalVolume(lv)
                               : DetectorConstruction::Phase::Unknown;
        if (phase == DetectorConstruction::Phase::BN)
            return "BN";
        if (phase == DetectorConstruction::Phase::ZnS)
            return "ZnS";
        if (phase == DetectorConstruction::Phase::Matrix)
            return "binder_void";
        if (phase == DetectorConstruction::Phase::World)
            return "outside";

        return "other";
    }

    G4bool IsOutsidePhase(const std::string &phase)
    {
        return (phase == "outside");
    }

    G4double PatchHalfZUm(const DetectorConstruction *det)
    {
        return det->GetBoxHalfZ() / um;
    }

    std::string ExitFaceLabel(const G4ThreeVector &position, const DetectorConstruction *det)
    {
        const G4double tolUm = 1.0e-3;
        const G4double xUm = position.x() / um;
        const G4double yUm = position.y() / um;
        const G4double zUm = position.z() / um;
        const G4double halfXUm = det->GetBoxHalfX() / um;
        const G4double halfYUm = det->GetBoxHalfY() / um;
        const G4double halfZUm = PatchHalfZUm(det);

        if (zUm >= halfZUm - tolUm)
            return "+Z";
        if (zUm <= -halfZUm + tolUm)
            return "-Z";
        if (xUm >= halfXUm - tolUm)
            return "+X";
        if (xUm <= -halfXUm + tolUm)
            return "-X";
        if (yUm >= halfYUm - tolUm)
            return "+Y";
        if (yUm <= -halfYUm + tolUm)
            return "-Y";
        return "unknown";
    }

    RunAction::BoundaryExitClass ClassifyBoundaryExit(
        const std::string &surfaceMode,
        const std::string &exitFace)
    {
        if (surfaceMode == "front_surface" && exitFace == "+Z")
        {
            return RunAction::BoundaryExitClass::PhysicalSurfaceExit;
        }
        if (surfaceMode == "back_surface" && exitFace == "-Z")
        {
            return RunAction::BoundaryExitClass::PhysicalSurfaceExit;
        }
        return RunAction::BoundaryExitClass::UnexpectedArtificialExit;
    }

    G4bool ContinueAcrossPeriodicBoundary(
        G4Track *track,
        const DetectorConstruction *detector,
        const G4ThreeVector &exitPosition,
        const G4ThreeVector &direction)
    {
        if (track == nullptr || detector == nullptr || direction.mag2() <= 0.0)
            return false;

        auto *stackManager = G4EventManager::GetEventManager()->GetStackManager();
        if (stackManager == nullptr)
            return false;

        auto *dynamicParticle = new G4DynamicParticle(*track->GetDynamicParticle());
        dynamicParticle->SetMomentumDirection(direction.unit());

        auto *continuationTrack = new G4Track(
            dynamicParticle,
            track->GetGlobalTime(),
            detector->WrapToPrimaryCellInside(exitPosition, direction));
        const auto *previousInfo = GetPeriodicTrackInfo(track);
        const G4int originalTrackId = previousInfo != nullptr
                                          ? previousInfo->OriginalTrackId()
                                          : track->GetTrackID();
        const G4int completedSteps =
            (previousInfo != nullptr ? previousInfo->CompletedSteps() : 0) +
            track->GetCurrentStepNumber();
        continuationTrack->SetTrackID(originalTrackId);
        continuationTrack->SetParentID(track->GetParentID());
        continuationTrack->SetLocalTime(track->GetLocalTime());
        continuationTrack->SetProperTime(track->GetProperTime());
        continuationTrack->SetWeight(track->GetWeight());
        continuationTrack->SetUserInformation(
            new PeriodicTrackInfo(originalTrackId, completedSteps));
        stackManager->PushOneTrack(continuationTrack);
        track->SetTrackStatus(fStopAndKill);
        return true;
    }
}

// --------------------------------------------------------------------

SteppingAction::SteppingAction(EventAction *eventAction,
                               const PrimaryGeneratorAction *primaryAction)
    : G4UserSteppingAction(),
      fEventAction(eventAction),
      fPrimaryAction(primaryAction)
{
}

// --------------------------------------------------------------------

SteppingAction::~SteppingAction() = default;

// --------------------------------------------------------------------

void SteppingAction::UserSteppingAction(const G4Step *step)
{
    if (!step)
        return;

    auto *track = step->GetTrack();
    if (!track)
        return;

    if (!IsTrackedHeavyParticle(track))
        return;

    const auto *prePoint = step->GetPreStepPoint();
    const auto *postPoint = step->GetPostStepPoint();
    if (!prePoint || !postPoint)
        return;

    const auto *prePV = prePoint->GetPhysicalVolume();
    const auto *postPV = postPoint->GetPhysicalVolume();

    const auto *detector = dynamic_cast<const DetectorConstruction *>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    const std::string phasePre = PhaseLabel(prePV, detector);
    const std::string phasePost = PhaseLabel(postPV, detector);

    // Do not keep tracking into world vacuum after leaving the patch.
    // But record the boundary-crossing step itself.
    const G4ThreeVector &xPre = prePoint->GetPosition();
    const G4ThreeVector &xPost = postPoint->GetPosition();

    const G4double stepLen = step->GetStepLength();
    const G4double edep = step->GetTotalEnergyDeposit();
    const G4double ekinPre = prePoint->GetKineticEnergy();
    const G4double ekinPost = postPoint->GetKineticEnergy();

    if (fEventAction)
    {
        fEventAction->AddEdep(edep);
    }

    auto *runAction = fEventAction ? fEventAction->GetRunAction() : nullptr;
    if (runAction && fPrimaryAction)
    {
        runAction->SwitchOutputCsvForInputPath(fPrimaryAction->GetCurrentRecordInputFile());
    }

    if (runAction && fPrimaryAction && runAction->IsFullMode() && runAction->IsFullStepCsvOpen())
    {
        std::ofstream &csv = runAction->GetFullStepCsv();

        const auto &rec = fPrimaryAction->GetCurrentRecord();
        const auto &capturePos = fPrimaryAction->GetCurrentLocalCapturePosition();
        const auto &bnCenter = fPrimaryAction->GetCurrentSelectedBNCenter();
        std::string placementFile = "unknown";

        const auto *det = dynamic_cast<const DetectorConstruction *>(
            G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        if (det)
        {
            placementFile = det->GetLoadedPlacementFileForRecord();
        }

        csv
            << fPrimaryAction->MakeCurrentPhysicalEventUid() << ","
            << rec.eventID << ","
            << rec.thickness_um << ","
            << rec.bn_wt << ","
            << rec.zns_wt << ","
            << rec.capture_x_um << ","
            << rec.capture_y_um << ","
            << rec.source_x_um << ","
            << rec.source_y_um << ","
            << rec.depth_um << ","
            << placementFile << ","
            << capturePos.x() / um << ","
            << capturePos.y() / um << ","
            << capturePos.z() / um << ","
            << fPrimaryAction->GetCurrentSurfaceMode() << ","
            << fPrimaryAction->GetCurrentTargetLocalZ() / um << ","
            << fPrimaryAction->GetCurrentUsedLocalZ() / um << ","
            << bnCenter.x() / um << ","
            << bnCenter.y() / um << ","
            << bnCenter.z() / um << ","
            << fPrimaryAction->GetCurrentAlphaLiReplayIndex() << ","
            << fPrimaryAction->GetCurrentAlphaLiReplayCount() << ","
            << OutputTrackId(track) << ","
            << OutputStepId(track) << ","
            << ParticleLabel(track) << ","
            << phasePre << ","
            << phasePost << ","
            << xPre.x() / um << ","
            << xPre.y() / um << ","
            << xPre.z() / um << ","
            << xPost.x() / um << ","
            << xPost.y() / um << ","
            << xPost.z() / um << ","
            << stepLen / um << ","
            << edep / keV << ","
            << ekinPre / keV << ","
            << ekinPost / keV << ","
            << fPrimaryAction->MakeCurrentSourceEventUid() << ","
            << rec.record_index << ","
            << fPrimaryAction->GetCurrentTrajectoryWeight()
            << "\n";
    }
    else if (runAction && fEventAction && runAction->IsSlimMode() &&
             runAction->IsSlimTrackCsvOpen() &&
             phasePre == "ZnS" && stepLen > 0.0)
    {
        std::ofstream &csv = runAction->GetSlimTrackCsv();
        const auto anchor = fEventAction->MakeCurrentCaptureAnchorRow();

        csv
            << anchor.physical_event_uid << ","
            << anchor.source_event_uid << ","
            << anchor.eventID << ","
            << anchor.record_index << ","
            << OutputTrackId(track) << ","
            << OutputStepId(track) << ","
            << ParticleLabel(track) << ","
            << phasePre << ","
            << phasePost << ","
            << xPre.x() / um << ","
            << xPre.y() / um << ","
            << xPre.z() / um << ","
            << xPost.x() / um << ","
            << xPost.y() / um << ","
            << xPost.z() / um << ","
            << stepLen / um << ","
            << edep / keV << ","
            << ekinPre / keV << ","
            << ekinPost / keV << ","
            << anchor.alphali_replay_index << ","
            << anchor.alphali_replay_count << ","
            << anchor.trajectory_weight
            << "\n";
    }

    // Preserve true screen-surface exits; wrap artificial RVE faces exactly.
    if (!IsOutsidePhase(phasePre) && IsOutsidePhase(phasePost))
    {
        const auto anchor = (fEventAction != nullptr)
                                ? fEventAction->MakeCurrentCaptureAnchorRow()
                                : RunAction::CaptureAnchorRow{};
        const std::string exitFace = detector != nullptr
                                         ? ExitFaceLabel(xPost, detector)
                                         : "unknown";
        const auto exitClass = ClassifyBoundaryExit(anchor.surface_mode, exitFace);
        const G4bool isPhysicalSurface =
            exitClass == RunAction::BoundaryExitClass::PhysicalSurfaceExit;

        if (!isPhysicalSurface && detector != nullptr &&
            detector->UsesPeriodicTransport() &&
            ContinueAcrossPeriodicBoundary(
                track, detector, xPost, prePoint->GetMomentumDirection()))
        {
            return;
        }

        if (runAction && fEventAction)
        {
            if (detector)
            {
                RunAction::UnexpectedBoundaryExitRow exitRow;
                exitRow.physical_event_uid = anchor.physical_event_uid;
                exitRow.source_event_uid = anchor.source_event_uid;
                exitRow.eventID = anchor.eventID;
                exitRow.record_index = anchor.record_index;
                exitRow.thickness_um = anchor.thickness_um;
                exitRow.bn_wt = anchor.bn_wt;
                exitRow.zns_wt = anchor.zns_wt;
                exitRow.placement_file = anchor.placement_file;
                exitRow.surface_mode = anchor.surface_mode;
                exitRow.particle = ParticleLabel(track);
                exitRow.trackID = OutputTrackId(track);
                exitRow.stepID = OutputStepId(track);
                exitRow.phase_pre = phasePre;
                exitRow.phase_post = phasePost;
                exitRow.exit_face = exitFace;
                exitRow.exit_class =
                    (exitClass == RunAction::BoundaryExitClass::PhysicalSurfaceExit)
                        ? "physical_surface_exit"
                        : "unexpected_artificial_exit";
                exitRow.ekin_post_keV = ekinPost / keV;
                exitRow.x_post_um = xPost.x() / um;
                exitRow.y_post_um = xPost.y() / um;
                exitRow.z_post_um = xPost.z() / um;
                exitRow.alphali_replay_index = anchor.alphali_replay_index;
                exitRow.alphali_replay_count = anchor.alphali_replay_count;
                exitRow.trajectory_weight = anchor.trajectory_weight;
                runAction->RecordBoundaryExit(exitRow, exitClass);
            }
        }
        track->SetTrackStatus(fStopAndKill);
    }
}
