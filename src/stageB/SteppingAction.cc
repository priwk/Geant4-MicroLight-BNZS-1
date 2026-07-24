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
        PeriodicTrackInfo(G4int originalTrackId,
                          G4int completedSteps,
                          G4int cellIx,
                          G4int cellIy,
                          G4int cellIz)
            : fOriginalTrackId(originalTrackId),
              fCompletedSteps(completedSteps),
              fCellIx(cellIx),
              fCellIy(cellIy),
              fCellIz(cellIz)
        {
        }

        void Print() const override {}
        G4int OriginalTrackId() const { return fOriginalTrackId; }
        G4int CompletedSteps() const { return fCompletedSteps; }
        G4int CellIx() const { return fCellIx; }
        G4int CellIy() const { return fCellIy; }
        G4int CellIz() const { return fCellIz; }

    private:
        G4int fOriginalTrackId;
        G4int fCompletedSteps;
        G4int fCellIx;
        G4int fCellIy;
        G4int fCellIz;
    };

    struct PeriodicCell
    {
        G4int ix = 0;
        G4int iy = 0;
        G4int iz = 0;
    };

    struct CoordinateFrame
    {
        PeriodicCell cell;
        G4ThreeVector unwrappedPre;
        G4ThreeVector unwrappedPost;
        G4ThreeVector screenPre;
        G4ThreeVector screenPost;
        G4double screenDepthPre = 0.0;
        G4double screenDepthPost = 0.0;
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
        return info != nullptr && info->OriginalTrackId() >= 0
                   ? info->OriginalTrackId()
                   : track->GetTrackID();
    }

    G4int OutputStepId(const G4Track *track)
    {
        const auto *info = GetPeriodicTrackInfo(track);
        return (info != nullptr ? info->CompletedSteps() : 0) +
               track->GetCurrentStepNumber();
    }

    PeriodicCell CurrentCell(const G4Track *track)
    {
        const auto *info = GetPeriodicTrackInfo(track);
        if (info == nullptr)
            return {};
        return {info->CellIx(), info->CellIy(), info->CellIz()};
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

    PeriodicCell ExitCellShift(const G4ThreeVector &position,
                               const G4ThreeVector &direction,
                               const DetectorConstruction *det)
    {
        const G4ThreeVector wrapped =
            det->WrapToPrimaryCellInside(position, direction);
        const auto imageShift = [](G4double before,
                                   G4double after,
                                   G4double length)
        {
            if (length <= 0.0)
                return 0;
            return static_cast<G4int>(std::llround((before - after) / length));
        };
        PeriodicCell shift;
        shift.ix = imageShift(position.x(), wrapped.x(), det->GetBoxXUm() * um);
        shift.iy = imageShift(position.y(), wrapped.y(), det->GetBoxYUm() * um);
        shift.iz = imageShift(position.z(), wrapped.z(), det->GetBoxZUm() * um);
        return shift;
    }

    std::string ExitFaceLabel(const PeriodicCell &shift)
    {
        if (shift.iz > 0)
            return "+Z";
        if (shift.iz < 0)
            return "-Z";
        if (shift.ix > 0)
            return "+X";
        if (shift.ix < 0)
            return "-X";
        if (shift.iy > 0)
            return "+Y";
        if (shift.iy < 0)
            return "-Y";
        return "unknown";
    }

    RunAction::BoundaryExitClass ClassifyBoundaryExit(
        G4double screenThicknessUm,
        G4double screenDepthPostUm,
        const PeriodicCell &exitShift,
        const DetectorConstruction *detector)
    {
        if (detector == nullptr || exitShift.iz == 0)
        {
            return RunAction::BoundaryExitClass::UnexpectedArtificialExit;
        }
        const G4double surfaceToleranceUm = 1.0e-3;
        if (screenDepthPostUm <= surfaceToleranceUm ||
            screenDepthPostUm >= screenThicknessUm - surfaceToleranceUm)
        {
            return RunAction::BoundaryExitClass::PhysicalSurfaceExit;
        }
        return RunAction::BoundaryExitClass::UnexpectedArtificialExit;
    }

    CoordinateFrame MakeCoordinateFrame(
        const G4Track *track,
        const DetectorConstruction *detector,
        const RunAction::CaptureAnchorRow &anchor,
        const G4ThreeVector &localPre,
        const G4ThreeVector &localPost)
    {
        CoordinateFrame frame;
        frame.cell = CurrentCell(track);
        const G4ThreeVector cellOffset(
            frame.cell.ix * detector->GetBoxXUm() * um,
            frame.cell.iy * detector->GetBoxYUm() * um,
            frame.cell.iz * detector->GetBoxZUm() * um);
        frame.unwrappedPre = localPre + cellOffset;
        frame.unwrappedPost = localPost + cellOffset;

        const G4ThreeVector localCapture(
            anchor.local_capture_x_um * um,
            anchor.local_capture_y_um * um,
            anchor.local_capture_z_um * um);
        const G4double captureScreenZ =
            (0.5 * anchor.thickness_um - anchor.depth_um) * um;
        const G4ThreeVector screenCapture(
            anchor.capture_x_um * um,
            anchor.capture_y_um * um,
            captureScreenZ);
        frame.screenPre = screenCapture + frame.unwrappedPre - localCapture;
        frame.screenPost = screenCapture + frame.unwrappedPost - localCapture;
        frame.screenDepthPre =
            anchor.depth_um * um - (frame.unwrappedPre.z() - localCapture.z());
        frame.screenDepthPost =
            anchor.depth_um * um - (frame.unwrappedPost.z() - localCapture.z());
        return frame;
    }

    void AppendCoordinateFrame(std::ofstream &csv, const CoordinateFrame &frame)
    {
        csv << "," << frame.cell.ix
            << "," << frame.cell.iy
            << "," << frame.cell.iz
            << "," << frame.unwrappedPre.x() / um
            << "," << frame.unwrappedPre.y() / um
            << "," << frame.unwrappedPre.z() / um
            << "," << frame.unwrappedPost.x() / um
            << "," << frame.unwrappedPost.y() / um
            << "," << frame.unwrappedPost.z() / um
            << "," << frame.screenPre.x() / um
            << "," << frame.screenPre.y() / um
            << "," << frame.screenPre.z() / um
            << "," << frame.screenDepthPre / um
            << "," << frame.screenPost.x() / um
            << "," << frame.screenPost.y() / um
            << "," << frame.screenPost.z() / um
            << "," << frame.screenDepthPost / um;
    }

    void PropagateCellToSecondaries(const G4Step *step, const PeriodicCell &cell)
    {
        const auto *secondaries = step->GetSecondaryInCurrentStep();
        if (secondaries == nullptr)
            return;
        for (const auto *secondary : *secondaries)
        {
            if (secondary == nullptr || secondary->GetUserInformation() != nullptr)
                continue;
            secondary->SetUserInformation(
                new PeriodicTrackInfo(-1, 0, cell.ix, cell.iy, cell.iz));
        }
    }

    G4bool ContinueAcrossPeriodicBoundary(
        G4Track *track,
        const DetectorConstruction *detector,
        const G4ThreeVector &exitPosition,
        const G4ThreeVector &direction,
        const PeriodicCell &exitShift)
    {
        if (track == nullptr || detector == nullptr || direction.mag2() <= 0.0)
            return false;

        auto *stackManager = G4EventManager::GetEventManager()->GetStackManager();
        if (stackManager == nullptr)
            return false;

        auto *dynamicParticle = new G4DynamicParticle(*track->GetDynamicParticle());
        dynamicParticle->SetMomentumDirection(direction.unit());

        const G4double inwardEpsilon = 1.0e-4 * um;
        const G4ThreeVector continuationPosition(
            exitPosition.x() - exitShift.ix * detector->GetBoxXUm() * um +
                exitShift.ix * inwardEpsilon,
            exitPosition.y() - exitShift.iy * detector->GetBoxYUm() * um +
                exitShift.iy * inwardEpsilon,
            exitPosition.z() - exitShift.iz * detector->GetBoxZUm() * um +
                exitShift.iz * inwardEpsilon);
        auto *continuationTrack = new G4Track(
            dynamicParticle,
            track->GetGlobalTime(),
            continuationPosition);
        const auto *previousInfo = GetPeriodicTrackInfo(track);
        const G4int originalTrackId = OutputTrackId(track);
        const G4int completedSteps =
            (previousInfo != nullptr ? previousInfo->CompletedSteps() : 0) +
            track->GetCurrentStepNumber();
        const PeriodicCell previousCell = CurrentCell(track);
        continuationTrack->SetTrackID(originalTrackId);
        continuationTrack->SetParentID(track->GetParentID());
        continuationTrack->SetLocalTime(track->GetLocalTime());
        continuationTrack->SetProperTime(track->GetProperTime());
        continuationTrack->SetWeight(track->GetWeight());
        continuationTrack->SetUserInformation(
            new PeriodicTrackInfo(
                originalTrackId,
                completedSteps,
                previousCell.ix + exitShift.ix,
                previousCell.iy + exitShift.iy,
                previousCell.iz + exitShift.iz));
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
    const G4bool isHeavyParticle = IsTrackedHeavyParticle(track);
    const PeriodicCell currentCell = CurrentCell(track);

    PropagateCellToSecondaries(step, currentCell);

    if (fEventAction && phasePre == "ZnS")
    {
        fEventAction->AddEdep(edep);
    }

    auto *runAction = fEventAction ? fEventAction->GetRunAction() : nullptr;

    const auto anchor = (fEventAction != nullptr)
                            ? fEventAction->MakeCurrentCaptureAnchorRow()
                            : RunAction::CaptureAnchorRow{};
    const CoordinateFrame coordinateFrame = detector != nullptr
                                                ? MakeCoordinateFrame(
                                                      track, detector, anchor, xPre, xPost)
                                                : CoordinateFrame{};

    if (isHeavyParticle && runAction && fPrimaryAction &&
        runAction->IsFullMode() && runAction->IsFullStepCsvOpen())
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
            << fPrimaryAction->GetCurrentTrajectoryWeight() << ","
            << track->GetParentID();
        AppendCoordinateFrame(csv, coordinateFrame);
        csv << "\n";
    }
    if (runAction && fEventAction &&
             runAction->IsSlimTrackCsvOpen() &&
             phasePre == "ZnS" &&
             ((isHeavyParticle && stepLen > 0.0) || edep > 0.0))
    {
        std::ofstream &csv = runAction->GetSlimTrackCsv();

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
            << anchor.trajectory_weight << ","
            << track->GetParentID();
        AppendCoordinateFrame(csv, coordinateFrame);
        csv << "\n";
    }

    // Preserve true screen-surface exits; wrap artificial RVE faces exactly.
    if (!IsOutsidePhase(phasePre) && IsOutsidePhase(phasePost))
    {
        const PeriodicCell exitShift = detector != nullptr
                                           ? ExitCellShift(
                                                 xPost,
                                                 prePoint->GetMomentumDirection(),
                                                 detector)
                                           : PeriodicCell{};
        const std::string exitFace = ExitFaceLabel(exitShift);
        const auto exitClass = ClassifyBoundaryExit(
            anchor.thickness_um,
            coordinateFrame.screenDepthPost / um,
            exitShift,
            detector);
        const G4bool isPhysicalSurface =
            exitClass == RunAction::BoundaryExitClass::PhysicalSurfaceExit;

        if (!isPhysicalSurface && detector != nullptr &&
            detector->UsesPeriodicTransport() &&
            ContinueAcrossPeriodicBoundary(
                track,
                detector,
                xPost,
                postPoint->GetMomentumDirection(),
                exitShift))
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
