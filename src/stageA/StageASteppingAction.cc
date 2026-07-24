#include "StageASteppingAction.hh"

#include "StageARunAction.hh"
#include "AnalysisConfig.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4VProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4Exception.hh"
#include "G4ios.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4TrackStatus.hh"
#include "G4RunManager.hh"

namespace
{
    G4String PhaseLabel(const G4VPhysicalVolume *pv, const DetectorConstruction *detector)
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

    const char *TrackStatusLabel(G4TrackStatus st)
    {
        switch (st)
        {
        case fAlive:
            return "fAlive";
        case fStopButAlive:
            return "fStopButAlive";
        case fStopAndKill:
            return "fStopAndKill";
        case fKillTrackAndSecondaries:
            return "fKillTrackAndSecondaries";
        case fSuspend:
            return "fSuspend";
        case fPostponeToNextEvent:
            return "fPostponeToNextEvent";
        default:
            return "unknown";
        }
    }
}

StageASteppingAction::StageASteppingAction(StageARunAction *runAction, AnalysisConfig *config)
    : G4UserSteppingAction(),
      fRunAction(runAction),
      fConfig(config)
{
    if (fRunAction == nullptr)
    {
        G4Exception("StageASteppingAction::StageASteppingAction",
                    "BNZS_A_STEP_001", FatalException,
                    "StageARunAction pointer is null.");
        return;
    }

    if (fConfig == nullptr)
    {
        G4Exception("StageASteppingAction::StageASteppingAction",
                    "BNZS_A_STEP_002", FatalException,
                    "AnalysisConfig pointer is null.");
        return;
    }
}

StageASteppingAction::~StageASteppingAction() = default;

void StageASteppingAction::UserSteppingAction(const G4Step *step)
{
    if (step == nullptr || fRunAction == nullptr || fConfig == nullptr)
    {
        return;
    }

    if (!IsPrimaryNeutronStep(step))
    {
        return;
    }

    const G4Track *track = step->GetTrack();
    const G4StepPoint *prePoint = step->GetPreStepPoint();
    const G4StepPoint *postPoint = step->GetPostStepPoint();

    if (track == nullptr || prePoint == nullptr || postPoint == nullptr)
    {
        return;
    }

    static G4int dbgStepPrintCount = 0;
    if (dbgStepPrintCount < 10)
    {
        ++dbgStepPrintCount;

        G4cout << "[StageA-step-debug]"
               << " trackID=" << track->GetTrackID()
               << " stepNo=" << track->GetCurrentStepNumber()
               << " preZ_um=" << prePoint->GetPosition().z() / um
               << " postZ_um=" << postPoint->GetPosition().z() / um
               << " stepLen_um=" << step->GetStepLength() / um
               << G4endl;
    }

    // ---- count one incident neutron only for primary neutron track ----
    if (track->GetParentID() == 0 && track->GetCurrentStepNumber() == 1)
    {
        fRunAction->AddIncident(1);
        fRunAction->MaybePrintProgress("incident");
    }

    // ---- accumulate total neutron track length inside fixed microstructure patch ----
    const G4ThreeVector midPoint =
        0.5 * (prePoint->GetPosition() + postPoint->GetPosition());

    const G4double stepLength = step->GetStepLength();
    if (stepLength > 0.0 && IsPointInsidePatch(midPoint))
    {
        fRunAction->AddTrackLength(stepLength);
    }

    const G4VProcess *proc = postPoint->GetProcessDefinedStep();
    const G4String procName = (proc != nullptr) ? proc->GetProcessName() : "none";
    const auto *detector = dynamic_cast<const DetectorConstruction *>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    const G4String prePhase = PhaseLabel(prePoint->GetPhysicalVolume(), detector);
    const G4String postPhase = PhaseLabel(postPoint->GetPhysicalVolume(), detector);

    // ---- count neutron absorption-like termination in BN ----
    // 对 BN 中的 10B(n,alpha)7Li，HP 下可能表现为 nCaptureHP，
    // 也可能表现为 neutronInelastic / hadInelastic。
    // 这里先用“发生在 BN 相内的吸收样终止”来统计。
    const G4bool isCaptureLike =
        (procName == "nCapture" || procName == "nCaptureHP") ||
        ((procName == "neutronInelastic" || procName == "hadInelastic") &&
         prePhase == "BN");

    if (isCaptureLike)
    {
        fRunAction->AddCapture(1);

        if (fRunAction->GetCaptureCount() <= 10)
        {
            G4cout << "[StageA-capture-debug]"
                   << " proc=" << procName
                   << " stepNo=" << track->GetCurrentStepNumber()
                   << " prePhase=" << prePhase
                   << " postPhase=" << postPhase
                   << " postZ_um=" << postPoint->GetPosition().z() / um
                   << G4endl;

            fRunAction->ForcePrintProgress("capture");
        }
    }

    // ---- diagnose how neutron tracks terminate / leave ----
    // 只打印前 30 个“neutron 终止/离开”案例，看看大多数到底是 capture、散射后出界，还是别的。
    static G4int dbgTerminalPrintCount = 0;
    const G4bool leftWorld =
        (postPoint->GetPhysicalVolume() == nullptr) ||
        (postPoint->GetStepStatus() == fWorldBoundary);

    const G4bool terminal =
        leftWorld || (track->GetTrackStatus() != fAlive);

    if (terminal && dbgTerminalPrintCount < 30)
    {
        ++dbgTerminalPrintCount;

        G4cout << "[StageA-terminal-debug]"
               << " proc=" << procName
               << " trackStatus=" << TrackStatusLabel(track->GetTrackStatus())
               << " stepNo=" << track->GetCurrentStepNumber()
               << " prePhase=" << prePhase
               << " postPhase=" << postPhase
               << " preZ_um=" << prePoint->GetPosition().z() / um
               << " postZ_um=" << postPoint->GetPosition().z() / um
               << " stepLen_um=" << step->GetStepLength() / um
               << G4endl;
    }
}

G4bool StageASteppingAction::IsPrimaryNeutronStep(const G4Step *step) const
{
    if (step == nullptr)
    {
        return false;
    }

    const G4Track *track = step->GetTrack();
    if (track == nullptr)
    {
        return false;
    }

    const G4ParticleDefinition *particle = track->GetParticleDefinition();
    if (particle == nullptr)
    {
        return false;
    }

    // 关键修改：
    // Stage A 统计所有 neutron track 的路径长度与 capture，
    // 不再把 secondary neutron 过滤掉。
    if (particle->GetParticleName() != "neutron")
    {
        return false;
    }

    return true;
}

G4bool StageASteppingAction::IsPointInsidePatch(const G4ThreeVector &p) const
{
    const auto *detector = dynamic_cast<const DetectorConstruction *>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    if (detector == nullptr)
        return false;
    return (p.x() >= -detector->GetBoxHalfX() && p.x() <= detector->GetBoxHalfX() &&
            p.y() >= -detector->GetBoxHalfY() && p.y() <= detector->GetBoxHalfY() &&
            p.z() >= -detector->GetBoxHalfZ() && p.z() <= detector->GetBoxHalfZ());
}
