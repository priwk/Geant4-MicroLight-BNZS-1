#include "StageDOpticalSteppingAction.hh"

#include "AnalysisConfig.hh"
#include "DetectorConstruction.hh"
#include "StageDOpticalEventAction.hh"
#include "StageDOpticalRunAction.hh"
#include "StageDReentrySampler.hh"

#include "G4EventManager.hh"
#include "G4Exception.hh"
#include "G4DynamicParticle.hh"
#include "G4GeometryTolerance.hh"
#include "G4OpticalPhoton.hh"
#include "G4LogicalVolume.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4RunManager.hh"
#include "G4StackManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4VProcess.hh"
#include "G4VPhysicalVolume.hh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace
{
  constexpr G4double kBoundaryEpsilon = 1.0e-4 * um;

  enum class RveFace
  {
    PosX,
    NegX,
    PosY,
    NegY,
    PosZ,
    NegZ,
    None
  };

  enum class OpticalBoundaryTransition
  {
    Reflection,
    Transmission,
    Absorption,
    Detection,
    StepTooSmall,
    Unavailable,
    Other
  };

  G4bool IsOpticalPhoton(const G4Track *track)
  {
    return track != nullptr &&
           track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition();
  }

  G4double AngleDeg(const G4ThreeVector &a, const G4ThreeVector &b)
  {
    const G4double dot = std::clamp(a.dot(b), -1.0, 1.0);
    return std::acos(dot) / deg;
  }

  G4double ClampCosTheta(const G4ThreeVector &a, const G4ThreeVector &b)
  {
    return std::clamp(a.dot(b), -1.0, 1.0);
  }

  G4bool UsesAngleThresholdEncounterMetric(const AnalysisConfig *config)
  {
    return config != nullptr &&
           config->stageD_scatter_metric == "particle_encounter_angle_threshold";
  }

  G4int PhaseFunctionBin(const G4double cosTheta);

  std::string ProcessName(const G4StepPoint *point)
  {
    if (point == nullptr)
      return "";
    const auto *process = point->GetProcessDefinedStep();
    if (process == nullptr)
      return "";
    return process->GetProcessName();
  }

  DetectorConstruction::Phase PhaseFromPhysicalVolume(
      const G4VPhysicalVolume *pv,
      const DetectorConstruction *detector)
  {
    if (pv == nullptr)
      return DetectorConstruction::Phase::World;

    const auto *lv = pv->GetLogicalVolume();
    if (lv == nullptr)
      return DetectorConstruction::Phase::World;

    return detector != nullptr
               ? detector->GetPhaseFromLogicalVolume(lv)
               : DetectorConstruction::Phase::Unknown;
  }

  G4bool IsOnRveSurface(const G4ThreeVector &position,
                        const DetectorConstruction *detector,
                        G4double tolerance,
                        RveFace *face)
  {
    if (detector == nullptr)
      return false;

    const G4double halfX = detector->GetPatchHalfXUm() * um;
    const G4double halfY = detector->GetPatchHalfYUm() * um;
    const G4double halfZ = detector->GetPatchHalfZUm() * um;
    const G4double tol = std::max(tolerance, 1.0e-9 * um);

    const G4bool withinX = std::abs(position.x()) <= halfX + tol;
    const G4bool withinY = std::abs(position.y()) <= halfY + tol;
    const G4bool withinZ = std::abs(position.z()) <= halfZ + tol;
    if (!withinX || !withinY || !withinZ)
      return false;

    struct FaceCandidate
    {
      RveFace face = RveFace::None;
      G4double delta = 0.0;
    };

    FaceCandidate best{RveFace::None, std::numeric_limits<G4double>::infinity()};
    auto consider = [&](RveFace candidate, G4double delta)
    {
      if (delta <= tol && delta < best.delta)
        best = FaceCandidate{candidate, delta};
    };

    consider((position.x() >= 0.0) ? RveFace::PosX : RveFace::NegX,
             std::abs(std::abs(position.x()) - halfX));
    consider((position.y() >= 0.0) ? RveFace::PosY : RveFace::NegY,
             std::abs(std::abs(position.y()) - halfY));
    consider((position.z() >= 0.0) ? RveFace::PosZ : RveFace::NegZ,
             std::abs(std::abs(position.z()) - halfZ));

    if (best.face == RveFace::None)
      return false;
    if (face != nullptr)
      *face = best.face;
    return true;
  }

  G4bool IsInsideRveBoxWithTolerance(const G4ThreeVector &position,
                                     const DetectorConstruction *detector,
                                     G4double tolerance)
  {
    if (detector == nullptr)
      return false;

    const G4double halfX = detector->GetPatchHalfXUm() * um;
    const G4double halfY = detector->GetPatchHalfYUm() * um;
    const G4double halfZ = detector->GetPatchHalfZUm() * um;
    const G4double tol = std::max(tolerance, 1.0e-9 * um);
    return std::abs(position.x()) <= halfX + tol &&
           std::abs(position.y()) <= halfY + tol &&
           std::abs(position.z()) <= halfZ + tol;
  }

  G4ThreeVector RveFaceNormal(RveFace face)
  {
    switch (face)
    {
    case RveFace::PosX:
      return G4ThreeVector(1.0, 0.0, 0.0);
    case RveFace::NegX:
      return G4ThreeVector(-1.0, 0.0, 0.0);
    case RveFace::PosY:
      return G4ThreeVector(0.0, 1.0, 0.0);
    case RveFace::NegY:
      return G4ThreeVector(0.0, -1.0, 0.0);
    case RveFace::PosZ:
      return G4ThreeVector(0.0, 0.0, 1.0);
    case RveFace::NegZ:
      return G4ThreeVector(0.0, 0.0, -1.0);
    case RveFace::None:
    default:
      return G4ThreeVector();
    }
  }

  G4bool IsRveOuterBoundaryStep(const G4Step *step,
                                const DetectorConstruction *detector,
                                RveFace *face)
  {
    if (step == nullptr || detector == nullptr ||
        step->GetPreStepPoint() == nullptr ||
        step->GetPostStepPoint() == nullptr)
      return false;

    const auto *prePoint = step->GetPreStepPoint();
    const auto *postPoint = step->GetPostStepPoint();
    if (postPoint->GetStepStatus() != fGeomBoundary)
      return false;

    const G4double tolerance =
        G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
    if (!IsInsideRveBoxWithTolerance(prePoint->GetPosition(), detector, tolerance))
      return false;

    RveFace candidateFace = RveFace::None;
    if (!IsOnRveSurface(postPoint->GetPosition(), detector, tolerance, &candidateFace))
      return false;

    const G4ThreeVector chord = postPoint->GetPosition() - prePoint->GetPosition();
    const G4ThreeVector normal = RveFaceNormal(candidateFace);
    if (chord.mag2() > 0.0 && normal.mag2() > 0.0 &&
        chord.dot(normal) < -tolerance)
    {
      return false;
    }
    const G4ThreeVector preDir = prePoint->GetMomentumDirection();
    if (preDir.mag2() > 0.0 && normal.mag2() > 0.0 &&
        preDir.unit().dot(normal) < -1.0e-9)
    {
      return false;
    }

    if (face != nullptr)
      *face = candidateFace;
    return true;
  }

  G4bool IsParticlePhase(DetectorConstruction::Phase phase)
  {
    return phase == DetectorConstruction::Phase::BN ||
           phase == DetectorConstruction::Phase::ZnS;
  }

  G4bool IsReentryPhase(DetectorConstruction::Phase phase)
  {
    return IsParticlePhase(phase) ||
           phase == DetectorConstruction::Phase::Matrix;
  }

  std::string PhaseLabel(DetectorConstruction::Phase phase)
  {
    if (phase == DetectorConstruction::Phase::Unknown)
      return "Unknown";
    return DetectorConstruction::PhaseName(phase);
  }

  DetectorConstruction::Phase PhaseFromLabel(const std::string &label)
  {
    if (label == "BN")
      return DetectorConstruction::Phase::BN;
    if (label == "ZnS")
      return DetectorConstruction::Phase::ZnS;
    if (label == "Matrix")
      return DetectorConstruction::Phase::Matrix;
    if (label == "World")
      return DetectorConstruction::Phase::World;
    return DetectorConstruction::Phase::Unknown;
  }

  G4bool IsBoundaryReflection(G4OpBoundaryProcessStatus status)
  {
    return status == FresnelReflection ||
           status == TotalInternalReflection ||
           status == LambertianReflection ||
           status == LobeReflection ||
           status == SpikeReflection ||
           status == BackScattering ||
           status == PolishedLumirrorAirReflection ||
           status == PolishedLumirrorGlueReflection ||
           status == PolishedAirReflection ||
           status == PolishedTeflonAirReflection ||
           status == PolishedTiOAirReflection ||
           status == PolishedTyvekAirReflection ||
           status == PolishedVM2000AirReflection ||
           status == PolishedVM2000GlueReflection ||
           status == EtchedLumirrorAirReflection ||
           status == EtchedLumirrorGlueReflection ||
           status == EtchedAirReflection ||
           status == EtchedTeflonAirReflection ||
           status == EtchedTiOAirReflection ||
           status == EtchedTyvekAirReflection ||
           status == EtchedVM2000AirReflection ||
           status == EtchedVM2000GlueReflection ||
           status == GroundLumirrorAirReflection ||
           status == GroundLumirrorGlueReflection ||
           status == GroundAirReflection ||
           status == GroundTeflonAirReflection ||
           status == GroundTiOAirReflection ||
           status == GroundTyvekAirReflection ||
           status == GroundVM2000AirReflection ||
           status == GroundVM2000GlueReflection ||
           status == CoatedDielectricReflection;
  }

  OpticalBoundaryTransition ClassifyBoundaryTransition(
      G4OpBoundaryProcessStatus status)
  {
    if (IsBoundaryReflection(status))
      return OpticalBoundaryTransition::Reflection;
    if (status == Transmission ||
        status == FresnelRefraction ||
        status == SameMaterial ||
        status == CoatedDielectricRefraction ||
        status == CoatedDielectricFrustratedTransmission)
    {
      return OpticalBoundaryTransition::Transmission;
    }
    if (status == Absorption)
      return OpticalBoundaryTransition::Absorption;
    if (status == Detection)
      return OpticalBoundaryTransition::Detection;
    if (status == StepTooSmall)
      return OpticalBoundaryTransition::StepTooSmall;
    if (status == Undefined || status == NotAtBoundary)
      return OpticalBoundaryTransition::Unavailable;
    return OpticalBoundaryTransition::Other;
  }

  G4bool IsPhaseTransmission(DetectorConstruction::Phase prePhase,
                             DetectorConstruction::Phase postPhase)
  {
    return prePhase != postPhase &&
           prePhase != DetectorConstruction::Phase::Unknown &&
           postPhase != DetectorConstruction::Phase::Unknown &&
           prePhase != DetectorConstruction::Phase::World &&
           postPhase != DetectorConstruction::Phase::World;
  }

  void AccumulateOuterBoundaryStatus(StageDPhotonEventRecord &event,
                                     G4OpBoundaryProcessStatus status)
  {
    if (status == FresnelReflection)
      ++event.num_outer_boundary_fresnel_reflection;
    else if (status == TotalInternalReflection)
      ++event.num_outer_boundary_total_internal_reflection;
    else if (status == FresnelRefraction || status == CoatedDielectricRefraction)
      ++event.num_outer_boundary_refraction;
    else if (status == Transmission ||
             status == SameMaterial ||
             status == CoatedDielectricFrustratedTransmission)
      ++event.num_outer_boundary_transmission;
    else
      ++event.num_outer_boundary_other_status;
  }

  G4bool ComputeExitPointOnRveBox(const G4ThreeVector &prePos,
                                  const G4ThreeVector &oldDir,
                                  const DetectorConstruction *detector,
                                  G4ThreeVector &exitPoint)
  {
    if (detector == nullptr || oldDir.mag2() <= 0.0)
      return false;

    const G4double halfX = detector->GetPatchHalfXUm() * um;
    const G4double halfY = detector->GetPatchHalfYUm() * um;
    const G4double halfZ = detector->GetPatchHalfZUm() * um;

    G4double tMin = -std::numeric_limits<G4double>::infinity();
    G4double tMax = std::numeric_limits<G4double>::infinity();

    auto updateAxis = [&](G4double pos, G4double dir, G4double half) -> G4bool
    {
      if (std::abs(dir) <= 1.0e-18)
        return std::abs(pos) <= half;

      G4double t1 = (-half - pos) / dir;
      G4double t2 = (+half - pos) / dir;
      if (t1 > t2)
        std::swap(t1, t2);
      tMin = std::max(tMin, t1);
      tMax = std::min(tMax, t2);
      return tMin <= tMax;
    };

    if (!updateAxis(prePos.x(), oldDir.x(), halfX) ||
        !updateAxis(prePos.y(), oldDir.y(), halfY) ||
        !updateAxis(prePos.z(), oldDir.z(), halfZ))
    {
      return false;
    }

    if (tMax <= 0.0 || !std::isfinite(tMax))
      return false;

    exitPoint = prePos + tMax * oldDir;
    return true;
  }

  DetectorConstruction::Phase ProbeForwardPhase(const G4Step *step,
                                                const DetectorConstruction *detector,
                                                const G4ThreeVector &direction)
  {
    if (step == nullptr || detector == nullptr || step->GetPostStepPoint() == nullptr ||
        direction.mag2() <= 0.0)
    {
      return DetectorConstruction::Phase::Unknown;
    }

    const G4double tolerance =
        G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
    const G4double offset = std::max(4.0 * tolerance, 4.0 * kBoundaryEpsilon);
    return detector->FindPeriodicPhaseAtPoint(step->GetPostStepPoint()->GetPosition() +
                                              offset * direction.unit());
  }

  void RecordCompleteEncounter(StageDPhotonEventRecord &event,
                               G4double cosTheta,
                               DetectorConstruction::Phase particlePhase,
                               const AnalysisConfig *config)
  {
    cosTheta = std::clamp(cosTheta, -1.0, 1.0);
    const G4double oneMinusCosTheta = 1.0 - cosTheta;
    const G4double cos2Theta = cosTheta * cosTheta;
    const G4double thetaDeg = std::acos(cosTheta) / deg;
    const G4bool useThresholdedEncounterMetric =
        UsesAngleThresholdEncounterMetric(config);
    const G4bool passesEncounterThreshold =
        (config != nullptr) ? (thetaDeg >= config->stageD_theta_threshold_deg) : true;

    ++event.num_complete_encounter_total;
    ++event.num_encounter_total;
    event.sum_cos_theta_encounter += cosTheta;
    event.sum_one_minus_cos_theta_encounter += oneMinusCosTheta;
    event.sum_cos2_theta_encounter += cos2Theta;
    if (passesEncounterThreshold)
    {
      ++event.num_encounter_effective_total;
      event.sum_cos_theta_encounter_effective += cosTheta;
      event.sum_one_minus_cos_theta_encounter_effective += oneMinusCosTheta;
      event.sum_cos2_theta_encounter_effective += cos2Theta;
    }
    if (!useThresholdedEncounterMetric || passesEncounterThreshold)
      ++event.phase_function_histogram[PhaseFunctionBin(cosTheta)];

    ++event.num_particle_scatter;
    event.sum_cos_theta_particle += cosTheta;

    if (particlePhase == DetectorConstruction::Phase::BN)
    {
      ++event.num_complete_encounter_BN;
      ++event.num_encounter_BN;
      event.sum_cos_theta_encounter_BN += cosTheta;
      event.sum_one_minus_cos_theta_encounter_BN += oneMinusCosTheta;
      event.sum_cos2_theta_encounter_BN += cos2Theta;
      if (passesEncounterThreshold)
      {
        ++event.num_encounter_effective_BN;
        event.sum_cos_theta_encounter_effective_BN += cosTheta;
        event.sum_one_minus_cos_theta_encounter_effective_BN += oneMinusCosTheta;
        event.sum_cos2_theta_encounter_effective_BN += cos2Theta;
      }
      ++event.num_particle_scatter_BN;
      event.sum_cos_theta_particle_BN += cosTheta;
    }
    else if (particlePhase == DetectorConstruction::Phase::ZnS)
    {
      ++event.num_complete_encounter_ZnS;
      ++event.num_encounter_ZnS;
      event.sum_cos_theta_encounter_ZnS += cosTheta;
      event.sum_one_minus_cos_theta_encounter_ZnS += oneMinusCosTheta;
      event.sum_cos2_theta_encounter_ZnS += cos2Theta;
      if (passesEncounterThreshold)
      {
        ++event.num_encounter_effective_ZnS;
        event.sum_cos_theta_encounter_effective_ZnS += cosTheta;
        event.sum_one_minus_cos_theta_encounter_effective_ZnS += oneMinusCosTheta;
        event.sum_cos2_theta_encounter_effective_ZnS += cos2Theta;
      }
      ++event.num_particle_scatter_ZnS;
      event.sum_cos_theta_particle_ZnS += cosTheta;
    }
  }

  StageDReentryDiagnosticRecord MakeReentryDiagnosticRecord(
      const StageDReentrySampler::ReentryContext &ctx,
      const StageDReentrySampler::ReentryDiagnostics &diag)
  {
    StageDReentryDiagnosticRecord record;
    record.event_id = ctx.eventID;
    record.reentry_index = ctx.reentryIndex;
    record.strategy = diag.strategy;
    record.fallback_level = diag.fallbackLevel;
    record.exit_phase = PhaseLabel(diag.exitPhase);
    record.entry_phase = PhaseLabel(diag.entryPhase);
    record.old_dir = ctx.oldDir;
    record.exit_point = ctx.exitPoint;
    record.entry_point = diag.entryPoint;
    record.particle_q_exit = diag.particleQExit;
    record.particle_q_entry = diag.particleQEntry;
    record.particle_mu_exit = diag.particleMuExit;
    record.particle_mu_entry = diag.particleMuEntry;
    record.matrix_clearance_exit_um = diag.matrixClearanceExitUm;
    record.matrix_clearance_entry_um = diag.matrixClearanceEntryUm;
    record.matrix_nearest_phase_exit = diag.matrixNearestPhaseExit;
    record.matrix_nearest_phase_entry = diag.matrixNearestPhaseEntry;
    record.matrix_clearance_bin_exit = diag.matrixClearanceBinExit;
    record.matrix_clearance_bin_entry = diag.matrixClearanceBinEntry;
    record.trials = diag.trials;
    return record;
  }

  void AccumulateReentryCounters(StageDPhotonEventRecord &event,
                                 const StageDReentrySampler::ReentryDiagnostics &diag)
  {
    if (diag.strategy == "particle_sphere_q_mu")
      ++event.num_reentry_particle_q_mu;
    else if (diag.strategy == "particle_sphere_q_only")
    {
      ++event.num_reentry_particle_q_only;
      if (diag.fallbackLevel == "q_only")
        ++event.num_reentry_particle_q_only_fallback;
    }
    else if (diag.strategy == "particle_same_phase_volume_random")
      ++event.num_reentry_particle_volume_random;
    else if (diag.strategy == "matrix_clearance_binned_portal")
      ++event.num_reentry_matrix_clearance_portal;
    else if (diag.strategy == "matrix_random_debug")
      ++event.num_reentry_random_matrix_debug;

    if (diag.fallbackLevel == "same_bin")
      ++event.num_reentry_fallback_same_bin;
    else if (diag.fallbackLevel == "adjacent_bin")
      ++event.num_reentry_fallback_adjacent_bin;
    else if (diag.fallbackLevel == "same_phase_any_bin")
      ++event.num_reentry_fallback_any_bin;
    else if (diag.fallbackLevel == "any_phase_same_bin")
      ++event.num_reentry_fallback_any_phase_same_bin;
    else if (diag.fallbackLevel == "any_portal")
      ++event.num_reentry_fallback_any_portal;
  }

  G4bool PushOpticalContinuation(
      G4Track *track,
      const G4ThreeVector &direction,
      const G4ThreeVector &polarization,
      G4double kineticEnergy,
      const G4ThreeVector &position)
  {
    auto *stackManager = G4EventManager::GetEventManager()->GetStackManager();
    if (track == nullptr || stackManager == nullptr || direction.mag2() <= 0.0)
      return false;

    auto *dynamicParticle = new G4DynamicParticle(
        G4OpticalPhoton::OpticalPhotonDefinition(),
        direction.unit(),
        kineticEnergy);
    dynamicParticle->SetPolarization(polarization);

    auto *continuationTrack = new G4Track(
        dynamicParticle,
        track->GetGlobalTime(),
        position);
    continuationTrack->SetTrackID(track->GetTrackID());
    continuationTrack->SetParentID(track->GetParentID());
    continuationTrack->SetLocalTime(track->GetLocalTime());
    continuationTrack->SetProperTime(track->GetProperTime());
    continuationTrack->SetWeight(track->GetWeight());
    stackManager->PushOneTrack(continuationTrack);
    return true;
  }

  G4int PhaseFunctionBin(const G4double cosTheta)
  {
    constexpr G4double kMin = -1.0;
    constexpr G4double kMax = 1.0;
    constexpr G4double kSpan = kMax - kMin;
    const G4double shifted = std::clamp(cosTheta, kMin, kMax) - kMin;
    const G4double normalized = shifted / kSpan;
    const G4int index = static_cast<G4int>(
        normalized * static_cast<G4double>(StageDPhotonEventRecord::kPhaseFunctionBins));
    return std::clamp(index, 0,
                      static_cast<G4int>(StageDPhotonEventRecord::kPhaseFunctionBins) - 1);
  }
}

StageDOpticalSteppingAction::StageDOpticalSteppingAction(
    StageDOpticalRunAction *runAction,
    StageDOpticalEventAction *eventAction,
    AnalysisConfig *config)
    : G4UserSteppingAction(),
      fConfig(config),
      fRunAction(runAction),
      fEventAction(eventAction),
      fDetector(nullptr),
      fBoundaryProcess(nullptr),
      fReentrySampler(nullptr)
{
}

StageDOpticalSteppingAction::~StageDOpticalSteppingAction()
{
  delete fReentrySampler;
}

void StageDOpticalSteppingAction::PrepareForNewRun()
{
  delete fReentrySampler;
  fReentrySampler = nullptr;
  fDetector = nullptr;
  fBoundaryProcess = nullptr;
}

const DetectorConstruction *StageDOpticalSteppingAction::ResolveDetector() const
{
  if (fDetector == nullptr)
  {
    fDetector = dynamic_cast<const DetectorConstruction *>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  }
  return fDetector;
}

G4OpBoundaryProcess *StageDOpticalSteppingAction::ResolveBoundaryProcess() const
{
  if (fBoundaryProcess != nullptr)
    return fBoundaryProcess;

  auto *definition = G4OpticalPhoton::OpticalPhotonDefinition();
  auto *processManager = definition != nullptr ? definition->GetProcessManager() : nullptr;
  auto *processList = processManager != nullptr ? processManager->GetProcessList() : nullptr;
  if (processList == nullptr)
    return nullptr;

  const G4int nProcesses = processManager->GetProcessListLength();
  for (G4int i = 0; i < nProcesses; ++i)
  {
    fBoundaryProcess = dynamic_cast<G4OpBoundaryProcess *>((*processList)[i]);
    if (fBoundaryProcess != nullptr)
      break;
  }
  return fBoundaryProcess;
}

G4bool StageDOpticalSteppingAction::HandleBoundaryReentry(
    const G4Step *step,
    G4Track *track,
    const DetectorConstruction *detector,
    DetectorConstruction::Phase prePhase)
{
  if (fConfig == nullptr || fEventAction == nullptr || detector == nullptr)
    return false;

  auto &event = fEventAction->MutableCurrentEvent();
  if (fConfig->stageD_boundary_mode == "escape")
  {
    fEventAction->MarkCensoredEncounterIfActive();
    fEventAction->SetFinalStatus("escaped_debug", false);
    track->SetTrackStatus(fStopAndKill);
    return true;
  }

  const G4bool isPeriodicWrap =
      fConfig->stageD_boundary_mode == "periodic_wrap";
  if (!isPeriodicWrap &&
      fConfig->stageD_boundary_mode != "same_phase_reentry")
    return false;

  if (event.num_reentry >= fConfig->stageD_max_reentry)
  {
    fEventAction->SetFinalStatus("max_reentry", false);
    track->SetTrackStatus(fStopAndKill);
    return true;
  }

  const G4ThreeVector prePos = step->GetPreStepPoint()->GetPosition();
  const G4ThreeVector oldDir = step->GetPreStepPoint()->GetMomentumDirection();
  const G4ThreeVector oldPolarization = step->GetPreStepPoint()->GetPolarization();
  const G4ThreeVector postPos = step->GetPostStepPoint()->GetPosition();
  if (oldDir.mag2() <= 0.0)
  {
    StageDReentrySampler::ReentryDiagnostics diag;
    diag.exitPhase = prePhase;
    diag.exitInsidePoint = postPos;
    diag.strategy = "invalid_exit_direction";
    diag.fallbackLevel = "unsupported_exit_direction";
    StageDReentrySampler::ReentryContext ctx;
    ctx.phase = prePhase;
    ctx.prePos = prePos;
    ctx.postPos = postPos;
    ctx.oldDir = oldDir;
    ctx.exitPoint = postPos;
    ctx.exitInsidePoint = postPos;
    ctx.wavelengthNm = event.wavelength_nm;
    ctx.eventID = event.photonID;
    ctx.reentryIndex = event.num_reentry + event.num_reentry_failed + 1;
    if (fRunAction != nullptr)
      fRunAction->RecordReentryDiagnostic(MakeReentryDiagnosticRecord(ctx, diag));
    ++event.num_reentry_failed;
    fEventAction->SetFinalStatus("reentry_failed", false);
    track->SetTrackStatus(fStopAndKill);
    return true;
  }

  G4ThreeVector exitPoint = postPos;
  ComputeExitPointOnRveBox(prePos, oldDir, detector, exitPoint);
  const G4ThreeVector exitInsidePoint = exitPoint - kBoundaryEpsilon * oldDir.unit();

  if (isPeriodicWrap)
  {
    StageDReentrySampler::ReentryContext ctx;
    ctx.phase = prePhase;
    ctx.prePos = prePos;
    ctx.postPos = postPos;
    ctx.oldDir = oldDir;
    ctx.exitPoint = exitPoint;
    ctx.exitInsidePoint = exitInsidePoint;
    ctx.wavelengthNm = event.wavelength_nm;
    ctx.eventID = event.photonID;
    ctx.reentryIndex = event.num_reentry + event.num_reentry_failed + 1;

    StageDReentrySampler::ReentryDiagnostics diag;
    diag.strategy = "periodic_wrap";
    diag.fallbackLevel = "none";
    diag.exitPhase = prePhase;
    diag.exitInsidePoint = exitInsidePoint;

    if (!detector->UsesPeriodicTransport())
    {
      diag.fallbackLevel = "periodic_geometry_required";
      if (fRunAction != nullptr)
        fRunAction->RecordReentryDiagnostic(MakeReentryDiagnosticRecord(ctx, diag));
      ++event.num_reentry_failed;
      fEventAction->SetFinalStatus("periodic_geometry_required", false);
      track->SetTrackStatus(fStopAndKill);
      return true;
    }

    const G4ThreeVector newPosition = detector->WrapToPrimaryCellInside(
        postPos, oldDir, kBoundaryEpsilon);
    const auto entryPhase = detector->FindPhaseAtPoint(newPosition);
    diag.entryPhase = entryPhase;
    diag.entryPoint = newPosition;

    if (!IsReentryPhase(prePhase) || entryPhase != prePhase)
    {
      diag.fallbackLevel = "phase_mismatch";
      if (fRunAction != nullptr)
        fRunAction->RecordReentryDiagnostic(MakeReentryDiagnosticRecord(ctx, diag));
      ++event.num_reentry_failed;
      fEventAction->SetFinalStatus("periodic_phase_mismatch", false);
      track->SetTrackStatus(fStopAndKill);
      return true;
    }

    if (!PushOpticalContinuation(
            track,
            oldDir,
            oldPolarization,
            step->GetPreStepPoint()->GetKineticEnergy(),
            newPosition))
    {
      diag.fallbackLevel = "stack_unavailable";
      if (fRunAction != nullptr)
        fRunAction->RecordReentryDiagnostic(MakeReentryDiagnosticRecord(ctx, diag));
      ++event.num_reentry_failed;
      fEventAction->SetFinalStatus("reentry_failed", false);
      track->SetTrackStatus(fStopAndKill);
      return true;
    }

    if (fRunAction != nullptr)
      fRunAction->RecordReentryDiagnostic(MakeReentryDiagnosticRecord(ctx, diag));
    ++event.num_reentry;
    if (prePhase == DetectorConstruction::Phase::BN)
      ++event.num_reentry_BN;
    else if (prePhase == DetectorConstruction::Phase::ZnS)
      ++event.num_reentry_ZnS;
    else if (prePhase == DetectorConstruction::Phase::Matrix)
      ++event.num_reentry_matrix;
    track->SetTrackStatus(fStopAndKill);
    return true;
  }

  if (fReentrySampler == nullptr)
  {
    fReentrySampler = new StageDReentrySampler(detector, fConfig);
    if (fRunAction != nullptr)
      fRunAction->SetReentryPortalSummary(fReentrySampler->GetPortalSummary());
  }

  auto phase = prePhase;
  if (!IsReentryPhase(phase))
  {
    phase = fReentrySampler->FastPhaseAtPointForReentry(exitInsidePoint);
    if (phase == DetectorConstruction::Phase::World ||
        phase == DetectorConstruction::Phase::Unknown)
    {
      phase = fReentrySampler->FastPhaseAtPointForReentry(prePos);
    }
  }

  if (event.encounter_active)
  {
    const auto encounterPhase = PhaseFromLabel(event.encounter_particle_phase);
    if (!IsParticlePhase(prePhase) || encounterPhase != prePhase)
    {
      ++event.num_inconsistent_encounter_state;
      fEventAction->MarkCensoredEncounterIfActive();
    }
  }
  else if (event.source_inside_particle_pending_exit && !IsParticlePhase(prePhase))
  {
    ++event.num_inconsistent_encounter_state;
    fEventAction->MarkCensoredEncounterIfActive();
  }

  StageDReentrySampler::ReentryContext ctx;
  ctx.phase = phase;
  ctx.prePos = prePos;
  ctx.postPos = postPos;
  ctx.oldDir = oldDir;
  ctx.exitPoint = exitPoint;
  ctx.exitInsidePoint = exitInsidePoint;
  ctx.wavelengthNm = event.wavelength_nm;
  ctx.eventID = event.photonID;
  ctx.reentryIndex = event.num_reentry + event.num_reentry_failed + 1;

  StageDReentrySampler::ReentryDiagnostics diag;
  diag.exitPhase = phase;
  diag.exitInsidePoint = exitInsidePoint;

  if (!IsReentryPhase(phase))
  {
    diag.strategy = "phase_resolve_failed";
    diag.fallbackLevel = "unsupported_exit_phase";
    if (fRunAction != nullptr)
      fRunAction->RecordReentryDiagnostic(MakeReentryDiagnosticRecord(ctx, diag));
    ++event.num_reentry_failed;
    fEventAction->SetFinalStatus("reentry_failed", false);
    track->SetTrackStatus(fStopAndKill);
    return true;
  }

  G4ThreeVector newPosition;
  const G4bool ok = fReentrySampler->SampleReentry(ctx, newPosition, diag);
  if (fRunAction != nullptr)
    fRunAction->RecordReentryDiagnostic(MakeReentryDiagnosticRecord(ctx, diag));

  if (!ok)
  {
    ++event.num_reentry_failed;
    fEventAction->SetFinalStatus("reentry_failed", false);
    track->SetTrackStatus(fStopAndKill);
    return true;
  }

  if (!PushOpticalContinuation(
          track,
          oldDir,
          oldPolarization,
          step->GetPreStepPoint()->GetKineticEnergy(),
          newPosition))
  {
    ++event.num_reentry_failed;
    fEventAction->SetFinalStatus("reentry_failed", false);
    track->SetTrackStatus(fStopAndKill);
    return true;
  }

  ++event.num_reentry;
  if (phase == DetectorConstruction::Phase::BN)
    ++event.num_reentry_BN;
  else if (phase == DetectorConstruction::Phase::ZnS)
    ++event.num_reentry_ZnS;
  else if (phase == DetectorConstruction::Phase::Matrix)
    ++event.num_reentry_matrix;
  AccumulateReentryCounters(event, diag);

  track->SetTrackStatus(fStopAndKill);
  return true;
}

G4bool StageDOpticalSteppingAction::HandleLimitKills(const G4Step *step, G4Track *track)
{
  if (fConfig == nullptr || fEventAction == nullptr)
    return false;

  auto &event = fEventAction->MutableCurrentEvent();
  const G4int primaryEncounterCount =
      UsesAngleThresholdEncounterMetric(fConfig)
          ? event.num_encounter_effective_total
          : event.num_encounter_total;

  if (fConfig->stageD_target_primary_scatter > 0 &&
      primaryEncounterCount >= fConfig->stageD_target_primary_scatter)
  {
    fEventAction->SetFinalStatus("target_primary_scatter", false);
    track->SetTrackStatus(fStopAndKill);
    return true;
  }

  if (HandleHardPathLimits(track))
    return true;

  if (track->GetTrackStatus() == fStopAndKill &&
      event.final_status == "in_progress")
  {
    fEventAction->SetFinalStatus("lost", false);
  }

  return false;
}

G4bool StageDOpticalSteppingAction::HandleHardPathLimits(G4Track *track)
{
  if (fConfig == nullptr || fEventAction == nullptr || track == nullptr)
    return false;

  auto &event = fEventAction->MutableCurrentEvent();
  if (event.num_steps >= fConfig->stageD_max_steps)
  {
    fEventAction->SetFinalStatus("max_steps", false);
    track->SetTrackStatus(fStopAndKill);
    return true;
  }

  if (event.total_path_length_um >= fConfig->stageD_max_path_length_um)
  {
    fEventAction->SetFinalStatus("max_path_length", false);
    track->SetTrackStatus(fStopAndKill);
    return true;
  }

  return false;
}

void StageDOpticalSteppingAction::UserSteppingAction(const G4Step *step)
{
  if (step == nullptr || fEventAction == nullptr)
    return;

  G4Track *track = step->GetTrack();
  if (!IsOpticalPhoton(track))
    return;

  auto &event = fEventAction->MutableCurrentEvent();
  ++event.num_steps;
  event.total_path_length_um += step->GetStepLength() / um;

  const auto *prePoint = step->GetPreStepPoint();
  const auto *postPoint = step->GetPostStepPoint();
  if (prePoint == nullptr || postPoint == nullptr)
    return;

  const auto *detector = ResolveDetector();
  RveFace outerFace = RveFace::None;
  const G4bool isOuterRveBoundary =
      (detector != nullptr && IsRveOuterBoundaryStep(step, detector, &outerFace));
  const auto *prePV = prePoint->GetPhysicalVolume();
  const auto *postPV = postPoint->GetPhysicalVolume();
  const auto prePhase = PhaseFromPhysicalVolume(prePV, detector);
  const auto postPhase =
      isOuterRveBoundary ? DetectorConstruction::Phase::World
                         : PhaseFromPhysicalVolume(postPV, detector);

  const G4double stepLengthUm = step->GetStepLength() / um;
  if (prePhase == DetectorConstruction::Phase::BN)
    event.path_length_bn_um += stepLengthUm;
  else if (prePhase == DetectorConstruction::Phase::ZnS)
    event.path_length_zns_um += stepLengthUm;
  else if (prePhase == DetectorConstruction::Phase::Matrix)
    event.path_length_matrix_um += stepLengthUm;
  else
    event.path_length_world_um += stepLengthUm;

  const std::string processName = ProcessName(postPoint);
  auto *boundaryProcess = ResolveBoundaryProcess();
  const G4OpBoundaryProcessStatus boundaryStatus =
      boundaryProcess != nullptr ? boundaryProcess->GetStatus() : Undefined;
  const OpticalBoundaryTransition boundaryTransition =
      ClassifyBoundaryTransition(boundaryStatus);
  const G4ThreeVector preDir = prePoint->GetMomentumDirection();
  const G4ThreeVector postDir = postPoint->GetMomentumDirection();
  const G4bool isGeomBoundary = (postPoint->GetStepStatus() == fGeomBoundary);

  if (isOuterRveBoundary)
  {
    if (HandleHardPathLimits(track))
      return;

    ++event.num_outer_boundary_hits;
    if (fConfig == nullptr || fConfig->stageD_boundary_mode != "periodic_wrap")
      AccumulateOuterBoundaryStatus(event, boundaryStatus);

    const G4int reentryBefore = event.num_reentry;
    if (HandleBoundaryReentry(step, track, detector, prePhase))
    {
      if (event.num_reentry > reentryBefore)
      {
        ++event.num_outer_boundary_reentry_success;
      }
      else if (event.final_status == "escaped_debug")
      {
        ++event.num_outer_boundary_escape;
      }
      else if (event.final_status == "max_reentry")
      {
        ++event.num_outer_boundary_max_reentry_stop;
      }
      else if (event.final_status == "periodic_phase_mismatch")
      {
        ++event.num_outer_boundary_periodic_phase_mismatch;
        ++event.num_outer_boundary_reentry_failed;
      }
      else
      {
        ++event.num_outer_boundary_reentry_algorithm_failed;
        ++event.num_outer_boundary_reentry_failed;
      }
      return;
    }

    ++event.num_outer_boundary_reentry_algorithm_failed;
    ++event.num_outer_boundary_reentry_failed;
    fEventAction->MarkCensoredEncounterIfActive();
    fEventAction->SetFinalStatus("reentry_failed", false);
    track->SetTrackStatus(fStopAndKill);
    return;
  }

  if (processName == "OpAbsorption" ||
      boundaryTransition == OpticalBoundaryTransition::Absorption)
  {
    if ((event.encounter_active ||
         event.source_inside_particle_pending_exit) &&
        !IsParticlePhase(prePhase))
    {
      ++event.num_inconsistent_encounter_state;
    }
    fEventAction->MarkCensoredEncounterIfActive();
    fEventAction->MarkAbsorbed(DetectorConstruction::PhaseName(prePhase));
    return;
  }

  if (boundaryTransition == OpticalBoundaryTransition::Detection)
  {
    fEventAction->MarkCensoredEncounterIfActive();
    fEventAction->SetFinalStatus("detected", false);
    return;
  }

  const G4bool phaseTransmissionFallback =
      boundaryTransition == OpticalBoundaryTransition::Unavailable &&
      IsPhaseTransmission(prePhase, postPhase);
  const G4bool isBoundaryTransmission =
      boundaryTransition == OpticalBoundaryTransition::Transmission ||
      phaseTransmissionFallback;
  const G4bool isBoundaryReflection =
      boundaryTransition == OpticalBoundaryTransition::Reflection;

  const G4bool isMaterialBoundary =
      boundaryTransition != OpticalBoundaryTransition::StepTooSmall &&
      (isGeomBoundary ||
       processName == "OpBoundary" ||
       phaseTransmissionFallback);
  if (isMaterialBoundary)
  {
    ++event.num_material_boundary;
  }

  if (prePhase == DetectorConstruction::Phase::Matrix &&
      isMaterialBoundary &&
      isBoundaryReflection)
  {
    DetectorConstruction::Phase particlePhase = postPhase;
    if (!IsParticlePhase(particlePhase) && detector != nullptr)
      particlePhase = ProbeForwardPhase(step, detector, preDir);
    if (IsParticlePhase(particlePhase))
    {
      RecordCompleteEncounter(event, ClampCosTheta(preDir, postDir), particlePhase, fConfig);
      ++event.num_surface_reflection_encounter;
    }
    else
    {
      ++event.num_unknown_particle_reflection;
    }
  }
  else if (prePhase == DetectorConstruction::Phase::Matrix &&
           IsParticlePhase(postPhase) &&
           isBoundaryTransmission)
  {
    if (event.encounter_active || event.source_inside_particle_pending_exit)
    {
      ++event.num_inconsistent_encounter_state;
      fEventAction->MarkCensoredEncounterIfActive();
    }
    event.encounter_active = true;
    event.encounter_has_matrix_entry = true;
    event.encounter_particle_phase = DetectorConstruction::PhaseName(postPhase);
    event.encounter_matrix_entry_direction = preDir;
  }
  else if (IsParticlePhase(prePhase) &&
           postPhase == DetectorConstruction::Phase::Matrix &&
           isBoundaryTransmission)
  {
    const auto encounterPhase = PhaseFromLabel(event.encounter_particle_phase);
    if (event.source_inside_particle_pending_exit &&
        (!event.encounter_active || !event.encounter_has_matrix_entry))
    {
      ++event.num_incomplete_initial_particle_exit;
      event.source_inside_particle_pending_exit = false;
      event.encounter_particle_phase.clear();
    }
    else if (event.encounter_active &&
             event.encounter_has_matrix_entry &&
             encounterPhase == prePhase)
    {
      RecordCompleteEncounter(
          event,
          ClampCosTheta(event.encounter_matrix_entry_direction, postDir),
          prePhase,
          fConfig);
      event.encounter_active = false;
      event.encounter_has_matrix_entry = false;
      event.encounter_particle_phase.clear();
    }
    else
    {
      ++event.num_inconsistent_encounter_state;
      fEventAction->MarkCensoredEncounterIfActive();
    }
  }
  else if (event.encounter_active &&
           IsParticlePhase(prePhase) &&
           prePhase == PhaseFromLabel(event.encounter_particle_phase) &&
           isBoundaryReflection)
  {
    // Internal particle reflections do not complete a matrix-to-matrix encounter.
  }
  else if (IsParticlePhase(prePhase) &&
           IsParticlePhase(postPhase) &&
           prePhase != postPhase &&
           isBoundaryTransmission)
  {
    ++event.num_particle_to_particle_boundary;
  }

  if (processName != "OpAbsorption")
  {
    const G4double thetaDeg =
        AngleDeg(preDir, postDir);
    const G4double thetaThresholdDeg =
        (fConfig != nullptr) ? fConfig->stageD_theta_threshold_deg : 0.0;
    if (thetaDeg > thetaThresholdDeg)
    {
      ++event.num_real_scatter;
      const G4double cosTheta = std::clamp(
          preDir.dot(postDir),
          -1.0, 1.0);
      event.sum_cos_theta += cosTheta;
      if (isMaterialBoundary)
      {
        ++event.num_boundary_scatter;
        event.sum_cos_theta_boundary += cosTheta;
        if (prePhase == DetectorConstruction::Phase::BN)
        {
          ++event.num_boundary_scatter_BN;
          event.sum_cos_theta_boundary_BN += cosTheta;
        }
        else if (prePhase == DetectorConstruction::Phase::ZnS)
        {
          ++event.num_boundary_scatter_ZnS;
          event.sum_cos_theta_boundary_ZnS += cosTheta;
        }
      }
      else
      {
        ++event.num_bulk_scatter;
        event.sum_cos_theta_bulk += cosTheta;
      }
    }
  }

  HandleLimitKills(step, track);
}
