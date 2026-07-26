#include "StageDOpticalEventAction.hh"

#include "AnalysisConfig.hh"
#include "DetectorConstruction.hh"
#include "StageDOpticalPrimaryGeneratorAction.hh"
#include "StageDOpticalRunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

#include <algorithm>
#include <filesystem>
#include <cmath>
#include <limits>
#include <sstream>

namespace
{
  std::string RatioTag(G4double bnWt, G4double znsWt)
  {
    auto toTag = [](G4double value)
    {
      const G4double rounded = std::round(value);
      std::ostringstream oss;
      if (std::abs(value - rounded) < 1.0e-9)
        oss << static_cast<long long>(rounded);
      else
        oss << value;
      return oss.str();
    };
    return toTag(bnWt) + "-" + toTag(znsWt);
  }
}

StageDOpticalEventAction::StageDOpticalEventAction(
    StageDOpticalRunAction *runAction,
    const StageDOpticalPrimaryGeneratorAction *primaryAction,
    AnalysisConfig *config)
    : G4UserEventAction(),
      fConfig(config),
      fRunAction(runAction),
      fPrimaryAction(primaryAction),
      fCurrentEvent()
{
}

StageDOpticalEventAction::~StageDOpticalEventAction() = default;

void StageDOpticalEventAction::BeginOfEventAction(const G4Event *event)
{
  fCurrentEvent = StageDPhotonEventRecord{};

  if (event == nullptr || fConfig == nullptr || fPrimaryAction == nullptr)
    return;

  const auto *detector = dynamic_cast<const DetectorConstruction *>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  const StageDPhotonLaunchRecord &launch = fPrimaryAction->GetCurrentPhotonRecord();

  fCurrentEvent.photonID = event->GetEventID();
  fCurrentEvent.ratio = RatioTag(fConfig->bnWt, fConfig->znsWt);
  fCurrentEvent.placement_file =
      detector ? detector->GetLoadedPlacementFileForRecord() : "unknown";
  fCurrentEvent.source_mode = fConfig->stageD_source_mode;
  fCurrentEvent.boundary_mode = fConfig->stageD_boundary_mode;
  fCurrentEvent.reentry_mode = fConfig->stageD_reentry_mode;
  fCurrentEvent.particle_reentry_mode = fConfig->stageD_particle_reentry_mode;
  fCurrentEvent.matrix_reentry_mode = fConfig->stageD_matrix_reentry_mode;
  fCurrentEvent.wavelength_nm = launch.wavelength_nm;
  fCurrentEvent.source_phase = launch.source_phase;
  fCurrentEvent.source_x_um = launch.source_position.x() / um;
  fCurrentEvent.source_y_um = launch.source_position.y() / um;
  fCurrentEvent.source_z_um = launch.source_position.z() / um;
  fCurrentEvent.final_status = "in_progress";
  fCurrentEvent.weight = launch.photon_weight;

  if (launch.source_phase == "BN" || launch.source_phase == "ZnS")
  {
    fCurrentEvent.source_inside_particle_pending_exit = true;
    fCurrentEvent.encounter_particle_phase = launch.source_phase;
  }
}

void StageDOpticalEventAction::EndOfEventAction(const G4Event *event)
{
  (void)event;

  const G4bool useThresholdedEncounterMetric =
      (fConfig != nullptr &&
       fConfig->stageD_scatter_metric == "particle_encounter_angle_threshold");
  const G4double nan = std::numeric_limits<G4double>::quiet_NaN();

  if (fCurrentEvent.final_status == "in_progress" ||
      fCurrentEvent.final_status == "continued_reentry")
    fCurrentEvent.final_status = "lost";

  if ((fCurrentEvent.encounter_active ||
       fCurrentEvent.source_inside_particle_pending_exit) &&
      !fCurrentEvent.absorbed)
  {
    MarkCensoredEncounterIfActive("event_end");
  }

  if (fCurrentEvent.num_real_scatter > 0)
  {
    fCurrentEvent.mean_cos_theta_for_this_photon =
        fCurrentEvent.sum_cos_theta /
        static_cast<G4double>(fCurrentEvent.num_real_scatter);
  }
  else
  {
    fCurrentEvent.mean_cos_theta_for_this_photon = 0.0;
  }

  if (fCurrentEvent.num_particle_scatter > 0)
  {
    fCurrentEvent.mean_cos_theta_particle_for_this_photon =
        fCurrentEvent.sum_cos_theta_particle /
        static_cast<G4double>(fCurrentEvent.num_particle_scatter);
  }
  else
  {
    fCurrentEvent.mean_cos_theta_particle_for_this_photon = 0.0;
  }

  if (fCurrentEvent.num_bulk_scatter > 0)
  {
    fCurrentEvent.mean_cos_theta_bulk_for_this_photon =
        fCurrentEvent.sum_cos_theta_bulk /
        static_cast<G4double>(fCurrentEvent.num_bulk_scatter);
  }
  else
  {
    fCurrentEvent.mean_cos_theta_bulk_for_this_photon = 0.0;
  }

  if (fCurrentEvent.num_boundary_scatter > 0)
  {
    fCurrentEvent.mean_cos_theta_boundary_for_this_photon =
        fCurrentEvent.sum_cos_theta_boundary /
        static_cast<G4double>(fCurrentEvent.num_boundary_scatter);
  }
  else
  {
    fCurrentEvent.mean_cos_theta_boundary_for_this_photon = 0.0;
  }

  if (fCurrentEvent.num_encounter_total > 0)
  {
    const G4double nEncounterRaw =
        static_cast<G4double>(fCurrentEvent.num_encounter_total);
    const G4double g1Raw =
        fCurrentEvent.sum_cos_theta_encounter / nEncounterRaw;
    const G4double meanCos2Raw =
        fCurrentEvent.sum_cos2_theta_encounter / nEncounterRaw;
    fCurrentEvent.g1_encounter_raw_for_this_photon = g1Raw;
    fCurrentEvent.g2_encounter_raw_for_this_photon =
        0.5 * (3.0 * meanCos2Raw - 1.0);
  }
  else
  {
    fCurrentEvent.g1_encounter_raw_for_this_photon = nan;
    fCurrentEvent.g2_encounter_raw_for_this_photon = nan;
  }

  if (useThresholdedEncounterMetric &&
      fCurrentEvent.num_encounter_effective_total > 0)
  {
    const G4double nEncounter =
        static_cast<G4double>(fCurrentEvent.num_encounter_effective_total);
    const G4double g1 =
        fCurrentEvent.sum_cos_theta_encounter_effective / nEncounter;
    const G4double meanCos2 =
        fCurrentEvent.sum_cos2_theta_encounter_effective / nEncounter;
    fCurrentEvent.g1_encounter_for_this_photon = g1;
    fCurrentEvent.g2_encounter_for_this_photon =
        0.5 * (3.0 * meanCos2 - 1.0);
  }
  else if (!useThresholdedEncounterMetric &&
           fCurrentEvent.num_encounter_total > 0)
  {
    fCurrentEvent.g1_encounter_for_this_photon =
        fCurrentEvent.g1_encounter_raw_for_this_photon;
    fCurrentEvent.g2_encounter_for_this_photon =
        fCurrentEvent.g2_encounter_raw_for_this_photon;
  }
  else
  {
    fCurrentEvent.g1_encounter_for_this_photon = nan;
    fCurrentEvent.g2_encounter_for_this_photon = nan;
  }

  const G4double mediumPathLengthUm =
      fCurrentEvent.path_length_bn_um +
      fCurrentEvent.path_length_zns_um +
      fCurrentEvent.path_length_matrix_um;
  fCurrentEvent.mu_s_prime_direct_encounter_raw_per_um_for_this_photon =
      (mediumPathLengthUm > 0.0)
          ? (fCurrentEvent.sum_one_minus_cos_theta_encounter / mediumPathLengthUm)
          : 0.0;
  fCurrentEvent.mu_s_prime_direct_encounter_per_um_for_this_photon =
      (mediumPathLengthUm > 0.0)
          ? ((useThresholdedEncounterMetric
                  ? fCurrentEvent.sum_one_minus_cos_theta_encounter_effective
                  : fCurrentEvent.sum_one_minus_cos_theta_encounter) /
             mediumPathLengthUm)
          : 0.0;

  if (fRunAction != nullptr)
    fRunAction->RecordPhotonEvent(fCurrentEvent);
}

void StageDOpticalEventAction::SetFinalStatus(
    const std::string &status,
    G4bool absorbed)
{
  if (fCurrentEvent.encounter_active ||
      fCurrentEvent.source_inside_particle_pending_exit)
  {
    if (status == "max_steps")
      MarkCensoredEncounterIfActive("max_steps");
    else if (status == "max_path_length")
      MarkCensoredEncounterIfActive("max_path");
    else if (status == "max_reentry")
      MarkCensoredEncounterIfActive("max_reentry");
    else if (status == "reentry_failed" ||
             status == "periodic_phase_mismatch" ||
             status == "periodic_geometry_required")
      MarkCensoredEncounterIfActive("reentry_failure");
    else if (status == "detected")
      MarkCensoredEncounterIfActive("detection");
    else if (status == "target_primary_scatter")
      MarkCensoredEncounterIfActive("target_scatter");
    else if (status == "escaped_debug")
      MarkCensoredEncounterIfActive("escape");
  }

  fCurrentEvent.final_status = status;
  fCurrentEvent.absorbed = absorbed;
}

void StageDOpticalEventAction::MarkAbsorbed(const std::string &phaseLabel)
{
  ++fCurrentEvent.num_absorbed_total;
  if (phaseLabel == "BN")
    ++fCurrentEvent.num_absorbed_BN;
  else if (phaseLabel == "ZnS")
    ++fCurrentEvent.num_absorbed_ZnS;
  else if (phaseLabel == "Matrix")
    ++fCurrentEvent.num_absorbed_Matrix;
  else
    ++fCurrentEvent.num_absorbed_World;

  SetFinalStatus("absorbed", true);
}

void StageDOpticalEventAction::MarkCensoredEncounterIfActive(
    const std::string &reason)
{
  if (!fCurrentEvent.encounter_active &&
      !fCurrentEvent.source_inside_particle_pending_exit)
  {
    return;
  }

  ++fCurrentEvent.num_censored_particle_encounter;
  if (reason == "absorption")
    ++fCurrentEvent.num_censored_by_absorption;
  else if (reason == "max_steps")
    ++fCurrentEvent.num_censored_by_max_steps;
  else if (reason == "max_path")
    ++fCurrentEvent.num_censored_by_max_path;
  else if (reason == "max_reentry")
    ++fCurrentEvent.num_censored_by_max_reentry;
  else if (reason == "reentry_failure")
    ++fCurrentEvent.num_censored_by_reentry_failure;
  else if (reason == "detection")
    ++fCurrentEvent.num_censored_by_detection;
  else if (reason == "target_scatter")
    ++fCurrentEvent.num_censored_by_target_scatter;
  else if (reason == "escape")
    ++fCurrentEvent.num_censored_by_escape;
  else if (reason == "event_end")
    ++fCurrentEvent.num_censored_by_event_end;
  else if (reason == "state_inconsistency")
    ++fCurrentEvent.num_censored_by_state_inconsistency;
  fCurrentEvent.encounter_active = false;
  fCurrentEvent.encounter_has_matrix_entry = false;
  fCurrentEvent.source_inside_particle_pending_exit = false;
  fCurrentEvent.encounter_particle_phase.clear();
}
