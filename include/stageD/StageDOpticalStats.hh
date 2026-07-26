#ifndef StageDOpticalStats_h
#define StageDOpticalStats_h 1

#include "G4ThreeVector.hh"
#include "globals.hh"

#include <array>
#include <string>

struct StageDPhotonLaunchRecord
{
  G4int geantEventID = -1;
  G4double wavelength_nm = 0.0;
  std::string source_mode;
  std::string source_phase;
  G4ThreeVector source_position;
  G4ThreeVector momentum_direction;
  G4ThreeVector polarization;
  G4double photon_weight = 1.0;
  G4bool is_continuation = false;
};

struct StageDReentryDiagnosticRecord
{
  G4int event_id = -1;
  G4int reentry_index = 0;
  std::string strategy;
  std::string fallback_level;
  std::string exit_phase;
  std::string entry_phase;

  G4ThreeVector old_dir;
  G4ThreeVector exit_point;
  G4ThreeVector entry_point;

  G4double particle_q_exit = -1.0;
  G4double particle_q_entry = -1.0;
  G4double particle_mu_exit = 0.0;
  G4double particle_mu_entry = 0.0;

  G4double matrix_clearance_exit_um = -1.0;
  G4double matrix_clearance_entry_um = -1.0;
  std::string matrix_nearest_phase_exit;
  std::string matrix_nearest_phase_entry;
  G4int matrix_clearance_bin_exit = -1;
  G4int matrix_clearance_bin_entry = -1;

  G4int trials = 0;
  G4double periodic_entry_offset_um = 0.0;
  G4int periodic_entry_search_trials = 0;
};

struct StageDReentryPortalSummary
{
  static constexpr std::size_t kFaceCount = 6;
  static constexpr std::size_t kBinCount = 4;

  G4int total_portal_count = 0;
  std::array<G4int, kFaceCount> portal_count_by_face{};
  std::array<G4int, kBinCount> portal_count_by_bin{};
};

struct StageDPhotonEventRecord
{
  static constexpr std::size_t kPhaseFunctionBins = 64;

  G4int photonID = -1;
  std::string ratio;
  std::string placement_file;
  std::string source_mode;
  std::string boundary_mode;
  std::string reentry_mode;
  std::string particle_reentry_mode;
  std::string matrix_reentry_mode;
  G4double wavelength_nm = 0.0;

  std::string source_phase;
  G4double source_x_um = 0.0;
  G4double source_y_um = 0.0;
  G4double source_z_um = 0.0;

  std::string final_status;
  G4bool absorbed = false;

  G4double total_path_length_um = 0.0;
  G4double path_length_bn_um = 0.0;
  G4double path_length_zns_um = 0.0;
  G4double path_length_matrix_um = 0.0;
  G4double path_length_world_um = 0.0;
  G4int num_steps = 0;

  G4int num_absorbed_total = 0;
  G4int num_absorbed_BN = 0;
  G4int num_absorbed_ZnS = 0;
  G4int num_absorbed_Matrix = 0;
  G4int num_absorbed_World = 0;

  G4int num_encounter_total = 0;
  G4int num_encounter_BN = 0;
  G4int num_encounter_ZnS = 0;
  G4int num_encounter_effective_total = 0;
  G4int num_encounter_effective_BN = 0;
  G4int num_encounter_effective_ZnS = 0;
  G4double sum_cos_theta_encounter = 0.0;
  G4double sum_cos_theta_encounter_BN = 0.0;
  G4double sum_cos_theta_encounter_ZnS = 0.0;
  G4double sum_cos_theta_encounter_effective = 0.0;
  G4double sum_cos_theta_encounter_effective_BN = 0.0;
  G4double sum_cos_theta_encounter_effective_ZnS = 0.0;
  G4double sum_one_minus_cos_theta_encounter = 0.0;
  G4double sum_one_minus_cos_theta_encounter_BN = 0.0;
  G4double sum_one_minus_cos_theta_encounter_ZnS = 0.0;
  G4double sum_one_minus_cos_theta_encounter_effective = 0.0;
  G4double sum_one_minus_cos_theta_encounter_effective_BN = 0.0;
  G4double sum_one_minus_cos_theta_encounter_effective_ZnS = 0.0;
  G4double sum_cos2_theta_encounter = 0.0;
  G4double sum_cos2_theta_encounter_BN = 0.0;
  G4double sum_cos2_theta_encounter_ZnS = 0.0;
  G4double sum_cos2_theta_encounter_effective = 0.0;
  G4double sum_cos2_theta_encounter_effective_BN = 0.0;
  G4double sum_cos2_theta_encounter_effective_ZnS = 0.0;

  G4int num_particle_scatter = 0;
  G4int num_particle_scatter_BN = 0;
  G4int num_particle_scatter_ZnS = 0;
  G4int num_real_scatter = 0;
  G4int num_bulk_scatter = 0;
  G4int num_boundary_scatter = 0;
  G4int num_boundary_scatter_BN = 0;
  G4int num_boundary_scatter_ZnS = 0;
  G4int num_material_boundary = 0;

  G4int num_complete_encounter_total = 0;
  G4int num_complete_encounter_BN = 0;
  G4int num_complete_encounter_ZnS = 0;
  G4int num_surface_reflection_encounter = 0;
  G4int num_incomplete_initial_particle_exit = 0;
  G4int num_censored_particle_encounter = 0;
  G4int num_inconsistent_encounter_state = 0;
  G4int num_particle_to_particle_boundary = 0;
  G4int num_unknown_particle_reflection = 0;

  G4int num_outer_boundary_hits = 0;
  G4int num_outer_boundary_reentry_success = 0;
  G4int num_outer_boundary_reentry_failed = 0;
  G4int num_outer_boundary_escape = 0;
  G4int num_outer_boundary_max_reentry_stop = 0;
  G4int num_outer_boundary_reentry_algorithm_failed = 0;
  G4int num_outer_boundary_periodic_phase_mismatch = 0;
  G4int num_outer_boundary_fresnel_reflection = 0;
  G4int num_outer_boundary_total_internal_reflection = 0;
  G4int num_outer_boundary_refraction = 0;
  G4int num_outer_boundary_transmission = 0;
  G4int num_outer_boundary_other_status = 0;

  G4int num_reentry = 0;
  G4int num_reentry_BN = 0;
  G4int num_reentry_ZnS = 0;
  G4int num_reentry_matrix = 0;
  G4int num_reentry_particle_q_mu = 0;
  G4int num_reentry_particle_q_only = 0;
  G4int num_reentry_particle_q_only_fallback = 0;
  G4int num_reentry_particle_volume_random = 0;
  G4int num_reentry_matrix_clearance_portal = 0;
  G4int num_reentry_fallback_same_bin = 0;
  G4int num_reentry_fallback_adjacent_bin = 0;
  G4int num_reentry_fallback_any_bin = 0;
  G4int num_reentry_fallback_any_phase_same_bin = 0;
  G4int num_reentry_fallback_any_portal = 0;
  G4int num_reentry_random_matrix_debug = 0;
  G4int num_reentry_failed = 0;

  G4double sum_cos_theta = 0.0;
  G4double sum_cos_theta_particle = 0.0;
  G4double sum_cos_theta_particle_BN = 0.0;
  G4double sum_cos_theta_particle_ZnS = 0.0;
  G4double sum_cos_theta_bulk = 0.0;
  G4double sum_cos_theta_boundary = 0.0;
  G4double sum_cos_theta_boundary_BN = 0.0;
  G4double sum_cos_theta_boundary_ZnS = 0.0;
  G4double g1_encounter_for_this_photon = 0.0;
  G4double g2_encounter_for_this_photon = 0.0;
  G4double mu_s_prime_direct_encounter_per_um_for_this_photon = 0.0;
  G4double g1_encounter_raw_for_this_photon = 0.0;
  G4double g2_encounter_raw_for_this_photon = 0.0;
  G4double mu_s_prime_direct_encounter_raw_per_um_for_this_photon = 0.0;
  G4double mean_cos_theta_particle_for_this_photon = 0.0;
  G4double mean_cos_theta_for_this_photon = 0.0;
  G4double mean_cos_theta_bulk_for_this_photon = 0.0;
  G4double mean_cos_theta_boundary_for_this_photon = 0.0;
  G4double weight = 1.0;

  std::array<G4int, kPhaseFunctionBins> phase_function_histogram_raw{};
  std::array<G4int, kPhaseFunctionBins> phase_function_histogram_thresholded{};

  G4bool first_complete_encounter_seen = false;
  G4double path_before_first_complete_encounter_um = 0.0;
  G4double post_first_encounter_path_bn_um = 0.0;
  G4double post_first_encounter_path_zns_um = 0.0;
  G4double post_first_encounter_path_matrix_um = 0.0;
  G4int post_first_encounter_count_raw = 0;
  G4double post_first_sum_cos_theta_raw = 0.0;
  G4double post_first_sum_one_minus_cos_theta_raw = 0.0;

  G4int num_censored_by_absorption = 0;
  G4int num_censored_by_max_steps = 0;
  G4int num_censored_by_max_path = 0;
  G4int num_censored_by_max_reentry = 0;
  G4int num_censored_by_reentry_failure = 0;
  G4int num_censored_by_event_end = 0;
  G4int num_censored_by_detection = 0;
  G4int num_censored_by_target_scatter = 0;
  G4int num_censored_by_escape = 0;
  G4int num_censored_by_state_inconsistency = 0;

  G4int num_periodic_entry_search_attempts = 0;
  G4int num_periodic_entry_search_success = 0;
  G4int num_periodic_entry_first_try_success = 0;
  G4int sum_periodic_entry_search_trials = 0;
  G4int max_periodic_entry_search_trials = 0;
  G4double max_periodic_entry_offset_um = 0.0;

  G4bool encounter_active = false;
  G4bool encounter_has_matrix_entry = false;
  G4bool source_inside_particle_pending_exit = false;
  std::string encounter_particle_phase;
  G4ThreeVector encounter_matrix_entry_direction;
};

struct StageDSourcePhaseAccumulator
{
  G4long nPhotons = 0;
  G4long nAbsorbed = 0;
  G4double totalMediumPathUm = 0.0;
  G4long totalEncounter = 0;
  G4double sumCosThetaEncounter = 0.0;
  G4double sumOneMinusCosThetaEncounter = 0.0;
  G4long nWithFirstCompleteEncounter = 0;
  G4long nPostFirstAbsorbed = 0;
  G4double postFirstPathBNUm = 0.0;
  G4double postFirstPathZnSUm = 0.0;
  G4double postFirstPathMatrixUm = 0.0;
  G4long postFirstEncounter = 0;
  G4double postFirstSumCosThetaEncounter = 0.0;
  G4double postFirstSumOneMinusCosThetaEncounter = 0.0;
};

struct StageDRunAccumulator
{
  G4long nPhotons = 0;
  G4long nAbsorbed = 0;
  G4long nAbsorbedBN = 0;
  G4long nAbsorbedZnS = 0;
  G4long nAbsorbedMatrix = 0;
  G4long nLost = 0;
  G4long nReentryFailed = 0;
  G4long nEscapedDebug = 0;
  G4long nPeriodicPhaseMismatch = 0;
  G4long nPeriodicGeometryRequired = 0;
  G4long nDetected = 0;
  G4long nMaxReentry = 0;
  G4long nMaxSteps = 0;
  G4long nMaxPathLength = 0;
  G4long nTargetPrimaryScatter = 0;
  G4double totalPathLengthUm = 0.0;
  G4double totalPathLengthBNUm = 0.0;
  G4double totalPathLengthZnSUm = 0.0;
  G4double totalPathLengthMatrixUm = 0.0;
  G4double totalPathLengthWorldUm = 0.0;
  G4long totalEncounter = 0;
  G4long totalEncounterBN = 0;
  G4long totalEncounterZnS = 0;
  G4long totalEncounterEffective = 0;
  G4long totalEncounterEffectiveBN = 0;
  G4long totalEncounterEffectiveZnS = 0;
  G4double sumCosThetaEncounter = 0.0;
  G4double sumCosThetaEncounterBN = 0.0;
  G4double sumCosThetaEncounterZnS = 0.0;
  G4double sumCosThetaEncounterEffective = 0.0;
  G4double sumCosThetaEncounterEffectiveBN = 0.0;
  G4double sumCosThetaEncounterEffectiveZnS = 0.0;
  G4double sumOneMinusCosThetaEncounter = 0.0;
  G4double sumOneMinusCosThetaEncounterBN = 0.0;
  G4double sumOneMinusCosThetaEncounterZnS = 0.0;
  G4double sumOneMinusCosThetaEncounterEffective = 0.0;
  G4double sumOneMinusCosThetaEncounterEffectiveBN = 0.0;
  G4double sumOneMinusCosThetaEncounterEffectiveZnS = 0.0;
  G4double sumCos2ThetaEncounter = 0.0;
  G4double sumCos2ThetaEncounterBN = 0.0;
  G4double sumCos2ThetaEncounterZnS = 0.0;
  G4double sumCos2ThetaEncounterEffective = 0.0;
  G4double sumCos2ThetaEncounterEffectiveBN = 0.0;
  G4double sumCos2ThetaEncounterEffectiveZnS = 0.0;
  std::array<G4long, StageDPhotonEventRecord::kPhaseFunctionBins> phaseFunctionCountsRaw{};
  std::array<G4long, StageDPhotonEventRecord::kPhaseFunctionBins> phaseFunctionCountsThresholded{};
  G4long nPhotonsWithFirstCompleteEncounter = 0;
  G4double totalPathBeforeFirstCompleteEncounterUm = 0.0;
  G4double totalPathBeforeFirstCompleteEncounterSeenUm = 0.0;
  G4double totalPostFirstEncounterPathBNUm = 0.0;
  G4double totalPostFirstEncounterPathZnSUm = 0.0;
  G4double totalPostFirstEncounterPathMatrixUm = 0.0;
  G4long totalPostFirstEncounter = 0;
  G4double sumPostFirstCosThetaEncounter = 0.0;
  G4double sumPostFirstOneMinusCosThetaEncounter = 0.0;
  G4long totalPostFirstAbsorbed = 0;
  G4long totalCensoredByAbsorption = 0;
  G4long totalCensoredByMaxSteps = 0;
  G4long totalCensoredByMaxPath = 0;
  G4long totalCensoredByMaxReentry = 0;
  G4long totalCensoredByReentryFailure = 0;
  G4long totalCensoredByEventEnd = 0;
  G4long totalCensoredByDetection = 0;
  G4long totalCensoredByTargetScatter = 0;
  G4long totalCensoredByEscape = 0;
  G4long totalCensoredByStateInconsistency = 0;
  std::array<StageDSourcePhaseAccumulator, 4> sourcePhase{};
  G4long totalPeriodicEntrySearchAttempts = 0;
  G4long totalPeriodicEntrySearchSuccess = 0;
  G4long totalPeriodicEntryFirstTrySuccess = 0;
  G4long totalPeriodicEntrySearchTrials = 0;
  G4int maxPeriodicEntrySearchTrials = 0;
  G4double maxPeriodicEntryOffsetUm = 0.0;
  G4long totalParticleScatter = 0;
  G4long totalParticleScatterBN = 0;
  G4long totalParticleScatterZnS = 0;
  G4long totalRealScatter = 0;
  G4long totalBulkScatter = 0;
  G4long totalBoundaryScatter = 0;
  G4long totalBoundaryScatterBN = 0;
  G4long totalBoundaryScatterZnS = 0;
  G4long totalMaterialBoundary = 0;
  G4long totalCompleteEncounter = 0;
  G4long totalCompleteEncounterBN = 0;
  G4long totalCompleteEncounterZnS = 0;
  G4long totalSurfaceReflectionEncounter = 0;
  G4long totalIncompleteInitialParticleExit = 0;
  G4long totalCensoredParticleEncounter = 0;
  G4long totalInconsistentEncounterState = 0;
  G4long totalParticleToParticleBoundary = 0;
  G4long totalUnknownParticleReflection = 0;
  G4long totalOuterBoundaryHits = 0;
  G4long totalOuterBoundaryReentrySuccess = 0;
  G4long totalOuterBoundaryReentryFailed = 0;
  G4long totalOuterBoundaryEscape = 0;
  G4long totalOuterBoundaryMaxReentryStop = 0;
  G4long totalOuterBoundaryReentryAlgorithmFailed = 0;
  G4long totalOuterBoundaryPeriodicPhaseMismatch = 0;
  G4long totalOuterBoundaryFresnelReflection = 0;
  G4long totalOuterBoundaryTotalInternalReflection = 0;
  G4long totalOuterBoundaryRefraction = 0;
  G4long totalOuterBoundaryTransmission = 0;
  G4long totalOuterBoundaryOtherStatus = 0;
  G4long totalReentry = 0;
  G4long totalReentryBN = 0;
  G4long totalReentryZnS = 0;
  G4long totalReentryMatrix = 0;
  G4long totalReentryParticleQMu = 0;
  G4long totalReentryParticleQOnly = 0;
  G4long totalReentryParticleQOnlyFallback = 0;
  G4long totalReentryParticleVolumeRandom = 0;
  G4long totalReentryMatrixClearancePortal = 0;
  G4long totalReentryFallbackSameBin = 0;
  G4long totalReentryFallbackAdjacentBin = 0;
  G4long totalReentryFallbackAnyBin = 0;
  G4long totalReentryFallbackAnyPhaseSameBin = 0;
  G4long totalReentryFallbackAnyPortal = 0;
  G4long totalReentryRandomMatrixDebug = 0;
  G4long totalReentryFailed = 0;
  G4double sumCosThetaParticleScatter = 0.0;
  G4double sumCosThetaParticleScatterBN = 0.0;
  G4double sumCosThetaParticleScatterZnS = 0.0;
  G4double sumCosThetaAllScatter = 0.0;
  G4double sumCosThetaBulkScatter = 0.0;
  G4double sumCosThetaBoundaryScatter = 0.0;
  G4double sumCosThetaBoundaryScatterBN = 0.0;
  G4double sumCosThetaBoundaryScatterZnS = 0.0;
};

#endif
