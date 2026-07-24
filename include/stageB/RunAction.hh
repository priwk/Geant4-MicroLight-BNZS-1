#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <cstdint>
#include <fstream>
#include <map>
#include <string>

class G4Run;
class PrimaryGeneratorAction;
class AnalysisConfig;

class RunAction : public G4UserRunAction
{
public:
  explicit RunAction(PrimaryGeneratorAction *primaryAction, AnalysisConfig *config = nullptr);
  ~RunAction() override;

  void BeginOfRunAction(const G4Run *) override;
  void EndOfRunAction(const G4Run *) override;

  void SetPrimaryAction(const PrimaryGeneratorAction *primaryAction);

  // For multi-input streaming:
  // switch output CSV according to the current input CSV path.
  void SwitchOutputCsvForInputPath(const std::string &inputPath);
  void SwitchOutputCsvForThicknessTag(const std::string &thicknessTag);

  enum class OutputMode
  {
    Full,
    Slim,
  };

  enum class BoundaryExitClass
  {
    PhysicalSurfaceExit,
    UnexpectedArtificialExit,
  };

  struct BoundarySummary
  {
    G4double thickness_um = 0.0;
    std::string placement_file;
    std::int64_t n_physical_surface_exit = 0;
    G4double sum_physical_surface_exit_ekin_post_keV = 0.0;
    std::int64_t n_unexpected_artificial_exit = 0;
    G4double sum_unexpected_artificial_exit_ekin_post_keV = 0.0;
    G4double max_unexpected_artificial_exit_ekin_post_keV = 0.0;
    std::int64_t n_unexpected_artificial_exit_alpha = 0;
    G4double sum_unexpected_artificial_exit_alpha_ekin_post_keV = 0.0;
    std::int64_t n_unexpected_artificial_exit_Li7 = 0;
    G4double sum_unexpected_artificial_exit_Li7_ekin_post_keV = 0.0;
    std::int64_t n_unexpected_bulk_exit = 0;
    G4double sum_unexpected_bulk_exit_ekin_post_keV = 0.0;
  };

  struct CaptureAnchorRow
  {
    std::string physical_event_uid;
    std::string source_event_uid;
    G4int eventID = -1;
    G4int record_index = -1;
    G4double thickness_um = 0.0;
    G4double bn_wt = 0.0;
    G4double zns_wt = 0.0;
    G4double capture_x_um = 0.0;
    G4double capture_y_um = 0.0;
    G4double source_x_um = 0.0;
    G4double source_y_um = 0.0;
    G4double depth_um = 0.0;
    std::string placement_file;
    G4double local_capture_x_um = 0.0;
    G4double local_capture_y_um = 0.0;
    G4double local_capture_z_um = 0.0;
    std::string surface_mode;
    G4double target_local_z_um = 0.0;
    G4double used_local_z_um = 0.0;
    G4double bn_center_x_um = 0.0;
    G4double bn_center_y_um = 0.0;
    G4double bn_center_z_um = 0.0;
    G4int alphali_replay_index = 0;
    G4int alphali_replay_count = 1;
    G4double trajectory_weight = 1.0;
  };

  struct UnexpectedBoundaryExitRow
  {
    std::string physical_event_uid;
    std::string source_event_uid;
    G4int eventID = -1;
    G4int record_index = -1;
    G4double thickness_um = 0.0;
    G4double bn_wt = 0.0;
    G4double zns_wt = 0.0;
    std::string placement_file;
    std::string surface_mode;
    std::string particle;
    G4int trackID = -1;
    G4int stepID = -1;
    std::string phase_pre;
    std::string phase_post;
    std::string exit_face;
    std::string exit_class;
    G4double ekin_post_keV = 0.0;
    G4double x_post_um = 0.0;
    G4double y_post_um = 0.0;
    G4double z_post_um = 0.0;
    G4int alphali_replay_index = 0;
    G4int alphali_replay_count = 1;
    G4double trajectory_weight = 1.0;
  };

  OutputMode GetOutputMode() const { return fOutputMode; }
  G4bool IsFullMode() const { return fOutputMode == OutputMode::Full; }
  G4bool IsSlimMode() const { return fOutputMode == OutputMode::Slim; }

  std::ofstream &GetFullStepCsv() { return fFullStepCsv; }
  std::ofstream &GetSlimTrackCsv() { return fSlimTrackCsv; }
  std::ofstream &GetUnexpectedBoundaryExitCsv() { return fUnexpectedBoundaryExitCsv; }
  const std::string &GetFullStepCsvPath() const { return fFullStepCsvPath; }
  const std::string &GetSlimTrackCsvPath() const { return fSlimTrackCsvPath; }
  const std::string &GetUnexpectedBoundaryExitCsvPath() const { return fUnexpectedBoundaryExitCsvPath; }
  G4bool IsFullStepCsvOpen() const { return fFullStepCsv.is_open(); }
  G4bool IsSlimTrackCsvOpen() const { return fSlimTrackCsv.is_open(); }
  G4bool IsUnexpectedBoundaryExitCsvOpen() const { return fUnexpectedBoundaryExitCsv.is_open(); }

  void AppendCaptureAnchor(const CaptureAnchorRow &row);
  void RecordBoundaryExit(const UnexpectedBoundaryExitRow &row,
                          BoundaryExitClass exitClass);

private:
  struct OutputPaths
  {
    std::string full_steps;
    std::string capture_anchors;
    std::string zns_track_steps;
    std::string unexpected_boundary_exits;
    std::string boundary_stop_summary;
  };

  std::string MakeOutputCsvPath() const;
  OutputPaths MakeOutputPathsFromInputPath(const std::string &inputPath) const;
  OutputPaths MakeOutputPathsFromThicknessTag(const std::string &thicknessTag) const;
  std::string RecordInputPathForSummary(const std::string &inputPath) const;
  std::string ExtractThicknessTagFromInputPath(const std::string &inputPath) const;
  std::string MakeThicknessTag(G4double thickness_um) const;
  void EnsureDataDirectory() const;
  OutputMode ReadOutputMode() const;
  void CloseOpenOutputs();
  void OpenOutputsForInputPath(const std::string &inputPath);
  void OpenOutputsForPaths(const OutputPaths &paths);
  void WriteFullStepCsvHeader();
  void WriteCaptureAnchorCsvHeader();
  void WriteSlimTrackCsvHeader();
  void WriteUnexpectedBoundaryExitCsvHeader();
  void WriteBoundarySummaryCsv();

private:
  using BoundarySummaryKey = std::pair<G4double, std::string>;

  const PrimaryGeneratorAction *fPrimaryAction;
  AnalysisConfig *fConfig;
  OutputMode fOutputMode;
  std::ofstream fFullStepCsv;
  std::string fFullStepCsvPath;
  std::ofstream fCaptureAnchorCsv;
  std::string fCaptureAnchorCsvPath;
  std::ofstream fSlimTrackCsv;
  std::string fSlimTrackCsvPath;
  std::ofstream fUnexpectedBoundaryExitCsv;
  std::string fUnexpectedBoundaryExitCsvPath;
  std::string fBoundaryStopSummaryCsvPath;
  std::string fCurrentOutputInputPath;
  std::map<BoundarySummaryKey, BoundarySummary> fBoundarySummaries;
};

#endif
