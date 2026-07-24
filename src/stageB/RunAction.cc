#include "RunAction.hh"

#include "AnalysisConfig.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <sstream>

namespace
{
std::string Trim(const std::string &value)
{
    const auto begin = value.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos)
        return "";
    const auto end = value.find_last_not_of(" \t\r\n");
    return value.substr(begin, end - begin + 1);
}
} // namespace

RunAction::RunAction(PrimaryGeneratorAction *primaryAction, AnalysisConfig *config)
    : G4UserRunAction(),
      fPrimaryAction(primaryAction),
      fConfig(config),
      fOutputMode(ReadOutputMode()),
      fFullStepCsv(),
      fFullStepCsvPath(""),
      fCaptureAnchorCsv(),
      fCaptureAnchorCsvPath(""),
      fSlimTrackCsv(),
      fSlimTrackCsvPath(""),
      fUnexpectedBoundaryExitCsv(),
      fUnexpectedBoundaryExitCsvPath(""),
      fBoundaryStopSummaryCsvPath(""),
      fCurrentOutputInputPath(""),
      fBoundarySummaries()
{
}

RunAction::~RunAction()
{
    CloseOpenOutputs();
}

void RunAction::SetPrimaryAction(const PrimaryGeneratorAction *primaryAction)
{
    fPrimaryAction = primaryAction;
}

void RunAction::EnsureDataDirectory() const
{
    namespace fs = std::filesystem;
    fs::path outputDir = fs::current_path().parent_path() / "Output" / "stageB";
    std::error_code ec;

    if (!fs::exists(outputDir))
    {
        fs::create_directories(outputDir, ec);
        if (ec)
        {
            G4cerr << "[RunAction] Warning: failed to create Output directory: "
                   << ec.message() << G4endl;
        }
    }
}

std::string RunAction::ExtractThicknessTagFromInputPath(const std::string &inputPath) const
{
    const std::string key = "_neutron_capture_positions.csv";

    std::size_t slashPos = inputPath.find_last_of("/\\");
    std::string fileName = (slashPos == std::string::npos)
                               ? inputPath
                               : inputPath.substr(slashPos + 1);

    std::size_t keyPos = fileName.find(key);
    if (keyPos == std::string::npos)
    {
        return "unknown";
    }

    return fileName.substr(0, keyPos);
}

std::string RunAction::MakeThicknessTag(G4double thickness_um) const
{
    std::ostringstream oss;
    const G4double rounded = std::round(thickness_um);
    if (std::abs(thickness_um - rounded) < 1.0e-9)
    {
        oss << static_cast<long long>(rounded);
    }
    else
    {
        oss << thickness_um;
    }
    return oss.str();
}

RunAction::OutputPaths RunAction::MakeOutputPathsFromThicknessTag(const std::string &tag) const
{
    namespace fs = std::filesystem;
    std::string ratioTag = "unknown";

    const char *ratioEnv = std::getenv("BNZS_OUTPUT_RATIO");
    if (ratioEnv != nullptr && std::string(ratioEnv).size() > 0)
    {
        ratioTag = ratioEnv;
    }

    const fs::path outDir = fs::current_path().parent_path() / "Output" / "stageB" / ratioTag;

    OutputPaths paths;
    paths.full_steps = (outDir / (tag + "_alpha_li_steps.csv")).string();
    paths.capture_anchors = (outDir / (tag + "_capture_anchors.csv")).string();
    paths.zns_track_steps = (outDir / (tag + "_zns_track_steps.csv")).string();
    paths.unexpected_boundary_exits = (outDir / (tag + "_unexpected_boundary_exits.csv")).string();
    paths.boundary_stop_summary = (outDir / (tag + "_boundary_stop_summary.csv")).string();
    return paths;
}

RunAction::OutputPaths RunAction::MakeOutputPathsFromInputPath(const std::string &inputPath) const
{
    const std::string tag = ExtractThicknessTagFromInputPath(inputPath);
    OutputPaths paths = MakeOutputPathsFromThicknessTag(tag);

    const char *ratioEnv = std::getenv("BNZS_OUTPUT_RATIO");
    if (ratioEnv != nullptr && std::string(ratioEnv).size() > 0)
    {
        return paths;
    }

    namespace fs = std::filesystem;
    const fs::path inPath(inputPath);
    if (!inPath.has_parent_path())
    {
        return paths;
    }

    const std::string parentName = inPath.parent_path().filename().string();
    if (parentName.empty())
    {
        return paths;
    }

    const fs::path outDir = fs::current_path().parent_path() / "Output" / "stageB" / parentName;
    paths.full_steps = (outDir / (tag + "_alpha_li_steps.csv")).string();
    paths.capture_anchors = (outDir / (tag + "_capture_anchors.csv")).string();
    paths.zns_track_steps = (outDir / (tag + "_zns_track_steps.csv")).string();
    paths.unexpected_boundary_exits = (outDir / (tag + "_unexpected_boundary_exits.csv")).string();
    paths.boundary_stop_summary = (outDir / (tag + "_boundary_stop_summary.csv")).string();
    return paths;
}

std::string RunAction::RecordInputPathForSummary(const std::string &inputPath) const
{
    return AnalysisConfig::PathForRecord(inputPath);
}

std::string RunAction::MakeOutputCsvPath() const
{
    if (fPrimaryAction)
    {
        return MakeOutputPathsFromThicknessTag(
                   MakeThicknessTag(fPrimaryAction->GetCurrentRecord().thickness_um))
            .full_steps;
    }

    return MakeOutputPathsFromThicknessTag("unknown").full_steps;
}

RunAction::OutputMode RunAction::ReadOutputMode() const
{
    const char *modeEnv = std::getenv("BNZS_STAGEB_OUTPUT_MODE");
    if (modeEnv == nullptr)
        return OutputMode::Full;

    const std::string mode = Trim(modeEnv);
    if (mode == "slim")
        return OutputMode::Slim;
    return OutputMode::Full;
}

void RunAction::CloseOpenOutputs()
{
    if (fFullStepCsv.is_open())
    {
        fFullStepCsv.flush();
        fFullStepCsv.close();
    }
    if (fCaptureAnchorCsv.is_open())
    {
        fCaptureAnchorCsv.flush();
        fCaptureAnchorCsv.close();
    }
    if (fSlimTrackCsv.is_open())
    {
        fSlimTrackCsv.flush();
        fSlimTrackCsv.close();
    }
    if (fUnexpectedBoundaryExitCsv.is_open())
    {
        fUnexpectedBoundaryExitCsv.flush();
        fUnexpectedBoundaryExitCsv.close();
    }
}

void RunAction::WriteFullStepCsvHeader()
{
    if (!fFullStepCsv.is_open())
        return;

    fFullStepCsv
        << "physical_event_uid,"
        << "eventID,"
        << "thickness_um,"
        << "bn_wt,"
        << "zns_wt,"
        << "capture_x_um,"
        << "capture_y_um,"
        << "source_x_um,"
        << "source_y_um,"
        << "depth_um,"
        << "placement_file,"
        << "local_capture_x_um,"
        << "local_capture_y_um,"
        << "local_capture_z_um,"
        << "surface_mode,"
        << "target_local_z_um,"
        << "used_local_z_um,"
        << "bn_center_x_um,"
        << "bn_center_y_um,"
        << "bn_center_z_um,"
        << "alphali_replay_index,"
        << "alphali_replay_count,"
        << "trackID,"
        << "stepID,"
        << "particle,"
        << "phase_pre,"
        << "phase_post,"
        << "x_pre_um,"
        << "y_pre_um,"
        << "z_pre_um,"
        << "x_post_um,"
        << "y_post_um,"
        << "z_post_um,"
        << "step_len_um,"
        << "edep_keV,"
        << "ekin_pre_keV,"
        << "ekin_post_keV,"
        << "source_event_uid,"
        << "record_index,"
        << "trajectory_weight,parentID,"
        << "cell_ix,cell_iy,cell_iz,"
        << "unwrapped_x_pre_um,unwrapped_y_pre_um,unwrapped_z_pre_um,"
        << "unwrapped_x_post_um,unwrapped_y_post_um,unwrapped_z_post_um,"
        << "screen_x_pre_um,screen_y_pre_um,screen_z_pre_um,screen_depth_pre_um,"
        << "screen_x_post_um,screen_y_post_um,screen_z_post_um,screen_depth_post_um"
        << "\n";
}

void RunAction::WriteCaptureAnchorCsvHeader()
{
    if (!fCaptureAnchorCsv.is_open())
        return;

    fCaptureAnchorCsv
        << "physical_event_uid,"
        << "source_event_uid,"
        << "eventID,"
        << "record_index,"
        << "thickness_um,"
        << "bn_wt,"
        << "zns_wt,"
        << "capture_x_um,"
        << "capture_y_um,"
        << "source_x_um,"
        << "source_y_um,"
        << "depth_um,"
        << "placement_file,"
        << "local_capture_x_um,"
        << "local_capture_y_um,"
        << "local_capture_z_um,"
        << "surface_mode,"
        << "target_local_z_um,"
        << "used_local_z_um,"
        << "bn_center_x_um,"
        << "bn_center_y_um,"
        << "bn_center_z_um,"
        << "alphali_replay_index,"
        << "alphali_replay_count,"
        << "trajectory_weight"
        << "\n";
}

void RunAction::WriteSlimTrackCsvHeader()
{
    if (!fSlimTrackCsv.is_open())
        return;

    fSlimTrackCsv
        << "physical_event_uid,"
        << "source_event_uid,"
        << "eventID,"
        << "record_index,"
        << "trackID,"
        << "stepID,"
        << "particle,"
        << "phase_pre,"
        << "phase_post,"
        << "x_pre_um,"
        << "y_pre_um,"
        << "z_pre_um,"
        << "x_post_um,"
        << "y_post_um,"
        << "z_post_um,"
        << "step_len_um,"
        << "edep_keV,"
        << "ekin_pre_keV,"
        << "ekin_post_keV,"
        << "alphali_replay_index,"
        << "alphali_replay_count,"
        << "trajectory_weight,parentID,"
        << "cell_ix,cell_iy,cell_iz,"
        << "unwrapped_x_pre_um,unwrapped_y_pre_um,unwrapped_z_pre_um,"
        << "unwrapped_x_post_um,unwrapped_y_post_um,unwrapped_z_post_um,"
        << "screen_x_pre_um,screen_y_pre_um,screen_z_pre_um,screen_depth_pre_um,"
        << "screen_x_post_um,screen_y_post_um,screen_z_post_um,screen_depth_post_um"
        << "\n";
}

void RunAction::WriteUnexpectedBoundaryExitCsvHeader()
{
    if (!fUnexpectedBoundaryExitCsv.is_open())
        return;

    fUnexpectedBoundaryExitCsv
        << "physical_event_uid,"
        << "source_event_uid,"
        << "eventID,"
        << "record_index,"
        << "thickness_um,"
        << "bn_wt,"
        << "zns_wt,"
        << "placement_file,"
        << "surface_mode,"
        << "particle,"
        << "trackID,"
        << "stepID,"
        << "phase_pre,"
        << "phase_post,"
        << "exit_face,"
        << "exit_class,"
        << "ekin_post_keV,"
        << "x_post_um,"
        << "y_post_um,"
        << "z_post_um,"
        << "alphali_replay_index,"
        << "alphali_replay_count,"
        << "trajectory_weight"
        << "\n";
}

void RunAction::OpenOutputsForInputPath(const std::string &inputPath)
{
    if (inputPath == fCurrentOutputInputPath)
        return;
    OpenOutputsForPaths(MakeOutputPathsFromInputPath(inputPath));
    fCurrentOutputInputPath = inputPath;
}

void RunAction::OpenOutputsForPaths(const OutputPaths &paths)
{
    const char *exceptionSource = "RunAction::OpenOutputsForPaths";

    if ((!IsFullMode() || (fFullStepCsv.is_open() && paths.full_steps == fFullStepCsvPath)) &&
        fCaptureAnchorCsv.is_open() &&
        paths.capture_anchors == fCaptureAnchorCsvPath &&
        fSlimTrackCsv.is_open() &&
        paths.zns_track_steps == fSlimTrackCsvPath &&
        fUnexpectedBoundaryExitCsv.is_open() &&
        paths.unexpected_boundary_exits == fUnexpectedBoundaryExitCsvPath)
    {
        return;
    }

    if (!fBoundaryStopSummaryCsvPath.empty())
    {
        WriteBoundarySummaryCsv();
    }

    CloseOpenOutputs();
    fBoundarySummaries.clear();

    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(fs::path(paths.full_steps).parent_path(), ec);
    if (ec)
    {
        G4cerr << "[RunAction] Warning: failed to create output directory: "
               << ec.message() << G4endl;
    }

    fFullStepCsvPath = paths.full_steps;
    fCaptureAnchorCsvPath = paths.capture_anchors;
    fSlimTrackCsvPath = paths.zns_track_steps;
    fUnexpectedBoundaryExitCsvPath = paths.unexpected_boundary_exits;
    fBoundaryStopSummaryCsvPath = paths.boundary_stop_summary;

    if (IsFullMode())
    {
        fFullStepCsv.open(fFullStepCsvPath.c_str(), std::ios::out);
        if (!fFullStepCsv.is_open())
        {
            G4Exception(exceptionSource,
                        "BNZS101", FatalException,
                        ("Failed to open output CSV: " + fFullStepCsvPath).c_str());
            return;
        }
        WriteFullStepCsvHeader();
    }

    fCaptureAnchorCsv.open(fCaptureAnchorCsvPath.c_str(), std::ios::out);
    if (!fCaptureAnchorCsv.is_open())
    {
        G4Exception(exceptionSource,
                    "BNZS102", FatalException,
                    ("Failed to open capture anchor CSV: " + fCaptureAnchorCsvPath).c_str());
        return;
    }
    WriteCaptureAnchorCsvHeader();

    fSlimTrackCsv.open(fSlimTrackCsvPath.c_str(), std::ios::out);
    if (!fSlimTrackCsv.is_open())
    {
        G4Exception(exceptionSource,
                    "BNZS103", FatalException,
                    ("Failed to open slim track CSV: " + fSlimTrackCsvPath).c_str());
        return;
    }
    WriteSlimTrackCsvHeader();

    fUnexpectedBoundaryExitCsv.open(fUnexpectedBoundaryExitCsvPath.c_str(), std::ios::out);
    if (!fUnexpectedBoundaryExitCsv.is_open())
    {
        G4Exception(exceptionSource,
                    "BNZS105", FatalException,
                    ("Failed to open unexpected boundary exit CSV: " + fUnexpectedBoundaryExitCsvPath).c_str());
        return;
    }
    WriteUnexpectedBoundaryExitCsvHeader();

    G4cout << "[RunAction] Switched Stage B output CSVs to:"
           << (IsFullMode() ? "\n  full alpha/Li tracks = " + fFullStepCsvPath : "")
           << "\n  anchors = " << fCaptureAnchorCsvPath
           << "\n  ZnS deposition steps = " << fSlimTrackCsvPath
           << "\n  unexpected exits = " << fUnexpectedBoundaryExitCsvPath
           << "\n  summary = " << fBoundaryStopSummaryCsvPath
           << G4endl;
}

void RunAction::SwitchOutputCsvForInputPath(const std::string &inputPath)
{
    OpenOutputsForInputPath(inputPath);
}

void RunAction::SwitchOutputCsvForThicknessTag(const std::string &thicknessTag)
{
    fCurrentOutputInputPath.clear();
    OpenOutputsForPaths(MakeOutputPathsFromThicknessTag(thicknessTag));
}

void RunAction::AppendCaptureAnchor(const CaptureAnchorRow &row)
{
    const BoundarySummaryKey key{row.thickness_um, row.placement_file};
    BoundarySummary &summary = fBoundarySummaries[key];
    summary.thickness_um = row.thickness_um;
    summary.placement_file = row.placement_file;

    if (!fCaptureAnchorCsv.is_open())
        return;

    fCaptureAnchorCsv
        << row.physical_event_uid << ","
        << row.source_event_uid << ","
        << row.eventID << ","
        << row.record_index << ","
        << row.thickness_um << ","
        << row.bn_wt << ","
        << row.zns_wt << ","
        << row.capture_x_um << ","
        << row.capture_y_um << ","
        << row.source_x_um << ","
        << row.source_y_um << ","
        << row.depth_um << ","
        << row.placement_file << ","
        << row.local_capture_x_um << ","
        << row.local_capture_y_um << ","
        << row.local_capture_z_um << ","
        << row.surface_mode << ","
        << row.target_local_z_um << ","
        << row.used_local_z_um << ","
        << row.bn_center_x_um << ","
        << row.bn_center_y_um << ","
        << row.bn_center_z_um << ","
        << row.alphali_replay_index << ","
        << row.alphali_replay_count << ","
        << row.trajectory_weight
        << "\n";
}

void RunAction::RecordBoundaryExit(const UnexpectedBoundaryExitRow &row,
                                   BoundaryExitClass exitClass)
{
    const BoundarySummaryKey key{row.thickness_um, row.placement_file};
    BoundarySummary &summary = fBoundarySummaries[key];
    summary.thickness_um = row.thickness_um;
    summary.placement_file = row.placement_file;

    if (exitClass == BoundaryExitClass::PhysicalSurfaceExit)
    {
        ++summary.n_physical_surface_exit;
        summary.sum_physical_surface_exit_ekin_post_keV += row.ekin_post_keV;
        return;
    }

    ++summary.n_unexpected_artificial_exit;
    summary.sum_unexpected_artificial_exit_ekin_post_keV += row.ekin_post_keV;
    summary.max_unexpected_artificial_exit_ekin_post_keV =
        std::max(summary.max_unexpected_artificial_exit_ekin_post_keV, row.ekin_post_keV);

    if (fUnexpectedBoundaryExitCsv.is_open())
    {
        fUnexpectedBoundaryExitCsv
            << row.physical_event_uid << ","
            << row.source_event_uid << ","
            << row.eventID << ","
            << row.record_index << ","
            << row.thickness_um << ","
            << row.bn_wt << ","
            << row.zns_wt << ","
            << row.placement_file << ","
            << row.surface_mode << ","
            << row.particle << ","
            << row.trackID << ","
            << row.stepID << ","
            << row.phase_pre << ","
            << row.phase_post << ","
            << row.exit_face << ","
            << row.exit_class << ","
            << row.ekin_post_keV << ","
            << row.x_post_um << ","
            << row.y_post_um << ","
            << row.z_post_um << ","
            << row.alphali_replay_index << ","
            << row.alphali_replay_count << ","
            << row.trajectory_weight
            << "\n";
    }

    if (row.particle == "alpha")
    {
        ++summary.n_unexpected_artificial_exit_alpha;
        summary.sum_unexpected_artificial_exit_alpha_ekin_post_keV += row.ekin_post_keV;
    }
    else if (row.particle == "Li7")
    {
        ++summary.n_unexpected_artificial_exit_Li7;
        summary.sum_unexpected_artificial_exit_Li7_ekin_post_keV += row.ekin_post_keV;
    }

    if (row.surface_mode == "bulk")
    {
        ++summary.n_unexpected_bulk_exit;
        summary.sum_unexpected_bulk_exit_ekin_post_keV += row.ekin_post_keV;
    }
}

void RunAction::WriteBoundarySummaryCsv()
{
    if (fBoundaryStopSummaryCsvPath.empty())
        return;

    std::ofstream out(fBoundaryStopSummaryCsvPath.c_str(), std::ios::out);
    if (!out.is_open())
    {
        G4Exception("RunAction::WriteBoundarySummaryCsv",
                    "BNZS104", FatalException,
                    ("Failed to open boundary summary CSV: " + fBoundaryStopSummaryCsvPath).c_str());
        return;
    }

    out << "thickness_um,"
        << "placement_file,"
        << "n_physical_surface_exit,"
        << "sum_physical_surface_exit_ekin_post_keV,"
        << "n_unexpected_artificial_exit,"
        << "sum_unexpected_artificial_exit_ekin_post_keV,"
        << "max_unexpected_artificial_exit_ekin_post_keV,"
        << "n_unexpected_artificial_exit_alpha,"
        << "sum_unexpected_artificial_exit_alpha_ekin_post_keV,"
        << "n_unexpected_artificial_exit_Li7,"
        << "sum_unexpected_artificial_exit_Li7_ekin_post_keV,"
        << "n_unexpected_bulk_exit,"
        << "sum_unexpected_bulk_exit_ekin_post_keV"
        << "\n";

    std::int64_t unexpectedCount = 0;
    G4double unexpectedSum = 0.0;
    G4double unexpectedMax = 0.0;

    for (const auto &[key, summary] : fBoundarySummaries)
    {
        (void)key;
        out << summary.thickness_um << ","
            << summary.placement_file << ","
            << summary.n_physical_surface_exit << ","
            << summary.sum_physical_surface_exit_ekin_post_keV << ","
            << summary.n_unexpected_artificial_exit << ","
            << summary.sum_unexpected_artificial_exit_ekin_post_keV << ","
            << summary.max_unexpected_artificial_exit_ekin_post_keV << ","
            << summary.n_unexpected_artificial_exit_alpha << ","
            << summary.sum_unexpected_artificial_exit_alpha_ekin_post_keV << ","
            << summary.n_unexpected_artificial_exit_Li7 << ","
            << summary.sum_unexpected_artificial_exit_Li7_ekin_post_keV << ","
            << summary.n_unexpected_bulk_exit << ","
            << summary.sum_unexpected_bulk_exit_ekin_post_keV
            << "\n";

        unexpectedCount += summary.n_unexpected_artificial_exit;
        unexpectedSum += summary.sum_unexpected_artificial_exit_ekin_post_keV;
        unexpectedMax = std::max(
            unexpectedMax,
            summary.max_unexpected_artificial_exit_ekin_post_keV);
    }

    G4cout << "[RunAction] unexpected artificial exits = " << unexpectedCount
           << "\n[RunAction] sum ekin_post = " << unexpectedSum << " keV"
           << "\n[RunAction] max ekin_post = " << unexpectedMax << " keV"
           << G4endl;
}

void RunAction::BeginOfRunAction(const G4Run *run)
{
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);

    EnsureDataDirectory();

    if (fPrimaryAction)
    {
        OpenOutputsForInputPath(fPrimaryAction->GetLoadedInputFile());
    }
    else
    {
        OpenOutputsForInputPath("unknown_neutron_capture_positions.csv");
    }

    G4cout << "\n[RunAction] Begin run " << run->GetRunID();
    if (IsFullMode())
        G4cout << "\n  alpha/Li full CSV = " << fFullStepCsvPath;
    G4cout << "\n  capture anchors = " << fCaptureAnchorCsvPath
           << "\n  ZnS deposition CSV = " << fSlimTrackCsvPath
           << "\n  unexpected exits = " << fUnexpectedBoundaryExitCsvPath
           << "\n  boundary summary = " << fBoundaryStopSummaryCsvPath;

    if (fPrimaryAction)
    {
        G4cout
            << "\n  input CSV  = " << RecordInputPathForSummary(fPrimaryAction->GetLoadedInputFile())
            << "\n  streamed input records so far = " << fPrimaryAction->GetTotalLoadedEvents();
    }

    G4cout << G4endl;
}

void RunAction::EndOfRunAction(const G4Run *run)
{
    WriteBoundarySummaryCsv();

    CloseOpenOutputs();

    G4cout << "\n[RunAction] End run " << run->GetRunID();
    if (IsFullMode())
        G4cout << "\n  alpha/Li full CSV = " << fFullStepCsvPath;
    G4cout << "\n  anchors CSV = " << fCaptureAnchorCsvPath
           << "\n  ZnS deposition CSV = " << fSlimTrackCsvPath
           << "\n  unexpected exits CSV = " << fUnexpectedBoundaryExitCsvPath
           << "\n  summary CSV = " << fBoundaryStopSummaryCsvPath;
    G4cout << G4endl;
}
