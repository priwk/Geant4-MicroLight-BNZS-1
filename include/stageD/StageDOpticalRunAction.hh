#ifndef StageDOpticalRunAction_h
#define StageDOpticalRunAction_h 1

#include "G4UserRunAction.hh"
#include "StageDOpticalStats.hh"

#include <array>
#include <fstream>
#include <string>
#include <vector>

class G4Run;
class AnalysisConfig;
class DetectorConstruction;

class StageDOpticalRunAction : public G4UserRunAction
{
public:
  explicit StageDOpticalRunAction(AnalysisConfig *config);
  ~StageDOpticalRunAction() override;

  void BeginOfRunAction(const G4Run *run) override;
  void EndOfRunAction(const G4Run *run) override;

  void RecordPhotonEvent(const StageDPhotonEventRecord &event);
  void RecordReentryDiagnostic(const StageDReentryDiagnosticRecord &record);
  void SetReentryPortalSummary(const StageDReentryPortalSummary &summary);

  const std::string &GetEventsCsvPath() const { return fEventsCsvPath; }
  const std::string &GetSummaryCsvPath() const { return fSummaryCsvPath; }
  const std::string &GetPhaseFunctionCsvPath() const { return fPhaseFunctionCsvPath; }
  const std::string &GetReentryDiagnosticsCsvPath() const { return fReentryDiagnosticsCsvPath; }
  const std::string &GetRatioTag() const { return fRatioTag; }
  const std::string &GetPlacementStem() const { return fPlacementStem; }

private:
  std::string MakeRatioTag() const;
  std::string MakePlacementStem() const;
  std::string MakeWavelengthTag() const;
  std::string ResolveOutputDirectory() const;
  void OpenOutputs();
  void WriteEventHeader();
  void WriteReentryDiagnosticsHeader();
  void WriteSummaryFile() const;
  void WritePhaseFunctionFile() const;
  void WriteGeometryMetadataFile(const DetectorConstruction *detector) const;

private:
  AnalysisConfig *fConfig;
  std::ofstream fEventsCsv;
  std::ofstream fReentryDiagnosticsCsv;
  std::string fEventsCsvPath;
  std::string fSummaryCsvPath;
  std::string fPhaseFunctionCsvPath;
  std::string fReentryDiagnosticsCsvPath;
  std::string fOutputDir;
  std::string fRatioTag;
  std::string fPlacementFile;
  std::string fPlacementStem;
  std::vector<StageDPhotonEventRecord> fEvents;
  StageDReentryPortalSummary fReentryPortalSummary;
};

#endif
