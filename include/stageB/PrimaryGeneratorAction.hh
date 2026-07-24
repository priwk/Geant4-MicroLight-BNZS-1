#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

class G4Event;
class AnalysisConfig;

struct CaptureRecord
{
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
};

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  explicit PrimaryGeneratorAction(AnalysisConfig *config = nullptr);
  ~PrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event *event) override;
  void RefreshInputSelectionFromConfig();

  G4ParticleGun *GetParticleGun() const { return fParticleGun; }

  // ---- current event metadata for later use in EventAction / SteppingAction ----
  const CaptureRecord &GetCurrentRecord() const { return fCurrentRecord; }
  const G4ThreeVector &GetCurrentLocalCapturePosition() const { return fCurrentLocalCapturePosition; }
  const G4ThreeVector &GetCurrentSelectedBNCenter() const { return fCurrentSelectedBNCenter; }
  const std::string &GetCurrentSurfaceMode() const { return fCurrentSurfaceMode; }
  G4double GetCurrentTargetLocalZ() const { return fCurrentTargetLocalZ; }
  G4double GetCurrentUsedLocalZ() const { return fCurrentUsedLocalZ; }
  G4int GetCurrentAlphaLiReplayIndex() const { return fCurrentAlphaLiReplayIndex; }
  G4int GetCurrentAlphaLiReplayCount() const { return fAlphaLiReplayPerCapture; }
  G4double GetCurrentTrajectoryWeight() const;
  std::string MakeCurrentPhysicalEventUid() const;
  std::string MakeCurrentSourceEventUid() const;
  G4int GetTotalLoadedEvents() const { return static_cast<G4int>(fTotalStreamedRecords); }
  const std::string &GetLoadedInputFile() const { return fCurrentInputFile; }
  const std::string &GetCurrentRecordInputFile() const { return fCurrentRecordInputFile; }

private:
  using HeaderIndex = std::unordered_map<std::string, std::size_t>;

  // ---- input handling ----
  void InitializeInputStreaming();
  std::vector<std::string> FindInputCsvFiles() const;
  G4bool OpenNextInputFile();
  G4bool ReadNextRecord(CaptureRecord &rec);
  G4bool ReadFirstValidRecordFromFile(const std::string &path, CaptureRecord &rec) const;
  G4bool ReadHeaderIndex(std::istream &input, HeaderIndex &headerIndex) const;
  G4bool ParseOneRecordLine(
      const std::string &line,
      const HeaderIndex &headerIndex,
      G4int fallbackRecordIndex,
      CaptureRecord &rec) const;
  void ConfigureDetectorFromInput();
  G4bool IsInputThicknessCompatible(G4double thickness_um, G4double localT_um) const;
  G4int ReadAlphaLiReplayPerCapture() const;
  G4bool PrepareCurrentCaptureReplayState();

  // ---- event classification ----
  std::string DetermineSurfaceMode(const CaptureRecord &rec) const;
  G4double DetermineTargetLocalZ(const CaptureRecord &rec, const std::string &mode) const;

  // ---- BN selection and point sampling ----
  G4bool SelectBNSphereForTargetZ(
      G4double targetZ,
      G4ThreeVector &chosenCenter,
      G4double &chosenRadius,
      G4double &usedZ,
      G4bool &usedFallback) const;

  G4ThreeVector SamplePointInSphereSlice(
      const G4ThreeVector &center,
      G4double zSlice,
      G4double sphereRadius) const;

  G4bool SampleSafePointInSphereSlice(
      const G4ThreeVector &center,
      G4double zSlice,
      G4double sphereRadius,
      G4ThreeVector &point) const;

  G4bool SampleBulkCapturePoint(
      G4ThreeVector &chosenCenter,
      G4ThreeVector &capturePoint,
      G4bool &usedFallback) const;

  G4ThreeVector SamplePointInSphereVolume(
      const G4ThreeVector &center,
      G4double sphereRadius) const;

  // ---- reaction generation ----
  void GenerateReactionProducts(
      G4Event *event,
      const G4ThreeVector &position,
      G4bool useGroundStateBranch) const;

private:
  AnalysisConfig *fConfig;
  G4ParticleGun *fParticleGun;

  // ---- multi-file streaming state ----
  std::vector<std::string> fInputFiles;
  std::size_t fCurrentFileIndex = 0;
  std::ifstream fCurrentInputStream;
  std::string fCurrentInputFile;
  std::string fCurrentRecordInputFile;
  HeaderIndex fCurrentHeaderIndex;
  G4int fCurrentInputRecordCounter = 0;

  CaptureRecord fFirstRecordForGeometry;
  G4bool fHasFirstRecordForGeometry = false;
  G4bool fNoMoreInput = false;
  std::size_t fTotalStreamedRecords = 0;

  G4int fAlphaLiReplayPerCapture = 1;
  G4int fCurrentAlphaLiReplayIndex = 0;
  G4int fRemainingReplaysForCurrentCapture = 0;
  std::string fInitializedCaptureCsvPath;
  std::string fInitializedCaptureInputDir;

  // current event cache
  CaptureRecord fCurrentRecord;
  G4ThreeVector fCurrentLocalCapturePosition;
  G4ThreeVector fCurrentSelectedBNCenter;
  std::string fCurrentSurfaceMode;
  G4double fCurrentTargetLocalZ;
  G4double fCurrentUsedLocalZ;
};

#endif
