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
class DetectorConstruction;

struct CaptureRecord
{
  G4int eventID = -1;
  G4int record_index = -1;
  std::string input_file_uid;
  G4int placement_replay_index = 0;
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
  void ValidateDetectorAgainstInput() const;

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
  G4int GetCurrentSelectedBNParticleId() const { return fCurrentSelectedBNParticleId; }
  G4int GetCurrentSelectedBNRadiusClassId() const { return fCurrentSelectedBNRadiusClassId; }
  G4double GetCurrentSelectedBNRadius() const { return fCurrentSelectedBNRadius; }
  G4int GetCurrentBNImageIx() const { return fCurrentBNImageIx; }
  G4int GetCurrentBNImageIy() const { return fCurrentBNImageIy; }
  G4int GetCurrentBNImageIz() const { return fCurrentBNImageIz; }
  const std::string &GetCurrentReactionBranch() const { return fCurrentReactionBranch; }
  const G4ThreeVector &GetCurrentLaunchDirection() const { return fCurrentLaunchDirection; }
  G4bool HasValidCurrentReplay() const { return fCurrentReplayValid; }
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
  G4bool IsInputThicknessCompatible(G4double thickness_um, G4double localT_um) const;
  G4int ReadAlphaLiReplayPerCapture() const;
  G4bool PrepareCurrentCaptureReplayState();
  std::string InputFileUid(const std::string &path);

  // ---- event classification ----
  std::string DetermineSurfaceMode(const CaptureRecord &rec) const;
  G4double DetermineTargetLocalZ(const CaptureRecord &rec, const std::string &mode) const;

  // ---- BN selection and point sampling ----
  G4bool SelectBNSphereForTargetZ(
      G4double targetZ,
      G4int &chosenSphereIndex,
      G4ThreeVector &chosenCenter,
      G4double &chosenRadius,
      G4double &usedZ,
      G4bool &usedFallback);

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
      G4int &chosenSphereIndex,
      G4ThreeVector &chosenCenter,
      G4ThreeVector &capturePoint,
      G4bool &usedFallback);

  void EnsureBNSamplingCache();
  G4int BNSamplingZBin(G4double z) const;

  G4ThreeVector SamplePointInSphereVolume(
      const G4ThreeVector &center,
      G4double sphereRadius) const;

  // ---- reaction generation ----
  void GenerateReactionProducts(
      G4Event *event,
      const G4ThreeVector &position,
      G4bool useGroundStateBranch,
      const G4ThreeVector &launchDirection) const;

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
  std::unordered_map<std::string, std::string> fInputFileUidByPath;

  CaptureRecord fFirstRecordForGeometry;
  G4bool fHasFirstRecordForGeometry = false;
  G4bool fNoMoreInput = false;
  std::size_t fTotalStreamedRecords = 0;

  G4int fAlphaLiReplayPerCapture = 1;
  G4int fCurrentAlphaLiReplayIndex = 0;
  G4int fRemainingReplaysForCurrentCapture = 0;
  G4bool fCurrentReplayValid = false;
  std::string fInitializedCaptureCsvPath;
  std::string fInitializedCaptureInputDir;

  // current event cache
  CaptureRecord fCurrentRecord;
  G4ThreeVector fCurrentLocalCapturePosition;
  G4ThreeVector fCurrentSelectedBNCenter;
  std::string fCurrentSurfaceMode;
  G4double fCurrentTargetLocalZ;
  G4double fCurrentUsedLocalZ;
  G4int fCurrentSelectedBNParticleId;
  G4int fCurrentSelectedBNRadiusClassId;
  G4double fCurrentSelectedBNRadius;
  G4int fCurrentBNImageIx;
  G4int fCurrentBNImageIy;
  G4int fCurrentBNImageIz;
  std::string fCurrentReactionBranch;
  G4ThreeVector fCurrentLaunchDirection;

  const DetectorConstruction *fBNSamplingDetectorIdentity;
  std::size_t fBNSamplingSphereCount;
  G4double fBNSamplingBoxZ;
  G4double fBNSamplingZBinWidth;
  G4int fBNSamplingZBinCount;
  std::vector<G4double> fBNVolumeCdf;
  std::vector<std::vector<G4int>> fBNZBins;
};

#endif
