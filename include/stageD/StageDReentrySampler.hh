#ifndef StageDReentrySampler_h
#define StageDReentrySampler_h 1

#include "DetectorConstruction.hh"
#include "StageDOpticalStats.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <array>
#include <string>
#include <vector>

class AnalysisConfig;

class StageDReentrySampler
{
public:
  struct ReentryContext
  {
    DetectorConstruction::Phase phase = DetectorConstruction::Phase::Unknown;

    G4ThreeVector prePos;
    G4ThreeVector postPos;
    G4ThreeVector oldDir;

    G4ThreeVector exitPoint;
    G4ThreeVector exitInsidePoint;

    G4double wavelengthNm = 0.0;

    G4int eventID = -1;
    G4int reentryIndex = 0;
  };

  struct ReentryDiagnostics
  {
    std::string strategy;
    std::string fallbackLevel;

    DetectorConstruction::Phase exitPhase = DetectorConstruction::Phase::Unknown;
    DetectorConstruction::Phase entryPhase = DetectorConstruction::Phase::Unknown;

    G4ThreeVector exitInsidePoint;
    G4ThreeVector entryPoint;

    G4double particleQExit = -1.0;
    G4double particleQEntry = -1.0;
    G4double particleMuExit = 0.0;
    G4double particleMuEntry = 0.0;

    G4double matrixClearanceExitUm = -1.0;
    G4double matrixClearanceEntryUm = -1.0;
    std::string matrixNearestPhaseExit;
    std::string matrixNearestPhaseEntry;
    G4int matrixClearanceBinExit = -1;
    G4int matrixClearanceBinEntry = -1;

    G4int trials = 0;
  };

  explicit StageDReentrySampler(const DetectorConstruction *detector,
                                const AnalysisConfig *config);

  G4bool SampleReentry(const ReentryContext &ctx,
                       G4ThreeVector &newPosition,
                       ReentryDiagnostics &diag) const;

  G4bool InsideRveBox(const G4ThreeVector &point) const;
  DetectorConstruction::Phase FastPhaseAtPointForReentry(const G4ThreeVector &point) const;
  const StageDReentryPortalSummary &GetPortalSummary() const { return fPortalSummary; }

private:
  struct MicroSphere
  {
    DetectorConstruction::Phase phase = DetectorConstruction::Phase::Unknown;
    G4ThreeVector center;
    G4double radius = 0.0;
    G4double radius2 = 0.0;
    G4double visibleVolume = 0.0;
    G4bool isClipped = false;
  };

  struct NearestSurface
  {
    G4double clearance = 1.0e9;
    DetectorConstruction::Phase nearestPhase = DetectorConstruction::Phase::Unknown;
    G4int sphereId = -1;
    G4ThreeVector surfaceNormal;
  };

  struct MatrixPortal
  {
    G4ThreeVector position;
    G4double clearance = 0.0;
    DetectorConstruction::Phase nearestPhase = DetectorConstruction::Phase::Unknown;
    G4int clearanceBin = 0;
  };

  enum class PortalFace : std::size_t
  {
    PosX = 0,
    NegX = 1,
    PosY = 2,
    NegY = 3,
    PosZ = 4,
    NegZ = 5,
    Count = 6
  };

  enum class PortalPhaseBucket : std::size_t
  {
    BN = 0,
    ZnS = 1,
    Any = 2,
    Count = 3
  };

  static constexpr std::size_t kPortalFaceCount = 6;
  static constexpr std::size_t kPortalPhaseBucketCount = 3;
  static constexpr std::size_t kClearanceBinCount = 4;

  using PortalBinPools = std::array<std::vector<MatrixPortal>, kClearanceBinCount>;
  using PortalPhasePools = std::array<PortalBinPools, kPortalPhaseBucketCount>;
  using PortalFacePools = std::array<PortalPhasePools, kPortalFaceCount>;

  void BuildSphereCache();
  void BuildSpatialGrid();
  void BuildMatrixPortalPool();

  G4bool SampleParticleSphereQMuReentry(const ReentryContext &ctx,
                                        G4ThreeVector &newPosition,
                                        ReentryDiagnostics &diag) const;
  G4bool SampleParticleSphereQOnlyReentry(const ReentryContext &ctx,
                                          const MicroSphere *oldSphere,
                                          G4double q,
                                          G4ThreeVector &newPosition,
                                          ReentryDiagnostics &diag) const;
  G4bool SampleParticleVolumeUniformReentry(const ReentryContext &ctx,
                                            G4ThreeVector &newPosition,
                                            ReentryDiagnostics &diag) const;
  G4bool SampleMatrixClearanceBinnedPortalReentry(const ReentryContext &ctx,
                                                  G4ThreeVector &newPosition,
                                                  ReentryDiagnostics &diag) const;
  G4bool SampleRandomMatrixDebugReentry(const ReentryContext &ctx,
                                        G4ThreeVector &newPosition,
                                        ReentryDiagnostics &diag) const;

  const MicroSphere *FindContainingSphereOnly(DetectorConstruction::Phase phase,
                                              const G4ThreeVector &position) const;
  NearestSurface FindNearestParticleSurface(const G4ThreeVector &position) const;

  G4bool ClampInsideSameSphereRoundoffOnly(G4ThreeVector &point,
                                           const MicroSphere &sphere) const;
  G4bool ValidatePhase(const G4ThreeVector &point,
                       DetectorConstruction::Phase phase) const;

  G4ThreeVector RandomUnitVector() const;
  G4ThreeVector RandomUnitVectorWithFixedDot(const G4ThreeVector &axis,
                                             G4double mu) const;
  G4ThreeVector RandomPointInMatrixBox() const;

  G4double SphereSelectionWeight(const MicroSphere &sphere) const;
  G4int SampleWeightedSamePhaseSphereId(DetectorConstruction::Phase phase) const;
  G4int ClearanceBin(G4double clearance) const;

  std::size_t FlatCellIndex(G4int ix, G4int iy, G4int iz) const;
  G4bool PointToCell(const G4ThreeVector &point,
                     G4int &ix,
                     G4int &iy,
                     G4int &iz) const;
  void CollectCandidateSphereIds(const G4ThreeVector &point,
                                 std::vector<G4int> &candidateIds) const;
  G4ThreeVector MinimumImageDelta(const G4ThreeVector &point,
                                  const G4ThreeVector &center) const;

  G4double HalfExtentX() const;
  G4double HalfExtentY() const;
  G4double HalfExtentZ() const;
  G4ThreeVector InwardNormal(PortalFace face) const;
  G4double FaceArea(PortalFace face) const;
  G4ThreeVector PointOnVirtualFace(PortalFace face,
                                   G4double u01,
                                   G4double v01) const;

  static std::size_t PhaseVectorIndex(DetectorConstruction::Phase phase);
  static std::size_t PortalPhaseBucketIndex(DetectorConstruction::Phase phase);
  static const char *PhaseNameOrUnknown(DetectorConstruction::Phase phase);

private:
  const DetectorConstruction *fDetector;
  const AnalysisConfig *fConfig;

  std::vector<MicroSphere> fSpheres;
  std::array<std::vector<G4int>, 2> fSphereIdsByPhase;
  std::array<std::vector<G4double>, 2> fSphereWeightCdfByPhase;

  G4double fHalfX;
  G4double fHalfY;
  G4double fHalfZ;
  G4double fMaxSphereRadius;
  G4double fGridCellSize;
  G4ThreeVector fGridMinCorner;
  G4int fGridNx;
  G4int fGridNy;
  G4int fGridNz;
  std::vector<std::vector<G4int>> fGridCells;
  mutable std::vector<G4int> fCandidateVisitStamps;
  mutable G4int fCandidateVisitToken;

  G4double fPortalHalfX;
  G4double fPortalHalfY;
  G4double fPortalHalfZ;
  G4double fPortalMargin;
  G4int fPortalNu;
  G4int fPortalNv;
  std::array<G4double, 3> fClearanceBinEdges;
  G4int fMaxParticleReentryTrials;
  G4int fMaxPortalFallbackLevel;

  PortalFacePools fMatrixPortals;
  StageDReentryPortalSummary fPortalSummary;
};

#endif
