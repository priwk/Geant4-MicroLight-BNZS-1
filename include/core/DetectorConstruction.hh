#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4Orb;
class AnalysisConfig;
struct RvePlacementData;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  enum class Phase
  {
    BN,
    ZnS,
    Matrix,
    World,
    Unknown
  };

  struct SphereInfo
  {
    G4int particleId = -1;
    G4int phaseId = -1;
    G4int radiusBinId = -1;
    G4int radiusClassId = -1;
    Phase phase = Phase::Unknown;
    G4String radiusClass;
    G4ThreeVector center;
    G4double radius = 0.0;
  };

  struct GeometryCopyInfo
  {
    G4int copyId = -1;
    G4int sourceParticleId = -1;
    G4bool isPrimary = false;
    G4int shiftIx = 0;
    G4int shiftIy = 0;
    G4int shiftIz = 0;
    G4int radiusClassId = -1;
    Phase phase = Phase::Unknown;
    G4double radius = 0.0;
  };

  explicit DetectorConstruction(AnalysisConfig *config);
  ~DetectorConstruction() override;

  G4VPhysicalVolume *Construct() override;

  // ----- 几何结构参数 -----
  void SetScreenThicknessUm(G4double thicknessUm); // CSV 中的实际屏幕厚度
  void SetMicroThicknessUm(G4double thicknessUm);  // CSV 中的局部微结构 z 窗口
  void SetPatchXYUm(G4double patchXYUm);           // CSV 中的局部 x/y 窗口
  void SetWeightRatio(G4double bnWt, G4double znsWt);
  void SetPresetRatio(G4int presetIndex);

  G4double GetScreenThicknessUm() const { return fScreenThickness / um; }
  G4double GetMicroThicknessUm() const { return fMicroThickness / um; }
  G4double GetPatchXYUm() const { return fPatchXY / um; }
  G4double GetBnWt() const { return fBnWt; }
  G4double GetZnsWt() const { return fZnsWt; }

  G4double GetEffectiveLocalThickness() const;
  G4double GetFrontSurfaceZ() const { return GetBoxHalfZ(); }
  G4double GetBackSurfaceZ() const { return -GetBoxHalfZ(); }

  G4double GetSafeMarginXY() const { return fSafeMarginXY; }
  G4double GetMinParticleRadius(Phase phase) const;
  G4double GetMaxParticleRadius(Phase phase) const;
  G4double GetGlobalMaxParticleRadius() const;

  G4LogicalVolume *GetScoringVolume() const { return fScoringVolume; }
  G4LogicalVolume *GetMatrixLogical() const { return fMatrixLogical; }
  Phase GetPhaseFromLogicalVolume(const G4LogicalVolume *logical) const;
  G4bool IsBNLogical(const G4LogicalVolume *logical) const;
  G4bool IsZnSLogical(const G4LogicalVolume *logical) const;
  const GeometryCopyInfo *GetGeometryCopyInfo(G4int copyId) const;
  const std::vector<SphereInfo> &GetBNSpheres() const { return fBNSpheres; }
  const std::vector<SphereInfo> &GetZnSSpheres() const { return fZnSSpheres; }

  const G4String &GetLoadedPlacementFile() const { return fLoadedPlacementFile; }
  const G4String &GetLoadedPlacementFileForRecord() const { return fLoadedPlacementFileForRecord; }
  const G4String &GetLoadedPeriodicImagesFile() const { return fLoadedPeriodicImagesFile; }
  const G4String &GetLoadedPeriodicImagesFileForRecord() const { return fLoadedPeriodicImagesFileForRecord; }
  const G4String &GetLoadedPlacementSeedBase() const { return fLoadedPlacementSeedBase; }
  G4int GetPlacementFormatVersion() const { return fPlacementFormatVersion; }
  std::uint64_t GetPlacementSeed() const { return fPlacementSeed; }
  std::size_t GetUniqueParticleCount() const { return fBNSpheres.size() + fZnSSpheres.size(); }
  std::size_t GetGeometryCopyCount() const { return fGeometryCopyById.size(); }
  std::size_t GetRadiusClassCount() const { return fRadiusClassCache.size(); }
  G4double GetPlacementPhiAchieved() const { return fPlacementPhiAchieved; }
  G4double GetPlacementZnSToBNMassRatio() const { return fPlacementZnSToBNMassRatio; }

  G4double GetBoxHalfX() const { return 0.5 * fBoxX; }
  G4double GetBoxHalfY() const { return 0.5 * fBoxY; }
  G4double GetBoxHalfZ() const { return 0.5 * fBoxZ; }
  G4double GetBoxXUm() const { return fBoxX / um; }
  G4double GetBoxYUm() const { return fBoxY / um; }
  G4double GetBoxZUm() const { return fBoxZ / um; }
  G4double GetPatchHalfXUm() const { return GetBoxHalfX() / um; }
  G4double GetPatchHalfYUm() const { return GetBoxHalfY() / um; }
  G4double GetPatchHalfZUm() const { return GetBoxHalfZ() / um; }
  G4double GetPatchXUm() const { return GetBoxXUm(); }
  G4double GetPatchYUm() const { return GetBoxYUm(); }
  G4double GetPatchZUm() const { return GetBoxZUm(); }

  Phase FindPhaseAtPointUm(G4double xUm, G4double yUm, G4double zUm) const;
  Phase FindPhaseAtPoint(const G4ThreeVector &point) const;
  Phase FindPeriodicPhaseAtPoint(const G4ThreeVector &point) const;
  G4ThreeVector WrapToPrimaryCell(const G4ThreeVector &point) const;
  G4ThreeVector WrapToPrimaryCellInside(
      const G4ThreeVector &point,
      const G4ThreeVector &direction,
      G4double inwardEpsilon = 1.0e-4 * um) const;
  G4bool UsesPeriodicTransport() const
  {
    return fPlacementFormatVersion == 3 && !fLoadedPeriodicImagesFile.empty();
  }
  static const char *PhaseName(Phase phase);

private:
  void DefineMaterials();
  void NotifyGeometryChanged();
  G4VPhysicalVolume *ConstructVersion3(const RvePlacementData &data);

  static G4double HashToUnit(G4int ix, G4int iy, G4int iz, G4int salt);

  G4int ComputeTargetZnSCount(G4int placedBNCount, G4int usableZnSCount) const;
  G4bool IsInsideSafeXY(const G4ThreeVector &pos) const;
  void BuildPeriodicSpatialIndex();
  Phase FindPhaseInsidePrimary(
      const G4ThreeVector &point,
      G4bool useMinimumImage) const;
  void RegisterLogicalPhase(G4LogicalVolume *logical, Phase phase);

  struct RadiusClassGeometry
  {
    G4Orb *solid = nullptr;
    G4LogicalVolume *logical = nullptr;
    Phase phase = Phase::Unknown;
    G4double radius = 0.0;
  };

private:
  // ---- screen-level thickness ----
  G4double fScreenThickness; // CSV 中的实际全屏厚度

  // ---- local microstructure window ----
  G4double fMicroThickness; // CSV 中的局部 z 窗口
  G4double fPatchXY;        // CSV 中的局部 x/y 窗口
  G4double fBoxX;
  G4double fBoxY;
  G4double fBoxZ;

  // ---- legacy particle geometry ----
  G4double fBNRadius;
  G4double fZnSRadius;
  G4double fBNPitch;      // BN lattice pitch
  G4double fZnSPitch;     // ZnS lattice pitch
  G4double fOverlapGap;   // avoid touching overlaps
  G4double fSafeMarginXY; // avoid choosing BN too close to x/y patch edge

  // ---- background ----
  G4double fVoidVolumeFraction;

  // ---- composition ----
  G4double fBnWt;
  G4double fZnsWt;

  // ---- materials ----
  G4Material *fVacuumMaterial;
  G4Material *fMatrixMaterial;
  G4Material *fBNMaterial;
  G4Material *fZnSMaterial;

  // ---- logical volumes ----
  G4LogicalVolume *fScoringVolume;
  G4LogicalVolume *fWorldLogical;
  G4LogicalVolume *fMatrixLogical;
  G4LogicalVolume *fBNLogical;
  G4LogicalVolume *fZnSLogical;

  // ---- unique periodic particles and geometry lookup ----
  std::vector<SphereInfo> fBNSpheres;
  std::vector<SphereInfo> fZnSSpheres;
  std::vector<G4ThreeVector> fPlacedBNCenters;
  std::vector<G4ThreeVector> fSafeBNCenters;
  std::vector<G4ThreeVector> fUsableZnSCandidateCenters;
  std::vector<G4ThreeVector> fPlacedZnSCenters;
  std::unordered_map<int, RadiusClassGeometry> fRadiusClassCache;
  std::unordered_map<const G4LogicalVolume *, Phase> fLogicalPhase;
  std::unordered_map<int, GeometryCopyInfo> fGeometryCopyById;
  std::vector<const SphereInfo *> fSpatialSpheres;
  std::vector<std::vector<std::size_t>> fSpatialCells;
  G4int fSpatialNx;
  G4int fSpatialNy;
  G4int fSpatialNz;
  G4double fSpatialCellX;
  G4double fSpatialCellY;
  G4double fSpatialCellZ;

  // ---- shared analysis config ----
  AnalysisConfig *fConfig;

  G4String fLoadedPlacementFile;
  G4String fLoadedPlacementFileForRecord;
  G4String fLoadedPeriodicImagesFile;
  G4String fLoadedPeriodicImagesFileForRecord;
  G4String fLoadedPlacementSeedBase;
  G4int fPlacementFormatVersion;
  std::uint64_t fPlacementSeed;
  G4double fPlacementPhiAchieved;
  G4double fPlacementZnSToBNMassRatio;
};

#endif
