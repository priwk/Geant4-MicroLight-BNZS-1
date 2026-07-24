#ifndef RvePlacementReader_h
#define RvePlacementReader_h 1

#include "G4ThreeVector.hh"

#include <cstdint>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

enum class RvePhase
{
  BN = 0,
  ZnS = 1,
  Unknown = -1
};

struct RveMetadata
{
  int formatVersion = 0;
  std::string generator;
  std::string boundary;
  std::string coordinateConvention;
  bool oneRowPerUniquePeriodicParticle = false;
  std::string geant4MaterialKey;
  std::string geant4SolidCacheKey;
  double boxXUm = 0.0;
  double boxYUm = 0.0;
  double boxZUm = 0.0;
  double phiTarget = 0.0;
  double phiAchieved = 0.0;
  double overlapGapUm = 0.0;
  double finalMaxOverlapViolationUm = 0.0;
  std::uint64_t seed = 0;
  std::string radiusStatsCsv;
  std::string periodicImagesCsv;
  std::unordered_map<std::string, std::string> values;
};

struct RveParticle
{
  int particleId = -1;
  RvePhase phase = RvePhase::Unknown;
  int phaseId = -1;
  int radiusBinId = -1;
  int radiusClassId = -1;
  std::string radiusClass;
  double radiusUm = 0.0;
  G4ThreeVector center;
};

struct RveGeometryCopy
{
  int copyId = -1;
  int sourceParticleId = -1;
  bool isPrimary = false;
  int shiftIx = 0;
  int shiftIy = 0;
  int shiftIz = 0;
  RvePhase phase = RvePhase::Unknown;
  int phaseId = -1;
  int radiusBinId = -1;
  int radiusClassId = -1;
  std::string radiusClass;
  double radiusUm = 0.0;
  G4ThreeVector center;
};

struct RveRadiusClass
{
  int radiusClassId = -1;
  RvePhase phase = RvePhase::Unknown;
  int phaseId = -1;
  std::string radiusClass;
  double radiusUm = 0.0;
};

struct RvePlacementData
{
  RveMetadata metadata;
  std::filesystem::path placementPath;
  std::filesystem::path periodicImagesPath;
  std::vector<RveParticle> particles;
  std::vector<RveGeometryCopy> geometryCopies;
  std::unordered_map<int, RveRadiusClass> radiusClasses;
  std::unordered_map<int, std::size_t> particleIndexById;
  std::unordered_map<int, std::size_t> copyIndexById;
  std::size_t bnCount = 0;
  std::size_t znsCount = 0;
  double bnVolumeUm3 = 0.0;
  double znsVolumeUm3 = 0.0;
  double phiAchieved = 0.0;
  double znsToBnMassRatio = 0.0;
};

class RvePlacementReader
{
public:
  static bool IsFormatVersion3(const std::filesystem::path &path);
  static bool IsMainPlacementCandidate(const std::filesystem::path &path);
  static RvePlacementData ReadVersion3(
      const std::filesystem::path &placementPath,
      const std::filesystem::path &periodicImagesOverride = {},
      bool readPeriodicImages = true);
};

#endif
