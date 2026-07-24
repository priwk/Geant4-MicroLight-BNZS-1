#include "RvePlacementReader.hh"

#include "G4Exception.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace
{
  constexpr double kBnDensityGPerCm3 = 2.10;
  constexpr double kZnSDensityGPerCm3 = 4.09;
  constexpr double kMetadataRelativeTolerance = 2.0e-9;
  constexpr double kCoordinateToleranceUm = 2.0e-9;
  constexpr double kRadiusToleranceUm = 2.0e-9;

  [[noreturn]] void Fail(const std::filesystem::path &path,
                         std::size_t lineNumber,
                         const std::string &message)
  {
    std::ostringstream out;
    out << path.string();
    if (lineNumber > 0)
      out << ':' << lineNumber;
    out << ": " << message;
    G4Exception("RvePlacementReader", "BNZS_RVE_001", FatalException, out.str().c_str());
    throw std::runtime_error(out.str());
  }

  std::string Trim(const std::string &value)
  {
    const auto begin = std::find_if_not(value.begin(), value.end(),
                                        [](unsigned char c) { return std::isspace(c); });
    const auto end = std::find_if_not(value.rbegin(), value.rend(),
                                      [](unsigned char c) { return std::isspace(c); }).base();
    return (begin < end) ? std::string(begin, end) : std::string();
  }

  std::vector<std::string> SplitCsv(const std::string &line)
  {
    std::vector<std::string> fields;
    std::string field;
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i)
    {
      const char c = line[i];
      if (c == '"')
      {
        if (quoted && i + 1 < line.size() && line[i + 1] == '"')
        {
          field.push_back('"');
          ++i;
        }
        else
        {
          quoted = !quoted;
        }
      }
      else if (c == ',' && !quoted)
      {
        fields.push_back(Trim(field));
        field.clear();
      }
      else
      {
        field.push_back(c);
      }
    }
    fields.push_back(Trim(field));
    return fields;
  }

  std::unordered_map<std::string, std::size_t> BuildColumnIndex(
      const std::filesystem::path &path,
      std::size_t lineNumber,
      const std::vector<std::string> &headers,
      const std::vector<std::string> &required)
  {
    std::unordered_map<std::string, std::size_t> index;
    for (std::size_t i = 0; i < headers.size(); ++i)
    {
      const std::string name = Trim(headers[i]);
      if (name.empty())
        Fail(path, lineNumber, "empty CSV column name");
      if (!index.emplace(name, i).second)
        Fail(path, lineNumber, "duplicate CSV column: " + name);
    }
    for (const auto &name : required)
    {
      if (index.find(name) == index.end())
        Fail(path, lineNumber, "missing required CSV column: " + name);
    }
    return index;
  }

  const std::string &Field(const std::filesystem::path &path,
                           std::size_t lineNumber,
                           const std::vector<std::string> &fields,
                           const std::unordered_map<std::string, std::size_t> &columns,
                           const std::string &name)
  {
    const std::size_t index = columns.at(name);
    if (index >= fields.size())
      Fail(path, lineNumber, "missing value for column: " + name);
    return fields[index];
  }

  int ParseInt(const std::filesystem::path &path, std::size_t lineNumber,
               const std::string &column, const std::string &value)
  {
    try
    {
      std::size_t used = 0;
      const long parsed = std::stol(value, &used);
      if (used != value.size() || parsed < std::numeric_limits<int>::min() ||
          parsed > std::numeric_limits<int>::max())
        throw std::invalid_argument("integer range");
      return static_cast<int>(parsed);
    }
    catch (...)
    {
      Fail(path, lineNumber, "invalid integer in " + column + ": " + value);
    }
  }

  std::uint64_t ParseUInt64(const std::filesystem::path &path,
                            const std::string &key,
                            const std::string &value)
  {
    try
    {
      std::size_t used = 0;
      const auto parsed = std::stoull(value, &used);
      if (used != value.size())
        throw std::invalid_argument("uint64 suffix");
      return parsed;
    }
    catch (...)
    {
      Fail(path, 0, "invalid unsigned integer metadata " + key + '=' + value);
    }
  }

  double ParseDouble(const std::filesystem::path &path, std::size_t lineNumber,
                     const std::string &name, const std::string &value)
  {
    try
    {
      std::size_t used = 0;
      const double parsed = std::stod(value, &used);
      if (used != value.size() || !std::isfinite(parsed))
        throw std::invalid_argument("non-finite or suffix");
      return parsed;
    }
    catch (...)
    {
      Fail(path, lineNumber, "invalid floating-point value in " + name + ": " + value);
    }
  }

  bool ParseBool(const std::filesystem::path &path, std::size_t lineNumber,
                 const std::string &name, std::string value)
  {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (value == "1" || value == "true")
      return true;
    if (value == "0" || value == "false")
      return false;
    Fail(path, lineNumber, "invalid boolean value in " + name + ": " + value);
  }

  bool NearlyEqual(double actual, double expected, double relativeTolerance)
  {
    return std::abs(actual - expected) <=
           relativeTolerance * std::max({1.0, std::abs(actual), std::abs(expected)});
  }

  const std::string &RequireMetadata(const std::filesystem::path &path,
                                     const RveMetadata &metadata,
                                     const std::string &key)
  {
    const auto found = metadata.values.find(key);
    if (found == metadata.values.end() || found->second.empty())
      Fail(path, 0, "missing required metadata: " + key);
    return found->second;
  }

  RvePhase ValidatePhase(const std::filesystem::path &path,
                         std::size_t lineNumber,
                         int phaseId,
                         const std::string &phaseName)
  {
    if (phaseId == 0 && phaseName == "BN")
      return RvePhase::BN;
    if (phaseId == 1 && phaseName == "ZnS")
      return RvePhase::ZnS;
    Fail(path, lineNumber,
         "phase/phase_id mismatch: phase=" + phaseName +
             ", phase_id=" + std::to_string(phaseId));
  }

  template <typename Row>
  void ValidateRadiusClass(const std::filesystem::path &path,
                           std::size_t lineNumber,
                           const Row &row,
                           std::unordered_map<int, RveRadiusClass> &classes)
  {
    const auto found = classes.find(row.radiusClassId);
    if (found == classes.end())
    {
      RveRadiusClass radiusClass;
      radiusClass.radiusClassId = row.radiusClassId;
      radiusClass.phase = row.phase;
      radiusClass.phaseId = row.phaseId;
      radiusClass.radiusClass = row.radiusClass;
      radiusClass.radiusUm = row.radiusUm;
      classes.emplace(row.radiusClassId, radiusClass);
      return;
    }
    const auto &expected = found->second;
    if (expected.phaseId != row.phaseId || expected.phase != row.phase ||
        expected.radiusClass != row.radiusClass ||
        !NearlyEqual(expected.radiusUm, row.radiusUm, kRadiusToleranceUm))
    {
      Fail(path, lineNumber,
           "inconsistent radius_class_id=" + std::to_string(row.radiusClassId));
    }
  }

  void ParseMetadataLine(const std::filesystem::path &path,
                         std::size_t lineNumber,
                         const std::string &line,
                         RveMetadata &metadata)
  {
    const std::string content = Trim(line.substr(1));
    const auto equals = content.find('=');
    if (equals == std::string::npos)
      return;
    const std::string key = Trim(content.substr(0, equals));
    const std::string value = Trim(content.substr(equals + 1));
    if (key.empty())
      Fail(path, lineNumber, "empty metadata key");
    metadata.values[key] = value;
  }

  void PopulateAndValidateMetadata(const std::filesystem::path &path,
                                   RveMetadata &metadata)
  {
    metadata.formatVersion = ParseInt(path, 0, "format_version",
                                      RequireMetadata(path, metadata, "format_version"));
    metadata.generator = RequireMetadata(path, metadata, "generator");
    metadata.boundary = RequireMetadata(path, metadata, "boundary");
    metadata.coordinateConvention = RequireMetadata(path, metadata, "coordinate_convention");
    metadata.oneRowPerUniquePeriodicParticle = ParseBool(
        path, 0, "one_row_per_unique_periodic_particle",
        RequireMetadata(path, metadata, "one_row_per_unique_periodic_particle"));
    metadata.geant4MaterialKey = RequireMetadata(path, metadata, "geant4_material_key");
    metadata.geant4SolidCacheKey = RequireMetadata(path, metadata, "geant4_solid_cache_key");
    metadata.boxXUm = ParseDouble(path, 0, "box_x_um", RequireMetadata(path, metadata, "box_x_um"));
    metadata.boxYUm = ParseDouble(path, 0, "box_y_um", RequireMetadata(path, metadata, "box_y_um"));
    metadata.boxZUm = ParseDouble(path, 0, "box_z_um", RequireMetadata(path, metadata, "box_z_um"));
    metadata.phiTarget = ParseDouble(path, 0, "phi_target", RequireMetadata(path, metadata, "phi_target"));
    metadata.phiAchieved = ParseDouble(path, 0, "phi_achieved", RequireMetadata(path, metadata, "phi_achieved"));
    metadata.overlapGapUm = ParseDouble(path, 0, "overlap_gap_um", RequireMetadata(path, metadata, "overlap_gap_um"));
    metadata.finalMaxOverlapViolationUm = ParseDouble(
        path, 0, "final_max_overlap_violation_um",
        RequireMetadata(path, metadata, "final_max_overlap_violation_um"));
    metadata.seed = ParseUInt64(path, "seed", RequireMetadata(path, metadata, "seed"));
    const auto radiusStats = metadata.values.find("radius_stats_csv");
    if (radiusStats != metadata.values.end())
      metadata.radiusStatsCsv = radiusStats->second;
    metadata.periodicImagesCsv = RequireMetadata(path, metadata, "periodic_images_csv");

    if (metadata.formatVersion != 3)
      Fail(path, 0, "unsupported format_version=" + std::to_string(metadata.formatVersion));
    if (metadata.boundary != "periodic_xyz")
      Fail(path, 0, "format_version=3 requires boundary=periodic_xyz");
    if (metadata.coordinateConvention != "centered_primary_cell[-L/2,L/2)")
      Fail(path, 0, "unsupported coordinate_convention=" + metadata.coordinateConvention);
    if (!metadata.oneRowPerUniquePeriodicParticle)
      Fail(path, 0, "one_row_per_unique_periodic_particle must be true");
    if (metadata.geant4MaterialKey != "phase_id")
      Fail(path, 0, "geant4_material_key must be phase_id");
    if (metadata.geant4SolidCacheKey != "radius_class_id")
      Fail(path, 0, "geant4_solid_cache_key must be radius_class_id");
    if (metadata.boxXUm <= 0.0 || metadata.boxYUm <= 0.0 || metadata.boxZUm <= 0.0)
      Fail(path, 0, "box_x_um, box_y_um, and box_z_um must be positive");
    if (metadata.overlapGapUm < 0.0 || metadata.finalMaxOverlapViolationUm < 0.0)
      Fail(path, 0, "overlap tolerances must be non-negative");
    if (RequireMetadata(path, metadata, "phase_0_name") != "BN" ||
        RequireMetadata(path, metadata, "phase_1_name") != "ZnS")
      Fail(path, 0, "phase_0_name/phase_1_name must be BN/ZnS");
    const double bnDensity = ParseDouble(
        path, 0, "phase_0_density_g_cm3",
        RequireMetadata(path, metadata, "phase_0_density_g_cm3"));
    const double znsDensity = ParseDouble(
        path, 0, "phase_1_density_g_cm3",
        RequireMetadata(path, metadata, "phase_1_density_g_cm3"));
    if (!NearlyEqual(bnDensity, kBnDensityGPerCm3, kMetadataRelativeTolerance) ||
        !NearlyEqual(znsDensity, kZnSDensityGPerCm3, kMetadataRelativeTolerance))
      Fail(path, 0, "phase densities do not match the Geant4 BN/ZnS materials");
  }

  void ValidatePrimaryCoordinates(const std::filesystem::path &path,
                                  std::size_t lineNumber,
                                  const RveParticle &particle,
                                  const RveMetadata &metadata)
  {
    const auto inside = [](double coordinateUm, double lengthUm)
    {
      return coordinateUm >= -0.5 * lengthUm && coordinateUm < +0.5 * lengthUm;
    };
    if (!inside(particle.center.x() / um, metadata.boxXUm) ||
        !inside(particle.center.y() / um, metadata.boxYUm) ||
        !inside(particle.center.z() / um, metadata.boxZUm))
    {
      Fail(path, lineNumber, "primary particle center lies outside [-L/2,L/2)");
    }
  }

  void ValidateMetadataCountsAndTotals(RvePlacementData &data)
  {
    const auto &path = data.placementPath;
    const auto &metadata = data.metadata;
    const auto requireCount = [&](const std::string &key, std::size_t actual)
    {
      const int expected = ParseInt(path, 0, key, RequireMetadata(path, metadata, key));
      if (expected < 0 || static_cast<std::size_t>(expected) != actual)
        Fail(path, 0, key + " does not match parsed rows");
    };
    requireCount("phase_0_count", data.bnCount);
    requireCount("phase_1_count", data.znsCount);

    const double expectedBnVolume = ParseDouble(
        path, 0, "phase_0_volume_um3", RequireMetadata(path, metadata, "phase_0_volume_um3"));
    const double expectedZnSVolume = ParseDouble(
        path, 0, "phase_1_volume_um3", RequireMetadata(path, metadata, "phase_1_volume_um3"));
    if (!NearlyEqual(data.bnVolumeUm3, expectedBnVolume, kMetadataRelativeTolerance) ||
        !NearlyEqual(data.znsVolumeUm3, expectedZnSVolume, kMetadataRelativeTolerance))
      Fail(path, 0, "recomputed phase volumes do not match metadata");

    const double boxVolume = metadata.boxXUm * metadata.boxYUm * metadata.boxZUm;
    data.phiAchieved = (data.bnVolumeUm3 + data.znsVolumeUm3) / boxVolume;
    if (!NearlyEqual(data.phiAchieved, metadata.phiAchieved, kMetadataRelativeTolerance))
      Fail(path, 0, "recomputed phi_achieved does not match metadata");

    data.znsToBnMassRatio =
        (data.znsVolumeUm3 * kZnSDensityGPerCm3) /
        (data.bnVolumeUm3 * kBnDensityGPerCm3);
    const std::string ratio = RequireMetadata(path, metadata, "actual_mass_ratio_BN_to_ZnS");
    const auto colon = ratio.find(':');
    if (colon == std::string::npos)
      Fail(path, 0, "invalid actual_mass_ratio_BN_to_ZnS=" + ratio);
    const double lhs = ParseDouble(path, 0, "actual_mass_ratio_BN_to_ZnS", ratio.substr(0, colon));
    const double rhs = ParseDouble(path, 0, "actual_mass_ratio_BN_to_ZnS", ratio.substr(colon + 1));
    if (lhs <= 0.0 || rhs <= 0.0 ||
        !NearlyEqual(data.znsToBnMassRatio, rhs / lhs, kMetadataRelativeTolerance))
      Fail(path, 0, "recomputed BN:ZnS mass ratio does not match metadata");
  }

  void ReadPeriodicImages(RvePlacementData &data,
                          const std::filesystem::path &path)
  {
    std::ifstream input(path);
    if (!input)
      Fail(path, 0, "failed to open periodic images CSV");

    const std::vector<std::string> required = {
        "copy_id", "source_particle_id", "is_primary", "shift_ix", "shift_iy", "shift_iz",
        "phase", "phase_id", "radius_bin_id", "radius_class_id", "radius_class",
        "radius_um", "x_um", "y_um", "z_um"};
    std::unordered_map<std::string, std::size_t> columns;
    std::string line;
    std::size_t lineNumber = 0;
    bool haveHeader = false;
    std::unordered_map<int, int> primaryCountByParticle;
    while (std::getline(input, line))
    {
      ++lineNumber;
      line = Trim(line);
      if (line.empty() || line[0] == '#')
        continue;
      const auto fields = SplitCsv(line);
      if (!haveHeader)
      {
        columns = BuildColumnIndex(path, lineNumber, fields, required);
        haveHeader = true;
        continue;
      }

      RveGeometryCopy copy;
      copy.copyId = ParseInt(path, lineNumber, "copy_id", Field(path, lineNumber, fields, columns, "copy_id"));
      copy.sourceParticleId = ParseInt(path, lineNumber, "source_particle_id", Field(path, lineNumber, fields, columns, "source_particle_id"));
      copy.isPrimary = ParseBool(path, lineNumber, "is_primary", Field(path, lineNumber, fields, columns, "is_primary"));
      copy.shiftIx = ParseInt(path, lineNumber, "shift_ix", Field(path, lineNumber, fields, columns, "shift_ix"));
      copy.shiftIy = ParseInt(path, lineNumber, "shift_iy", Field(path, lineNumber, fields, columns, "shift_iy"));
      copy.shiftIz = ParseInt(path, lineNumber, "shift_iz", Field(path, lineNumber, fields, columns, "shift_iz"));
      const std::string phaseName = Field(path, lineNumber, fields, columns, "phase");
      copy.phaseId = ParseInt(path, lineNumber, "phase_id", Field(path, lineNumber, fields, columns, "phase_id"));
      copy.phase = ValidatePhase(path, lineNumber, copy.phaseId, phaseName);
      copy.radiusBinId = ParseInt(path, lineNumber, "radius_bin_id", Field(path, lineNumber, fields, columns, "radius_bin_id"));
      copy.radiusClassId = ParseInt(path, lineNumber, "radius_class_id", Field(path, lineNumber, fields, columns, "radius_class_id"));
      copy.radiusClass = Field(path, lineNumber, fields, columns, "radius_class");
      copy.radiusUm = ParseDouble(path, lineNumber, "radius_um", Field(path, lineNumber, fields, columns, "radius_um"));
      const double xUm = ParseDouble(path, lineNumber, "x_um", Field(path, lineNumber, fields, columns, "x_um"));
      const double yUm = ParseDouble(path, lineNumber, "y_um", Field(path, lineNumber, fields, columns, "y_um"));
      const double zUm = ParseDouble(path, lineNumber, "z_um", Field(path, lineNumber, fields, columns, "z_um"));
      copy.center = G4ThreeVector(xUm * um, yUm * um, zUm * um);
      if (copy.copyId < 0 || copy.sourceParticleId < 0 || copy.radiusClassId < 0 || copy.radiusUm <= 0.0)
        Fail(path, lineNumber, "copy IDs and radius must be non-negative/positive");
      if (!data.copyIndexById.emplace(copy.copyId, data.geometryCopies.size()).second)
        Fail(path, lineNumber, "duplicate copy_id=" + std::to_string(copy.copyId));
      const auto sourceFound = data.particleIndexById.find(copy.sourceParticleId);
      if (sourceFound == data.particleIndexById.end())
        Fail(path, lineNumber, "unknown source_particle_id=" + std::to_string(copy.sourceParticleId));
      const auto &source = data.particles[sourceFound->second];
      if (copy.phaseId != source.phaseId || copy.radiusBinId != source.radiusBinId ||
          copy.radiusClassId != source.radiusClassId || copy.radiusClass != source.radiusClass ||
          !NearlyEqual(copy.radiusUm, source.radiusUm, kRadiusToleranceUm))
        Fail(path, lineNumber, "PBC copy attributes do not match source particle");
      const G4ThreeVector expected = source.center +
                                     G4ThreeVector(copy.shiftIx * data.metadata.boxXUm * um,
                                                   copy.shiftIy * data.metadata.boxYUm * um,
                                                   copy.shiftIz * data.metadata.boxZUm * um);
      if ((copy.center - expected).mag() / um > 5.0e-9)
        Fail(path, lineNumber, "PBC copy coordinate does not match source plus periodic shift");
      if (copy.isPrimary != (copy.shiftIx == 0 && copy.shiftIy == 0 && copy.shiftIz == 0))
        Fail(path, lineNumber, "is_primary is inconsistent with shift indices");
      if (copy.isPrimary)
        ++primaryCountByParticle[copy.sourceParticleId];
      ValidateRadiusClass(path, lineNumber, copy, data.radiusClasses);
      data.geometryCopies.push_back(copy);
    }
    if (!haveHeader || data.geometryCopies.empty())
      Fail(path, 0, "periodic images CSV has no data rows");
    if (primaryCountByParticle.size() != data.particles.size())
      Fail(path, 0, "periodic images CSV does not contain every primary particle exactly once");
    for (const auto &entry : primaryCountByParticle)
    {
      if (entry.second != 1)
        Fail(path, 0, "source_particle_id=" + std::to_string(entry.first) +
                          " has multiple primary copies");
    }
  }
}

bool RvePlacementReader::IsFormatVersion3(const std::filesystem::path &path)
{
  std::ifstream input(path);
  if (!input)
    return false;
  std::string line;
  while (std::getline(input, line))
  {
    line = Trim(line);
    if (line.empty())
      continue;
    if (line == "# format_version=3" || line == "#format_version=3")
      return true;
    if (line[0] != '#')
      return false;
  }
  return false;
}

bool RvePlacementReader::IsMainPlacementCandidate(const std::filesystem::path &path)
{
  if (path.extension() != ".csv")
    return false;
  const std::string name = path.filename().string();
  if (name.size() >= 15 && name.rfind("_pbc_images.csv") == name.size() - 15)
    return false;
  if (name.size() >= 17 && name.rfind("_radius_stats.csv") == name.size() - 17)
    return false;
  return true;
}

RvePlacementData RvePlacementReader::ReadVersion3(
    const std::filesystem::path &placementPath,
    const std::filesystem::path &periodicImagesOverride,
    bool readPeriodicImages)
{
  RvePlacementData data;
  std::error_code canonicalError;
  data.placementPath = std::filesystem::weakly_canonical(placementPath, canonicalError);
  if (canonicalError)
    data.placementPath = placementPath.lexically_normal();
  std::ifstream input(data.placementPath);
  if (!input)
    Fail(data.placementPath, 0, "failed to open placement CSV");

  const std::vector<std::string> required = {
      "particle_id", "phase", "phase_id", "radius_bin_id", "radius_class_id",
      "radius_class", "radius_um", "x_um", "y_um", "z_um"};
  std::unordered_map<std::string, std::size_t> columns;
  std::string line;
  std::size_t lineNumber = 0;
  bool haveHeader = false;
  while (std::getline(input, line))
  {
    ++lineNumber;
    line = Trim(line);
    if (line.empty())
      continue;
    if (line[0] == '#')
    {
      if (haveHeader)
        Fail(data.placementPath, lineNumber, "metadata found after CSV header");
      ParseMetadataLine(data.placementPath, lineNumber, line, data.metadata);
      continue;
    }
    const auto fields = SplitCsv(line);
    if (!haveHeader)
    {
      PopulateAndValidateMetadata(data.placementPath, data.metadata);
      columns = BuildColumnIndex(data.placementPath, lineNumber, fields, required);
      haveHeader = true;
      continue;
    }

    RveParticle particle;
    particle.particleId = ParseInt(data.placementPath, lineNumber, "particle_id",
                                   Field(data.placementPath, lineNumber, fields, columns, "particle_id"));
    const std::string phaseName = Field(data.placementPath, lineNumber, fields, columns, "phase");
    particle.phaseId = ParseInt(data.placementPath, lineNumber, "phase_id",
                                Field(data.placementPath, lineNumber, fields, columns, "phase_id"));
    particle.phase = ValidatePhase(data.placementPath, lineNumber, particle.phaseId, phaseName);
    particle.radiusBinId = ParseInt(data.placementPath, lineNumber, "radius_bin_id",
                                    Field(data.placementPath, lineNumber, fields, columns, "radius_bin_id"));
    particle.radiusClassId = ParseInt(data.placementPath, lineNumber, "radius_class_id",
                                      Field(data.placementPath, lineNumber, fields, columns, "radius_class_id"));
    particle.radiusClass = Field(data.placementPath, lineNumber, fields, columns, "radius_class");
    particle.radiusUm = ParseDouble(data.placementPath, lineNumber, "radius_um",
                                    Field(data.placementPath, lineNumber, fields, columns, "radius_um"));
    const double xUm = ParseDouble(data.placementPath, lineNumber, "x_um",
                                   Field(data.placementPath, lineNumber, fields, columns, "x_um"));
    const double yUm = ParseDouble(data.placementPath, lineNumber, "y_um",
                                   Field(data.placementPath, lineNumber, fields, columns, "y_um"));
    const double zUm = ParseDouble(data.placementPath, lineNumber, "z_um",
                                   Field(data.placementPath, lineNumber, fields, columns, "z_um"));
    particle.center = G4ThreeVector(xUm * um, yUm * um, zUm * um);
    if (particle.particleId < 0 || particle.radiusBinId < 0 ||
        particle.radiusClassId < 0 || particle.radiusUm <= 0.0 || particle.radiusClass.empty())
      Fail(data.placementPath, lineNumber, "particle IDs and radius must be non-negative/positive");
    if (!data.particleIndexById.emplace(particle.particleId, data.particles.size()).second)
      Fail(data.placementPath, lineNumber, "duplicate particle_id=" + std::to_string(particle.particleId));
    ValidatePrimaryCoordinates(data.placementPath, lineNumber, particle, data.metadata);
    ValidateRadiusClass(data.placementPath, lineNumber, particle, data.radiusClasses);
    const double sphereVolume = (4.0 / 3.0) * CLHEP::pi *
                                particle.radiusUm * particle.radiusUm * particle.radiusUm;
    if (particle.phase == RvePhase::BN)
    {
      ++data.bnCount;
      data.bnVolumeUm3 += sphereVolume;
    }
    else
    {
      ++data.znsCount;
      data.znsVolumeUm3 += sphereVolume;
    }
    data.particles.push_back(particle);
  }
  if (!haveHeader || data.particles.empty())
    Fail(data.placementPath, 0, "placement CSV has no data rows");
  ValidateMetadataCountsAndTotals(data);

  if (!readPeriodicImages)
    return data;

  data.periodicImagesPath = periodicImagesOverride.empty()
                                ? data.placementPath.parent_path() / data.metadata.periodicImagesCsv
                                : periodicImagesOverride;
  canonicalError.clear();
  const auto canonicalPbc = std::filesystem::weakly_canonical(data.periodicImagesPath, canonicalError);
  if (!canonicalError)
    data.periodicImagesPath = canonicalPbc;
  ReadPeriodicImages(data, data.periodicImagesPath);
  return data;
}
