#include "DetectorConstruction.hh"
#include "AnalysisConfig.hh"
#include "RvePlacementReader.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4MaterialPropertiesTable.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
  constexpr G4double kPlacementOverlapTolerance = 5.0e-5 * um;

  struct PlacementData
  {
    std::vector<G4ThreeVector> bnCenters;
    std::vector<G4ThreeVector> znsCenters;
    std::string seedBase;
    G4double patchXYUm = -1.0;
    G4double localThicknessUm = -1.0;
    G4double bnWt = -1.0;
    G4double znsWt = -1.0;
    G4double bnRadiusUm = -1.0;
    G4double znsRadiusUm = -1.0;
    G4double overlapGapUm = -1.0;
  };

  // --------------------------------------------------------------------
  // Overall composition parts for the whole micro-patch
  // Keep consistent with the homogeneous model:
  //   solid  = BN + ZnS particles
  //   binder = binder in matrix background
  //   air    = void in matrix background
  // Current background split: binder:air = 21.6:14.4 -> void fraction = 0.40
  // --------------------------------------------------------------------
  G4double gSolidPart = 64.0;
  G4double gBinderPart = 21.6;
  G4double gAirPart = 14.4;
  constexpr G4double kBinderPlexiglassRIndex = 1.49;
  constexpr G4double kAirRIndex = 1.0;

  inline G4double ComputeMatrixVoidFractionFromOverallParts()
  {
    const G4double bgPart = gBinderPart + gAirPart;
    if (bgPart <= 0.0)
      return 0.0;

    return gAirPart / bgPart; // 14.4 / (21.6 + 14.4) = 0.40
  }

  inline G4double LorentzLorenzTerm(G4double rindex)
  {
    const G4double n2 = rindex * rindex;
    return (n2 - 1.0) / (n2 + 2.0);
  }

  inline G4double RIndexFromLorentzLorenzTerm(G4double term)
  {
    const G4double clamped = std::clamp(term, -0.499999, 0.999999);
    const G4double n2 = (1.0 + 2.0 * clamped) / (1.0 - clamped);
    return (n2 > 0.0) ? std::sqrt(n2) : 1.0;
  }

  G4double ComputeMatrixRIndexLorentzLorenz()
  {
    const G4double bgPart = gBinderPart + gAirPart;
    if (bgPart <= 0.0)
      return kBinderPlexiglassRIndex;

    const G4double binderFrac = gBinderPart / bgPart;
    const G4double airFrac = gAirPart / bgPart;
    const G4double mixedTerm =
        binderFrac * LorentzLorenzTerm(kBinderPlexiglassRIndex) +
        airFrac * LorentzLorenzTerm(kAirRIndex);
    return RIndexFromLorentzLorenzTerm(mixedTerm);
  }

  std::string Trim(const std::string &s)
  {
    std::size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b])))
      ++b;

    std::size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1])))
      --e;

    return s.substr(b, e - b);
  }

  bool IsTruthyEnvValue(const char *value)
  {
    if (value == nullptr)
      return false;

    std::string lower = Trim(value);
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c)
                   { return static_cast<char>(std::tolower(c)); });

    return lower == "1" ||
           lower == "true" ||
           lower == "yes" ||
           lower == "on" ||
           lower == "y";
  }

  bool IsEnvFlagEnabled(const char *name)
  {
    return IsTruthyEnvValue(std::getenv(name));
  }

  std::vector<std::string> SplitFlexible(const std::string &line)
  {
    std::string s = line;
    for (char &c : s)
    {
      if (c == ',' || c == '\t')
        c = ' ';
    }

    std::vector<std::string> tokens;
    std::stringstream ss(s);
    std::string token;
    while (ss >> token)
    {
      tokens.push_back(token);
    }
    return tokens;
  }

  std::string ParseMetadataValue(const std::string &line, const std::string &key)
  {
    std::string s = Trim(line);
    if (s.empty() || s[0] != '#')
      return "";

    s = Trim(s.substr(1));
    const std::string prefix = key + "=";
    if (s.rfind(prefix, 0) == 0)
    {
      return Trim(s.substr(prefix.size()));
    }

    const auto tokens = SplitFlexible(s);
    if (tokens.size() >= 2 && tokens[0] == key)
    {
      return tokens[1];
    }

    return "";
  }

  void TrySetPlacementMetadataDouble(const std::string &line,
                                     const std::string &key,
                                     G4double &target)
  {
    const std::string value = ParseMetadataValue(line, key);
    if (value.empty())
      return;

    try
    {
      target = std::stod(value);
    }
    catch (...)
    {
      G4Exception("DetectorConstruction::Construct",
                  "BNZS205", FatalException,
                  ("Invalid placement metadata value: " + key + "=" + value).c_str());
    }
  }

  std::filesystem::path MakePlacementRootDirectoryPath()
  {
    return AnalysisConfig::ProjectRootPath() / "Input" / "placements";
  }

  std::filesystem::path MakeVersion3PlacementRootDirectoryPath()
  {
    return AnalysisConfig::ProjectRootPath() / "Input" / "output_pbc";
  }

  std::string WeightPartToTagString(double v)
  {
    const double rounded = std::round(v);
    if (std::fabs(v - rounded) < 1.0e-9)
    {
      std::ostringstream oss;
      oss << static_cast<long long>(rounded);
      return oss.str();
    }

    std::ostringstream oss;
    oss << v;
    return oss.str();
  }

  std::string MakeRatioFolderName(G4double bnWt, G4double znsWt)
  {
    if (bnWt <= 0.0 || znsWt <= 0.0)
    {
      G4Exception("DetectorConstruction::Construct",
                  "BNZS201", FatalException,
                  "Invalid BN/ZnS weight ratio.");
      return "";
    }

    return WeightPartToTagString(bnWt) + "-" +
           WeightPartToTagString(znsWt);
  }

  std::filesystem::path ResolvePlacementPathForCurrentProject(const std::string &rawPath)
  {
    namespace fs = std::filesystem;

    const fs::path input(rawPath);
    std::error_code ec;

    if (input.is_absolute() && fs::exists(input, ec))
      return fs::weakly_canonical(input, ec);

    if (fs::exists(input, ec))
      return fs::weakly_canonical(input, ec);

    const fs::path placementsRoot = MakePlacementRootDirectoryPath();
    const fs::path version3Root = MakeVersion3PlacementRootDirectoryPath();
    std::vector<fs::path> candidates;
    const auto returnExistingCandidate =
        [](const std::vector<fs::path> &paths) -> fs::path
    {
      for (const auto &candidate : paths)
      {
        std::error_code existsEc;
        if (fs::exists(candidate, existsEc) && !existsEc)
          return fs::weakly_canonical(candidate, existsEc);
      }
      return {};
    };

    const auto findUniqueRecursiveMatch =
        [&](const fs::path &searchRoot,
            const fs::path &fileName) -> fs::path
    {
      std::error_code searchEc;
      if (fileName.empty() ||
          !fs::exists(searchRoot, searchEc) ||
          !fs::is_directory(searchRoot, searchEc))
      {
        return {};
      }

      std::vector<fs::path> matches;
      for (const auto &entry : fs::recursive_directory_iterator(searchRoot, searchEc))
      {
        if (searchEc)
          break;
        if (!entry.is_regular_file())
          continue;
        if (entry.path().filename() == fileName)
          matches.push_back(entry.path());
      }

      if (matches.empty())
        return {};

      std::sort(matches.begin(), matches.end());
      if (matches.size() > 1)
      {
        G4Exception("DetectorConstruction::Construct",
                    "BNZS202", FatalException,
                    ("Ambiguous placement filename under " + searchRoot.string() + ": "
                     + fileName.string())
                        .c_str());
      }

      return fs::weakly_canonical(matches.front(), searchEc);
    };

    const auto parts = input.lexically_normal();
    if (!parts.empty())
    {
      std::vector<std::string> strings;
      strings.reserve(std::distance(parts.begin(), parts.end()));
      for (const auto &part : parts)
      {
        strings.push_back(part.string());
      }

      for (std::size_t i = 0; i + 1 < strings.size(); ++i)
      {
        if (strings[i] == "Input" && strings[i + 1] == "placements")
        {
          fs::path tail;
          for (std::size_t j = i + 2; j < strings.size(); ++j)
            tail /= strings[j];
          if (!tail.empty())
            candidates.push_back(placementsRoot / tail);
          break;
        }
        if (strings[i] == "Input" && strings[i + 1] == "output_pbc")
        {
          fs::path tail;
          for (std::size_t j = i + 2; j < strings.size(); ++j)
            tail /= strings[j];
          if (!tail.empty())
            candidates.push_back(version3Root / tail);
          break;
        }
      }
    }

    // Prefer explicit project-relative placement paths such as
    // Prefer an exact ratio directory match over filename-only fallback matching.
    if (const fs::path directMatch = returnExistingCandidate(candidates);
        !directMatch.empty())
    {
      return directMatch;
    }

    if (!input.parent_path().filename().empty())
    {
      candidates.push_back(version3Root / input.parent_path().filename() / input.filename());
      candidates.push_back(placementsRoot / input.parent_path().filename() / input.filename());

      if (const fs::path parentScopedMatch = returnExistingCandidate(candidates);
          !parentScopedMatch.empty())
      {
        return parentScopedMatch;
      }

      const fs::path nestedMatch =
          findUniqueRecursiveMatch(version3Root / input.parent_path().filename(), input.filename());
      if (!nestedMatch.empty())
        return nestedMatch;

      const fs::path legacyNestedMatch =
          findUniqueRecursiveMatch(placementsRoot / input.parent_path().filename(), input.filename());
      if (!legacyNestedMatch.empty())
        return legacyNestedMatch;
    }

    if (!input.filename().empty())
    {
      const fs::path version3Match = findUniqueRecursiveMatch(version3Root, input.filename());
      if (!version3Match.empty())
        return version3Match;
      for (const auto &entry : fs::directory_iterator(placementsRoot, ec))
      {
        if (ec || !entry.is_directory())
          continue;
        const fs::path nestedMatch = findUniqueRecursiveMatch(entry.path(), input.filename());
        if (!nestedMatch.empty())
          return nestedMatch;
      }
    }

    return input;
  }

  std::string PickRandomPlacementFilePath(G4double bnWt, G4double znsWt)
  {
    namespace fs = std::filesystem;

    const std::string ratioFolderName = MakeRatioFolderName(bnWt, znsWt);
    const fs::path version3RatioDir = MakeVersion3PlacementRootDirectoryPath() / ratioFolderName;
    const fs::path legacyRatioDir = MakePlacementRootDirectoryPath() / ratioFolderName;
    fs::path ratioDir;
    if (fs::exists(version3RatioDir) && fs::is_directory(version3RatioDir))
      ratioDir = version3RatioDir;
    else if (fs::exists(legacyRatioDir) && fs::is_directory(legacyRatioDir))
      ratioDir = legacyRatioDir;

    if (ratioDir.empty())
    {
      std::ostringstream oss;
      oss << "No placement folder for BN:ZnS wt = "
          << bnWt << ":" << znsWt
          << ", checked: " << version3RatioDir.string()
          << " and " << legacyRatioDir.string();
      G4Exception("DetectorConstruction::Construct",
                  "BNZS203", FatalException,
                  oss.str().c_str());
      return "";
    }

    std::vector<fs::path> files;
    const bool selectingVersion3 = (ratioDir == version3RatioDir);
    for (const auto &entry : fs::recursive_directory_iterator(ratioDir))
    {
      if (!entry.is_regular_file())
        continue;

      const std::string ext = entry.path().extension().string();
      const bool legacyCandidate = !selectingVersion3 && ext == ".txt";
      const bool csvCandidate = RvePlacementReader::IsMainPlacementCandidate(entry.path()) &&
                                (!selectingVersion3 ||
                                 RvePlacementReader::IsFormatVersion3(entry.path()));
      if (legacyCandidate || csvCandidate)
        files.push_back(entry.path());
    }

    if (files.empty())
    {
      G4Exception("DetectorConstruction::Construct",
                  "BNZS204", FatalException,
                  ("No placement CSV/TXT files found in: " + ratioDir.string()).c_str());
      return "";
    }

    std::sort(files.begin(), files.end());
    static thread_local std::mt19937 rng{std::random_device{}()};
    std::uniform_int_distribution<std::size_t> dist(0, files.size() - 1);

    return files[dist(rng)].string();
  }

  PlacementData ReadPlacementCSV(const std::string &filePath)
  {
    PlacementData data;

    const auto resolvedPath = ResolvePlacementPathForCurrentProject(filePath);
    std::ifstream fin(resolvedPath.string().c_str(), std::ios::in);
    if (!fin.is_open())
    {
      G4Exception("DetectorConstruction::Construct",
                  "BNZS205", FatalException,
                  ("Failed to open placement CSV: " + resolvedPath.string()).c_str());
      return data;
    }

    std::string line;
    bool headerHandled = false;

    while (std::getline(fin, line))
    {
      line = Trim(line);
      if (line.empty())
        continue;

      if (line[0] == '#')
      {
        const std::string seedBase = ParseMetadataValue(line, "seedBase");
        if (!seedBase.empty())
          data.seedBase = seedBase;
        TrySetPlacementMetadataDouble(line, "patchXY_um", data.patchXYUm);
        TrySetPlacementMetadataDouble(line, "localThickness_um", data.localThicknessUm);
        TrySetPlacementMetadataDouble(line, "bnWt", data.bnWt);
        TrySetPlacementMetadataDouble(line, "znsWt", data.znsWt);
        TrySetPlacementMetadataDouble(line, "bnRadius_um", data.bnRadiusUm);
        TrySetPlacementMetadataDouble(line, "znsRadius_um", data.znsRadiusUm);
        TrySetPlacementMetadataDouble(line, "overlapGap_um", data.overlapGapUm);
        continue;
      }

      if (!headerHandled)
      {
        headerHandled = true;

        std::string lower = line;
        for (char &ch : lower)
          ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));

        if (lower.find("type") != std::string::npos &&
            lower.find("x") != std::string::npos &&
            lower.find("y") != std::string::npos &&
            lower.find("z") != std::string::npos)
        {
          continue;
        }
      }

      const auto tokens = SplitFlexible(line);
      if (tokens.size() < 4)
        continue;

      const std::string type = Trim(tokens[0]);
      const std::string sx = Trim(tokens[1]);
      const std::string sy = Trim(tokens[2]);
      const std::string sz = Trim(tokens[3]);

      const G4double x = std::stod(sx) * um;
      const G4double y = std::stod(sy) * um;
      const G4double z = std::stod(sz) * um;

      G4ThreeVector pos(x, y, z);

      std::string typeLower = type;
      for (char &ch : typeLower)
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));

      if (typeLower == "bn")
      {
        data.bnCenters.push_back(pos);
      }
      else if (typeLower == "zns")
      {
        data.znsCenters.push_back(pos);
      }
      else
      {
        G4Exception("DetectorConstruction::Construct",
                    "BNZS206", FatalException,
                    ("Unknown particle type in placement CSV: " + type).c_str());
        return data;
      }
    }

    if (data.bnCenters.empty() && data.znsCenters.empty())
    {
      G4Exception("DetectorConstruction::Construct",
                  "BNZS207", FatalException,
                  ("Placement CSV is empty: " + resolvedPath.string()).c_str());
      return data;
    }

    return data;
  }

  G4bool IsFullyInsideBox(const G4ThreeVector &pos,
                          G4double radius,
                          G4double patchXY,
                          G4double localThickness)
  {
    return (std::abs(pos.x()) + radius <= 0.5 * patchXY &&
            std::abs(pos.y()) + radius <= 0.5 * patchXY &&
            std::abs(pos.z()) + radius <= 0.5 * localThickness);
  }

  void ValidatePlacementData(const PlacementData &data,
                             G4double patchXY,
                             G4double localThickness,
                             G4double bnRadius,
                             G4double znsRadius,
                             G4double overlapGap)
  {
    auto checkInside = [&](const G4ThreeVector &p, const char *name)
    {
      const G4double xMin = -0.5 * patchXY;
      const G4double xMax = 0.5 * patchXY;
      const G4double yMin = -0.5 * patchXY;
      const G4double yMax = 0.5 * patchXY;
      const G4double zMin = -0.5 * localThickness;
      const G4double zMax = 0.5 * localThickness;

      if (p.x() < xMin || p.x() > xMax ||
          p.y() < yMin || p.y() > yMax ||
          p.z() < zMin || p.z() > zMax)
      {
        std::ostringstream oss;
        oss << name << " center outside allowed region: ("
            << p.x() / um << ", " << p.y() / um << ", " << p.z() / um << ") um";
        G4Exception("DetectorConstruction::Construct",
                    "BNZS206", FatalException,
                    oss.str().c_str());
      }
    };

    for (const auto &p : data.bnCenters)
      checkInside(p, "BN");

    for (const auto &p : data.znsCenters)
      checkInside(p, "ZnS");

    // --- O(N) 网格空间索引代码 ---
    struct Part
    {
      G4ThreeVector pos;
      G4double radius;
    };
    std::vector<Part> allParts;
    allParts.reserve(data.bnCenters.size() + data.znsCenters.size());
    for (const auto &p : data.bnCenters)
      allParts.push_back({p, bnRadius});
    for (const auto &p : data.znsCenters)
      allParts.push_back({p, znsRadius});

    if (allParts.empty())
      return;

    G4int reducedGapPairCount = 0;
    G4double minPositiveClearance = std::numeric_limits<G4double>::infinity();
    G4double maxGapShortfall = 0.0;

    // 1. 构建网格
    const G4double maxRadius = std::max(bnRadius, znsRadius);
    const G4double cellSize = 2.0 * maxRadius + overlapGap;

    const int nx = std::max(1, static_cast<int>(std::ceil(patchXY / cellSize)));
    const int ny = std::max(1, static_cast<int>(std::ceil(patchXY / cellSize)));
    const int nz = std::max(1, static_cast<int>(std::ceil(localThickness / cellSize)));

    std::vector<std::vector<int>> grid(nx * ny * nz);

    auto getCellIndex = [&](const G4ThreeVector &p)
    {
      int cx = std::clamp(static_cast<int>((p.x() + 0.5 * patchXY) / cellSize), 0, nx - 1);
      int cy = std::clamp(static_cast<int>((p.y() + 0.5 * patchXY) / cellSize), 0, ny - 1);
      int cz = std::clamp(static_cast<int>((p.z() + 0.5 * localThickness) / cellSize), 0, nz - 1);
      return cx + cy * nx + cz * nx * ny;
    };

    // 2. 将球体分配到网格中
    for (std::size_t i = 0; i < allParts.size(); ++i)
    {
      grid[getCellIndex(allParts[i].pos)].push_back(static_cast<int>(i));
    }

    // 3. 仅检查相邻网格
    for (int cz = 0; cz < nz; ++cz)
    {
      for (int cy = 0; cy < ny; ++cy)
      {
        for (int cx = 0; cx < nx; ++cx)
        {
          const auto &bucket = grid[cx + cy * nx + cz * nx * ny];
          for (int i : bucket)
          {
            for (int dz = -1; dz <= 1; ++dz)
            {
              for (int dy = -1; dy <= 1; ++dy)
              {
                for (int dx = -1; dx <= 1; ++dx)
                {
                  int ncx = cx + dx, ncy = cy + dy, ncz = cz + dz;
                  if (ncx < 0 || ncx >= nx || ncy < 0 || ncy >= ny || ncz < 0 || ncz >= nz)
                    continue;

                  for (int j : grid[ncx + ncy * nx + ncz * nx * ny])
                  {
                    if (i >= j)
                      continue; // 避免重复比较和自己比自己

                    const G4double touchingDist =
                        allParts[i].radius + allParts[j].radius;
                    const G4double preferredDist = touchingDist + overlapGap;
                    const G4double dist2 =
                        (allParts[i].pos - allParts[j].pos).mag2();
                    const G4double hardMinDist =
                        std::max(0.0, touchingDist - kPlacementOverlapTolerance);

                    if (dist2 < hardMinDist * hardMinDist)
                    {
                      const G4double dist = std::sqrt(dist2);
                      const G4double clearance = dist - touchingDist;
                      std::ostringstream oss;
                      oss << "Overlap detected in placement CSV! "
                          << "surface clearance = " << clearance / um
                          << " um, but must be >= 0 within tolerance.";
                      G4Exception("DetectorConstruction::Construct",
                                  "BNZS207", FatalException,
                                  oss.str().c_str());
                    }

                    const G4double softMinDist =
                        std::max(0.0, preferredDist - kPlacementOverlapTolerance);
                    if (overlapGap > 0.0 && dist2 < softMinDist * softMinDist)
                    {
                      const G4double dist = std::sqrt(dist2);
                      const G4double clearance = dist - touchingDist;
                      if (clearance >= -kPlacementOverlapTolerance)
                      {
                        ++reducedGapPairCount;
                        minPositiveClearance = std::min(minPositiveClearance, clearance);
                        maxGapShortfall = std::max(maxGapShortfall, overlapGap - clearance);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    if (reducedGapPairCount > 0)
    {
      G4cerr << "[DetectorConstruction] Warning: placement contains "
             << reducedGapPairCount
             << " particle pairs whose surface clearance stays non-negative but is below the "
             << "configured overlapGap."
             << " min_clearance=" << minPositiveClearance / um
             << " um"
             << " overlapGap=" << overlapGap / um
             << " um"
             << " max_shortfall=" << maxGapShortfall / um
             << " um"
             << ". Continuing because only true overlaps are fatal."
             << G4endl;
    }
  }

  // --- 数学解析球冠体积函数 ---

  G4double CalculateClippedVolume(const G4ThreeVector &c,
                                  G4double r,
                                  G4double patchXY,
                                  G4double localThickness,
                                  G4int nz = 40,
                                  G4int nx = 40)
  {
    const G4double pi = CLHEP::pi;

    // 捷径 1：如果球体完全在箱子内部，直接返回完整体积
    if (std::abs(c.x()) + r <= 0.5 * patchXY &&
        std::abs(c.y()) + r <= 0.5 * patchXY &&
        std::abs(c.z()) + r <= 0.5 * localThickness)
    {
      return (4.0 / 3.0) * pi * r * r * r;
    }

    // 捷径 2：确定 Z 轴有效积分范围
    const G4double zMinBox = -0.5 * localThickness;
    const G4double zMaxBox = 0.5 * localThickness;
    const G4double z0 = std::max(zMinBox, c.z() - r);
    const G4double z1 = std::min(zMaxBox, c.z() + r);

    // 如果完全在箱子上方或下方
    if (z1 <= z0)
      return 0.0;

    G4double volume = 0.0;
    const G4double dz = (z1 - z0) / static_cast<G4double>(nz);

    // 第一层：沿 Z 轴切成圆盘
    for (G4int iz = 0; iz < nz; ++iz)
    {
      G4double z = z0 + (iz + 0.5) * dz;
      G4double rz2 = r * r - (z - c.z()) * (z - c.z());
      if (rz2 <= 0.0)
        continue; // 理论上不会发生，防御性判断

      G4double rho = std::sqrt(rz2); // 当前高度的二维圆半径

      // 确定 X 轴有效积分范围
      G4double x0 = std::max(-0.5 * patchXY, c.x() - rho);
      G4double x1 = std::min(0.5 * patchXY, c.x() + rho);

      if (x1 <= x0)
        continue; // 当前切片在这个高度上完全超出了 X 边界

      G4double dx = (x1 - x0) / static_cast<G4double>(nx);
      G4double area = 0.0;

      // 第二层：沿 X 轴切成条带
      for (G4int ix = 0; ix < nx; ++ix)
      {
        G4double x = x0 + (ix + 0.5) * dx;
        G4double ry2 = rho * rho - (x - c.x()) * (x - c.x());
        if (ry2 <= 0.0)
          continue;

        G4double ySpan = std::sqrt(ry2); // 当前条带在 Y 轴的半长

        // 解析求解 Y 轴落在箱内的实际长度 (极其优雅的处理了角落问题)
        G4double y0 = std::max(-0.5 * patchXY, c.y() - ySpan);
        G4double y1 = std::min(0.5 * patchXY, c.y() + ySpan);

        if (y1 > y0)
          area += (y1 - y0) * dx;
      }
      volume += area * dz;
    }

    return volume;
  }
}

// --------------------------------------------------------------------

DetectorConstruction::DetectorConstruction(AnalysisConfig *config)
    : G4VUserDetectorConstruction(),
      fScreenThickness(125.0 * um),
      fMicroThickness(30.0 * um),
      fPatchXY(50.0 * um),
      fBoxX(50.0 * um),
      fBoxY(50.0 * um),
      fBoxZ(30.0 * um),
      fBNRadius(2.5 * um),  // BN diameter = 5 um
      fZnSRadius(1.0 * um), // ZnS diameter = 2 um
      fBNPitch(5.0 * um),
      fZnSPitch(2.0 * um),
      fOverlapGap(0.01 * um),
      fSafeMarginXY(8.0 * um),
      fVoidVolumeFraction(ComputeMatrixVoidFractionFromOverallParts()),
      fBnWt(1.0),
      fZnsWt(2.0),
      fVacuumMaterial(nullptr),
      fMatrixMaterial(nullptr),
      fBNMaterial(nullptr),
      fZnSMaterial(nullptr),
      fScoringVolume(nullptr),
      fWorldLogical(nullptr),
      fMatrixLogical(nullptr),
      fBNLogical(nullptr),
      fZnSLogical(nullptr),
      fSpatialNx(1),
      fSpatialNy(1),
      fSpatialNz(1),
      fSpatialCellX(50.0 * um),
      fSpatialCellY(50.0 * um),
      fSpatialCellZ(30.0 * um),
      fConfig(config),
      fLoadedPlacementFile(""),
      fLoadedPlacementFileForRecord(""),
      fLoadedPeriodicImagesFile(""),
      fLoadedPeriodicImagesFileForRecord(""),
      fLoadedPlacementSeedBase(""),
      fPlacementFormatVersion(0),
      fPlacementSeed(0),
      fPlacementPhiAchieved(0.0),
      fPlacementZnSToBNMassRatio(0.0)
{
  if (fConfig != nullptr)
  {
    fPatchXY = fConfig->patchXY_um * um;
    fMicroThickness = fConfig->microThickness_um * um;
    fBoxX = fPatchXY;
    fBoxY = fPatchXY;
    fBoxZ = fMicroThickness;
    fBnWt = fConfig->bnWt;
    fZnsWt = fConfig->znsWt;

    G4cout << "[DetectorConstruction] AnalysisConfig attached:"
           << " mode=" << AnalysisConfig::RunModeName(fConfig->runMode)
           << " patchXY=" << fConfig->patchXY_um << " um"
           << " microThickness=" << fConfig->microThickness_um << " um"
           << " bnWt=" << fConfig->bnWt
           << " znsWt=" << fConfig->znsWt
           << G4endl;
  }
}

// --------------------------------------------------------------------

DetectorConstruction::~DetectorConstruction() = default;

DetectorConstruction::Phase DetectorConstruction::FindPhaseAtPointUm(
    G4double xUm,
    G4double yUm,
    G4double zUm) const
{
  return FindPhaseAtPoint(G4ThreeVector(xUm * um, yUm * um, zUm * um));
}

DetectorConstruction::Phase DetectorConstruction::FindPhaseAtPoint(
    const G4ThreeVector &point) const
{
  if (point.x() < -0.5 * fBoxX || point.x() >= 0.5 * fBoxX ||
      point.y() < -0.5 * fBoxY || point.y() >= 0.5 * fBoxY ||
      point.z() < -0.5 * fBoxZ || point.z() >= 0.5 * fBoxZ)
    return Phase::World;
  return FindPhaseInsidePrimary(point, fPlacementFormatVersion == 3);
}

DetectorConstruction::Phase DetectorConstruction::FindPeriodicPhaseAtPoint(
    const G4ThreeVector &point) const
{
  return FindPhaseInsidePrimary(WrapToPrimaryCell(point), true);
}

G4ThreeVector DetectorConstruction::WrapToPrimaryCell(const G4ThreeVector &point) const
{
  const auto wrap = [](G4double value, G4double length)
  {
    if (length <= 0.0)
      return value;
    G4double wrapped = value - length * std::floor((value + 0.5 * length) / length);
    if (wrapped >= 0.5 * length)
      wrapped -= length;
    return wrapped;
  };
  return G4ThreeVector(wrap(point.x(), fBoxX),
                       wrap(point.y(), fBoxY),
                       wrap(point.z(), fBoxZ));
}

G4ThreeVector DetectorConstruction::WrapToPrimaryCellInside(
    const G4ThreeVector &point,
    const G4ThreeVector &direction,
    G4double inwardEpsilon) const
{
  if (direction.mag2() <= 0.0)
    return WrapToPrimaryCell(point);
  return WrapToPrimaryCell(
      point + std::max(0.0, inwardEpsilon) * direction.unit());
}

DetectorConstruction::Phase DetectorConstruction::FindPhaseInsidePrimary(
    const G4ThreeVector &point,
    G4bool useMinimumImage) const
{
  if (fSpatialCells.empty() || fSpatialSpheres.empty())
    return Phase::Matrix;
  const auto cell = [](G4double coordinate, G4double length, G4int count)
  {
    return std::clamp(static_cast<G4int>(std::floor(
                          ((coordinate + 0.5 * length) / length) * count)),
                      0, count - 1);
  };
  const auto modulo = [](G4int value, G4int count)
  {
    const G4int result = value % count;
    return result < 0 ? result + count : result;
  };
  const auto flat = [&](G4int ix, G4int iy, G4int iz)
  {
    return static_cast<std::size_t>(ix + fSpatialNx * (iy + fSpatialNy * iz));
  };
  const G4int ix = cell(point.x(), fBoxX, fSpatialNx);
  const G4int iy = cell(point.y(), fBoxY, fSpatialNy);
  const G4int iz = cell(point.z(), fBoxZ, fSpatialNz);
  std::vector<std::size_t> visitedCells;
  visitedCells.reserve(27);
  for (G4int dz = -1; dz <= 1; ++dz)
  {
    for (G4int dy = -1; dy <= 1; ++dy)
    {
      for (G4int dx = -1; dx <= 1; ++dx)
      {
        const std::size_t cellIndex = flat(modulo(ix + dx, fSpatialNx),
                                           modulo(iy + dy, fSpatialNy),
                                           modulo(iz + dz, fSpatialNz));
        if (std::find(visitedCells.begin(), visitedCells.end(), cellIndex) != visitedCells.end())
          continue;
        visitedCells.push_back(cellIndex);
        for (std::size_t sphereIndex : fSpatialCells[cellIndex])
        {
          const SphereInfo &sphere = *fSpatialSpheres[sphereIndex];
          G4ThreeVector delta = point - sphere.center;
          if (useMinimumImage)
          {
            delta.setX(delta.x() - fBoxX * std::nearbyint(delta.x() / fBoxX));
            delta.setY(delta.y() - fBoxY * std::nearbyint(delta.y() / fBoxY));
            delta.setZ(delta.z() - fBoxZ * std::nearbyint(delta.z() / fBoxZ));
          }
          if (delta.mag2() <= sphere.radius * sphere.radius)
            return sphere.phase;
        }
      }
    }
  }
  return Phase::Matrix;
}

void DetectorConstruction::BuildPeriodicSpatialIndex()
{
  fSpatialSpheres.clear();
  fSpatialSpheres.reserve(fBNSpheres.size() + fZnSSpheres.size());
  for (const auto &sphere : fBNSpheres)
    fSpatialSpheres.push_back(&sphere);
  for (const auto &sphere : fZnSSpheres)
    fSpatialSpheres.push_back(&sphere);
  const G4double targetCellSize = std::max(2.0 * GetGlobalMaxParticleRadius() + fOverlapGap,
                                           1.0e-6 * um);
  fSpatialNx = std::max(1, static_cast<G4int>(std::floor(fBoxX / targetCellSize)));
  fSpatialNy = std::max(1, static_cast<G4int>(std::floor(fBoxY / targetCellSize)));
  fSpatialNz = std::max(1, static_cast<G4int>(std::floor(fBoxZ / targetCellSize)));
  fSpatialCellX = fBoxX / fSpatialNx;
  fSpatialCellY = fBoxY / fSpatialNy;
  fSpatialCellZ = fBoxZ / fSpatialNz;
  fSpatialCells.assign(static_cast<std::size_t>(fSpatialNx * fSpatialNy * fSpatialNz), {});
  const auto cell = [](G4double coordinate, G4double length, G4int count)
  {
    return std::clamp(static_cast<G4int>(std::floor(
                          ((coordinate + 0.5 * length) / length) * count)),
                      0, count - 1);
  };
  for (std::size_t sphereIndex = 0; sphereIndex < fSpatialSpheres.size(); ++sphereIndex)
  {
    const auto &center = fSpatialSpheres[sphereIndex]->center;
    const G4int ix = cell(center.x(), fBoxX, fSpatialNx);
    const G4int iy = cell(center.y(), fBoxY, fSpatialNy);
    const G4int iz = cell(center.z(), fBoxZ, fSpatialNz);
    const std::size_t flat = static_cast<std::size_t>(ix + fSpatialNx * (iy + fSpatialNy * iz));
    fSpatialCells[flat].push_back(sphereIndex);
  }
}

G4double DetectorConstruction::GetMinParticleRadius(Phase phase) const
{
  const auto &spheres = phase == Phase::BN ? fBNSpheres : fZnSSpheres;
  G4double result = std::numeric_limits<G4double>::infinity();
  for (const auto &sphere : spheres)
    result = std::min(result, sphere.radius);
  return std::isfinite(result) ? result : 0.0;
}

G4double DetectorConstruction::GetMaxParticleRadius(Phase phase) const
{
  const auto &spheres = phase == Phase::BN ? fBNSpheres : fZnSSpheres;
  G4double result = 0.0;
  for (const auto &sphere : spheres)
    result = std::max(result, sphere.radius);
  return result;
}

G4double DetectorConstruction::GetGlobalMaxParticleRadius() const
{
  return std::max(GetMaxParticleRadius(Phase::BN), GetMaxParticleRadius(Phase::ZnS));
}

void DetectorConstruction::RegisterLogicalPhase(G4LogicalVolume *logical, Phase phase)
{
  if (logical != nullptr)
    fLogicalPhase[logical] = phase;
}

DetectorConstruction::Phase DetectorConstruction::GetPhaseFromLogicalVolume(
    const G4LogicalVolume *logical) const
{
  if (logical == nullptr)
    return Phase::Unknown;
  if (logical == fMatrixLogical)
    return Phase::Matrix;
  if (logical == fWorldLogical)
    return Phase::World;
  const auto found = fLogicalPhase.find(logical);
  if (found != fLogicalPhase.end())
    return found->second;
  if (logical->GetMaterial() == fBNMaterial)
    return Phase::BN;
  if (logical->GetMaterial() == fZnSMaterial)
    return Phase::ZnS;
  return Phase::Unknown;
}

G4bool DetectorConstruction::IsBNLogical(const G4LogicalVolume *logical) const
{
  return GetPhaseFromLogicalVolume(logical) == Phase::BN;
}

G4bool DetectorConstruction::IsZnSLogical(const G4LogicalVolume *logical) const
{
  return GetPhaseFromLogicalVolume(logical) == Phase::ZnS;
}

const DetectorConstruction::GeometryCopyInfo *DetectorConstruction::GetGeometryCopyInfo(
    G4int copyId) const
{
  const auto found = fGeometryCopyById.find(copyId);
  return found == fGeometryCopyById.end() ? nullptr : &found->second;
}

const char *DetectorConstruction::PhaseName(Phase phase)
{
  switch (phase)
  {
  case Phase::BN:
    return "BN";
  case Phase::ZnS:
    return "ZnS";
  case Phase::Matrix:
    return "Matrix";
  case Phase::World:
    return "World";
  case Phase::Unknown:
  default:
    return "Unknown";
  }
}

// --------------------------------------------------------------------

void DetectorConstruction::NotifyGeometryChanged()
{
  auto *rm = G4RunManager::GetRunManager();
  if (rm)
  {
    rm->GeometryHasBeenModified();
  }
}

// --------------------------------------------------------------------

void DetectorConstruction::SetScreenThicknessUm(G4double thicknessUm)
{
  if (thicknessUm <= 0.0)
    return;
  fScreenThickness = thicknessUm * um;
  NotifyGeometryChanged();
}

// --------------------------------------------------------------------

void DetectorConstruction::SetMicroThicknessUm(G4double thicknessUm)
{
  if (thicknessUm <= 0.0)
    return;
  fMicroThickness = thicknessUm * um;
  NotifyGeometryChanged();
}

// --------------------------------------------------------------------

void DetectorConstruction::SetPatchXYUm(G4double patchXYUm)
{
  if (patchXYUm <= 0.0)
    return;
  fPatchXY = patchXYUm * um;
  NotifyGeometryChanged();
}

// --------------------------------------------------------------------

void DetectorConstruction::SetWeightRatio(G4double bnWt, G4double znsWt)
{
  if (bnWt <= 0.0 || znsWt <= 0.0)
    return;
  fBnWt = bnWt;
  fZnsWt = znsWt;
  NotifyGeometryChanged();
}

// --------------------------------------------------------------------

void DetectorConstruction::SetPresetRatio(G4int presetIndex)
{
  switch (presetIndex)
  {
  case 0:
    SetWeightRatio(2.0, 1.0);
    break;
  case 1:
    SetWeightRatio(1.0, 1.0);
    break;
  case 2:
    SetWeightRatio(1.0, 1.5);
    break;
  case 3:
    SetWeightRatio(1.0, 2.0);
    break;
  case 4:
    SetWeightRatio(1.0, 2.5);
    break;
  case 5:
    SetWeightRatio(1.0, 3.0);
    break;
  default:
    SetWeightRatio(1.0, 1.0);
    break;
  }
}

// --------------------------------------------------------------------

G4double DetectorConstruction::GetEffectiveLocalThickness() const
{
  if (fPlacementFormatVersion == 3)
    return fBoxZ;
  return std::min(fScreenThickness, fMicroThickness);
}

// --------------------------------------------------------------------

G4double DetectorConstruction::HashToUnit(G4int ix, G4int iy, G4int iz, G4int salt)
{
  std::uint32_t h = 2166136261u;
  h = (h ^ static_cast<std::uint32_t>(ix + 10007 + 13 * salt)) * 16777619u;
  h = (h ^ static_cast<std::uint32_t>(iy + 10009 + 17 * salt)) * 16777619u;
  h = (h ^ static_cast<std::uint32_t>(iz + 10037 + 19 * salt)) * 16777619u;

  return static_cast<G4double>(h & 0x00FFFFFFu) /
         static_cast<G4double>(0x01000000u);
}

// --------------------------------------------------------------------

G4bool DetectorConstruction::IsInsideSafeXY(const G4ThreeVector &pos) const
{
  const G4double halfXY = 0.5 * fPatchXY;
  return (std::abs(pos.x()) < (halfXY - fSafeMarginXY) &&
          std::abs(pos.y()) < (halfXY - fSafeMarginXY));
}

// --------------------------------------------------------------------

void DetectorConstruction::DefineMaterials()
{
  auto *nist = G4NistManager::Instance();

  // ---------- vacuum ----------
  fVacuumMaterial = nist->FindOrBuildMaterial("G4_Galactic");

  // ---------- common materials ----------
  auto *air = nist->FindOrBuildMaterial("G4_AIR");
  auto *binder = nist->FindOrBuildMaterial("G4_PLEXIGLASS");

  // ---------- enriched boron ----------
  auto *enrichedB = G4Element::GetElement("EnrichedBoron_10B_99p15", false);
  if (!enrichedB)
  {
    auto *B10 = new G4Isotope("B10_iso", 5, 10, 10.0129370 * g / mole);
    auto *B11 = new G4Isotope("B11_iso", 5, 11, 11.0093050 * g / mole);

    enrichedB = new G4Element("EnrichedBoron_10B_99p15", "B_enr", 2);
    enrichedB->AddIsotope(B10, 19.78 * perCent);
    enrichedB->AddIsotope(B11, 80.22 * perCent);
  }

  auto *N = nist->FindOrBuildElement("N");
  auto *Zn = nist->FindOrBuildElement("Zn");
  auto *S = nist->FindOrBuildElement("S");
  auto *Ag = nist->FindOrBuildElement("Ag");

  // ---------- BN ----------
  fBNMaterial = G4Material::GetMaterial("BN_10B99p15", false);
  if (!fBNMaterial)
  {
    fBNMaterial = new G4Material("BN_10B99p15", 2.10 * g / cm3, 2);
    fBNMaterial->AddElement(enrichedB, 1);
    fBNMaterial->AddElement(N, 1);
  }

  // ---------- ZnS base ----------
  auto *znsBase = G4Material::GetMaterial("ZnS_Base", false);
  if (!znsBase)
  {
    znsBase = new G4Material("ZnS_Base", 4.09 * g / cm3, 2);
    znsBase->AddElement(Zn, 1);
    znsBase->AddElement(S, 1);
  }

  // ---------- ZnS:Ag 0.1 wt% ----------
  fZnSMaterial = G4Material::GetMaterial("ZnS_Ag_0p1wt", false);
  if (!fZnSMaterial)
  {
    fZnSMaterial = new G4Material("ZnS_Ag_0p1wt", 4.09 * g / cm3, 2);
    fZnSMaterial->AddMaterial(znsBase, 99.9 * perCent);
    fZnSMaterial->AddElement(Ag, 0.1 * perCent);
  }

  // ---------- Binder + air voids background ----------
  fMatrixMaterial = G4Material::GetMaterial("BinderVoidMix", false);
  if (!fMatrixMaterial)
  {
    const G4double rhoBinder = binder->GetDensity();
    const G4double rhoAir = air->GetDensity();
    const G4double vfVoid = fVoidVolumeFraction;

    const G4double rhoMix =
        (1.0 - vfVoid) * rhoBinder + vfVoid * rhoAir;

    const G4double binderMass = (1.0 - vfVoid) * rhoBinder;
    const G4double airMass = vfVoid * rhoAir;

    const G4double binderMassFrac = binderMass / (binderMass + airMass);
    const G4double airMassFrac = airMass / (binderMass + airMass);

    fMatrixMaterial = new G4Material("BinderVoidMix", rhoMix, 2);
    fMatrixMaterial->AddMaterial(binder, binderMassFrac);
    fMatrixMaterial->AddMaterial(air, airMassFrac);
  }

  // =========================================================
  // 光学属性定义 (Optical Properties)
  // =========================================================

  // 定义光子能量范围：2.0 eV (~620 nm) 到 3.0 eV (~413 nm)
  // ZnS:Ag 的发光峰值通常在 450 nm 左右，约等于 2.75 eV
  G4double photonEnergy[] = {2.0 * eV, 2.75 * eV, 3.0 * eV};
  const G4int nEntries = sizeof(photonEnergy) / sizeof(G4double);

  // ---------------------------------------------------------
  // 1. 真空 (Vacuum)
  // ---------------------------------------------------------
  G4double rindexVacuum[] = {1.0, 1.0, 1.0};
  auto *mptVacuum = new G4MaterialPropertiesTable();
  mptVacuum->AddProperty("RINDEX", photonEnergy, rindexVacuum, nEntries);
  fVacuumMaterial->SetMaterialPropertiesTable(mptVacuum);

  // ---------------------------------------------------------
  // 2. 基质/空气孔隙混合物 (Binder + Voids Matrix)
  // ---------------------------------------------------------
  const G4double matrixRIndex =
      (fConfig != nullptr && fConfig->opticalMatrixRIndex > 0.0)
          ? fConfig->opticalMatrixRIndex
          : ComputeMatrixRIndexLorentzLorenz();
  const G4double matrixAbsLength = (fConfig ? fConfig->opticalMatrixAbsLengthUm : 1.0e6) * um;
  if (fConfig != nullptr && fConfig->opticalMatrixRIndex <= 0.0)
    fConfig->opticalMatrixRIndex = matrixRIndex;
  G4double rindexMatrix[] = {matrixRIndex, matrixRIndex, matrixRIndex};
  G4double absMatrix[] = {matrixAbsLength, matrixAbsLength, matrixAbsLength};
  auto *mptMatrix = new G4MaterialPropertiesTable();
  mptMatrix->AddProperty("RINDEX", photonEnergy, rindexMatrix, nEntries);
  mptMatrix->AddProperty("ABSLENGTH", photonEnergy, absMatrix, nEntries);
  fMatrixMaterial->SetMaterialPropertiesTable(mptMatrix);

  // ---------------------------------------------------------
  // 3. 氮化硼 (BN)
  // ---------------------------------------------------------
  const G4double bnRIndex = fConfig ? fConfig->opticalBnRIndex : 2.1;
  const G4double bnAbsLength = (fConfig ? fConfig->opticalBnAbsLengthUm : 1.0e6) * um;
  G4double rindexBN[] = {bnRIndex, bnRIndex, bnRIndex};
  G4double absBN[] = {bnAbsLength, bnAbsLength, bnAbsLength};
  auto *mptBN = new G4MaterialPropertiesTable();
  mptBN->AddProperty("RINDEX", photonEnergy, rindexBN, nEntries);
  mptBN->AddProperty("ABSLENGTH", photonEnergy, absBN, nEntries);
  fBNMaterial->SetMaterialPropertiesTable(mptBN);

  // ---------------------------------------------------------
  // 4. 闪烁体 (ZnS:Ag 0.1 wt%)
  // ---------------------------------------------------------
  const G4double znsRIndex = fConfig ? fConfig->opticalZnsRIndex : 2.36;
  const G4double znsAbsLength = (fConfig ? fConfig->opticalZnsAbsLengthUm : 1.0e4) * um;
  G4double rindexZnS[] = {znsRIndex, znsRIndex, znsRIndex};
  G4double absZnS[] = {znsAbsLength, znsAbsLength, znsAbsLength};

  // 发光光谱相对强度 (在 2.75 eV 时达到峰值 1.0)
  G4double emissionZnS[] = {0.1, 1.0, 0.2};

  auto *mptZnS = new G4MaterialPropertiesTable();
  mptZnS->AddProperty("RINDEX", photonEnergy, rindexZnS, nEntries);
  mptZnS->AddProperty("ABSLENGTH", photonEnergy, absZnS, nEntries);

  // ----- 核心发光属性 -----
  // 发射谱分布 (Geant4 10.7 及以上版本使用 SCINTILLATIONCOMPONENT1)
  // 注意：如果你的 G4 版本低于 10.7，请将 SCINTILLATIONCOMPONENT1 替换为 FASTCOMPONENT
  mptZnS->AddProperty("SCINTILLATIONCOMPONENT1", photonEnergy, emissionZnS, nEntries);

  const G4bool enableZnSScintillation =
      (fConfig != nullptr) &&
      fConfig->runMode == RunMode::StageD_OpticalHomogenization;

  // Stage A/B 不需要生成 optical photons；将 scintillation yield 设为 0
  // 以避免误注册光学物理时引入灾难性的 photon tracking 开销。
  mptZnS->AddConstProperty(
      "SCINTILLATIONYIELD",
      enableZnSScintillation ? (75000. / MeV) : 0.0);

  // 展宽系数，设为 1.0 表现出符合泊松分布的光子涨落
  mptZnS->AddConstProperty("RESOLUTIONSCALE", 1.0);

  // ZnS:Ag 衰减时间较长，主峰大约 200 ns
  // 注意：低版本 G4 请将 SCINTILLATIONTIMECONSTANT1 替换为 FASTTIMECONSTANT
  mptZnS->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 200. * ns);
  mptZnS->AddConstProperty("SCINTILLATIONYIELD1", 1.0);

  fZnSMaterial->SetMaterialPropertiesTable(mptZnS);

  if (fConfig != nullptr)
  {
    if (!fConfig->opticalParamsProvided &&
        fConfig->runMode == RunMode::StageD_OpticalHomogenization)
    {
      G4cout << "[DetectorConstruction] Warning: optical parameters are using built-in defaults/estimates for "
             << AnalysisConfig::RunModeName(fConfig->runMode) << ". "
             << "Built-in absorption defaults keep matrix/BN near-transparent and assign a shorter ZnS absorption "
             << "length as a fallback, not as phase-resolved material data. "
             << "Set /cfg/setOpticalParams matrix_n matrix_abs_um bn_n bn_abs_um zns_n zns_abs_um for physics runs."
             << G4endl;
    }

    if (!fConfig->opticalParamsProvided)
    {
      G4cout << "[DetectorConstruction] Auto-estimated matrix n with Lorentz-Lorenz:"
             << " binder_n=" << kBinderPlexiglassRIndex
             << " binder_frac=" << (gBinderPart / (gBinderPart + gAirPart))
             << " air_frac=" << (gAirPart / (gBinderPart + gAirPart))
             << " -> matrix_n=" << matrixRIndex
             << G4endl;
    }

    G4cout << "[DetectorConstruction] Optical params:"
           << " matrix n/abs_um=" << matrixRIndex << "/" << matrixAbsLength / um
           << " BN n/abs_um=" << bnRIndex << "/" << bnAbsLength / um
           << " ZnS n/abs_um=" << znsRIndex << "/" << znsAbsLength / um
           << G4endl;
  }
}

// --------------------------------------------------------------------

G4int DetectorConstruction::ComputeTargetZnSCount(G4int placedBNCount, G4int usableZnSCount) const
{
  if (placedBNCount <= 0 || usableZnSCount <= 0)
    return 0;

  const G4double vBN = (4.0 / 3.0) * pi * std::pow(fBNRadius, 3);
  const G4double vZnS = (4.0 / 3.0) * pi * std::pow(fZnSRadius, 3);

  const G4double rhoBN = 2.10 * g / cm3;
  const G4double rhoZnS = 4.09 * g / cm3;

  const G4double mBN = rhoBN * vBN;
  const G4double mZnS = rhoZnS * vZnS;

  const G4double targetZnSCountReal =
      static_cast<G4double>(placedBNCount) *
      (mBN / mZnS) *
      (fZnsWt / fBnWt);

  G4int targetZnSCount = static_cast<G4int>(std::llround(targetZnSCountReal));
  targetZnSCount = std::max(0, std::min(targetZnSCount, usableZnSCount));

  return targetZnSCount;
}

// --------------------------------------------------------------------

G4VPhysicalVolume *DetectorConstruction::ConstructVersion3(const RvePlacementData &data)
{
  fPlacementFormatVersion = data.metadata.formatVersion;
  fPlacementSeed = data.metadata.seed;
  fPlacementPhiAchieved = data.phiAchieved;
  fPlacementZnSToBNMassRatio = data.znsToBnMassRatio;
  fLoadedPlacementSeedBase = std::to_string(data.metadata.seed);
  fLoadedPeriodicImagesFile = data.periodicImagesPath.string();
  fLoadedPeriodicImagesFileForRecord = AnalysisConfig::PathForRecord(data.periodicImagesPath);
  fBoxX = data.metadata.boxXUm * um;
  fBoxY = data.metadata.boxYUm * um;
  fBoxZ = data.metadata.boxZUm * um;
  fPatchXY = fBoxX;
  fMicroThickness = fBoxZ;
  fOverlapGap = data.metadata.overlapGapUm * um;
  if (fConfig != nullptr)
  {
    fConfig->patchXY_um = data.metadata.boxXUm;
    fConfig->microThickness_um = data.metadata.boxZUm;
  }

  const auto phaseFromRve = [](RvePhase phase)
  {
    return phase == RvePhase::BN ? Phase::BN :
           phase == RvePhase::ZnS ? Phase::ZnS : Phase::Unknown;
  };
  fBNSpheres.reserve(data.bnCount);
  fZnSSpheres.reserve(data.znsCount);
  fPlacedBNCenters.reserve(data.bnCount);
  fPlacedZnSCenters.reserve(data.znsCount);
  for (const auto &particle : data.particles)
  {
    SphereInfo sphere;
    sphere.particleId = particle.particleId;
    sphere.phaseId = particle.phaseId;
    sphere.radiusBinId = particle.radiusBinId;
    sphere.radiusClassId = particle.radiusClassId;
    sphere.phase = phaseFromRve(particle.phase);
    sphere.radiusClass = particle.radiusClass;
    sphere.center = particle.center;
    sphere.radius = particle.radiusUm * um;
    if (sphere.phase == Phase::BN)
    {
      fBNSpheres.push_back(sphere);
      fPlacedBNCenters.push_back(sphere.center);
    }
    else
    {
      fZnSSpheres.push_back(sphere);
      fPlacedZnSCenters.push_back(sphere.center);
    }
  }
  BuildPeriodicSpatialIndex();
  DefineMaterials();

  const G4double margin = 50.0 * um;
  auto *solidWorld = new G4Box("WorldSolid",
                               0.5 * fBoxX + margin,
                               0.5 * fBoxY + margin,
                               0.5 * fBoxZ + margin);
  fWorldLogical = new G4LogicalVolume(solidWorld, fVacuumMaterial, "WorldLV");
  auto *physWorld = new G4PVPlacement(
      nullptr, G4ThreeVector(), fWorldLogical, "WorldPV", nullptr, false, 0, false);

  auto *solidMatrix = new G4Box("MicroPatchSolid", 0.5 * fBoxX, 0.5 * fBoxY, 0.5 * fBoxZ);
  fMatrixLogical = new G4LogicalVolume(solidMatrix, fMatrixMaterial, "MatrixLV");
  new G4PVPlacement(nullptr, G4ThreeVector(), fMatrixLogical, "MatrixPV",
                    fWorldLogical, false, 0, false);
  fScoringVolume = fMatrixLogical;

  auto *worldVis = new G4VisAttributes();
  worldVis->SetVisibility(false);
  fWorldLogical->SetVisAttributes(worldVis);
  auto *matrixVis = new G4VisAttributes();
  matrixVis->SetVisibility(false);
  fMatrixLogical->SetVisAttributes(matrixVis);
  auto *bnVis = new G4VisAttributes(G4Colour(0.20, 0.60, 1.00));
  bnVis->SetForceSolid(true);
  auto *znsVis = new G4VisAttributes(G4Colour(1.00, 0.85, 0.20));
  znsVis->SetForceSolid(true);
  auto *bnClipVis = new G4VisAttributes(G4Colour(0.20, 0.60, 1.00, 0.20));
  bnClipVis->SetForceSolid(true);
  auto *znsClipVis = new G4VisAttributes(G4Colour(1.00, 0.85, 0.20, 0.20));
  znsClipVis->SetForceSolid(true);

  for (const auto &item : data.radiusClasses)
  {
    const auto &radiusClass = item.second;
    RadiusClassGeometry entry;
    entry.phase = phaseFromRve(radiusClass.phase);
    entry.radius = radiusClass.radiusUm * um;
    const std::string suffix = std::to_string(radiusClass.radiusClassId);
    entry.solid = new G4Orb("RVE_Sphere_RC_" + suffix, entry.radius);
    entry.logical = new G4LogicalVolume(
        entry.solid,
        entry.phase == Phase::BN ? fBNMaterial : fZnSMaterial,
        "RVE_LV_RC_" + suffix);
    entry.logical->SetVisAttributes(entry.phase == Phase::BN ? bnVis : znsVis);
    RegisterLogicalPhase(entry.logical, entry.phase);
    fRadiusClassCache.emplace(radiusClass.radiusClassId, entry);
  }

  auto *solidClipBox = new G4Box("ClipBoxSolid", 0.5 * fBoxX, 0.5 * fBoxY, 0.5 * fBoxZ);
  const bool checkOverlaps = IsEnvFlagEnabled("BNZS_CHECK_OVERLAPS");
  const bool primaryOnly = fConfig != nullptr && fConfig->placementGeometryMode == "primary_only";
  if (fConfig != nullptr && !primaryOnly && fConfig->placementGeometryMode != "pbc_clipped")
  {
    G4Exception("DetectorConstruction::ConstructVersion3", "BNZS_RVE_020", FatalException,
                ("Unsupported placementGeometryMode=" + fConfig->placementGeometryMode).c_str());
  }

  G4int clippedCount = 0;
  G4int placedCount = 0;
  const auto placeCopy = [&](G4int copyId,
                             G4int sourceParticleId,
                             G4bool isPrimary,
                             G4int shiftIx,
                             G4int shiftIy,
                             G4int shiftIz,
                             G4int radiusClassId,
                             Phase phase,
                             G4double radius,
                             const G4ThreeVector &center)
  {
    const G4double halfX = 0.5 * fBoxX;
    const G4double halfY = 0.5 * fBoxY;
    const G4double halfZ = 0.5 * fBoxZ;
    const bool intersects = center.x() + radius > -halfX && center.x() - radius < halfX &&
                            center.y() + radius > -halfY && center.y() - radius < halfY &&
                            center.z() + radius > -halfZ && center.z() - radius < halfZ;
    if (!intersects)
    {
      G4Exception("DetectorConstruction::ConstructVersion3", "BNZS_RVE_021", JustWarning,
                  ("Skipping non-intersecting geometry copy_id=" + std::to_string(copyId)).c_str());
      return;
    }
    const auto cache = fRadiusClassCache.find(radiusClassId);
    if (cache == fRadiusClassCache.end() || cache->second.phase != phase ||
        std::abs(cache->second.radius - radius) > 1.0e-9 * um)
    {
      G4Exception("DetectorConstruction::ConstructVersion3", "BNZS_RVE_022", FatalException,
                  ("Invalid radius class for copy_id=" + std::to_string(copyId)).c_str());
      return;
    }
    const bool fullyInside = std::abs(center.x()) + radius <= halfX &&
                             std::abs(center.y()) + radius <= halfY &&
                             std::abs(center.z()) + radius <= halfZ;
    G4LogicalVolume *logical = cache->second.logical;
    std::string physicalName = "RVE_PV_" + std::to_string(copyId);
    if (!fullyInside)
    {
      const std::string suffix = std::to_string(copyId);
      auto *clippedSolid = new G4IntersectionSolid(
          "RVE_ClipSolid_" + suffix, cache->second.solid, solidClipBox,
          nullptr, G4ThreeVector(-center.x(), -center.y(), -center.z()));
      logical = new G4LogicalVolume(
          clippedSolid, phase == Phase::BN ? fBNMaterial : fZnSMaterial,
          "RVE_ClipLV_" + suffix);
      logical->SetVisAttributes(phase == Phase::BN ? bnClipVis : znsClipVis);
      RegisterLogicalPhase(logical, phase);
      ++clippedCount;
    }
    new G4PVPlacement(nullptr, center, logical, physicalName,
                      fMatrixLogical, false, copyId, checkOverlaps);
    GeometryCopyInfo info;
    info.copyId = copyId;
    info.sourceParticleId = sourceParticleId;
    info.isPrimary = isPrimary;
    info.shiftIx = shiftIx;
    info.shiftIy = shiftIy;
    info.shiftIz = shiftIz;
    info.radiusClassId = radiusClassId;
    info.phase = phase;
    info.radius = radius;
    if (!fGeometryCopyById.emplace(copyId, info).second)
    {
      G4Exception("DetectorConstruction::ConstructVersion3", "BNZS_RVE_023", FatalException,
                  ("Duplicate physical copy number=" + std::to_string(copyId)).c_str());
    }
    ++placedCount;
  };

  if (primaryOnly)
  {
    for (const auto &particle : data.particles)
    {
      placeCopy(particle.particleId, particle.particleId, true, 0, 0, 0,
                particle.radiusClassId, phaseFromRve(particle.phase),
                particle.radiusUm * um, particle.center);
    }
    fLoadedPeriodicImagesFile = "";
    fLoadedPeriodicImagesFileForRecord = "";
  }
  else
  {
    for (const auto &copy : data.geometryCopies)
    {
      placeCopy(copy.copyId, copy.sourceParticleId, copy.isPrimary,
                copy.shiftIx, copy.shiftIy, copy.shiftIz,
                copy.radiusClassId, phaseFromRve(copy.phase),
                copy.radiusUm * um, copy.center);
    }
  }

  const G4double bnMassG = data.bnVolumeUm3 * 2.10e-12;
  const G4double znsMassG = data.znsVolumeUm3 * 4.09e-12;
  G4cout << "\n[DetectorConstruction] format_version=3 periodic RVE loaded"
         << "\n  placement file            = " << fLoadedPlacementFileForRecord
         << "\n  periodic images file      = "
         << (fLoadedPeriodicImagesFileForRecord.empty() ? "disabled" : fLoadedPeriodicImagesFileForRecord)
         << "\n  geometry mode             = " << (primaryOnly ? "primary_only" : "pbc_clipped")
         << "\n  format version            = " << fPlacementFormatVersion
         << "\n  seed                      = " << fPlacementSeed
         << "\n  box (um)                  = " << GetBoxXUm() << " x " << GetBoxYUm() << " x " << GetBoxZUm()
         << "\n  unique particles          = " << data.particles.size()
         << "\n  BN / ZnS unique count     = " << data.bnCount << " / " << data.znsCount
         << "\n  radius classes            = " << data.radiusClasses.size()
         << "\n  geometry copies           = " << placedCount
         << "\n  clipped geometry copies   = " << clippedCount
         << std::setprecision(15)
         << "\n  phi achieved              = " << data.phiAchieved
         << "\n  BN:ZnS mass               = 1:" << (bnMassG > 0.0 ? znsMassG / bnMassG : 0.0)
         << std::setprecision(6)
         << G4endl;
  return physWorld;
}

// --------------------------------------------------------------------

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  fPlacedBNCenters.clear();
  fSafeBNCenters.clear();
  fUsableZnSCandidateCenters.clear();
  fPlacedZnSCenters.clear();
  fBNSpheres.clear();
  fZnSSpheres.clear();
  fRadiusClassCache.clear();
  fLogicalPhase.clear();
  fGeometryCopyById.clear();
  fSpatialSpheres.clear();
  fSpatialCells.clear();
  fLoadedPeriodicImagesFile = "";
  fLoadedPeriodicImagesFileForRecord = "";
  fLoadedPlacementSeedBase = "";
  fPlacementFormatVersion = 0;
  fPlacementSeed = 0;
  fPlacementPhiAchieved = 0.0;
  fPlacementZnSToBNMassRatio = 0.0;

  // ---- sync latest config before geometry construction ----
  if (fConfig != nullptr)
  {
    fPatchXY = fConfig->patchXY_um * um;
    fMicroThickness = fConfig->microThickness_um * um;
    fBoxX = fPatchXY;
    fBoxY = fPatchXY;
    fBoxZ = fMicroThickness;
    fBnWt = fConfig->bnWt;
    fZnsWt = fConfig->znsWt;
  }

  // =========================================================
  // Load external placement file instead of generating packing internally
  // =========================================================
  G4cout << "[DetectorConstruction] Using BN:ZnS weight ratio = "
         << fBnWt << ":" << fZnsWt << G4endl;

  if (fConfig != nullptr &&
      !fConfig->useRandomPlacement &&
      !fConfig->placementFilePath.empty())
  {
    fLoadedPlacementFile = ResolvePlacementPathForCurrentProject(fConfig->placementFilePath).string();
  }
  else
  {
    fLoadedPlacementFile = PickRandomPlacementFilePath(fBnWt, fZnsWt);
  }

  G4cout << "[DetectorConstruction] Loading placement file: "
         << fLoadedPlacementFile << G4endl;

  fLoadedPlacementFileForRecord = AnalysisConfig::PathForRecord(fLoadedPlacementFile);

  const std::filesystem::path loadedPlacementPath(fLoadedPlacementFile.c_str());
  if (RvePlacementReader::IsFormatVersion3(loadedPlacementPath))
  {
    const bool readPeriodicImages =
        fConfig == nullptr || fConfig->placementGeometryMode != "primary_only";
    std::filesystem::path periodicOverride;
    if (readPeriodicImages && fConfig != nullptr && !fConfig->periodicImagesFilePath.empty())
      periodicOverride = ResolvePlacementPathForCurrentProject(fConfig->periodicImagesFilePath);
    const RvePlacementData version3 =
        RvePlacementReader::ReadVersion3(
            loadedPlacementPath, periodicOverride, readPeriodicImages);
    return ConstructVersion3(version3);
  }

  DefineMaterials();

  const G4double localThickness = GetEffectiveLocalThickness();
  fBoxX = fPatchXY;
  fBoxY = fPatchXY;
  fBoxZ = localThickness;

  // ---------- legacy world and local patch ----------
  const G4double worldXY = 2 * fPatchXY;
  const G4double worldZ = localThickness + 100.0 * um;
  auto *solidWorld = new G4Box("WorldSolid", 0.5 * worldXY, 0.5 * worldXY, 0.5 * worldZ);
  fWorldLogical = new G4LogicalVolume(solidWorld, fVacuumMaterial, "WorldLV");
  auto *physWorld = new G4PVPlacement(
      nullptr, G4ThreeVector(), fWorldLogical, "WorldPV", nullptr, false, 0, false);
  auto *solidMatrix = new G4Box("MicroPatchSolid", 0.5 * fPatchXY, 0.5 * fPatchXY, 0.5 * localThickness);
  fMatrixLogical = new G4LogicalVolume(solidMatrix, fMatrixMaterial, "MatrixLV");
  new G4PVPlacement(nullptr, G4ThreeVector(), fMatrixLogical, "MatrixPV",
                    fWorldLogical, false, 0, false);
  fScoringVolume = fMatrixLogical;

  PlacementData loaded = ReadPlacementCSV(fLoadedPlacementFile);
  fLoadedPlacementSeedBase = loaded.seedBase.empty() ? "unknown" : loaded.seedBase;

  const G4double metadataTol = 1.0e-9;

  if (loaded.patchXYUm > 0.0 && std::abs(loaded.patchXYUm - fPatchXY / um) > metadataTol)
  {
    G4Exception("DetectorConstruction::Construct",
                "BNZS208", FatalException,
                ("Placement patchXY_um does not match current geometry: "
                 + std::to_string(loaded.patchXYUm) + " vs " + std::to_string(fPatchXY / um))
                    .c_str());
  }

  if (loaded.localThicknessUm > 0.0 &&
      std::abs(loaded.localThicknessUm - localThickness / um) > metadataTol)
  {
    G4Exception("DetectorConstruction::Construct",
                "BNZS209", FatalException,
                ("Placement localThickness_um does not match current geometry: "
                 + std::to_string(loaded.localThicknessUm) + " vs " + std::to_string(localThickness / um))
                    .c_str());
  }

  if (loaded.bnWt > 0.0 && std::abs(loaded.bnWt - fBnWt) > metadataTol)
  {
    G4Exception("DetectorConstruction::Construct",
                "BNZS210", FatalException,
                ("Placement bnWt does not match current geometry ratio: "
                 + std::to_string(loaded.bnWt) + " vs " + std::to_string(fBnWt))
                    .c_str());
  }

  if (loaded.znsWt > 0.0 && std::abs(loaded.znsWt - fZnsWt) > metadataTol)
  {
    G4Exception("DetectorConstruction::Construct",
                "BNZS211", FatalException,
                ("Placement znsWt does not match current geometry ratio: "
                 + std::to_string(loaded.znsWt) + " vs " + std::to_string(fZnsWt))
                    .c_str());
  }

  if (loaded.bnRadiusUm > 0.0)
    fBNRadius = loaded.bnRadiusUm * um;

  if (loaded.znsRadiusUm > 0.0)
    fZnSRadius = loaded.znsRadiusUm * um;

  if (loaded.overlapGapUm >= 0.0)
    fOverlapGap = loaded.overlapGapUm * um;

  // ---------- particle solids ----------
  auto *solidBN = new G4Orb("BNSolid", fBNRadius);
  auto *solidZnS = new G4Orb("ZnSSolid", fZnSRadius);

  fBNLogical = new G4LogicalVolume(solidBN, fBNMaterial, "BN_LV");
  fZnSLogical = new G4LogicalVolume(solidZnS, fZnSMaterial, "ZnS_LV");
  RegisterLogicalPhase(fBNLogical, Phase::BN);
  RegisterLogicalPhase(fZnSLogical, Phase::ZnS);

  ValidatePlacementData(loaded,
                        fPatchXY,
                        localThickness,
                        fBNRadius,
                        fZnSRadius,
                        fOverlapGap);

  fPlacedBNCenters = loaded.bnCenters;
  fPlacedZnSCenters = loaded.znsCenters;

  fBNSpheres.reserve(fPlacedBNCenters.size());
  for (const auto &center : fPlacedBNCenters)
  {
    SphereInfo sphere;
    sphere.phase = Phase::BN;
    sphere.center = center;
    sphere.radius = fBNRadius;
    fBNSpheres.push_back(sphere);
  }

  fZnSSpheres.reserve(fPlacedZnSCenters.size());
  for (const auto &center : fPlacedZnSCenters)
  {
    SphereInfo sphere;
    sphere.phase = Phase::ZnS;
    sphere.center = center;
    sphere.radius = fZnSRadius;
    fZnSSpheres.push_back(sphere);
  }
  BuildPeriodicSpatialIndex();

  fSafeBNCenters.clear();
  for (const auto &pos : fPlacedBNCenters)
  {
    if (IsInsideSafeXY(pos))
      fSafeBNCenters.push_back(pos);
  }

  // ---------- quantities for report ----------
  const G4double vBN = (4.0 / 3.0) * pi * std::pow(fBNRadius, 3);
  const G4double vZnS = (4.0 / 3.0) * pi * std::pow(fZnSRadius, 3);

  const G4double rhoBN = 2.10 * g / cm3;
  const G4double rhoZnS = 4.09 * g / cm3;

  const G4double volumeUnit_um3 = um * um * um;

  // whole matrix-box volume
  const G4double Vpatch = fPatchXY * fPatchXY * localThickness;

  // actual particle volume kept inside the box
  G4double insideBNVolume = 0.0;
  G4double insideZnSVolume = 0.0;

  // debug counters for clipped particles
  G4int clippedBNCount = 0;
  G4int clippedZnSCount = 0;

  G4double minClippedBNVolume = -1.0;
  G4double minClippedZnSVolume = -1.0;

  // ---------- visualization ----------
  auto *worldVis = new G4VisAttributes();
  worldVis->SetVisibility(false);
  fWorldLogical->SetVisAttributes(worldVis);

  // 调试阶段：先把 matrix 隐掉，避免挡住裁剪体
  auto *matrixVis = new G4VisAttributes();
  matrixVis->SetVisibility(false);
  fMatrixLogical->SetVisAttributes(matrixVis);

  // 正常完整球
  auto *bnVis = new G4VisAttributes(G4Colour(0.20, 0.60, 1.00)); // 蓝
  bnVis->SetForceSolid(true);
  fBNLogical->SetVisAttributes(bnVis);

  auto *znsVis = new G4VisAttributes(G4Colour(1.00, 0.85, 0.20)); // 黄
  znsVis->SetForceSolid(true);
  fZnSLogical->SetVisAttributes(znsVis);

  // 裁剪球专用显示
  auto *bnClipVis = new G4VisAttributes(G4Colour(0.20, 0.60, 1.00, 0.20));
  bnClipVis->SetForceSolid(true);

  auto *znsClipVis = new G4VisAttributes(G4Colour(1.00, 0.85, 0.20, 0.20));
  znsClipVis->SetForceSolid(true);

  // =========================================================
  // 4) Build physical placements only after final centers are ready
  // =========================================================
  auto *solidClipBox = new G4Box("ClipBoxSolid",
                                 0.5 * fPatchXY,
                                 0.5 * fPatchXY,
                                 0.5 * localThickness);
  const bool checkOverlaps = IsEnvFlagEnabled("BNZS_CHECK_OVERLAPS");
  const bool geometryVerbose = IsEnvFlagEnabled("BNZS_GEOMETRY_VERBOSE");

  G4int copyBN = 0;
  for (const auto &pos : fPlacedBNCenters)
  {
    if (IsFullyInsideBox(pos, fBNRadius, fPatchXY, localThickness))
    {
      new G4PVPlacement(
          nullptr, pos, fBNLogical, "BN_PV",
          fMatrixLogical, false, copyBN++, false);

      insideBNVolume += vBN;
    }
    else
    {
      const std::string solidName = "BN_ClipSolid_" + std::to_string(copyBN);
      const std::string lvName = "BN_ClipLV_" + std::to_string(copyBN);
      const std::string pvName = "BN_ClipPV_" + std::to_string(copyBN);

      auto *clippedSolid = new G4IntersectionSolid(
          solidName,
          solidBN,
          solidClipBox,
          nullptr,
          G4ThreeVector(-pos.x(), -pos.y(), -pos.z()));

      const G4double clipVol = CalculateClippedVolume(pos, fBNRadius, fPatchXY, localThickness);

      auto *clippedLV = new G4LogicalVolume(
          clippedSolid, fBNMaterial, lvName);
      RegisterLogicalPhase(clippedLV, Phase::BN);
      clippedLV->SetVisAttributes(bnClipVis);

      new G4PVPlacement(
          nullptr, pos, clippedLV, pvName,
          fMatrixLogical, false, copyBN++, checkOverlaps);

      insideBNVolume += clipVol;
      ++clippedBNCount;

      if (minClippedBNVolume < 0.0 || clipVol < minClippedBNVolume)
        minClippedBNVolume = clipVol;

      if (geometryVerbose)
      {
        G4cout
            << "[Clip BN] center = ("
            << pos.x() / um << ", "
            << pos.y() / um << ", "
            << pos.z() / um << ") um"
            << "  volume = " << clipVol / volumeUnit_um3 << " um^3"
            << G4endl;
      }
    }
  }

  G4int copyZnS = 0;
  for (const auto &pos : fPlacedZnSCenters)
  {
    if (IsFullyInsideBox(pos, fZnSRadius, fPatchXY, localThickness))
    {
      new G4PVPlacement(
          nullptr, pos, fZnSLogical, "ZnS_PV",
          fMatrixLogical, false, copyZnS++, false);

      insideZnSVolume += vZnS;
    }
    else
    {
      const std::string solidName = "ZnS_ClipSolid_" + std::to_string(copyZnS);
      const std::string lvName = "ZnS_ClipLV_" + std::to_string(copyZnS);
      const std::string pvName = "ZnS_ClipPV_" + std::to_string(copyZnS);

      auto *clippedSolid = new G4IntersectionSolid(
          solidName,
          solidZnS,
          solidClipBox,
          nullptr,
          G4ThreeVector(-pos.x(), -pos.y(), -pos.z()));

      const G4double clipVol = CalculateClippedVolume(pos, fZnSRadius, fPatchXY, localThickness);

      auto *clippedLV = new G4LogicalVolume(
          clippedSolid, fZnSMaterial, lvName);
      RegisterLogicalPhase(clippedLV, Phase::ZnS);
      clippedLV->SetVisAttributes(znsClipVis);

      new G4PVPlacement(
          nullptr, pos, clippedLV, pvName,
          fMatrixLogical, false, copyZnS++, checkOverlaps);

      insideZnSVolume += clipVol;
      ++clippedZnSCount;

      if (minClippedZnSVolume < 0.0 || clipVol < minClippedZnSVolume)
        minClippedZnSVolume = clipVol;

      if (geometryVerbose)
      {
        G4cout
            << "[Clip ZnS] center = ("
            << pos.x() / um << ", "
            << pos.y() / um << ", "
            << pos.z() / um << ") um"
            << "  volume = " << clipVol / volumeUnit_um3 << " um^3"
            << G4endl;
      }
    }
  }

  // ---------- achieved ratio report ----------
  const G4double totalBNMass = rhoBN * insideBNVolume;
  const G4double totalZnSMass = rhoZnS * insideZnSVolume;

  const G4double particleVolume = insideBNVolume + insideZnSVolume;

  const G4double phiAchieved =
      (Vpatch > 0.0) ? (particleVolume / Vpatch) : 0.0;

  const G4double targetTotalPart = gSolidPart + gBinderPart + gAirPart;

  const G4double targetSolidFrac =
      (targetTotalPart > 0.0) ? (gSolidPart / targetTotalPart) : 0.0;
  const G4double targetBinderFrac =
      (targetTotalPart > 0.0) ? (gBinderPart / targetTotalPart) : 0.0;
  const G4double targetAirFrac =
      (targetTotalPart > 0.0) ? (gAirPart / targetTotalPart) : 0.0;

  // actual overall fractions in the whole patch
  const G4double achievedSolidFrac = phiAchieved;
  const G4double achievedBinderFrac = (1.0 - phiAchieved) * (1.0 - fVoidVolumeFraction);
  const G4double achievedAirFrac = (1.0 - phiAchieved) * fVoidVolumeFraction;

  const G4double achievedBNToZnSMass =
      (totalZnSMass > 0.0) ? (totalBNMass / totalZnSMass) : 0.0;

  const G4double achievedZnSToBNMass =
      (totalBNMass > 0.0) ? (totalZnSMass / totalBNMass) : 0.0;

  G4cout
      << "\n[DetectorConstruction] External placement loaded"
      << "\n  placement file            = " << fLoadedPlacementFileForRecord
      << "\n  screen thickness          = " << fScreenThickness / um << " um"
      << "\n  local thickness           = " << localThickness / um << " um"
      << "\n  patch XY                  = " << fPatchXY / um << " um"
      << "\n  target BN:ZnS (wt)        = " << fBnWt << " : " << fZnsWt
      << "\n  BN radius                 = " << fBNRadius / um << " um"
      << "\n  ZnS radius                = " << fZnSRadius / um << " um"
      << "\n  BN pitch                  = " << fBNPitch / um << " um"
      << "\n  ZnS pitch                 = " << fZnSPitch / um << " um"
      << "\n  safe margin XY            = " << fSafeMarginXY / um << " um"
      << "\n  void fraction             = " << 100.0 * fVoidVolumeFraction << " %"
      << "\n  placed BN spheres         = " << fPlacedBNCenters.size()
      << "\n  safe BN spheres           = " << fSafeBNCenters.size()
      << "\n  placed ZnS spheres        = " << fPlacedZnSCenters.size()
      << "\n  clipped BN spheres        = " << clippedBNCount
      << "\n  clipped ZnS spheres       = " << clippedZnSCount
      << "\n  achieved BN mass          = " << totalBNMass / g << " g"
      << "\n  achieved ZnS mass         = " << totalZnSMass / g << " g"
      << "\n  achieved BN/ZnS mass      = " << achievedBNToZnSMass
      << "\n  achieved ZnS/BN mass      = " << achievedZnSToBNMass
      << "\n  phi achieved              = " << phiAchieved
      << "\n  [solid | binder | air] target   = "
      << gSolidPart << " | " << gBinderPart << " | " << gAirPart
      << "\n  [solid | binder | air] achieved = "
      << 100.0 * achievedSolidFrac << " | "
      << 100.0 * achievedBinderFrac << " | "
      << 100.0 * achievedAirFrac;

  if (clippedBNCount > 0)
  {
    G4cout << "\n  min clipped BN volume     = "
           << minClippedBNVolume / volumeUnit_um3 << " um^3";
  }

  if (clippedZnSCount > 0)
  {
    G4cout << "\n  min clipped ZnS volume    = "
           << minClippedZnSVolume / volumeUnit_um3 << " um^3";
  }

  G4cout << G4endl;

  return physWorld;
}
