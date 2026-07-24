#include "StageDReentrySampler.hh"

#include "AnalysisConfig.hh"

#include "Randomize.hh"
#include "G4Exception.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>

namespace
{
  constexpr G4double kBoundaryEpsilon = 1.0e-4 * um;
  constexpr G4double kNearZero = 1.0e-18;
  constexpr G4double kLargeClearance = 1.0e9 * um;

  G4double Uniform01FromHash(std::uint64_t key)
  {
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdULL;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53ULL;
    key ^= key >> 33;
    constexpr long double kScale =
        1.0L / static_cast<long double>(std::numeric_limits<std::uint64_t>::max());
    return static_cast<G4double>(static_cast<long double>(key) * kScale);
  }
}

StageDReentrySampler::StageDReentrySampler(const DetectorConstruction *detector,
                                           const AnalysisConfig *config)
    : fDetector(detector),
      fConfig(config),
      fSpheres(),
      fSphereIdsByPhase(),
      fSphereWeightCdfByPhase(),
      fHalfX(detector ? detector->GetPatchHalfXUm() * um : 0.0),
      fHalfY(detector ? detector->GetPatchHalfYUm() * um : 0.0),
      fHalfZ(detector ? detector->GetPatchHalfZUm() * um : 0.0),
      fMaxSphereRadius(0.0),
      fGridCellSize(0.0),
      fGridMinCorner(-fHalfX, -fHalfY, -fHalfZ),
      fGridNx(1),
      fGridNy(1),
      fGridNz(1),
      fGridCells(),
      fCandidateVisitStamps(),
      fCandidateVisitToken(1),
      fPortalHalfX(0.0),
      fPortalHalfY(0.0),
      fPortalHalfZ(0.0),
      fPortalMargin((config != nullptr && config->stageD_portal_margin_um > 0.0)
                        ? config->stageD_portal_margin_um * um
                        : 0.0),
      fPortalNu((config != nullptr) ? config->stageD_portal_nu : 128),
      fPortalNv((config != nullptr) ? config->stageD_portal_nv : 128),
      fClearanceBinEdges{
          (config != nullptr ? config->stageD_clearance_bin0_um : 0.02) * um,
          (config != nullptr ? config->stageD_clearance_bin1_um : 0.10) * um,
          (config != nullptr ? config->stageD_clearance_bin2_um : 0.50) * um},
      fMaxParticleReentryTrials((config != nullptr) ? config->stageD_max_particle_reentry_trials : 64),
      fMaxPortalFallbackLevel((config != nullptr) ? config->stageD_max_portal_fallback_level : 4),
      fMatrixPortals(),
      fPortalSummary()
{
  BuildSphereCache();
  BuildSpatialGrid();
  BuildMatrixPortalPool();
}

G4bool StageDReentrySampler::SampleReentry(const ReentryContext &ctx,
                                           G4ThreeVector &newPosition,
                                           ReentryDiagnostics &diag) const
{
  diag = ReentryDiagnostics{};
  diag.exitPhase = ctx.phase;
  diag.exitInsidePoint = ctx.exitInsidePoint;

  if (ctx.phase == DetectorConstruction::Phase::BN ||
      ctx.phase == DetectorConstruction::Phase::ZnS)
  {
    std::string mode =
        (fConfig != nullptr) ? fConfig->stageD_particle_reentry_mode : std::string();
    if (mode.empty())
    {
      mode = (fConfig != nullptr && fConfig->stageD_reentry_mode == "state_matched")
                 ? "sphere_q_mu"
                 : "same_phase_random";
    }

    if (mode == "sphere_q_mu")
    {
      return SampleParticleSphereQMuReentry(ctx, newPosition, diag);
    }

    if (mode == "same_phase_rho_over_R" || mode == "same_phase_random")
    {
      const MicroSphere *oldSphere = FindContainingSphereOnly(ctx.phase, ctx.exitInsidePoint);
      G4ThreeVector oldPosition = ctx.exitInsidePoint;
      if (oldSphere == nullptr)
      {
        oldSphere = FindContainingSphereOnly(ctx.phase, ctx.prePos);
        oldPosition = ctx.prePos;
      }
      if (oldSphere == nullptr || oldSphere->radius <= 0.0)
        return false;

      const G4double rho = MinimumImageDelta(oldPosition, oldSphere->center).mag();
      const G4double qUpper =
          std::max(0.0, 1.0 - (kBoundaryEpsilon / std::max(oldSphere->radius, kBoundaryEpsilon)));
      const G4double q =
          std::clamp(rho / std::max(oldSphere->radius, kBoundaryEpsilon), 0.0, qUpper);

      diag.particleQExit = q;
      diag.particleMuExit = 0.0;
      const G4bool ok = SampleParticleSphereQOnlyReentry(ctx, oldSphere, q, newPosition, diag);
      if (ok)
      {
        if (mode == "same_phase_rho_over_R")
        {
          diag.strategy = "particle_sphere_q_only";
          diag.fallbackLevel = "legacy_q_mode";
        }
        else
        {
          diag.strategy = "particle_same_phase_random_debug";
          diag.fallbackLevel = "legacy_random_mode";
        }
      }
      return ok;
    }

    G4Exception("StageDReentrySampler::SampleReentry",
                "BNZS_D_REENTRY_000", FatalException,
                ("Unsupported Stage D particle re-entry mode: " + mode).c_str());
    return false;
  }

  if (ctx.phase == DetectorConstruction::Phase::Matrix)
  {
    const std::string &matrixMode =
        (fConfig != nullptr) ? fConfig->stageD_matrix_reentry_mode : std::string();
    if (matrixMode == "clearance_binned_portal" || matrixMode.empty())
      return SampleMatrixClearanceBinnedPortalReentry(ctx, newPosition, diag);

    if (matrixMode == "random_matrix_debug")
      return SampleRandomMatrixDebugReentry(ctx, newPosition, diag);

    if (matrixMode == "distance_matched_matrix")
    {
      G4Exception("StageDReentrySampler::SampleReentry",
                  "BNZS_D_REENTRY_001", FatalException,
                  "Stage D matrix re-entry mode distance_matched_matrix is not implemented. "
                  "Use clearance_binned_portal or random_matrix_debug.");
      return false;
    }

    if (matrixMode == "random_matrix")
    {
      G4Exception("StageDReentrySampler::SampleReentry",
                  "BNZS_D_REENTRY_002", FatalException,
                  "Stage D matrix re-entry mode random_matrix was renamed to random_matrix_debug. "
                  "It is debug-only and cannot be the production default.");
      return false;
    }

    G4Exception("StageDReentrySampler::SampleReentry",
                "BNZS_D_REENTRY_003", FatalException,
                ("Unsupported Stage D matrix re-entry mode: " + matrixMode).c_str());
    return false;
  }

  return false;
}

G4bool StageDReentrySampler::InsideRveBox(const G4ThreeVector &point) const
{
  return std::abs(point.x()) <= fHalfX &&
         std::abs(point.y()) <= fHalfY &&
         std::abs(point.z()) <= fHalfZ;
}

DetectorConstruction::Phase StageDReentrySampler::FastPhaseAtPointForReentry(
    const G4ThreeVector &point) const
{
  if (!InsideRveBox(point))
    return DetectorConstruction::Phase::World;
  return fDetector != nullptr ? fDetector->FindPhaseAtPoint(point)
                              : DetectorConstruction::Phase::Unknown;
}

void StageDReentrySampler::BuildSphereCache()
{
  if (fDetector == nullptr)
    return;

  auto appendPhase = [&](const std::vector<DetectorConstruction::SphereInfo> &src,
                         DetectorConstruction::Phase phase)
  {
    const std::size_t phaseIndex = PhaseVectorIndex(phase);
    for (const auto &sphere : src)
    {
      MicroSphere micro;
      micro.phase = phase;
      micro.center = sphere.center;
      micro.radius = sphere.radius;
      micro.radius2 = sphere.radius * sphere.radius;
      micro.visibleVolume = (4.0 / 3.0) * CLHEP::pi * sphere.radius * sphere.radius * sphere.radius;
      micro.isClipped =
          (std::abs(micro.center.x()) + micro.radius > fHalfX ||
           std::abs(micro.center.y()) + micro.radius > fHalfY ||
           std::abs(micro.center.z()) + micro.radius > fHalfZ);

      const G4int sphereId = static_cast<G4int>(fSpheres.size());
      fSpheres.push_back(micro);
      fSphereIdsByPhase[phaseIndex].push_back(sphereId);
      fMaxSphereRadius = std::max(fMaxSphereRadius, sphere.radius);
    }
  };

  appendPhase(fDetector->GetBNSpheres(), DetectorConstruction::Phase::BN);
  appendPhase(fDetector->GetZnSSpheres(), DetectorConstruction::Phase::ZnS);

  for (std::size_t phaseIndex = 0; phaseIndex < fSphereIdsByPhase.size(); ++phaseIndex)
  {
    G4double running = 0.0;
    auto &cdf = fSphereWeightCdfByPhase[phaseIndex];
    cdf.reserve(fSphereIdsByPhase[phaseIndex].size());
    for (G4int sphereId : fSphereIdsByPhase[phaseIndex])
    {
      running += SphereSelectionWeight(fSpheres[static_cast<std::size_t>(sphereId)]);
      cdf.push_back(running);
    }
  }

  fCandidateVisitStamps.assign(fSpheres.size(), 0);
}

void StageDReentrySampler::BuildSpatialGrid()
{
  fGridCellSize = std::max(fMaxSphereRadius, 2.0 * kBoundaryEpsilon);
  if (fGridCellSize <= 0.0)
  {
    fGridCellSize = std::max({fHalfX, fHalfY, fHalfZ, 1.0 * um});
  }

  const auto computeDim = [&](G4double halfExtent)
  {
    return std::max(1, static_cast<G4int>(std::floor((2.0 * halfExtent) / fGridCellSize)));
  };

  fGridNx = computeDim(fHalfX);
  fGridNy = computeDim(fHalfY);
  fGridNz = computeDim(fHalfZ);

  fGridCells.clear();
  fGridCells.resize(static_cast<std::size_t>(fGridNx) *
                    static_cast<std::size_t>(fGridNy) *
                    static_cast<std::size_t>(fGridNz));

  for (std::size_t sphereId = 0; sphereId < fSpheres.size(); ++sphereId)
  {
    const MicroSphere &sphere = fSpheres[sphereId];
    auto coordToCell = [&](G4double coord, G4double halfExtent, G4int dim)
    {
      const G4double normalized = (coord + halfExtent) / (2.0 * halfExtent);
      return std::clamp(static_cast<G4int>(std::floor(normalized * dim)), 0, dim - 1);
    };
    const G4int ix = coordToCell(sphere.center.x(), fHalfX, fGridNx);
    const G4int iy = coordToCell(sphere.center.y(), fHalfY, fGridNy);
    const G4int iz = coordToCell(sphere.center.z(), fHalfZ, fGridNz);
    fGridCells[FlatCellIndex(ix, iy, iz)].push_back(static_cast<G4int>(sphereId));
  }
}

void StageDReentrySampler::BuildMatrixPortalPool()
{
  fPortalSummary = StageDReentryPortalSummary{};

  const G4double defaultMargin = (fPortalMargin > 0.0) ? fPortalMargin : fMaxSphereRadius;
  if (fHalfX > defaultMargin && fHalfY > defaultMargin && fHalfZ > defaultMargin)
  {
    fPortalHalfX = fHalfX - defaultMargin;
    fPortalHalfY = fHalfY - defaultMargin;
    fPortalHalfZ = fHalfZ - defaultMargin;
    fPortalMargin = defaultMargin;
  }
  else
  {
    fPortalHalfX = std::max(kBoundaryEpsilon, fHalfX - kBoundaryEpsilon);
    fPortalHalfY = std::max(kBoundaryEpsilon, fHalfY - kBoundaryEpsilon);
    fPortalHalfZ = std::max(kBoundaryEpsilon, fHalfZ - kBoundaryEpsilon);
    fPortalMargin = kBoundaryEpsilon;
  }

  fPortalNu = std::max(1, fPortalNu);
  fPortalNv = std::max(1, fPortalNv);

  for (std::size_t faceIndex = 0; faceIndex < kPortalFaceCount; ++faceIndex)
  {
    const PortalFace face = static_cast<PortalFace>(faceIndex);
    const G4ThreeVector inward = InwardNormal(face);
    for (G4int iu = 0; iu < fPortalNu; ++iu)
    {
      for (G4int iv = 0; iv < fPortalNv; ++iv)
      {
        const std::uint64_t base =
            (static_cast<std::uint64_t>(faceIndex) << 48) ^
            (static_cast<std::uint64_t>(iu) << 24) ^
            static_cast<std::uint64_t>(iv);
        const G4double u01 = (static_cast<G4double>(iu) + Uniform01FromHash(base ^ 0x9e3779b97f4a7c15ULL)) /
                             static_cast<G4double>(fPortalNu);
        const G4double v01 = (static_cast<G4double>(iv) + Uniform01FromHash(base ^ 0xc2b2ae3d27d4eb4fULL)) /
                             static_cast<G4double>(fPortalNv);

        const G4ThreeVector surfacePoint = PointOnVirtualFace(face, u01, v01);
        const G4ThreeVector insidePoint = surfacePoint + kBoundaryEpsilon * inward;
        if (FastPhaseAtPointForReentry(insidePoint) != DetectorConstruction::Phase::Matrix)
          continue;

        const NearestSurface nearest = FindNearestParticleSurface(insidePoint);
        const G4int bin = ClearanceBin(nearest.clearance);
        const std::size_t bucket = PortalPhaseBucketIndex(nearest.nearestPhase);

        MatrixPortal portal;
        portal.position = insidePoint;
        portal.clearance = nearest.clearance;
        portal.nearestPhase = nearest.nearestPhase;
        portal.clearanceBin = bin;

        fMatrixPortals[faceIndex][bucket][static_cast<std::size_t>(bin)].push_back(portal);
        ++fPortalSummary.total_portal_count;
        ++fPortalSummary.portal_count_by_face[faceIndex];
        ++fPortalSummary.portal_count_by_bin[static_cast<std::size_t>(bin)];
      }
    }
  }
}

G4bool StageDReentrySampler::SampleParticleSphereQMuReentry(
    const ReentryContext &ctx,
    G4ThreeVector &newPosition,
    ReentryDiagnostics &diag) const
{
  diag.strategy = "particle_sphere_q_mu";
  diag.fallbackLevel = "none";

  const MicroSphere *oldSphere = FindContainingSphereOnly(ctx.phase, ctx.exitInsidePoint);
  G4ThreeVector oldPosition = ctx.exitInsidePoint;
  if (oldSphere == nullptr)
  {
    oldSphere = FindContainingSphereOnly(ctx.phase, ctx.prePos);
    oldPosition = ctx.prePos;
  }
  if (oldSphere == nullptr || oldSphere->radius <= 0.0)
    return false;

  const G4ThreeVector oldOffset = MinimumImageDelta(oldPosition, oldSphere->center);
  const G4double rho = oldOffset.mag();
  const G4double qUpper = std::max(0.0, 1.0 - (kBoundaryEpsilon / std::max(oldSphere->radius, kBoundaryEpsilon)));
  const G4double q = std::clamp(rho / std::max(oldSphere->radius, kBoundaryEpsilon), 0.0, qUpper);

  const G4ThreeVector oldDir = ctx.oldDir.unit();
  G4double mu = 0.0;
  if (rho > kNearZero)
  {
    const G4ThreeVector nOld = oldOffset / rho;
    mu = std::clamp(oldDir.dot(nOld), -1.0, 1.0);
  }

  diag.particleQExit = q;
  diag.particleMuExit = mu;

  const G4int maxTrials = std::max(1, fMaxParticleReentryTrials);
  for (G4int trial = 0; trial < maxTrials; ++trial)
  {
    ++diag.trials;
    const G4int sphereId = SampleWeightedSamePhaseSphereId(ctx.phase);
    if (sphereId < 0)
      break;

    const MicroSphere &newSphere = fSpheres[static_cast<std::size_t>(sphereId)];
    const G4ThreeVector nNew = RandomUnitVectorWithFixedDot(oldDir, mu);
    G4ThreeVector candidate = newSphere.center + q * newSphere.radius * nNew;

    if (!ClampInsideSameSphereRoundoffOnly(candidate, newSphere))
      continue;
    const G4double rhoEntry = (candidate - newSphere.center).mag();
    candidate = fDetector->WrapToPrimaryCell(candidate);
    if (!ValidatePhase(candidate, ctx.phase))
      continue;

    G4ThreeVector nEntry = nNew;

    newPosition = candidate;
    diag.entryPoint = candidate;
    diag.entryPhase = ctx.phase;
    diag.particleQEntry = rhoEntry / std::max(newSphere.radius, kBoundaryEpsilon);
    diag.particleMuEntry = std::clamp(oldDir.dot(nEntry), -1.0, 1.0);
    return true;
  }

  return SampleParticleSphereQOnlyReentry(ctx, oldSphere, q, newPosition, diag);
}

G4bool StageDReentrySampler::SampleParticleSphereQOnlyReentry(
    const ReentryContext &ctx,
    const MicroSphere *oldSphere,
    G4double q,
    G4ThreeVector &newPosition,
    ReentryDiagnostics &diag) const
{
  (void)oldSphere;

  const G4int maxTrials = std::max(1, fMaxParticleReentryTrials / 2);
  for (G4int trial = 0; trial < maxTrials; ++trial)
  {
    ++diag.trials;
    const G4int sphereId = SampleWeightedSamePhaseSphereId(ctx.phase);
    if (sphereId < 0)
      break;

    const MicroSphere &newSphere = fSpheres[static_cast<std::size_t>(sphereId)];
    const G4ThreeVector nNew = RandomUnitVector();
    G4ThreeVector candidate = newSphere.center + q * newSphere.radius * nNew;

    if (!ClampInsideSameSphereRoundoffOnly(candidate, newSphere))
      continue;
    const G4double rhoEntry = (candidate - newSphere.center).mag();
    const G4ThreeVector radial = (rhoEntry > kNearZero)
                                     ? ((candidate - newSphere.center) / rhoEntry)
                                     : nNew;
    candidate = fDetector->WrapToPrimaryCell(candidate);
    if (!ValidatePhase(candidate, ctx.phase))
      continue;

    newPosition = candidate;
    diag.strategy = "particle_sphere_q_only";
    diag.fallbackLevel = "q_only";
    diag.entryPoint = candidate;
    diag.entryPhase = ctx.phase;
    diag.particleQEntry = rhoEntry / std::max(newSphere.radius, kBoundaryEpsilon);
    diag.particleMuEntry = std::clamp(ctx.oldDir.unit().dot(radial), -1.0, 1.0);
    return true;
  }

  return false;
}

G4bool StageDReentrySampler::SampleMatrixClearanceBinnedPortalReentry(
    const ReentryContext &ctx,
    G4ThreeVector &newPosition,
    ReentryDiagnostics &diag) const
{
  diag.strategy = "matrix_clearance_binned_portal";
  diag.fallbackLevel = "none";

  const NearestSurface exitNearest = FindNearestParticleSurface(ctx.exitInsidePoint);
  const G4int exitBin = ClearanceBin(exitNearest.clearance);
  const DetectorConstruction::Phase exitNearestPhase = exitNearest.nearestPhase;
  diag.matrixClearanceExitUm = exitNearest.clearance / um;
  diag.matrixNearestPhaseExit = PhaseNameOrUnknown(exitNearestPhase);
  diag.matrixClearanceBinExit = exitBin;

  struct PoolRef
  {
    PortalFace face;
    const std::vector<MatrixPortal> *pool;
  };

  const G4ThreeVector oldDir = ctx.oldDir.unit();

  auto isInflowFace = [&](PortalFace face)
  {
    return oldDir.dot(InwardNormal(face)) > 0.0;
  };

  auto addPoolRef = [&](std::vector<PoolRef> &refs,
                        PortalFace face,
                        std::size_t bucketIndex,
                        std::size_t binIndex)
  {
    const auto &pool =
        fMatrixPortals[static_cast<std::size_t>(face)][bucketIndex][binIndex];
    if (!pool.empty())
      refs.push_back(PoolRef{face, &pool});
  };

  auto attemptPools = [&](const std::vector<std::size_t> &bucketIndices,
                          const std::vector<std::size_t> &binIndices,
                          const std::string &fallbackLabel) -> G4bool
  {
    struct FaceChoice
    {
      PortalFace face;
      G4double weight = 0.0;
      std::vector<PoolRef> pools;
    };

    std::vector<FaceChoice> faceChoices;
    for (std::size_t faceIndex = 0; faceIndex < kPortalFaceCount; ++faceIndex)
    {
      const PortalFace face = static_cast<PortalFace>(faceIndex);
      if (!isInflowFace(face))
        continue;

      FaceChoice choice;
      choice.face = face;
      for (std::size_t bucketIndex : bucketIndices)
      {
        for (std::size_t binIndex : binIndices)
        {
          addPoolRef(choice.pools, face, bucketIndex, binIndex);
        }
      }

      if (choice.pools.empty())
        continue;

      std::size_t totalPoolSize = 0;
      for (const PoolRef &pool : choice.pools)
        totalPoolSize += pool.pool->size();

      if (totalPoolSize == 0)
        continue;

      choice.weight = oldDir.dot(InwardNormal(face)) *
                      FaceArea(face) *
                      static_cast<G4double>(totalPoolSize);
      if (choice.weight > 0.0)
        faceChoices.push_back(choice);
    }

    if (faceChoices.empty())
      return false;

    G4double totalWeight = 0.0;
    for (const auto &choice : faceChoices)
      totalWeight += choice.weight;
    if (totalWeight <= 0.0)
      return false;

    G4double pick = G4UniformRand() * totalWeight;
    const FaceChoice *selectedFace = &faceChoices.front();
    for (const auto &choice : faceChoices)
    {
      pick -= choice.weight;
      if (pick <= 0.0)
      {
        selectedFace = &choice;
        break;
      }
    }

    std::size_t totalPortalCount = 0;
    for (const PoolRef &poolRef : selectedFace->pools)
      totalPortalCount += poolRef.pool->size();
    if (totalPortalCount == 0)
      return false;

    std::size_t portalPick =
        static_cast<std::size_t>(G4UniformRand() * static_cast<G4double>(totalPortalCount));
    if (portalPick >= totalPortalCount)
      portalPick = totalPortalCount - 1;

    const MatrixPortal *selectedPortal = nullptr;
    for (const PoolRef &poolRef : selectedFace->pools)
    {
      if (portalPick < poolRef.pool->size())
      {
        selectedPortal = &(*poolRef.pool)[portalPick];
        break;
      }
      portalPick -= poolRef.pool->size();
    }
    if (selectedPortal == nullptr)
      return false;

    if (!ValidatePhase(selectedPortal->position, DetectorConstruction::Phase::Matrix))
      return false;

    const NearestSurface entryNearest = FindNearestParticleSurface(selectedPortal->position);
    const G4int entryBin = ClearanceBin(entryNearest.clearance);

    newPosition = selectedPortal->position;
    diag.fallbackLevel = fallbackLabel;
    diag.entryPoint = selectedPortal->position;
    diag.entryPhase = DetectorConstruction::Phase::Matrix;
    diag.matrixClearanceEntryUm = entryNearest.clearance / um;
    diag.matrixNearestPhaseEntry = PhaseNameOrUnknown(entryNearest.nearestPhase);
    diag.matrixClearanceBinEntry = entryBin;
    return true;
  };

  const std::size_t sameBucket = PortalPhaseBucketIndex(exitNearestPhase);
  const std::array<std::size_t, 3> anyBuckets{0, 1, 2};
  const std::array<std::size_t, 1> sameBucketOnly{sameBucket};
  const std::array<std::size_t, 1> sameBinOnly{static_cast<std::size_t>(exitBin)};

  if (attemptPools(std::vector<std::size_t>(sameBucketOnly.begin(), sameBucketOnly.end()),
                   std::vector<std::size_t>(sameBinOnly.begin(), sameBinOnly.end()),
                   "same_bin"))
  {
    return true;
  }

  if (fMaxPortalFallbackLevel >= 1)
  {
    std::vector<std::size_t> adjacentBins;
    if (exitBin > 0)
      adjacentBins.push_back(static_cast<std::size_t>(exitBin - 1));
    if (exitBin + 1 < static_cast<G4int>(kClearanceBinCount))
      adjacentBins.push_back(static_cast<std::size_t>(exitBin + 1));

    if (!adjacentBins.empty() &&
        attemptPools(std::vector<std::size_t>(sameBucketOnly.begin(), sameBucketOnly.end()),
                     adjacentBins,
                     "adjacent_bin"))
    {
      return true;
    }
  }

  if (fMaxPortalFallbackLevel >= 2)
  {
    std::vector<std::size_t> allBins{0, 1, 2, 3};
    if (attemptPools(std::vector<std::size_t>(sameBucketOnly.begin(), sameBucketOnly.end()),
                     allBins,
                     "same_phase_any_bin"))
    {
      return true;
    }
  }

  if (fMaxPortalFallbackLevel >= 3)
  {
    if (attemptPools(std::vector<std::size_t>(anyBuckets.begin(), anyBuckets.end()),
                     std::vector<std::size_t>(sameBinOnly.begin(), sameBinOnly.end()),
                     "any_phase_same_bin"))
    {
      return true;
    }
  }

  if (fMaxPortalFallbackLevel >= 4)
  {
    std::vector<std::size_t> allBins{0, 1, 2, 3};
    if (attemptPools(std::vector<std::size_t>(anyBuckets.begin(), anyBuckets.end()),
                     allBins,
                     "any_portal"))
    {
      return true;
    }
  }

  return false;
}

G4bool StageDReentrySampler::SampleRandomMatrixDebugReentry(
    const ReentryContext &ctx,
    G4ThreeVector &newPosition,
    ReentryDiagnostics &diag) const
{
  (void)ctx;

  constexpr G4int kMaxRandomMatrixTrials = 10000;
  for (G4int trial = 0; trial < kMaxRandomMatrixTrials; ++trial)
  {
    ++diag.trials;
    const G4ThreeVector candidate = RandomPointInMatrixBox();
    if (!ValidatePhase(candidate, DetectorConstruction::Phase::Matrix))
      continue;

    const NearestSurface nearest = FindNearestParticleSurface(candidate);
    newPosition = candidate;
    diag.strategy = "matrix_random_debug";
    diag.fallbackLevel = "random_matrix_debug";
    diag.entryPoint = candidate;
    diag.entryPhase = DetectorConstruction::Phase::Matrix;
    diag.matrixClearanceEntryUm = nearest.clearance / um;
    diag.matrixNearestPhaseEntry = PhaseNameOrUnknown(nearest.nearestPhase);
    diag.matrixClearanceBinEntry = ClearanceBin(nearest.clearance);
    return true;
  }

  return false;
}

const StageDReentrySampler::MicroSphere *StageDReentrySampler::FindContainingSphereOnly(
    DetectorConstruction::Phase phase,
    const G4ThreeVector &position) const
{
  if (!InsideRveBox(position))
    return nullptr;

  std::vector<G4int> candidateIds;
  CollectCandidateSphereIds(position, candidateIds);

  for (G4int sphereId : candidateIds)
  {
    const MicroSphere &sphere = fSpheres[static_cast<std::size_t>(sphereId)];
    if (sphere.phase != phase)
      continue;
    if (MinimumImageDelta(position, sphere.center).mag2() <= sphere.radius2)
      return &sphere;
  }

  return nullptr;
}

StageDReentrySampler::NearestSurface StageDReentrySampler::FindNearestParticleSurface(
    const G4ThreeVector &position) const
{
  NearestSurface best;
  best.clearance = kLargeClearance;

  G4int ix = 0;
  G4int iy = 0;
  G4int iz = 0;
  if (!PointToCell(position, ix, iy, iz))
    return best;

  if (fSpheres.empty())
    return best;

  if (fCandidateVisitToken == std::numeric_limits<G4int>::max())
  {
    std::fill(fCandidateVisitStamps.begin(), fCandidateVisitStamps.end(), 0);
    fCandidateVisitToken = 1;
  }
  ++fCandidateVisitToken;

  const G4double maxSearchDistance =
      fClearanceBinEdges[2] + fMaxSphereRadius + fGridCellSize;
  const G4int maxLayer =
      std::max(1, static_cast<G4int>(std::ceil(maxSearchDistance / fGridCellSize)));

  for (G4int layer = 0; layer <= maxLayer; ++layer)
  {
    const auto modulo = [](G4int value, G4int count)
    {
      const G4int result = value % count;
      return result < 0 ? result + count : result;
    };
    for (G4int dx = -layer; dx <= layer; ++dx)
    {
      const G4int cx = modulo(ix + dx, fGridNx);
      for (G4int dy = -layer; dy <= layer; ++dy)
      {
        const G4int cy = modulo(iy + dy, fGridNy);
        for (G4int dz = -layer; dz <= layer; ++dz)
        {
          const G4int cz = modulo(iz + dz, fGridNz);
          const auto &cell = fGridCells[FlatCellIndex(cx, cy, cz)];
          for (G4int sphereId : cell)
          {
            const std::size_t sphereIndex = static_cast<std::size_t>(sphereId);
            if (fCandidateVisitStamps[sphereIndex] == fCandidateVisitToken)
              continue;
            fCandidateVisitStamps[sphereIndex] = fCandidateVisitToken;

            const MicroSphere &sphere = fSpheres[sphereIndex];
            const G4ThreeVector offset = MinimumImageDelta(position, sphere.center);
            const G4double distance = offset.mag();
            const G4double clearance = distance - sphere.radius;
            if (clearance < best.clearance)
            {
              best.clearance = clearance;
              best.nearestPhase = sphere.phase;
              best.sphereId = sphereId;
              best.surfaceNormal = (distance > kNearZero) ? (offset / distance) : RandomUnitVector();
            }
          }
        }
      }
    }
  }

  if (best.clearance > fClearanceBinEdges[2])
  {
    best.nearestPhase = DetectorConstruction::Phase::Unknown;
    best.sphereId = -1;
  }

  return best;
}

G4bool StageDReentrySampler::ClampInsideSameSphereRoundoffOnly(
    G4ThreeVector &point,
    const MicroSphere &sphere) const
{
  G4ThreeVector offset = point - sphere.center;
  const G4double rho = offset.mag();
  const G4double maxInsideRadius = std::max(0.0, sphere.radius - kBoundaryEpsilon);

  if (rho <= maxInsideRadius)
    return true;

  if (rho > sphere.radius + 8.0 * kBoundaryEpsilon)
    return false;

  if (rho <= kNearZero)
  {
    point = sphere.center;
    return true;
  }

  offset /= rho;
  point = sphere.center + maxInsideRadius * offset;
  return true;
}

G4bool StageDReentrySampler::ValidatePhase(const G4ThreeVector &point,
                                           DetectorConstruction::Phase phase) const
{
  if (phase == DetectorConstruction::Phase::Matrix)
  {
    return InsideRveBox(point) &&
           FastPhaseAtPointForReentry(point) == DetectorConstruction::Phase::Matrix;
  }

  if (phase == DetectorConstruction::Phase::BN ||
      phase == DetectorConstruction::Phase::ZnS)
  {
    return InsideRveBox(point) &&
           FastPhaseAtPointForReentry(point) == phase;
  }

  return false;
}

G4ThreeVector StageDReentrySampler::RandomUnitVector() const
{
  const G4double cosTheta = 2.0 * G4UniformRand() - 1.0;
  const G4double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
  const G4double phi = CLHEP::twopi * G4UniformRand();
  return G4ThreeVector(sinTheta * std::cos(phi),
                       sinTheta * std::sin(phi),
                       cosTheta);
}

G4ThreeVector StageDReentrySampler::RandomUnitVectorWithFixedDot(
    const G4ThreeVector &axis,
    G4double mu) const
{
  if (axis.mag2() <= kNearZero)
    return RandomUnitVector();

  const G4ThreeVector e3 = axis.unit();
  G4ThreeVector ref(0.0, 0.0, 1.0);
  if (std::abs(e3.dot(ref)) > 0.9)
    ref = G4ThreeVector(0.0, 1.0, 0.0);

  G4ThreeVector e1 = e3.cross(ref);
  if (e1.mag2() <= kNearZero)
    return RandomUnitVector();
  e1 = e1.unit();
  G4ThreeVector e2 = e3.cross(e1).unit();

  const G4double phi = CLHEP::twopi * G4UniformRand();
  const G4double muClamped = std::clamp(mu, -1.0, 1.0);
  const G4double s = std::sqrt(std::max(0.0, 1.0 - muClamped * muClamped));
  return muClamped * e3 + s * (std::cos(phi) * e1 + std::sin(phi) * e2);
}

G4ThreeVector StageDReentrySampler::RandomPointInMatrixBox() const
{
  return G4ThreeVector((2.0 * G4UniformRand() - 1.0) * fHalfX,
                       (2.0 * G4UniformRand() - 1.0) * fHalfY,
                       (2.0 * G4UniformRand() - 1.0) * fHalfZ);
}

G4double StageDReentrySampler::SphereSelectionWeight(const MicroSphere &sphere) const
{
  if (sphere.visibleVolume > 0.0)
    return sphere.visibleVolume;
  return std::max(kNearZero, sphere.radius * sphere.radius * sphere.radius);
}

G4int StageDReentrySampler::SampleWeightedSamePhaseSphereId(
    DetectorConstruction::Phase phase) const
{
  const std::size_t phaseIndex = PhaseVectorIndex(phase);
  const auto &sphereIds = fSphereIdsByPhase[phaseIndex];
  const auto &cdf = fSphereWeightCdfByPhase[phaseIndex];
  if (sphereIds.empty() || cdf.empty() || cdf.back() <= 0.0)
    return -1;

  const G4double target = G4UniformRand() * cdf.back();
  const auto it = std::lower_bound(cdf.begin(), cdf.end(), target);
  const std::size_t index = static_cast<std::size_t>(std::distance(cdf.begin(), it));
  return sphereIds[std::min(index, sphereIds.size() - 1)];
}

G4int StageDReentrySampler::ClearanceBin(G4double clearance) const
{
  if (clearance < fClearanceBinEdges[0])
    return 0;
  if (clearance < fClearanceBinEdges[1])
    return 1;
  if (clearance < fClearanceBinEdges[2])
    return 2;
  return 3;
}

std::size_t StageDReentrySampler::FlatCellIndex(G4int ix, G4int iy, G4int iz) const
{
  return (static_cast<std::size_t>(ix) * static_cast<std::size_t>(fGridNy) +
          static_cast<std::size_t>(iy)) *
             static_cast<std::size_t>(fGridNz) +
         static_cast<std::size_t>(iz);
}

G4bool StageDReentrySampler::PointToCell(const G4ThreeVector &point,
                                         G4int &ix,
                                         G4int &iy,
                                         G4int &iz) const
{
  if (!InsideRveBox(point))
    return false;

  auto coordToCell = [&](G4double coord, G4double halfExtent, G4int dim)
  {
    const G4double normalized = (coord + halfExtent) / (2.0 * halfExtent);
    return std::clamp(static_cast<G4int>(std::floor(normalized * dim)), 0, dim - 1);
  };

  ix = coordToCell(point.x(), fHalfX, fGridNx);
  iy = coordToCell(point.y(), fHalfY, fGridNy);
  iz = coordToCell(point.z(), fHalfZ, fGridNz);
  return true;
}

void StageDReentrySampler::CollectCandidateSphereIds(
    const G4ThreeVector &point,
    std::vector<G4int> &candidateIds) const
{
  candidateIds.clear();

  G4int ix = 0;
  G4int iy = 0;
  G4int iz = 0;
  if (!PointToCell(point, ix, iy, iz))
    return;

  if (fSpheres.empty())
    return;

  if (fCandidateVisitToken == std::numeric_limits<G4int>::max())
  {
    std::fill(fCandidateVisitStamps.begin(), fCandidateVisitStamps.end(), 0);
    fCandidateVisitToken = 1;
  }
  ++fCandidateVisitToken;

  const auto modulo = [](G4int value, G4int count)
  {
    const G4int result = value % count;
    return result < 0 ? result + count : result;
  };

  for (G4int dx = -1; dx <= 1; ++dx)
  {
    const G4int cx = modulo(ix + dx, fGridNx);
    for (G4int dy = -1; dy <= 1; ++dy)
    {
      const G4int cy = modulo(iy + dy, fGridNy);
      for (G4int dz = -1; dz <= 1; ++dz)
      {
        const G4int cz = modulo(iz + dz, fGridNz);

        const auto &cell = fGridCells[FlatCellIndex(cx, cy, cz)];
        for (G4int sphereId : cell)
        {
          const std::size_t sphereIndex = static_cast<std::size_t>(sphereId);
          if (fCandidateVisitStamps[sphereIndex] == fCandidateVisitToken)
            continue;
          fCandidateVisitStamps[sphereIndex] = fCandidateVisitToken;
          candidateIds.push_back(sphereId);
        }
      }
    }
  }
}

G4ThreeVector StageDReentrySampler::MinimumImageDelta(
    const G4ThreeVector &point,
    const G4ThreeVector &center) const
{
  G4ThreeVector delta = point - center;
  const G4double lengthX = 2.0 * fHalfX;
  const G4double lengthY = 2.0 * fHalfY;
  const G4double lengthZ = 2.0 * fHalfZ;
  delta.setX(delta.x() - lengthX * std::nearbyint(delta.x() / lengthX));
  delta.setY(delta.y() - lengthY * std::nearbyint(delta.y() / lengthY));
  delta.setZ(delta.z() - lengthZ * std::nearbyint(delta.z() / lengthZ));
  return delta;
}

G4double StageDReentrySampler::HalfExtentX() const
{
  return fPortalHalfX;
}

G4double StageDReentrySampler::HalfExtentY() const
{
  return fPortalHalfY;
}

G4double StageDReentrySampler::HalfExtentZ() const
{
  return fPortalHalfZ;
}

G4ThreeVector StageDReentrySampler::InwardNormal(PortalFace face) const
{
  switch (face)
  {
  case PortalFace::PosX:
    return G4ThreeVector(-1.0, 0.0, 0.0);
  case PortalFace::NegX:
    return G4ThreeVector(+1.0, 0.0, 0.0);
  case PortalFace::PosY:
    return G4ThreeVector(0.0, -1.0, 0.0);
  case PortalFace::NegY:
    return G4ThreeVector(0.0, +1.0, 0.0);
  case PortalFace::PosZ:
    return G4ThreeVector(0.0, 0.0, -1.0);
  case PortalFace::NegZ:
  default:
    return G4ThreeVector(0.0, 0.0, +1.0);
  }
}

G4double StageDReentrySampler::FaceArea(PortalFace face) const
{
  switch (face)
  {
  case PortalFace::PosX:
  case PortalFace::NegX:
    return 4.0 * HalfExtentY() * HalfExtentZ();
  case PortalFace::PosY:
  case PortalFace::NegY:
    return 4.0 * HalfExtentX() * HalfExtentZ();
  case PortalFace::PosZ:
  case PortalFace::NegZ:
  default:
    return 4.0 * HalfExtentX() * HalfExtentY();
  }
}

G4ThreeVector StageDReentrySampler::PointOnVirtualFace(PortalFace face,
                                                       G4double u01,
                                                       G4double v01) const
{
  const G4double x = -HalfExtentX() + 2.0 * HalfExtentX() * u01;
  const G4double y = -HalfExtentY() + 2.0 * HalfExtentY() * u01;
  const G4double z = -HalfExtentZ() + 2.0 * HalfExtentZ() * v01;

  switch (face)
  {
  case PortalFace::PosX:
    return G4ThreeVector(+HalfExtentX(), -HalfExtentY() + 2.0 * HalfExtentY() * u01, -HalfExtentZ() + 2.0 * HalfExtentZ() * v01);
  case PortalFace::NegX:
    return G4ThreeVector(-HalfExtentX(), -HalfExtentY() + 2.0 * HalfExtentY() * u01, -HalfExtentZ() + 2.0 * HalfExtentZ() * v01);
  case PortalFace::PosY:
    return G4ThreeVector(-HalfExtentX() + 2.0 * HalfExtentX() * u01, +HalfExtentY(), -HalfExtentZ() + 2.0 * HalfExtentZ() * v01);
  case PortalFace::NegY:
    return G4ThreeVector(-HalfExtentX() + 2.0 * HalfExtentX() * u01, -HalfExtentY(), -HalfExtentZ() + 2.0 * HalfExtentZ() * v01);
  case PortalFace::PosZ:
    return G4ThreeVector(-HalfExtentX() + 2.0 * HalfExtentX() * u01, -HalfExtentY() + 2.0 * HalfExtentY() * v01, +HalfExtentZ());
  case PortalFace::NegZ:
  default:
    return G4ThreeVector(-HalfExtentX() + 2.0 * HalfExtentX() * u01, -HalfExtentY() + 2.0 * HalfExtentY() * v01, -HalfExtentZ());
  }
}

std::size_t StageDReentrySampler::PhaseVectorIndex(DetectorConstruction::Phase phase)
{
  return (phase == DetectorConstruction::Phase::ZnS) ? 1u : 0u;
}

std::size_t StageDReentrySampler::PortalPhaseBucketIndex(
    DetectorConstruction::Phase phase)
{
  if (phase == DetectorConstruction::Phase::BN)
    return 0u;
  if (phase == DetectorConstruction::Phase::ZnS)
    return 1u;
  return 2u;
}

const char *StageDReentrySampler::PhaseNameOrUnknown(
    DetectorConstruction::Phase phase)
{
  if (phase == DetectorConstruction::Phase::Unknown)
    return "Unknown";
  return DetectorConstruction::PhaseName(phase);
}
