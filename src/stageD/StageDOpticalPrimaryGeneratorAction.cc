#include "StageDOpticalPrimaryGeneratorAction.hh"

#include "AnalysisConfig.hh"
#include "DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4Exception.hh"
#include "G4OpticalPhoton.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cmath>

namespace
{
  G4double EnergyFromWavelengthNm(G4double wavelengthNm)
  {
    return (1239.841984 / wavelengthNm) * eV;
  }

  const char *PhaseToString(DetectorConstruction::Phase phase)
  {
    return DetectorConstruction::PhaseName(phase);
  }
}

StageDOpticalPrimaryGeneratorAction::StageDOpticalPrimaryGeneratorAction(
    AnalysisConfig *config)
    : G4VUserPrimaryGeneratorAction(),
      fConfig(config),
      fParticleGun(new G4ParticleGun(1)),
      fCurrentPhoton(),
      fZnSCdfDetector(nullptr),
      fZnSCdfPlacementSeed(0),
      fZnSCdfSphereCount(0),
      fZnSSphereCdf(),
      fZnSSphereTotalWeight(0.0)
{
  fParticleGun->SetParticleDefinition(
      G4OpticalPhoton::OpticalPhotonDefinition());
  fParticleGun->SetParticleTime(0.0);
}

StageDOpticalPrimaryGeneratorAction::~StageDOpticalPrimaryGeneratorAction()
{
  delete fParticleGun;
}

const DetectorConstruction *StageDOpticalPrimaryGeneratorAction::ResolveDetector() const
{
  return dynamic_cast<const DetectorConstruction *>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
}

G4ThreeVector StageDOpticalPrimaryGeneratorAction::RandomUnitVector() const
{
  const G4double cosTheta = 2.0 * G4UniformRand() - 1.0;
  const G4double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
  const G4double phi = CLHEP::twopi * G4UniformRand();
  return G4ThreeVector(sinTheta * std::cos(phi),
                       sinTheta * std::sin(phi),
                       cosTheta);
}

G4ThreeVector StageDOpticalPrimaryGeneratorAction::RandomPolarization(
    const G4ThreeVector &direction) const
{
  G4ThreeVector reference(0.0, 0.0, 1.0);
  if (std::abs(direction.dot(reference)) > 0.9)
    reference = G4ThreeVector(0.0, 1.0, 0.0);

  G4ThreeVector e1 = direction.cross(reference).unit();
  G4ThreeVector e2 = direction.cross(e1).unit();
  const G4double phi = CLHEP::twopi * G4UniformRand();
  return std::cos(phi) * e1 + std::sin(phi) * e2;
}

G4ThreeVector StageDOpticalPrimaryGeneratorAction::SampleUniformPointInZnSSphere() const
{
  const auto *detector = ResolveDetector();
  if (detector == nullptr)
  {
    G4Exception("StageDOpticalPrimaryGeneratorAction::SampleUniformPointInZnSSphere",
                "BNZS_D_PRI_001", FatalException,
                "DetectorConstruction is null.");
    return G4ThreeVector();
  }

  const auto &spheres = detector->GetZnSSpheres();
  if (spheres.empty())
  {
    G4Exception("StageDOpticalPrimaryGeneratorAction::SampleUniformPointInZnSSphere",
                "BNZS_D_PRI_002", FatalException,
                "No ZnS spheres are available for Stage D uniform_ZnS source mode.");
    return G4ThreeVector();
  }

  PrepareZnSSphereCdf(detector);
  if (fZnSSphereTotalWeight <= 0.0 || fZnSSphereCdf.size() != spheres.size())
  {
    G4Exception("StageDOpticalPrimaryGeneratorAction::SampleUniformPointInZnSSphere",
                "BNZS_D_PRI_003", FatalException,
                "Invalid ZnS sphere weights.");
    return G4ThreeVector();
  }

  const G4double target = G4UniformRand() * fZnSSphereTotalWeight;
  const auto selectedIt = std::lower_bound(
      fZnSSphereCdf.begin(), fZnSSphereCdf.end(), target);
  const std::size_t selectedIndex = selectedIt == fZnSSphereCdf.end()
                                        ? spheres.size() - 1
                                        : static_cast<std::size_t>(
                                              selectedIt - fZnSSphereCdf.begin());
  const auto *selected = &spheres[selectedIndex];

  const G4double radius = selected->radius * std::cbrt(G4UniformRand());
  return detector->WrapToPrimaryCell(selected->center + radius * RandomUnitVector());
}

void StageDOpticalPrimaryGeneratorAction::PrepareZnSSphereCdf(
    const DetectorConstruction *detector) const
{
  if (detector == nullptr)
    return;

  const auto &spheres = detector->GetZnSSpheres();
  const std::uint64_t placementSeed = detector->GetPlacementSeed();
  if (fZnSCdfDetector == detector &&
      fZnSCdfPlacementSeed == placementSeed &&
      fZnSCdfSphereCount == spheres.size() &&
      fZnSSphereCdf.size() == spheres.size())
    return;

  fZnSCdfDetector = detector;
  fZnSCdfPlacementSeed = placementSeed;
  fZnSCdfSphereCount = spheres.size();
  fZnSSphereCdf.clear();
  fZnSSphereCdf.reserve(spheres.size());
  fZnSSphereTotalWeight = 0.0;
  for (const auto &sphere : spheres)
  {
    fZnSSphereTotalWeight += std::pow(sphere.radius, 3);
    fZnSSphereCdf.push_back(fZnSSphereTotalWeight);
  }
}

G4ThreeVector StageDOpticalPrimaryGeneratorAction::SampleUniformPointInWholeRve(
    std::string &phaseName) const
{
  const auto *detector = ResolveDetector();
  if (detector == nullptr)
  {
    G4Exception("StageDOpticalPrimaryGeneratorAction::SampleUniformPointInWholeRve",
                "BNZS_D_PRI_006", FatalException,
                "DetectorConstruction is null.");
    phaseName = "Unknown";
    return G4ThreeVector();
  }

  const G4double halfX = detector->GetPatchHalfXUm() * um;
  const G4double halfY = detector->GetPatchHalfYUm() * um;
  const G4double halfZ = detector->GetPatchHalfZUm() * um;
  const G4ThreeVector point(
      (2.0 * G4UniformRand() - 1.0) * halfX,
      (2.0 * G4UniformRand() - 1.0) * halfY,
      (2.0 * G4UniformRand() - 1.0) * halfZ);
  phaseName = PhaseToString(detector->FindPhaseAtPoint(point));
  if (phaseName != "Matrix" && phaseName != "BN" && phaseName != "ZnS")
  {
    G4Exception("StageDOpticalPrimaryGeneratorAction::SampleUniformPointInWholeRve",
                "BNZS_D_PRI_007", FatalException,
                ("Uniform RVE source resolved to invalid phase: " + phaseName).c_str());
  }
  return point;
}

void StageDOpticalPrimaryGeneratorAction::GeneratePrimaries(G4Event *event)
{
  if (fConfig == nullptr)
  {
    G4Exception("StageDOpticalPrimaryGeneratorAction::GeneratePrimaries",
                "BNZS_D_PRI_004", FatalException,
                "AnalysisConfig pointer is null.");
    return;
  }

  if (fConfig->stageD_source_mode != "uniform_ZnS" &&
      fConfig->stageD_source_mode != "uniform_all_phase")
  {
    G4Exception("StageDOpticalPrimaryGeneratorAction::GeneratePrimaries",
                "BNZS_D_PRI_005", FatalException,
                ("Stage D source mode is not implemented in v1: " +
                 fConfig->stageD_source_mode)
                    .c_str());
    return;
  }

  std::string sourcePhase = "ZnS";
  G4ThreeVector sourcePosition;
  if (fConfig->stageD_source_mode == "uniform_ZnS")
  {
    sourcePosition = SampleUniformPointInZnSSphere();
    const auto *detector = ResolveDetector();
    if (detector == nullptr ||
        detector->FindPhaseAtPoint(sourcePosition) != DetectorConstruction::Phase::ZnS)
    {
      G4Exception("StageDOpticalPrimaryGeneratorAction::GeneratePrimaries",
                  "BNZS_D_PRI_008", FatalException,
                  "uniform_ZnS source sampling did not resolve to the ZnS phase.");
      return;
    }
    sourcePhase = "ZnS";
  }
  else
  {
    sourcePosition = SampleUniformPointInWholeRve(sourcePhase);
  }
  const G4ThreeVector direction = RandomUnitVector();
  const G4ThreeVector polarization = RandomPolarization(direction);
  const G4double energy = EnergyFromWavelengthNm(fConfig->stageD_wavelength_nm);

  fCurrentPhoton = StageDPhotonLaunchRecord{};
  fCurrentPhoton.geantEventID = event->GetEventID();
  fCurrentPhoton.wavelength_nm = fConfig->stageD_wavelength_nm;
  fCurrentPhoton.source_mode = fConfig->stageD_source_mode;
  fCurrentPhoton.source_phase = sourcePhase;
  fCurrentPhoton.source_position = sourcePosition;
  fCurrentPhoton.momentum_direction = direction;
  fCurrentPhoton.polarization = polarization;
  fCurrentPhoton.photon_weight = 1.0;
  fCurrentPhoton.is_continuation = false;

  fParticleGun->SetParticleEnergy(energy);
  fParticleGun->SetParticlePosition(sourcePosition);
  fParticleGun->SetParticleMomentumDirection(direction);
  fParticleGun->SetParticlePolarization(polarization);

  const G4int vertexCountBefore = event->GetNumberOfPrimaryVertex();
  fParticleGun->GeneratePrimaryVertex(event);
  if (event->GetNumberOfPrimaryVertex() > vertexCountBefore)
  {
    auto *vertex = event->GetPrimaryVertex(vertexCountBefore);
    if (vertex != nullptr && vertex->GetPrimary() != nullptr)
      vertex->GetPrimary()->SetWeight(1.0);
  }
}
