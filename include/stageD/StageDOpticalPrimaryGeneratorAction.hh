#ifndef StageDOpticalPrimaryGeneratorAction_h
#define StageDOpticalPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "StageDOpticalStats.hh"

#include <cstdint>
#include <vector>

class G4Event;
class AnalysisConfig;
class DetectorConstruction;

class StageDOpticalPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  explicit StageDOpticalPrimaryGeneratorAction(AnalysisConfig *config);
  ~StageDOpticalPrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event *event) override;

  const StageDPhotonLaunchRecord &GetCurrentPhotonRecord() const { return fCurrentPhoton; }

private:
  const DetectorConstruction *ResolveDetector() const;
  void PrepareZnSSphereCdf(const DetectorConstruction *detector) const;
  G4ThreeVector SampleUniformPointInZnSSphere() const;
  G4ThreeVector SampleUniformPointInWholeRve(std::string &phaseName) const;
  G4ThreeVector RandomUnitVector() const;
  G4ThreeVector RandomPolarization(const G4ThreeVector &direction) const;

private:
  AnalysisConfig *fConfig;
  G4ParticleGun *fParticleGun;
  StageDPhotonLaunchRecord fCurrentPhoton;
  mutable const DetectorConstruction *fZnSCdfDetector;
  mutable std::uint64_t fZnSCdfPlacementSeed;
  mutable std::size_t fZnSCdfSphereCount;
  mutable std::vector<G4double> fZnSSphereCdf;
  mutable G4double fZnSSphereTotalWeight;
};

#endif
