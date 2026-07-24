#include "StageAPrimaryGeneratorAction.hh"

#include "AnalysisConfig.hh"
#include "DetectorConstruction.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Event.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Exception.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"

StageAPrimaryGeneratorAction::StageAPrimaryGeneratorAction(AnalysisConfig *config)
    : G4VUserPrimaryGeneratorAction(),
      fConfig(config),
      fParticleGun(nullptr)
{
    if (fConfig == nullptr)
    {
        G4Exception("StageAPrimaryGeneratorAction::StageAPrimaryGeneratorAction",
                    "BNZS_A_001", FatalException,
                    "AnalysisConfig pointer is null.");
        return;
    }

    fParticleGun = new G4ParticleGun(1);

    auto *particleTable = G4ParticleTable::GetParticleTable();
    auto *neutron = particleTable->FindParticle("neutron");
    if (neutron == nullptr)
    {
        G4Exception("StageAPrimaryGeneratorAction::StageAPrimaryGeneratorAction",
                    "BNZS_A_002", FatalException,
                    "Cannot find neutron definition.");
        return;
    }

    fParticleGun->SetParticleDefinition(neutron);
    fParticleGun->SetParticleEnergy(0.0253 * eV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));

    G4cout << "[StageAPrimaryGeneratorAction] Initialized:"
           << " particle=neutron"
           << " energy=0.0253 eV"
           << " direction=(0,0,-1)"
           << G4endl;
}

StageAPrimaryGeneratorAction::~StageAPrimaryGeneratorAction()
{
    delete fParticleGun;
}

void StageAPrimaryGeneratorAction::GeneratePrimaries(G4Event *event)
{
    if (fConfig == nullptr || fParticleGun == nullptr)
    {
        G4Exception("StageAPrimaryGeneratorAction::GeneratePrimaries",
                    "BNZS_A_003", FatalException,
                    "Generator is not properly initialized.");
        return;
    }

    const auto *detector = dynamic_cast<const DetectorConstruction *>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    if (detector == nullptr)
        return;

    // 在 patch 前表面上方 1 um 处，从上往下入射
    const G4double sourceZ = detector->GetBoxHalfZ() + 1.0 * um;

    // --- 修改这里 ---
    // 假设你想把边长缩小为原来的 90%（即面积缩小到 81%）
    const G4double scaleFactor = 0.92;
    const G4double halfX = detector->GetBoxHalfX() * scaleFactor;
    const G4double halfY = detector->GetBoxHalfY() * scaleFactor;
    // ---------------

    const G4double x = (2.0 * G4UniformRand() - 1.0) * halfX;
    const G4double y = (2.0 * G4UniformRand() - 1.0) * halfY;

    fParticleGun->SetParticlePosition(G4ThreeVector(x, y, sourceZ));

    static G4int dbgPrimaryPrintCount = 0;
    if (dbgPrimaryPrintCount < 5)
    {
        ++dbgPrimaryPrintCount;

        G4cout << "[StageA-primary-debug]"
               << " eventID=" << event->GetEventID()
               << " x_um=" << x / um
               << " y_um=" << y / um
               << " z_um=" << sourceZ / um
               << G4endl;
    }

    fParticleGun->GeneratePrimaryVertex(event);
}
