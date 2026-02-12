// ============================================================================
// PrimaryGeneratorAction.cc
//
// Shoots a 2.5 MeV neutron beam along +Z.
//
// BEAM PROFILE – controlled by the flag below:
//
//   USE_POINT_SOURCE = true   → pencil beam at (0, 0)      ← current setting
//   USE_POINT_SOURCE = false  → Gaussian profile σ = 5 mm  ← switch when ready
//
// Change the one flag to switch; nothing else needs editing.
// ============================================================================

#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4AnalysisManager.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

namespace B4
{

// ── Beam-profile switch ──────────────────────────────────────────────────────
// true  = pencil beam at (0,0) — use this to compare algorithms
// false = Gaussian beam σ = 5 mm — use for imaging
static constexpr bool USE_POINT_SOURCE = true;
static constexpr double BEAM_SIGMA     = 5.0;  // mm (only used when false)
// ────────────────────────────────────────────────────────────────────────────

// ── Constructor ──────────────────────────────────────────────────────────────
PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    fParticleGun = new G4ParticleGun(1);

    auto* particleDef =
        G4ParticleTable::GetParticleTable()->FindParticle("neutron");
    fParticleGun->SetParticleDefinition(particleDef);

    //  direction: straight along +Z
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));

    G4cout << "PrimaryGeneratorAction: "
           << (USE_POINT_SOURCE ? "POINT SOURCE at (0,0)"
                                : "GAUSSIAN beam  sigma=" +
                                  std::to_string(BEAM_SIGMA) + " mm")
           << G4endl;
}

// ── Destructor ───────────────────────────────────────────────────────────────
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

// ── GeneratePrimaries ────────────────────────────────────────────────────────
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    // Energy: 2.5 MeV neutrons
    fParticleGun->SetParticleEnergy(2.5 * MeV);

    // Transverse position
    G4double xpos = 0.0;
    G4double ypos = 0.0;

    if (!USE_POINT_SOURCE) {
        // Gaussian beam 
        xpos = CLHEP::RandGauss::shoot(0.0, BEAM_SIGMA * mm);
        ypos = CLHEP::RandGauss::shoot(0.0, BEAM_SIGMA * mm);
    }
    // else: both remain 0.0  (pencil beam)

    // Starting position: -0.5 cm upstream of the valve (valve at z = -15 cm)
    G4double zpos = -0.5 * cm;

    fParticleGun->SetParticlePosition(G4ThreeVector(xpos, ypos, zpos));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, 1.0));

    fParticleGun->GeneratePrimaryVertex(anEvent);
}

} // namespace B4
