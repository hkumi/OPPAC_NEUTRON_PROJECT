
#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "NeutronHPphysics.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4EmStandardPhysics.hh"
#include "G4OpticalPhysics.hh"


//......................................


PhysicsList::PhysicsList()
	: G4VModularPhysicsList()
{
	SetVerboseLevel(1);

	// Registrar física electromagnética
	RegisterPhysics(new G4EmStandardPhysics());

	// Registrar física óptica
	auto opticalPhysics = new G4OpticalPhysics();
	RegisterPhysics(opticalPhysics);

	// Registrar tu física de neutrones HP
	RegisterPhysics(new NeutronHPphysics("neutronHP"));
}


PhysicsList::~PhysicsList()
{}

//......................................

void PhysicsList::ConstructParticle()
{
	G4BosonConstructor pBosonConstructor;
	pBosonConstructor.ConstructParticle();

	G4LeptonConstructor pLeptonConstructor;
	pLeptonConstructor.ConstructParticle();

	G4MesonConstructor pMesonConstructor;
	pMesonConstructor.ConstructParticle();

	G4BaryonConstructor pBaryonConstructor;
	pBaryonConstructor.ConstructParticle();

	G4IonConstructor pIonConstructor;
	pIonConstructor.ConstructParticle();

	G4ShortLivedConstructor pShortLivedConstructor;
	pShortLivedConstructor.ConstructParticle();
}

