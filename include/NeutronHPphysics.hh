
#ifndef NeutronHPphysics_h
#define NeutronHPphysics_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"

//..............................................

class NeutronHPphysics : public G4VPhysicsConstructor
{
public:
	NeutronHPphysics(const G4String& name = "neutron");
	~NeutronHPphysics();

public:
	void ConstructProcess() override;
	void ConstructParticle() override {};

};

#endif