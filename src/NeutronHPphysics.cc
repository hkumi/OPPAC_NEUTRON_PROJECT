#include "NeutronHPphysics.hh"

#include "G4GenericMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessTable.hh"

#include "G4HadronElasticProcess.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4ParticleHPThermalScatteringData.hh"
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPThermalScattering.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4ParticleHPInelasticData.hh"
#include "G4ParticleHPInelastic.hh"

#include "G4NeutronCaptureProcess.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4ParticleHPCapture.hh"

#include "G4NeutronFissionProcess.hh"
#include "G4ParticleHPFissionData.hh"
#include "G4ParticleHPFission.hh"

#include "G4SystemOfUnits.hh"

//...............................................

NeutronHPphysics::NeutronHPphysics(const G4String& name)
	: G4VPhysicsConstructor(name)
{
	//ConstructProcess();
}

NeutronHPphysics::~NeutronHPphysics()
{
}

//................................................

void NeutronHPphysics::ConstructProcess()
{
	G4ParticleDefinition* neutron = G4Neutron::Neutron();
	G4ProcessManager* pManager = neutron->GetProcessManager();

	// delete previous neutron processes
	G4VProcess* process = 0;
	process = pManager->GetProcess("hadElastic");
	if (process) pManager->RemoveProcess(process);

	process = pManager->GetProcess("neutronInelastic");
	if (process) pManager->RemoveProcess(process);

	process = pManager->GetProcess("nCapture");
	if (process) pManager->RemoveProcess(process);

	process = pManager->GetProcess("nFission");
	if (process) pManager->RemoveProcess(process);


	// re-create elastic process
	G4HadronElasticProcess* processE = new G4HadronElasticProcess();
	pManager->AddDiscreteProcess(processE);
	G4ParticleHPElastic* modelE_fast = new G4ParticleHPElastic();
	processE->RegisterMe(modelE_fast);
	processE->AddDataSet(new G4ParticleHPElasticData());
	modelE_fast->SetMinEnergy(4 * eV);
	G4ParticleHPThermalScattering* modelE_th = new G4ParticleHPThermalScattering();
	processE->RegisterMe(modelE_th);
	processE->AddDataSet(new G4ParticleHPThermalScatteringData());

	// re-create inelastic process
	G4HadronInelasticProcess* processI = new G4HadronInelasticProcess(
		"neutronInelastic", G4Neutron::Definition());
	pManager->AddDiscreteProcess(processI);
	G4ParticleHPInelastic* modelI = new G4ParticleHPInelastic();
	processI->RegisterMe(modelI);
	processI->AddDataSet(new G4ParticleHPInelasticData());

	// re-create neutron capture process
	G4NeutronCaptureProcess* processC = new G4NeutronCaptureProcess();
	pManager->AddDiscreteProcess(processC);
	G4ParticleHPCapture* modelC = new G4ParticleHPCapture();
	processC->RegisterMe(modelC);
	processC->AddDataSet(new G4ParticleHPCaptureData());

	// re-create neutron fission
	G4NeutronFissionProcess* processF = new G4NeutronFissionProcess();
	pManager->AddDiscreteProcess(processF);
	G4ParticleHPFission* modelF = new G4ParticleHPFission();
	processF->RegisterMe(modelF);
	processF->AddDataSet(new G4ParticleHPFissionData());

}