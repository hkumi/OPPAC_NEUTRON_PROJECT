#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

#include "G4SystemOfUnits.hh"  
#include "G4PhysicalConstants.hh"

using namespace B4;

namespace B4a
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(const DetectorConstruction* detConstruction,
                               EventAction* eventAction)
  : fDetConstruction(detConstruction),
    fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    if (!step) return;
    
    auto preStepPoint = step->GetPreStepPoint();
    if (!preStepPoint) return;
    
    auto volume = preStepPoint->GetTouchableHandle()->GetVolume();
    if (!volume) return;
    
    G4String volName = volume->GetName();
    auto particle = step->GetTrack()->GetDefinition();
    G4String particleName = particle->GetParticleName();
    G4double energy = preStepPoint->GetKineticEnergy();
    
    auto analysisManager = G4AnalysisManager::Instance();
    
    // Get current event ID
    G4int eventID = 0;
    G4Event* currentEvent = G4EventManager::GetEventManager()->GetNonconstCurrentEvent();
    if (currentEvent) {
        eventID = currentEvent->GetEventID();
    }

    // Debug: Track neutron interactions
    if (particleName == "neutron") {
        const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
        if (process) {
            G4String processName = process->GetProcessName();
            if (processName != "Transportation") {
                G4cout << "NEUTRON INTERACTION: " << processName 
                       << " in " << volName 
                       << " at E=" << energy/MeV << " MeV" << G4endl;
            }
        }
    }

    // Optical photon detection 
    if (volName == "SiPM" && particleName == "opticalphoton") {
        G4int copyNo = volume->GetCopyNo();
        G4double photonEnergy = preStepPoint->GetKineticEnergy();
        G4double time = preStepPoint->GetGlobalTime();
        
        // Get sensor position
        G4ThreeVector position = volume->GetTranslation();
        
        G4cout << "OPTICAL PHOTON DETECTED: Event " << eventID 
               << ", SiPM " << copyNo 
               << " E=" << photonEnergy/eV << " eV" 
               << " Pos(" << position.x()/mm << ", " << position.y()/mm << ") mm" << G4endl;

        // Fill detailed sensor data ntuple
        analysisManager->FillNtupleIColumn(1, 0, eventID);              // EventID
        analysisManager->FillNtupleIColumn(1, 1, copyNo);               // SensorID
        analysisManager->FillNtupleDColumn(1, 2, position.x()/mm);      // PosX (mm)
        analysisManager->FillNtupleDColumn(1, 3, position.y()/mm);      // PosY (mm)
        analysisManager->FillNtupleDColumn(1, 4, photonEnergy/eV);      // Energy (eV)
        analysisManager->FillNtupleDColumn(1, 5, time/ns);              // Time (ns)
        analysisManager->FillNtupleIColumn(1, 6, 1);                    // PhotonCount (1 per step)
        analysisManager->AddNtupleRow(1);
        
        // Also update histograms for quick look
        analysisManager->FillH1(0, photonEnergy); // Energy_SiPM histogram

        // Fill the appropriate SiPM array histogram (for 54 sensors per array)
        if (copyNo < 54) { 
            analysisManager->FillH1(1, copyNo); // SiPM_bottom
        }
        else if (copyNo < 108) { 
            analysisManager->FillH1(2, copyNo - 54); // SiPM_left  
        }
        else if (copyNo < 162) { 
            analysisManager->FillH1(3, copyNo - 108); // SiPM_up
        }
        else if (copyNo < 216) { 
            analysisManager->FillH1(4, copyNo - 162); // SiPM_right
        }

        // Kill the photon after detection
        step->GetTrack()->SetTrackStatus(fStopAndKill); 
    }

    // Proton tracking
    if (particleName == "proton") {
        const G4VProcess* preProcess = step->GetPreStepPoint()->GetProcessDefinedStep();
        G4String preProcessName = preProcess ? preProcess->GetProcessName() : "none";
        
        if (preProcessName == "none") {
            G4cout << "PROTON CREATED: Event " << eventID 
                   << " E=" << energy/MeV << " MeV in " << volName << G4endl;
            analysisManager->FillH1(5, energy); // proton_conv histogram
        }
        
        if (volName == "MylA") {
            G4cout << "Proton reached mylar sheet: E=" << energy/MeV << " MeV" << G4endl;
            analysisManager->FillH1(6, energy); // proton_myl histogram
        }
        
        if (volName == "World") {
            G4cout << "Proton entered gas volume: E=" << energy/MeV << " MeV" << G4endl;
            analysisManager->FillH1(7, energy); // proton_gas histogram
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
