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

void SteppingAction::SaveEventData(G4int eventID)
{
    // Save to a text file for macro processing
    std::ofstream outfile;
    outfile.open("position_data.txt", std::ios_base::app);
    
    // Format: EventID nBottom nTop nLeft nRight bottom_positions... top_positions... left_positions... right_positions...
    outfile << eventID << " ";
    outfile << fXBottom.size() << " ";
    outfile << fXTop.size() << " ";
    outfile << fYLeft.size() << " ";
    outfile << fYRight.size() << " ";
    
    // Write all positions
    for (double x : fXBottom) outfile << x << " ";
    for (double x : fXTop) outfile << x << " ";
    for (double y : fYLeft) outfile << y << " ";
    for (double y : fYRight) outfile << y << " ";
    
    outfile << "\n";
    outfile.close();
}

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
    //G4ThreeVector posPhotons = postStepPoint->GetPosition();//accessing the position
    
    auto analysisManager = G4AnalysisManager::Instance();
    
    // Get current event ID
    G4int eventID = 0;
    G4Event* currentEvent = G4EventManager::GetEventManager()->GetNonconstCurrentEvent();
    if (currentEvent) {
        eventID = currentEvent->GetEventID();
    }

    // Optical photon detection 
    if (volName == "SiPM" && particleName == "opticalphoton") {
        G4int copyNo = volume->GetCopyNo();
        G4double photonEnergy = preStepPoint->GetKineticEnergy();
        G4double time = preStepPoint->GetGlobalTime();
        
        // Get sensor position
        G4ThreeVector position = volume->GetTranslation();
        G4double posX = position.x()/mm;
        G4double posY = position.y()/mm;
        
        G4double x = step->GetPostStepPoint()->GetPosition().x()/mm;
        G4double y = step->GetPostStepPoint()->GetPosition().y()/mm;
        G4double z = step->GetPostStepPoint()->GetPosition().z()/mm;

        // ===== DEBUG OUTPUT =====
    G4cout << "\n=== PHOTON HIT DEBUG ===" << G4endl;
    G4cout << "Event ID: " << eventID << G4endl;
    G4cout << "Sensor copyNo: " << copyNo << G4endl;
    G4cout << "Sensor (volume) position: (" << posX << ", " << posY << ", 0) mm" << G4endl;
    G4cout << "Hit (post-step) position: (" << x << ", " << y << ", " << z << ") mm" << G4endl;
    G4cout << "Distance from sensor center: " 
           << std::sqrt((x-posX)*(x-posX) + (y-posY)*(y-posY)) << " mm" << G4endl;
        // STORE POSITIONS FOR MACRO PROCESSING (NEW CODE)
       // Check if this is a new event
       if (eventID != fCurrentEventID) {
        // Save previous event data if we have any
           if (!fXBottom.empty() || !fXTop.empty() || !fYLeft.empty() || !fYRight.empty()) {
               SaveEventData(fCurrentEventID);
           }
        // Reset for new event
        fCurrentEventID = eventID;
        fXBottom.clear();
        fXTop.clear();
        fYLeft.clear();
        fYRight.clear();
       }
        // Fill detailed sensor data ntuple
        analysisManager->FillNtupleIColumn(0, 0, eventID);              // EventID
        analysisManager->FillNtupleIColumn(0, 1, copyNo);               // SensorID
        analysisManager->FillNtupleDColumn(0, 2, x);      // PosX (mm)
        analysisManager->FillNtupleDColumn(0, 3, y);      // PosY (mm)
        analysisManager->FillNtupleDColumn(0, 4, photonEnergy/eV);      // Energy (eV)
        analysisManager->FillNtupleDColumn(0, 5, time/ns);              // Time (ns)
        analysisManager->FillNtupleIColumn(0, 6, 1);                    // PhotonCount (1 per step)
        analysisManager->AddNtupleRow(0);
        
        // Also update histograms for quick look
        analysisManager->FillH1(0, photonEnergy); // Energy_SiPM histogram

        // Fill the appropriate SiPM array histogram (for 25 sensors per array)
        if (copyNo < 25) { 
            fXBottom.push_back(x);
            analysisManager->FillH1(1, copyNo); // SiPM_bottom
            analysisManager->FillH1(8, x);      // Xbottom_pos (X position) -



        }
        else if (copyNo < 50 && copyNo >= 25) {
             fYLeft.push_back(y); 
            analysisManager->FillH1(2, copyNo - 25); // SiPM_left
             analysisManager->FillH1(9, y);           // Yleft_pos (Y position)
  
        }
        else if (copyNo < 75 && copyNo >= 50) {
             fXTop.push_back(x);  
            analysisManager->FillH1(3, copyNo - 50); // SiPM_up
             analysisManager->FillH1(10, x);           // Xtop_pos (X position)

        }
        else if (copyNo < 100 && copyNo >= 75) {
             fYRight.push_back(y); 
            analysisManager->FillH1(4, copyNo - 75); // SiPM_right
             analysisManager->FillH1(11, y);           // Yright_pos (Y position)

        }

        // Kill the photon after detection
        step->GetTrack()->SetTrackStatus(fStopAndKill); 
    }

    // Proton tracking
    if (particleName == "proton") {
        const G4VProcess* preProcess = step->GetPreStepPoint()->GetProcessDefinedStep();
        G4String preProcessName = preProcess ? preProcess->GetProcessName() : "none";
        
        if (preProcessName == "none") {
            //G4cout << "PROTON CREATED: Event " << eventID 
                   //<< " E=" << energy/MeV << " MeV in " << volName << G4endl;
            analysisManager->FillH1(5, energy); // proton_conv histogram
        }
        
        if (volName == "MylA") {
            //G4cout << "Proton reached mylar sheet: E=" << energy/MeV << " MeV" << G4endl;
            analysisManager->FillH1(6, energy); // proton_myl histogram
        }
        
        if (volName == "World") {
            //G4cout << "Proton entered gas volume: E=" << energy/MeV << " MeV" << G4endl;
            analysisManager->FillH1(7, energy); // proton_gas histogram
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
