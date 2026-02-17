// ********************************************************************
// SteppingAction.cc
// ********************************************************************

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"

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
    auto volume   = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    auto particle = step->GetTrack()->GetDefinition();
    G4Track* track = step->GetTrack();
    G4String volName = volume->GetName();

    auto analysisManager = G4AnalysisManager::Instance();
    G4double particleEnergy = step->GetPreStepPoint()->GetKineticEnergy();

    // ========================================================================
    // CAPTURE FIRST PROTON ENTERING GAS VOLUME
    // 
    // ========================================================================
    if (!fEventAction->IsPositionSet()     &&
        particle->GetPDGCharge() != 0      &&
        volName == "DetectorBox")          
    {
         //where proton currently is in gas
        G4ThreeVector pos = step->GetPreStepPoint()->GetPosition();  
        G4ThreeVector dir = track->GetMomentumDirection();           

        fEventAction->SetInteractionPosition(pos);
        fEventAction->SetPrimaryDirection(dir);

        // Keep debug output but comment out for production runs
         G4cout << "First charged particle in gas: ("
               << pos.x()/mm << ", "
                << pos.y()/mm << ", "
                << pos.z()/mm << ") mm" << G4endl;
    }

    // ========================================================================
    // PROTON ENERGY TRACKING
    // ========================================================================
    if (particle->GetParticleName() == "proton") {

        // Proton entering converter 
        if (volName == "Conv"&& particle->GetParticleName() == "proton")
            analysisManager->FillH1(10, particleEnergy);  // proton_conv

        // Proton entering cathode mylar
        if (volName == "MylK"&& particle->GetParticleName() == "proton")
            analysisManager->FillH1(11, particleEnergy);  // proton_cathode

        // Proton inside gas 
        if (volName == "DetectorBox" && particle->GetParticleName() == "proton")
            analysisManager->FillH1(12, particleEnergy);  // proton_gas
    }

    // ========================================================================
    // OPTICAL PHOTON HANDLING
    // Records which SiPM was hit and photon energy
    // ========================================================================
    if (volName == "SiPM" && particle == G4OpticalPhoton::Definition()) {
        G4int copyNo = volume->GetCopyNo();

        analysisManager->FillH1(0, particleEnergy);  // Photon_Energy

        if      (copyNo <  25) analysisManager->FillH1(1, copyNo);       // SiPM_top
        else if (copyNo <  50) analysisManager->FillH1(2, copyNo - 25);  // SiPM_bottom
        else if (copyNo <  75) analysisManager->FillH1(3, copyNo - 50);  // SiPM_left
        else if (copyNo < 100) analysisManager->FillH1(4, copyNo - 75);  // SiPM_right

        step->GetTrack()->SetTrackStatus(fStopAndKill);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
