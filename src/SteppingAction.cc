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
  // Get volume of the current step
  auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  auto particle = step->GetTrack()->GetDefinition();
  G4String volName = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
  
  // Get track
  G4Track* track = step->GetTrack();
  
  // ========================================================================
  //  CAPTURE FIRST CHARGED PARTICLE INFORMATION
  // ========================================================================
  // This captures the interaction position and direction for angle analysis
  
  if (!fEventAction->IsPositionSet() && 
      particle->GetPDGCharge() != 0)  // Any charged particle
  {
      // Get vertex position (where particle was created)
      G4ThreeVector pos = track->GetVertexPosition();
      
      // Get vertex momentum direction (initial direction)
      G4ThreeVector dir = track->GetVertexMomentumDirection();
      
      // Store in EventAction
      fEventAction->SetInteractionPosition(pos);
      fEventAction->SetPrimaryDirection(dir);
      
      // Debug output
      G4cout << "\n=== CAPTURED IN STEPPING ACTION ===" << G4endl;
      G4cout << "Particle: " << particle->GetParticleName() << G4endl;
      G4cout << "Volume: " << volName << G4endl;
      G4cout << "Parent ID: " << track->GetParentID() << G4endl;
      G4cout << "Position: (" << pos.x()/mm << ", " << pos.y()/mm << ", " << pos.z()/mm << ") mm" << G4endl;
      G4cout << "Direction: (" << dir.x() << ", " << dir.y() << ", " << dir.z() << ")" << G4endl;
      
      // Calculate angles
      double theta = std::acos(std::abs(dir.z())) * 180.0 / CLHEP::pi;
      double phi = std::atan2(dir.y(), dir.x()) * 180.0 / CLHEP::pi;
      G4cout << "Theta: " << theta << " degrees" << G4endl;
      G4cout << "Phi: " << phi << " degrees" << G4endl;
      G4cout << "====================================\n" << G4endl;
  }
  
  // ========================================================================
 //========================================================================
  
  // Process of post and pre step points
  const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
  G4String processName = process->GetProcessName();
  const G4VProcess* preProcess = step->GetPreStepPoint()->GetProcessDefinedStep();
  G4String preProcessName;
  if (preProcess) { preProcessName = preProcess->GetProcessName(); }
  else { preProcessName = "none"; }

  auto analysisManager = G4AnalysisManager::Instance();

  G4double particleEnergy = step->GetPreStepPoint()->GetKineticEnergy();
  G4int pdgCode = particle->GetPDGEncoding();

  // Proton energy tracking 
  if (particle->GetParticleName() == "proton" && preProcessName == "none") {
      analysisManager->FillH1(10, particleEnergy);  // proton_conv
  }

  if (volName == "MylA" && particle->GetParticleName() == "proton") {
      analysisManager->FillH1(11, particleEnergy);  // proton_myl
  }

  if (volName == "World" && particle->GetParticleName() == "proton") {
      analysisManager->FillH1(12, particleEnergy);  // proton_gas
  }

  // OPTICAL PHOTON HANDLING 
  if (volName == "SiPM" && particle == G4OpticalPhoton::Definition()) {
      G4int copyNo = volume->GetCopyNo();
      
      // FILL PHOTON ENERGY HISTOGRAM 
      analysisManager->FillH1(0, particleEnergy);  // Photon_Energy histogram
      
      // Fill SiPM histograms (count photons per sensor)
      if (copyNo < 25) { 
          analysisManager->FillH1(1, copyNo);  // SiPM_top
      }
      else if (copyNo < 50) { 
          analysisManager->FillH1(2, copyNo - 25);  // SiPM_bottom
      }
      else if (copyNo < 75) { 
          analysisManager->FillH1(3, copyNo - 50);  // SiPM_left
      }
      else if (copyNo < 100) { 
          analysisManager->FillH1(4, copyNo - 75);  // SiPM_right
      }
      
      step->GetTrack()->SetTrackStatus(fStopAndKill);  
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
