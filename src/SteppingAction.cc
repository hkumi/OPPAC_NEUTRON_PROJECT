// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B4/B4a/src/SteppingAction.cc
/// \brief Implementation of the B4a::SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

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
// Collect energy and track length step by step

  // get volume of the current step
  auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();
  auto particle = step->GetTrack()->GetDefinition();
  G4String volName = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
  
  // process of post and pre step points
  const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
  G4String processName = process->GetProcessName();
  const G4VProcess* preProcess = step->GetPreStepPoint()->GetProcessDefinedStep();
  G4String preProcessName;
  if (preProcess) { preProcessName = preProcess->GetProcessName(); }
  else { preProcessName = "none"; }

  auto analysisManager = G4AnalysisManager::Instance();

  G4double photonEnergy = step->GetPreStepPoint()->GetKineticEnergy();
  G4int pdgCode = particle->GetPDGEncoding();
  G4double x = step->GetPostStepPoint()->GetPosition().x();
  G4double y = step->GetPostStepPoint()->GetPosition().y();
  G4double z = step->GetPostStepPoint()->GetPosition().z();

  G4double prex = step->GetPreStepPoint()->GetPosition().x();
  G4double prey = step->GetPreStepPoint()->GetPosition().y();
  G4double prez = step->GetPreStepPoint()->GetPosition().z();

  // Print relevant information
  //if (particle != G4OpticalPhoton::Definition()) {
  //    G4cout << particle->GetParticleName() << " event: " << preProcessName << " - " << processName << 
  //        "; Initial energy: " << photonEnergy << " MeV" << G4endl;
  //}

  if (particle->GetParticleName() == "proton" && preProcessName == "none") {
	  //G4cout << "Recoil proton produced with energy: " << photonEnergy * 1e3 << " keV " << G4endl;
	  analysisManager->FillH1(8, photonEnergy);  // Changed from 5 to 8 for proton_conv
  }

  if (volName == "MylA" && particle->GetParticleName() == "proton") {
	 // G4cout << "Proton reached mylar sheet with energy: " << photonEnergy * 1e3 << " keV " << G4endl;
	  analysisManager->FillH1(9, photonEnergy);  // Changed from 6 to 9 for proton_myl
  }

  if (volName == "World" && particle->GetParticleName() == "proton") {
      //G4cout << "Proton entered field region with energy: " << photonEnergy * 1e3 << " keV " << G4endl;
	  analysisManager->FillH1(10, photonEnergy); // Changed from 7 to 10 for proton_gas
  }

  if (volName == "SiPM" && pdgCode == -22) {
	 G4int copyNo = volume->GetCopyNo();
      //G4cout << "SiPM " << copyNo << " reached, PDG: " << pdgCode << G4endl;

      // Fill SiPM histograms (indices 1-4 as before)
      if (copyNo < 25) { analysisManager->FillH1(1, copyNo); }
      else if (copyNo < 50 && copyNo >= 25) { analysisManager->FillH1(2, copyNo - 25); }
      else if (copyNo < 75 && copyNo >= 50) { analysisManager->FillH1(3, copyNo - 50); }
	  else if (copyNo < 100 && copyNo >= 75) { analysisManager->FillH1(4, copyNo - 75); }

      step->GetTrack()->SetTrackStatus(fStopAndKill); 
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
