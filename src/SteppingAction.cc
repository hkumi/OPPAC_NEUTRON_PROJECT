//
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
// * use.  Please see the license in the file  LICENSE  and URL above *
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
/// 
/*
#include "SteppingAction.hh"
#include "G4Step.hh"
#include "G4LogicalVolume.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4StepPoint.hh"

#include "G4ios.hh" // Para usar G4cout
using namespace B4;
namespace B4a {
    SteppingAction::SteppingAction(const DetectorConstruction* detConstruction,
        EventAction* eventAction)
        : fDetConstruction(detConstruction),
        fEventAction(eventAction)
    {
    }

    void SteppingAction::UserSteppingAction(const G4Step* step)
    {
        // Obtener el punto donde termina el paso
        G4StepPoint* postPoint = step->GetPostStepPoint();
        if (!postPoint) return;

        // Obtener volumen físico donde ocurre el paso
        G4VPhysicalVolume* volume = postPoint->GetTouchableHandle()->GetVolume();
        if (!volume) return;

        G4String volumeName = volume->GetName();

        // Imprimir nombre del volumen donde ocurre el paso
        G4cout << "Paso en volumen: " << volumeName << G4endl;
    }

}
*/
#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"

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

  //optical photon detection 
  if (volName == "SiPM" && particleName == "opticalphoton") {
    G4int copyNo = volume->GetCopyNo();
    G4cout << "OPTICAL PHOTON DETECTED: SiPM " << copyNo 
           << " E=" << energy/eV << " eV" << G4endl;

    analysisManager->FillH1(0, energy); // Energy_SiPM histogram

    // Fill the appropriate SiPM array histogram
    if (copyNo < 25) { 
        analysisManager->FillH1(1, copyNo); // SiPM_bottom
    }
    else if (copyNo < 50) { 
        analysisManager->FillH1(2, copyNo - 25); // SiPM_left  
    }
    else if (copyNo < 75) { 
        analysisManager->FillH1(3, copyNo - 50); // SiPM_up
    }
    else if (copyNo < 100) { 
        analysisManager->FillH1(4, copyNo - 75); // SiPM_right
    }

    // Kill the photon after detection
    step->GetTrack()->SetTrackStatus(fStopAndKill); 
  }

  //  proton tracking ...
  if (particleName == "proton") {
    const G4VProcess* preProcess = step->GetPreStepPoint()->GetProcessDefinedStep();
    G4String preProcessName = preProcess ? preProcess->GetProcessName() : "none";
    
    if (preProcessName == "none") {
      G4cout << "PROTON CREATED: E=" << energy/MeV << " MeV in " << volName << G4endl;
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

}
