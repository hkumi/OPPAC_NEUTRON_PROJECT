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
/// \file B4/B4a/src/RunAction.cc
/// \brief Implementation of the B4::RunAction class

#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    RunAction::RunAction()
    {
      // set printing event number per each event
      G4RunManager::GetRunManager()->SetPrintProgress(1);

      // Create analysis manager
      // The choice of the output format is done via the specified
      // file extension.
      auto analysisManager = G4AnalysisManager::Instance();

      // Create directories
      //analysisManager->SetHistoDirectoryName("histograms");
      //analysisManager->SetNtupleDirectoryName("ntuple");
      analysisManager->SetVerboseLevel(1);
      analysisManager->SetNtupleMerging(true);

      // Book histograms, ntuple
      //

      // Creating histograms
      
      analysisManager->CreateH1("Energy_SiPM", "Photon spectrum in SiPMs", 100, 1.5 * eV, 6 * eV, "eV", "none");
	  analysisManager->CreateH1("SiPM_bottom", "Number of photons detected in top array", 25, 0, 25);
      analysisManager->CreateH1("SiPM_left", "Number of photons detected in left array", 25, 0, 25);
      analysisManager->CreateH1("SiPM_up", "Number of photons detected in bottom array", 25, 0, 25);
      analysisManager->CreateH1("SiPM_right", "Number of photons detected in right array", 25, 0, 25);

	  //analysisManager->CreateH1("proton_energy", "Recoil proton energy", 100, pow(10, -8) * MeV, 1 * MeV, "MeV", "none", "log");
	  analysisManager->CreateH1("proton_conv", "Recoil proton energy", 1000, 0, 2 * MeV, "keV", "none");
	  analysisManager->CreateH1("proton_myl", "Protons reaching mylar sheet", 1000, 0, 2 * MeV, "keV", "none");
	  analysisManager->CreateH1("proton_gas", "Protons reaching gas volume", 1000, 0, 2 * MeV, "keV", "none");

      analysisManager->CreateNtuple("B4", "Energy");
      analysisManager->CreateNtupleDColumn("WSF");
      analysisManager->CreateNtupleDColumn("SiPM");
      analysisManager->CreateNtupleDColumn("InitialEnergy");
      
      analysisManager->FinishNtuple();
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void RunAction::BeginOfRunAction(const G4Run* /*run*/)
    {
      //inform the runManager to save random number seed
      //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

      // Get analysis manager
      auto analysisManager = G4AnalysisManager::Instance();

      // Open an output file
      //
      G4String fileName = "B4.root";
      // Other supported output types:
      // G4String fileName = "B4.csv";
      // G4String fileName = "B4.hdf5";
      // G4String fileName = "B4.xml";
      analysisManager->OpenFile(fileName);
      G4cout << "Using " << analysisManager->GetType() << G4endl;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void RunAction::EndOfRunAction(const G4Run* /*run*/)
    {
      // print histogram statistics
      //
    
      auto analysisManager = G4AnalysisManager::Instance();


      // save histograms & ntuple
      //
      analysisManager->Write();
      analysisManager->CloseFile();
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
