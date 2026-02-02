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
        auto analysisManager = G4AnalysisManager::Instance();

        analysisManager->SetVerboseLevel(1);
        analysisManager->SetNtupleMerging(true);

        // Create ONLY the histograms and ntuples you're actually using in detector.cc:

        // 1D Histograms:
        analysisManager->CreateH1("Photon_Energy", "Photon energy spectrum in SiPMs", 
                                  100, 1.5 * eV, 6 * eV);               // Histogram ID 0
        
        // Note: Histograms 1-4 are used in SteppingAction.cc, not detector.cc
        // but you're filling them in SteppingAction so we keep them
        analysisManager->CreateH1("SiPM_top", "Number of photons detected in top array", 
                                  25, 0, 25);                           // Histogram ID 1
        analysisManager->CreateH1("SiPM_bottom", "Number of photons detected in bottom array", 
                                  25, 0, 25);                           // Histogram ID 2
        analysisManager->CreateH1("SiPM_left", "Number of photons detected in left array", 
                                  25, 0, 25);                           // Histogram ID 3
        analysisManager->CreateH1("SiPM_right", "Number of photons detected in right array", 
                                  25, 0, 25);                           // Histogram ID 4
        
        analysisManager->CreateH1("Total_Hits", "Total hits per event", 
                                  100, 0, 100);                         // Histogram ID 5
        
        analysisManager->CreateH1("Weighted_X", "Weighted X position", 
                                  100, -50, 50, "mm");                  // Histogram ID 6
        analysisManager->CreateH1("Weighted_Y", "Weighted Y position", 
                                  100, -50, 50, "mm");                  // Histogram ID 7
        
        // Resolution histograms (used in ReconstructFullEvent):
        analysisManager->CreateH1("Resolution_X", "X resolution", 
                                  100, -10, 10, "mm");                  // Histogram ID 8
        analysisManager->CreateH1("Resolution_Y", "Y resolution", 
                                  100, -10, 10, "mm");                  // Histogram ID 9
        
        // Proton energy histograms (used in SteppingAction):
        analysisManager->CreateH1("proton_conv", "Recoil proton energy", 
                                  1000, 0, 2 * MeV, "keV");             // Histogram ID 10
        analysisManager->CreateH1("proton_myl", "Protons reaching mylar sheet", 
                                  1000, 0, 2 * MeV, "keV");             // Histogram ID 11
        analysisManager->CreateH1("proton_gas", "Protons reaching gas volume", 
                                  1000, 0, 2 * MeV, "keV");             // Histogram ID 12
        
        // Resolution magnitude histogram (used in ReconstructFullEvent):
        analysisManager->CreateH1("Resolution_Mag", "Position resolution magnitude", 
                                  100, 0, 20, "mm");                    // Histogram ID 13
        
        // Border event histograms (used in ReconstructBorderEvent):
        analysisManager->CreateH1("Border_Res_X", "Border event X resolution", 
                                  100, -10, 10, "mm");                  // Histogram ID 14
        analysisManager->CreateH1("Border_Res_Y", "Border event Y resolution", 
                                  100, -10, 10, "mm");                  // Histogram ID 15

        // 2D Histograms:
        analysisManager->CreateH2("Position_2D_Full", "Full event reconstructed position", 
                                  100, -50, 50, 100, -50, 50, 
                                  "mm", "mm");                           // Histogram ID 0 (2D)
        analysisManager->CreateH2("Position_2D_Border", "Border event reconstructed position", 
                                  100, -50, 50, 100, -50, 50, 
                                  "mm", "mm");                           // Histogram ID 1 (2D)

        // Ntuples (4 arrays for sensor data - used in RecordSensorData):
        for (int i = 0; i < 4; i++) {
            std::string ntupleName = "Array_" + std::to_string(i);
            std::string ntupleTitle = "Array " + std::to_string(i) + " data";
            
            analysisManager->CreateNtuple(ntupleName, ntupleTitle);     // Ntuple IDs 0-3
            analysisManager->CreateNtupleDColumn("posX");
            analysisManager->CreateNtupleDColumn("posY");
            analysisManager->CreateNtupleIColumn("event");
            analysisManager->CreateNtupleIColumn("copyNo");
            analysisManager->FinishNtuple();
        }

        // REMOVED: Weighted_X and Weighted_Y ntuples (not being filled in detector.cc)
        // REMOVED: Original B4 ntuple (not being used)

        // Add a single reconstruction ntuple that includes both full and border events:
        analysisManager->CreateNtuple("Reconstruction", "Event reconstruction data");  // Ntuple ID 4
        analysisManager->CreateNtupleDColumn("recX");
        analysisManager->CreateNtupleDColumn("recY");
        analysisManager->CreateNtupleDColumn("trueX");
        analysisManager->CreateNtupleDColumn("trueY");
        analysisManager->CreateNtupleIColumn("borderFlag");  // 0=full, 1=border
        analysisManager->CreateNtupleIColumn("nArrays");     // Number of arrays with signal
        analysisManager->FinishNtuple();
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void RunAction::BeginOfRunAction(const G4Run* /*run*/)
    {
        // Get analysis manager
        auto analysisManager = G4AnalysisManager::Instance();

        // Open an output file
        G4String fileName = "B4.root";
        analysisManager->OpenFile(fileName);
        G4cout << "Using " << analysisManager->GetType() << G4endl;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void RunAction::EndOfRunAction(const G4Run* /*run*/)
    {
        auto analysisManager = G4AnalysisManager::Instance();

        // save histograms & ntuple
        analysisManager->Write();
        analysisManager->CloseFile();
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
