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

        // ====================================================================
        // 1D HISTOGRAMS
        // ====================================================================
        
        analysisManager->CreateH1("Photon_Energy", "Photon energy spectrum in SiPMs", 
                                  100, 1.5 * eV, 6 * eV,"eV", "none");               // Histogram ID 0
        
        // Note: Histograms 1-4 are used in SteppingAction.cc
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
                                  1000, 0, 1 * MeV, "keV", "none");             // Histogram ID 10
        analysisManager->CreateH1("proton_myl", "Protons reaching mylar sheet", 
                                  1000, 0, 1 * MeV, "keV", "none");             // Histogram ID 11
        analysisManager->CreateH1("proton_gas", "Protons reaching gas volume", 
                                  1000, 0, 1 * MeV, "keV", "none");             // Histogram ID 12
        
        // Resolution magnitude histogram (used in ReconstructFullEvent):
        analysisManager->CreateH1("Resolution_Mag", "Position resolution magnitude", 
                                  100, 0, 20, "mm");                    // Histogram ID 13
        
        // Border event histograms (used in ReconstructBorderEvent):
        analysisManager->CreateH1("Border_Res_X", "Border event X resolution", 
                                  100, -10, 10, "mm");                  // Histogram ID 14
        analysisManager->CreateH1("Border_Res_Y", "Border event Y resolution", 
                                  100, -10, 10, "mm");                  // Histogram ID 15

        // ====================================================================
        // 2D HISTOGRAMS - ALL CREATED BEFORE NTUPLES
        // ====================================================================
        
        analysisManager->CreateH2("Position_2D_Full", "Full event reconstructed position", 
                                  100, -50, 50, 100, -50, 50, 
                                  "mm", "mm");                           // Histogram ID 0 (2D)
        
        analysisManager->CreateH2("Position_2D_Border", "Border event reconstructed position", 
                                  100, -50, 50, 100, -50, 50, 
                                  "mm", "mm");                           // Histogram ID 1 (2D)

        // H2(2): Theta (polar angle) vs ΔX
        analysisManager->CreateH2("Theta_vs_ResX", 
            "Residual X vs Polar Angle;#theta (degrees);#DeltaX (mm)",
            90, 0, 90,      // theta: 0-90 degrees
            100, -2, 2);    // ΔX: ±2 mm                                // Histogram ID 2 (2D)
        
        // H2(3): Theta (polar angle) vs ΔY  
        analysisManager->CreateH2("Theta_vs_ResY",
            "Residual Y vs Polar Angle;#theta (degrees);#DeltaY (mm)",
            90, 0, 90,      // theta: 0-90 degrees
            100, -2, 2);    // ΔY: ±2 mm                                // Histogram ID 3 (2D)
        
        // H2(4): Theta vs Combined Resolution
        analysisManager->CreateH2("Theta_vs_Resolution",
            "Combined Resolution vs Polar Angle;#theta (degrees);Resolution (mm)",
            90, 0, 90,      // theta: 0-90 degrees
            100, 0, 2);     // resolution: 0-2 mm                       // Histogram ID 4 (2D)
        
        // H2(5): Phi (azimuthal angle) vs ΔX
        analysisManager->CreateH2("Phi_vs_ResX",
            "Residual X vs Azimuthal Angle;#phi (degrees);#DeltaX (mm)",
            180, -180, 180, // phi: -180 to +180 degrees
            100, -2, 2);    // ΔX: ±2 mm                                // Histogram ID 5 (2D)
        
        // H2(6): Phi (azimuthal angle) vs ΔY
        analysisManager->CreateH2("Phi_vs_ResY",
            "Residual Y vs Azimuthal Angle;#phi (degrees);#DeltaY (mm)",
            180, -180, 180, // phi: -180 to +180 degrees
            100, -2, 2);    // ΔY: ±2 mm                                // Histogram ID 6 (2D)

        // ====================================================================
        // NTUPLES
        // ====================================================================
        
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

        // Reconstruction ntuple that includes both full and border events:
        analysisManager->CreateNtuple("Reconstruction", "Event reconstruction data");  // Ntuple ID 4
        analysisManager->CreateNtupleDColumn("recX");
        analysisManager->CreateNtupleDColumn("recY");
        analysisManager->CreateNtupleDColumn("trueX");
        analysisManager->CreateNtupleDColumn("trueY");
        analysisManager->CreateNtupleIColumn("borderFlag");  // 0=full, 1=border
        analysisManager->CreateNtupleIColumn("nArrays");     // Number of arrays with signal
        analysisManager->FinishNtuple();
        
        // Ntuple 5: Complete reconstruction data with angles
        analysisManager->CreateNtuple("AngleAnalysis", "Position and Angle Data");  // Ntuple ID 5
        analysisManager->CreateNtupleDColumn("x_rec");       // 0: Reconstructed X (mm)
        analysisManager->CreateNtupleDColumn("y_rec");       // 1: Reconstructed Y (mm)
        analysisManager->CreateNtupleDColumn("x_true");      // 2: True X (mm)
        analysisManager->CreateNtupleDColumn("y_true");      // 3: True Y (mm)
        analysisManager->CreateNtupleDColumn("dx");          // 4: ΔX = x_rec - x_true (mm)
        analysisManager->CreateNtupleDColumn("dy");          // 5: ΔY = y_rec - y_true (mm)
        analysisManager->CreateNtupleDColumn("theta");       // 6: Polar angle (degrees)
        analysisManager->CreateNtupleDColumn("phi");         // 7: Azimuthal angle (degrees)
        analysisManager->CreateNtupleDColumn("resolution");  // 8: Combined resolution (mm)
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
