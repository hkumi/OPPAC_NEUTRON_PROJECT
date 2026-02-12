#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

namespace B4
{

RunAction::RunAction()
{
    G4RunManager::GetRunManager()->SetPrintProgress(1);

    auto* man = G4AnalysisManager::Instance();
    man->SetVerboseLevel(1);
    man->SetNtupleMerging(true);

    // 1D HISTOGRAMS
    // IDs 0-15 : algorithm 1
    // IDs 16-20: PAPER algorithm

    // Photon / sensor diagnostics
    man->CreateH1("Photon_Energy",
                  "Photon energy spectrum in SiPMs",
                  100, 1.5*eV, 6*eV, "eV", "none");              // ID 0

    man->CreateH1("SiPM_top",
                  "Photons detected – top array",
                  25, 0, 25);                                      // ID 1
    man->CreateH1("SiPM_bottom",
                  "Photons detected – bottom array",
                  25, 0, 25);                                      // ID 2
    man->CreateH1("SiPM_left",
                  "Photons detected – left array",
                  25, 0, 25);                                      // ID 3
    man->CreateH1("SiPM_right",
                  "Photons detected – right array",
                  25, 0, 25);                                      // ID 4

    man->CreateH1("Total_Hits",
                  "Total photon hits per event",
                  100, 0, 100);                                    // ID 5

    // Algorithm 1 (weight = N/σ²)
    man->CreateH1("Weighted_X",
                  "Algorithm 1: Reconstructed X position",
                  100, -14, 14, "mm");                             // ID 6
    man->CreateH1("Weighted_Y",
                  "Algorithm 1: Reconstructed Y position",
                  100, -14, 14, "mm");                             // ID 7
    man->CreateH1("Resolution_X",
                  "Algorithm 1: X residual (x_rec - x_true)",
                  200, -3, 3, "mm");                             // ID 8
    man->CreateH1("Resolution_Y",
                  "Algorithm 1: Y residual (y_rec - y_true)",
                  200, -3, 3, "mm");                             // ID 9

    // Proton energy histograms (SteppingAction)
    man->CreateH1("proton_conv",
                  "Recoil proton energy at conversion",
                  1000, 0, 1*MeV, "keV", "none");                  // ID 10
    man->CreateH1("proton_myl",
                  "Proton energy at mylar",
                  1000, 0, 1*MeV, "keV", "none");                  // ID 11
    man->CreateH1("proton_gas",
                  "Proton energy in gas volume",
                  1000, 0, 1*MeV, "keV", "none");                  // ID 12

    man->CreateH1("Resolution_Mag",
                  "Algorithm 1: Combined position resolution sqrt(dx²+dy²)",
                  200, 0, 20, "mm");                               // ID 13

    // IDs 14-15 were border-event histograms; kept as empty placeholders
    man->CreateH1("Border_Res_X_UNUSED",
                  "(unused – border reconstruction removed)",
                  1, -1, 1);                                       // ID 14
    man->CreateH1("Border_Res_Y_UNUSED",
                  "(unused – border reconstruction removed)",
                  1, -1, 1);                                       // ID 15

    // PAPER algorithm (weight = N/σ)
    man->CreateH1("Weighted_X_Paper",
                  "PAPER: Reconstructed X position",
                  100, -14, 14, "mm");                             // ID 16
    man->CreateH1("Weighted_Y_Paper",
                  "PAPER: Reconstructed Y position",
                  100, -14, 14, "mm");                             // ID 17
    man->CreateH1("Resolution_X_Paper",
                  "PAPER: X residual (x_rec - x_true)",
                  200, -3, 3, "mm");                             // ID 18
    man->CreateH1("Resolution_Y_Paper",
                  "PAPER: Y residual (y_rec - y_true)",
                  200, -3, 3, "mm");                             // ID 19
    man->CreateH1("Resolution_Mag_Paper",
                  "PAPER: Combined position resolution sqrt(dx²+dy²)",
                  200, 0, 20, "mm");                               // ID 20

    // 2D HISTOGRAMS
    // IDs 0-6  : algorithm 1
    // IDs 7-12 : PAPER algorithm

    // Algorithm 1
    man->CreateH2("Position_2D_Full",
                  "Algorithm 1: Full-event reconstructed XY position",
                  100, -14, 14, 100, -14, 14,
                  "mm", "mm");                                     // 2D ID 0

    man->CreateH2("Position_2D_Border_UNUSED",
                  "(unused – border reconstruction removed)",
                  1, -1, 1, 1, -1, 1);                            // 2D ID 1

    man->CreateH2("Theta_vs_ResX",
                  "Algorithm 1: Polar angle vs X residual;"
                  "#theta (deg);#Delta X (mm)",
                  90, 0, 90, 200, -5, 5);                         // 2D ID 2
    man->CreateH2("Theta_vs_ResY",
                  "Algorithm 1: Polar angle vs Y residual;"
                  "#theta (deg);#Delta Y (mm)",
                  90, 0, 90, 200, -5, 5);                         // 2D ID 3
    man->CreateH2("Theta_vs_Resolution",
                  "Algorithm 1: Polar angle vs combined resolution;"
                  "#theta (deg);Resolution (mm)",
                  90, 0, 90, 200, 0, 10);                         // 2D ID 4
    man->CreateH2("Phi_vs_ResX",
                  "Algorithm 1: Azimuthal angle vs X residual;"
                  "#phi (deg);#Delta X (mm)",
                  180, -180, 180, 200, -5, 5);                    // 2D ID 5
    man->CreateH2("Phi_vs_ResY",
                  "Algorithm 1: Azimuthal angle vs Y residual;"
                  "#phi (deg);#Delta Y (mm)",
                  180, -180, 180, 200, -5, 5);                    // 2D ID 6

    // PAPER algorithm
    man->CreateH2("Position_2D_Full_Paper",
                  "PAPER: Full-event reconstructed XY position",
                  100, -14, 14, 100, -14, 14,
                  "mm", "mm");                                     // 2D ID 7

    man->CreateH2("Theta_vs_ResX_Paper",
                  "PAPER: Polar angle vs X residual;"
                  "#theta (deg);#Delta X (mm)",
                  90, 0, 90, 200, -5, 5);                         // 2D ID 8
    man->CreateH2("Theta_vs_ResY_Paper",
                  "PAPER: Polar angle vs Y residual;"
                  "#theta (deg);#Delta Y (mm)",
                  90, 0, 90, 200, -5, 5);                         // 2D ID 9
    man->CreateH2("Theta_vs_Resolution_Paper",
                  "PAPER: Polar angle vs combined resolution;"
                  "#theta (deg);Resolution (mm)",
                  90, 0, 90, 200, 0, 10);                         // 2D ID 10
    man->CreateH2("Phi_vs_ResX_Paper",
                  "PAPER: Azimuthal angle vs X residual;"
                  "#phi (deg);#Delta X (mm)",
                  180, -180, 180, 200, -5, 5);                    // 2D ID 11
    man->CreateH2("Phi_vs_ResY_Paper",
                  "PAPER: Azimuthal angle vs Y residual;"
                  "#phi (deg);#Delta Y (mm)",
                  180, -180, 180, 200, -5, 5);                    // 2D ID 12

    // NTUPLES
    // IDs 0-3 : per-array sensor data
    // ID  4   : Reconstruction ntuple (placeholder, not filled)
    // ID  5   : AngleAnalysis        – algorithm 1
    // ID  6   : AngleAnalysis_Paper  – PAPER algorithm

    // Per-array sensor data
    for (int i = 0; i < 4; ++i) {
        std::string name  = "Array_" + std::to_string(i);
        std::string title = "Array "  + std::to_string(i) + " sensor data";
        man->CreateNtuple(name, title);                            // IDs 0-3
        man->CreateNtupleDColumn("posX");
        man->CreateNtupleDColumn("posY");
        man->CreateNtupleIColumn("event");
        man->CreateNtupleIColumn("copyNo");
        man->FinishNtuple();
    }

    // Reconstruction placeholder (ID 4, not filled)
    man->CreateNtuple("Reconstruction", "Placeholder – see AngleAnalysis ntuples");
    man->CreateNtupleDColumn("recX");
    man->CreateNtupleDColumn("recY");
    man->CreateNtupleDColumn("trueX");
    man->CreateNtupleDColumn("trueY");
    man->CreateNtupleIColumn("borderFlag");
    man->CreateNtupleIColumn("nArrays");
    man->FinishNtuple();                                           // ID 4

    // Algorithm 1 ntuple (ID 5)
    man->CreateNtuple("AngleAnalysis", "Algorithm 1 (weight = N/sigma^2)");
    man->CreateNtupleDColumn("x_rec");        // 0
    man->CreateNtupleDColumn("y_rec");        // 1
    man->CreateNtupleDColumn("x_true");       // 2
    man->CreateNtupleDColumn("y_true");       // 3
    man->CreateNtupleDColumn("dx");           // 4
    man->CreateNtupleDColumn("dy");           // 5
    man->CreateNtupleDColumn("theta");        // 6
    man->CreateNtupleDColumn("phi");          // 7
    man->CreateNtupleDColumn("resolution");   // 8
    man->FinishNtuple();                                           // ID 5

    // PAPER algorithm ntuple (ID 6)
    man->CreateNtuple("AngleAnalysis_Paper", "PAPER algorithm (weight = N/sigma)");
    man->CreateNtupleDColumn("x_rec");        // 0
    man->CreateNtupleDColumn("y_rec");        // 1
    man->CreateNtupleDColumn("x_true");       // 2
    man->CreateNtupleDColumn("y_true");       // 3
    man->CreateNtupleDColumn("dx");           // 4
    man->CreateNtupleDColumn("dy");           // 5
    man->CreateNtupleDColumn("theta");        // 6
    man->CreateNtupleDColumn("phi");          // 7
    man->CreateNtupleDColumn("resolution");   // 8
    man->FinishNtuple();                                           // ID 6
}

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
    auto* man = G4AnalysisManager::Instance();
    man->OpenFile("B4.root");
    G4cout << "RunAction: output file B4.root opened." << G4endl;
    G4cout << "Using analysis manager type: " << man->GetType() << G4endl;
}

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
    auto* man = G4AnalysisManager::Instance();
    man->Write();
    man->CloseFile();
    G4cout << "RunAction: B4.root written and closed." << G4endl;
}

}
