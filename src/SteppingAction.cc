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

// For Gaussian fitting
#include "TH1D.h"
#include "TF1.h"
#include <algorithm>
#include <cmath>
#include <fstream>


namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(const B4::DetectorConstruction* detConstruction,
                               EventAction* eventAction)
  : fDetConstruction(detConstruction),
    fEventAction(eventAction),
    fCurrentEventID(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::SaveEventData(G4int eventID)
{
    std::ofstream outfile("position_data.txt", std::ios_base::app);
    if (!outfile) return;
    
    outfile << eventID << " "
            << fXBottom.size() << " "
            << fXTop.size() << " "
            << fYLeft.size() << " "
            << fYRight.size() << " ";
    
    for (double x : fXBottom) outfile << x << " ";
    for (double x : fXTop) outfile << x << " ";
    for (double y : fYLeft) outfile << y << " ";
    for (double y : fYRight) outfile << y << " ";
    
    outfile << "\n";
    outfile.close();
}

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

    // Check if this is a new event
    if (eventID != fCurrentEventID) {
        // Save previous event data if we have any
        if (fCurrentEventID >= 0 && (!fXBottom.empty() || !fXTop.empty() || !fYLeft.empty() || !fYRight.empty())) {
            SaveEventData(fCurrentEventID);
        }
        // Reset for new event
        fCurrentEventID = eventID;
        fXBottom.clear();
        fXTop.clear();
        fYLeft.clear();
        fYRight.clear();
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

        // DEBUG OUTPUT (optional)
        // G4cout << "\n=== PHOTON HIT ===" << G4endl;
        // G4cout << "Event ID: " << eventID << G4endl;
        // G4cout << "Sensor copyNo: " << copyNo << G4endl;
        
        // Fill detailed sensor data ntuple (ntuple 0)
        analysisManager->FillNtupleIColumn(0, 0, eventID);              // EventID
        analysisManager->FillNtupleIColumn(0, 1, copyNo);               // SensorID
        analysisManager->FillNtupleDColumn(0, 2, posX);                 // PosX (mm)
        analysisManager->FillNtupleDColumn(0, 3, posY);                 // PosY (mm)
        analysisManager->FillNtupleDColumn(0, 4, photonEnergy/eV);      // Energy (eV)
        analysisManager->FillNtupleDColumn(0, 5, time/ns);              // Time (ns)
        analysisManager->FillNtupleIColumn(0, 6, 1);                    // PhotonCount
        analysisManager->AddNtupleRow(0);
        
        // Update histograms
        analysisManager->FillH1(0, photonEnergy);

        // Store positions based on sensor location
        if (copyNo < 25) { 
            fXBottom.push_back(posX);
            analysisManager->FillH1(1, copyNo);
        }
        else if (copyNo < 50) {
            fYLeft.push_back(posY); 
            analysisManager->FillH1(2, copyNo - 25);
        }
        else if (copyNo < 75) {
            fXTop.push_back(posX);  
            analysisManager->FillH1(3, copyNo - 50);
        }
        else if (copyNo < 100) {
            fYRight.push_back(posY); 
            analysisManager->FillH1(4, copyNo - 75);
        }

        // Kill the photon
        step->GetTrack()->SetTrackStatus(fStopAndKill); 
    }

    // Proton tracking
    if (particleName == "proton") {
        const G4VProcess* preProcess = step->GetPreStepPoint()->GetProcessDefinedStep();
        G4String preProcessName = preProcess ? preProcess->GetProcessName() : "none";
        
        if (preProcessName == "none") {
            analysisManager->FillH1(5, energy);
        }
        
        if (volName == "MylA") {
            analysisManager->FillH1(6, energy);
        }
        
        if (volName == "World") {
            analysisManager->FillH1(7, energy);
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::FillReconstructedPosition(G4int eventID)
{
    auto analysisManager = G4AnalysisManager::Instance();
    
    // Get sizes
    G4int nBottom = fXBottom.size();
    G4int nTop = fXTop.size();
    G4int nLeft = fYLeft.size();
    G4int nRight = fYRight.size();
    G4int nTotal = nBottom + nTop + nLeft + nRight;
    
    // Skip events with no photons
    if (nTotal == 0) return;
    
    // Function to calculate Gaussian stats
    auto getGaussianStats = [eventID](const std::vector<G4double>& positions,
                                     G4double& mean, G4double& sigma, G4double& N) {
        N = positions.size();
        if (N == 0) {
            mean = 0.0;
            sigma = 1.0;
            return;
        }
        
        if (N < 5) {
            // Simple stats
            G4double sum = 0.0;
            for (G4double pos : positions) sum += pos;
            mean = sum / N;
            
            G4double sumSq = 0.0;
            for (G4double pos : positions) sumSq += (pos - mean) * (pos - mean);
            sigma = std::sqrt(sumSq / N);
            if (sigma < 0.1) sigma = 0.1;
            return;
        }
        
        // Gaussian fit
        G4double minVal = *std::min_element(positions.begin(), positions.end());
        G4double maxVal = *std::max_element(positions.begin(), positions.end());
        G4double range = maxVal - minVal;
        
        minVal -= 0.1 * range;
        maxVal += 0.1 * range;
        
        TH1D* hist = new TH1D(("h_temp_" + std::to_string(eventID)).c_str(), 
                              "Temp histogram", 100, minVal, maxVal);
        
        for (G4double pos : positions) hist->Fill(pos);
        
        TF1* gauss = new TF1("gauss", "gaus", minVal, maxVal);
        gauss->SetParameters(N/10.0, (minVal + maxVal)/2.0, range/4.0);
        hist->Fit(gauss, "QN");
        
        mean = gauss->GetParameter(1);
        sigma = gauss->GetParameter(2);
        
        delete hist;
        delete gauss;
        
        if (sigma <= 0 || !std::isfinite(sigma)) {
            // Fallback
            G4double sum = 0.0;
            for (G4double pos : positions) sum += pos;
            mean = sum / N;
            
            G4double sumSq = 0.0;
            for (G4double pos : positions) sumSq += (pos - mean) * (pos - mean);
            sigma = std::sqrt(sumSq / N);
            if (sigma < 0.1) sigma = 0.1;
        }
    };
    
    // Calculate stats
    G4double Px1, sigmaX1, Nx1;  // Bottom
    G4double Px2, sigmaX2, Nx2;  // Top
    G4double Py1, sigmaY1, Ny1;  // Left
    G4double Py2, sigmaY2, Ny2;  // Right
    
    getGaussianStats(fXBottom, Px1, sigmaX1, Nx1);
    getGaussianStats(fXTop, Px2, sigmaX2, Nx2);
    getGaussianStats(fYLeft, Py1, sigmaY1, Ny1);
    getGaussianStats(fYRight, Py2, sigmaY2, Ny2);
    
    // Apply paper formula
    G4double xWeighted = 0.0, yWeighted = 0.0;
    
    // X reconstruction
    G4double numeratorX = 0.0;
    G4double denominatorX = 0.0;
    
    if (Nx1 > 0 && sigmaX1 > 0) {
        numeratorX += Px1 * Nx1 / sigmaX1;
        denominatorX += Nx1 / sigmaX1;
    }
    if (Nx2 > 0 && sigmaX2 > 0) {
        numeratorX += Px2 * Nx2 / sigmaX2;
        denominatorX += Nx2 / sigmaX2;
    }
    if (denominatorX > 0) {
        xWeighted = numeratorX / denominatorX;
    }
    
    // Y reconstruction
    G4double numeratorY = 0.0;
    G4double denominatorY = 0.0;
    
    if (Ny1 > 0 && sigmaY1 > 0) {
        numeratorY += Py1 * Ny1 / sigmaY1;
        denominatorY += Ny1 / sigmaY1;
    }
    if (Ny2 > 0 && sigmaY2 > 0) {
        numeratorY += Py2 * Ny2 / sigmaY2;
        denominatorY += Ny2 / sigmaY2;
    }
    if (denominatorY > 0) {
        yWeighted = numeratorY / denominatorY;
    }
    
    // Fill reconstructed positions ntuple (ntuple 1)
    analysisManager->FillNtupleIColumn(1, 0, eventID);
    analysisManager->FillNtupleDColumn(1, 1, xWeighted);
    analysisManager->FillNtupleDColumn(1, 2, yWeighted);
    analysisManager->FillNtupleDColumn(1, 3, nTotal);
    analysisManager->AddNtupleRow(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
