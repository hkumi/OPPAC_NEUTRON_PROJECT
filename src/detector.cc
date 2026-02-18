// ============================================================================
// detector.cc - COMPLETE IMPLEMENTATION WITH ANGLE CALCULATIONS 
// ============================================================================
// 
// KEY FEATURES:
// - Records INTERACTION position (where proton is created)
// - NOT entrance position (where neutron enters)
// - Calculates incident angle from particle direction
// - Fills angle-dependent resolution histograms
// - ENHANCED DEBUG OUTPUT
//
// ============================================================================

#include "detector.hh"
#include "EventAction.hh"
#include <vector>
#include <map>
#include <mutex>
#include <cmath>
#include "G4Poisson.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include "Randomize.hh"
#include "SensorHit.hh"

// --------------------------------------------------------------------
std::mutex MySensitiveDetector::analysisMutex;

// Thread-local variables for position tracking
G4ThreadLocal G4ThreeVector interactionPos;    // Where (n,p) reaction occurs
G4ThreadLocal G4ThreeVector primaryDirection;   // Particle initial direction  
G4ThreadLocal bool truePosSet = false;

// Debug counters (thread-local)
G4ThreadLocal int debugChargedParticleCount = 0;
G4ThreadLocal int debugCaptureAttempts = 0;

// --------------------------------------------------------------------
MySensitiveDetector::MySensitiveDetector(G4String name)
    : G4VSensitiveDetector(name),
      fPDE(1),           
      fNoiseLambda(6),  
      fMinSigma(0.3 * mm),
      fSensorPitch(1.10 * mm),
      fNSensorsPerArray(25),
      fEventAction(nullptr)  // Initialize to nullptr
{
    collectionName.insert("SensorCollection");
    
    G4cout << "==================================================" << G4endl;
    G4cout << "MySensitiveDetector Configuration:" << G4endl;
    G4cout << "  PDE: " << fPDE << G4endl;
    G4cout << "  Dark noise (λ): " << fNoiseLambda << " counts/event" << G4endl;
    G4cout << "  Min sigma: " << fMinSigma/mm << " mm" << G4endl;
    G4cout << "  Sensor pitch: " << fSensorPitch/mm << " mm" << G4endl;
    G4cout << "  Sensors per array: " << fNSensorsPerArray << G4endl;
    G4cout << "==================================================" << G4endl;
}

MySensitiveDetector::~MySensitiveDetector() {}

// --------------------------------------------------------------------
void MySensitiveDetector::Initialize(G4HCofThisEvent* HCE)
{
    SensorCollection =
        new SensorHitsCollection(SensitiveDetectorName, collectionName[0]);

    static G4int HCID = -1;
    if (HCID < 0)
        HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);

    HCE->AddHitsCollection(HCID, SensorCollection);
    
    // Reset debug counters at start of event
    debugChargedParticleCount = 0;
    debugCaptureAttempts = 0;
}

// --------------------------------------------------------------------
double MySensitiveDetector::IndexToPosition(double index)
{
    // Convert sensor index to physical position
    // For 25 sensors with 1.1mm pitch centered at 0
    return (index + 0.5) * fSensorPitch
           - (fNSensorsPerArray * fSensorPitch) / 2.0;
}

// --------------------------------------------------------------------
void MySensitiveDetector::RecordSensorData(int ntupleIndex,
                                           double x, double y,
                                           int event, int copyNo)
{
    std::lock_guard<std::mutex> lock(analysisMutex);
    auto* man = G4AnalysisManager::Instance();

    man->FillNtupleDColumn(ntupleIndex, 0, x / mm);
    man->FillNtupleDColumn(ntupleIndex, 1, y / mm);
    man->FillNtupleIColumn(ntupleIndex, 2, event);
    man->FillNtupleIColumn(ntupleIndex, 3, copyNo);
    man->AddNtupleRow(ntupleIndex);
}

// --------------------------------------------------------------------
void MySensitiveDetector::ReconstructFullEvent(
        const std::vector<double>& xTop,
        const std::vector<double>& xBottom,
        const std::vector<double>& yLeft,
        const std::vector<double>& yRight)
{
    std::lock_guard<std::mutex> lock(analysisMutex);
    auto* man = G4AnalysisManager::Instance();

    // Check if EventAction is available
    if (!fEventAction) {
        G4cout << "ERROR: fEventAction is null in ReconstructFullEvent!" << G4endl;
        return;
    }

    // Gaussian MLE for each array
    auto GaussianMLE = [&](const std::vector<double>& v,
                           double& mean, double& sigma, double& N)
    {
        N = v.size();
        if (N == 0) {
            mean = 0.0;
            sigma = fMinSigma;
            return;
        }

        double sum = 0.0;
        for (double i : v) sum += IndexToPosition(i);
        mean = sum / N;

        if (N < 2) {
            sigma = fMinSigma;
            return;
        }

        double var = 0.0;
        for (double i : v) {
            double x = IndexToPosition(i);
            var += (x - mean) * (x - mean);
        }

        sigma = std::sqrt(var / N);
        if (sigma < fMinSigma) sigma = fMinSigma;
    };

    double mxT, sxT, nT;
    double mxB, sxB, nB;
    double myL, syL, nL;
    double myR, syR, nR;

    GaussianMLE(xTop,    mxT, sxT, nT);
    GaussianMLE(xBottom, mxB, sxB, nB);
    GaussianMLE(yLeft,   myL, syL, nL);
    GaussianMLE(yRight,  myR, syR, nR);

    // Weighted MLE combining arrays
    auto WeightedMLE = [](double p1, double n1, double s1,
                          double p2, double n2, double s2)
    {
        if (n1 > 0 && n2 > 0)
            return (p1 * n1 / (s1) + p2 * n2 / (s2)) /
                   (n1 / (s1)     + n2 / (s2));
        else if (n1 > 0) return p1;
        else if (n2 > 0) return p2;
        return 0.0;
    };

    double x_rec = WeightedMLE(mxT, nT, sxT, mxB, nB, sxB);
    double y_rec = WeightedMLE(myL, nL, syL, myR, nR, syR);

    // Fill reconstructed position histograms
    man->FillH1(6, x_rec);
    man->FillH1(7, y_rec);
    man->FillH2(0, x_rec, y_rec);

    // ========================================================================
    // CALCULATE RESIDUALS using INTERACTION position from EventAction
    // ========================================================================
    G4ThreeVector trueInteractionPos = fEventAction->GetInteractionPosition();
    G4ThreeVector truePrimaryDirection = fEventAction->GetPrimaryDirection();
    
    double x_true = trueInteractionPos.x();
    double y_true = trueInteractionPos.y();

    double dx = x_rec - x_true;
    double dy = y_rec - y_true;
    double resolution = std::sqrt(dx*dx + dy*dy);

    // ========================================================================
// DEBUG: Print reconstruction for ALL events with non-zero true positions
// ========================================================================
// DEBUG: Always print
    G4cout << "\n=== RECONSTRUCTION CALLED ===" << G4endl;
    G4cout << "TRUE position: (" << x_true/mm << ", " << y_true/mm << ") mm" << G4endl;
    G4cout << "Sensors hit - Top:" << xTop.size() 
           << ", Bottom:" << xBottom.size()
           << ", Left:" << yLeft.size() 
           << ", Right:" << yRight.size() << G4endl;
    
    // Check if we have enough sensors
    int arrays = (!xTop.empty()) + (!xBottom.empty()) +
                 (!yLeft.empty()) + (!yRight.empty());
    
    G4cout << "Arrays with signal: " << arrays << G4endl;
    
    if (arrays >= 4) {
        G4cout << "Will do FULL reconstruction" << G4endl;
    } else if (arrays >= 2) {
        G4cout << "Will do BORDER reconstruction" << G4endl;
    } else {
        G4cout << "NOT ENOUGH ARRAYS - NO RECONSTRUCTION" << G4endl;
    }
    
    
    // Fill residual histograms
    man->FillH1(8, dx / mm);   // ΔX
    man->FillH1(9, dy / mm);   // ΔY
    man->FillH1(13, resolution / mm);  // Combined resolution

    // ========================================================================
    // CALCULATE INCIDENT ANGLE from EventAction
    // ========================================================================
    double px = truePrimaryDirection.x();
    double py = truePrimaryDirection.y();
    double pz = truePrimaryDirection.z();
    
    // Polar angle from Z-axis (detector normal)
    double theta = std::acos(std::abs(pz)) * 180.0 / M_PI;
    
    // Azimuthal angle in X-Y plane
    double phi = std::atan2(py, px) * 180.0 / M_PI;
    
    // ========================================================================
    // DEBUG: Print angle details
    // ========================================================================
   /* G4cout << "Primary direction (from EventAction): (" << px << ", " << py << ", " << pz << ")" << G4endl;
    G4cout << "Theta: " << theta << " degrees" << G4endl;
    G4cout << "Phi: " << phi << " degrees" << G4endl;
    G4cout << "****************************\n" << G4endl;*/
    
    // Fill angle-dependent histograms
    man->FillH2(2, theta, dx );          // H2(2): theta vs ΔX
    man->FillH2(3, theta, dy );          // H2(3): theta vs ΔY
    man->FillH2(4, theta, resolution );  // H2(4): theta vs resolution
    man->FillH2(5, phi, dx );            // H2(5): phi vs ΔX
    man->FillH2(6, phi, dy );            // H2(6): phi vs ΔY

    // ========================================================================
    // FILL NTUPLES
    // ========================================================================
    
    /*// Fill Reconstruction ntuple (ID 4) - Full events
    man->FillNtupleDColumn(4, 0, x_rec );
    man->FillNtupleDColumn(4, 1, y_rec );
    man->FillNtupleDColumn(4, 2, x_true );
    man->FillNtupleDColumn(4, 3, y_true);
    man->FillNtupleIColumn(4, 4, 0);  // borderFlag = 0 (full event)
    man->FillNtupleIColumn(4, 5, 4);  // nArrays = 4
    man->AddNtupleRow(4);*/

    // Fill AngleAnalysis ntuple (ID 5)
    man->FillNtupleDColumn(5, 0, x_rec );
    man->FillNtupleDColumn(5, 1, y_rec );
    man->FillNtupleDColumn(5, 2, x_true );
    man->FillNtupleDColumn(5, 3, y_true );
    man->FillNtupleDColumn(5, 4, dx );
    man->FillNtupleDColumn(5, 5, dy );
    man->FillNtupleDColumn(5, 6, theta);
    man->FillNtupleDColumn(5, 7, phi);
    man->FillNtupleDColumn(5, 8, resolution );
    man->AddNtupleRow(5);
}

// --------------------------------------------------------------------
void MySensitiveDetector::ReconstructBorderEvent(
        const std::vector<double>& xTop,
        const std::vector<double>& xBottom,
        const std::vector<double>& yLeft,
        const std::vector<double>& yRight,
        int nArrays)
{
    std::lock_guard<std::mutex> lock(analysisMutex);
    auto* man = G4AnalysisManager::Instance();

    // Check if EventAction is available
    if (!fEventAction) {
        G4cout << "ERROR: fEventAction is null in ReconstructBorderEvent!" << G4endl;
        return;
    }

    auto GaussianMLE = [&](const std::vector<double>& v,
                           double& mean, double& sigma, double& N)
    {
        N = v.size();
        if (N == 0) {
            mean = 0.0;
            sigma = fMinSigma;
            return;
        }

        double sum = 0.0;
        for (double i : v) sum += IndexToPosition(i);
        mean = sum / N;

        if (N < 2) {
            sigma = fMinSigma;
            return;
        }

        double var = 0.0;
        for (double i : v) {
            double x = IndexToPosition(i);
            var += (x - mean) * (x - mean);
        }

        sigma = std::sqrt(var / N);
        if (sigma < fMinSigma) sigma = fMinSigma;
    };

    double mxT, sxT, nT;
    double mxB, sxB, nB;
    double myL, syL, nL;
    double myR, syR, nR;

    GaussianMLE(xTop,    mxT, sxT, nT);
    GaussianMLE(xBottom, mxB, sxB, nB);
    GaussianMLE(yLeft,   myL, syL, nL);
    GaussianMLE(yRight,  myR, syR, nR);

    auto WeightedMLE = [](double p1, double n1, double s1,
                          double p2, double n2, double s2)
    {
        if (n1 > 0 && n2 > 0)
            return (p1 * n1 / (s1) + p2 * n2 / (s2)) /
                   (n1 / (s1)     + n2 / (s2));
        else if (n1 > 0) return p1;
        else if (n2 > 0) return p2;
        return 0.0;
    };

    double x_rec = WeightedMLE(mxT, nT, sxT, mxB, nB, sxB);
    double y_rec = WeightedMLE(myL, nL, syL, myR, nR, syR);

    man->FillH2(1, x_rec, y_rec);

    // Get true position from EventAction
    G4ThreeVector trueInteractionPos = fEventAction->GetInteractionPosition();
    double x_true = trueInteractionPos.x();
    double y_true = trueInteractionPos.y();
    
    double dx = x_rec - x_true;
    double dy = y_rec - y_true;
    
    man->FillH1(14, dx / mm);  // Border_Res_X
    man->FillH1(15, dy / mm);  // Border_Res_Y

   /* // Fill Reconstruction ntuple (ID 4) - Border events
    man->FillNtupleDColumn(4, 0, x_rec / mm);
    man->FillNtupleDColumn(4, 1, y_rec / mm);
    man->FillNtupleDColumn(4, 2, x_true / mm);
    man->FillNtupleDColumn(4, 3, y_true / mm);
    man->FillNtupleIColumn(4, 4, 1);  // borderFlag = 1 (border event)
    man->FillNtupleIColumn(4, 5, nArrays);
    man->AddNtupleRow(4);*/
}

// --------------------------------------------------------------------
// --------------------------------------------------------------------
G4bool MySensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    G4Track* track = step->GetTrack();

    // Only process optical photons
    if (track->GetDefinition() != G4OpticalPhoton::Definition())
        return false;

    auto* touchable = step->GetPreStepPoint()->GetTouchable();
    if (!touchable) return false;

    G4int copyNo = touchable->GetCopyNumber();
    G4ThreeVector pos = step->GetPostStepPoint()->GetPosition();

    // Validate copy number
    if (copyNo < 0 || copyNo >= 4 * fNSensorsPerArray) {
        track->SetTrackStatus(fStopAndKill);
        return false;
    }

    // Create sensor hit
    auto* hit = new SensorHit();
    hit->SetSensorPosition(pos);
    hit->SetSensorEnergy(track->GetKineticEnergy());
    hit->SetCopyNumber(copyNo);
    SensorCollection->insert(hit);

    // Record in ntuple
    int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    int N = fNSensorsPerArray;
    
    if (copyNo < N)          RecordSensorData(0, pos.x(), pos.y(), evt, copyNo);
    else if (copyNo < 2*N)   RecordSensorData(1, pos.x(), pos.y(), evt, copyNo);
    else if (copyNo < 3*N)   RecordSensorData(2, pos.x(), pos.y(), evt, copyNo);
    else if (copyNo < 4*N)   RecordSensorData(3, pos.x(), pos.y(), evt, copyNo);

    track->SetTrackStatus(fStopAndKill);
    return true;
}

// --------------------------------------------------------------------
void MySensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
    if (!SensorCollection || SensorCollection->entries() == 0)
        return;

    // Count photons per array
    std::map<int,int> topC, bottomC, leftC, rightC;
    int N = fNSensorsPerArray;

    for (int i = 0; i < SensorCollection->entries(); ++i) {
        int c = (*SensorCollection)[i]->GetCopyNumber();
        if (c < N)        topC[c]++;
        else if (c < 2*N) bottomC[c - N]++;
        else if (c < 3*N) leftC[c - 2*N]++;
        else if (c < 4*N) rightC[c - 3*N]++;
    }

    // Apply PDE and dark noise
    auto build = [&](const std::map<int,int>& C, std::vector<double>& V)
    {
        for (auto& [idx, ph] : C) {
            int npe = G4Poisson(ph * fPDE);
            for (int i = 0; i < npe; ++i) V.push_back(idx);
        }
        int noise = G4Poisson(fNoiseLambda);
        for (int i = 0; i < noise; ++i)
            V.push_back(static_cast<int>(G4UniformRand() * N));
    };

    std::vector<double> topI, botI, leftI, rightI;
    build(topC, topI);
    build(bottomC, botI);
    build(leftC, leftI);
    build(rightC, rightI);

    // Count arrays with signal
    int arrays = (!topI.empty()) + (!botI.empty()) +
                 (!leftI.empty()) + (!rightI.empty());

    // Reconstruct position
    if (arrays >= 4) {
        ReconstructFullEvent(topI, botI, leftI, rightI);
    }
   /* else if (arrays >= 2) {
        ReconstructBorderEvent(topI, botI, leftI, rightI, arrays);
    }*/
}
