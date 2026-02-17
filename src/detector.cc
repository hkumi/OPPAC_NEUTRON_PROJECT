// ============================================================================
// detector.cc - COMPLETE IMPLEMENTATION WITH ANGLE CALCULATIONS
// Records TRUE INTERACTION position from EventAction
// Calculates angle-dependent resolution
// ============================================================================

#include "detector.hh"
#include "EventAction.hh"
#include "SensorHit.hh"

#include <vector>
#include <map>
#include <mutex>
#include <cmath>

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"

// ============================================================================
// Thread Safety
// ============================================================================
std::mutex MySensitiveDetector::analysisMutex;

// ============================================================================
// Constructor
// ============================================================================
MySensitiveDetector::MySensitiveDetector(G4String name)
    : G4VSensitiveDetector(name),
      fPDE(1.0),
      fNoiseLambda(6),
      fMinSigma(0.3 * mm),
      fSensorPitch(1.10 * mm),
      fNSensorsPerArray(25),
      fEventAction(nullptr)
{
    collectionName.insert("SensorCollection");

    G4cout << "\n==== MySensitiveDetector Initialized ====" << G4endl;
    G4cout << "PDE: " << fPDE << G4endl;
    G4cout << "Dark Noise: " << fNoiseLambda << " counts/event" << G4endl;
    G4cout << "Sensor pitch: " << fSensorPitch/mm << " mm" << G4endl;
    G4cout << "Sensors/array: " << fNSensorsPerArray << G4endl;
    G4cout << "=========================================\n" << G4endl;
}

MySensitiveDetector::~MySensitiveDetector() {}

// ============================================================================
// Initialize
// ============================================================================
void MySensitiveDetector::Initialize(G4HCofThisEvent* HCE)
{
    SensorCollection =
        new SensorHitsCollection(SensitiveDetectorName, collectionName[0]);

    static G4int HCID = -1;
    if (HCID < 0)
        HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);

    HCE->AddHitsCollection(HCID, SensorCollection);
}

// ============================================================================
// Convert sensor index to position
// ============================================================================
double MySensitiveDetector::IndexToPosition(double index)
{
    return (index + 0.5) * fSensorPitch
           - (fNSensorsPerArray * fSensorPitch) / 2.0;
}

// ============================================================================
// Record per-sensor hit in ntuple
// ============================================================================
void MySensitiveDetector::RecordSensorData(int ntupleIndex,
                                           double x, double y,
                                           int event, int copyNo)
{
    std::lock_guard<std::mutex> lock(analysisMutex);
    auto* man = G4AnalysisManager::Instance();

    man->FillNtupleDColumn(ntupleIndex, 0, x/mm);
    man->FillNtupleDColumn(ntupleIndex, 1, y/mm);
    man->FillNtupleIColumn(ntupleIndex, 2, event);
    man->FillNtupleIColumn(ntupleIndex, 3, copyNo);
    man->AddNtupleRow(ntupleIndex);
}

// ============================================================================
// FULL EVENT RECONSTRUCTION
// ============================================================================
void MySensitiveDetector::ReconstructFullEvent(
        const std::vector<double>& xTop,
        const std::vector<double>& xBottom,
        const std::vector<double>& yLeft,
        const std::vector<double>& yRight)
{
    if (!fEventAction) return;

    std::lock_guard<std::mutex> lock(analysisMutex);
    auto* man = G4AnalysisManager::Instance();

    // ---------- Gaussian MLE ----------
    auto GaussianMLE = [&](const std::vector<double>& v,
                           double& mean, double& sigma, double& N)
    {
        N = v.size();
        if (N == 0) { mean = 0; sigma = fMinSigma; return; }

        double sum = 0;
        for (double i : v) sum += IndexToPosition(i);
        mean = sum / N;

        if (N < 2) { sigma = fMinSigma; return; }

        double var = 0;
        for (double i : v) {
            double x = IndexToPosition(i);
            var += (x - mean)*(x - mean);
        }
        sigma = std::sqrt(var/N);
        if (sigma < fMinSigma) sigma = fMinSigma;
    };

    double mxT,sxT,nT, mxB,sxB,nB;
    double myL,syL,nL, myR,syR,nR;

    GaussianMLE(xTop,mxT,sxT,nT);
    GaussianMLE(xBottom,mxB,sxB,nB);
    GaussianMLE(yLeft,myL,syL,nL);
    GaussianMLE(yRight,myR,syR,nR);

    auto WeightedMLE = [](double p1,double n1,double s1,
                          double p2,double n2,double s2)
    {
        if (n1>0 && n2>0)
            return (p1*n1/s1 + p2*n2/s2)/(n1/s1 + n2/s2);
        if (n1>0) return p1;
        if (n2>0) return p2;
        return 0.0;
    };

    double x_rec = WeightedMLE(mxT,nT,sxT,mxB,nB,sxB);
    double y_rec = WeightedMLE(myL,nL,syL,myR,nR,syR);

    // ---------- TRUE INTERACTION ----------
    G4ThreeVector truePos = fEventAction->GetInteractionPosition();
    G4ThreeVector trueDir = fEventAction->GetPrimaryDirection();

    double dx = x_rec - truePos.x();
    double dy = y_rec - truePos.y();
    double resolution = std::sqrt(dx*dx + dy*dy);

    // ---------- Histograms ----------
    man->FillH1(6, x_rec/mm);
    man->FillH1(7, y_rec/mm);
    man->FillH2(0, x_rec/mm, y_rec/mm);

    man->FillH1(8, dx/mm);
    man->FillH1(9, dy/mm);
    man->FillH1(13, resolution/mm);

    // ---------- ANGLE ----------
    double theta = std::acos(std::abs(trueDir.z())) * 180.0/M_PI;
    double phi   = std::atan2(trueDir.y(), trueDir.x()) * 180.0/M_PI;

    // Fill angle-dependent histograms
    man->FillH2(2, theta, dx/mm );          // H2(2): theta vs ΔX
    man->FillH2(3, theta, dy/mm );          // H2(3): theta vs ΔY
    man->FillH2(4, theta, resolution/mm );  // H2(4): theta vs resolution
    man->FillH2(5, phi, dx/mm );            // H2(5): phi vs ΔX
    man->FillH2(6, phi, dy/mm );            // H2(6): phi vs ΔY

    // ---------- Angle ntuple (ID 5) ----------
    man->FillNtupleDColumn(5,0,x_rec/mm);
    man->FillNtupleDColumn(5,1,y_rec/mm);
    man->FillNtupleDColumn(5,2,truePos.x()/mm);
    man->FillNtupleDColumn(5,3,truePos.y()/mm);
    man->FillNtupleDColumn(5,4,dx/mm);
    man->FillNtupleDColumn(5,5,dy/mm);
    man->FillNtupleDColumn(5, 6, theta);
    man->FillNtupleDColumn(5, 7, phi);
    man->FillNtupleDColumn(5,8,resolution/mm);
    man->AddNtupleRow(5);
}

// ============================================================================
// BORDER RECONSTRUCTION
// ============================================================================
void MySensitiveDetector::ReconstructBorderEvent(
        const std::vector<double>& xTop,
        const std::vector<double>& xBottom,
        const std::vector<double>& yLeft,
        const std::vector<double>& yRight,
        int)
{
    /*if (!fEventAction) return;

    std::lock_guard<std::mutex> lock(analysisMutex);
    auto* man = G4AnalysisManager::Instance();

    double x_rec = 0, y_rec = 0;

    if (!xTop.empty())   x_rec = IndexToPosition(xTop.front());
    if (!yLeft.empty())  y_rec = IndexToPosition(yLeft.front());

    G4ThreeVector truePos = fEventAction->GetInteractionPosition();

    man->FillH2(1, x_rec/mm, y_rec/mm);
    man->FillH1(14, (x_rec-truePos.x())/mm);
    man->FillH1(15, (y_rec-truePos.y())/mm);*/
}

// ============================================================================
// PROCESS HITS (Optical photons only)
// ============================================================================
G4bool MySensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    G4Track* track = step->GetTrack();

    if (track->GetDefinition() != G4OpticalPhoton::Definition())
        return false;

    auto* touchable = step->GetPreStepPoint()->GetTouchable();
    if (!touchable) return false;

    int copyNo = touchable->GetCopyNumber();
    if (copyNo < 0 || copyNo >= 4*fNSensorsPerArray) {
        track->SetTrackStatus(fStopAndKill);
        return false;
    }

    auto* hit = new SensorHit();
    hit->SetSensorPosition(step->GetPostStepPoint()->GetPosition());
    hit->SetSensorEnergy(track->GetKineticEnergy());
    hit->SetCopyNumber(copyNo);
    SensorCollection->insert(hit);

    track->SetTrackStatus(fStopAndKill);
    return true;
}

// ============================================================================
// END OF EVENT
// ============================================================================
void MySensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
    if (!SensorCollection || SensorCollection->entries()==0)
        return;

    std::map<int,int> top,bottom,left,right;
    int N = fNSensorsPerArray;

    for (int i=0;i<SensorCollection->entries();++i) {
        int c = (*SensorCollection)[i]->GetCopyNumber();
        if (c<N) top[c]++;
        else if (c<2*N) bottom[c-N]++;
        else if (c<3*N) left[c-2*N]++;
        else right[c-3*N]++;
    }

    auto build = [&](const std::map<int,int>& C, std::vector<double>& V){
        for (auto& [idx,ph] : C) {
            int npe = G4Poisson(ph*fPDE);
            for(int i=0;i<npe;++i) V.push_back(idx);
        }
        int noise = G4Poisson(fNoiseLambda);
        for(int i=0;i<noise;++i)
            V.push_back((int)(G4UniformRand()*N));
    };

    std::vector<double> t,b,l,r;
    build(top,t); build(bottom,b);
    build(left,l); build(right,r);

    int arrays = (!t.empty())+(!b.empty())+(!l.empty())+(!r.empty());

    if (arrays==4)
        ReconstructFullEvent(t,b,l,r);
   // else if (arrays>=2)
        //ReconstructBorderEvent(t,b,l,r,arrays);
}

