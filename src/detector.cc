// ============================================================================
// detector.cc
//
// Two reconstruction algorithms run every event for direct comparison:
//
//   ReconstructFullEvent       → weight = N / σ²
//   ReconstructFullEventPaper  → PAPER Eq.1 weight = N / σ
//
// Histogram / ntuple layout:
//
//   H1 IDs 0-15   : algorithm uses IDs 6-9, 13
//   H1 IDs 16-20  : paper algorithm equivalents
//     16 = Weighted_X_Paper
//     17 = Weighted_Y_Paper
//     18 = Resolution_X_Paper
//     19 = Resolution_Y_Paper
//     20 = Resolution_Mag_Paper
//
//   H2 IDs 0-6    : algorithm
//   H2 IDs 7-12   : paper algorithm equivalents
//     7  = Position_2D_Full_Paper
//     8  = Theta_vs_ResX_Paper
//     9  = Theta_vs_ResY_Paper
//     10 = Theta_vs_Resolution_Paper
//     11 = Phi_vs_ResX_Paper
//     12 = Phi_vs_ResY_Paper
//
//   Ntuple ID 5   = AngleAnalysis
//   Ntuple ID 6   = AngleAnalysis_Paper
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

std::mutex MySensitiveDetector::analysisMutex;

MySensitiveDetector::MySensitiveDetector(G4String name)
    : G4VSensitiveDetector(name),
      fPDE(1.0),
      fNoiseLambda(0.5),
      fSensorPitch(1.10 * mm),
      fMinSigma(0.3),
      fNSensorsPerArray(25),
      fEventAction(nullptr)
{
    collectionName.insert("SensorCollection");

    G4cout << "==================================================" << G4endl;
    G4cout << "MySensitiveDetector – dual-algorithm mode" << G4endl;
    G4cout << "  Algorithm 1 : weight = N / sigma^2" << G4endl;
    G4cout << "  PAPER Eq.1  : weight = N / sigma" << G4endl;
    G4cout << "  PDE            : " << fPDE << G4endl;
    G4cout << "  Dark noise (λ) : " << fNoiseLambda << G4endl;
    G4cout << "  Min sigma      : " << fMinSigma / mm << " mm" << G4endl;
    G4cout << "  Sensor pitch   : " << fSensorPitch / mm << " mm" << G4endl;
    G4cout << "  Sensors/array  : " << fNSensorsPerArray << G4endl;
    G4cout << "==================================================" << G4endl;
}

MySensitiveDetector::~MySensitiveDetector() {}

void MySensitiveDetector::Initialize(G4HCofThisEvent* HCE)
{
    SensorCollection =
        new SensorHitsCollection(SensitiveDetectorName, collectionName[0]);

    static G4int HCID = -1;
    if (HCID < 0)
        HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);

    HCE->AddHitsCollection(HCID, SensorCollection);
}

double MySensitiveDetector::IndexToPosition(double index)
{
    // index 0 → leftmost sensor centre, index N-1 → rightmost
    return (index + 0.5) * fSensorPitch
           - (fNSensorsPerArray * fSensorPitch) / 2.0;
}

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

// Extracts P (mean), N (total pe count), σ (std dev) from the list of
// photoelectron sensor-index hits. Maps directly to quantities in Eq.1.
void MySensitiveDetector::ComputeGaussianParams(
        const std::vector<double>& hits,
        double& P, double& N, double& sigma)
{
    N = static_cast<double>(hits.size());

    if (N == 0) {
        P     = 0.0;
        sigma = fMinSigma;
        return;
    }

    double sum = 0.0;
    for (double idx : hits) sum += IndexToPosition(idx);
    P = sum / N;

    if (N < 2) {
        sigma = fMinSigma;
        return;
    }

    double var = 0.0;
    for (double idx : hits) {
        double pos = IndexToPosition(idx);
        var += (pos - P) * (pos - P);
    }
    sigma = std::sqrt(var / N);
    if (sigma < fMinSigma) sigma = fMinSigma;
}

// Algorithm 1: weight = N / σ² (inverse-variance MLE)
double MySensitiveDetector::WeightedMLE(double p1, double n1, double s1,
                                        double p2, double n2, double s2)
{
    if (n1 > 0 && n2 > 0)
        return (p1 * n1 / (s1*s1) + p2 * n2 / (s2*s2)) /
               (      n1 / (s1*s1) +       n2 / (s2*s2));
    else if (n1 > 0) return p1;
    else if (n2 > 0) return p2;
    return 0.0;
}

// PAPER Eq.1: weight = N / σ
double MySensitiveDetector::PaperEq1(double p1, double n1, double s1,
                                     double p2, double n2, double s2)
{
    double saf1 = (s1 > 0.0) ? s1 : fMinSigma;
    double saf2 = (s2 > 0.0) ? s2 : fMinSigma;

    double w1 = (n1 > 0.0) ? n1 / saf1 : 0.0;
    double w2 = (n2 > 0.0) ? n2 / saf2 : 0.0;

    double denom = w1 + w2;
    if (denom > 0.0)  return (w1 * p1 + w2 * p2) / denom;
    else if (w1 > 0.0) return p1;
    else if (w2 > 0.0) return p2;
    return 0.0;
}

static void ComputeAngles(const G4ThreeVector& dir,
                          double& theta, double& phi)
{
    double px = dir.x(), py = dir.y(), pz = dir.z();
    theta = std::acos(std::abs(pz)) * 180.0 / M_PI;
    phi   = std::atan2(py, px)      * 180.0 / M_PI;
}

// ReconstructFullEvent (Algorithm 1: weight = N/σ²)
void MySensitiveDetector::ReconstructFullEvent(
        const std::vector<double>& xTop,
        const std::vector<double>& xBottom,
        const std::vector<double>& yLeft,
        const std::vector<double>& yRight)
{
    std::lock_guard<std::mutex> lock(analysisMutex);
    auto* man = G4AnalysisManager::Instance();

    if (!fEventAction) {
        G4cout << "ERROR: fEventAction is null in ReconstructFullEvent!" << G4endl;
        return;
    }

    double PxT, NxT, sxT;
    double PxB, NxB, sxB;
    double PyL, NyL, syL;
    double PyR, NyR, syR;

    ComputeGaussianParams(xTop,    PxT, NxT, sxT);
    ComputeGaussianParams(xBottom, PxB, NxB, sxB);
    ComputeGaussianParams(yLeft,   PyL, NyL, syL);
    ComputeGaussianParams(yRight,  PyR, NyR, syR);

    double x_rec = WeightedMLE(PxT, NxT, sxT, PxB, NxB, sxB);
    double y_rec = WeightedMLE(PyL, NyL, syL, PyR, NyR, syR);

    G4ThreeVector truePos = fEventAction->GetInteractionPosition();
    double x_true = truePos.x();
    double y_true = truePos.y();
    double dx         = x_rec - x_true;
    double dy         = y_rec - y_true;
    double resolution = std::sqrt(dx*dx + dy*dy);

    double theta, phi;
    ComputeAngles(fEventAction->GetPrimaryDirection(), theta, phi);

    G4cout << "[Alg1] Rec=(" << x_rec/mm << "," << y_rec/mm
           << ") True=(" << x_true/mm << "," << y_true/mm
           << ") dx=" << dx/mm << " dy=" << dy/mm << " mm" << G4endl;

    man->FillH1(6,  x_rec / mm);
    man->FillH1(7,  y_rec / mm);
    man->FillH1(8,  dx / mm);
    man->FillH1(9,  dy / mm);
    man->FillH1(13, resolution / mm);

    man->FillH2(0, x_rec / mm, y_rec / mm);
    man->FillH2(2, theta, dx / mm);
    man->FillH2(3, theta, dy / mm);
    man->FillH2(4, theta, resolution / mm);
    man->FillH2(5, phi,   dx / mm);
    man->FillH2(6, phi,   dy / mm);

    man->FillNtupleDColumn(5, 0, x_rec / mm);
    man->FillNtupleDColumn(5, 1, y_rec / mm);
    man->FillNtupleDColumn(5, 2, x_true / mm);
    man->FillNtupleDColumn(5, 3, y_true / mm);
    man->FillNtupleDColumn(5, 4, dx / mm);
    man->FillNtupleDColumn(5, 5, dy / mm);
    man->FillNtupleDColumn(5, 6, theta);
    man->FillNtupleDColumn(5, 7, phi);
    man->FillNtupleDColumn(5, 8, resolution / mm);
    man->AddNtupleRow(5);
}

// ReconstructFullEventPaper (PAPER Eq.1: weight = N/σ)
void MySensitiveDetector::ReconstructFullEventPaper(
        const std::vector<double>& xTop,
        const std::vector<double>& xBottom,
        const std::vector<double>& yLeft,
        const std::vector<double>& yRight)
{
    std::lock_guard<std::mutex> lock(analysisMutex);
    auto* man = G4AnalysisManager::Instance();

    if (!fEventAction) {
        G4cout << "ERROR: fEventAction is null in ReconstructFullEventPaper!" << G4endl;
        return;
    }

    double PxT, NxT, sxT;
    double PxB, NxB, sxB;
    double PyL, NyL, syL;
    double PyR, NyR, syR;

    ComputeGaussianParams(xTop,    PxT, NxT, sxT);
    ComputeGaussianParams(xBottom, PxB, NxB, sxB);
    ComputeGaussianParams(yLeft,   PyL, NyL, syL);
    ComputeGaussianParams(yRight,  PyR, NyR, syR);

    double x_rec = PaperEq1(PxT, NxT, sxT, PxB, NxB, sxB);
    double y_rec = PaperEq1(PyL, NyL, syL, PyR, NyR, syR);

    G4ThreeVector truePos = fEventAction->GetInteractionPosition();
    double x_true = truePos.x();
    double y_true = truePos.y();
    double dx         = x_rec - x_true;
    double dy         = y_rec - y_true;
    double resolution = std::sqrt(dx*dx + dy*dy);

    double theta, phi;
    ComputeAngles(fEventAction->GetPrimaryDirection(), theta, phi);

    G4cout << "[PAPER] Rec=(" << x_rec/mm << "," << y_rec/mm
           << ") True=(" << x_true/mm << "," << y_true/mm
           << ") dx=" << dx/mm << " dy=" << dy/mm << " mm" << G4endl;

    man->FillH1(16, x_rec / mm);
    man->FillH1(17, y_rec / mm);
    man->FillH1(18, dx / mm);
    man->FillH1(19, dy / mm);
    man->FillH1(20, resolution / mm);

    man->FillH2(7,  x_rec / mm, y_rec / mm);
    man->FillH2(8,  theta, dx / mm);
    man->FillH2(9,  theta, dy / mm);
    man->FillH2(10, theta, resolution / mm);
    man->FillH2(11, phi,   dx / mm);
    man->FillH2(12, phi,   dy / mm);

    man->FillNtupleDColumn(6, 0, x_rec / mm);
    man->FillNtupleDColumn(6, 1, y_rec / mm);
    man->FillNtupleDColumn(6, 2, x_true / mm);
    man->FillNtupleDColumn(6, 3, y_true / mm);
    man->FillNtupleDColumn(6, 4, dx / mm);
    man->FillNtupleDColumn(6, 5, dy / mm);
    man->FillNtupleDColumn(6, 6, theta);
    man->FillNtupleDColumn(6, 7, phi);
    man->FillNtupleDColumn(6, 8, resolution / mm);
    man->AddNtupleRow(6);
}

G4bool MySensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    G4Track* track = step->GetTrack();

    if (track->GetDefinition() != G4OpticalPhoton::Definition())
        return false;

    auto* touchable = step->GetPreStepPoint()->GetTouchable();
    if (!touchable) return false;

    G4int copyNo = touchable->GetCopyNumber();
    G4ThreeVector pos = step->GetPostStepPoint()->GetPosition();

    if (copyNo < 0 || copyNo >= 4 * fNSensorsPerArray) {
        track->SetTrackStatus(fStopAndKill);
        return false;
    }

    auto* hit = new SensorHit();
    hit->SetSensorPosition(pos);
    hit->SetSensorEnergy(track->GetKineticEnergy());
    hit->SetCopyNumber(copyNo);
    SensorCollection->insert(hit);

    int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    int N   = fNSensorsPerArray;

    if      (copyNo < N)    RecordSensorData(0, pos.x(), pos.y(), evt, copyNo);
    else if (copyNo < 2*N)  RecordSensorData(1, pos.x(), pos.y(), evt, copyNo);
    else if (copyNo < 3*N)  RecordSensorData(2, pos.x(), pos.y(), evt, copyNo);
    else if (copyNo < 4*N)  RecordSensorData(3, pos.x(), pos.y(), evt, copyNo);

    track->SetTrackStatus(fStopAndKill);
    return true;
}

void MySensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
    if (!SensorCollection || SensorCollection->entries() == 0) {
        G4cout << "Event: no photon detections." << G4endl;
        return;
    }

    std::map<int,int> topC, bottomC, leftC, rightC;
    int N = fNSensorsPerArray;

    for (int i = 0; i < SensorCollection->entries(); ++i) {
        int c = (*SensorCollection)[i]->GetCopyNumber();
        if      (c < N)    topC   [c]++;
        else if (c < 2*N)  bottomC[c - N]++;
        else if (c < 3*N)  leftC  [c - 2*N]++;
        else if (c < 4*N)  rightC [c - 3*N]++;
    }

    auto buildPeVector = [&](const std::map<int,int>& counts,
                              std::vector<double>& peVec)
    {
        for (auto& [idx, nPhotons] : counts) {
            int npe = G4Poisson(nPhotons * fPDE);
            for (int i = 0; i < npe; ++i)
                peVec.push_back(static_cast<double>(idx));
        }
        int noise = G4Poisson(fNoiseLambda);
        for (int i = 0; i < noise; ++i)
            peVec.push_back(G4UniformRand() * N);
    };

    std::vector<double> topPE, botPE, leftPE, rightPE;
    buildPeVector(topC,    topPE);
    buildPeVector(bottomC, botPE);
    buildPeVector(leftC,   leftPE);
    buildPeVector(rightC,  rightPE);

    int nArrays = (!topPE.empty()) + (!botPE.empty()) +
                  (!leftPE.empty()) + (!rightPE.empty());

    G4cout << "=== EndOfEvent: arrays with signal = " << nArrays << "/4 ===" << G4endl;

    if (nArrays == 4) {
        ReconstructFullEvent     (topPE, botPE, leftPE, rightPE);
        ReconstructFullEventPaper(topPE, botPE, leftPE, rightPE);
    } else {
        G4cout << "Border/partial event: skipped (need all 4 arrays)." << G4endl;
    }
}
