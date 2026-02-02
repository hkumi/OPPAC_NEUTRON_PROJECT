#include "detector.hh"

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

G4ThreadLocal G4ThreeVector truePos;
G4ThreadLocal bool truePosSet = false;

// --------------------------------------------------------------------
MySensitiveDetector::MySensitiveDetector(G4String name)
    : G4VSensitiveDetector(name),
      fPDE(0.25),
      fNoiseLambda(6.0),
      fMinSigma(0.5 * mm),
      fSensorPitch(1.0 * mm),
      fNSensorsPerArray(25)
{
    collectionName.insert("SensorCollection");
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
}

// --------------------------------------------------------------------
double MySensitiveDetector::IndexToPosition(double index)
{
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
void MySensitiveDetector::ReconstructFullEvent(const std::vector<double>& xTop,
                                       const std::vector<double>& xBottom,
                                       const std::vector<double>& yLeft,
                                       const std::vector<double>& yRight)
{
    std::lock_guard<std::mutex> lock(analysisMutex);
    auto* man = G4AnalysisManager::Instance();

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
            return (p1 * n1 / (s1*s1) + p2 * n2 / (s2*s2)) /
                   (n1 / (s1*s1)     + n2 / (s2*s2));
        else if (n1 > 0) return p1;
        else if (n2 > 0) return p2;
        return 0.0;
    };

    double x_rec = WeightedMLE(mxT, nT, sxT, mxB, nB, sxB);
    double y_rec = WeightedMLE(myL, nL, syL, myR, nR, syR);

    man->FillH1(6, x_rec);
    man->FillH1(7, y_rec);
    man->FillH2(0, x_rec, y_rec);

    double dx = x_rec - truePos.x();
    double dy = y_rec - truePos.y();
    double resolution = std::sqrt(dx*dx + dy*dy);

    man->FillH1(8, dx / mm);
    man->FillH1(9, dy / mm);
    man->FillH1(13, resolution / mm);  // Changed from 14 to 13
}

// --------------------------------------------------------------------
void MySensitiveDetector::ReconstructBorderEvent(
        const std::vector<double>& xTop,
        const std::vector<double>& xBottom,
        const std::vector<double>& yLeft,
        const std::vector<double>& yRight)
{
    std::lock_guard<std::mutex> lock(analysisMutex);
    auto* man = G4AnalysisManager::Instance();

    auto SimpleMean = [&](const std::vector<double>& v) {
        if (v.empty()) return 0.0;
        double sum = 0.0;
        for (double i : v) sum += IndexToPosition(i);
        return sum / v.size();
    };

    double x_rec = 0.0;
    if (!xTop.empty() && !xBottom.empty()) {
        x_rec = 0.5 * (SimpleMean(xTop) + SimpleMean(xBottom));
    } else if (!xTop.empty()) {
        x_rec = SimpleMean(xTop);
    } else if (!xBottom.empty()) {
        x_rec = SimpleMean(xBottom);
    }

    double y_rec = 0.0;
    if (!yLeft.empty() && !yRight.empty()) {
        y_rec = 0.5 * (SimpleMean(yLeft) + SimpleMean(yRight));
    } else if (!yLeft.empty()) {
        y_rec = SimpleMean(yLeft);
    } else if (!yRight.empty()) {
        y_rec = SimpleMean(yRight);
    }

    man->FillH2(1, x_rec, y_rec);

    double dx = x_rec - truePos.x();
    double dy = y_rec - truePos.y();
    man->FillH1(14, dx / mm);  // Border_Res_X
    man->FillH1(15, dy / mm);  // Border_Res_Y
}

// --------------------------------------------------------------------
G4bool MySensitiveDetector::ProcessHits(G4Step* step,
                                       G4TouchableHistory*)
{
    G4Track* track = step->GetTrack();

    if (!truePosSet &&
        track->GetDefinition()->GetPDGCharge() != 0 &&
        track->GetParentID() == 0)
    {
        truePos = step->GetPreStepPoint()->GetPosition();
        truePosSet = true;
    }

    if (track->GetDefinition() != G4OpticalPhoton::Definition())
        return false;

    auto* touchable = step->GetPreStepPoint()->GetTouchable();
    if (!touchable) return false;

    G4int copyNo = touchable->GetCopyNumber();
    G4ThreeVector pos = step->GetPostStepPoint()->GetPosition();

    auto* hit = new SensorHit();
    hit->SetSensorPosition(pos);
    hit->SetSensorEnergy(track->GetKineticEnergy());
    hit->SetCopyNumber(copyNo);
    SensorCollection->insert(hit);

    int evt = G4RunManager::GetRunManager()
                  ->GetCurrentEvent()->GetEventID();

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
    truePosSet = false;
    if (!SensorCollection || SensorCollection->entries() == 0) return;

    std::map<int,int> topC, bottomC, leftC, rightC;
    int N = fNSensorsPerArray;

    for (int i = 0; i < SensorCollection->entries(); ++i) {
        int c = (*SensorCollection)[i]->GetCopyNumber();
        if (c < N)        topC[c]++;
        else if (c < 2*N) bottomC[c - N]++;
        else if (c < 3*N) leftC[c - 2*N]++;
        else if (c < 4*N) rightC[c - 3*N]++;
    }

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

    int arrays = (!topI.empty()) + (!botI.empty()) +
                 (!leftI.empty()) + (!rightI.empty());

    if (arrays >= 2) {
        if (arrays < 4)
            ReconstructBorderEvent(topI, botI, leftI, rightI);
        else
            ReconstructFullEvent(topI, botI, leftI, rightI);
    }

    //SensorCollection->clear();
}
