#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include <mutex>
#include <vector>
#include "SensorHit.hh"

namespace B4a {
    class EventAction;  
}

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class MySensitiveDetector : public G4VSensitiveDetector
{
public:
    MySensitiveDetector(G4String name);
    ~MySensitiveDetector();
    
    void SetEventAction(B4a::EventAction* evtAction) { fEventAction = evtAction; } 

    void Initialize(G4HCofThisEvent* HCE) override;
    G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
    void EndOfEvent(G4HCofThisEvent* HCE) override;

    void SetPDE(G4double pde) { fPDE = pde; }
    void SetNoiseLambda(G4double lambda) { fNoiseLambda = lambda; }
    void SetNSensorsPerArray(G4int n) { fNSensorsPerArray = n; }

    void SetSensorPitch(G4double pitch)
    {
        fSensorPitch = pitch;
        fMinSigma = pitch / 2.0;
    }

private:
    B4a::EventAction* fEventAction = nullptr; 
    
    void RecordSensorData(int ntupleIndex,
                          double x, double y,
                          int event, int copyNo);

    void ReconstructFullEvent(const std::vector<double>& xTop,
                              const std::vector<double>& xBottom,
                              const std::vector<double>& yLeft,
                              const std::vector<double>& yRight);

    void ReconstructBorderEvent(const std::vector<double>& xTop,
                                const std::vector<double>& xBottom,
                                const std::vector<double>& yLeft,
                                const std::vector<double>& yRight,
                                int nArrays); 

    double IndexToPosition(double index);

    SensorHitsCollection* SensorCollection;

    G4double fPDE;
    G4double fNoiseLambda;
    G4double fMinSigma;
    G4double fSensorPitch;
    G4int    fNSensorsPerArray;

    static std::mutex analysisMutex;
};

#endif
