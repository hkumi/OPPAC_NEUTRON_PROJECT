#ifndef SensorData_h
#define SensorData_h 1

#include <vector>
#include <map>

struct SensorHit {
    int sensorID;
    double posX;
    double posY; 
    double energy;
    double time;
    int photonCount;
};

struct EventData {
    int eventID;
    double trueX;
    double trueY;
    std::vector<SensorHit> hits;
    int totalPhotons;
    int sensorsHit;
};

class PositionReconstructor {
public:
    PositionReconstructor();
    ~PositionReconstructor();
    
    void LoadData(const std::string& filename);
    void ReconstructPositions();
    void CalculateResolution();
    void GeneratePlots();
    
private:
    std::vector<EventData> events;
    std::map<int, std::pair<double, double>> sensorPositions; // sensorID -> (x,y)
    
    void InitializeSensorMap();
    std::pair<double, double> CalculateCenterOfGravity(const EventData& event);
    std::pair<double, double> CalculateWeightedPosition(const EventData& event);
    double CalculateFWHM(const std::vector<double>& residuals);
};

#endif
