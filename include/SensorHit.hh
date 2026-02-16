#ifndef SENSORHIT_HH
#define SENSORHIT_HH 1

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"  

class SensorHit : public G4VHit {
public:
    SensorHit();
    ~SensorHit();
    SensorHit(const SensorHit&);

    const SensorHit& operator=(const SensorHit&);
    int operator==(const SensorHit&) const;

    // operator new/delete for G4
    static void* operator new(size_t);
    static void operator delete(void*);

    // Setters and Getters
    inline void SetSensorPosition(G4ThreeVector position) { positionSensor = position; };
    inline G4ThreeVector GetSensorPosition() const { return positionSensor; };

    inline void SetSensorNumber(G4int SensorNumber) { SensorNumberID = SensorNumber; };
    inline G4int GetSensorNumber() const { return SensorNumberID; };

    inline void SetVolumeName(G4String VolumeName) { VolumeNameID = VolumeName; };
    inline G4String GetVolumeName() const { return VolumeNameID; };

    inline void SetSensorEnergy(G4double energy) { energySensor = energy; };
    inline G4double GetSensorEnergy() const { return energySensor; };

    // Copy number methods
    inline void SetCopyNumber(G4int copyNo) { copyNumber = copyNo; }
    inline G4int GetCopyNumber() const { return copyNumber; }

private:
    G4ThreeVector positionSensor;
    G4double energySensor;
    G4int SensorNumberID;
    G4String VolumeNameID;
    G4int copyNumber;
};

// Vector collection of hits
typedef G4THitsCollection<SensorHit> SensorHitsCollection;

// Thread-local allocator
extern G4ThreadLocal G4Allocator<SensorHit>* SensorHitAllocator;

#endif
