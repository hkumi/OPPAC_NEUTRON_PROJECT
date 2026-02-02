#include "SensorHit.hh"
#include "G4UnitsTable.hh"
#include <iomanip>

// Initialize thread-local allocator
G4ThreadLocal G4Allocator<SensorHit>* SensorHitAllocator = nullptr;

SensorHit::SensorHit()
    : positionSensor(G4ThreeVector(0., 0., 0.)),
      energySensor(0.),
      SensorNumberID(-1),
      VolumeNameID(""),
      copyNumber(-1)
{ }

SensorHit::~SensorHit() { }

SensorHit::SensorHit(const SensorHit& right)
    : positionSensor(right.positionSensor),
      energySensor(right.energySensor),
      SensorNumberID(right.SensorNumberID),
      VolumeNameID(right.VolumeNameID),
      copyNumber(right.copyNumber)
{ }

const SensorHit& SensorHit::operator=(const SensorHit& right) {
    positionSensor = right.positionSensor;
    energySensor   = right.energySensor;
    SensorNumberID = right.SensorNumberID;
    VolumeNameID   = right.VolumeNameID;
    copyNumber     = right.copyNumber;
    return *this;
}

int SensorHit::operator==(const SensorHit& /*right*/) const {
    return 0;
}

// Thread-safe operator new
void* SensorHit::operator new(size_t) {
    if (!SensorHitAllocator) {
        SensorHitAllocator = new G4Allocator<SensorHit>;
    }
    return (void*)SensorHitAllocator->MallocSingle();
}

// Thread-safe operator delete
void SensorHit::operator delete(void* aSensorHit) {
    if (!SensorHitAllocator) {
        SensorHitAllocator = new G4Allocator<SensorHit>;
    }
    SensorHitAllocator->FreeSingle((SensorHit*)aSensorHit);
}
