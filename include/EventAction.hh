// ********************************************************************
// EventAction.hh
// ********************************************************************

#ifndef B4aEventAction_h
#define B4aEventAction_h 1

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

namespace B4a
{

class EventAction : public G4UserEventAction
{
  public:
    EventAction() = default;
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

    // NEW: Methods to store and retrieve particle information
    void SetInteractionPosition(const G4ThreeVector& pos) { fInteractionPos = pos; fPosSet = true; }
    void SetPrimaryDirection(const G4ThreeVector& dir) { fPrimaryDirection = dir; }
    
    G4ThreeVector GetInteractionPosition() const { return fInteractionPos; }
    G4ThreeVector GetPrimaryDirection() const { return fPrimaryDirection; }
    G4bool IsPositionSet() const { return fPosSet; }

  private:
    G4double fEnergySiPM = 0.;
    G4double fEnergySiPM2 = 0.;
    
    // NEW: Store particle tracking information
    G4ThreeVector fInteractionPos;
    G4ThreeVector fPrimaryDirection;
    G4bool fPosSet = false;
};

}

#endif
