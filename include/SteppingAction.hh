#ifndef STEPPINGACTION_HH
#define STEPPINGACTION_HH

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <vector>

namespace B4
{

// Forward declarations within the same namespace
class DetectorConstruction;
class EventAction;

class SteppingAction : public G4UserSteppingAction
{
public:
    SteppingAction(const DetectorConstruction* detConstruction,
                   EventAction* eventAction);
    virtual ~SteppingAction();

    virtual void UserSteppingAction(const G4Step* step);
    
    void FillReconstructedPosition(G4int eventID);
    
    // Arrays for storing positions
    std::vector<G4double> fXBottom;
    std::vector<G4double> fXTop;
    std::vector<G4double> fYLeft;
    std::vector<G4double> fYRight;

private:
    const DetectorConstruction* fDetConstruction;
    EventAction* fEventAction;
    
    G4int fCurrentEventID;
    
    void SaveEventData(G4int eventID);
};

}
#endif
