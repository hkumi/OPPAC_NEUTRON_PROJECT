#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include "G4UserEventAction.hh"
#include "globals.hh"

namespace B4
{

// Forward declaration
class SteppingAction;

class EventAction : public G4UserEventAction
{
public:
    EventAction();
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);
    
    void SetSteppingAction(SteppingAction* steppingAction);

private:
    SteppingAction* fSteppingAction;
};

}
#endif
