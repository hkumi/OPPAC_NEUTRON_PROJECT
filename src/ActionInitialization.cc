#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "DetectorConstruction.hh"

namespace B4
{

ActionInitialization::ActionInitialization(DetectorConstruction* detConstruction)
    : G4VUserActionInitialization(),
      fDetConstruction(detConstruction)
{}

ActionInitialization::~ActionInitialization()
{}

void ActionInitialization::BuildForMaster() const
{
    SetUserAction(new RunAction());
}

void ActionInitialization::Build() const
{
    SetUserAction(new PrimaryGeneratorAction());
    SetUserAction(new RunAction());
    
    // Create event action
    auto eventAction = new EventAction();
    SetUserAction(eventAction);
    
    // Create stepping action
    auto steppingAction = new SteppingAction(fDetConstruction, eventAction);
    SetUserAction(steppingAction);
    
    // Link them together
    eventAction->SetSteppingAction(steppingAction);
}

}
