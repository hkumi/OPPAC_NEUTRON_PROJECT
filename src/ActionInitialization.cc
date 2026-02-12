// ActionInitialization.cc
#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "DetectorConstruction.hh"

namespace B4a
{

ActionInitialization::ActionInitialization(B4::DetectorConstruction* detConstruction)
 : fDetConstruction(detConstruction)
{}

void ActionInitialization::BuildForMaster() const
{
  SetUserAction(new B4::RunAction);
}

void ActionInitialization::Build() const
{
  SetUserAction(new B4::PrimaryGeneratorAction);
  SetUserAction(new B4::RunAction);
  
  auto eventAction = new EventAction;
  SetUserAction(eventAction);
  
  SetUserAction(new SteppingAction(fDetConstruction, eventAction));
}

}
