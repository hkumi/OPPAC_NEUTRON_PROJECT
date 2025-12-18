#include "EventAction.hh"
#include "SteppingAction.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"

namespace B4  // <-- ADD THIS NAMESPACE
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
    : fSteppingAction(nullptr)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::SetSteppingAction(SteppingAction* steppingAction)
{
    fSteppingAction = steppingAction;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{
    // Optional initialization
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
    if (!fSteppingAction) return;
    
    // Tell stepping action to process and store the reconstructed position
    G4int eventID = event->GetEventID();
    fSteppingAction->FillReconstructedPosition(eventID);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // <-- END NAMESPACE B4
