// ********************************************************************
// EventAction.cc
// ********************************************************************

#include "EventAction.hh"
#include "RunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include <iomanip>

namespace B4a
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{
  // Initialisation per event
  fEnergySiPM = 0.;
  fEnergySiPM2 = 0.;
  
  //  Reset particle tracking flags
  fPosSet = false;
  fInteractionPos = G4ThreeVector(0, 0, 0);
  fPrimaryDirection = G4ThreeVector(0, 0, 0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  G4int nPrimaries = event->GetNumberOfPrimaryVertex();
  for (G4int iVertex = 0; iVertex < nPrimaries; ++iVertex) {
      G4PrimaryVertex* vertex = event->GetPrimaryVertex(iVertex);
      G4PrimaryParticle* primary = vertex->GetPrimary();
      
      if (primary->GetPDGcode() == 2112) { // neutron=2112; gamma=22
          auto analysisManager = G4AnalysisManager::Instance();
          auto eventID = event->GetEventID();
          auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
      }
  }
  
  //  Debug output at end of event
  if (fPosSet) {
      G4cout << "EventAction: Captured position: (" 
             << fInteractionPos.x()/CLHEP::mm << ", " 
             << fInteractionPos.y()/CLHEP::mm << ", " 
             << fInteractionPos.z()/CLHEP::mm << ") mm" << G4endl;
      G4cout << "EventAction: Direction: (" 
             << fPrimaryDirection.x() << ", " 
             << fPrimaryDirection.y() << ", " 
             << fPrimaryDirection.z() << ")" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

}
