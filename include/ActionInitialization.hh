#ifndef ACTIONINITIALIZATION_HH
#define ACTIONINITIALIZATION_HH

#include "G4VUserActionInitialization.hh"

namespace B4
{

class DetectorConstruction;

class ActionInitialization : public G4VUserActionInitialization
{
public:
    ActionInitialization(DetectorConstruction* detConstruction);
    virtual ~ActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

private:
    DetectorConstruction* fDetConstruction;
};

}
#endif
