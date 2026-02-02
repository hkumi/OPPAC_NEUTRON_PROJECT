#ifndef ACTIONINITIALIZATION_HH
#define ACTIONINITIALIZATION_HH

#include "G4VUserActionInitialization.hh"

namespace B4
{
  class DetectorConstruction;
}

namespace B4a
{

class ActionInitialization : public G4VUserActionInitialization
{
public:
  ActionInitialization(B4::DetectorConstruction* detConstruction);
  ~ActionInitialization() override = default;

  void BuildForMaster() const override;
  void Build() const override;

private:
  B4::DetectorConstruction* fDetConstruction;
};

}

#endif
