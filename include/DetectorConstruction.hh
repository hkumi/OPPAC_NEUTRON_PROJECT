#ifndef B4DetectorConstruction_h
#define B4DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4GlobalMagFieldMessenger;

namespace B4
{

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction() = default;
    ~DetectorConstruction() override = default;

  public:
    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    // get methods
    //
    const G4VPhysicalVolume* GetSiPMPV() const;
    const G4VPhysicalVolume* GetSiPM2PV() const;
	G4double GetPitch() const { return pitch; }
	G4double GetSize() const { return size; }


  private:
    // methods
    //
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes(G4double, G4double);

    // data members
    //
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;
                                      // magnetic field messenger

    G4VPhysicalVolume* SiPMPV = nullptr; // the absorber physical volume
    G4VPhysicalVolume* SiPM2PV = nullptr;
    G4double pitch;
	G4double size;

    G4bool fCheckOverlaps = true; // option to activate checking of volumes overlaps
};

// inline functions

inline const G4VPhysicalVolume* DetectorConstruction::GetSiPMPV() const {
  return SiPMPV;
}

inline const G4VPhysicalVolume* DetectorConstruction::GetSiPM2PV() const {
  return SiPM2PV;
}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

