// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B4/B4a/src/DetectorConstruction.cc
/// \brief Implementation of the B4::DetectorConstruction class
/// CORRECTED VERSION - Optimized for OPPAC detector

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVDivision.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SubtractionSolid.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"

#include "G4SDManager.hh"
#include "detector.hh" 

namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  
    pitch = 0.5 * cm;
    size = 5.0 * cm;
    
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes(pitch, size);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{

    G4int ncomponents, natoms;
    G4double massfraction;

    G4double Vdens = 1.e-25 * g / cm3;
    G4double Vpres = 1.e-19 * pascal;
    G4double Vtemp = 0.1 * kelvin;

    G4double a, z, dens;

    // G4 materials 
    auto nistManager = G4NistManager::Instance();
    nistManager->FindOrBuildElement("B");
    nistManager->FindOrBuildElement("Mn");
    nistManager->FindOrBuildElement("Cr");
    nistManager->FindOrBuildElement("Ni");
    nistManager->FindOrBuildElement("Al");
    nistManager->FindOrBuildElement("K");
    nistManager->FindOrBuildElement("Mg");
    nistManager->FindOrBuildElement("Na");
    nistManager->FindOrBuildElement("Ca");
    nistManager->FindOrBuildElement("Cs");
    nistManager->FindOrBuildElement("I");
    nistManager->FindOrBuildElement("H");
    nistManager->FindOrBuildElement("C");
    nistManager->FindOrBuildElement("Cl");
    nistManager->FindOrBuildElement("O");
    nistManager->FindOrBuildElement("Zn");
    nistManager->FindOrBuildElement("S");
    nistManager->FindOrBuildElement("N");
    nistManager->FindOrBuildElement("Ag");
    nistManager->FindOrBuildElement("Si");

    nistManager->FindOrBuildElement("Ar");
    nistManager->FindOrBuildElement("C");
    nistManager->FindOrBuildElement("F");

    G4Material* air = nistManager->FindOrBuildMaterial("G4_AIR");
    G4Material* concrete = nistManager->FindOrBuildMaterial("G4_CONCRETE");
    G4Material* lead = nistManager->FindOrBuildMaterial("G4_Pb");
    G4Material* tungsten = nistManager->FindOrBuildMaterial("G4_W");
    G4Material* SiO2 = nistManager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    
    // vacuum
    G4Material* vacuum = new G4Material("vacuum", z = 1, a = 1.01 * g / mole, Vdens, kStateGas, Vtemp, Vpres);

    // CF4
    auto C = G4Element::GetElement("C");
    auto F = G4Element::GetElement("F");
    G4Material* CF4 = new G4Material("CF4", 3.72 * mg / cm3, ncomponents = 2, kStateGas, 293.15 * kelvin, 1.0 * atmosphere);
    CF4->AddElement(C, natoms = 1);
    CF4->AddElement(F, natoms = 4);

    // Ar (gas)
    auto Ar = G4Element::GetElement("Ar");
    G4Material* Ar_gas = new G4Material("Ar_gas", 1.782 * mg / cm3, ncomponents = 1, kStateGas, 293.15 * kelvin, 1.0 * atmosphere);
    Ar_gas->AddElement(Ar, natoms = 1);

    // Ar:CF4 (90:10)
    G4Material* Ar_CF4 = new G4Material("Ar_CF4", 0.061 * mg / cm3, ncomponents = 2, kStateGas, 293.15 * kelvin, 30.e-3 * bar); //2.06mg/cm3 @1atm ; 0.061mg/cm3 @30mbar
    Ar_CF4->AddMaterial(Ar_gas, massfraction = 0.9);
    Ar_CF4->AddMaterial(CF4, massfraction = 0.1);
    
    // WSL (Wavelength Shifting Fiber)
        // Polystyrene core
    auto H = G4Element::GetElement("H");
    G4Material* polystyrene = new G4Material("polystyrene", 1.05 * g / cm3, 2);
    polystyrene->AddElement(C, natoms = 8); // 8 C
    polystyrene->AddElement(H, natoms = 8); // 8 H
        // PMMA cladding
    auto O = G4Element::GetElement("O");
    G4Material* PMMA = new G4Material("PMMA", 1.20 * g / cm3, 3);
    PMMA->AddElement(C, natoms = 5);
    PMMA->AddElement(H, natoms = 8);
    PMMA->AddElement(O, natoms = 2);

    // mylar
    G4Material* mylar = new G4Material("mylar", 1.39 * g / cm3, 3);
    mylar->AddElement(C, natoms = 5);
    mylar->AddElement(H, natoms = 4);
    mylar->AddElement(O, natoms = 2);

    // sensor (SiPM)
    auto Si = G4Element::GetElement("Si");
    G4Material* sensor = new G4Material("SiPM", 2.33 * g / cm3, 2);
    sensor->AddElement(Si, natoms = 1); // 1 Si
    sensor->AddElement(O, natoms = 2); // 2 O

    // Aluminium
    auto Al = G4Element::GetElement("Al");
    G4Material* alum = new G4Material("alum", 2.7 * g / cm3, 1);
    alum->AddElement(Al, natoms = 1);

    // HDPE
    G4Element* H_hp = new G4Element("TS_H_of_Polyethylene", "H", 1, 1.0079 * g / mole);
    G4Material* HDPE = new G4Material("HDPE", 0.93 * g / cm3, 2);
    HDPE->AddElement(C, 2); // 2 C
    HDPE->AddElement(H_hp, 4); // 4 H

    // Teflon - collimator material (BEST CHOICE for optical reflectivity)
    G4Material* Teflon = new G4Material("Teflon", 2.2 * g / cm3, 2, kStateSolid);
    Teflon->AddElement(C, 0.240183);
    Teflon->AddElement(F, 0.759817);

    // Teflon optical properties (highly reflective)
    const G4int TNbEntries = 4;
    G4double pdTeflonPhotonMomentum[TNbEntries] = {
        2.07 * eV, 2.34 * eV, 2.64 * eV, 3.10 * eV
    };
    G4double pdTeflonRefractiveIndex[TNbEntries] = {1.34, 1.34, 1.34, 1.34};
    G4double pdTeflonReflectivity[TNbEntries] = {0.95, 0.95, 0.95, 0.95};
    G4double pdTeflonAbsLength[TNbEntries] = {
        1.0 * m, 1.0 * m, 1.0 * m, 1.0 * m
    };

    G4MaterialPropertiesTable* pTeflonPropertiesTable = new G4MaterialPropertiesTable();
    pTeflonPropertiesTable->AddProperty("RINDEX", pdTeflonPhotonMomentum, 
                                        pdTeflonRefractiveIndex, TNbEntries);
    pTeflonPropertiesTable->AddProperty("ABSLENGTH", pdTeflonPhotonMomentum, 
                                        pdTeflonAbsLength, TNbEntries);
    pTeflonPropertiesTable->AddProperty("REFLECTIVITY", pdTeflonPhotonMomentum, 
                                        pdTeflonReflectivity, TNbEntries);
    Teflon->SetMaterialPropertiesTable(pTeflonPropertiesTable);

    // Optical properties // --------------------------------------------------------

    // Air & vacuum
    auto air_mpt = new G4MaterialPropertiesTable();
    air_mpt->AddProperty("RINDEX", "Air");
    air->SetMaterialPropertiesTable(air_mpt);
    vacuum->SetMaterialPropertiesTable(air_mpt);

    // Scintillating gas -  scintillation yield
    const G4int nScint = 18;
    G4double gasEnergy[nScint] = {
      1.55 * eV, 1.61 * eV, 1.70 * eV, 1.78 * eV, 1.88 * eV,
      2.00 * eV, 2.11 * eV, 2.22 * eV, 2.54 * eV, 2.87 * eV,
      3.42 * eV, 3.89 * eV, 4.35 * eV, 4.72 * eV, 4.88 * eV, 
      5.40 * eV, 5.83 * eV, 6.16 * eV
    };
    G4double gasScintSp[nScint] = {
      0.01, 0.01, 0.04, 0.10, 0.20,
      0.26, 0.18, 0.04, 0.00, 0.02,
      0.04, 0.13, 0.30, 0.19, 0.25,
      0.03, 0.00, 0.00
    };
    G4double gasAbsLength[nScint] = {
      4 * m, 4 * m, 4 * m, 4 * m, 4 * m,
      4 * m, 4 * m, 4 * m, 4 * m, 4 * m,
      4 * m, 4 * m, 4 * m, 4 * m, 4 * m,
      4 * m, 4 * m, 4 * m
    };
     auto ArCF4_mpt = new G4MaterialPropertiesTable();
    ArCF4_mpt->AddProperty("SCINTILLATIONCOMPONENT1", gasEnergy, gasScintSp, nScint);
    ArCF4_mpt->AddProperty("RINDEX", "Air");
    // CRITICAL: Use realistic scintillation yield
    // M.Cortesi et al. (2018): 2500/keV = 2500000/MeV for CF4
    ArCF4_mpt->AddConstProperty("SCINTILLATIONYIELD",2500000 / MeV); 
    ArCF4_mpt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
    ArCF4_mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 15 * ns);
    ArCF4_mpt->AddProperty("ABSLENGTH", gasEnergy, gasAbsLength, nScint);
    ArCF4_mpt->AddConstProperty("RESOLUTIONSCALE", 1.0);
    Ar_CF4->SetMaterialPropertiesTable(ArCF4_mpt);
    Ar_gas->SetMaterialPropertiesTable(ArCF4_mpt);

    // SiPM sensors
    const G4int nSiPM = 6;
    G4double SiPMEnergy[nSiPM] = {
      2.07 * eV, // 600 nm
      2.34 * eV, // 530 nm
      2.64 * eV, // 470 nm
      2.75 * eV, // 450 nm
      2.95 * eV, // 420 nm
      3.10 * eV  // 400 nm
    };
    G4double SiPMQE[nSiPM] = {
      0.25,  // 600 nm
      0.35,  // 530 nm
      0.45,  // 470 nm
      0.48,  // 450 nm
      0.50,  // 420 nm (máximo)
      0.45   // 400 nm
    };
    G4double SiPMrIndex[nSiPM] = {
      3.4, 3.4, 3.4, 3.4, 3.4, 3.4
    };
    auto SiPM_mpt = new G4MaterialPropertiesTable();
    SiPM_mpt->AddProperty("RINDEX", SiPMEnergy, SiPMrIndex, nSiPM);
    SiPM_mpt->AddProperty("EFFICIENCY", SiPMEnergy, SiPMQE, nSiPM);
    sensor->SetMaterialPropertiesTable(SiPM_mpt);

    // mylar & Al
    const G4int nAl = 6;
    G4double AlEnergy[nAl] = {
      2.07 * eV, // 600 nm
      2.34 * eV, // 530 nm
      2.64 * eV, // 470 nm
      2.75 * eV, // 450 nm
      2.95 * eV, // 420 nm
      3.10 * eV  // 400 nm
    };
    G4double AlrIndex[nAl] = {
      1.373, 1.373, 1.373, 1.373, 1.373, 1.373
    };
    G4double AlAbsLength[nAl] = {
        4 * cm, 4 * cm, 4 * cm, 4 * cm, 4 * cm, 4 * cm
    };
    auto Al_mtp = new G4MaterialPropertiesTable();
    Al_mtp->AddProperty("RINDEX", AlEnergy, AlrIndex, nAl);
    Al_mtp->AddProperty("ABSLENGTH", AlEnergy, AlAbsLength, nAl);
    alum->SetMaterialPropertiesTable(Al_mtp);
    mylar->SetMaterialPropertiesTable(Al_mtp);

    G4cout << *(G4Material::GetMaterialTable()) << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes(G4double pitch, G4double size)
{

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("vacuum");
  auto air = G4Material::GetMaterial("G4_AIR");
  auto atmosphere = G4Material::GetMaterial("Ar_CF4");
  auto mirror = G4Material::GetMaterial("alum");
  auto sensor = G4Material::GetMaterial("SiPM");
  auto PMMA = G4Material::GetMaterial("PMMA");
  auto Teflon = G4Material::GetMaterial("Teflon"); // Use Teflon for collimator

  auto mylar = G4Material::GetMaterial("mylar");
  auto HDPE = G4Material::GetMaterial("HDPE");
  auto polystyrene = G4Material::GetMaterial("polystyrene");

  G4RotationMatrix* rot = nullptr;


  // -------------------------- // 
  // world - ArCF4 (90/10)
  // -------------------------- //
  G4double fBoxSize = 10 * cm;
  G4double fBoxDepth = 1 * cm;

  G4Box* sBox = new G4Box("world",
      fBoxSize, fBoxSize, fBoxDepth/2);

  G4LogicalVolume* fLBox = new G4LogicalVolume(sBox,
      atmosphere, "World");

  G4VPhysicalVolume* fPBox = new G4PVPlacement(0,
      G4ThreeVector(), fLBox, "World", 0, false, 0, fCheckOverlaps);

 
  G4double sensPitch = pitch;
  G4double sensGap = 1.1 * mm;

  // OPTIMIZED: 1mm sensor size 
  G4double sensSize = 1.0 * mm; 
  G4double sipmSize = 1.0 * mm; // Match the sensor size to hole size
  
  // -------------------------- //
  // collimator cell + array
  // -------------------------- //
  G4double cellSize = sipmSize + 0.10 * mm; // Total cell size including walls
  
  // OPTIMIZED: Shorter collimator for better light collection efficiency
  // Trade-off: 10mm gives good angular resolution while maintaining efficiency
  G4double cellLength = 3.0 * mm; // OPTIMUM POINT
  
  G4double nCells = 25;
  G4double collSize = nCells * cellSize;

  G4Box* sBlock = new G4Box("cellBlock",
      cellLength / 2, cellSize / 2, cellSize / 2);

  G4Box* sHole = new G4Box("cellHole",
      cellLength / 2, sipmSize / 2, sipmSize / 2);

  G4SubtractionSolid* sCell = new G4SubtractionSolid("cell",
      sBlock, sHole, 0, G4ThreeVector(0, 0, 0));

  //  Use Teflon for highly reflective collimator
  G4LogicalVolume* fLCell = new G4LogicalVolume(sCell,
      Teflon, "Cell");

  G4int copyNo = 0;
  G4double xPos, yPos;
  
  // Store cell physical volumes for optical surface definition
  std::vector<G4VPhysicalVolume*> cellVolumes;
  
  for (G4int i = 0; i < 4; i++) {

      G4double angle = i * 90. * deg;
      auto cellRot = new G4RotationMatrix();
      cellRot->rotateZ(angle); 

      for (G4int j = 0; j < nCells; j++) {

          if (i == 0 || i == 2) {
            xPos = (collSize / 2 + cellLength / 2) * std::cos(angle);
            yPos = ((-1.0 * nCells / 2 + 1.0 * j + 1.0/2) * cellSize) + (collSize / 2 + cellLength / 2) * std::sin(angle);
          }
          else {
            xPos = ((-1.0 * nCells / 2 + 1.0 * j + 1.0/2) * cellSize) + (collSize / 2 + cellLength / 2) * std::cos(angle);
            yPos = (collSize / 2 + cellLength / 2) * std::sin(angle);
          }

          G4VPhysicalVolume* fPCell = new G4PVPlacement(cellRot,
              G4ThreeVector(xPos, yPos, 0), fLCell, "Cell", fLBox, false, copyNo, true);
          
          cellVolumes.push_back(fPCell);
          copyNo++;
      }
  }

  // Optical properties // --------------------------------------------------------
  auto sipmSurf = new G4OpticalSurface("sipmSurf");
  sipmSurf->SetType(dielectric_dielectric);
  sipmSurf->SetFinish(polished);
  sipmSurf->SetModel(unified);

  auto mylarSurf = new G4OpticalSurface("mylSurf");
  mylarSurf->SetType(dielectric_metal);
  mylarSurf->SetFinish(polished);
  mylarSurf->SetModel(unified);

  const G4int nOptic = 6;
  G4double photonEnergy[nOptic] = {
      2.07 * eV, 2.34 * eV, 2.64 * eV, 2.75 * eV, 2.95 * eV, 3.10 * eV
  };
  G4double mylarRefl[nOptic] = { 1, 1, 1, 1, 1, 1 };
  G4double mylarTrans[nOptic] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  G4double SiPMeff[nOptic] = { 0.25, 0.35, 0.45, 0.48, 0.50, 0.45 };

  auto sipm_mtp = new G4MaterialPropertiesTable();
  sipm_mtp->AddProperty("EFFICIENCY", photonEnergy, SiPMeff, nOptic);
  sipmSurf->SetMaterialPropertiesTable(sipm_mtp);

  auto mylarS_mtp = new G4MaterialPropertiesTable();
  mylarS_mtp->AddProperty("REFLECTIVITY", photonEnergy, mylarRefl, nOptic);
  mylarS_mtp->AddProperty("TRANSMITTANCE", photonEnergy, mylarTrans, nOptic);
  mylarSurf->SetMaterialPropertiesTable(mylarS_mtp);
  // ------------------------------------------------------------------------------


  // -------------------------- //
  // mylar foils - electrodes
  // -------------------------- //
  G4double mylThickness = 0.012 * mm;
  G4double collThickness = collSize;
  G4double mylPos = cellSize;

  G4Box* sMyl = new G4Box("myl",
      collSize/2, collSize/2, mylThickness/2);

  G4LogicalVolume* fLMyl = new G4LogicalVolume(sMyl,
      mirror, "Myl");

  G4VPhysicalVolume* fPMylA = new G4PVPlacement(0,
      G4ThreeVector(0, 0, mylPos/2), fLMyl, "MylA", fLBox, false, 0, true);

  G4VPhysicalVolume* fPMylK = new G4PVPlacement(0,
      G4ThreeVector(0, 0, -mylPos/2), fLMyl, "MylK", fLBox, false, 0, true);


  // -------------------------- //
  // (n,p) HDPE conversor
  // -------------------------- //
  G4double convThickness = 0.20 * mm;
  G4double convPos = mylPos / 2 + mylThickness / 2 + convThickness / 2;

  G4Box* sConv = new G4Box("conv",
      collSize / 2, collSize / 2, convThickness / 2);

  G4LogicalVolume* fLConv = new G4LogicalVolume(sConv,
      HDPE, "Conv");

  G4VPhysicalVolume* fPConv = new G4PVPlacement(0,
      G4ThreeVector(0, 0, convPos), fLConv, "Conv", fLBox, false, 0, true);


  // -------------------------- //
  // SiPMs - OPTIMIZED POSITIONING
  // -------------------------- //
  // STRATEGY: Place sensors OUTSIDE collimator but very close
  // This allows photons to exit the collimator channel and hit the sensor directly
  // with minimal air gap, maximizing collection efficiency
  
  G4double sipmWidth = 1.0 * mm; // Sensor thickness in beam direction
  G4double airGap = 0.1 * mm; // Minimal gap between collimator exit and sensor

  G4Box* sSiPM = new G4Box("sipm",
      sipmWidth/2, sipmSize/2, sipmSize/2);

  // Make sensors solid red and visible
  G4VisAttributes* sipmVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
  sipmVisAtt->SetVisibility(true);
  sipmVisAtt->SetForceSolid(true);

  G4LogicalVolume* fLSiPM = new G4LogicalVolume(sSiPM,
      sensor, "SiPM");
  fLSiPM->SetVisAttributes(sipmVisAtt);

  copyNo = 0;
  for (G4int i = 0; i < 4; i++) {

      G4double angle = i * 90. * deg;
      auto sipmRot = new G4RotationMatrix();
      sipmRot->rotateZ(angle);

      for (G4int j = 0; j < nCells; j++) {

          // Place sensors just outside collimator exit
          // Position = collimator_center + collimator_half_length  + sensor_half_width
          if (i == 0 || i == 2) {
              xPos = (collSize / 2 + cellLength  + sipmWidth/2) * std::cos(angle);
              yPos = ((-1.0 * nCells / 2 + 1.0 * j + 1.0 / 2) * cellSize) + 
                     (collSize / 2 + cellLength / 2) * std::sin(angle);
          }
          else {
              xPos = ((-1.0 * nCells / 2 + 1.0 * j + 1.0 / 2) * cellSize) + 
                     (collSize / 2 + cellLength / 2) * std::cos(angle);
              yPos = (collSize / 2 + cellLength + sipmWidth / 2) * std::sin(angle);
          }

          G4VPhysicalVolume* fPSiPM = new G4PVPlacement(sipmRot,
              G4ThreeVector(xPos, yPos, 0), fLSiPM, "SiPM", fLBox, false, copyNo, true);

          // Define optical border surface between SiPM and world
          //new G4LogicalBorderSurface("SiPM-Air", fPSiPM, fPBox, sipmSurf);
          
          copyNo++;
      }
  }

  // Define optical surfaces for mylar electrodes
  new G4LogicalBorderSurface("MylarA-Air_in", fPBox, fPMylA, mylarSurf);
  new G4LogicalBorderSurface("MylarA-Air_out", fPMylA, fPBox, mylarSurf);

  new G4LogicalBorderSurface("MylarK-Air_in", fPBox, fPMylK, mylarSurf);
  new G4LogicalBorderSurface("MylarK-Air_out", fPMylK, fPBox, mylarSurf);

  // Set visual attributes
  G4VisAttributes* cellVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.3)); // Transparent yellow
  cellVisAtt->SetVisibility(true);
  fLCell->SetVisAttributes(cellVisAtt);
  
  fLMyl->SetVisAttributes(G4VisAttributes(G4Colour(0.8, 0.8, 0.8, 0.9))); // Grey
  fLConv->SetVisAttributes(G4VisAttributes(G4Colour(0.0, 1.0, 0.0))); // Green
  
  fLBox->SetVisAttributes(G4VisAttributes::GetInvisible());
  
  return fPBox;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);

  // **********************************************
  // Register Sensitive Detectors for SiPMs
  // **********************************************
  
  auto SDmanager = G4SDManager::GetSDMpointer();
  SDmanager->SetVerboseLevel(1);
  
  // Create sensitive detector for SiPMs
  auto sipmSD = new MySensitiveDetector("SiPM_SD");
  
  // Register the sensitive detector
  SDmanager->AddNewDetector(sipmSD);
  
  // Find all SiPM logical volumes and assign the sensitive detector
  auto logicalVolumeStore = G4LogicalVolumeStore::GetInstance();
  
  for (auto logicalVolume : *logicalVolumeStore) {
    G4String volumeName = logicalVolume->GetName();
    
    // Assign SD to all SiPM logical volumes
    if (volumeName.find("SiPM") != std::string::npos) {
      logicalVolume->SetSensitiveDetector(sipmSD);
      G4cout << "Assigned sensitive detector to volume: " << volumeName << G4endl;
    }
  }
  
  G4cout << "Sensitive detectors registered successfully!" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
