//
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

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

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


    // Optical properties // --------------------------------------------------------

    // Air & vacuum
    auto air_mpt = new G4MaterialPropertiesTable();
    air_mpt->AddProperty("RINDEX", "Air");
    air->SetMaterialPropertiesTable(air_mpt);
    vacuum->SetMaterialPropertiesTable(air_mpt);
	//Ar_CF4->SetMaterialPropertiesTable(air_mpt);

    // Scintillating gas
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
    //G4double rIndexScint[nScint] = {
    //  2.356, 2.356, 2.356, 2.356, 2.356,
    //  2.356, 2.356, 2.356, 2.356, 2.356
    //};
    G4double gasAbsLength[nScint] = {
	  4 * m, 4 * m, 4 * m, 4 * m, 4 * m,
      4 * m, 4 * m, 4 * m, 4 * m, 4 * m,
      4 * m, 4 * m, 4 * m, 4 * m, 4 * m,
      4 * m, 4 * m, 4 * m
    };

    G4double gasScatLength[nScint];
     for (int i = 0; i < nScint; i++) {
         gasScatLength[i] = 5.0 * cm;  // Start with 5 cm scattering length
    }

    G4double gasRIndex[nScint];
     for (int i = 0; i < nScint; i++) {
     gasRIndex[i] = 1.0005;  // Refractive index for Ar:CF4 gas
    }
    auto ArCF4_mpt = new G4MaterialPropertiesTable();
    ArCF4_mpt->AddProperty("SCINTILLATIONCOMPONENT1", gasEnergy, gasScintSp, nScint);
    ArCF4_mpt->AddProperty("RAYLEIGH", gasEnergy, gasScatLength, nScint);
    //ArCF4_mpt->AddProperty("RINDEX", gasScintEnergy, rIndexScint, nScint);
    ArCF4_mpt->AddProperty("RINDEX", "Air");
    ArCF4_mpt->AddConstProperty("SCINTILLATIONYIELD", 5000 / MeV); // Based on nOPPAC notes for Ar:CF4 (90/10)
    ArCF4_mpt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
    ArCF4_mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 15 * ns);
    ArCF4_mpt->AddProperty("ABSLENGTH", gasEnergy, gasAbsLength, nScint);
    ArCF4_mpt->AddConstProperty("RESOLUTIONSCALE", 1.0);
    Ar_CF4->SetMaterialPropertiesTable(ArCF4_mpt);
	Ar_gas->SetMaterialPropertiesTable(ArCF4_mpt);

    // WSF (general)
    /*
    const G4int nWSF = 6;
    G4double photonEnergyWSF[nWSF] = {
      2.34 * eV, 2.43 * eV, 2.53 * eV, 2.70 * eV, 2.75 * eV, 2.82 * eV
    };
    G4double emitWSF[nWSF] = {
      0.8, 1.0, 0.9, 0.6, 0.3, 0.1
    };
    G4double absLengthWSF[nWSF] = {
      4.0 * m, 4.0 * m, 4.0 * m, 4.0 * m, 4.0 * m, 4.0 * m
    };
    G4double rIndexWSF[nWSF] = {
      1.60, 1.60, 1.60, 1.60, 1.60, 1.60
    };
    
    const G4int nWSF = 26;
    G4double photonEnergyWSF[nWSF] = {
        2.07 * eV, 2.10 * eV, 2.14 * eV, 2.17 * eV, 2.21 * eV, 2.26 * eV, 2.29 * eV, 2.34 * eV,
        2.38 * eV, 2.43 * eV, 2.48 * eV, 2.53 * eV, 2.58 * eV, 2.63 * eV, 2.70 * eV,
        2.77 * eV, 2.82 * eV, 2.88 * eV, 2.94 * eV, 3.02 * eV, 3.09 * eV, 3.18 * eV,
        3.26 * eV, 3.35 * eV, 3.44 * eV, 3.55 * eV 
    };
    G4double emitWSF[nWSF] = { 
        0.01, 0.01, 0.02, 0.05, 0.1, 0.16, 0.24,
        0.29, 0.35, 0.44, 0.67, 0.99, 0.88, 0.56,
        0.06, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
        0.00, 0.00, 0.00, 0.00, 0.00
    };
    G4double absLengthWSF[nWSF] = {
        4.00 * m, 4.00 * m, 4.00 * m, 4.00 * m, 4.00 * m, 4.00 * m, 4.00 * m,
        4.00 * m, 4.00 * m, 4.00 * m, 4.00 * m, 4.00 * m, 4.00 * m, 4.00 * m,
        4.00 * m, 1.00 * mm, 1.00 * mm, 1.00 * cm, 1.00 * mm, 1.00 * cm, 5.00 * cm,
        10.0 * cm, 1.00 * m, 1.00 * m, 4.00 * m, 4.00 * m
    };
    G4double rIndexWSF[nWSF] = {
      1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60,
      1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60,
      1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60,
      1.60, 1.60, 1.60, 1.60, 1.60
    };
    G4double wlsTimeConstant = 12. * ns;
    auto WSF_mpt = new G4MaterialPropertiesTable();
    WSF_mpt->AddProperty("WLSCOMPONENT", photonEnergyWSF, emitWSF, nWSF);
    WSF_mpt->AddProperty("WLSABSLENGTH", photonEnergyWSF, absLengthWSF, nWSF);
    WSF_mpt->AddConstProperty("WLSTIMECONSTANT", wlsTimeConstant);
    WSF_mpt->AddProperty("RINDEX", photonEnergyWSF, rIndexWSF, nWSF);
    polystyrene->SetMaterialPropertiesTable(WSF_mpt);
    */


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
      0.50,  // 420 nm (maximo)
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


    //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
      //------------------------------------ teflon -----------------------------------

  const G4int TNbEntries = 4;
  G4Material* Teflon = new G4Material("Teflon", 2.2*g/cm3, 2, kStateSolid);
  Teflon->AddElement(C, 0.240183);
  Teflon->AddElement(F, 0.759817);

  G4double pdTeflonPhotonMomentum[TNbEntries]  = {1.0*eV,6.91*eV, 6.98*eV, 7.05*eV};
  G4double pdTeflonRefractiveIndex[TNbEntries] = {1.34,1.34,1.34,1.34};
  G4double pdTeflonReflectivity[TNbEntries]    = {0.9,0.9,0.9,0.9};
  G4double pdTeflonAbsLength[TNbEntries]       = {0.2*mm,0.2*mm,0.2*mm,0.2*mm};
  G4double pdTeflonScatteringLength[TNbEntries] = {30*cm,30*cm,30*cm,30*cm};

  G4MaterialPropertiesTable *pTeflonPropertiesTable = new G4MaterialPropertiesTable();

  pTeflonPropertiesTable->AddProperty("RINDEX", pdTeflonPhotonMomentum, pdTeflonRefractiveIndex, TNbEntries);

  pTeflonPropertiesTable->AddProperty("ABSLENGTH", pdTeflonPhotonMomentum, pdTeflonAbsLength, TNbEntries);  
  pTeflonPropertiesTable->AddProperty("RAYLEIGH", pdTeflonPhotonMomentum, pdTeflonScatteringLength, TNbEntries);
  pTeflonPropertiesTable->AddProperty("REFLECTIVITY", pdTeflonPhotonMomentum, pdTeflonReflectivity, TNbEntries);
  Teflon->SetMaterialPropertiesTable(pTeflonPropertiesTable);

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

  auto mylar = G4Material::GetMaterial("mylar");
  auto HDPE = G4Material::GetMaterial("HDPE");
  auto polystyrene = G4Material::GetMaterial("polystyrene");
  auto Teflon = G4Material::GetMaterial("Teflon");

  
  
 
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

  G4double sensSize = size;
  //G4double sipmSize = 3.16 * mm;
  // *** NEW: High-resolution SiPM ***
  G4double sipmSize = 1.0 * mm;  

  // -------------------------- //
  // collimator cell + array
  // -------------------------- //
  G4double cellSize = sipmSize + 0.84 * mm;     // = 1.84 mm
  // Collimator parameters based on paper optimization (15-20mm optimal)
  G4double collimatorLength = 10.0 * mm;  // Optimal from paper: 15-20mm
  G4double cellLength = collimatorLength;  // Now using optimized length
  //G4double cellLength = 30.0 * mm;
  G4double nCells = 25;
  // *** NEW: number of cells that fit into 100 mm ***
  //G4double nCells = std::floor(fBoxSize/cellSize);        //54;   // 54 Ã— 1.84 = 99.36 mm â‰ˆ 100 mm
  G4double collSize = nCells * cellSize;

  G4Box* sBlock = new G4Box("cellBlock",
      cellLength / 2, cellSize / 2, cellSize / 2);

  G4Box* sHole = new G4Box("cellHole",
      cellLength / 2, sipmSize / 2, sipmSize / 2);

  G4SubtractionSolid* sCell = new G4SubtractionSolid("cell",
	  sBlock, sHole, 0, G4ThreeVector(0, 0, 0));

  G4LogicalVolume* fLCell = new G4LogicalVolume(sCell,
	  Teflon, "Cell");

    G4int copyNo = 0;
  G4double xPos, yPos;
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
              G4ThreeVector(xPos, yPos, 0), fLCell, "Cell", fLBox, false, copyNo++, true);

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
  mylarSurf->SetModel(unified); //unified / glisur

  const G4int nOptic = 6;
  G4double photonEnergy[nOptic] = {
      2.07 * eV, 2.34 * eV, 2.64 * eV, 2.75 * eV, 2.95 * eV, 3.10 * eV
      //0.1 * eV, 2.34 * eV, 2.64 * eV, 2.75 * eV, 2.95 * eV, 100.0 * eV
  };
  //G4double mylarRefl[nOptic] = { 0.90, 0.90, 0.89, 0.88, 0.87, 0.86 };
  G4double mylarRefl[nOptic] = { 1, 1, 1, 1, 1, 1 };
  G4double mylarTrans[nOptic] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  G4double SiPMeff[nOptic] = { 0.25, 0.35, 0.45, 0.48, 0.50, 0.45 };

  auto sipm_mtp = new G4MaterialPropertiesTable();
  sipm_mtp->AddProperty("EFFICIENCY", photonEnergy, SiPMeff, nOptic);
  sipmSurf->SetMaterialPropertiesTable(sipm_mtp);

  auto mylarS_mtp = new G4MaterialPropertiesTable();
  mylarS_mtp->AddProperty("REFLECTIVITY", photonEnergy, mylarRefl, nOptic);
  //mylarS_mtp->AddProperty("TRANSMITTANCE", photonEnergy, mylarTrans, nOptic);
  mylarS_mtp->AddProperty("TRANSMITTANCE", photonEnergy, mylarTrans, nOptic);
  mylarSurf->SetMaterialPropertiesTable(mylarS_mtp);
  // ------------------------------------------------------------------------------


  // -------------------------- //
  // mylar foils - electrodes
  // -------------------------- //
  G4double mylThickness = 0.012 * mm; // 0.012 * mm; 0.1 * mm
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
  G4double convThickness = 0.01 * mm;
  G4double convPos = mylPos / 2 + mylThickness / 2 + convThickness / 2;

  G4Box* sConv = new G4Box("conv",
	  collSize / 2, collSize / 2, convThickness / 2);

  G4LogicalVolume* fLConv = new G4LogicalVolume(sConv,
	  HDPE, "Conv");

  G4VPhysicalVolume* fPConv = new G4PVPlacement(0,
	  G4ThreeVector(0, 0, convPos), fLConv, "Conv", fLBox, false, 0, true);


   // -------------------------- //
  // SiPMs
  // -------------------------- //
  G4double sipmWidth = 1.0 * mm;

  G4Box* sSiPM = new G4Box("sipm",
      sipmWidth/2, sipmSize/2, sipmSize/2);

  G4LogicalVolume* fLSiPM = new G4LogicalVolume(sSiPM,
      sensor, "SiPM");

  copyNo = 0;
  for (G4int i = 0; i < 4; i++) {

      G4double angle = i * 90. * deg;
      auto sipmRot = new G4RotationMatrix();
      sipmRot->rotateZ(angle);

      for (G4int j = 0; j < nCells; j++) {

          if (i == 0 || i == 2) {
              xPos = (collSize / 2 + cellLength + sipmWidth/2) * std::cos(angle);
              yPos = ((-1.0 * nCells / 2 + 1.0 * j + 1.0 / 2) * cellSize) + (collSize / 2 + cellLength / 2) * std::sin(angle);
          }
          else {
              xPos = ((-1.0 * nCells / 2 + 1.0 * j + 1.0 / 2) * cellSize) + (collSize / 2 + cellLength / 2) * std::cos(angle);
              yPos = (collSize / 2 + cellLength + sipmWidth / 2) * std::sin(angle);
          }

          G4VPhysicalVolume* fPSiPM = new G4PVPlacement(sipmRot,
              G4ThreeVector(xPos, yPos, 0), fLSiPM, "SiPM", fLBox, false, copyNo++, true);

          auto fLSiPM_surf = new G4LogicalBorderSurface("SiPM-Air", fPSiPM, fPBox, sipmSurf);

      }
  }
  new G4LogicalBorderSurface("MylarA-Air_in", fPBox, fPMylA, mylarSurf);
  new G4LogicalBorderSurface("MylarA-Air_out", fPMylA, fPBox, mylarSurf);

  new G4LogicalBorderSurface("MylarK-Air_in", fPBox, fPMylK, mylarSurf);
  new G4LogicalBorderSurface("MylarK-Air_out", fPMylK, fPBox, mylarSurf);
  


  /*
  G4int copyNo = 0;
  for (G4int i = 0; i < nFibers; i++) {
      for (G4int j = 0; j < 2; j++) {

          G4double angle = (1.0*i + 1.0*j/2) * wsfAngle;
          G4double setRadius;
          if (j == 0) { setRadius = detRadius; }
		  else { setRadius = detRadius2; }
          G4double xPos = setRadius * std::cos(angle);
          G4double yPos = setRadius * std::sin(angle);

          // Placement of ring of fibers 
          // (must be here so fPWSF is defined for defining surface optical properties)

          G4VPhysicalVolume* fPWSF = new G4PVPlacement(WSFRot,
              G4ThreeVector(xPos, yPos, 0), fLWSF, "WSF", fLBox, false, copyNo++, true);

          auto fLwsfSurf = new G4LogicalBorderSurface("WSF_det", fPWSF, fPBox, wsfSurf);

      }   
  }
 
  //auto fLwsfSurf_b = new G4LogicalBorderSurface("WSF_air", fPWSF, fPBox, wsfSurf);

  //auto fLMylarScint_insurf =  new G4LogicalBorderSurface("MylarScint_in", fPDet, fPMyl, mylarSurf);
  //auto fLMylarScint_outsurf =  new G4LogicalBorderSurface("MylarScint_out", fPMyl, fPDet, mylarSurf);

  //auto fLMylarAir_insurf = new G4LogicalBorderSurface("MylarAir_in", fPBox, fPMyl, mylarSurf);
  //auto fLMylarAir_outsurf = new G4LogicalBorderSurface("MylarAir_out", fPMyl, fPBox, mylarSurf);

  new G4LogicalBorderSurface("Mylar1Scint_in", fPDet1, fPMyl1, mylarSurf);
  new G4LogicalBorderSurface("Mylar1Scint_out", fPMyl1, fPDet1, mylarSurf);

  new G4LogicalBorderSurface("Mylar2Scint_in", fPDet2, fPMyl2, mylarSurf);
  new G4LogicalBorderSurface("Mylar2Scint_out", fPMyl2, fPDet2, mylarSurf);

  new G4LogicalBorderSurface("Mylar3Scint_in", fPDet3, fPMyl3, mylarSurf);
  new G4LogicalBorderSurface("Mylar3Scint_out", fPMyl3, fPDet3, mylarSurf);

  new G4LogicalBorderSurface("Mylar4Scint_in", fPDet4, fPMyl4, mylarSurf);
  new G4LogicalBorderSurface("Mylar4Scint_out", fPMyl4, fPDet4, mylarSurf);
  
  new G4LogicalBorderSurface("Mylar1Air_in", fPBox, fPMyl1, mylarSurf);
  new G4LogicalBorderSurface("Mylar1Air_out", fPMyl1, fPBox, mylarSurf);

  new G4LogicalBorderSurface("Mylar4Air_in", fPBox, fPMyl4, mylarSurf);
  new G4LogicalBorderSurface("Mylar4Air_out", fPMyl4, fPBox, mylarSurf);
  */

  /*
  // Cilindro externo de aluminio
  G4double outer_radius = 10. * cm;
  G4double inner_radius = 9.98 * cm;
  G4double height = 55./2 * cm;
  auto detectorS = new G4Tubs("Detector", inner_radius, outer_radius, height, 0. * deg, 360. * deg);
  auto detectorLV = new G4LogicalVolume(detectorS, aluminio, "Detector");
  auto rotation = new G4RotationMatrix();
  rotation->rotateX(90. * deg); // Rotar el cilindro 90 grados en el eje X
  auto detectorPV = new G4PVPlacement(rotation, G4ThreeVector(0.0, 0.0, 0.0), detectorLV, "Detector", worldLV, false, 0);

  
  G4double out_radius = 9.98 * cm;
  G4double in_radius = 0.0 * cm;
  G4double h = 55./2 * cm;
  auto interior_detectorS = new G4Tubs("InteriorDetector", in_radius, out_radius, h, 0. * deg, 360. * deg);
  auto interior_detectorLV = new G4LogicalVolume(interior_detectorS, air, "InteriorDetector");
  auto interior_detectorPV = new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, 0), interior_detectorLV, "InteriorDetector", detectorLV, false, 0);
  

  auto mylarS = new G4Tubs("Mylar", 55 * mm, 55.25 * mm, 50.0/2 * cm, 0. * deg, 360. * deg);
  auto mylarLV = new G4LogicalVolume(mylarS, mylar, "Mylar");
  //auto mylarS = new G4Box("Mylar", 4.2 * cm, 92.25 * cm, 0.0125 * mm);
  auto mylarPV = new G4PVPlacement(rotation, G4ThreeVector(0, 0*mm, -0*cm), mylarLV, "Mylar", worldLV, false, 0);

  auto centelladorS = new G4Tubs("Centellador", 54.5 * mm, 55. * mm, 50.0/2 * cm, 0. * deg, 360. * deg);
  auto centelladorLV = new G4LogicalVolume(centelladorS, centellador, "Centellador");
  auto centelladorPV = new G4PVPlacement(rotation, G4ThreeVector(0, 0 * mm, -0 * cm), centelladorLV, "Centellador", worldLV, false, 0);

  // Definición del semi-largo de la fibra 
  G4double fiberHalfLength = 27.2 * cm; // semi-largo de la fibra (importante: no el largo total)

  // SiPMs 
  G4double sipmHalfThickness = 0.5 * mm; // semi-largo del SiPM
  G4double margin = 0.05 * mm; // margen de seguridad recomendado

  G4double SiPM_z = fiberHalfLength + sipmHalfThickness + margin;

  // Definición y colocación SiPM superior
  auto SiPMS = new G4Tubs("SiPM", 53.35 * mm, 54.35 * mm, sipmHalfThickness, 0. * deg, 360. * deg);
  auto SiPMLV = new G4LogicalVolume(SiPMS, silicon, "SiPM");
  auto SiPMPV = new G4PVPlacement(rotation, G4ThreeVector(0, +SiPM_z, 0), SiPMLV, "SiPM", worldLV, false, 0);

  // Definición y colocación SiPM inferior
  auto SiPM2S = new G4Tubs("SiPM2", 53.35 * mm, 54.35 * mm, sipmHalfThickness, 0. * deg, 360. * deg);
  auto SiPM2LV = new G4LogicalVolume(SiPM2S, silicon, "SiPM2");
  auto SiPM2PV = new G4PVPlacement(rotation, G4ThreeVector(0, -SiPM_z, 0), SiPM2LV, "SiPM2", worldLV, false, 0);

  // Fibras
  auto fiberS = new G4Tubs("Fiber", 0.0 * cm, 0.5 * mm, fiberHalfLength, 0. * deg, 360. * deg);
  auto fiberLV = new G4LogicalVolume(fiberS, polystyrene, "Fiber");

  G4int nFibers = 18;
  G4double fiberRadius = 0.5 * mm;
  // G4double fiberHalfLength = 27.2 * cm; // ya definido arriba
  G4double fiberLengthH = 10 * cm; // longitud horizontal que quieres


  G4double placementRadius = (53.35 * mm + 54.35 * mm) / 2.0; // = 53.85 mm

  //Fibra horizontal 
  auto fiberHorizontalS = new G4Tubs("FiberHorizontal", 0., fiberRadius, fiberLengthH / 2., 0. * deg, 360. * deg);
  auto fiberHorizontalLV = new G4LogicalVolume(fiberHorizontalS, polystyrene, "FiberHorizontal");
  *//*

  // Colocación de las fibras verticales
  for (int i = 0; i < nFibers; i++) {
      double angle = 2.0 * CLHEP::pi * i / nFibers;
      double x = placementRadius * std::cos(angle);
      double z = placementRadius * std::sin(angle);

      // Colocar fibra vertical (centrada en z=0)
      G4ThreeVector posVertical(x, 0.0, z);
      auto fiberPV = new G4PVPlacement(rotation, posVertical, fiberLV, "Fiber", worldLV, false, i, fCheckOverlaps);

      // ... el resto de tu código (óptica, etc) ...
  

      // Superficie óptica entre centellador y fibras
      auto opticalSurface = new G4OpticalSurface("ScintillatorToSiPM");
      opticalSurface->SetType(dielectric_dielectric);
      opticalSurface->SetModel(unified);
      opticalSurface->SetFinish(polished);  // también puedes probar ground o groundbackpainted

      const G4int nEntries = 6;
      G4double photonEnergy[nEntries] = {
        2.07 * eV, 2.34 * eV, 2.64 * eV, 2.75 * eV, 2.95 * eV, 3.10 * eV
      };

      // 100% transmisión (puedes ajustar esto)
      G4double reflectivity[nEntries] = { 0., 0., 0., 0., 0., 0. }; // No refleja
      G4double efficiency[nEntries] = { 0.25, 0.35, 0.45, 0.48, 0.50, 0.45 }; // PDE del SiPM

      auto surfaceMPT = new G4MaterialPropertiesTable();
      surfaceMPT->AddProperty("REFLECTIVITY", photonEnergy, reflectivity, nEntries);
      surfaceMPT->AddProperty("EFFICIENCY", photonEnergy, efficiency, nEntries); // Activa la detección

      opticalSurface->SetMaterialPropertiesTable(surfaceMPT);
      auto opticalSurfacePV = new G4LogicalBorderSurface("ScintillatorToSiPM_Surface", fiberPV, SiPMPV, opticalSurface);

      auto opticalSurface2PV = new G4LogicalBorderSurface("ScintillatorToSiPM_Surface", fiberPV, SiPM2PV, opticalSurface);



  }



  auto centellador2S = new G4Tubs("Centellador2", 52.7 * mm, 53.2 * mm, 50.0/2 * cm, 0. * deg, 360. * deg);
  auto centellador2LV = new G4LogicalVolume(centellador2S, centellador, "Centellador2");
  auto centellador2PV = new G4PVPlacement(rotation, G4ThreeVector(0, 0 * mm, -0 * cm), centellador2LV, "Centellador2", worldLV, false, 0);

  auto mylar2S = new G4Tubs("Mylar2", 52.45 * mm, 52.7 * mm, 50.0/2 * cm, 0. * deg, 360. * deg);
  auto mylar2LV = new G4LogicalVolume(mylar2S, mylar, "Mylar2");
  auto mylar2PV = new G4PVPlacement(rotation, G4ThreeVector(0, 0 * mm, -0 * cm), mylar2LV, "Mylar", worldLV, false, 0);




  */
  /*
  auto mylar2S = new G4Box("Mylar2", 4.2 * cm, 92.25 * cm, 0.0125 * mm);
  auto mylar2LV = new G4LogicalVolume(mylar2S, mylar, "Mylar2");
  auto mylar2PV = new G4PVPlacement(rotation, G4ThreeVector(0, +0.7625*mm, -12.75*cm), mylar2LV, "Mylar2", interior_detectorLV, false, 0);

  auto centelladorS = new G4Box("Centellador", 4.2 * cm, 92.25 * cm, 0.125 * mm);
  auto centelladorLV = new G4LogicalVolume(centelladorS, centellador, "Centellador");
  auto centelladorPV = new G4PVPlacement(rotation, G4ThreeVector(0, -0.625*mm, -12.75*cm), centelladorLV, "Centellador", interior_detectorLV, false, 0);

  auto centellador2S = new G4Box("Centellador2", 4.2 * cm, 92.25 * cm, 0.125 * mm);
  auto centellador2LV = new G4LogicalVolume(centellador2S, centellador, "Centellador2");
  auto centellador2PV = new G4PVPlacement(rotation, G4ThreeVector(0, +0.625*mm, -12.75*cm), centellador2LV, "Centellador2", interior_detectorLV, false, 0);

  auto fiberS = new G4Tubs("Fiber", 0.0 * cm, 0.5 * mm, 100.0 * cm, 0. * deg, 360. * deg);
  auto fiberLV = new G4LogicalVolume(fiberS, polystyrene, "Fiber");
  //auto fiberPV = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -5.0*cm), fiberLV, "Fiber", interior_detectorLV, false, 0);

  G4int nfibers = 12;
  G4double spacing = 7.0 * mm;
  G4double startY = -((nfibers - 1) / 2.0) * spacing; // centra el arreglo

  for (G4int i = 0; i < nfibers; ++i) {
      G4double posY = startY + i * spacing;
      new G4PVPlacement(
          nullptr,
          G4ThreeVector(posY, 0, -5*cm), // fibras alineadas en X, separadas en Y
          fiberLV,
          "Fiber",
          interior_detectorLV,
          false,
          i
      );
  }

  auto sipmS = new G4Box("SiPM", 42 * mm, 1 * mm, 0.5 * mm);
  auto sipmLV = new G4LogicalVolume(sipmS, silicon, "SiPM");
  auto sipmPV = new G4PVPlacement(rotation, G4ThreeVector(0, 0, +0.95*m), sipmLV, "SiPM", interior_detectorLV, false, 0);

  auto optic_sipmS = new G4Box("OpticSiPM", 42 * mm, 1 * mm, 0.1 * mm);
  auto optic_sipmLV = new G4LogicalVolume(optic_sipmS, sio2, "OpticSiPM");
  auto optic_sipmPV = new G4PVPlacement(rotation, G4ThreeVector(0, 0, 0.9503*m), optic_sipmLV, "OpticSiPM", interior_detectorLV, false, 0);
  */

  fLCell->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0.0))); // Yellow
  fLMyl->SetVisAttributes(G4VisAttributes(G4Colour(0.8, 0.8, 0.8, 0.9))); // Grey
  fLConv->SetVisAttributes(G4VisAttributes(G4Colour(0.0, 1.0, 0.0))); // Green
  fLSiPM->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 0.0, 0.0))); // Red
  
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
}

} 
