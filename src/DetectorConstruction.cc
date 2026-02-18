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
// * use.  Please see the license in the file  LICENSE and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  agree  to acknowledge its     *
// * use  in  resulting  scientific  publications,  and indicate      *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file B4/B4a/src/DetectorConstruction.cc
/// \brief Implementation of the B4::DetectorConstruction class

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

#include "EventAction.hh"
#include "G4RunManager.hh"

namespace B4
{

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

// ============================================================================
// Main Construction
// ============================================================================
G4VPhysicalVolume* DetectorConstruction::Construct()
{
    pitch = 0.5 * cm;
    size = 5.0 * cm;
    
    DefineMaterials();
    return DefineVolumes(pitch, size);
}

// ============================================================================
// Material Definitions
// ============================================================================
void DetectorConstruction::DefineMaterials()
{
    G4int ncomponents, natoms;
    G4double massfraction;
    G4double a, z, dens;
    
    G4double Vdens = 1.e-25 * g / cm3;
    G4double Vpres = 1.e-19 * pascal;
    G4double Vtemp = 0.1 * kelvin;
    
    auto nistManager = G4NistManager::Instance();
    
    nistManager->FindOrBuildElement("H");
    nistManager->FindOrBuildElement("He");
    nistManager->FindOrBuildElement("C");
    nistManager->FindOrBuildElement("N");
    nistManager->FindOrBuildElement("O");
    nistManager->FindOrBuildElement("F");
    nistManager->FindOrBuildElement("Na");
    nistManager->FindOrBuildElement("Mg");
    nistManager->FindOrBuildElement("Al");
    nistManager->FindOrBuildElement("Si");
    nistManager->FindOrBuildElement("Ar");
    nistManager->FindOrBuildElement("K");
    nistManager->FindOrBuildElement("Ca");
    nistManager->FindOrBuildElement("Cr");
    nistManager->FindOrBuildElement("Mn");
    nistManager->FindOrBuildElement("Fe");
    nistManager->FindOrBuildElement("Ni");
    nistManager->FindOrBuildElement("Zn");
    nistManager->FindOrBuildElement("Ag");
    nistManager->FindOrBuildElement("I");
    nistManager->FindOrBuildElement("Cs");
    nistManager->FindOrBuildElement("Pb");
    nistManager->FindOrBuildElement("W");
    
    G4Material* air = nistManager->FindOrBuildMaterial("G4_AIR");
    G4Material* steel = nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    G4Material* water = nistManager->FindOrBuildMaterial("G4_WATER");
    G4Material* aluminum = nistManager->FindOrBuildMaterial("G4_Al");
    
    G4Material* vacuum = new G4Material("vacuum", z = 1, a = 1.01 * g / mole, 
                                        Vdens, kStateGas, Vtemp, Vpres);
    
    // CF4 gas
    auto C = G4Element::GetElement("C");
    auto F = G4Element::GetElement("F");
    G4Material* CF4 = new G4Material("CF4", 3.72 * mg / cm3, ncomponents = 2, 
                                     kStateGas, 293.15 * kelvin, 1.0 * atmosphere);
    CF4->AddElement(C, natoms = 1);
    CF4->AddElement(F, natoms = 4);
    
    // Argon gas
    auto Ar = G4Element::GetElement("Ar");
    G4Material* Ar_gas = new G4Material("Ar_gas", 1.782 * mg / cm3, ncomponents = 1, 
                                        kStateGas, 293.15 * kelvin, 1.0 * atmosphere);
    Ar_gas->AddElement(Ar, natoms = 1);
    
    // Ar:CF4 detector gas (90:10)
    G4Material* Ar_CF4 = new G4Material("Ar_CF4", 0.061 * mg / cm3, ncomponents = 2, 
                                        kStateGas, 293.15 * kelvin, 30.e-3 * bar);
    Ar_CF4->AddMaterial(Ar_gas, massfraction = 0.9);
    Ar_CF4->AddMaterial(CF4, massfraction = 0.1);
    
    // Polystyrene
    auto H = G4Element::GetElement("H");
    G4Material* polystyrene = new G4Material("polystyrene", 1.05 * g / cm3, 2);
    polystyrene->AddElement(C, natoms = 8);
    polystyrene->AddElement(H, natoms = 8);
    
    // PMMA
    auto O = G4Element::GetElement("O");
    G4Material* PMMA = new G4Material("PMMA", 1.20 * g / cm3, 3);
    PMMA->AddElement(C, natoms = 5);
    PMMA->AddElement(H, natoms = 8);
    PMMA->AddElement(O, natoms = 2);
    
    // Mylar
    G4Material* mylar = new G4Material("mylar", 1.39 * g / cm3, 3);
    mylar->AddElement(C, natoms = 5);
    mylar->AddElement(H, natoms = 4);
    mylar->AddElement(O, natoms = 2);
    
    // SiPM sensor
    auto Si = G4Element::GetElement("Si");
    G4Material* sensor = new G4Material("SiPM", 2.33 * g / cm3, 2);
    sensor->AddElement(Si, natoms = 1);
    sensor->AddElement(O, natoms = 2);
    
    // Aluminum
    auto Al = G4Element::GetElement("Al");
    G4Material* alum = new G4Material("alum", 2.7 * g / cm3, 1);
    alum->AddElement(Al, natoms = 1);
    
    // HDPE converter
    G4Element* H_hp = new G4Element("TS_H_of_Polyethylene", "H", 1, 1.0079 * g / mole);
    G4Material* HDPE = new G4Material("HDPE", 0.93 * g / cm3, 2);
    HDPE->AddElement(C, 2);
    HDPE->AddElement(H_hp, 4);
    
    // Teflon collimator
    G4Material* Teflon = new G4Material("Teflon", 2.2 * g / cm3, 2, kStateSolid);
    Teflon->AddElement(C, 0.240183);
    Teflon->AddElement(F, 0.759817);
    
    // Teflon optical properties
    const G4int TNbEntries = 4;
    G4double pdTeflonPhotonMomentum[TNbEntries] = {
        2.07 * eV, 2.34 * eV, 2.64 * eV, 3.10 * eV
    };
    G4double pdTeflonRefractiveIndex[TNbEntries] = {1.34, 1.34, 1.34, 1.34};
    G4double pdTeflonReflectivity[TNbEntries] = {0.95, 0.95, 0.95, 0.95};
    G4double pdTeflonAbsLength[TNbEntries] = {1.0 * m, 1.0 * m, 1.0 * m, 1.0 * m};
    
    G4MaterialPropertiesTable* pTeflonPropertiesTable = new G4MaterialPropertiesTable();
    pTeflonPropertiesTable->AddProperty("RINDEX", pdTeflonPhotonMomentum, 
                                        pdTeflonRefractiveIndex, TNbEntries);
    pTeflonPropertiesTable->AddProperty("ABSLENGTH", pdTeflonPhotonMomentum, 
                                        pdTeflonAbsLength, TNbEntries);
    pTeflonPropertiesTable->AddProperty("REFLECTIVITY", pdTeflonPhotonMomentum, 
                                        pdTeflonReflectivity, TNbEntries);
    Teflon->SetMaterialPropertiesTable(pTeflonPropertiesTable);
    
    // ========================================================================
    // Optical Properties
    // ========================================================================
    
    // Air & vacuum
    auto air_mpt = new G4MaterialPropertiesTable();
    air_mpt->AddProperty("RINDEX", "Air");
    air->SetMaterialPropertiesTable(air_mpt);
    vacuum->SetMaterialPropertiesTable(air_mpt);
    
    // Scintillating gas
    const G4int nScint = 18;
    G4double gasEnergy[nScint] = {
        1.55*eV, 1.61*eV, 1.70*eV, 1.78*eV, 1.88*eV,
        2.00*eV, 2.11*eV, 2.22*eV, 2.54*eV, 2.87*eV,
        3.42*eV, 3.89*eV, 4.35*eV, 4.72*eV, 4.88*eV,
        5.40*eV, 5.83*eV, 6.16*eV
    };
    G4double gasScintSp[nScint] = {
        0.01, 0.01, 0.04, 0.10, 0.20,
        0.26, 0.18, 0.04, 0.00, 0.02,
        0.04, 0.13, 0.30, 0.19, 0.25,
        0.03, 0.00, 0.00
    };
    G4double gasAbsLength[nScint] = {
        4*m, 4*m, 4*m, 4*m, 4*m,
        4*m, 4*m, 4*m, 4*m, 4*m,
        4*m, 4*m, 4*m, 4*m, 4*m,
        4*m, 4*m, 4*m
    };
    
    auto ArCF4_mpt = new G4MaterialPropertiesTable();
    ArCF4_mpt->AddProperty("SCINTILLATIONCOMPONENT1", gasEnergy, gasScintSp, nScint);
    ArCF4_mpt->AddProperty("RINDEX", "Air");
    ArCF4_mpt->AddConstProperty("SCINTILLATIONYIELD", 2500000 / MeV);
    ArCF4_mpt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
    ArCF4_mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 15 * ns);
    ArCF4_mpt->AddProperty("ABSLENGTH", gasEnergy, gasAbsLength, nScint);
    ArCF4_mpt->AddConstProperty("RESOLUTIONSCALE", 1.0);
    Ar_CF4->SetMaterialPropertiesTable(ArCF4_mpt);
    Ar_gas->SetMaterialPropertiesTable(ArCF4_mpt);
    
    // SiPM sensors
    const G4int nSiPM = 6;
    G4double SiPMEnergy[nSiPM] = {
        2.07*eV, 2.34*eV, 2.64*eV, 2.75*eV, 2.95*eV, 3.10*eV
    };
    G4double SiPMQE[nSiPM] = {
        0.25, 0.35, 0.45, 0.48, 0.50, 0.45
    };
    G4double SiPMrIndex[nSiPM] = {
        3.4, 3.4, 3.4, 3.4, 3.4, 3.4
    };
    auto SiPM_mpt = new G4MaterialPropertiesTable();
    SiPM_mpt->AddProperty("RINDEX", SiPMEnergy, SiPMrIndex, nSiPM);
    SiPM_mpt->AddProperty("EFFICIENCY", SiPMEnergy, SiPMQE, nSiPM);
    sensor->SetMaterialPropertiesTable(SiPM_mpt);
    
    // Mylar & Aluminum
    const G4int nAl = 6;
    G4double AlEnergy[nAl] = {
        2.07*eV, 2.34*eV, 2.64*eV, 2.75*eV, 2.95*eV, 3.10*eV
    };
    G4double AlrIndex[nAl] = {
        1.373, 1.373, 1.373, 1.373, 1.373, 1.373
    };
    G4double AlAbsLength[nAl] = {
        4*cm, 4*cm, 4*cm, 4*cm, 4*cm, 4*cm
    };
    auto Al_mtp = new G4MaterialPropertiesTable();
    Al_mtp->AddProperty("RINDEX", AlEnergy, AlrIndex, nAl);
    Al_mtp->AddProperty("ABSLENGTH", AlEnergy, AlAbsLength, nAl);
    alum->SetMaterialPropertiesTable(Al_mtp);
    mylar->SetMaterialPropertiesTable(Al_mtp);
    
    G4cout << "Materials defined successfully!" << G4endl;
}

// ============================================================================
// Volume Definitions
// ============================================================================
G4VPhysicalVolume* DetectorConstruction::DefineVolumes(G4double pitch, G4double size)
{
    auto vacuum = G4Material::GetMaterial("vacuum");
    auto atmosphere = G4Material::GetMaterial("Ar_CF4");
    auto mirror = G4Material::GetMaterial("alum");
    auto sensor = G4Material::GetMaterial("SiPM");
    auto Teflon = G4Material::GetMaterial("Teflon");
    auto mylar = G4Material::GetMaterial("mylar");
    auto HDPE = G4Material::GetMaterial("HDPE");
    auto steel = G4Material::GetMaterial("G4_STAINLESS-STEEL");
    auto water = G4Material::GetMaterial("G4_WATER");
    auto aluminum = G4Material::GetMaterial("G4_Al");
    
    // ========================================================================
    // World Volume
    // ========================================================================
    G4double worldSize = 60 * cm;
    G4double worldZ = 100 * cm;
    
    G4Box* worldBox = new G4Box("World", worldSize/2, worldSize/2, worldZ/2);
    G4LogicalVolume* logicWorld = new G4LogicalVolume(worldBox, vacuum, "World");
    G4VPhysicalVolume* physWorld = new G4PVPlacement(0,
                                                     G4ThreeVector(),
                                                     logicWorld,
                                                     "World",
                                                     0,
                                                     false,
                                                     0,
                                                     fCheckOverlaps);
    
    // ========================================================================
    // Valve Phantom
    // ========================================================================
    G4double valveRadius = 8.0 * mm;
    G4double valveLength = 4.0 * mm;
    G4double cavityRadius = 5.0 * mm - 0.01*mm;
    G4double cavityLength = 4.0 * mm - 0.02*mm;
    G4double stemRadius = 1.95 * mm;
    G4double stemLength = 3.5 * mm;
    
    // Valve body (hollow steel cylinder)
    G4Tubs* solidValveOuter = new G4Tubs("ValveOuter",
                                         0,
                                         valveRadius,
                                         valveLength/2,
                                         0,
                                         360*deg);
    
    G4Tubs* solidValveHole = new G4Tubs("ValveHole",
                                        0,
                                        cavityRadius,
                                        cavityLength/2,
                                        0,
                                        360*deg);
    
    G4SubtractionSolid* solidValveBody = new G4SubtractionSolid("ValveBody",
                                                                solidValveOuter,
                                                                solidValveHole);
    
    G4LogicalVolume* logicValveBody = new G4LogicalVolume(solidValveBody,
                                                          steel,
                                                          "ValveBody");
    
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, -15*cm),
                      logicValveBody,
                      "ValveBody",
                      logicWorld,
                      false,
                      0,
                      fCheckOverlaps);
    
    // Water cavity
    G4Tubs* solidCavity = new G4Tubs("Cavity",
                                     0,
                                     cavityRadius,
                                     cavityLength/2,
                                     0,
                                     360*deg);
    
    G4LogicalVolume* logicCavity = new G4LogicalVolume(solidCavity,
                                                       water,
                                                       "Cavity");
    
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, 0),
                      logicCavity,
                      "Cavity",
                      logicValveBody,
                      false,
                      0,
                      fCheckOverlaps);
    
    // Aluminum valve stem
    G4Tubs* solidStem = new G4Tubs("ValveStem",
                                   0,
                                   stemRadius,
                                   stemLength/2,
                                   0,
                                   360*deg);
    
    G4LogicalVolume* logicStem = new G4LogicalVolume(solidStem,
                                                     aluminum,
                                                     "ValveStem");
    
    G4RotationMatrix* stemRot = new G4RotationMatrix();
    stemRot->rotateY(90*deg);
    
    new G4PVPlacement(stemRot,
                      G4ThreeVector(0, 0, 0),
                      logicStem,
                      "ValveStem",
                      logicCavity,
                      false,
                      0,
                      fCheckOverlaps);
    
    // ========================================================================
    // Detector Box
    // ========================================================================
    G4double fBoxSize = 10 * cm;
    G4double fBoxDepth = 1 * cm;
    
    G4Box* sBox = new G4Box("DetectorBox", fBoxSize, fBoxSize, fBoxDepth/2);
    G4LogicalVolume* fLBox = new G4LogicalVolume(sBox, atmosphere, "DetectorBox");
    G4VPhysicalVolume* fPBox = new G4PVPlacement(0,
                                                 G4ThreeVector(0, 0, 0),
                                                 fLBox,
                                                 "DetectorBox",
                                                 logicWorld,
                                                 false,
                                                 0,
                                                 fCheckOverlaps);
    
    // ========================================================================
    // Collimator Array
    // ========================================================================
    G4double sipmSize = 1.0 * mm;
    G4double cellSize = sipmSize + 0.10 * mm;
    G4double cellLength = 3.0 * mm;
    G4int nCells = 25;
    G4double collSize = nCells * cellSize;
    
    G4Box* sBlock = new G4Box("cellBlock", cellLength/2, cellSize/2, cellSize/2);
    G4Box* sHole = new G4Box("cellHole", cellLength/2, sipmSize/2, sipmSize/2);
    
    G4SubtractionSolid* sCell = new G4SubtractionSolid("cell",
                                                       sBlock,
                                                       sHole,
                                                       0,
                                                       G4ThreeVector(0, 0, 0));
    
    G4LogicalVolume* fLCell = new G4LogicalVolume(sCell, Teflon, "Cell");
    
    G4int copyNo = 0;
    G4double xPos, yPos;
    
    for (G4int i = 0; i < 4; i++) {
        G4double angle = i * 90. * deg;
        auto cellRot = new G4RotationMatrix();
        cellRot->rotateZ(angle);
        
        for (G4int j = 0; j < nCells; j++) {
            if (i == 0 || i == 2) {
                xPos = (collSize/2 + cellLength/2) * std::cos(angle);
                yPos = ((-1.0*nCells/2 + 1.0*j + 1.0/2) * cellSize) +
                       (collSize/2 + cellLength/2) * std::sin(angle);
            } else {
                xPos = ((-1.0*nCells/2 + 1.0*j + 1.0/2) * cellSize) +
                       (collSize/2 + cellLength/2) * std::cos(angle);
                yPos = (collSize/2 + cellLength/2) * std::sin(angle);
            }
            
            new G4PVPlacement(cellRot,
                             G4ThreeVector(xPos, yPos, 0),
                             fLCell,
                             "Cell",
                             fLBox,
                             false,
                             copyNo++,
                             fCheckOverlaps);
        }
    }
    
    // ========================================================================
    // Mylar Electrodes
    // ========================================================================
    G4double mylThickness = 0.012 * mm;
    G4double mylPos = cellSize;
    
    G4Box* sMyl = new G4Box("myl", collSize/2, collSize/2, mylThickness/2);
    G4LogicalVolume* fLMyl = new G4LogicalVolume(sMyl, mirror, "Myl");
    
    G4VPhysicalVolume* fPMylA = new G4PVPlacement(0,
                                                  G4ThreeVector(0, 0, mylPos/2),
                                                  fLMyl,
                                                  "MylA",
                                                  fLBox,
                                                  false,
                                                  0,
                                                  fCheckOverlaps);
    
    G4VPhysicalVolume* fPMylK = new G4PVPlacement(0,
                                                  G4ThreeVector(0, 0, -mylPos/2),
                                                  fLMyl,
                                                  "MylK",
                                                  fLBox,
                                                  false,
                                                  0,
                                                  fCheckOverlaps);
    
    // ========================================================================
    // HDPE Converter
    // ========================================================================
    G4double convThickness = 1.5 * mm;
    G4double convPos = mylPos/2 + mylThickness/2 + convThickness/2;
    
    G4Box* sConv = new G4Box("conv", collSize/2, collSize/2, convThickness/2);
    G4LogicalVolume* fLConv = new G4LogicalVolume(sConv, HDPE, "Conv");
    
    new G4PVPlacement(0,
                     G4ThreeVector(0, 0, convPos),
                     fLConv,
                     "Conv",
                     fLBox,
                     false,
                     0,
                     fCheckOverlaps);
    
    // ========================================================================
    // SiPM Arrays
    // ========================================================================
    G4double sipmWidth = 1.0 * mm;
    
    G4Box* sSiPM = new G4Box("sipm", sipmWidth/2, sipmSize/2, sipmSize/2);
    G4LogicalVolume* fLSiPM = new G4LogicalVolume(sSiPM, sensor, "SiPM");
    
    copyNo = 0;
    for (G4int i = 0; i < 4; i++) {
        G4double angle = i * 90. * deg;
        auto sipmRot = new G4RotationMatrix();
        sipmRot->rotateZ(angle);
        
        for (G4int j = 0; j < nCells; j++) {
            if (i == 0 || i == 2) {
                xPos = (collSize/2 + cellLength + sipmWidth/2) * std::cos(angle);
                yPos = ((-1.0*nCells/2 + 1.0*j + 1.0/2) * cellSize) +
                       (collSize/2 + cellLength/2) * std::sin(angle);
            } else {
                xPos = ((-1.0*nCells/2 + 1.0*j + 1.0/2) * cellSize) +
                       (collSize/2 + cellLength/2) * std::cos(angle);
                yPos = (collSize/2 + cellLength + sipmWidth/2) * std::sin(angle);
            }
            
            new G4PVPlacement(sipmRot,
                             G4ThreeVector(xPos, yPos, 0),
                             fLSiPM,
                             "SiPM",
                             fLBox,
                             false,
                             copyNo++,
                             fCheckOverlaps);
        }
    }
    
    // ========================================================================
    // Optical Surfaces
    // ========================================================================
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
        2.07*eV, 2.34*eV, 2.64*eV, 2.75*eV, 2.95*eV, 3.10*eV
    };
    G4double mylarRefl[nOptic] = {1, 1, 1, 1, 1, 1};
    G4double mylarTrans[nOptic] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    G4double SiPMeff[nOptic] = {0.25, 0.35, 0.45, 0.48, 0.50, 0.45};
    
    auto sipm_mtp = new G4MaterialPropertiesTable();
    sipm_mtp->AddProperty("EFFICIENCY", photonEnergy, SiPMeff, nOptic);
    sipmSurf->SetMaterialPropertiesTable(sipm_mtp);
    
    auto mylarS_mtp = new G4MaterialPropertiesTable();
    mylarS_mtp->AddProperty("REFLECTIVITY", photonEnergy, mylarRefl, nOptic);
    mylarS_mtp->AddProperty("TRANSMITTANCE", photonEnergy, mylarTrans, nOptic);
    mylarSurf->SetMaterialPropertiesTable(mylarS_mtp);
    
    new G4LogicalBorderSurface("MylarA-Air_in", fPBox, fPMylA, mylarSurf);
    new G4LogicalBorderSurface("MylarA-Air_out", fPMylA, fPBox, mylarSurf);
    new G4LogicalBorderSurface("MylarK-Air_in", fPBox, fPMylK, mylarSurf);
    new G4LogicalBorderSurface("MylarK-Air_out", fPMylK, fPBox, mylarSurf);
    
    // ========================================================================
    // Visualization Attributes
    // ========================================================================
    G4VisAttributes* cellVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.3));
    cellVisAtt->SetVisibility(true);
    fLCell->SetVisAttributes(cellVisAtt);
    
    fLMyl->SetVisAttributes(G4VisAttributes(G4Colour(0.8, 0.8, 0.8, 0.9)));
    fLConv->SetVisAttributes(G4VisAttributes(G4Colour(0.0, 1.0, 0.0)));
    
    G4VisAttributes* sipmVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
    sipmVisAtt->SetVisibility(true);
    sipmVisAtt->SetForceSolid(true);
    fLSiPM->SetVisAttributes(sipmVisAtt);
    
    G4VisAttributes* valveVis = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.7));
    logicValveBody->SetVisAttributes(valveVis);
    
    G4VisAttributes* cavityVis = new G4VisAttributes(G4Colour(0.0, 0.5, 1.0, 0.5));
    logicCavity->SetVisAttributes(cavityVis);
    
    G4VisAttributes* stemVis = new G4VisAttributes(G4Colour(0.8, 0.8, 0.0, 0.8));
    logicStem->SetVisAttributes(stemVis);
    
    fLBox->SetVisAttributes(G4VisAttributes(G4Colour(0.0, 1.0, 1.0, 0.2)));
    logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());
    
    return physWorld;
}

// ============================================================================
// Sensitive Detector Construction
// ============================================================================
void DetectorConstruction::ConstructSDandField()
{
    G4ThreeVector fieldValue;
    fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
    fMagFieldMessenger->SetVerboseLevel(1);
    G4AutoDelete::Register(fMagFieldMessenger);
    
    auto SDmanager = G4SDManager::GetSDMpointer();
    SDmanager->SetVerboseLevel(1);
    
    auto runManager = G4RunManager::GetRunManager();
    auto eventAction = const_cast<B4a::EventAction*>(
        dynamic_cast<const B4a::EventAction*>(runManager->GetUserEventAction()));
    
    auto sipmSD = new MySensitiveDetector("SiPM_SD");
    
    if (eventAction) {
        sipmSD->SetEventAction(eventAction);
        G4cout << "EventAction set on MySensitiveDetector" << G4endl;
    }
    
    SDmanager->AddNewDetector(sipmSD);
    
    auto logicalVolumeStore = G4LogicalVolumeStore::GetInstance();
    for (auto logicalVolume : *logicalVolumeStore) {
        if (logicalVolume->GetName().find("SiPM") != std::string::npos) {
            logicalVolume->SetSensitiveDetector(sipmSD);
            G4cout << "Assigned sensitive detector to: " 
                   << logicalVolume->GetName() << G4endl;
        }
    }
    
    G4cout << "Sensitive detectors registered successfully!" << G4endl;
}

} // namespace B4

