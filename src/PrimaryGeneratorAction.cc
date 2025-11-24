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
/// \file B4/B4a/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the B4::PrimaryGeneratorAction class
#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4AnalysisManager.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

namespace B4
{
    PrimaryGeneratorAction::PrimaryGeneratorAction()
    {
        G4int nofParticles = 1;
        fParticleGun = new G4ParticleGun(nofParticles);

        auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
        fParticleGun->SetParticleDefinition(particleDefinition);
        
        // Fixed energy: 2.5 MeV
        fParticleGun->SetParticleEnergy(2.5 * MeV);
        
        // Default direction (can be overridden in GeneratePrimaries)
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    }

    PrimaryGeneratorAction::~PrimaryGeneratorAction()
    {
        delete fParticleGun;
    }

    void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
    {
        // Energy is already set to 2.5 MeV in constructor
        
        // Set gun position - neutron source positioned to hit the detector
        G4double xpos = 0.0 * cm;    // Center in X
        G4double ypos = 0.0 * cm;    // Center in Y  
        G4double zpos = 0.6 * cm;    // Positioned to hit the HDPE converter at z ≈ 0.5cm

        fParticleGun->SetParticlePosition(G4ThreeVector(xpos, ypos, zpos));

        // Set direction - pointing into the detector (negative Z direction)
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, -1.0));

        // Generate the primary vertex
        fParticleGun->GeneratePrimaryVertex(anEvent);
        
        // Optional: Log the event for debugging
        G4cout << "Generated 2.5 MeV neutron at (" << xpos/cm << ", " << ypos/cm << ", " << zpos/cm << ") cm" << G4endl;
    }
}
