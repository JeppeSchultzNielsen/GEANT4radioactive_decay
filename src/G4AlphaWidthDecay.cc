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
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4NeutronWidthDecay.cc                                            //
//  Author: L.G. Sarmiento (Lund)                                             //
//  Date:   10 October 2015                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4AlphaWidthDecay.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <iomanip>

G4AlphaWidthDecay::G4AlphaWidthDecay(const G4ParticleDefinition* theParentNucleus,
                                         const G4double& branch, const G4double& Qvalue,
                                         const G4double& excitationE,
                                         const G4Ions::G4FloatLevelBase& flb)
        : G4NuclearDecay("alpha decay", AlphaWidth, excitationE, flb), transitionQ(Qvalue)
{
    SetParent(theParentNucleus);  // Store name of parent nucleus, delete G4MT_parent
    SetBR(branch);

    SetNumberOfDaughters(2);
    theIonTable = (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
    daughterZ = theParentNucleus->GetAtomicNumber()-2;
    daughterA = theParentNucleus->GetAtomicMass() - 4;
    nominalExcitation = excitationE;
    nominalFlb = flb;
    widthReadSucces = ReadWidthFile(daughterZ,daughterA,excitationE*MeV, transitionQ*MeV);
}


G4AlphaWidthDecay::~G4AlphaWidthDecay()
{}


G4DecayProducts* G4AlphaWidthDecay::DecayIt(G4double)
{
    G4DecayProducts* products;
    if(widthReadSucces){
        int chosenLevel = ChooseDecaySublevel();

        double sublevelExcitation = sublevelExs[chosenLevel]/1000;
        double sublevelQvalue = sublevelQvalues[chosenLevel];

        SetDaughter(0, theIonTable->GetIon(daughterZ, daughterA, sublevelExcitation, nominalFlb, nominalExcitation) );
        SetDaughter(1, "alpha");

        // Fill G4MT_parent with theParentNucleus (stored by SetParent in ctor)
        CheckAndFillParent();

        // Fill G4MT_daughters with alpha and residual nucleus (stored by SetDaughter)
        //CheckAndFillDaughters();
        CheckAndFillDaughters();

        G4double alphaMass = G4MT_daughters[1]->GetPDGMass();
        // Excitation energy included in PDG mass
        G4double nucleusMass = G4MT_daughters[0]->GetPDGMass();

        // Q value was calculated from atomic masses.
        // Use it to get correct neutron energy.
        G4double cmMomentum = std::sqrt(sublevelQvalue*(sublevelQvalue + 2.*alphaMass)*
                                        (sublevelQvalue + 2.*nucleusMass)*
                                        (sublevelQvalue + 2.*alphaMass + 2.*nucleusMass) )/
                              (sublevelQvalue + alphaMass + nucleusMass)/2.;

        // Set up final state
        // parentParticle is set at rest here because boost with correct momentum
        // is done later
        G4DynamicParticle parentParticle(G4MT_parent, G4ThreeVector(0,0,0), 0.0);
        products = new G4DecayProducts(parentParticle);

        G4double costheta = 2.*G4UniformRand()-1.0;
        G4double sintheta = std::sqrt(1.0 - costheta*costheta);
        G4double phi  = twopi*G4UniformRand()*rad;
        G4ThreeVector direction(sintheta*std::cos(phi),sintheta*std::sin(phi),
                                costheta);

        G4double KE = std::sqrt(cmMomentum*cmMomentum + alphaMass*alphaMass)
                      - alphaMass;

        /*G4cout << "Decay from " << G4MT_parent -> GetParticleName() << G4endl;
        G4cout << "To alpha" << "with " << KE << G4endl;*/

        G4DynamicParticle* daughterparticle =
                new G4DynamicParticle(G4MT_daughters[1], direction, KE, alphaMass);
        products->PushProducts(daughterparticle);

        KE = std::sqrt(cmMomentum*cmMomentum + nucleusMass*nucleusMass) - nucleusMass;
        //G4cout << "And " << G4MT_daughters[0]->GetParticleName() << " with " << KE << G4endl;
        daughterparticle =
                new G4DynamicParticle(G4MT_daughters[0], -1.0*direction, KE, nucleusMass);
        products->PushProducts(daughterparticle);
    }
    else{
        SetDaughter(0, theIonTable->GetIon(daughterZ, daughterA, nominalExcitation, nominalFlb) );
        SetDaughter(1, "alpha");

        // Fill G4MT_parent with theParentNucleus (stored by SetParent in ctor)
        CheckAndFillParent();

        // Fill G4MT_daughters with neutron and residual nucleus (stored by SetDaughter)
        CheckAndFillDaughters();

        G4double alphaMass = G4MT_daughters[1]->GetPDGMass();
        // Excitation energy included in PDG mass
        G4double nucleusMass = G4MT_daughters[0]->GetPDGMass();

        // Q value was calculated from atomic masses.
        // Use it to get correct neutron energy.
        G4double cmMomentum = std::sqrt(transitionQ*(transitionQ + 2.*alphaMass)*
                                        (transitionQ + 2.*nucleusMass)*
                                        (transitionQ + 2.*alphaMass + 2.*nucleusMass) )/
                              (transitionQ + alphaMass + nucleusMass)/2.;

        // Set up final state
        // parentParticle is set at rest here because boost with correct momentum
        // is done later
        G4DynamicParticle parentParticle(G4MT_parent, G4ThreeVector(0,0,0), 0.0);
        products = new G4DecayProducts(parentParticle);
        G4double costheta = 2.*G4UniformRand()-1.0;
        G4double sintheta = std::sqrt(1.0 - costheta*costheta);
        G4double phi  = twopi*G4UniformRand()*rad;
        G4ThreeVector direction(sintheta*std::cos(phi),sintheta*std::sin(phi),
                                costheta);

        G4double KE = std::sqrt(cmMomentum*cmMomentum + alphaMass*alphaMass)
                      - alphaMass;
        G4DynamicParticle* daughterparticle =
                new G4DynamicParticle(G4MT_daughters[1], direction, KE, alphaMass);
        products->PushProducts(daughterparticle);

        KE = std::sqrt(cmMomentum*cmMomentum + nucleusMass*nucleusMass) - nucleusMass;
        daughterparticle =
                new G4DynamicParticle(G4MT_daughters[0], -1.0*direction, KE, nucleusMass);
        products->PushProducts(daughterparticle);
    }
    return products;
}


void G4AlphaWidthDecay::DumpNuclearInfo()
{
    G4cout << " G4AlphaDecay for parent nucleus " << GetParentName() << G4endl;
    G4cout << " decays to " << GetDaughterName(0) << " + " << GetDaughterName(1)
           << " with branching ratio " << GetBR() << "% and Q value "
           << transitionQ << G4endl;
}

