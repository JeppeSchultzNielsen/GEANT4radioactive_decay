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
//ANT decay controls 3-body breakup of He-8.

#include "G4ANTDecay.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <iomanip>

G4ANTDecay::G4ANTDecay(const G4ParticleDefinition* theParentNucleus,
                             const G4double& branch, const G4double& Qvalue,
                             const G4double& excitationE,
                             const G4Ions::G4FloatLevelBase& flb, DalitzHandler *newDalitzHandler)
        : G4NuclearDecay("ANT decay", ANT, excitationE, flb), transitionQ(Qvalue), dalitzHandler(newDalitzHandler)
{
    SetParent(theParentNucleus);  // Store name of parent nucleus, delete G4MT_parent
    SetBR(branch);

    SetNumberOfDaughters(3);
    SetDaughter(0, "neutron");
    SetDaughter(1, "triton");
    SetDaughter(2, "alpha");
}


G4ANTDecay::~G4ANTDecay()
{}


G4DecayProducts* G4ANTDecay::DecayIt(G4double)
{
    // Fill G4MT_parent with theParentNucleus (stored by SetParent in ctor)
    CheckAndFillParent();

    // Fill G4MT_daughters with triton and residual nucleus (stored by SetDaughter)
    CheckAndFillDaughters();

    //naming convention stems from file in which DalitzPlot was defined; 1 is neutron, 2 is triton, 3 is alpha.
    G4double m1 = G4MT_daughters[0]->GetPDGMass();
    G4double m2 = G4MT_daughters[1]->GetPDGMass();
    G4double m3 = G4MT_daughters[2]->GetPDGMass();
    G4double M = G4MT_parent -> GetPDGMass();
    G4double M = G4MT_parent -> GetPDGMass();

    vector<G4double> s1s2 = dalitzHandler -> GetCorrectedS1S2(M);
    G4double s1 = s1s2[0];
    G4double s2 = s1s2[1];
    G4double s3 =

    //kinematics are fully contained within s1 and s2 parameters.

    // Q value was calculated from atomic masses.
    // Use it to get correct triton energy.
    G4double cmMomentum = std::sqrt(transitionQ*(transitionQ + 2.*tritonMass)*
                                    (transitionQ + 2.*nucleusMass)*
                                    (transitionQ + 2.*tritonMass + 2.*nucleusMass) )/
                          (transitionQ + tritonMass + nucleusMass)/2.;

    // Set up final state
    // parentParticle is set at rest here because boost with correct momentum
    // is done later
    G4DynamicParticle parentParticle(G4MT_parent, G4ThreeVector(0,0,0), 0.0);
    G4DecayProducts* products = new G4DecayProducts(parentParticle);

    G4double costheta = 2.*G4UniformRand()-1.0;
    G4double sintheta = std::sqrt(1.0 - costheta*costheta);
    G4double phi  = twopi*G4UniformRand()*rad;
    G4ThreeVector direction(sintheta*std::cos(phi),sintheta*std::sin(phi),
                            costheta);

    G4double KE = std::sqrt(cmMomentum*cmMomentum + tritonMass*tritonMass)
                  - tritonMass;
    G4DynamicParticle* daughterparticle =
            new G4DynamicParticle(G4MT_daughters[1], direction, KE, tritonMass);
    products->PushProducts(daughterparticle);

    KE = std::sqrt(cmMomentum*cmMomentum + nucleusMass*nucleusMass) - nucleusMass;
    daughterparticle =
            new G4DynamicParticle(G4MT_daughters[0], -1.0*direction, KE, nucleusMass);
    products->PushProducts(daughterparticle);

    return products;
}


void G4ANTDecay::DumpNuclearInfo()
{
    G4cout << " G4TritonDecay for parent nucleus " << GetParentName() << G4endl;
    G4cout << " decays to " << GetDaughterName(0) << " + " << GetDaughterName(1) << " + " << GetDaughterName(2)
           << " with branching ratio " << GetBR() << "% and Q value "
           << transitionQ << G4endl;
}

