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
                             const G4Ions::G4FloatLevelBase& flb, G4DalitzHandler *newDalitzHandler)
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
    G4double m1 = G4MT_daughters[0]->GetPDGMass() * 1000;
    G4double m2 = G4MT_daughters[1]->GetPDGMass() * 1000;
    G4double m3 = G4MT_daughters[2]->GetPDGMass() * 1000;
    G4double M = G4MT_parent -> GetPDGMass() * 1000;
    G4double nomM = dalitzHandler -> GetNomMass() * 1000;

    //for now, all Decays are treated like nominal decays. For this reason s3 is calculated from nominal mass.
    std::vector<G4double> s1s2 = dalitzHandler -> GetCorrectedS1S2(M);
    G4double s = nomM*nomM;
    G4double s1 = s1s2[0];
    G4double s2 = s1s2[1];
    G4double s3 = s+ m1*m1 + m2*m2 + m3*m3 - s1 - s2;

    //kinematics are fully contained within s1 and s2 parameters (which also constrain s3). Reference to Byckling Particle Kinematics.
    //first, place neutron along z-axis.
    G4double neutronMomentum = std::sqrt(tri(s,m1*m1,s2)/s)/2.*keV;
    G4ThreeVector neutronDirection(0, 0, 1);
    G4ThreeVector neutronVectorMomentum = neutronMomentum*neutronDirection;

    //theta12 is given by kinematics. When neutron momentum is placed along z-axis, there is cylindrical symmetry and phi
    //of triton can be chosen arbitrarily. Only relative angles matter at this point; reaction plane will be rotated
    //randomly at later point.
    G4double costheta12 = ( (s+m1*m1-s2)*(s+m2*m2-s3) + 2*s*(m1*m1 + m2*m2 - s1) )/std::sqrt(tri(s,m1*m1,s2)*tri(s,m2*m2,s3));
    G4double tritonMomentum = std::sqrt(tri(s,m2*m2,s3)/s)/2.*keV;
    G4ThreeVector tritonDirection(std::sin(std::acos(costheta12)), 0, costheta12);
    G4ThreeVector tritonVectorMomentum = tritonMomentum*tritonDirection;

    //by conservation of momentum, the momentum of the alpha particle is now decided:
    G4ThreeVector alphaVectorMomentum = - neutronVectorMomentum - tritonVectorMomentum;

    //now, rotate the three vectors by random rotation to ensure isotropic radiation.
    //phi must be chosen randomly between 0 and 2pi. sin(theta) must be chosen between 0 and pi for isotropic distb.
    //ThreeVector's rotate function is given wrt angle and rotation normal vector. Find these through https://en.wikipedia.org/wiki/Rotation_matrix
    G4double alpha = std::asin(2.*G4UniformRand()-1.0);
    G4double beta = twopi*G4UniformRand()*rad;
    G4double gamma = twopi*G4UniformRand()*rad;

    G4ThreeVector rotationVector(std::sin(alpha)*std::cos(beta) - std::cos(alpha)*std::sin(beta)*std::sin(gamma)+
        std::sin(alpha)*std::cos(gamma),
        std::cos(alpha)*std::sin(beta)*std::cos(gamma)+std::sin(alpha)*std::sin(gamma)+std::sin(beta),
        std::cos(beta)*std::sin(gamma)-std::sin(alpha)*std::sin(beta)*std::sin(gamma)+std::cos(alpha)*std::sin(gamma));

    G4double rotationAngle = std::acos( (std::cos(beta)*std::cos(gamma)+std::sin(alpha)*std::sin(beta)*std::sin(gamma)
            + std::cos(alpha)*std::cos(gamma)+std::cos(alpha)*cos(beta)-1.)/2. );

    //rotate all vectors.
    neutronVectorMomentum.rotate(rotationAngle,rotationVector);
    tritonVectorMomentum.rotate(rotationAngle,rotationVector);
    alphaVectorMomentum.rotate(rotationAngle,rotationVector);

    G4DynamicParticle parentParticle(G4MT_parent, G4ThreeVector(0,0,0), 0.0);
    G4DecayProducts* products = new G4DecayProducts(parentParticle);

    G4DynamicParticle* dynamicNeutron
            = new G4DynamicParticle(G4MT_daughters[0], neutronVectorMomentum);
    products->PushProducts(dynamicNeutron);

    G4DynamicParticle* dynamicTriton
            = new G4DynamicParticle(G4MT_daughters[1], tritonVectorMomentum);
    products->PushProducts(dynamicTriton);

    G4DynamicParticle* dynamicAlpha
            = new G4DynamicParticle(G4MT_daughters[2], alphaVectorMomentum);
    products->PushProducts(dynamicAlpha);

    return products;
}


void G4ANTDecay::DumpNuclearInfo()
{
    G4cout << " G4TritonDecay for parent nucleus " << GetParentName() << G4endl;
    G4cout << " decays to " << GetDaughterName(0) << " + " << GetDaughterName(1) << " + " << GetDaughterName(2)
           << " with branching ratio " << GetBR() << "% and Q value "
           << transitionQ << G4endl;
}

