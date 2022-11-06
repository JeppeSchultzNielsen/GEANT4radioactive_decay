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
//  File:   G4NuclearDecay.cc                                                 //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   11 December 2014                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4NuclearDecay.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>
#include "Randomize.hh"

G4NuclearDecay::G4NuclearDecay(const G4String& channelName,
                               const G4RadioactiveDecayMode& aMode,
                               const G4double& excitationE,
                               const G4Ions::G4FloatLevelBase& flb)
 : G4VDecayChannel(channelName), theMode(aMode), daughterEx(excitationE),
   floatingLevel(flb), halflifeThreshold(nanosecond) 
{}

G4NuclearDecay::~G4NuclearDecay()
{}

//helper method for ChooseDecaySublevel and ChooseDalitzConf
G4int G4NuclearDecay::BinarySearch(std::vector<G4double> *cumProbs, G4double x, G4int low, G4int high){
    if(low > high){return -1;} //level was not found
    else{
        G4int mid = (low + high)/2;
        if(mid == 0){
            if(x < cumProbs->at(1)) return mid; //should have reached bottom; hopefully x is smaller than 1st entry (0th is 0)
            else{return -1;}
        }
        if(low == high){ //in case we reached top
            if(x > cumProbs->at(high) && x <= cumProbs->at(high)+1) return mid;
            else{return -1;}
        }
        if(x > cumProbs->at(mid)){
            if(x <= cumProbs->at(mid+1)) return mid; //found the level
            else return BinarySearch(cumProbs,x,mid+1,high); //x is in upper half of the array
        }
        else return BinarySearch(cumProbs,x,low,mid-1); //x is in lower half of the array
    }
}

G4int G4NuclearDecay::ChooseDecaySublevel(){
    //in case probability distribution is not properly normalized
    double probRoof = probabilitySum * G4UniformRand();
    double probSum = 0;

    //finding by weighted average. Can be optimized by look-up in binary tree somehow?
    for(int i = 0; i < sublevelBRs.size(); i++){
        probSum += sublevelBRs[i];
        if(probSum > probRoof){
            return i;
        }
    }

    //find by binary search.
    //return binarySearch(&sublevelCumBRs,probRoof,0,sublevelCumBRs.size()-1);

    //error has occured, as probSum never crossed probRoof (should be impossible)
    return -1;
}

G4bool G4NuclearDecay::ReadWidthFile(G4int daughterZ, G4int daughterA, G4double nominalDaughterEx, G4double nominalQvalue){
    //initialize vectors for storing output
    sublevelCumBRs = {};
    sublevelBRs = {};
    sublevelExs = {};
    sublevelQvalues = {};
    probabilitySum = 0;

    //nominalDaughterEx is apparently in MeV; i need it in keV:
    nominalDaughterEx = 1000 * nominalDaughterEx;

    //first check if the daughterNominalEx exists in the decayfile. I don't know how G4 does this, but i can find the
    //correct directory by doing this:
    G4String radDirPath = "";
    char* path_var = std::getenv("G4RADIOACTIVEDATA");
    if (!path_var) {
        G4Exception("G4RadioactiveDecay()","HAD_RDM_200",FatalException,
                    "Environment variable G4RADIOACTIVEDATA is not set");
    } else {
        radDirPath = path_var;   // convert to string
    }
    std::ostringstream os;
    os << radDirPath << "/Widths/z" << daughterZ << "WIDTH.a" << daughterA << '\0';

    char inputChars[120]={' '};
    G4String inputLine;
    G4bool found(false);
    G4String dump;
    G4double nominalRead;
    G4double sublevelBR;
    G4double sublevel;

    std::ifstream file;
    file.open(os.str());
    if (!file.is_open() ){
        return false;
    }
    else {
        //read file to see if nominalDaughter ex is present:
        while (!file.getline(inputChars, 120).eof()) {
            inputLine = inputChars;
            std::istringstream tmpStream(inputLine);
            if(inputChars[0] == 'N'){
                //if the correct level has already been read, just return; at this point the values have been read.
                if(found) return true;
                //else
                else{
                    tmpStream >> dump >> nominalRead;
                    if(abs(nominalRead-nominalDaughterEx) < 0.01){
                        //if the read nominallevel is less than 10eV from the daughterlevel, we consider it found.
                        found = true;
                    }
                }
            }
            else{
                //when the line does not start with N, we are reading sublevels and sublevelBRs. If the nominalLevel is
                //found, these should be loaded into the vectors.
                if(found){
                    //if Q value is less than 0, sharply set BR to 0 - this should be implemented in a smoother way
                    //later
                    if((nominalQvalue*1000 - (sublevel - nominalDaughterEx))/1000 < 0){
                        sublevelBR = 0;
                    }

                    tmpStream >> sublevel >> sublevelBR;
                    sublevelExs.push_back(sublevel);
                    sublevelBRs.push_back(sublevelBR);
                    probabilitySum += sublevelBR;
                    sublevelCumBRs.push_back(probabilitySum);
                    //Q-value of decay to the sublevel is the nominal Q-value minus the relative excitation of the
                    //sublevelnominalQvalue*1000
                    sublevelQvalues.push_back((nominalQvalue*1000 - (sublevel - nominalDaughterEx))/1000);
                }
            }
        }
        //file has been read. If the nominalLevel wasnt found, return false. If it was found, it was the last level,
        //return true.
        return found;
    }
}