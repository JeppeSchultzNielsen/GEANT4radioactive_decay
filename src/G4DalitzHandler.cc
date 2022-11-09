//
// Created by jeppe on 11/6/22.
//

#include "G4DalitzHandler.hh"
#include "Randomize.hh"
#include <fstream>

G4DalitzHandler::G4DalitzHandler(const G4int daughterZ, const G4int daughterA, const G4double nomEx, const G4double nomMass){
    readSucces = ReadDalitzFile(daughterZ, daughterA, nomEx);
    this -> nomMass = nomMass;
}

G4DalitzHandler::~G4DalitzHandler()
{
    s1s.clear();
    s2s.clear();
    dalitzConfs.clear();
    dalitzProbs.clear();
}

G4bool G4DalitzHandler::ReadDalitzFile(G4int daughterZ, G4int daughterA, G4double nominalParentEx){
    //initialize vectors for storing output
    s1s = {};
    s2s = {};
    dalitzConfs = {};
    dalitzProbs = {};
    probabilitySum = 0;
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
    os << radDirPath << "/DalitzFiles/z" << daughterZ << "DALITZ.a" << daughterA << '\0';

    char inputChars[25000]={' '};
    G4String inputLine;
    G4bool found(false);
    G4String dump;
    G4double nominalRead;
    G4double s1;
    G4double s2;
    G4double prob;
    std::vector<G4int> currentConf;

    G4bool hasNext = true;

    std::ifstream file;
    file.open(os.str());
    if (!file.is_open() ){
        G4String error ="Did not find DalitzFile for A " + std::to_string(daughterA) + " Z " + std::to_string(daughterZ);
        G4Exception("G4DalitzHandler()","HAD_RDM_200",FatalException,
                    error);
        return false;
    }
    else {
        //read file to see if nominalParentEx is present:
        while (!file.getline(inputChars, 25000).eof()) {
            inputLine = inputChars;
            std::istringstream tmpStream(inputLine);
            if(inputChars[0] == 'N'){
                //if the correct level has already been read, just return; at this point the values have been read.
                if(found) return true;
                    //else
                else{
                    tmpStream >> dump >> nominalRead;
                    if(abs(nominalRead/1000.-nominalParentEx) < 0.01){
                        //if the read nominallevel is less than 10eV from the read level, we consider it found.
                        found = true;
                        //the entries following the nominallevel are the s1 values.
                        while(tmpStream >> s1) {
                            s1s.push_back(s1);
                        }
                    }
                }
            }
            else{
                //when the line does not start with N, we are reading s2-values and probabilities. If the correct parentlevel
                //found, these should be loaded into the vectors.
                if(found){
                    tmpStream >> s2;
                    s2s.push_back(s2);
                    G4int s1no = 0;
                    G4int s2no = s2s.size();
                    while(tmpStream >> prob){
                        //now the probabilites come
                        if(prob > 0){
                            probabilitySum += prob;
                            dalitzProbs.push_back(probabilitySum);
                            currentConf = {s1no,s2no};
                            dalitzConfs.push_back(currentConf);
                        }
                        s1no++;
                    }
                }
            }
        }
        //file has been read. If the nominalLevel wasnt found, return false. If it was found, it was the last level,
        //return true.
        return found;
    }
}

G4int G4DalitzHandler::BinarySearch(std::vector<G4double> *cumProbs, G4double x, G4int low, G4int high){
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

std::vector<G4int> G4DalitzHandler::ChooseDalitzConfiguration(){
    G4double probRoof  = probabilitySum * G4UniformRand();
    return dalitzConfs[BinarySearch(&dalitzProbs,probRoof,0,dalitzProbs.size()-1)];
}

std::vector<G4double> G4DalitzHandler::GetCorrectedS1S2(G4double sublevelMass){
    //for now, just return uncorrected S1s and S2s.
    std::vector<G4int> chosenIndeces = ChooseDalitzConfiguration();
    /*G4String radDirPath = "";
    char* path_var = std::getenv("G4RADIOACTIVEDATA");
    if (!path_var) {
        G4Exception("G4RadioactiveDecay()","HAD_RDM_200",FatalException,
                    "Environment variable G4RADIOACTIVEDATA is not set");
    } else {
        radDirPath = path_var;   // convert to string
    }

    G4String path = radDirPath + "/DalitzFiles/check.txt";

    std::ofstream outfile;

    outfile.open(path, std::ios_base::app); // append instead of overwrite
    outfile << s1s[chosenIndeces[0]] << "\t" << s2s[chosenIndeces[1]] << "\n";

    outfile.close();*/
    return {s1s[chosenIndeces[0]],s2s[chosenIndeces[1]]};
}