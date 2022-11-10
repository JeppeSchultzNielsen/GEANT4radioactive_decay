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

std::vector<G4double> G4DalitzHandler::GetCorrectedS1S2(G4double sublevelMass, G4double m1, G4double m2, G4double m3){

    std::vector<G4int> chosenIndeces = ChooseDalitzConfiguration();
    G4double s1i = s1s[chosenIndeces[0]];
    G4double s2i = s2s[chosenIndeces[1]];
    G4double si = nomMass*nomMass;
    G4double s3i = si + m1*m1 + m2*m2 + m3*m3 - s1i - s2i;

    G4double sf = sublevelMass*sublevelMass;
    G4double ds = sf-si;

    /*G4cout << "nomMass " << nomMass << G4endl;
    G4cout << "subMass " << sublevelMass << G4endl;
    G4cout << "ds " << ds << G4endl;
    G4cout << "s1i " << s1i << G4endl;
    G4cout << "s2i " << s2i << G4endl;*/

    //since s = s1+s2+s3 - m1² - m2² - m3², sf = si+ds = ds+s1+s2+s3 - m1² - m2² - m3², preferred situation is that
    //new s1, s2 and s3 are chosen such that ds is distributed in a weighted way. The hope is that this would correspond
    //to rescaling the dalitz plot without changing it significantly. Need a more correct way to do this.

    G4double s1f = s1i + ds*s1i/(s1i+s2i+s3i);
    G4double s2f = s2i + ds*s2i/(s1i+s2i+s3i);
    G4double s3f = s3i + ds*s3i/(s1i+s2i+s3i);

    //check if they obey their lower bounds. We can only enter these situations if ds is negative.
    if(ds < 0){
        if(s1f < (m1+m2)*(m1+m2)){
            if(s2f < (m2+m3)*(m2+m3)){
                //both s1 and s2 crossed their lower bounds in the correction. Set them equal to their lower bounds and s3
                //takes it all.
                s1f = (m1+m2)*(m1+m2);
                s2f = (m2+m3)*(m2+m3);
                s3f = sf + m1*m1 + m2*m2 + m3*m3 - s1f - s2f;
                //this must be an allowed configuration.
                return {s1f,s2f};
            }
            if(s3f < (m1+m3)*(m1+m3)){
                //both s1 and s3 crossed their lower bounds in the correction. Set them equal to their lower bounds and s2
                //takes it all.
                s1f = (m1+m2)*(m1+m2);
                s3f = (m1+m3)*(m1+m3);
                s2f = sf + m1*m1 + m2*m2 + m3*m3 - s1f - s3f;
                //this must be an allowed configuration.
                return {s1f,s2f};
            }
            //only s1f is below lower bound. Set it to lower bound and recalculate s2f and s3f.
            s1f = (m1+m2)*(m1+m2);
            //I've already used som of ds to correct s1
            s2f = s2i + (ds-(s1f-s1i))*s2i/(s2i+s3i);
            s3f = s3i + (ds-(s1f-s1i))*s3i/(s2i+s3i);
            //check if upper bounds are now obeyed - maybe this is not neceasary?
            if(s2f > std::pow((sublevelMass - m1),2)){
                G4cout << "s1 was at lowerbound, checking s2s upper bound" << G4endl;
                //set s2f to upper limit and recalculate s3f:
                s2f = std::pow(sublevelMass - m1,2);
                s3f = sf + m1*m1 + m2*m2 + m3*m3 - s1f - s2f;
                return {s1f,s2f};
            }
            if(s3f > std::pow(sublevelMass - m2,2)){
                //set s3f to upper limit and recalculate s2f:
                s3f = std::pow(sublevelMass - m2,2);
                s2f = sf + m1*m1 + m2*m2 + m3*m3 - s1f - s3f;
                return {s1f,s2f};
            }
            //all limits are upheld
            return {s1f,s2f};
        }
        if(s2f < (m3+m2)*(m3+m2)){
            if(s3f < (m1+m3)*(m1+m3)){
                //both s2 and s3 crossed their lower bounds in the correction. Set them equal to their lower bounds and s1
                //takes it all.
                s2f = (m3+m2)*(m3+m2);
                s3f = (m1+m3)*(m1+m3);
                s1f = sf + m1*m1 + m2*m2 + m3*m3 - s2f - s3f;
                //this must be an allowed configuration.
                return {s1f,s2f};
            }
            //only s2f is below lower bound. Set it to lower bound and recalculate s1f and s3f.
            s2f = (m3+m2)*(m3+m2);
            //I've already used som of ds to correct s2
            s1f = s1i + (ds-(s2f-s2i))*s1i/(s1i+s3i);
            s3f = s3i + (ds-(s2f-s2i))*s3i/(s1i+s3i);
            //check if upper bounds are now obeyed
            if(s1f > std::pow(sublevelMass - m3,2)){
                //set s1f to upper limit and recalculate s3f:
                s1f = std::pow(sublevelMass - m3,2);
                s3f = sf + m1*m1 + m2*m2 + m3*m3 - s1f - s2f;
                return {s1f,s2f};
            }
            if(s3f > std::pow(sublevelMass - m2,2)){
                //set s3f to upper limit and recalculate s1f:
                s3f = std::pow(sublevelMass - m2,2);
                s1f = sf + m1*m1 + m2*m2 + m3*m3 - s2f - s3f;
                return {s1f,s2f};
            }
            //all limits are upheld
            return {s1f,s2f};
        }
        if(s3f < (m3+m1)*(m3+m1)){
            //only s2f is below lower bound. Set it to lower bound and recalculate s1f and s2f.
            s3f = (m3+m1)*(m3+m1);
            //I've already used som of ds to correct s3
            s1f = s1i + (ds-(s3f-s3i))*s1i/(s1i+s2i);
            s2f = s2i + (ds-(s3f-s3i))*s2i/(s1i+s2i);
            //check if upper bounds are now obeyed
            if(s1f > std::pow(sublevelMass - m3,2)){
                //set s1f to upper limit and recalculate s2f:
                s1f = std::pow(sublevelMass - m3,2);
                s2f = sf + m1*m1 + m2*m2 + m3*m3 - s3f - s1f;
                return {s1f,s2f};
            }
            if(s2f > std::pow(sublevelMass - m1,2)){
                //set s2f to upper limit and recalculate s1f:
                s2f = std::pow(sublevelMass - m1,2);
                s1f = sf + m1*m1 + m2*m2 + m3*m3 - s2f - s3f;
                return {s1f,s2f};
            }
            //all limits are upheld
            return {s1f,s2f};
        }
        return {s1f,s2f};
    }

    //case ds > 0. Do opposite of what i did before.
    if(ds > 0){
        if(s1f > std::pow(sublevelMass - m3,2)){
            if(s2f > std::pow(sublevelMass - m1,2)){
                s1f = std::pow(sublevelMass - m3,2);
                s2f = std::pow(sublevelMass - m1,2);
                return {s1f,s2f};
            }
            if(s3f > std::pow(sublevelMass - m2,2)){
                s1f = std::pow(sublevelMass - m3,2);
                s3f = std::pow(sublevelMass - m2,2);
                s2f = s2f = sf + m1*m1 + m2*m2 + m3*m3 - s3f - s1f;
                return {s1f,s2f};
            }
            //only s1f is above upper bound. Set it to upper bound and recalculate s2f and s3f.
            s1f = (m1+m2)*(m1+m2);
            //I've already used som of ds to correct s1
            s2f = s2i + (ds-(s1f-s1i))*s2i/(s2i+s3i);
            s3f = s3i + (ds-(s1f-s1i))*s3i/(s2i+s3i);
            return {s1f, s2f}
        }
        if(s2f > std::pow(sublevelMass - m1,2)){

        }
        if(s3f > std::pow(sublevelMass - m2,2)){

        }
        //all limits are upheld
        return {s1f,s2f};
    }

    return {s1i,s2i};
}