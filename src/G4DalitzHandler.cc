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
    probabilitySum = 0.;
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
    G4cout << std::setprecision (std::numeric_limits<double>::digits10 + 1);
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
                            /*if(prob > 1){
                                G4cout << prob << G4endl;
                                G4cout << probabilitySum << G4endl;
                            }*/
                            dalitzProbs.push_back(probabilitySum);
                            currentConf = {s1no,s2no};
                            dalitzConfs.push_back(currentConf);
                        }
                        s1no++;
                    }
                }
            }
        }
        G4cout << std::setprecision (std::numeric_limits<double>::digits10 + 1);
        /*for(int i = 0; i < dalitzProbs.size(); i++){
            G4cout << dalitzProbs[i]+0.1 << G4endl;
        }*/
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
            if(x > cumProbs->at(high) && x <= cumProbs->at(high+1)) return mid;
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
    G4double probRoof  = probabilitySum *1.* G4UniformRand();
    G4int chosenIndex = BinarySearch(&dalitzProbs,probRoof,0,dalitzProbs.size()-1);
    return dalitzConfs[chosenIndex];
}

std::vector<G4double> G4DalitzHandler::GetCorrectedS1S2(G4double sublevelMass, G4double m1, G4double m2, G4double m3){

    std::vector<G4int> chosenIndeces = ChooseDalitzConfiguration();
    G4double s1i = s1s[chosenIndeces[0]]/1000000; //change units from keV to MeV
    G4double s2i = s2s[chosenIndeces[1]]/1000000; //change units from keV to MeV
    G4double si = nomMass*nomMass;
    G4double s3i = si + m1*m1 + m2*m2 + m3*m3 - s1i - s2i;

    /*G4cout << s1i << G4endl;
    G4cout << s2i << G4endl;*/
    G4double sf = sublevelMass*sublevelMass;
    G4double ds = sf-si;

    //G4double s1f = s1i + ds*s1i/(si + m1*m1 + m2*m2 + m3*m3);
    //G4double s2f = s2i + ds*s2i/(si + m1*m1 + m2*m2 + m3*m3);

    G4double s1weight = s1i*G4UniformRand();
    G4double s2weight = s2i*G4UniformRand();
    G4double s3weight = s3i*G4UniformRand();

    s1weight *= 1/(s1weight + s2weight + s3weight);
    s2weight *= 1/(s1weight + s2weight + s3weight);
    s3weight *= 1/(s1weight + s2weight + s3weight);

    G4double s1f = s1i + ds * s1weight;
    G4double s2f = s2i + ds * s2weight;
    G4double s3f = sf + m1*m1 + m2*m2 + m3*m3 - s1f - s2f;

    G4double upBs1 = UpperBoundS1(sf,s1f,s2f,m1,m2,m3);
    G4double lowBs1 = LowerBoundS1(sf,s1f,s2f,m1,m2,m3);
    G4double upBs3 = std::pow((std::sqrt(sf)-m2),2);
    G4double lowBs3 = (m1+m3)*(m1+m3);
    G4double lowBs2 = (m2+m3)*(m2+m3);
    G4double upBs2 = (std::sqrt(sf)-m1)*(std::sqrt(sf)-m1);


    int i = 0;
    //try new random numbers until we find some that work.
    while(s1f > upBs1 || s1f < lowBs1 || s3f > upBs3 || s3f < lowBs3 || std::isnan(upBs1) || std::isnan(lowBs1)){
        s1weight = s1i*G4UniformRand();
        s2weight = s2i*G4UniformRand();
        s3weight = s3i*G4UniformRand();
        s1weight *= 1/(s1weight + s2weight + s3weight);
        s2weight *= 1/(s1weight + s2weight + s3weight);
        s3weight *= 1/(s1weight + s2weight + s3weight);
        s1f = s1i + ds * s1weight;
        s2f = s2i + ds * s2weight;
        s3f = sf + m1*m1 + m2*m2 + m3*m3 - s1f - s2f;

        upBs1 = UpperBoundS1(sf,s1f,s2f,m1,m2,m3);
        lowBs1 = LowerBoundS1(sf,s1f,s2f,m1,m2,m3);
        i++;
        if(i > 20){
            //at edges of Dalitz plot, s1 and s2 can be at their minimum values (can happen for any pair of invariant
            //masses). At this point, the algorithm for choosing new values will not progress, because two of the weights
            //would have to be zero. If the algorithm doesnt progress, we should therefore choose the most recent value
            //and force it to be legal.
            //or really just pull a new value
            chosenIndeces = ChooseDalitzConfiguration();
            s1i = s1s[chosenIndeces[0]]/1000000; //change units from keV to MeV
            s2i = s2s[chosenIndeces[1]]/1000000; //change units from keV to MeV
            G4cout << i << G4endl;
            i = 0;
            /*
            if(s2f < lowBs2){
                G4double ds = ds - (s2f - lowBs2);
                s2f = lowBs2;
                //as default, share remaining ds among s1 and s3 in weighted manner.
                s1f = s1i + ds * s1i /(s1i+s3i);
                s3f = s3i + ds * s3i /(s1i+s3i);
                upBs1 = UpperBoundS1(sf,s1f,s2f,m1,m2,m3);
                lowBs1 = LowerBoundS1(sf,s1f,s2f,m1,m2,m3);
                if(s1f < lowBs1){
                    s1f = lowBs1;
                    s3f = sf + m1*m1 + m2*m2 + m3*m3 - s1f - s2f;
                }
                if(s1f > upBs1){
                    s1f = upBs1;
                    s3f = sf + m1*m1 + m2*m2 + m3*m3 - s1f - s2f;
                }
                if(s3f < lowBs3){
                    s3f = lowBs3;
                    s1f = sf + m1*m1 + m2*m2 + m3*m3 - s3f - s2f;
                }
                if(s3f > upBs3){
                    s3f = upBs3;
                    s1f = sf + m1*m1 + m2*m2 + m3*m3 - s3f - s2f;
                }
            }
            if(s2f > upBs2){
                G4double ds = ds - (s2f - lowBs2);
                s2f = upBs2;
                //as default, share remaining ds among s1 and s3 in weighted manner.
                s1f = s1i + ds * s1i /(s1i+s3i);
                s3f = s3i + ds * s3i /(s1i+s3i);
                upBs1 = UpperBoundS1(sf,s1f,s2f,m1,m2,m3);
                lowBs1 = LowerBoundS1(sf,s1f,s2f,m1,m2,m3);
                if(s1f < lowBs1){
                    s1f = lowBs1;
                    s3f = sf + m1*m1 + m2*m2 + m3*m3 - s1f - s2f;
                }
                if(s1f > upBs1){
                    s1f = upBs1;
                    s3f = sf + m1*m1 + m2*m2 + m3*m3 - s1f - s2f;
                }
                if(s3f < lowBs3){
                    s3f = lowBs3;
                    s1f = sf + m1*m1 + m2*m2 + m3*m3 - s3f - s2f;
                }
                if(s3f > upBs3){
                    s3f = upBs3;
                    s1f = sf + m1*m1 + m2*m2 + m3*m3 - s3f - s2f;
                }
            }
            else{
                //s2f is legal. now choose s1f so that its also legal.
                if(s1f < lowBs1){
                    s1f = lowBs1;
                    s3f = sf + m1*m1 + m2*m2 + m3*m3 - s1f - s2f;
                }
                if(s1f > upBs1){
                    s1f = upBs1;
                    s3f = sf + m1*m1 + m2*m2 + m3*m3 - s1f - s2f;
                }
            }*/
        }
    }

    return {s1f,s2f};
}

G4double G4DalitzHandler::UpperBoundS1(G4double s, G4double s1, G4double s2, G4double m1, G4double m2, G4double m3){
    return m1*m1 + m2*m2 - 1/(2*s2)*((s2 - s + m1*m1)*(s2+m2*m2-m3*m3)-std::sqrt(tri(s2,s,m1*m1)*tri(s2,m2*m2,m3*m3)));
};
G4double G4DalitzHandler::LowerBoundS1(G4double s, G4double s1, G4double s2, G4double m1, G4double m2, G4double m3){
    return m1*m1 + m2*m2 - 1/(2*s2)*((s2 - s + m1*m1)*(s2+m2*m2-m3*m3)+std::sqrt(tri(s2,s,m1*m1)*tri(s2,m2*m2,m3*m3)));
}