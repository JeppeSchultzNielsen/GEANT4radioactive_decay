//
// Created by jeppe on 11/6/22.
//
// class for handling loading of Dalitz plots and giving s1 and s2 values.
#ifndef G4DalitzHandler_h
#define G4DalitzHandler_h 1

#include "vector"
#include "globals.hh"

class G4DalitzHandler{
public:
    G4DalitzHandler(const G4int daughterZ, const G4int daughterA, const G4double nomEx, const G4double nomMass);
    ~G4DalitzHandler();

    G4double GetNomMass(){return nomMass;};

    std::vector<G4double> GetCorrectedS1S2(G4double sublevelMass);
private:
    //called in construction
    G4bool ReadDalitzFile(G4int daughterZ, G4int daughterA, G4double nominalParentEx);

    //helper function for making weighted random choice of dalitzconfiguration
    G4int BinarySearch(std::vector<G4double> *cumProbs, G4double x, G4int low, G4int high);

    //chooses a dalitz configuration
    std::vector<G4int> ChooseDalitzConfiguration();

    G4double probabilitySum;
    G4double nomMass;
    G4bool readSucces;

    //data for Dalitz configurations
    std::vector<G4double> s1s;
    std::vector<G4double> s2s;
    std::vector<G4double> dalitzProbs;
    std::vector<std::vector<G4int>> dalitzConfs;
};

#endif