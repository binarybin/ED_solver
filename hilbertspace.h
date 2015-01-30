//
//  hilbertspace.h
//  ED_solver
//
//  Created by Bin Xu on 12/19/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#ifndef __ED_solver__hilbertspace__
#define __ED_solver__hilbertspace__

#include <iostream>
#include <unordered_map>
#include <vector>
#include "hamiltonian.h"
#include "miscellaneous.h"
using namespace std;

class Interaction;
class Solver;
class Measurement;

class HilbertSpace
{
    friend class Interaction;
    friend class Solver;
    friend class Measurement;
    Hamiltonian &ham;
    vector<Orbital> orbital_list;
    unordered_map<Orbital, size_t, Orbital_hasher> reverse_orbital_map;
    unordered_map<Orbital, vector<OrbPair>, Orbital_hasher> pair_map;
    unordered_map<Orbital, vector<CompactState>, Orbital_hasher> state_map;
    
public:
    HilbertSpace(Hamiltonian &o_ham): ham(o_ham){}
    void buildOrbitalList();
    void buildReverseOrbitalMap();
    void buildPairMap();
    void buildStateMap();
    unordered_map<Orbital, vector<CompactState>, Orbital_hasher>& getStateMap() {return state_map;}
    
private: // some auxiliary functions
    void generateHaldaneOrbList();
    void computeFacHaldane();
    double last_discarded_r2;
    
};

Pxy computeMomentum(CompactState &cstate, vector<Orbital>& orbital_list);
Orbital computeMomentumAndSpin(CompactState &cstate, vector<Orbital>& orbital_list);
void csPlusplus(CompactState& cstate, int norb);
bool csNotHighest(CompactState& cstate, int norb, int nele);

#endif /* defined(__ED_solver__hilbertspace__) */
