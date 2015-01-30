//
//  interaction.h
//  ED_solver
//
//  Created by Bin Xu on 12/19/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#ifndef __ED_solver__interaction__
#define __ED_solver__interaction__

#include "hilbertspace.h"

class Solver;
class Measurement;

class Interaction
{
    friend class Solver;
    friend class Measurement;
    vector<State> state_list;
    unordered_map<unsigned long, unsigned> reverse_state_map;
    
    HilbertSpace &hilbert_space;
    Orbital orb_sector;
public:
    Interaction(HilbertSpace &o_hilbert_space, Orbital o_orb): hilbert_space(o_hilbert_space), orb_sector(o_orb){}
    vector<MatEle> matrix;
    void decorateState();
    void build2bodyMatrix();
};

#endif /* defined(__ED_solver__interaction__) */
