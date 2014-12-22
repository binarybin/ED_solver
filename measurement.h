//
//  measurement.h
//  ED_solver
//
//  Created by Bin Xu on 12/21/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#ifndef __ED_solver__measurement__
#define __ED_solver__measurement__

#include "solver.h"

class Measurement
{
    HilbertSpace &hilbert_space;
    Solver &solver;
    Hamiltonian &ham;
    
public:
    Measurement(Solver &o_solver, HilbertSpace &o_hilbert_space, Hamiltonian &o_ham) : hilbert_space(o_hilbert_space), solver(o_solver), ham(o_ham) {}
    void measure(){}
};

#endif /* defined(__ED_solver__measurement__) */
