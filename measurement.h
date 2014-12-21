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
    vector<rs_eigenvector> &rs_result;
    vector<ch_eigenvector> &ch_result;
    HilbertSpace &hilbert_space;
    
public:
    Measurement(Solver &solver, HilbertSpace &hilbert_space, Hamiltonian &ham);
    void measure();
};

#endif /* defined(__ED_solver__measurement__) */
