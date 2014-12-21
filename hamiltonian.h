//
//  hamiltonian.h
//  ED_solver
//
//  Created by Bin Xu on 12/19/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#ifndef __ED_solver__hamiltonian__
#define __ED_solver__hamiltonian__

#include <iostream>

#include "constants.h"
#include "output.h"

class HilbertSpace;
class Interaction;
class Solver;

class Hamiltonian
{
    friend class HilbertSpace;
    friend class Interaction;
    friend class Solver;
    double U, V1, vf, k0, delta; // interaction parameters
    size_t norb, nele, max_norb, pbc_type; // geometry parameters
    double epsilon_kx, epsilon_ky, aspect2; // geomatry parameters
    int lanczos_ne; // problem parameter
    
public:
    void getParameters(); // This method can be overloaded
    size_t getNorb() {return norb;}
};

#endif /* defined(__ED_solver__hamiltonian__) */
