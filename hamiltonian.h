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
class Measurement;


class Hamiltonian
{
    friend class HilbertSpace;
    friend class Interaction;
    friend class Solver;
    friend class Measurement;
    double U, V1, vf, k0, delta; // interaction parameters
    size_t norb, nele, max_norb, pbc_type; // geometry parameters
    double epsilon_Ix, epsilon_Iy, epsilon_Jx, epsilon_Jy, center_x, center_y; // geomatry parameters
    int lanczos_ne; // problem parameter
    int MaxLapackSize;
    
public:
    void getParameters(); // This method can be overloaded
    size_t getNorb() {return norb;}
    pair<double, double> computeK(int i, int j);
};

#endif /* defined(__ED_solver__hamiltonian__) */
