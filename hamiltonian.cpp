//
//  hamiltonian.cpp
//  ED_solver
//
//  Created by Bin Xu on 12/19/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#include "hamiltonian.h"

void Hamiltonian::getParameters()
{
    std::cout<<"Please input parameters:\n U, V1, Vf, delta, Max Norb, PBC, aspect^2, Ne, Nev" <<std::endl;

    std::cin>>U>>V1>>vf>>delta>>max_norb>>pbc_type>>aspect2>>nele>>lanczos_ne;
}