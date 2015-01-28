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
    std::cout<<"Please input parameters:\n U, V1, Vf, delta, Max Norb, PBC, Ne (999 for half-filling), Nev, I_x, I_y, J_x, J_y, C_x, C_y, MaxLapackSize" <<std::endl;

    std::cin>>U>>V1>>vf>>delta>>max_norb>>pbc_type>>nele>>lanczos_ne>>epsilon_Ix>>epsilon_Iy>>epsilon_Jx>>epsilon_Jy>>center_x>>center_y>>MaxLapackSize;
}

pair<double, double> Hamiltonian::computeK(int i, int j)
{
    return pair<double, double>((i-center_x)*epsilon_Ix + (j-center_x)*epsilon_Jx, (i-center_y)*epsilon_Iy + (j-center_y)*epsilon_Jy);
}
