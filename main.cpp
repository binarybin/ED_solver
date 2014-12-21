//
//  main.cpp
//  ED_solver
//
//  Created by Bin Xu on 12/19/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#include <iostream>
#include "measurement.h"

int main(int argc, const char * argv[])
{
    Hamiltonian ham;
    ham.getParameters(); // This will invoke the keyboard input
    HilbertSpace hilbert_space(ham);
    hilbert_space.buildOrbitalList();
    hilbert_space.buildReverseOrbitalMap();
    hilbert_space.buildPairMap();
    hilbert_space.buildStateMap();
    
    for (auto k : hilbert_space.getStateMap())
    {
        Interaction interaction(hilbert_space, k.first); //This is the interaction matrix
        interaction.decorateState();
        interaction.buildKineticMatrix();
        interaction.build2bodyMatrix();
        Solver diagonalize(interaction);
        diagonalize.diagonalize();
        Measurement measurement(diagonalize, hilbert_space, ham);
        measurement.measure();
    }
    
    return 0;
}
