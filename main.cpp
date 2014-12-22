//
//  main.cpp
//  ED_solver
//
//  Created by Bin Xu on 12/19/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#include <iostream>
#include "measurement.h"
#include "output.h"

int main(int argc, const char * argv[])
{
    Hamiltonian ham;
    ham.getParameters(); // This will invoke the keyboard input
    HilbertSpace hilbert_space(ham);
    hilbert_space.buildOrbitalList();
    hilbert_space.buildReverseOrbitalMap();
    hilbert_space.buildPairMap();
    hilbert_space.buildStateMap();
    
    vector<Solver> solver_list;
    for (auto k : hilbert_space.getStateMap())
    {
        cout<<"Sector "<<k.first<<" Dimension: "<<k.second.size()<<endl;
        Interaction interaction(hilbert_space, k.first); //This is the interaction matrix
        interaction.decorateState();
        interaction.buildKineticMatrix();
        interaction.build2bodyMatrix();
        if (k.second.size() > 100) 
        {
            Solver diagonalize(interaction, k.first);
            diagonalize.diagonalize();
            Measurement measurement(diagonalize, hilbert_space, ham);
            measurement.measure();
            solver_list.push_back(diagonalize);
        }
        
        
    }
    
    for (auto it : solver_list)
    {
        cout<<it.p<<endl;
        for (auto it2 : it.eigenvalues)
        {
            cout<<it2<<" ";
        }
        cout<<endl;
    }
    
    /*   Pxy p(0, 0);
     Interaction interaction(hilbert_space, p);
     interaction.decorateState();
     interaction.buildKineticMatrix();
     interaction.build2bodyMatrix();
     Solver diagonalize(interaction);
     diagonalize.diagonalize();
     Measurement measurement(diagonalize, hilbert_space, ham);
     measurement.measure();*/
    
    
    return 0;
}
