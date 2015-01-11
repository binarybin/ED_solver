//
//  main.cpp
//  ED_solver
//
//  Created by Bin Xu on 12/19/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include "measurement.h"
#include "output.h"

bool interest(const Pxy & p) // describe the interesting sectors
{
    if (p.px == 0 || p.py == 0) return true;
    else return false;
}

int main(int argc, const char * argv[])
{
    Hamiltonian ham;
    ham.getParameters(); // This will invoke the keyboard input
    HilbertSpace hilbert_space(ham);
    hilbert_space.buildOrbitalList();
    hilbert_space.buildReverseOrbitalMap();
    hilbert_space.buildPairMap();
    hilbert_space.buildStateMap();
    
    vector<Measurement> measure_list;
    for (auto k : hilbert_space.getStateMap())
    if (interest(k.first))
    {
        cout<<"Sector "<<k.first<<" Dimension: "<<k.second.size()<<endl;
        Interaction interaction(hilbert_space, k.first); //This is the interaction matrix
        interaction.decorateState();
        interaction.buildKineticMatrix();
        interaction.build2bodyMatrix();
        Solver diagonalize(interaction, k.first);
        diagonalize.diagonalize();
        Measurement measurement(diagonalize, hilbert_space, ham);
        measurement.measure();
        measure_list.push_back(measurement);
    }
    
    for (auto it : measure_list)
    {
        cout<<it.p<<endl;
        for (int i = 0; i < it.eigenvalues.size(); i++)
        {
            cout<<"Eigenvalue: "<<it.eigenvalues[i]<<endl;
            cout<<"Density: "<<endl;
            for(auto it3 : it.one_state_results[i].occupation_nbr) cout<<it3<<" ";
            cout<<endl;
            cout<<"Form"<<endl;
            for (int k = 0; k < it.one_state_results[i].form.data[0][0].size(); k++)
            {
                for (int l = 0; l < it.one_state_results[i].form.data[0][0][0].size(); l++)
                {
                    cout<<setw(10)<<fixed<<real(it.one_state_results[i].form.data[0][0][k][l]);
                }
                cout<<endl;
            }
            cout<<endl;
        }
        cout<<endl;
    }
    
    return 0;
}
