//
//  main.cpp
//  ED_solver
//
//  Created by Bin Xu on 12/19/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include "measurement.h"
#include "output.h"

bool interest(const Pxy & p) // describe the interesting sectors
{
    if (p.px == 0 &&  p.py == 8) return true;
    else return false;
}

int main(int argc, const char * argv[])
{
    Hamiltonian ham;
    ham.getParameters(); // This will invoke the keyboard input
    HilbertSpace hilbert_space(ham);
    hilbert_space.buildOrbitalList();
    if (ham.nele == 999) ham.nele = ham.norb/2;
    cout<<ham.norb<<" Orbitals with "<<ham.nele<<" Electrons"<<endl;
    hilbert_space.buildReverseOrbitalMap();
    hilbert_space.buildPairMap();
    hilbert_space.buildStateMap();
    
    vector<Measurement> measure_list;
    fstream fout(argv[1], fstream::out|fstream::app);
    for (auto p : hilbert_space.getStateMap())
//    if (interest(p.first))
    {
        cout<<"Sector "<<p.first<<" Dimension: "<<p.second.size()<<endl;
        pair<double, double> k = ham.computeK(p.first.px, p.first.py);
        cout<<"Momentum: ("<<k.first<<", "<<k.second<<")"<<endl;
        Interaction interaction(hilbert_space, p.first); //This is the interaction matrix
        interaction.decorateState();
        interaction.buildKineticMatrix();
        interaction.build2bodyMatrix();
        Solver diagonalize(interaction, p.first);
        diagonalize.diagonalize();
        Measurement measurement(diagonalize, hilbert_space, ham);
        measurement.measure();
//        measure_list.push_back(measurement);
        fout<<measurement.p<<endl;
        cout<<measurement.p<<endl;
        pair<double, double> kk = ham.computeK(measurement.p.px, measurement.p.py);
        fout<<"Momentum: ("<<kk.first<<", "<<kk.second<<")"<<endl;
        cout<<"Momentum: ("<<kk.first<<", "<<kk.second<<")"<<endl;
        for (int i = 0; i < measurement.eigenvalues.size(); i++)
        {
            fout<<"Eigenvalues: "<<measurement.eigenvalues[i]<<endl;    
            cout<<"Eigenvalues: "<<measurement.eigenvalues[i]<<endl;    
        }
        fout<<endl;
        cout<<endl;
    }
    
    fout.close();
    

    for (auto it : measure_list)
    {
        cout<<it.p<<endl;
        pair<double, double> k = ham.computeK(it.p.px, it.p.py);
        cout<<"Momentum: ("<<k.first<<", "<<k.second<<")"<<endl;
        for (int i = 0; i < it.eigenvalues.size(); i++)
        {
            cout<<"Eigenvalue: "<<it.eigenvalues[i]<<endl;
//            cout<<"Density: "<<endl;
//            for(auto it3 : it.one_state_results[i].occupation_nbr) cout<<it3<<" ";
//            cout<<endl;
        }
        cout<<endl;
    }

    
    return 0;
}
