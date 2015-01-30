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

bool interest(const Orbital & orb) // describe the interesting sectors
{
    if (orb.p.px == 0 &&  orb.p.py == 8) return true;
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
    for (auto orb : hilbert_space.getStateMap())
//    if (interest(orb.first))
    {
        cout<<"Sector "<<orb.first<<" Dimension: "<<orb.second.size()<<endl;
        pair<double, double> k = ham.computeK(orb.first.p.px, orb.first.p.py);
        cout<<"Momentum: ("<<k.first<<", "<<k.second<<")"<<endl;
        Interaction interaction(hilbert_space, orb.first); //This is the interaction matrix
        interaction.decorateState();
        interaction.build2bodyMatrix();
        Solver diagonalize(interaction, orb.first);
        diagonalize.diagonalize();
        Measurement measurement(diagonalize, hilbert_space, ham);
        measurement.measure();
//        measure_list.push_back(measurement);
        fout<<measurement.orb_sector<<endl;
        cout<<measurement.orb_sector<<endl;
        pair<double, double> kk = ham.computeK(measurement.orb_sector.p.px, measurement.orb_sector.p.py);
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
        cout<<it.orb_sector<<endl;
        pair<double, double> k = ham.computeK(it.orb_sector.p.px, it.orb_sector.p.py);
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
