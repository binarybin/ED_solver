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

struct BcsSc
{
    vector<vector<complex<double> > > corr_mat;
    BcsSc(int norb)
    {
        int nk = norb/2;
        corr_mat.resize(nk);
        for (auto & it : corr_mat) it.resize(nk);
    }
};

struct PwaveSc
{
    vector<vector<vector<vector<complex<double> > > > > corr_mat;
    PwaveSc(int norb)
    {
        int nk = norb/2;
        corr_mat.resize(2);
        for (auto & it : corr_mat)
        {
            it.resize(2);
            for (auto & it2 : corr_mat)
            {
                it2.resize(nk);
                for (auto & it3 : corr_mat)
                    it3.resize(nk);
            }
        }
    }
};

struct FormFactor
{
    vector<vector<double> > form_mat;
    FormFactor(int norb)
    {
        int nk = norb/2;
        form_mat.resize(nk);
        for (auto & it : form_mat)
            it.resize(nk);
    }
};

struct OneStateMeasurement
{
    int state_id;
    vector<double> eigenvalues;
    vector<double> occupation_nbr;
    vector<double> spin_nbr;
    BcsSc bcs_results;
    FormFactor form_results;
};

class Measurement
{
    HilbertSpace &hilbert_space;
    Solver &solver;
    Hamiltonian &ham;
    
    vector<OneStateMeasurement> one_state_results;
    PwaveSc pwave_results;
    
    
public:
    Measurement(Solver &o_solver, HilbertSpace &o_hilbert_space, Hamiltonian &o_ham) : hilbert_space(o_hilbert_space), solver(o_solver), ham(o_ham), pwave_results((int)solver.interaction.state_list.size())
    {
        one_state_results.resize(ham.lanczos_ne);
    }
    void measureDensity(int num);
    void measureSz(int num);
    void measureBCS(int num);
    void measurePwave();
    void measureFormFactor(int num);
    
    void measure()//This function is what we modify the most
    {}
};

#endif /* defined(__ED_solver__measurement__) */
