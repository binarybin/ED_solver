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
    vector<vector<vector<vector<complex<double> > > > > data;
    FormFactor(int norb, int max_qx, int max_qy)
    {
        data.resize(2);
        for (auto & it : data)
        {
            it.resize(2);
            for (auto & it2 : it)
            {
                it2.resize(2*max_qx + 1);
                for (auto & it3 : it2)
                    it3.resize(2 * max_qy + 1);
            }
        }
    }
};

struct OneStateMeasurement
{
    OneStateMeasurement(int norb, int max_qx, int max_qy): form(norb, max_qx, max_qy){}
    int state_id;
    vector<double> occupation_nbr;
    vector<double> spin_nbr;
    FormFactor form;
};

class Measurement
{
    HilbertSpace &hilbert_space;
    Solver &solver;
    Hamiltonian &ham;
    const int max_qx, max_qy;
public:
    vector<OneStateMeasurement> one_state_results;
    PwaveSc pwave_results;
    vector<double> eigenvalues;
    Pxy p;
    

    Measurement(Solver &o_solver, HilbertSpace &o_hilbert_space, Hamiltonian &o_ham) : hilbert_space(o_hilbert_space), solver(o_solver), ham(o_ham), pwave_results((int)solver.interaction.state_list.size()), p(o_solver.p), eigenvalues(solver.eigenvalues), max_qx(4), max_qy(4)
    {
        one_state_results.resize(solver.Nev, OneStateMeasurement((int)ham.norb, max_qx, max_qy));
    }
    void measureDensity(int num);
    void measureSz(int num);
    void measureBCS(int num){};
    void measurePwave(){};
    void measureForm(int num);
    
    
    void measure();
    
private: //auxiliary functions
    complex<double> formOneStep(int num, int state1, int state2, int qx, int qy, int k1, int k2);
    int formEle(int num, int state1, int state2, int qx, int qy, int &i, int j, int k1, int k2);
    complex<double> formWrapper(int num, int state1, int state2, int qx, int qy);
};

#endif /* defined(__ED_solver__measurement__) */
