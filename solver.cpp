//
//  solver.cpp
//  ED_solver
//
//  Created by Bin Xu on 12/19/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#include "solver.h"

int *fast_bra_list;
int *fast_ket_list;
complex<double> *fast_amp_list;
int fast_size;

void matvec(int *size, complex<double> *vec_in, complex<double> *vec_out, bool *add);

void matvec(int *size, complex<double> *vec_in, complex<double> *vec_out, bool *add)
{
    if (!(*add))
        for (int i = 0; i < (*size); i++)
            vec_out[i] = 0.0;
    
    for (int i =0; i < fast_size; i++)
        vec_out[fast_bra_list[i]] += fast_amp_list[i] * vec_in[fast_ket_list[i]];
}

void Solver::chDiagonalize()
{
    fast_size = interaction.matrix.size();
    cout<<"fast size: "<<fast_size<<endl;
    fast_bra_list = new int [fast_size];
    fast_ket_list = new int [fast_size];
    fast_amp_list = new complex<double> [fast_size];
    
    int i = 0;
    for (auto it = interaction.matrix.begin(); it != interaction.matrix.end(); it++)
    {
        fast_bra_list[i] = it->bra;
        fast_ket_list[i] = it->ket;
        fast_amp_list[i] = it->amplitude;
        i++;
    }
    
    vector<double> variance;
    
    lanczos_diag((int)interaction.state_list.size(), interaction.hilbert_space.ham.lanczos_ne, matvec, eigenvalues, variance, ch_result);

    delete [] fast_bra_list;
    delete [] fast_ket_list;
    delete [] fast_amp_list;
}
