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

void Solver::chLanczosDiagonalize()
{
    fast_size = (int)interaction.matrix.size();
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
    
    lanczos_diag((int)interaction.state_list.size(), (int)Nev, matvec, eigenvalues, variance, ch_result);

    delete [] fast_bra_list;
    delete [] fast_ket_list;
    delete [] fast_amp_list;
}

//bool zheevpp(vector<vector<complex<double> > >& matrix, vector<double>& evals, vector<ch_eigenvector>& evecs, int nevecs);
void Solver::chLapackDiagonalize()
{
    int size = (int) interaction.state_list.size();
    vector<vector<complex<double> > > ham_mat(size);
    for (auto & it : ham_mat)
        it.resize(size);
    
    for (auto it : interaction.matrix)
    {
        ham_mat[it.bra][it.ket] += it.amplitude;
    }
    
    
    if (!zheevpp(ham_mat, eigenvalues, ch_result, (int)Nev))
    {
        cout<<"Error while diagonalizing the Hamiltonian using lapack."<<endl;
    }

}

void Solver::testEvec(int num) // This function is not supposed to be used...
{
    complex<double> *ev = new complex<double>[ch_result[num].eigenvector.size()];
    complex<double> *evout = new complex<double>[ch_result[num].eigenvector.size()];
    
    for (int i = 0; i < ch_result[num].eigenvector.size(); i++)
    {
        ev[i] = ch_result[num].eigenvector[i];
    }
    
    void matvec(int *size, complex<double> *vec_in, complex<double> *vec_out, bool *add);
    
    int size = interaction.state_list.size();
    bool FF = false;
    
    matvec(&size, ev, evout, &FF);
    
    cout<<"testing evec"<<endl;
    
    for(int i = 0; i < 100; i++) //row
    {
        cout<<"before:\t"<<ev[i];
        cout<<"\tafter:\t"<<evout[i];
        cout<<"\tratio:\t"<<evout[i]/ev[i];
        cout<<"\tEigenvalue:\t"<<eigenvalues[num]<<endl;
    }
    cout<<endl;

}

void Solver::testNormalization(int num)
{
    double Norm = 0;
    for(auto it : ch_result[num].eigenvector)
    {
        Norm += norm(it);
    }
    cout<<"Norm = "<<Norm<<endl;
}

ostream& operator<<(ostream &os, ch_eigenvector &ch_result)
{
    os<<"Eigenvalue: "<<ch_result.eigenvalue<<endl;
    os<<"Eigenvector: "<<endl;
    for (int i = 0; i < ch_result.eigenvector.size(); ++i)
    {
        os<<fixed<<ch_result.eigenvector[i]<<"\t";
    }
    os<<endl;
    return os;
}
