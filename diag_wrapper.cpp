//
//  diag_wrapper.cpp
//  SPT
//
//  Created by Bin Xu on 9/5/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#include "diag_wrapper.h"

#include <complex.h>
#include <Accelerate/Accelerate.h>

bool zheevpp(vector<vector<complex<double> > >& matrix, vector<double>& evals, vector<ch_eigenvector>& evecs, int nevecs)
// The interface should be fixed for different implementations.
// nevecs = 0 gives no eigenvectors
// evec[n] gives the nth eigenvector
// This version is the gnu lapacke implementation
{
    __CLPK_integer size = matrix.size();
    __CLPK_doublecomplex *mat = new __CLPK_doublecomplex[size*size];
    
    for (int row = 0; row < size; row++)
    for (int col = 0; col < size; col++)
    {
        mat[row * size + col].r = real(matrix[row][col]);
        mat[row * size + col].i = imag(matrix[row][col]);
    }
    
    __CLPK_integer lwork = 2 * size;
    double *Evals = new double[lwork];
    __CLPK_doublecomplex *work = new  __CLPK_doublecomplex[2*size];
    double *rwork = new double[3 * size];
    int info;
    char tempV = 'V';
    char tempU = 'U';
    
    zheev_(&tempV, &tempU, &size, mat, &size, Evals, work, &lwork, rwork, &info);
    
    if (info > 0)
    {
        cout << "Error in Lapack while computing eigenvalues." <<endl;
        return false;
    }
    evecs.resize(nevecs);
    evals.resize(nevecs);
    for (auto& it : evecs) it.eigenvector.resize(size);
    
    for(int i = 0; i < nevecs; i++)
    {
        for (int j = 0; j < size; j++)
        {
            evecs[i].eigenvector[j] = mat[i*size + j].r + complex<double>(0, 1) * mat[j*size + i].i ; //I may have mixed up i and j
        }
        evals[i] = Evals[i];
    }
    
    for (int i = 0; i < nevecs; i++)
    {
    
        double Norm = 0;
        for(auto it : evecs[i].eigenvector)
        {
            Norm += norm(it);
        }
        
        Norm = sqrt(Norm);
        
        for (int j = 0; j < size; j++)
        {
            evecs[i].eigenvector[j] /= Norm;
        }
    }
    
    delete [] Evals;
    delete [] mat;

    return true;
}