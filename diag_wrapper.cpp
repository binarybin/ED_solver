//
//  diag_wrapper.cpp
//  SPT
//
//  Created by Bin Xu on 9/5/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#include "diag_wrapper.h"
#include "mkl_lapacke.h"
#include <complex.h>

bool zheevpp(vector<vector<complex<double> > >& matrix, vector<double>& evals, vector<ch_eigenvector>& evecs, int nevecs)
// The interface should be fixed for different implementations.
// nevecs = 0 gives no eigenvectors
// evec[n] gives the nth eigenvector
// This version is the Intel Lapacke implementation
// Tested to be working on della.princeton.edu
{
    size_t size = matrix.size();
    MKL_Complex16 *mat = new MKL_Complex16[size*size];
    
    for (int row = 0; row < size; row++)
    for (int col = 0; col < size; col++)
    {
        mat[col * size + row].real = real(matrix[row][col]);
        mat[col * size + row].imag = imag(matrix[row][col]);
    }
    
    double *Evals = new double[size];
    int info;
    
    info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'L', size, mat, size, Evals);
    
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
            evecs[i].eigenvector[j] = mat[j*size + i].real + complex<double>(0, 1) * mat[i*size + j].imag ; //I may have mixed up i and j
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
