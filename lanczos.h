//-------------------------------------------------------------------------//
//              C++ header file for the Lanczos library                    //
//     Written by Bin Xu                                                   //
//     Version 0.2                                                         //
//     Last modified: 6 Aug 2014                                           //
//-------------------------------------------------------------------------//


#ifndef LANCZOSDIAG_H
#define LANCZOSDIAG_H

#include <complex.h>
#include <vector>
using namespace std;
// Structure to store eigenvectors
struct rs_eigenvector
{
	vector<double> eigenvector;
	double eigenvalue;
};

struct ch_eigenvector
{
	vector<complex<double> > eigenvector;
	double eigenvalue;
};

// Lanczos methods with a function pointer input as matrix vector multiplication
void lanczos_diag(int dim, int nevec, void (*matvec) (int*, double*, double*, bool*), vector<double>& eigenvalues, vector<double>& variance, vector<rs_eigenvector>& eigenvectors, int maxstep=1000, int report=1000, int seed=123456);
void lanczos_diag(int dim, int nevec, void (*matvec) (int*, double*, double*, bool*), vector<double>& eigenvalues, vector<double>& variance, int maxstep=1000, int report=1000, int seed=123456);
void lanczos_diag(int dim, int nevec, void (*matvec) (int*, complex<double> *, complex<double> *, bool*), vector<double>& eigenvalues, vector<double>& variance, vector<ch_eigenvector>& eigenvectors, int maxstep=1000, int report=1000, int seed=123456);
void lanczos_diag(int dim, int nevec, void (*matvec) (int*, complex<double> *, complex<double> *, bool*), vector<double>& eigenvalues, vector<double>& variance, int maxstep=1000, int report=1000, int seed=123456);
// Lanczos methods for sparse matrices stored in a sparse matrix structure 
struct rs_sparse_matrix
{
	int bra, ket;
	double element;
};

struct ch_sparse_matrix
{
	int bra, ket;
	complex<double>  element;
};

void lanczos_diag(int dim, int nevec, vector<rs_sparse_matrix>& matrix, vector<double>& eigenvalues, vector<double>& variance, vector<rs_eigenvector>& eigenvectors, int maxstep=1000, int report=1000, int seed=123456);
void lanczos_diag(int dim, int nevec, vector<rs_sparse_matrix>& matrix, vector<double>& eigenvalues, vector<double>& variance, int maxstep=1000, int report=1000, int seed=123456);
void lanczos_diag(int dim, int nevec, vector<ch_sparse_matrix>& matrix, vector<double>& eigenvalues, vector<double>& variance, vector<ch_eigenvector>& eigenvectors, int maxstep=1000, int report=1000, int seed=123456);
void lanczos_diag(int dim, int nevec, vector<ch_sparse_matrix>& matrix, vector<double>& eigenvalues, vector<double>& variance, int maxstep=1000, int report=1000, int seed=123456);

// Declaration of the interface to the corresponding Fortran routine
extern "C" void lanczos_diag_rs_(int *n, int *nevec, void (*matvec) (int*, double*, double*, bool*), double *eval, double *evec, double *variance, int *number, double *resolution, int *maxstep, int *report, int *seed);

extern "C" void lanczos_diag_ch_(int *n, int *nevec, void (*matvec) (int*, complex<double> *, complex<double> *, bool*), double *eval, complex<double>  *evec, double *variance, int *number, double *resolution, int *maxstep, int *report, int *seed);



#endif /* LANCZOSDIAG_H */
