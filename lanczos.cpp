//
//  lanczos.cpp
//  SPT
//
//  Created by Bin Xu on 9/6/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#include "lanczos.h"
#include <iostream>
using namespace std;

void lanczos_diag(int dim, int nevec, void (*matvec) (int*, double*, double*, bool*), vector<double>& eigenvalues, vector<double>& variance, vector<rs_eigenvector>& eigenvectors, int maxstep, int report, int seed)
{
	double *eval = new double[nevec];
	double *vari = new double[nevec];
	int *number = nullptr;
	double *reso = new double;
	double *evec = new double[nevec * dim];
	
	lanczos_diag_rs_(&dim, &nevec, matvec, eval, evec, vari, number, reso, &maxstep, &report, &seed);
	
	eigenvalues.resize(nevec);
	eigenvectors.resize(nevec);
	variance.resize(nevec);
	
	int k = 0;
	
	for (int i = 0; i < nevec; i++)
	{
		eigenvalues[i] = eval[i];
		variance[i] = vari[i];
		eigenvectors[i].eigenvector.resize(dim);
		for(int j = 0; j < dim; j++)
		{
			eigenvectors[i].eigenvector[j] = evec[k];
			k++;
		}
		eigenvectors[i].eigenvalue = eval[i];
	}
	
	delete [] eval;
	delete [] vari;
	delete  reso;
	delete [] evec;
    
}

void lanczos_diag(int dim, int nevec, void (*matvec) (int*, double*, double*, bool*), vector<double>& eigenvalues, vector<double>& variance, int maxstep, int report, int seed)
{
    double *eval = new double[nevec];
    double *vari = new double[nevec];
    int *number = nullptr;
    double *reso = new double;
    double *evec = new double[nevec * dim];
    lanczos_diag_rs_(&dim, &nevec, matvec, eval, evec, vari, number, reso, &maxstep, &report, &seed);
    
    eigenvalues.resize(nevec);
    variance.resize(nevec);
    
    for (int i = 0; i < nevec; i++)
    {
        eigenvalues[i] = eval[i];
        variance[i] = vari[i];
    }
    
    delete [] eval;
    delete [] vari;
    delete  reso;
    delete [] evec;
}

void lanczos_diag(int dim, int nevec, void (*matvec) (int*, complex<double>  *, complex<double>  *, bool*), vector<double>& eigenvalues, vector<double>& variance, vector<ch_eigenvector>& eigenvectors, int maxstep, int report, int seed)
{
    double *eval = new double[nevec];
    double *vari = new double[nevec];
    int *number = nullptr;
    double *reso = new double;
    complex<double>  *evec = new complex<double> [nevec * dim];
    lanczos_diag_ch_(&dim, &nevec, matvec, eval, evec, vari, number, reso, &maxstep, &report, &seed);
    cout<<nevec<<endl;
    eigenvalues.resize(nevec);
    eigenvectors.resize(nevec);
    variance.resize(nevec);
    
    int k = 0;
    
    for (int i = 0; i < nevec; i++)
    {
        eigenvalues[i] = eval[i];
        variance[i] = vari[i];
        eigenvectors[i].eigenvector.resize(dim);
        for(int j = 0; j < dim; j++)
        {
            eigenvectors[i].eigenvector[j] = evec[k];
            k++;
        }
        eigenvectors[i].eigenvalue = eval[i];
    }
    
    delete [] eval;
    delete [] vari;
    delete  reso;
    delete [] evec;
    
}

void lanczos_diag(int dim, int nevec, void (*matvec) (int*, complex<double>  *, complex<double>  *, bool*), vector<double>& eigenvalues, vector<double>& variance, int maxstep, int report, int seed)
{
    double *eval = new double[nevec];
    double *vari = new double[nevec];
    int *number = nullptr;
    double *reso = new double;
    
    complex<double>  *evec = new complex<double> [nevec * dim];
    
    lanczos_diag_ch_(&dim, &nevec, matvec, eval, evec, vari, number, reso, &maxstep, &report, &seed);
    
    eigenvalues.resize(nevec);
    variance.resize(nevec);
    
    for (int i = 0; i < nevec; i++)
    {
        eigenvalues[i] = eval[i];
        variance[i] = vari[i];
    }
    
    delete [] eval;
    delete [] vari;
    delete  reso;
    delete [] evec;
    
}
