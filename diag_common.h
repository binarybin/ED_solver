//
//  diag_common.h
//  ED_solver
//
//  Created by Bin Xu on 12/22/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#ifndef ED_solver_diag_common_h
#define ED_solver_diag_common_h

#include <vector>
#include <complex>
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

#endif
