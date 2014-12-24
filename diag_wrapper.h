//
//  diag_wrapper.h
//  SPT
//
//  Created by Bin Xu on 9/5/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#ifndef __SPT__diag_wrapper__
#define __SPT__diag_wrapper__

#include <iostream>
#include <vector>
#include <complex>
#include "diag_common.h"
using namespace std;

bool zheevpp(vector<vector<complex<double> > >& matrix, vector<double>& evals, vector<ch_eigenvector>& evecs, int nevecs);
#endif /* defined(__SPT__diag_wrapper__) */
