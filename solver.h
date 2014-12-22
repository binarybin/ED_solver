//
//  solver.h
//  ED_solver
//
//  Created by Bin Xu on 12/19/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#ifndef __ED_solver__solver__
#define __ED_solver__solver__

#include "interaction.h"
#include "lanczos.h"

class Solver
{
    vector<rs_eigenvector> rs_result;
    Interaction &interaction;
    
    
public:
    Solver(Interaction &o_interaction, Pxy o_p): interaction(o_interaction), p(o_p){}
    Pxy p;
    vector<ch_eigenvector> ch_result;
    vector<double> eigenvalues;
    void diagonalize()
    {
//        if(interaction.hilbert_space.ham.vf == 0)
//            this->rsDiagonalize();
//        else
//            this->chDiagonalize();
        chDiagonalize();
        cout<<"Result: "<<endl;
        for (auto it : ch_result)
        {
            cout<<it.eigenvalue<<" ";
        }
        cout<<endl;
    }
    void rsDiagonalize();
    void chDiagonalize();
};

void matvec(int *size, complex<double> *vec_in, complex<double> *vec_out, bool *add);

#endif /* defined(__ED_solver__solver__) */
