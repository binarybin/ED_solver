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
    vector<ch_eigenvector> ch_result;
    Interaction &interaction;
    
public:
    Solver(Interaction &o_interaction): interaction(o_interaction){}
    void diagonalize()
    {
        if(interaction.hilbert_space.ham.vf == 0)
            this->rsDiagonalize();
        else
            this->chDiagonalize();
    }
    void rsDiagonalize();
    void chDiagonalize();
};

#endif /* defined(__ED_solver__solver__) */
