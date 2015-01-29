//
//  solver.h
//  ED_solver
//
//  Created by Bin Xu on 12/19/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#ifndef __ED_solver__solver__
#define __ED_solver__solver__

#include <algorithm>
#include "interaction.h"
#include "lanczos.h"
#include "diag_wrapper.h"

class Measurement;

class Solver
{
    friend class Measurement;
    Interaction &interaction;
    bool ch, lapack;
    size_t Nev;
    
public:
    Solver(Interaction &o_interaction, Orbital o_orb): interaction(o_interaction), orb_sector(o_orb)
    {
        Nev = std::min((size_t)interaction.hilbert_space.ham.lanczos_ne, interaction.state_list.size());
    }
    Orbital orb_sector;
    vector<ch_eigenvector> ch_result;
    vector<rs_eigenvector> rs_result;
    vector<double> eigenvalues;

    void diagonalize()
    {
        if (interaction.hilbert_space.ham.vf == -100)
        {
            ch = false;
            if (interaction.state_list.size() < interaction.hilbert_space.ham.MaxLapackSize)
            {
                rsLapackDiagonalize();
                lapack = true;
            }
            else
            {
                rsLanczosDiagonalize();
                lapack = false;
            }
        }
        else
        {
            ch = true;
            if (interaction.state_list.size() < interaction.hilbert_space.ham.MaxLapackSize)
            {
                chLapackDiagonalize();
                lapack = true;
            }
            else
            {
                chLanczosDiagonalize();
                lapack = false;
            }
        }
    }
    void rsLanczosDiagonalize(){};
    void chLanczosDiagonalize();
    void rsLapackDiagonalize(){};
    void chLapackDiagonalize();
    
    void testEvec(int num);
    void testNormalization(int num);
    
};

ostream& operator<<(ostream &os, ch_eigenvector &ch_result);

void matvec(int *size, complex<double> *vec_in, complex<double> *vec_out, bool *add);

#endif /* defined(__ED_solver__solver__) */
