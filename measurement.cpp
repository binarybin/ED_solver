//
//  measurement.cpp
//  ED_solver
//
//  Created by Bin Xu on 12/21/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#include "measurement.h"

void Measurement::measureDensity(int num)
{
    one_state_results[num].occupation_nbr.resize(ham.norb);
    for(int i = 0; i < solver.ch_result[num].eigenvector.size(); i++) //state coefficient
    {
        for (int j = 0; j < ham.norb; j++) //check occupancy of orbital
        {
            if (solver.interaction.state_list[i].cstate[j])
                one_state_results[num].occupation_nbr[j] += norm(solver.ch_result[num].eigenvector[i]);
        }
    }
}

void Measurement::measureSz(int num)
{
    one_state_results[num].spin_nbr.resize(ham.norb/2);
    for (int i = 0; i < ham.norb; i+=2)
    {
        one_state_results[num].spin_nbr[i/2] += one_state_results[num].occupation_nbr[i]*0.5;
        one_state_results[num].spin_nbr[i/2] -= one_state_results[num].occupation_nbr[i+1]*0.5;
    }
}

void Measurement::measureFormFactor(int num)
{
    one_state_results[num].form_results;
    
}
