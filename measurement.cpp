//
//  measurement.cpp
//  ED_solver
//
//  Created by Bin Xu on 12/21/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//
#include <stdexcept>

#include "measurement.h"


void Measurement::measure()//This function is what we modify the most
{
    for (int i = 0; i < solver.Nev; i++)
    {
        measureDensity(i);
        measureSz(i);
//        measureForm(i);
    }
}

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

int Measurement::formEle(int num, int state1, int state2, int qx, int qy, int &i, int j, int k1, int k2)
{
    //<i|c^(k+q) c(k) c^(k'-q) c(k')|j>
    
    int sign_counter = 0;
    int kill2 = 2 * k2 + state2; // The first killed orbital
    if (!solver.interaction.state_list[j].cstate[kill2]) return 0;
    
    auto temp_state = solver.interaction.state_list[j].cstate;
    for (int i = 0; i < kill2; i++)
        if (temp_state[i])
            sign_counter++;
    
    temp_state[kill2] = false;
    
    Orbital cre_orb_2(hilbert_space.orbital_list[kill2].p.px - qx, hilbert_space.orbital_list[kill2].p.py - qy, hilbert_space.orbital_list[kill2].spin);
    size_t create_2 = 10000;
    
    try
    {
        create_2 = hilbert_space.reverse_orbital_map.at(cre_orb_2);
    }
    catch (const std::out_of_range& oor) // This orbital does not exist
    {
        return 0;
    }
    
    if (temp_state[create_2]) return 0;

    for (int i = 0; i < create_2; i++)
        if (temp_state[i])
            sign_counter++;
    
    temp_state[create_2] = true;
    int kill_1 = 2 * k1 + state1;
    if (!temp_state[kill_1]) return 0;
    
    for (int i = 0; i < kill_1; i++)
        if (temp_state[i])
            sign_counter++;
    
    temp_state[kill_1] = false;
    size_t create_1 = 10000;
    Orbital cre_orb_1(hilbert_space.orbital_list[kill_1].p.px + qx, hilbert_space.orbital_list[kill_1].p.py + qy, hilbert_space.orbital_list[kill_1].spin);
    try
    {
        create_1 = hilbert_space.reverse_orbital_map.at(cre_orb_1);
    }
    catch (const std::out_of_range& oor)
    {
        return 0;
    }
    
    if (temp_state[create_1]) return 0;
    
    for (int i = 0; i < create_1; i++)
        if (temp_state[i])
            sign_counter++;
    
    temp_state[create_1] = true;
    
    i = solver.interaction.reverse_state_map.at(temp_state.to_ulong());
    
    if (sign_counter % 2 == 0) {
        return 1;
    }
    else
        return -1;

}

complex<double> Measurement::formOneStep(int num, int state1, int state2, int qx, int qy, int k1, int k2)
{
    // \sum_{i, j}<i|c^(k+q) c(k) c^(k'-q) c(k')|j> A*(i) A(j)
    complex<double> result = 0;
    auto state_vec = solver.ch_result[num].eigenvector;
    for (int j = 0; j < state_vec.size(); j++)
    {
        int i = 0;
        int sign_counter = formEle(num, state1, state2, qx, qy, i, j, k1, k2);
        result += conj(state_vec[i]) * state_vec[j] * (complex<double>)sign_counter;
    }
    return result;
}

complex<double> Measurement::formWrapper(int num, int state1, int state2, int qx, int qy)
{
    // sum_{k1, k2}\sum_{i, j}<i|c^(k+q) c(k) c^(k'-q) c(k')|j> A*(i) A(j)
    complex<double> result = 0;
    for (int k1 = 0; k1 < hilbert_space.orbital_list.size(); k1++)
        for (int k2 = 0; k2 < hilbert_space.orbital_list.size(); k2++)
            result += formOneStep(num, state1, state2, qx, qy, k1, k2);
    
    return result;
}

void Measurement::measureForm(int num)
{
    one_state_results[num].form.data[0][0][0][0] = 0;
    for (int state1 = 0; state1 < 2; state1++) // spin index, out
    for (int state2 = 0; state2 < 2; state2++) // spin index, in
    for (int qx = -max_qx; qx < max_qx; qx++) // momentum index
    for (int qy = -max_qy; qy < max_qy; qy++) // momentum index
        one_state_results[num].form.data[state1][state2][qx+max_qx][qy+max_qy] = formWrapper(num, state1, state2, qx, qy);
}