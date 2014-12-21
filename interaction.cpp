//
//  interaction.cpp
//  ED_solver
//
//  Created by Bin Xu on 12/19/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#include "interaction.h"

void Interaction::decorateState()
{
    auto cstates = hilbert_space.state_map.at(p);
    for(int i = 0; i < cstates.size(); i++)
    {
        State temp(cstates[i]);
        temp.state_id = i;
        reverse_state_map[cstates[i].to_ulong()] = i;
        state_list.push_back(temp);
    }
}

void Interaction::buildKineticMatrix()
{
    MatEle mat_ele;
    Hamiltonian& ham = hilbert_space.ham;
    for (auto it : state_list)
    {
        for (int i = 0; i < ham.norb; i += 2)
        {
            if (it.cstate[i] != it.cstate[i + 1])
            {
                CompactState temp_cstate;
                temp_cstate = it.cstate;
                if (temp_cstate[i])
                    temp_cstate[i] = 0, temp_cstate[i + 1] = 1;
                else
                    temp_cstate[i] = 1, temp_cstate[i + 1] = 0;
                mat_ele.bra = it.state_id;
                mat_ele.ket = reverse_state_map.at(temp_cstate.to_ulong());

                // Attention: 2 states that we operated on are neighbours with 0
                // and 1 filling, which means the 2 signs in fermionic second
                // quantized form are the same and they multiply to 1
                Pxy tempp = hilbert_space.orbital_list[i].p;
                if (it.cstate[i])
                    mat_ele.amplitude = complex<double>(ham.vf*ham.epsilon_kx*tempp.px,-ham.vf*ham.epsilon_ky*tempp.py)
                    + complex<double>(ham.delta * ham.epsilon_ky*tempp.py,  ham.delta * ham.epsilon_kx * tempp.px);
                else
                    mat_ele.amplitude = complex<double>(ham.vf*ham.epsilon_kx*tempp.px, ham.vf*ham.epsilon_ky*tempp.py)
                    + complex<double>(ham.delta * ham.epsilon_ky*tempp.py, -ham.delta * ham.epsilon_kx * tempp.px);
                if(abs(real(mat_ele.amplitude))>SmallDouble || abs(imag(mat_ele.amplitude))>SmallDouble)
                    matrix.push_back(mat_ele);
            }
        }
    }

}

void Interaction::build2bodyMatrix()
{
    Hamiltonian &ham = hilbert_space.ham;
    auto &orblist = hilbert_space.orbital_list;
    int count = 0;
    for (auto it : state_list)//go through all states in the list of "states"
    {
        //for a specific state, take a look at the first element of a pair
        for (int pos1 = 0; pos1 < ham.norb; pos1++) if (it.cstate[pos1])
            for (int pos2 = 0; pos2 < ham.norb; pos2++)if (it.cstate[pos2])
                if (pos2 != pos1)
                {
                    Orbital totalprop = hilbert_space.orbital_list[pos1] + hilbert_space.orbital_list[pos2];
                    //get pairs with same spin and momentum quantum numbers
                    vector<OrbPair> possible_pairs = hilbert_space.pair_map.at(totalprop);
                    for (auto it2 : possible_pairs)
                    {
                        int new1 = hilbert_space.reverse_orbital_map.at(it2.orb1);
                        int new2 = hilbert_space.reverse_orbital_map.at(it2.orb2);
                        CompactState tempcstate = it.cstate;
                        tempcstate[pos1] = 0;
                        tempcstate[pos2] = 0;
                        
                        if ((!tempcstate[new1]) && (!tempcstate[new2]))
                        {
                            if (it2.orb1.spin == it2.orb2.spin)
                            {
                                //dirty block
                                CompactState tempcstate(it.cstate);
                                
                                int sign_counter = 0;
                                
                                //applying the operators
                                for (int i = 0; i < pos1; i++) if (tempcstate[i]) sign_counter++;
                                tempcstate[pos1] = 0;
                                for (int i = 0; i < pos2; i++) if (tempcstate[i]) sign_counter++;
                                tempcstate[pos2] = 0;
                                for (int i = 0; i < new2; i++) if (tempcstate[i]) sign_counter++;
                                tempcstate[new2] = 1;
                                for (int i = 0; i < new1; i++) if (tempcstate[i]) sign_counter++;
                                tempcstate[new1] = 1;
                                
                                int newid = reverse_state_map.at(tempcstate.to_ulong());
                                
                                MatEle mat_ele;
                                mat_ele.bra = it.state_id;
                                mat_ele.ket = newid;
                                double amplitude=ham.norb/8.0*ham.V1*((orblist[new1].p.px*ham.epsilon_kx-orblist[pos1].p.px*ham.epsilon_kx)*\
                                                                      (orblist[new1].p.px*ham.epsilon_kx-orblist[pos1].p.px*ham.epsilon_kx)+\
                                                                      (orblist[new1].p.py*ham.epsilon_ky-orblist[pos1].p.py*ham.epsilon_ky)*\
                                                                      (orblist[new1].p.py*ham.epsilon_ky - orblist[pos1].p.py*ham.epsilon_ky))\
                                *orblist[new1].fac*orblist[new2].fac*orblist[pos1].fac*orblist[pos2].fac;
                                
                                if (sign_counter % 2 == 0) {
                                    mat_ele.amplitude = amplitude;
                                } else {
                                    mat_ele.amplitude = -amplitude;
                                }
                                
                                
                                if(abs(real(amplitude))>SmallDouble || abs(imag(amplitude))>SmallDouble) {
                                    matrix.push_back(mat_ele);
                                }
                                //dirty block ends
                            }
                            else if (it2.orb1.spin == orblist[pos1].spin && it2.orb2.spin == orblist[pos2].spin)
                            {
                                //dirty block
                                CompactState tempcstate(it.cstate);
                                
                                int sign_counter = 0;
                                
                                //applying the operators
                                for (int i = 0; i < pos1; i++) if (tempcstate[i]) sign_counter++;
                                tempcstate[pos1] = 0;
                                for (int i = 0; i < pos2; i++) if (tempcstate[i]) sign_counter++;
                                tempcstate[pos2] = 0;
                                for (int i = 0; i < new2; i++) if (tempcstate[i]) sign_counter++;
                                tempcstate[new2] = 1;
                                for (int i = 0; i < new1; i++) if (tempcstate[i]) sign_counter++;
                                tempcstate[new1] = 1;
                                
                                int newid = reverse_state_map.at(tempcstate.to_ulong());
                                
                                MatEle mat_ele;
                                mat_ele.bra = it.state_id;
                                mat_ele.ket = newid;
                                double amplitude=ham.norb/8.0*ham.V1*((orblist[new1].p.px*ham.epsilon_kx-orblist[pos1].p.px*ham.epsilon_kx)*\
                                                                      (orblist[new1].p.px*ham.epsilon_kx-orblist[pos1].p.px*ham.epsilon_kx)+\
                                                                      (orblist[new1].p.py*ham.epsilon_ky-orblist[pos1].p.py*ham.epsilon_ky)*\
                                                                      (orblist[new1].p.py*ham.epsilon_ky - orblist[pos1].p.py*ham.epsilon_ky))\
                                *orblist[new1].fac*orblist[new2].fac*orblist[pos1].fac*orblist[pos2].fac+\
                                ham.U*ham.norb/8.0*orblist[pos1].fac*orblist[pos2].fac*orblist[new1].fac*orblist[new2].fac;
                                
                                if (sign_counter % 2 == 0) {
                                    mat_ele.amplitude = amplitude;
                                } else {
                                    mat_ele.amplitude = -amplitude;
                                }
                                
                                
                                if(abs(real(amplitude))>SmallDouble || abs(imag(amplitude))>SmallDouble) {
                                    matrix.push_back(mat_ele);
                                }
                                //dirty block ends
                            }
                        }
                    }
                }
        count++;
        if(count %1000 == 0) cout<<count<<" states finished in part 1"<<endl;
    }
    count = 0;
    for (auto it : state_list)
    {
        int NbrOccupied = 0;
        for (int i = 0; i < ham.norb; i++) {
            if (it.cstate[i]) {
                NbrOccupied++;
            }
        }
        if (NbrOccupied != ham.nele) {
            cout<<"Warning: a state with wrong occupation number came into build_interaction_mat"<<endl;
            continue;
        }
        for (int pos1 = 0; pos1 < ham.norb; pos1++) if(!it.cstate[pos1])
            for (int pos2 = 0; pos2 < ham.norb; pos2++) if(!it.cstate[pos2])
                if(pos2 != pos1)
                {
                    Orbital totalprop = orblist[pos1] + orblist[pos2];
                    vector<OrbPair> possible_pairs = hilbert_space.pair_map.at(totalprop);
                    
                    for (auto it2 : possible_pairs)
                    {
                        int new1 = hilbert_space.reverse_orbital_map.at(it2.orb1);
                        int new2 = hilbert_space.reverse_orbital_map.at(it2.orb2);
                        CompactState tempcstate = it.cstate;
                        tempcstate[pos1] = 1;
                        tempcstate[pos2] = 1;
                        
                        if ((tempcstate[new1] && tempcstate[new2]))
                        {
                            if (it2.orb1.spin == it2.orb2.spin)
                            {
                                //dirty block
                                CompactState tempcstate = it.cstate;
                                int sign_counter = 0;
                                
                                //applying the operators
                                for (int i = 0; i < pos1; i++) if (tempcstate[i]) sign_counter++;
                                tempcstate[pos1] = 1;
                                for (int i = 0; i < pos2; i++) if (tempcstate[i]) sign_counter++;
                                tempcstate[pos2] = 1;
                                for (int i = 0; i < new2; i++) if (tempcstate[i]) sign_counter++;
                                tempcstate[new2] = 0;
                                for (int i = 0; i < new1; i++) if (tempcstate[i]) sign_counter++;
                                tempcstate[new1] = 0;
                                
                                int newid = reverse_state_map.at(tempcstate.to_ulong());
                                
                                MatEle mat_ele;
                                mat_ele.bra = it.state_id;
                                mat_ele.ket = newid;
                                double amplitude=ham.V1*ham.norb/8.0*((orblist[pos1].p.px*ham.epsilon_kx-orblist[new1].p.px*ham.epsilon_kx)*\
                                                                      (orblist[pos1].p.px*ham.epsilon_kx-orblist[new1].p.px*ham.epsilon_kx)+\
                                                                      (orblist[pos1].p.py*ham.epsilon_ky-orblist[new1].p.py*ham.epsilon_ky)*\
                                                                      (orblist[pos1].p.py*ham.epsilon_ky - orblist[new1].p.py*ham.epsilon_ky))\
                                *orblist[pos1].fac*orblist[pos2].fac*orblist[new1].fac*orblist[new2].fac;
                                
                                if (sign_counter % 2 == 0) {
                                    mat_ele.amplitude = amplitude;
                                } else {
                                    mat_ele.amplitude = -amplitude;
                                }
                                
                                if(abs(real(amplitude))>SmallDouble || abs(imag(amplitude))>SmallDouble) {
                                    matrix.push_back(mat_ele);
                                }
                                //dirty block ends
                            }
                            else if (it2.orb1.spin == orblist[pos1].spin && it2.orb2.spin == orblist[pos2].spin)
                            {
                                //dirty block
                                CompactState tempcstate = it.cstate;
                                int sign_counter = 0;
                                
                                //applying the operators
                                for (int i = 0; i < pos1; i++) if (tempcstate[i]) sign_counter++;
                                tempcstate[pos1] = 1;
                                for (int i = 0; i < pos2; i++) if (tempcstate[i]) sign_counter++;
                                tempcstate[pos2] = 1;
                                for (int i = 0; i < new2; i++) if (tempcstate[i]) sign_counter++;
                                tempcstate[new2] = 0;
                                for (int i = 0; i < new1; i++) if (tempcstate[i]) sign_counter++;
                                tempcstate[new1] = 0;
                                
                                int newid = reverse_state_map.at(tempcstate.to_ulong());
                                
                                MatEle mat_ele;
                                mat_ele.bra = it.state_id;
                                mat_ele.ket = newid;
                                double amplitude=ham.V1*ham.norb/8.0*((orblist[pos1].p.px*ham.epsilon_kx-orblist[new1].p.px*ham.epsilon_kx)*\
                                                                      (orblist[pos1].p.px*ham.epsilon_kx-orblist[new1].p.px*ham.epsilon_kx)+\
                                                                      (orblist[pos1].p.py*ham.epsilon_ky-orblist[new1].p.py*ham.epsilon_ky)*\
                                                                      (orblist[pos1].p.py*ham.epsilon_ky - orblist[new1].p.py*ham.epsilon_ky))\
                                *orblist[pos1].fac*orblist[pos2].fac*orblist[new1].fac*orblist[new2].fac+\
                                ham.U*ham.norb/8.0*orblist[pos1].fac*orblist[pos2].fac*orblist[new1].fac*orblist[new2].fac;
                                
                                if (sign_counter % 2 == 0) {
                                    mat_ele.amplitude = amplitude;
                                } else {
                                    mat_ele.amplitude = -amplitude;
                                }
                                
                                if(abs(real(amplitude))>SmallDouble || abs(imag(amplitude))>SmallDouble) {
                                    matrix.push_back(mat_ele);
                                }
                                //dirty block ends
                            }
                        }
                    }
                }
        count++;
        if(count %1000 == 0) cout<<count<<" states finished in part 2"<<endl;
    }

}

