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
    auto cstates = hilbert_space.state_map.at(orb_sector);
    for(int i = 0; i < cstates.size(); i++)
    {
        State temp(cstates[i]);
        temp.state_id = i;
        reverse_state_map[cstates[i].to_ulong()] = i;
        state_list.push_back(temp);
    }
}

void Interaction::build2bodyMatrix()
{
    Hamiltonian &ham = hilbert_space.ham;
    auto &orblist = hilbert_space.orbital_list;
    int count = 0;
    for (auto it : state_list)//go through all states in the list of "states"
    {
        unordered_map<int, complex<double> > temp_mat_ele_list;
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
                        unsigned long new1 = hilbert_space.reverse_orbital_map.at(it2.orb1);
                        unsigned long new2 = hilbert_space.reverse_orbital_map.at(it2.orb2);
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

                                pair<double, double> q = ham.computeK(orblist[new1].p.px - orblist[pos1].p.px, orblist[new1].p.py - orblist[pos1].p.py); // exchange momentum
                                double amplitude = ham.norb/8.0*ham.V1 * (pow(q.first, 2) + pow(q.second, 2)) * orblist[new1].fac*orblist[new2].fac*orblist[pos1].fac*orblist[pos2].fac;
                                
                                if (sign_counter % 2 == 0)
                                    mat_ele.amplitude = amplitude;
                                else
                                    mat_ele.amplitude = -amplitude;
                                
                                
                                if(abs(real(amplitude))>SmallDouble || abs(imag(amplitude))>SmallDouble)
                                {
                                    temp_mat_ele_list[mat_ele.ket] += mat_ele.amplitude;
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
                                
                                pair<double, double> q = ham.computeK(orblist[new1].p.px - orblist[pos1].p.px, orblist[new1].p.py - orblist[pos1].p.py);
                                double amplitude = ham.norb/8.0*(ham.U + ham.V1 * (pow(q.first, 2) + pow(q.second, 2))) * orblist[new1].fac*orblist[new2].fac*orblist[pos1].fac*orblist[pos2].fac;
                                
                                if (sign_counter % 2 == 0)
                                    mat_ele.amplitude = amplitude;
                                else
                                    mat_ele.amplitude = -amplitude;
                                
                                if(abs(real(amplitude))>SmallDouble || abs(imag(amplitude))>SmallDouble)
                                {
                                    temp_mat_ele_list[mat_ele.ket] += mat_ele.amplitude;
                                }
                                //dirty block ends
                            }
                        }
                    }
                }
        for (auto mat_ele_it : temp_mat_ele_list)
        {
            MatEle temp(it.state_id, mat_ele_it.first, mat_ele_it.second);
            matrix.push_back(temp);
        }
        count++;
        if(count %1000 == 0) cout<<count<<" states finished in part 1"<<endl;
    }
    count = 0;
    for (auto it : state_list)
    {
        unordered_map<int, complex<double> > temp_mat_ele_list;
        for (int pos1 = 0; pos1 < ham.norb; pos1++) if(!it.cstate[pos1])
            for (int pos2 = 0; pos2 < ham.norb; pos2++) if(!it.cstate[pos2])
                if(pos2 != pos1)
                {
                    Orbital totalprop = orblist[pos1] + orblist[pos2];
                    vector<OrbPair> possible_pairs = hilbert_space.pair_map.at(totalprop);
                    
                    for (auto it2 : possible_pairs)
                    {
                        unsigned long new1 = hilbert_space.reverse_orbital_map.at(it2.orb1);
                        unsigned long new2 = hilbert_space.reverse_orbital_map.at(it2.orb2);
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
                                pair<double, double> q = ham.computeK(orblist[new1].p.px - orblist[pos1].p.px, orblist[new1].p.py - orblist[pos1].p.py);
                                double amplitude = ham.norb/8.0*ham.V1 * (pow(q.first, 2) + pow(q.second, 2)) * orblist[new1].fac*orblist[new2].fac*orblist[pos1].fac*orblist[pos2].fac;
                                
                                if (sign_counter % 2 == 0)
                                    mat_ele.amplitude = amplitude;
                                else
                                    mat_ele.amplitude = -amplitude;
                                
                                if(abs(real(amplitude))>SmallDouble || abs(imag(amplitude))>SmallDouble) {
                                    temp_mat_ele_list[mat_ele.ket] += mat_ele.amplitude;
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
                                pair<double, double> q = ham.computeK(orblist[new1].p.px - orblist[pos1].p.px, orblist[new1].p.py - orblist[pos1].p.py);
                                double amplitude = ham.norb/8.0*(ham.U + ham.V1 * (pow(q.first, 2) + pow(q.second, 2))) * orblist[new1].fac*orblist[new2].fac*orblist[pos1].fac*orblist[pos2].fac;
                                
                                if (sign_counter % 2 == 0)
                                    mat_ele.amplitude = amplitude;
                                else
                                    mat_ele.amplitude = -amplitude;
                                
                                if(abs(real(amplitude))>SmallDouble || abs(imag(amplitude))>SmallDouble) {
                                    temp_mat_ele_list[mat_ele.ket] += mat_ele.amplitude;
                                }
                                //dirty block ends
                            }
                        }
                    }
                }
        for (auto mat_ele_it : temp_mat_ele_list)
        {
            MatEle temp(it.state_id, mat_ele_it.first, mat_ele_it.second);
            matrix.push_back(temp);
        }
        count++;
        if(count %1000 == 0) cout<<count<<" states finished in part 2"<<endl;
    }

}