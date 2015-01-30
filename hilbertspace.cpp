//
//  hilbertspace.cpp
//  ED_solver
//
//  Created by Bin Xu on 12/19/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//
#include <algorithm>
#include "hilbertspace.h"

void HilbertSpace::buildOrbitalList()
{
    generateHaldaneOrbList();
    computeFacHaldane();
    ham.norb = orbital_list.size();
}

void HilbertSpace::buildReverseOrbitalMap()
{
    for(int i = 0; i < orbital_list.size(); i++)
    {
        reverse_orbital_map[orbital_list[i]] = i;
    }
}

void HilbertSpace::buildPairMap()
{
    for (int id1 = 0; id1 < orbital_list.size(); id1++)
        for (int id2 = 0; id2 < orbital_list.size(); id2++)
        {
            if (id1 != id2)
            {
                pair_map[orbital_list[id1]+orbital_list[id2]].push_back(OrbPair(orbital_list[id1], orbital_list[id2]));
            }
        }
}

void HilbertSpace::buildStateMap()
{
    CompactState cstate;
    for (int i = 0; i < MaxOrbital; i++)
    {
        cstate[i] = 0;
    }
    for (int i = 0; i < ham.nele; i++) cstate[i] = 1;
    Orbital orb_prop = computeMomentumAndSpin(cstate, orbital_list);
    state_map[orb_prop].push_back(cstate);
    while (csNotHighest(cstate, (int)ham.norb, (int)ham.nele))
    {
        csPlusplus(cstate, (int)ham.norb);
        orb_prop = computeMomentumAndSpin(cstate, orbital_list);
        state_map[orb_prop].push_back(cstate);
    }
}




/////////////////////////////////////////////////////////////////////////////////////////
// below are some auxiliary functions, don't touch if you are not sure you understand it
/////////////////////////////////////////////////////////////////////////////////////////



void HilbertSpace::generateHaldaneOrbList()
{
    int norb_max = (int)ham.max_norb;
    size_t pbc_type = ham.pbc_type;
    bool even1, even2;
    
    vector<Orbital> raworblist;
    
    for(int i = -norb_max; i < norb_max; i++)
        for(int j = -norb_max; j< norb_max; j++)
        {
            even1 = (i % 2 == 0);
            even2 = (j % 2 == 0);
            switch (pbc_type)
            {
                case 1:
                    if ((!even1)||(!even2)) continue;
                    break;
                case 2:
                    if (even1 != even2) continue;
                    break;
                case 3:
                    if (even1 || (!even2)) continue;
                    break;
                case 4:
                    if ((!even1) || even2) continue;
                    break;
                case 5:
                    if (even1 || even2) continue;
                    break;
                case 6:
                    if (even1 == even2) continue;
                    break;
                default:
                    cout<<"Error in orbital generation: PBC type not supported."<<endl;
                    abort();
            }
            
            Orbital temporb1(i, j, -1);
            
            temporb1.k = ham.computeK(i, j);
            temporb1.radius2 = pow(temporb1.k.first, 2) + pow(temporb1.k.second, 2);
            raworblist.push_back(temporb1);
        }
    
    
    sort(raworblist.begin(), raworblist.end(),
         [](const Orbital & a, const Orbital & b) -> bool
         {
             return a.radius2 < b.radius2;
         });
    cout<<"=================================================="<<endl;
    
    
    vector<Orbital> semiraworblist;
    
    if (raworblist.size() < norb_max + 1)
    {
        cout << "Error in orbital generation: not enough orbitals generated."<<endl;
        abort();
    }
    
    for(int i = 0; i < norb_max + 1; i++)
    {
        semiraworblist.push_back(raworblist[i]);
    }
    
    bool keep_poping = true;
    
    while(keep_poping)
    {
        cout<<"One orbital omitted..."<<endl;
        last_discarded_r2 = semiraworblist[semiraworblist.size()-1].radius2;
        keep_poping=(semiraworblist[semiraworblist.size() - 1].radius2 - semiraworblist[semiraworblist.size() - 2].radius2 < SmallDouble);
        semiraworblist.pop_back();
    }
    
    if (semiraworblist.size() <= 0)
    {
        cout <<"Error in orbital generation: "<<semiraworblist.size()<<" orbitals are kept."<<endl;
        abort();
    }
    
    for(auto it = semiraworblist.begin(); it != semiraworblist.end(); it++)
    {
        orbital_list.push_back(*it);
        Orbital theotherorbital(*it);
        theotherorbital.spin = 1;
        orbital_list.push_back(theotherorbital);
    }
    
}

void HilbertSpace::computeFacHaldane()
{
    double r2_max;
    switch(orbital_cut_mode)
    {
        case 1:
            r2_max = pow((sqrt(last_discarded_r2) + sqrt(orbital_list.back().radius2))/2.0, 2);
            break;
        case 2:
            r2_max = last_discarded_r2;
            break;
        default:
            cout<<"Orbital cut mode error"<<endl;
            abort();
    }
    for (auto& it : orbital_list)
        it.fac = r2_max - it.radius2;
    double fac_sum = 0;
    for (auto it = orbital_list.begin(); it != orbital_list.end(); it++)
        fac_sum += it->fac/2.0;
    for (auto& it : orbital_list)
        it.fac = sqrt(it.fac/fac_sum);
    while (orbital_list.back().fac < Tiny)
    {
        cout<<"Insignificant orbital omitted: "<<orbital_list.back()<<" fac = "<<orbital_list.back().fac<<endl;
        orbital_list.pop_back();
    }
    
    fac_sum = 0;
    for (auto it = orbital_list.begin(); it != orbital_list.end(); it++)
        fac_sum += pow(it->fac, 2)/2.0;
    fac_sum = sqrt(fac_sum);
    
    for (auto& it : orbital_list)
        it.fac = it.fac/fac_sum;
    
    // Check that we don't screw up computing the fac
    double check_fac_sum = 0;
    for (auto it : orbital_list) if (it.spin == 1)
        check_fac_sum += pow(it.fac, 2);
    cout<<"Orbital list: "<<endl;
    for (auto it : orbital_list)
    {
        cout<<fixed<<it<<"\tfac: "<<it.fac<<"\t"<<"kx = "<<it.k.first<<"\tky = "<<it.k.second<<endl;
        
    }
    cout<<"Sum of fac's: "<<check_fac_sum<<endl;
}

Pxy computeMomentum(CompactState &cstate, vector<Orbital>& orbital_list)
{
    Pxy p(0, 0);
    for (int i = 0; i < orbital_list.size(); i++)
    {
        if (cstate[i])
        {
            p.px += orbital_list[i].p.px;
            p.py += orbital_list[i].p.py;
        }
    }
    return p;
}

Orbital computeMomentumAndSpin(CompactState &cstate, vector<Orbital>& orblist)
{
    Orbital orb(0, 0, 0);
    for (int i = 0; i < orblist.size(); i++)
    {
        if (cstate[i])
        {
            orb.p.px += orblist[i].p.px;
            orb.p.py += orblist[i].p.py;
            orb.spin += orblist[i].spin;
        }
    }
    return orb;
}

void csPlusplus(CompactState& cstate, int norb)
{
    if (!cstate[norb - 1]) //no carry case
    {
        int i = norb - 1;
        while (!cstate[i]) i--;
        cstate[i] = 0, cstate[i + 1] = 1;
    }
    else           //with carry case
    {
        int i = norb - 1;
        int i2 = 0;
        while (cstate[i])
        {
            i--;
            i2++;
        }
        while (!cstate[i]) i--;
        cstate[i] = 0, cstate[i + 1] = 1;
        i += 2;
        for (int j = i; j < i + i2; j++) cstate[j] = 1;
        if (i + i2 < norb)
            for (int j = i + i2; j < norb; j++) cstate[j] = 0;
    }
}

bool csNotHighest(CompactState& cstate, int norb, int nele)
{
    for (int i = norb - 1; i > norb - nele - 1; i--)
    {
        if (!cstate[i]) return true;
    }
    return false;
}
