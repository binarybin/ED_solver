//
//  miscellaneous.h
//  ED_solver
//
//  Created by Bin Xu on 12/19/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#ifndef ED_solver_miscellaneous_h
#define ED_solver_miscellaneous_h

#include <bitset>
#include <complex>
#include "constants.h"
using namespace std;

struct Pxy
{
    int px, py;
    bool operator== (const Pxy &rhs) const
    { return px == rhs.px && py == rhs.py; }
    Pxy operator+ (const Pxy &rhs) const
    { Pxy temp(px + rhs.px, py + rhs.py); return temp; }
    
    Pxy(){}
    Pxy(int o_px, int o_py){ px = o_px, py = o_py; }
};

struct Orbital
{
    char spin;	//In spin 1/2 case, one sets 1 for up spin and 2 for down spin.
    Pxy p;
    double fac;
    double radius2;
    
    bool operator== (const Orbital &rhs) const
    { return p == rhs.p && spin == rhs.spin; }
    
    Orbital operator+ (const Orbital &rhs) const
    {Orbital temp(p.px+rhs.p.px, p.py+rhs.p.py, spin+rhs.spin); return temp;}
    
    Orbital(){fac = 0; radius2 = 0;}
    
    Orbital(int px, int py, char spin)
    { p.px = px, p.py = py, this->spin = spin; fac = 0; radius2 = 0;}
    
    
};

struct OrbPair
{
    Orbital orb1, orb2;
    char spin;
    Pxy p;
    
    OrbPair(){}
    OrbPair(const Orbital& orba, const Orbital& orbb)
    {
        orb1 = orba, orb2 = orbb; p = orb1.p + orb2.p;
        spin = orb1.spin + orb2.spin;
    }
};

struct OrbTriplet
{
    Orbital orb1, orb2, orb3;
    char spins;
    Pxy p;
    
    OrbTriplet(){}
    OrbTriplet(const Orbital& orba, const Orbital& orbb, const Orbital& orbc)
    : orb1(orba), orb2(orbb), orb3(orbc)
    {
        spins = orb1.spin + orb2.spin + orb3.spin;
    }
};

typedef bitset<MaxOrbital> CompactState;

struct State
{
    CompactState cstate;
    unsigned state_id;
    State(CompactState &o_cstate): cstate(o_cstate), state_id(0){}
    State(): state_id(0){}
};

////////////////////////////////////////////////////////////////////////////////
// Hashers
////////////////////////////////////////////////////////////////////////////////
struct Pxy_hasher
{
    size_t operator()(const Pxy& p) const
    {
        return (hash<int>()(p.px) ^ (hash<int>()(p.py) << 1));
    }
};

struct Orbital_hasher
{
    size_t operator()(const Orbital& orb) const
    {
        return (hash<char>()(orb.spin)^((hash<int>()(orb.p.px) ^ \
                                         (hash<int>()(orb.p.py) << 1))<<1));
    }
};

struct MatEle
{
    unsigned bra, ket;
    complex<double> amplitude;
    MatEle()
    {
        bra = 0, ket = 0, amplitude = 0;
    }
    MatEle(int b, int k, complex<double> amp) : bra(b), ket(k), amplitude(amp) {}
};


#endif
