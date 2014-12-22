//
//  output.cpp
//  SPT
//
//  Created by Bin Xu on 9/5/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#include "output.h"

////////////////////////////////////////////////////////////////////////////////
// Output functions
////////////////////////////////////////////////////////////////////////////////

ostream& operator<<(ostream& os, Pxy& p)
{ os << "Px = " << p.px << "\tPy = " << p.py; return os; }
ostream& operator<<(ostream& os, const Pxy& p)
{ os << "Px = " << p.px << "\tPy = " << p.py; return os; }
ostream& operator<<(ostream& os, Orbital& orb)
{
    os << "Px = " << orb.p.px << "\tPy = " << orb.p.py\
    << "\t Spin: " << (int)orb.spin;
    return os;
}
ostream& operator<<(ostream& os, const Orbital& orb)
{
    os << "Px = " << orb.p.px << "\tPy = " << orb.p.py\
    << "\t Spin is " <<(int)orb.spin;
    return os;
}
ostream& operator<<(ostream& os, CompactState& cstate)
{
    
    for (int i = 0; i < cstate.size(); i++)
    {
        os<<cstate[i];
    }
	return os;
}

ostream& operator<<(ostream& os, OrbPair& pair)
{
    const char up[] = "\xe2\x88\x86";
    const char down[] = "\xe2\x88\x87";
    os << "Orbital "<<" : ( " << pair.orb1.p.px << ", " << pair.orb1.p.py <<\
    " )\t" << (pair.orb1.spin == 1 ? up : down) << "\tOrbital "<<\
    " : ( " << pair.orb2.p.px << ", " << pair.orb2.p.py << " )\t" <<\
    (pair.orb2.spin == 1 ? up : down);
    return os;
}
ostream& operator<<(ostream & os, MatEle & mat_ele)
{
    os << "<" << mat_ele.bra << "|H|" << mat_ele.ket << "> = "\
    << mat_ele.amplitude;
    return os;
}