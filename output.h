//
//  output.h
//  SPT
//
//  Created by Bin Xu on 9/5/14.
//  Copyright (c) 2014 Bin Xu. All rights reserved.
//

#ifndef __SPT__output__
#define __SPT__output__

#include <iostream>
#include "miscellaneous.h"
using namespace std;

ostream& operator<<(ostream& os, Pxy& p);
ostream& operator<<(ostream& os, const Pxy& p);
ostream& operator<<(ostream& os, Orbital& orb);
ostream& operator<<(ostream& os, const Orbital& orb);
ostream& operator<<(ostream& os, CompactState& cstate);
ostream& operator<<(ostream& os, OrbPair& pair);
ostream& operator<<(ostream& os, MatEle & mat_ele);


#endif /* defined(__SPT__output__) */
