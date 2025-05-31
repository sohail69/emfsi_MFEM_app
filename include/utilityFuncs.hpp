#pragma once
#include "mfem.hpp"
#include <functional>

using namespace std;
using namespace mfem;

//
// Kroncecker delta function
//
real_t kDelta(int I, int J){return ( (I==J)? 1.00:0.00);};

//
// Copy one vector to another
//
void copyVec(const Vector & x, Vector & y){y = x;};

//
// Apply Dirchelet Boundary conditions
//

