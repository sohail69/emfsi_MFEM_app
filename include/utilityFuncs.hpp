#pragma once
#include "mfem.hpp"
#include <functional>

using namespace std;
using namespace mfem;

real_t kDelta(int I, int J){return ( (I==J)? 1.00:0.00);};
