#ifndef _GEODETIC_
#define _GEODETIC_

#include "SAT_Const.hpp"
#include <tuple>
#include <iostream>
#include "matrix.hpp"
using namespace std;

tuple<double, double, double> Geodetic(Matrix& r);

#endif
