#ifndef _IERS_
#define _IERS_

#include "SAT_Const.hpp"
#include "matrix.hpp"
#include "iostream"
#include <tuple>
using namespace std;

tuple<double,double,double,double,double,double,double,double,double> IERS(Matrix eop, double Mjd_UTC, char interp='n');

#endif
