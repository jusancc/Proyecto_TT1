#ifndef _DEINTEG_
#define _DEINTEG_

#include "matrix.hpp"
#include "SAT_Const.hpp"
#include "Sign.hpp"
#include <cmath>
#include <iostream>
#include <tuple>
#include <cfloat>

Matrix &DEInteg(Matrix& func(double, Matrix&),double &t,double tout,double &relerr,double &abserr,int n_eqn,Matrix &y);

#endif
