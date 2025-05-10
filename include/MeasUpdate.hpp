#ifndef _MEAS_UPDATE_
#define _MEAS_UPDATE_

#include <tuple>
#include "matrix.hpp"

tuple<Matrix&,Matrix&,Matrix&> MeasUpdate(Matrix &x, double z, double g, double s, Matrix &G, Matrix &P, int n);

#endif 
