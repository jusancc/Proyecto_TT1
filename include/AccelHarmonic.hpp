#ifndef _ACCEL_HARMONIC_
#define _ACCEL_HARMONIC_

#include "matrix.hpp"
#include <cmath>
#include "Legendre.hpp"
#include "SAT_Const.hpp"
#include "global.hpp"

Matrix& AccelHarmonic(Matrix &r, Matrix &E, int n_max, int m_max);

#endif 
