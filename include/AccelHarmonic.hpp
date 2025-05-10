#ifndef _ACCEL_HARMONIC_
#define _ACCEL_HARMONIC_

#include "matrix.hpp"
#include <cmath>
#include "Legendre.hpp"
#include "SAT_Const.hpp"
#include "global.hpp"

Matrix& AccelHarmonic(Matrix &r, Matrix &E, double n_max, double m_max);

#endif 
