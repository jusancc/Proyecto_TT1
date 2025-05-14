#ifndef _ACCEL_
#define _ACCEL_

#include "global.hpp"
#include "SAT_Const.hpp"
#include "matrix.hpp"
#include "IERS.hpp"
#include "Timediff.hpp"
#include "PrecMatrix.hpp"
#include "NutMatrix.hpp"
#include "PoleMatrix.hpp"
#include "GHAMatrix.hpp"
#include "Mjday.hpp"
#include "AccelHarmonic.hpp"
#include "AccelPointMass.hpp"
#include "Mjday_TDB.hpp"
#include "JPL_Eph_DE430.hpp"

Matrix& accel(double x, Matrix &Y);

#endif
