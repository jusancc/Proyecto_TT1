#ifndef _VAR_EQN_
#define _VAR_EQN_

#include "matrix.hpp"
#include "IERS.hpp"
#include "global.hpp"
#include "SAT_Const.hpp"
#include "Timediff.hpp"
#include "PrecMatrix.hpp"
#include "NutMatrix.hpp"
#include "PoleMatrix.hpp"
#include "GHAMatrix.hpp"
#include "AccelHarmonic.hpp"
#include "G_AccelHarmonic.hpp"

Matrix& varEqn(double x, Matrix &yPhi);

#endif
