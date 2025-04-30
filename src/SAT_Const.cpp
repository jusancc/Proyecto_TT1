#include "../include/SAT_Const.hpp"

// Mathematical constants
const double SAT_Const::pi       = 3.14159265358979324;
const double SAT_Const::pi2      = 2 * SAT_Const::pi;
const double SAT_Const::Rad      = SAT_Const::pi / 180;
const double SAT_Const::Deg      = 180 / SAT_Const::pi;
const double SAT_Const::Arcs     = 3600 * 180 / SAT_Const::pi;

// General constants
const double SAT_Const::MJD_J2000 = 51544.5;
const double SAT_Const::T_B1950   = -0.500002108;
const double SAT_Const::c_light   = 299792458.000000000;
const double SAT_Const::AU        = 149597870700.000000;

// Physical parameters
const double SAT_Const::R_Earth   = 6378.1363e3;
const double SAT_Const::f_Earth   = 1.0 / 298.257223563;
const double SAT_Const::R_Sun     = 696000e3;
const double SAT_Const::R_Moon    = 1738e3;
const double SAT_Const::omega_Earth = 15.04106717866910 / 3600 * SAT_Const::Rad;

// Gravitational coefficients
const double SAT_Const::GM_Earth    = 398600.435436e9;
const double SAT_Const::GM_Sun      = 132712440041.939400e9;
const double SAT_Const::GM_Moon     = SAT_Const::GM_Earth / 81.30056907419062;
const double SAT_Const::GM_Mercury  = 22031.780000e9;
const double SAT_Const::GM_Venus    = 324858.592000e9;
const double SAT_Const::GM_Mars     = 42828.375214e9;
const double SAT_Const::GM_Jupiter  = 126712764.800000e9;
const double SAT_Const::GM_Saturn   = 37940585.200000e9;
const double SAT_Const::GM_Uranus   = 5794548.600000e9;
const double SAT_Const::GM_Neptune  = 6836527.100580e9;
const double SAT_Const::GM_Pluto    = 977.0000000000009e9;

// Solar radiation
const double SAT_Const::P_Sol = 1367 / SAT_Const::c_light;