#include "../include/EqnEquinox.hpp"

double EqnEquinox(double Mjd_TT){
    double EqE, dpsi, deps;
    tie(dpsi, deps) = NutAngles(Mjd_TT);
    EqE = dpsi * cos(MeanObliquity(Mjd_TT));
    return EqE;
}