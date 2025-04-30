#ifndef _SAT_CONST_
#define _SAT_CONST_

class SAT_Const {
public:
    // Mathematical constants
    static const double pi;
    static const double pi2;
    static const double Rad;
    static const double Deg;
    static const double Arcs;

    // General
    static const double MJD_J2000;
    static const double T_B1950;
    static const double c_light;
    static const double AU;

    // Physical parameters
    static const double R_Earth;
    static const double f_Earth;
    static const double R_Sun;
    static const double R_Moon;
    static const double omega_Earth;

    // Gravitational coefficients
    static const double GM_Earth;
    static const double GM_Sun;
    static const double GM_Moon;
    static const double GM_Mercury;
    static const double GM_Venus;
    static const double GM_Mars;
    static const double GM_Jupiter;
    static const double GM_Saturn;
    static const double GM_Uranus;
    static const double GM_Neptune;
    static const double GM_Pluto;

    // Solar radiation
    static const double P_Sol;
};

#endif // SAT_CONST_HPP