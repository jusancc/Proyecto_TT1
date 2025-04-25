#include "../include/MeanObliquity.hpp"


namespace Const {

    const double MJD_J2000 = 51544.5;  
    const double Rad = M_PI / 180.0;  
}

double MeanObliquity(double Mjd_TT){
    double T, MOblq;

    T = (Mjd_TT-Const::MJD_J2000)/36525;
    MOblq = Const::Rad*(84381.448/3600-(46.8150+(0.00059-0.001813*T)*T)*T/3600);

    return MOblq;
}