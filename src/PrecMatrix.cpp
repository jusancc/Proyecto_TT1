#include "../include/PrecMatrix.hpp"

Matrix& PrecMatrix(double Mjd_1, double Mjd_2){
    double T = (Mjd_1-SAT_Const::MJD_J2000)/36525;
    double dT = (Mjd_2-Mjd_1)/36525;

    double zeta = ( (2306.2181+(1.39656-0.000139*T)*T)+((0.30188-0.000344*T)+0.017998*dT)*dT )*dT/SAT_Const::Arcs;
    double z = zeta + ( (0.79280+0.000411*T)+0.000205*dT)*dT*dT/SAT_Const::Arcs;
    double theta = ( (2004.3109-(0.85330+0.000217*T)*T)-((0.42665+0.000217*T)+0.041833*dT)*dT )*dT/SAT_Const::Arcs;

    Matrix &rz1 = R_z(-z);
    Matrix &ry = R_y(theta);
    Matrix &rz2 = R_z(-zeta);

    Matrix &PrecMat = rz1*ry*rz2;
    return PrecMat;
}