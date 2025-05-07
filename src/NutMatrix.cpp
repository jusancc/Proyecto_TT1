#include "../include/NutMatrix.hpp"

Matrix& NutMatrix(double Mjd_TT){
    double eps, dpsi, deps;

    eps = MeanObliquity(Mjd_TT);
    tie(dpsi,deps) = NutAngles(Mjd_TT);

    Matrix rx1 = R_x(-eps-deps);
    Matrix rz = R_z(-dpsi);
    Matrix rx2 = R_x(+eps);

    Matrix &NutMat = rx1 * rz * rx2;
    return NutMat;
}