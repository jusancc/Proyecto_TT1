#include "../include/Geodetic.hpp"

tuple<double, double, double> Geodetic(Matrix& r){
    double R_equ = SAT_Const::R_Earth;
    double f = SAT_Const::f_Earth;

    double epsRequ = SAT_Const::eps * R_equ;
    double e2 = f * (2.0 - f);

    double X = r(1);
    double Y = r(2);
    double Z = r(3);

    double rho2 = X*X + Y*Y;

    if (r.norm() == 0){
        cout << "ERROR: invalid input in Geodetic constructor\n";
        double lon = 0.0;
        double lat = 0.0;
        double h   = - SAT_Const::R_Earth;
        return tie(lon, lat, h);
    }

    double dZ = e2 * Z;
    double ZdZ, Nh, SinPhi, N, dZ_new;

    while (1)
    {
        ZdZ =  Z + dZ;
        Nh  =  sqrt ( rho2 + ZdZ*ZdZ ); 
        SinPhi =  ZdZ / Nh;                    
        N =  R_equ / sqrt(1.0-e2*SinPhi*SinPhi);
        dZ_new =  N*e2*SinPhi;
        if ( abs(dZ-dZ_new) < epsRequ ){
            break;
        }
        dZ = dZ_new;
    }

    double lon = atan2(Y,X);
    double lat = atan2(ZdZ,sqrt(rho2));
    double h = Nh - N;
    return tie(lon, lat, h);
}