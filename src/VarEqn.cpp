#include "../include/VarEqn.hpp"

Matrix& varEqn(double x, Matrix &yPhi){
    auto[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,AuxParam.Mjd_UTC,'l');
    auto[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
    double Mjd_UT1 = AuxParam.Mjd_TT+(UT1_UTC-TT_UTC)/86400;

    Matrix &P = PrecMatrix(SAT_Const::MJD_J2000,AuxParam.Mjd_TT+x/86400);
    Matrix &N = NutMatrix(AuxParam.Mjd_TT+x/86400);
    Matrix &T = N * P;
    Matrix &E = PoleMatrix(x_pole,y_pole)*GHAMatrix(Mjd_UT1)*T;

    Matrix &r = yPhi.extract_vector(1,3);
    Matrix &v = yPhi.extract_vector(4,6);
    Matrix &Phi = zeros(6);

    for (int i = 1; i <= 6; i++)
    {
        Phi.assign_column(i, yPhi.extract_vector(6*i+1,6*i+6));
    }

    Matrix &a = AccelHarmonic(r,E,AuxParam.n,AuxParam.m);
    Matrix &G = G_AccelHarmonic(r,E,AuxParam.n,AuxParam.m);

    Matrix &yPhip = zeros(42,1);
    Matrix &dfdy = zeros(6);
    
    for (int i = 1; i <= 3; i++)
    {
        for (int j = 1; j <= 3; j++)
        {
            dfdy(i,j) = 0.0;
            dfdy(i+3,j) = G(i,j);
            if (i == j)
            {
                dfdy(i,j+3) = 1;
            } else {
                dfdy(i,j+3) = 0;
            }

            dfdy(i+3,j+3) = 0.0;
        }
    }

    Matrix &Phip = dfdy*Phi;
    
    for (int i = 1; i <= 3; i++)
    {
        yPhip(i) = v(i);
        yPhip(i+3) = a(i);
    }
    for (int i = 1; i <= 6; i++)
    {
        for (int j = 1; j <= 6; j++)
        {
            yPhip(6*j+i) = Phip(i,j);
        }
    }
    
    return yPhip;
}