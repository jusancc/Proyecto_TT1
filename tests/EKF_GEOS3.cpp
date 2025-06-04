#include "../include/matrix.hpp"
#include "../include/global.hpp"
#include "../include/SAT_Const.hpp"
#include "../include/Mjday.hpp"
#include "../include/Position.hpp"
#include "../include/DEInteg.hpp"
#include "../include/Accel.hpp"
#include "../include/R_z.hpp"
#include "../include/TimeUpdate.hpp"
#include "../include/VarEqn.hpp"
#include "../include/AzElPa.hpp"
#include "../include/MeasUpdate.hpp"
#include "../include/LTC.hpp"
#include <iostream>
#include <cstdio>
#include <cstring>
#include <tuple>

using namespace std;

int main() {

    DE430Coeff(2285, 1020);
    GGM03S(181);
    eop19620101(21413);

    int num_obs = 46;
    Matrix& observations = zeros(num_obs, 4);

    FILE* file = fopen("../data/GEOS3.txt", "r");
    if (!file) {
        cerr << "Error al abrir GEOS3.txt" << endl;
        return 1;
    }

    char buffer[256], temp[12];
    for (int i = 1; i <= num_obs; i++) {
        fgets(buffer, sizeof(buffer), file);

        strncpy(temp, &buffer[0], 4); temp[4] = '\0';
        double YY = atof(temp);
        strncpy(temp, &buffer[5], 2); temp[2] = '\0';
        double MM = atof(temp);
        strncpy(temp, &buffer[8], 2); temp[2] = '\0';
        double DD = atof(temp);
        strncpy(temp, &buffer[12], 2); temp[2] = '\0';
        double hh = atof(temp);
        strncpy(temp, &buffer[15], 2); temp[2] = '\0';
        double mm = atof(temp);
        strncpy(temp, &buffer[18], 6); temp[6] = '\0';
        double ss = atof(temp);
        strncpy(temp, &buffer[25], 8); temp[8] = '\0';
        double az = atof(temp);
        strncpy(temp, &buffer[35], 7); temp[7] = '\0';
        double el = atof(temp);
        strncpy(temp, &buffer[44], 10); temp[10] = '\0';
        double dist = atof(temp);

        observations(i, 1) = Mjday(YY, MM, DD, hh, mm, ss);
        observations(i, 2) = SAT_Const::Rad * az;
        observations(i, 3) = SAT_Const::Rad * el;
        observations(i, 4) = 1000.0 * dist;
    }

    fclose(file);

    const double sig_range = 92.5;
    const double sig_az = 0.0224 * SAT_Const::Rad;
    const double sig_el = 0.0139 * SAT_Const::Rad;

    double latitude = 21.5748 * SAT_Const::Rad;
    double longitude = -158.2706 * SAT_Const::Rad;
    double altitude = 300.2;

    Matrix& Rs = transpose(Position(longitude, latitude, altitude));

    double Mjd1 = observations(1, 1);
    double Mjd2 = observations(9, 1);
    double Mjd3 = observations(18, 1);

    Matrix& r2 = zeros(1);
    Matrix& v2 = zeros(1);
    //tie(r2, v2) = anglesg(observations(1, 2), observations(9, 2), observations(18, 2),observations(1, 3), observations(9, 3), observations(18, 3),Mjd1, Mjd2, Mjd3,transpose(Rs), transpose(Rs), transpose(Rs));

    Matrix& Y0_apr = transpose(union_vector(r2,v2));
    double Mjd0 = Mjday(1995, 1, 29, 2, 38, 0);
    double Mjd_UTC = observations(9, 1);

    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;

    int n_eqn = 6;
    double t = 0;
    double relerr = 1e-13;
    double abserr = 1e-6;
    Matrix& Y = DEInteg(accel, t, -(Mjd_UTC - Mjd0) * 86400.0, relerr, abserr, n_eqn, transpose(Y0_apr));

    Matrix P = eye(n_eqn) * 1e3;
    for (int i = 1; i <= 3; i++) P(i, i) = 1e8;

    Matrix& LT = LTC(longitude, latitude);
    Matrix& yPhi = zeros(42, 1);
    Matrix& Phi = zeros(6, 6);

    double t = 0.0, t_old;
    Matrix& Y_old = zeros(6, 1), U, r, s, dAds, dEds, dAdY, dEdY, dDds, dDdY, K;
    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC, Mjd_TT, Mjd_UT1, theta, Azim, Elev, Dist;

    for (int i = 1; i <= num_obs; i++) {
        t_old = t;
        Y_old = Y;
        Mjd_UTC = observations(i, 1);
        t = (Mjd_UTC - Mjd0) * 86400.0;

        tie(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC) = IERS(eopdata,AuxParam.Mjd_UTC,'l');
        tie(UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC) = timediff(UT1_UTC, TAI_UTC);
        Mjd_TT = Mjd_UTC + TT_UTC / 86400.0;
        Mjd_UT1 = Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;
        AuxParam.Mjd_UTC = Mjd_UTC;
        AuxParam.Mjd_TT = Mjd_TT;

        for (int ii = 1; ii <= 6; ii++) {
            yPhi(ii) = Y_old(ii);
            for (int j = 1; j <= 6; j++)
                yPhi(6 * j + ii) = (ii == j) ? 1.0 : 0.0;
        }


        yPhi = DEInteg(varEqn, t, t - t_old, relerr, abserr, 42, yPhi);
        for (int j = 1; j <= 6; j++)
            Phi.assign_column(j, yPhi.extract_vector(6 * j + 1, 6 * j + 6));

        Y = DEInteg(accel, t, t - t_old, relerr, abserr, 6, Y_old);

        theta = gmst(Mjd_UT1);
        U = R_z(theta);
        r = transpose(Y.extract_vector(1, 3));
        s = LT * (U * r - Rs);

        P = TimeUpdate(P, Phi,0.0);

        tie(Azim, Elev, dAds, dEds) = AzElPa(transpose(s));
        dAdY = union_vector((dAds * LT * U),zeros(3));
        tie(K, Y, P) = MeasUpdate(Y, observations(i, 2), Azim, sig_az, dAdY, P, 6);

        r = transpose(Y.extract_vector(1, 3));
        s = LT * (U * r - Rs);
        tie(Azim, Elev, dAds, dEds) = AzElPa(transpose(s));
        dEdY = union_vector((dEds * LT * U),zeros(3));
        tie(K, Y, P) = MeasUpdate(Y, observations(i, 3), Elev, sig_el, dEdY, P, 6);

        r = transpose(Y.extract_vector(1, 3));
        s = LT * (U * r - Rs);
        Dist = transpose(s).norm();
        dDds = transpose(s / Dist);
        dDdY = union_vector((dDds * LT * U),zeros(3));
        tie(K, Y, P) = MeasUpdate(Y, observations(i, 4), Dist, sig_range, dDdY, P, 6);
    }

    tie(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC) = IERS(eopdata,observations(46, 1), 'l');
    tie(UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC) = timediff(UT1_UTC, TAI_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC / 86400.0;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;

    Matrix& Y0 = DEInteg(accel, t, -(observations(46, 1) - observations(1, 1)) * 86400.0, relerr, abserr, 6, Y);

    Matrix& Y_ref = zeros(6, 1);
    Y_ref(1) = 5753.173e3;
    Y_ref(2) = 2673.361e3;
    Y_ref(3) = 3440.304e3;
    Y_ref(4) = 4324.207;
    Y_ref(5) = -1924.299;
    Y_ref(6) = -5728.216;

    printf("\nError of Position Estimation\n");
    printf("dX%10.1f [m]\n", Y0(1) - Y_ref(1));
    printf("dY%10.1f [m]\n", Y0(2) - Y_ref(2));
    printf("dZ%10.1f [m]\n", Y0(3) - Y_ref(3));
    printf("\nError of Velocity Estimation\n");
    printf("dVx%8.1f [m/s]\n", Y0(4) - Y_ref(4));
    printf("dVy%8.1f [m/s]\n", Y0(5) - Y_ref(5));
    printf("dVz%8.1f [m/s]\n", Y0(6) - Y_ref(6));

    return 0;
}
