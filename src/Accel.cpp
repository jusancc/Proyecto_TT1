#include "../include/Accel.hpp"

Matrix& accel(double x, Matrix &Y){
    auto[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,AuxParam.Mjd_UTC + x / 86400, 'l');
    auto[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC, TAI_UTC);

    double Mjd_UT1 = AuxParam.Mjd_UTC + (x / 86400) + UT1_UTC / 86400;
    double Mjd_TT = AuxParam.Mjd_UTC + (x / 86400) + TT_UTC / 86400;

    Matrix &P = PrecMatrix(SAT_Const::MJD_J2000, Mjd_TT);
    Matrix &N = NutMatrix(Mjd_TT);

    Matrix &T = N * P;


    Matrix &Pole = PoleMatrix(x_pole, y_pole);
    Matrix &GHA = GHAMatrix(Mjd_UT1);

    Matrix &EG = Pole * GHA;

    Matrix &E = EG * T;

    double MJD_TDB = Mjday_TDB(Mjd_TT);
    auto[r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun] = JPL_Eph_DE430(MJD_TDB);

    Matrix &pos = Y.extract_vector(1, 3);

    Matrix &a = AccelHarmonic(transpose(pos), E, AuxParam.n, AuxParam.m);

    if (AuxParam.sun) {
        a = a + transpose(AccelPointMass(pos, r_Sun, SAT_Const::GM_Sun));
    }
    if (AuxParam.moon) {
        a = a + transpose(AccelPointMass(pos, r_Moon, SAT_Const::GM_Moon));
    }
    if (AuxParam.planets) {
        a = a + transpose(AccelPointMass(pos, r_Mercury, SAT_Const::GM_Mercury));
        a = a + transpose(AccelPointMass(pos, r_Venus, SAT_Const::GM_Venus));
        a = a + transpose(AccelPointMass(pos, r_Mars, SAT_Const::GM_Mars));
        a = a + transpose(AccelPointMass(pos, r_Jupiter, SAT_Const::GM_Jupiter));
        a = a + transpose(AccelPointMass(pos, r_Saturn, SAT_Const::GM_Saturn));
        a = a + transpose(AccelPointMass(pos, r_Uranus, SAT_Const::GM_Uranus));
        a = a + transpose(AccelPointMass(pos, r_Neptune, SAT_Const::GM_Neptune));
        a = a + transpose(AccelPointMass(pos, r_Pluto, SAT_Const::GM_Pluto));
    }

    Matrix &vel = Y.extract_vector(4, 6);

    Matrix &dY = union_vector(vel, transpose(a));

    return dY;
}
