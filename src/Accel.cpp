//$Header$
//------------------------------------------------------------------------------
//                                accel
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file accel.cpp
 * @brief Implementación de la función que calcula la aceleración total que actúa 
 *        sobre un satélite, considerando perturbaciones armónicas y gravitatorias
 *        de cuerpos celestes usando efemérides JPL DE430.
 */
//------------------------------------------------------------------------------

#include "../include/Accel.hpp"

//------------------------------------------------------------------------------
//  Matrix& accel(double x, Matrix &Y)
//------------------------------------------------------------------------------
/**
 * @brief Calcula la aceleración total sobre un cuerpo a partir de su estado y tiempo.
 * 
 * Esta función computa las fuerzas que actúan sobre un satélite, incluyendo:
 * - Aceleración por armónicos esféricos del campo gravitacional terrestre.
 * - Atracción gravitatoria del Sol, la Luna y planetas (opcional según configuración).
 * 
 * El sistema de referencia se transforma aplicando matrices de precesión, nutación,
 * movimiento del polo y ángulo horario de Greenwich.
 *
 * @param x Tiempo transcurrido (en segundos) desde la época de referencia (MJD_UTC).
 * @param Y Vector de estado (6x1) con posición y velocidad cartesianas [km, km/s].
 * @return Referencia a un vector de derivadas del estado (6x1): [vx,vy,vz,ax,ay,az].
 */
//------------------------------------------------------------------------------
Matrix& accel(double x, Matrix &Y){
    // Obtener parámetros de interpolación IERS para la época actual
    auto [x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC] =
        IERS(eopdata, AuxParam.Mjd_UTC + x / 86400, 'l');
    
    // Diferencias entre escalas de tiempo
    auto [UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);

    double Mjd_UT1 = AuxParam.Mjd_UTC + (x / 86400) + UT1_UTC / 86400;
    double Mjd_TT = AuxParam.Mjd_UTC + (x / 86400) + TT_UTC / 86400;

    // Matrices de transformación: precesión, nutación, polo, ángulo horario
    Matrix &P = PrecMatrix(SAT_Const::MJD_J2000, Mjd_TT);
    Matrix &N = NutMatrix(Mjd_TT);
    Matrix &T = N * P;

    Matrix &Pole = PoleMatrix(x_pole, y_pole);
    Matrix &GHA = GHAMatrix(Mjd_UT1);
    Matrix &EG = Pole * GHA;
    Matrix &E = EG * T;

    // Cálculo de efemérides JPL para cuerpos del sistema solar
    double MJD_TDB = Mjday_TDB(Mjd_TT);
    auto [r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn,
          r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun] = JPL_Eph_DE430(MJD_TDB);

    Matrix &pos = Y.extract_vector(1, 3);  // Vector posición

    // Aceleración por armónicos esféricos
    Matrix &a = AccelHarmonic(transpose(pos), E, AuxParam.n, AuxParam.m);

    // Aceleración por atracción de cuerpos puntuales (opcional)
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

    Matrix &vel = Y.extract_vector(4, 6); // Vector velocidad

    // Concatenar velocidad y aceleración como derivadas del estado
    Matrix &dY = union_vector(vel, transpose(a));

    return dY;
}
