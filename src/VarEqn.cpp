//$Header$
//------------------------------------------------------------------------------
//                                 VarEqn.cpp
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file VarEqn.cpp
 * @brief Implementación de la función de ecuaciones variacionales para
 * propagación orbital con perturbaciones armónicas.
 *
 * Esta función calcula la derivada del vector de estado extendido (posición,
 * velocidad y matriz Phi de variación del estado), teniendo en cuenta efectos
 * de precesión, nutación, movimiento del polo y rotación de la Tierra, así como
 * el campo gravitatorio armónico.
 */
//------------------------------------------------------------------------------

#include "../include/VarEqn.hpp"
#include <iostream>
#include <iomanip>
using namespace std;

/**
 * @brief Calcula la derivada del vector extendido (estado + matriz de variación).
 * 
 * @param x Tiempo en segundos desde la época inicial (MJD).
 * @param yPhi Vector 42x1 que contiene:
 *             - posición (elementos 1-3)
 *             - velocidad (elementos 4-6)
 *             - matriz de variación del estado Phi (6x6), almacenada por columnas
 * @return Referencia a un vector 42x1 con la derivada temporal del estado extendido.
 */
Matrix& varEqn(double x, Matrix &yPhi) {

    // Parámetros de orientación de la Tierra
    auto [x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC] =
        IERS(eopdata, AuxParam.Mjd_UTC, 'l');

    // Diferencias temporales
    auto [UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);
    double Mjd_UT1 = AuxParam.Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;

    // Transformación de coordenadas: precesión, nutación, polo, GHA
    Matrix &P = PrecMatrix(SAT_Const::MJD_J2000, AuxParam.Mjd_TT + x / 86400.0);
    Matrix &N = NutMatrix(AuxParam.Mjd_TT + x / 86400.0);
    Matrix &T = N * P;
    Matrix &E = PoleMatrix(x_pole, y_pole) * GHAMatrix(Mjd_UT1) * T;

    // Extraer posición y velocidad
    Matrix &r = zeros(3, 1);
    Matrix &v = zeros(3, 1);
    for (int i = 1; i <= 3; i++) {
        r(i, 1) = yPhi(i, 1);
        v(i, 1) = yPhi(i + 3, 1);
    }

    // Reconstruir matriz Phi (6x6)
    Matrix &Phi = zeros(6, 6);
    for (int j = 1; j <= 6; j++) {
        Matrix &col = yPhi.extract_vector(6 * j + 1, 6 * j + 6); // vector fila 1x6
        Phi.assign_column(j, transpose(col));
    }

    // Aceleración y gradiente armónico
    Matrix &a = AccelHarmonic(r, E, AuxParam.n, AuxParam.m);
    Matrix &G = G_AccelHarmonic(r, E, AuxParam.n, AuxParam.m);

    // Jacobiano del sistema
    Matrix &dfdy = zeros(6, 6);
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            dfdy(i, j)       = 0.0;
            dfdy(i + 3, j)   = G(i, j);
            dfdy(i, j + 3)   = (i == j) ? 1.0 : 0.0;
            dfdy(i + 3, j + 3) = 0.0;
        }
    }

    // Derivada de la matriz Phi
    Matrix &Phip = dfdy * Phi;

    // Construir derivada del vector yPhi
    Matrix &yPhip = zeros(42, 1);
    for (int i = 1; i <= 3; i++) {
        yPhip(i, 1)     = v(i, 1);
        yPhip(i + 3, 1) = a(i, 1);
    }
    for (int j = 1; j <= 6; j++) {
        for (int i = 1; i <= 6; i++) {
            yPhip(6 * j + i, 1) = Phip(i, j);
        }
    }

    return yPhip;
}
