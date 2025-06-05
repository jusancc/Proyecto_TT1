//$Header$
//------------------------------------------------------------------------------
//                                AzElPa
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file AzElPa.cpp
 * @brief Implementación de la función que calcula el azimut, la elevación y las
 *        derivadas parciales respecto a la posición.
 */
//------------------------------------------------------------------------------

#include "../include/AzElPa.hpp"

//------------------------------------------------------------------------------
//  tuple<double, double, Matrix&, Matrix&> AzElPa(Matrix &s)
//------------------------------------------------------------------------------
/**
 * @brief Calcula los ángulos de observación y las derivadas parciales respecto a la posición.
 *
 * A partir de un vector de posición `s` en coordenadas cartesianas, esta función devuelve:
 * - El azimut (medido desde el norte hacia el este).
 * - La elevación (ángulo sobre el horizonte).
 * - El vector derivada parcial del azimut respecto a `s`.
 * - El vector derivada parcial de la elevación respecto a `s`.
 *
 * @param s Vector de posición en coordenadas cartesianas (1x3) [km].
 * @return Tupla con:
 *         - Azimut [rad]
 *         - Elevación [rad]
 *         - dAz/ds (1x3)
 *         - dEl/ds (1x3)
 */
//------------------------------------------------------------------------------
tuple<double,double,Matrix&,Matrix&> AzElPa(Matrix &s){
    double pi2 = SAT_Const::pi2;
    double rho = sqrt(s(1,1)*s(1,1) + s(1,2)*s(1,2));  // distancia proyectada en plano XY

    // Azimut en radianes (desde el norte hacia el este)
    double Az = atan2(s(1,1), s(1,2));
    if (Az < 0.0) {
        Az += pi2;
    }

    // Elevación sobre el horizonte
    double El = atan(s(1,3) / rho);

    // Derivadas parciales
    Matrix &dAds = zeros(1, 3);  // dAz/ds
    Matrix &dEds = zeros(1, 3);  // dEl/ds

    dAds(1,1) = s(1,2) / (rho * rho);
    dAds(1,2) = -s(1,1) / (rho * rho);
    dAds(1,3) = 0.0;

    dEds(1,1) = -s(1,1) * s(1,3) / rho;
    dEds(1,2) = -s(1,2) * s(1,3) / rho;
    dEds(1,3) = rho;

    dEds = dEds / s.dot(s);  // normalización con la norma cuadrada

    return tie(Az, El, dAds, dEds);
}
