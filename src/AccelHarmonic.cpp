//$Header$
//------------------------------------------------------------------------------
//                            AccelHarmonic
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file AccelHarmonic.cpp
 * @brief Implementación de la función que calcula la aceleración gravitacional
 *        perturbativa debida a los armónicos esféricos del campo gravitacional terrestre.
 */
//------------------------------------------------------------------------------

#include "../include/AccelHarmonic.hpp"

//------------------------------------------------------------------------------
//  Matrix& AccelHarmonic(Matrix &r, Matrix &E, int n_max, int m_max)
//------------------------------------------------------------------------------
/**
 * @brief Calcula la aceleración debido a los armónicos esféricos del campo terrestre.
 *
 * Esta función transforma la posición del satélite al sistema de referencia del cuerpo,
 * computa los potenciales asociados a los armónicos esféricos de orden y grado dados,
 * y deriva la aceleración resultante en el sistema inercial.
 * 
 * @param r Vector de posición (3x1) del satélite en coordenadas cartesianas [km].
 * @param E Matriz de transformación (de inercial a cuerpo fijo).
 * @param n_max Orden máximo del desarrollo armónico.
 * @param m_max Grado máximo del desarrollo armónico.
 * @return Referencia a una matriz (3x1) con la aceleración resultante [km/s²].
 */
//------------------------------------------------------------------------------
Matrix& AccelHarmonic(Matrix &r, Matrix &E, int n_max, int m_max) {
    double r_ref = 6378.1363e3;            // radio de referencia en metros
    double gm = 398600.4415e9;             // constante gravitacional de la Tierra [m³/s²]

    Matrix &r_bf = E * r;                  // posición en sistema del cuerpo fijo
    double d = transpose(r_bf).norm();    // distancia radial
    double latgc = asin(r_bf(3, 1) / d);   // latitud geocéntrica
    double lon = atan2(r_bf(2, 1), r_bf(1, 1)); // longitud

    // Cálculo de funciones de Legendre y derivadas
    auto [pnm, dpnm] = Legendre(n_max, m_max, latgc);

    // Inicialización de derivadas del potencial
    double dUdr = 0.0;
    double dUdlatgc = 0.0;
    double dUdlon = 0.0;

    // Bucle sobre orden y grado
    for (int n = 0; n <= n_max; n++) {
        double b1 = (-gm / (d * d)) * pow(r_ref / d, n) * (n + 1);
        double b2 = (gm / d) * pow(r_ref / d, n);
        double b3 = (gm / d) * pow(r_ref / d, n);

        double q1 = 0.0, q2 = 0.0, q3 = 0.0;

        for (int m = 0; m <= m_max; m++) {
            double cos_ml = cos(m * lon);
            double sin_ml = sin(m * lon);

            double term1 = pnm(n + 1, m + 1) * (Cnm(n + 1, m + 1) * cos_ml + Snm(n + 1, m + 1) * sin_ml);
            double term2 = dpnm(n + 1, m + 1) * (Cnm(n + 1, m + 1) * cos_ml + Snm(n + 1, m + 1) * sin_ml);
            double term3 = m * pnm(n + 1, m + 1) * (Snm(n + 1, m + 1) * cos_ml - Cnm(n + 1, m + 1) * sin_ml);

            q1 += term1;
            q2 += term2;
            q3 += term3;
        }

        dUdr += q1 * b1;
        dUdlatgc += q2 * b2;
        dUdlon += q3 * b3;
    }

    // Componentes en el sistema del cuerpo
    double r2xy = pow(r_bf(1, 1), 2) + pow(r_bf(2, 1), 2);
    double ax = (1.0 / d * dUdr - r_bf(3, 1) / (d * d * sqrt(r2xy)) * dUdlatgc) * r_bf(1, 1) - (1.0 / r2xy * dUdlon) * r_bf(2, 1);
    double ay = (1.0 / d * dUdr - r_bf(3, 1) / (d * d * sqrt(r2xy)) * dUdlatgc) * r_bf(2, 1) + (1.0 / r2xy * dUdlon) * r_bf(1, 1);
    double az = 1.0 / d * dUdr * r_bf(3, 1) + sqrt(r2xy) / (d * d) * dUdlatgc;

    // Vector de aceleración en sistema del cuerpo
    Matrix a_bf(3, 1);
    a_bf(1, 1) = ax;
    a_bf(2, 1) = ay;
    a_bf(3, 1) = az;

    // Transformar de vuelta al sistema inercial
    Matrix &a = transpose(E) * a_bf;

    return a;
}
