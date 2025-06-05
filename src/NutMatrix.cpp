//$Header$
//------------------------------------------------------------------------------
//                                NutMatrix
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file NutMatrix.cpp
 * @brief Implementación de la función que calcula la matriz de nutación.
 *
 * Esta función construye la matriz de transformación de nutación, que representa
 * los pequeños movimientos periódicos del eje de la Tierra en respuesta a las
 * fuerzas gravitacionales de la Luna y el Sol. Se basa en la oblicuidad media
 * de la eclíptica y en los ángulos de nutación en longitud y oblicuidad.
 */
//------------------------------------------------------------------------------

#include "../include/NutMatrix.hpp"

/**
 * @brief Calcula la matriz de rotación de nutación para un instante dado.
 *
 * @param Mjd_TT Tiempo Terrestre en días julianos modificados.
 * @return Matriz 3x3 de rotación asociada a la nutación.
 *
 * La matriz resultante es el producto de tres rotaciones:
 * 1. Una rotación alrededor del eje X por el ángulo `-eps - deps`,
 * 2. Una rotación alrededor del eje Z por `-dpsi`,
 * 3. Una rotación alrededor del eje X por `+eps`.
 */
Matrix& NutMatrix(double Mjd_TT) {
    double eps, dpsi, deps;

    // Oblicuidad media de la eclíptica
    eps = MeanObliquity(Mjd_TT);

    // Ángulos de nutación en longitud y oblicuidad
    tie(dpsi, deps) = NutAngles(Mjd_TT);

    // Rotaciones intermedias para componer la matriz de nutación
    Matrix rx1 = R_x(-eps - deps);
    Matrix rz = R_z(-dpsi);
    Matrix rx2 = R_x(+eps);

    // Matriz de nutación final
    Matrix &NutMat = rx1 * rz * rx2;

    return NutMat;
}
