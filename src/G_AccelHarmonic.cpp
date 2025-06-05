//$Header$
//------------------------------------------------------------------------------
//                          G_AccelHarmonic
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file G_AccelHarmonic.cpp
 * @brief Implementación de la función que calcula el gradiente de la aceleración
 *        armónica respecto a la posición.
 */
//------------------------------------------------------------------------------

#include "../include/G_AccelHarmonic.hpp"
#include <iostream>

//------------------------------------------------------------------------------
//  Matrix& G_AccelHarmonic(Matrix &r, Matrix &U, double n_max, double m_max)
//------------------------------------------------------------------------------
/**
 * @brief Calcula numéricamente el gradiente de la aceleración armónica ∂a/∂r.
 *
 * Esta función evalúa la derivada parcial de la aceleración producida por el
 * campo gravitacional armónico respecto a la posición usando diferencias finitas
 * centradas. Se perturba cada componente de `r` una cantidad `d/2` y se calcula
 * la diferencia entre aceleraciones.
 *
 * @param r Vector de posición (3x1).
 * @param U Matriz de transformación (3x3) del sistema inercial al cuerpo.
 * @param n_max Orden máximo del modelo armónico.
 * @param m_max Grado máximo del modelo armónico.
 * @return Referencia a una matriz (3x3) con el gradiente ∂a/∂r.
 *
 * @note Esta aproximación es numérica y puede ser costosa para modelos grandes.
 */
//------------------------------------------------------------------------------
Matrix& G_AccelHarmonic(Matrix &r, Matrix &U, double n_max, double m_max){
    double d = 1.0;  // Paso para diferencias finitas

    Matrix &G = zeros(3, 3);   // Matriz resultante del gradiente
    Matrix &dr = zeros(3, 1);  // Vector de perturbación

    for (int i = 1; i <= 3; i++) {
        dr = zeros(3, 1);
        dr(i) = d;

        Matrix &a_plus = AccelHarmonic(r + dr / 2, U, n_max, m_max);
        Matrix &a_minus = AccelHarmonic(r - dr / 2, U, n_max, m_max);

        Matrix &da = a_plus - a_minus;

        G.assign_column(i, da / d);
    }

    return G;
}
