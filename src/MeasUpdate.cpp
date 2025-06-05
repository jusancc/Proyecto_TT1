//$Header$
//------------------------------------------------------------------------------
//                              MeasUpdate.cpp
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file MeasUpdate.cpp
 * @brief Implementación de la función MeasUpdate para el filtro de Kalman extendido.
 *
 * Esta función aplica la corrección del filtro de Kalman extendido sobre el estado estimado `x`
 * y su covarianza `P`, utilizando una nueva observación `z`. Se calcula la ganancia de Kalman
 * `K`, el residuo entre la observación y su estimación `g`, y se actualizan las variables.
 *
 * La ganancia de Kalman se calcula mediante:
 * \f[
 * K = P G^T (G P G^T + R)^{-1}
 * \f]
 * Donde `R = s²` es la varianza de la observación.
 * 
 * Luego se actualiza:
 * \f[
 * x := x + K(z - g) \\
 * P := (I - K G) P
 * \f]
 */
//------------------------------------------------------------------------------

#include "../include/MeasUpdate.hpp"

/**
 * @brief Realiza la actualización del estado y la covarianza del filtro de Kalman extendido.
 * 
 * @param x Vector de estado estimado antes de la medición.
 * @param z Valor observado de la medición (escalar).
 * @param g Valor estimado de la medición (modelo).
 * @param s Varianza de la medición (desviación estándar).
 * @param G Gradiente de la función de observación (1 x n).
 * @param P Matriz de covarianza del estado estimado.
 * @param n Dimensión del vector de estado.
 * @return Tupla con:
 *         - K: ganancia de Kalman
 *         - x: estado actualizado
 *         - P: nueva matriz de covarianza
 */
tuple<Matrix&, Matrix&, Matrix&> MeasUpdate(Matrix &x, double z, double g, double s, Matrix &G, Matrix &P, int n) {
    double Inv_W = s * s;

    // K = P * G^T * (G * P * G^T + R)^-1
    Matrix& K = transpose((P * transpose(G)) * (G * P * transpose(G) + Inv_W).inv());

    // x = x + K * (z - g)
    x = x + (K * (z - g));

    // P = (I - K * G) * P
    P = (eye(n) - transpose(K) * G) * P;

    return tie(K, x, P);
}
