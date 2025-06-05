//$Header$
//------------------------------------------------------------------------------
//                              MeasUpdate.hpp
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file MeasUpdate.hpp
 * @brief Declaración de la función MeasUpdate, utilizada en el filtro de Kalman extendido.
 *
 * Esta función realiza la actualización del estado estimado, la matriz de covarianza,
 * y el residuo innovador en un filtro de Kalman extendido, a partir de una medición
 * escalar y su predicción.
 *
 * La medición puede ser de tipo rango, velocidad radial, o ángulo, y está modelada
 * mediante una función no lineal linealizada alrededor del estado actual.
 */
//------------------------------------------------------------------------------

#ifndef _MEAS_UPDATE_
#define _MEAS_UPDATE_

#include <tuple>
#include "matrix.hpp"

/**
 * @brief Actualiza el estado estimado, la matriz de covarianza y el residuo en un filtro de Kalman.
 * 
 * @param x Vector de estado estimado antes de la medición.
 * @param z Valor observado de la medición (escalar).
 * @param g Valor estimado de la medición (modelo).
 * @param s Varianza de la medición.
 * @param G Gradiente (jacobiano) de la medición con respecto al estado (1xn).
 * @param P Matriz de covarianza del estado estimado.
 * @param n Dimensión del vector de estado.
 * @return Tupla con:
 *         - x_new: estado actualizado
 *         - P_new: nueva matriz de covarianza
 *         - K: vector ganancia de Kalman
 */
tuple<Matrix&, Matrix&, Matrix&> MeasUpdate(Matrix &x, double z, double g, double s, Matrix &G, Matrix &P, int n);

#endif
