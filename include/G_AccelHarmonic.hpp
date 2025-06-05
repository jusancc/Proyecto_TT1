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
 * @file G_AccelHarmonic.hpp
 * @brief Declaración de la función que calcula el gradiente de la aceleración
 *        producida por el campo gravitacional armónico.
 *
 * Esta función proporciona la derivada parcial de la aceleración respecto a la posición,
 * lo cual es útil en algoritmos de integración, propagación de errores o filtros de Kalman.
 */
//------------------------------------------------------------------------------

#ifndef _G_ACCEL_HARMONIC_
#define _G_ACCEL_HARMONIC_

#include "matrix.hpp"          ///< Clase Matrix personalizada
#include "AccelHarmonic.hpp"   ///< Aceleración por armónicos esféricos

/**
 * @brief Calcula el gradiente de la aceleración armónica respecto a la posición.
 * 
 * @param r Vector de posición (3x1).
 * @param U Matriz de rotación de sistema inercial a cuerpo fijo (3x3).
 * @param n_max Orden máximo del desarrollo armónico.
 * @param m_max Grado máximo del desarrollo armónico.
 * @return Referencia a una matriz (3x3) con el gradiente ∂a/∂r.
 */
Matrix& G_AccelHarmonic(Matrix &r, Matrix &U, double n_max, double m_max);

#endif  // _G_ACCEL_HARMONIC_
