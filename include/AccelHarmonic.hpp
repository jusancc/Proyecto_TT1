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
 * @file AccelHarmonic.hpp
 * @brief Declaración de la función que calcula la aceleración debida al
 *        campo gravitacional terrestre modelado mediante armónicos esféricos.
 *
 * Esta cabecera define la interfaz de la función `AccelHarmonic`, que utiliza
 * coeficientes de armónicos hasta orden y grado especificado para calcular
 * perturbaciones gravitacionales sobre un satélite.
 */
//------------------------------------------------------------------------------

#ifndef _ACCEL_HARMONIC_
#define _ACCEL_HARMONIC_

#include "matrix.hpp"         ///< Clase Matrix personalizada para álgebra lineal
#include <cmath>              ///< Funciones matemáticas estándar
#include "Legendre.hpp"       ///< Funciones de Legendre asociadas
#include "SAT_Const.hpp"      ///< Constantes físicas y astronómicas
#include "global.hpp"         ///< Parámetros globales y auxiliares

/**
 * @brief Calcula la aceleración perturbativa debido a armónicos esféricos terrestres.
 * 
 * @param r Vector de posición (3x1) en coordenadas cartesianas.
 * @param E Matriz de transformación de coordenadas (MOD -> ECEF).
 * @param n_max Orden máximo del desarrollo armónico.
 * @param m_max Grado máximo del desarrollo armónico.
 * @return Referencia a un vector (3x1) con la aceleración [km/s²].
 */
Matrix& AccelHarmonic(Matrix &r, Matrix &E, int n_max, int m_max);

#endif  // _ACCEL_HARMONIC_
