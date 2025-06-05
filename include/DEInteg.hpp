//$Header$
//------------------------------------------------------------------------------
//                                DEInteg
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file DEInteg.hpp
 * @brief Declaración de la función que implementa el integrador multistep de
 *        Shampine y Gordon para resolver sistemas de EDOs con control de paso.
 *
 * Este integrador se usa para resolver ecuaciones diferenciales ordinarias
 * mediante un método predictor-corrector con orden variable y control de errores.
 */
//------------------------------------------------------------------------------

#ifndef _DEINTEG_
#define _DEINTEG_

#include "matrix.hpp"        ///< Clase Matrix personalizada
#include "SAT_Const.hpp"     ///< Constantes físicas y astronómicas
#include "Sign.hpp"          ///< Función sign_ para signo real
#include <cmath>             ///< Funciones matemáticas estándar
#include <iostream>          ///< Entrada/salida estándar
#include <tuple>             ///< Uso de tuplas en funciones internas
#include <cfloat>            ///< Límites numéricos y precisión

/**
 * @brief Integrador multistep con control de paso y errores de Shampine y Gordon.
 *
 * @param func Función que define el sistema EDO a resolver: dy/dt = f(t, y).
 * @param t Tiempo actual (modificado al finalizar la integración).
 * @param tout Tiempo objetivo al que se desea integrar.
 * @param relerr Tolerancia relativa permitida.
 * @param abserr Tolerancia absoluta permitida.
 * @param n_eqn Número de ecuaciones del sistema.
 * @param y Vector de estado (modificado al final con la solución).
 * @return Referencia a un `Matrix` con el estado en `tout`.
 */
Matrix &DEInteg(Matrix& func(double, Matrix&), double &t, double tout, double &relerr, double &abserr, int n_eqn, Matrix &y);

#endif  // _DEINTEG_
