//$Header$
//------------------------------------------------------------------------------
//                                EccAnom
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file EccAnom.hpp
 * @brief Declaración de la función que resuelve la ecuación de Kepler para obtener
 *        la anomalía excéntrica a partir de la anomalía media y la excentricidad.
 *
 * La función emplea un método iterativo de Newton-Raphson para resolver:
 *     M = E - e·sin(E)
 */
//------------------------------------------------------------------------------

#ifndef _ECC_ANOM_
#define _ECC_ANOM_

#include "SAT_Const.hpp"   ///< Constantes físicas y astronómicas
#include <iostream>        ///< Entrada/salida estándar
#include <cmath>           ///< Funciones matemáticas estándar

using namespace std;

/**
 * @brief Calcula la anomalía excéntrica E a partir de la anomalía media M y la excentricidad e.
 *
 * @param M Anomalía media (rad).
 * @param e Excentricidad (0 ≤ e < 1).
 * @return Anomalía excéntrica (rad).
 */
double EccAnom(double M, double e);

#endif  // _ECC_ANOM_
