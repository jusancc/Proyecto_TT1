//$Header$
//------------------------------------------------------------------------------
//                                  Frac
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file Frac.hpp
 * @brief Declaración de la función que extrae la parte fraccionaria de un número real.
 *
 * La función devuelve la fracción decimal de un número real, es decir:
 *     Frac(x) = x - floor(x)
 * útil para normalizar ángulos o ciclos periódicos.
 */
//------------------------------------------------------------------------------

#ifndef _FRAC_
#define _FRAC_

#include <cmath>  ///< Función floor()

/**
 * @brief Devuelve la parte fraccionaria de un número real.
 *
 * @param x Número real.
 * @return Parte fraccionaria de x (en [0, 1)).
 */
double Frac(double x);

#endif  // _FRAC_
