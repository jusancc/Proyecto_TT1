//$Header$
//------------------------------------------------------------------------------
//                                   sign_
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file Sign.hpp
 * @brief Declaración de la función `sign_` que devuelve el valor absoluto de `a` con el signo de `b`.
 *
 * Esta función es útil en contextos numéricos y de dinámica orbital donde se necesita 
 * garantizar el signo correcto de una magnitud sin alterar su módulo.
 */
//------------------------------------------------------------------------------

#ifndef _SIGN_
#define _SIGN_

#include <cmath>

/**
 * @brief Devuelve el valor absoluto de `a` con el signo de `b`.
 *
 * @param a Valor cuya magnitud se conservará.
 * @param b Valor cuyo signo se tomará.
 * @return double |a| si b ≥ 0, -|a| si b < 0.
 */
double sign_(double a, double b);

#endif
