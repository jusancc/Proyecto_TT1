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
 * @file Sign.cpp
 * @brief Implementación de la función `sign_`, que devuelve el valor absoluto de `a` con el signo de `b`.
 *
 * Esta función es comúnmente usada en algoritmos numéricos para garantizar que
 * una magnitud siga la dirección (signo) de una variable de referencia.
 */
//------------------------------------------------------------------------------

#include "../include/Sign.hpp"

/**
 * @brief Devuelve el valor absoluto de `a` con el signo de `b`.
 *
 * @param a Valor cuya magnitud se conservará.
 * @param b Valor cuyo signo se tomará.
 * @return double |a| si b ≥ 0, -|a| si b < 0.
 */
double sign_(double a, double b){
    double result;
    if (b >= 0.0)
    {
        result = abs(a);
    } else {
        result = -abs(a);
    }

    return result;
}
