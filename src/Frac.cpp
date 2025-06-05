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
 * @file Frac.cpp
 * @brief Implementación de la función que devuelve la parte fraccionaria de un número real.
 */
//------------------------------------------------------------------------------

#include "../include/Frac.hpp"

//------------------------------------------------------------------------------
//  double Frac(double x)
//------------------------------------------------------------------------------
/**
 * @brief Devuelve la parte fraccionaria de un número real.
 *
 * Calcula: Frac(x) = x - floor(x), de forma que siempre queda en el intervalo [0, 1).
 * 
 * Útil para operaciones de normalización periódica como ángulos o fases.
 *
 * @param x Número real de entrada.
 * @return Parte fraccionaria de `x` (en el rango [0, 1)).
 */
//------------------------------------------------------------------------------
double Frac(double x){
    double res;
    res = x - floor(x);
    return res;
}
