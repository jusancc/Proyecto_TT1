//$Header$
//------------------------------------------------------------------------------
//                                Legendre
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file Legendre.hpp
 * @brief Prototipo de la función que calcula los polinomios de Legendre y sus derivadas.
 *
 * Esta función genera los valores de los polinomios asociados de Legendre y sus derivadas
 * evaluadas en un ángulo geocéntrico `fi`. Se utiliza comúnmente en el cálculo del potencial
 * gravitacional terrestre para modelos esférico-armónicos.
 */
//------------------------------------------------------------------------------

#ifndef _LEGENDRE_
#define _LEGENDRE_

#include "matrix.hpp"
#include <iostream>
#include <tuple>
#include <cmath>
using namespace std;

/**
 * @brief Calcula los polinomios de Legendre y sus derivadas.
 * 
 * @param n Grado máximo.
 * @param m Orden máximo.
 * @param fi Ángulo geocéntrico (latitud geocéntrica), en radianes.
 * @return Tupla con:
 *         - pnm: Matriz con los polinomios asociados de Legendre (n+1)x(m+1)
 *         - dpnm: Matriz con las derivadas de los polinomios (n+1)x(m+1)
 */
tuple<Matrix&, Matrix&> Legendre(int n, int m, double fi);

#endif
