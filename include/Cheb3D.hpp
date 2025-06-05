//$Header$
//------------------------------------------------------------------------------
//                                Cheb3D
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file Cheb3D.hpp
 * @brief Declaración de la función que evalúa una serie de interpolación de
 *        Chebyshev en tres dimensiones.
 *
 * Esta función evalúa una serie de polinomios de Chebyshev para obtener una
 * posición 3D interpolada en el intervalo [Ta, Tb] a partir de coeficientes dados.
 */
//------------------------------------------------------------------------------

#ifndef _CHEB_3D_
#define _CHEB_3D_

#include "matrix.hpp"     ///< Clase Matrix personalizada para álgebra lineal
#include <iostream>       ///< Entrada/salida estándar

using namespace std;

/**
 * @brief Evalúa una serie de interpolación de Chebyshev en 3D.
 * 
 * @param t Tiempo en el que se evalúa la serie.
 * @param N Número de coeficientes (orden del polinomio + 1).
 * @param Ta Límite inferior del intervalo de interpolación.
 * @param Tb Límite superior del intervalo de interpolación.
 * @param Cx Coeficientes Chebyshev en X (1xN o Nx1).
 * @param Cy Coeficientes Chebyshev en Y (1xN o Nx1).
 * @param Cz Coeficientes Chebyshev en Z (1xN o Nx1).
 * @return Referencia a un vector (1x3) con la posición interpolada [x, y, z].
 */
Matrix& Cheb3D(double t, int N, double Ta, double Tb, Matrix Cx, Matrix Cy, Matrix Cz);

#endif  // _CHEB_3D_
