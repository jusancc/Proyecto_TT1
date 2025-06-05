//$Header$
//------------------------------------------------------------------------------
//                                   R_y
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file R_y.hpp
 * @brief Declaración de la función R_y que genera una matriz de rotación sobre el eje Y.
 *
 * Esta función retorna una matriz de rotación 3x3 correspondiente a un ángulo
 * dado en radianes alrededor del eje Y. Es útil para transformaciones espaciales
 * y rotaciones en dinámica orbital.
 */
//------------------------------------------------------------------------------

#ifndef _R_Y_
#define _R_Y_

#include "matrix.hpp"
#include <cmath>

/**
 * @brief Devuelve la matriz de rotación sobre el eje Y.
 * @param angle Ángulo de rotación en radianes.
 * @return Referencia a una matriz 3x3 que representa la rotación.
 */
Matrix& R_y(double angle);

#endif
