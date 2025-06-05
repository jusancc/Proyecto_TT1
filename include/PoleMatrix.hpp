//$Header$
//------------------------------------------------------------------------------
//                               PoleMatrix
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file PoleMatrix.hpp
 * @brief Declaración de la función que calcula la matriz de transformación por movimiento polar.
 *
 * Esta función construye una matriz de rotación 3x3 que representa el desplazamiento de los polos 
 * (xp, yp) con respecto al eje de rotación de la Tierra, y se utiliza para corregir el sistema
 * de coordenadas terrestres.
 */
//------------------------------------------------------------------------------

#ifndef _POLE_MATRIX_
#define _POLE_MATRIX_

#include "R_y.hpp"
#include "R_x.hpp"
#include "matrix.hpp"

/**
 * @brief Calcula la matriz de rotación por efecto del movimiento polar.
 * 
 * @param xp Desplazamiento del polo en el eje X [rad].
 * @param yp Desplazamiento del polo en el eje Y [rad].
 * @return Referencia a una matriz 3x3 que representa la transformación por movimiento polar.
 */
Matrix& PoleMatrix(double xp, double yp);

#endif
