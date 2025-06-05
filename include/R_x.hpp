//$Header$
//------------------------------------------------------------------------------
//                                   R_x
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file R_x.cpp
 * @brief Implementación de la matriz de rotación respecto al eje X.
 *
 * Esta función genera una matriz de rotación 3x3 para un ángulo dado,
 * que rota un vector en sentido antihorario alrededor del eje X.
 */
//------------------------------------------------------------------------------


#include <iostream>
#include <iomanip>
#include "../include/Matrix.hpp"

/**
 * @brief Genera la matriz de rotación respecto al eje X.
 * 
 * @param alpha Ángulo de rotación en radianes.
 * @return Matriz 3x3 que representa la rotación en torno al eje X.
 * 
 * La forma de la matriz es:
 * 
 *     [ 1     0          0     ]
 *     [ 0   cos(α)   sin(α) ]
 *     [ 0  -sin(α)   cos(α) ]
 */
Matrix& R_x(double alpha);