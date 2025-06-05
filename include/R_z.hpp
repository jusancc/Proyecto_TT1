
//$Header$
//------------------------------------------------------------------------------
//                                   R_z
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file R_z.cpp
 * @brief Implementación de la función R_z para generar una matriz de rotación sobre el eje Z.
 *
 * Esta función construye una matriz 3x3 de rotación alrededor del eje Z para un
 * ángulo especificado en radianes. Se utiliza en transformaciones espaciales,
 * especialmente en contextos de mecánica celeste y navegación orbital.
 */
//------------------------------------------------------------------------------

#ifndef _R_Z_
#define _R_Z_

#include "matrix.hpp"
#include <cmath>


/**
 * @brief Genera la matriz de rotación sobre el eje Z.
 * 
 * @param angle Ángulo de rotación en radianes.
 * @return Referencia a una matriz 3x3 correspondiente a la rotación.
 */
Matrix& R_z(double angle);


#endif 
