//$Header$
//------------------------------------------------------------------------------
//                                  angl
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file angl.hpp
 * @brief Declaración de la función que calcula el ángulo entre dos vectores.
 *
 * Esta función devuelve el ángulo (en radianes) entre dos vectores 3D usando
 * el producto escalar y la norma de los vectores.
 */
//------------------------------------------------------------------------------

#ifndef _ANGL_
#define _ANGL_

#include "matrix.hpp"  ///< Clase Matrix personalizada para operaciones vectoriales

/**
 * @brief Calcula el ángulo entre dos vectores tridimensionales.
 *
 * @param vec1 Primer vector (3x1 o 1x3).
 * @param vec2 Segundo vector (3x1 o 1x3).
 * @return Ángulo entre ambos vectores en radianes.
 */
double angl(Matrix& vec1, Matrix& vec2);

#endif  // _ANGL_
