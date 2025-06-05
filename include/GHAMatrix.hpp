//$Header$
//------------------------------------------------------------------------------
//                              GHAMatrix
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file GHAMatrix.hpp
 * @brief Declaración de la función que genera la matriz de rotación para aplicar
 *        el Ángulo Horario de Greenwich (GHA).
 *
 * Esta matriz transforma coordenadas desde el marco inercial al marco rotante
 * de la Tierra, aplicando la rotación sobre el eje Z usando el GAST.
 */
//------------------------------------------------------------------------------

#ifndef _GHA_MATRIX_
#define _GHA_MATRIX_

#include "matrix.hpp"   ///< Clase Matrix personalizada
#include "R_z.hpp"      ///< Rotación alrededor del eje Z
#include "gast.hpp"     ///< Cálculo del ángulo sidéreo aparente de Greenwich

/**
 * @brief Genera la matriz de rotación de Greenwich a partir del GAST.
 * 
 * @param Mjd_UT1 Fecha en tiempo universal UT1 (Modified Julian Date).
 * @return Referencia a una matriz de rotación (3x3) que aplica el GHA.
 */
Matrix& GHAMatrix(double Mjd_UT1);

#endif  // _GHA_MATRIX_
