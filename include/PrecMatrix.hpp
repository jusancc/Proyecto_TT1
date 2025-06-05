//$Header$
//------------------------------------------------------------------------------
//                                PrecMatrix
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file PrecMatrix.hpp
 * @brief Declaración de la función que calcula la matriz de precesión entre dos fechas.
 *
 * Esta función genera la matriz de rotación debida a la precesión de la Tierra entre
 * dos épocas dadas en tiempo juliano modificado (MJD), según el modelo IAU 1976.
 */
//------------------------------------------------------------------------------

#ifndef _PREC_MATRIX_
#define _PREC_MATRIX_

#include "matrix.hpp"
#include "R_y.hpp"
#include "R_z.hpp"
#include "SAT_Const.hpp"

/**
 * @brief Calcula la matriz de precesión entre dos fechas dadas.
 *
 * La precesión es una rotación lenta del eje de la Tierra que afecta a las coordenadas
 * astronómicas a lo largo del tiempo. Esta función construye la matriz de transformación
 * que lleva coordenadas referidas a la época Mjd_1 a la época Mjd_2.
 *
 * @param Mjd_1 Fecha inicial en tiempo juliano modificado.
 * @param Mjd_2 Fecha final en tiempo juliano modificado.
 * @return Matriz de rotación 3x3 que representa la precesión entre ambas fechas.
 */
Matrix& PrecMatrix(double Mjd_1, double Mjd_2);

#endif
