//$Header$
//------------------------------------------------------------------------------
//                                NutMatrix
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file NutMatrix.hpp
 * @brief Declaración de la función que calcula la matriz de nutación.
 *
 * Esta función construye la matriz de transformación de nutación que corrige
 * la orientación del eje de la Tierra debida a los efectos periódicos de la
 * nutación en longitud y oblicuidad. La matriz se genera utilizando la oblicuidad
 * media de la eclíptica y los ángulos de nutación calculados para el tiempo dado.
 */
//------------------------------------------------------------------------------

#ifndef _NUT_MATRIX_
#define _NUT_MATRIX_

#include "MeanObliquity.hpp"  ///< Cálculo de la oblicuidad media
#include "NutAngles.hpp"      ///< Cálculo de ángulos de nutación
#include "R_x.hpp"            ///< Rotación alrededor del eje X
#include "R_z.hpp"            ///< Rotación alrededor del eje Z
#include "matrix.hpp"         ///< Clase Matrix

/**
 * @brief Construye la matriz de nutación para un instante dado (Mjd_TT).
 *
 * @param Mjd_TT Tiempo Terrestre en días julianos modificados.
 * @return Matriz de rotación 3x3 que representa la nutación.
 */
Matrix& NutMatrix(double Mjd_TT);

#endif
