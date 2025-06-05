//$Header$
//------------------------------------------------------------------------------
//                                accel
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file accel.hpp
 * @brief Declaración de la función que calcula la aceleración total sobre un cuerpo
 *        celeste, considerando perturbaciones armónicas y de cuerpos puntuales.
 *
 * Esta cabecera define la interfaz de la función `accel`, que computa la aceleración
 * resultante sobre un cuerpo debido a efectos gravitatorios como armónicos de la Tierra
 * y atracciones de cuerpos puntuales usando efemérides y modelos de precesión/nutación.
 */
//------------------------------------------------------------------------------

#ifndef _ACCEL_
#define _ACCEL_

#include "global.hpp"            ///< Parámetros globales
#include "SAT_Const.hpp"         ///< Constantes astronómicas
#include "matrix.hpp"            ///< Clase Matrix personalizada
#include "IERS.hpp"              ///< Interpolación de parámetros IERS
#include "Timediff.hpp"          ///< Diferencias de tiempo (TT, UTC, etc.)
#include "PrecMatrix.hpp"        ///< Matriz de precesión
#include "NutMatrix.hpp"         ///< Matriz de nutación
#include "PoleMatrix.hpp"        ///< Matriz de movimiento del polo
#include "GHAMatrix.hpp"         ///< Matriz de ángulo horario de Greenwich
#include "Mjday.hpp"             ///< Conversión a Modified Julian Date
#include "AccelHarmonic.hpp"     ///< Aceleración por armónicos esféricos
#include "AccelPointMass.hpp"    ///< Aceleración por cuerpos puntuales
#include "Mjday_TDB.hpp"         ///< Tiempo TDB (dínamico baricéntrico)
#include "JPL_Eph_DE430.hpp"     ///< Efemérides JPL DE430

/**
 * @brief Calcula la aceleración total sobre un cuerpo a partir del estado y tiempo.
 * 
 * @param x Tiempo transcurrido en segundos desde un instante de referencia.
 * @param Y Vector de estado (6x1): posición (x,y,z) y velocidad (vx,vy,vz).
 * @return Referencia a un objeto Matrix que contiene la aceleración total (3x1).
 */
Matrix& accel(double x, Matrix &Y);

#endif  // _ACCEL_
