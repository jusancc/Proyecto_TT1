//$Header$
//------------------------------------------------------------------------------
//                                  gmst
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file GMST.hpp
 * @brief Declaración de la función que calcula el Tiempo Sidéreo Medio de Greenwich (GMST).
 *
 * Esta función permite obtener la rotación de la Tierra respecto al marco inercial,
 * necesaria para transformar coordenadas entre marcos de referencia.
 */
//------------------------------------------------------------------------------

#ifndef _GMST_
#define _GMST_

#include "SAT_Const.hpp"  ///< Constantes astronómicas
#include "Frac.hpp"       ///< Función para obtener la parte fraccionaria
#include <cmath>          ///< Funciones matemáticas estándar

/**
 * @brief Calcula el Tiempo Sidéreo Medio de Greenwich en radianes.
 *
 * @param Mjd_UT1 Fecha en tiempo universal UT1 (Modified Julian Date).
 * @return GMST en radianes, en el rango [0, 2π).
 */
double gmst(double Mjd_UT1);

#endif  // _GMST_
