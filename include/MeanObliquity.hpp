//$Header$
//------------------------------------------------------------------------------
//                           MeanObliquity.hpp
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file MeanObliquity.hpp
 * @brief Declaración de la función MeanObliquity.
 *
 * Esta función calcula la oblicuidad media de la eclíptica (IAU 1980) en función
 * de la fecha juliana modificada en Tiempo Terrestre (Mjd_TT). Es usada
 * en transformaciones astronómicas como parte del modelo de precesión-nutación.
 */
//------------------------------------------------------------------------------

#ifndef _MEAN_OBLIQUITY_
#define _MEAN_OBLIQUITY_

#include "SAT_Const.hpp"
#include <cmath> 

/**
 * @brief Calcula la oblicuidad media de la eclíptica.
 * 
 * @param Mjd_TT Fecha juliana modificada en Tiempo Terrestre.
 * @return Oblicuidad media (en radianes).
 */
double MeanObliquity(double Mjd_TT);

#endif
