//$Header$
//------------------------------------------------------------------------------
//                                 NutAngles
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file NutAngles.hpp
 * @brief Declaración de la función que calcula los ángulos de nutación.
 *
 * Esta función devuelve los ángulos de nutación en longitud (`dpsi`) y 
 * en oblicuidad (`deps`) para una fecha dada en Tiempo Terrestre (TT).
 * Se basa en el modelo de la IAU para las perturbaciones del movimiento
 * de la Tierra, necesarias para transformaciones precisas entre marcos 
 * de referencia astronómicos.
 */
//------------------------------------------------------------------------------

#ifndef _NUT_ANGLES_
#define _NUT_ANGLES_

#include "SAT_Const.hpp"
#include <iostream>
#include <cmath>
#include <tuple>
using namespace std;

/**
 * @brief Calcula los ángulos de nutación para una fecha en Tiempo Terrestre.
 *
 * @param Mjd_TT Fecha en MJD (Modified Julian Date) en Tiempo Terrestre.
 * @return Tupla con:
 *         - dpsi: Nutación en longitud [rad]
 *         - deps: Nutación en oblicuidad [rad]
 */
tuple<double, double> NutAngles(double Mjd_TT);

#endif
