//$Header$
//------------------------------------------------------------------------------
//                                  gast
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file gast.hpp
 * @brief Declaración de la función que calcula el ángulo sidéreo aparente de Greenwich (GAST).
 *
 * El GAST se obtiene sumando el GMST (ángulo sidéreo medio) y la ecuación de los equinoccios.
 */
//------------------------------------------------------------------------------

#ifndef _GAST_
#define _GAST_

#include "GMST.hpp"          ///< Ángulo sidéreo medio de Greenwich
#include "EqnEquinox.hpp"    ///< Ecuación de los equinoccios
#include "SAT_Const.hpp"     ///< Constantes astronómicas
#include <cmath>             ///< Función fmod()

/**
 * @brief Calcula el ángulo sidéreo aparente de Greenwich (GAST).
 *
 * @param Mjd_UT1 Fecha en tiempo universal UT1 (Modified Julian Date).
 * @return Ángulo sidéreo aparente de Greenwich en radianes (en [0, 2π)).
 */
double gast(double Mjd_UT1);

#endif  // _GAST_
