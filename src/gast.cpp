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
 * @file gast.cpp
 * @brief Implementación de la función que calcula el ángulo sidéreo aparente de Greenwich.
 */
//------------------------------------------------------------------------------

#include "../include/gast.hpp"

//------------------------------------------------------------------------------
//  double gast(double Mjd_UT1)
//------------------------------------------------------------------------------
/**
 * @brief Calcula el ángulo sidéreo aparente de Greenwich (GAST).
 *
 * Se obtiene sumando el ángulo sidéreo medio (GMST) con la ecuación de los equinoccios.
 * El resultado se normaliza al intervalo [0, 2π).
 *
 * @param Mjd_UT1 Fecha en tiempo universal UT1 (Modified Julian Date).
 * @return Ángulo sidéreo aparente de Greenwich (radianes).
 */
//------------------------------------------------------------------------------
double gast(double Mjd_UT1){
    return fmod(gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2 * SAT_Const::pi);
}
