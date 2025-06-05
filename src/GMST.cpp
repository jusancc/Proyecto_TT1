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
 * @file GMST.cpp
 * @brief Implementación de la función que calcula el Tiempo Sidéreo Medio de Greenwich.
 */
//------------------------------------------------------------------------------

#include "../include/GMST.hpp"

//------------------------------------------------------------------------------
//  double gmst(double Mjd_UT1)
//------------------------------------------------------------------------------
/**
 * @brief Calcula el Tiempo Sidéreo Medio de Greenwich (GMST) en radianes.
 *
 * Esta función sigue la fórmula del IAU 1982 para calcular el GMST, que representa
 * el ángulo entre el meridiano de Greenwich y el punto Aries, medido en el plano ecuatorial.
 *
 * @param Mjd_UT1 Tiempo universal UT1 (Modified Julian Date).
 * @return GMST en radianes, normalizado al intervalo [0, 2π).
 */
//------------------------------------------------------------------------------
double gmst(double Mjd_UT1){
    double Secs = 86400;                 // Segundos por día
    double MJD_J2000 = 51544.5;          // Época J2000.0

    double Mjd_0 = floor(Mjd_UT1);       // Día entero
    double UT1 = Secs * (Mjd_UT1 - Mjd_0);      // Segundos desde medianoche
    double T_0 = (Mjd_0 - MJD_J2000) / 36525.0; // Tiempo en siglos desde J2000
    double T = (Mjd_UT1 - MJD_J2000) / 36525.0; // Tiempo continuo en siglos

    // Fórmula para GMST en segundos (IAU 1982)
    double gmst = 24110.54841
                + 8640184.812866 * T_0
                + 1.002737909350795 * UT1
                + (0.093104 - 6.2e-6 * T) * T * T;

    // Conversión a radianes y normalización
    double gmstime = 2 * SAT_Const::pi * Frac(gmst / Secs);

    return gmstime;
}
