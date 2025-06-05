//$Header$
//------------------------------------------------------------------------------
//                           MeanObliquity.cpp
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file MeanObliquity.cpp
 * @brief Implementación de la función MeanObliquity.
 *
 * Esta función calcula la oblicuidad media de la eclíptica (IAU 1980) a partir
 * del Tiempo Terrestre (Mjd_TT), usando una fórmula polinómica que considera
 * el tiempo juliano desde J2000.0 en siglos.
 *
 * Referencia: "The Explanatory Supplement to the Astronomical Almanac"
 */
//------------------------------------------------------------------------------

#include "../include/MeanObliquity.hpp"

/**
 * @brief Calcula la oblicuidad media de la eclíptica.
 * 
 * @param Mjd_TT Fecha juliana modificada en Tiempo Terrestre.
 * @return Oblicuidad media en radianes.
 */
double MeanObliquity(double Mjd_TT){
    double T, MOblq;

    // Tiempo juliano en siglos desde J2000.0
    T = (Mjd_TT - SAT_Const::MJD_J2000) / 36525.0;

    // Oblicuidad media en arcosegundos, convertida a radianes
    MOblq = SAT_Const::Rad * (
        84381.448 / 3600.0 -
        (46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600.0
    );

    return MOblq;
}
