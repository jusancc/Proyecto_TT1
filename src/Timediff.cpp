//$Header$
//------------------------------------------------------------------------------
//                                 timediff
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file timediff.cpp
 * @brief Implementación de la función `timediff` para diferencias entre escalas de tiempo.
 *
 * Esta función permite obtener relaciones entre distintas escalas temporales como UT1, TAI, TT, GPS y UTC,
 * a partir de dos entradas: UT1 - UTC y TAI - UTC.
 */
//------------------------------------------------------------------------------

#include "../include/Timediff.hpp"

/**
 * @brief Calcula diferencias entre escalas temporales astronómicas.
 * 
 * @param UT1_UTC Diferencia entre UT1 y UTC [segundos].
 * @param TAI_UTC Diferencia entre TAI y UTC [segundos].
 * @return Tupla con:
 *         - UT1_TAI: UT1 - TAI
 *         - UTC_GPS: UTC - GPS
 *         - UT1_GPS: UT1 - GPS
 *         - TT_UTC:  TT - UTC
 *         - GPS_UTC: GPS - UTC
 */
tuple<double, double, double, double, double> timediff(double UT1_UTC, double TAI_UTC) {
    double TT_TAI  = 32.184;                   // TT - TAI
    double GPS_TAI = -19.0;                    // GPS - TAI
    double TT_GPS  = TT_TAI - GPS_TAI;         // TT - GPS
    double TAI_GPS = -GPS_TAI;                 // TAI - GPS
    double UT1_TAI = UT1_UTC - TAI_UTC;        // UT1 - TAI
    double UTC_TAI = -TAI_UTC;                 // UTC - TAI
    double UTC_GPS = UTC_TAI - GPS_TAI;        // UTC - GPS
    double UT1_GPS = UT1_TAI - GPS_TAI;        // UT1 - GPS
    double TT_UTC  = TT_TAI - UTC_TAI;         // TT - UTC
    double GPS_UTC = GPS_TAI - UTC_TAI;        // GPS - UTC

    return tie(UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
}
