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
 * @file timediff.hpp
 * @brief Prototipo de la función `timediff`, que calcula diferencias entre escalas temporales.
 *
 * Esta función convierte entre distintas escalas de tiempo astronómicas:
 * UT1, UTC, TAI, TT y GPS, a partir de las diferencias UT1-UTC y TAI-UTC.
 */
//------------------------------------------------------------------------------

#ifndef _TIMEDIFF_
#define _TIMEDIFF_

#include <iostream>
#include <tuple>
using namespace std;

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
tuple<double,double,double,double,double> timediff(double UT1_UTC, double TAI_UTC);

#endif
