//$Header$
//------------------------------------------------------------------------------
//                                  LTC
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file LTC.hpp
 * @brief Declaración de la función que genera la matriz de transformación del sistema geodésico local.
 *
 * La matriz LTC (Local Tangent Coordinate) transforma de coordenadas cartesianas 
 * ECEF a un sistema local tangente al elipsoide terrestre en una ubicación dada
 * por latitud y longitud.
 */
//------------------------------------------------------------------------------

#ifndef _LTC_
#define _LTC_

#include "R_y.hpp"       ///< Rotación alrededor del eje Y
#include "R_z.hpp"       ///< Rotación alrededor del eje Z
#include "matrix.hpp"    ///< Clase de matrices personalizada

/**
 * @brief Genera la matriz de transformación LTC (Local Tangent Coordinate).
 *
 * @param lon Longitud geodésica [rad]
 * @param lat Latitud geodésica [rad]
 * @return Matriz de transformación (3x3) del sistema ECEF al sistema local tangente (topocéntrico).
 */
Matrix& LTC(double lon, double lat);

#endif  // _LTC_
