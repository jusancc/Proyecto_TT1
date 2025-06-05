//$Header$
//------------------------------------------------------------------------------
//                                 Position
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file Position.hpp
 * @brief Prototipo de la función que convierte coordenadas geodésicas a posición en coordenadas cartesianas.
 *
 * Esta función transforma las coordenadas geodésicas (longitud, latitud, altitud)
 * a coordenadas cartesianas ECEF (Earth-Centered, Earth-Fixed) usando el modelo elipsoidal de la Tierra.
 */
//------------------------------------------------------------------------------

#ifndef _POSITION_
#define _POSITION_

#include "matrix.hpp"
#include "SAT_Const.hpp"
#include <cmath>

/**
 * @brief Calcula la posición en coordenadas cartesianas ECEF a partir de coordenadas geodésicas.
 *
 * @param lon Longitud [rad].
 * @param lat Latitud geodésica [rad].
 * @param h Altitud sobre el elipsoide [m].
 * @return Referencia a una matriz columna (3x1) con las coordenadas X, Y, Z en el sistema ECEF.
 */
Matrix& Position(double lon, double lat, double h);

#endif
