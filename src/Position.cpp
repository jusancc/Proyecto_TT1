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
 * @file Position.cpp
 * @brief Implementación de la función que convierte coordenadas geodésicas a cartesianas ECEF.
 *
 * Esta función transforma una ubicación dada por longitud, latitud y altitud (geodésicas)
 * en coordenadas cartesianas del sistema ECEF (Earth-Centered, Earth-Fixed), utilizando
 * el modelo elipsoidal del WGS84.
 */
//------------------------------------------------------------------------------

#include "../include/Position.hpp"

/**
 * @brief Transforma coordenadas geodésicas (lon, lat, h) en cartesianas ECEF (X, Y, Z).
 *
 * @param lon Longitud [rad].
 * @param lat Latitud geodésica [rad].
 * @param h Altitud sobre el elipsoide [m].
 * @return Matriz columna (3x1) con coordenadas cartesianas [X, Y, Z] en metros.
 */
Matrix& Position(double lon, double lat, double h){
    double R_equ = SAT_Const::R_Earth;   // Radio ecuatorial [m]
    double f = SAT_Const::f_Earth;       // Aplanamiento de la Tierra

    double e2 = f * (2.0 - f);           // Excentricidad cuadrada
    double CosLat = cos(lat);
    double SinLat = sin(lat);

    double N = R_equ / sqrt(1.0 - e2 * SinLat * SinLat);  // Radio de curvatura en la vertical

    Matrix &r = zeros(1, 3);  // Vector de posición [X, Y, Z]
    r(1, 1) = (N + h) * CosLat * cos(lon);
    r(1, 2) = (N + h) * CosLat * sin(lon);
    r(1, 3) = ((1.0 - e2) * N + h) * SinLat;

    return r;
}
