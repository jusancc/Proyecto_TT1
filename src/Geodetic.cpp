//$Header$
//------------------------------------------------------------------------------
//                                Geodetic
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file Geodetic.cpp
 * @brief Implementación de la función que convierte coordenadas cartesianas
 *        (ECEF) en coordenadas geodésicas (latitud, longitud y altura).
 */
//------------------------------------------------------------------------------

#include "../include/Geodetic.hpp"

//------------------------------------------------------------------------------
//  tuple<double, double, double> Geodetic(Matrix& r)
//------------------------------------------------------------------------------
/**
 * @brief Convierte coordenadas cartesianas a geodésicas usando el modelo elipsoidal.
 *
 * A partir de un vector de posición `r` en coordenadas cartesianas (ECEF), esta función
 * calcula la latitud geodésica, la longitud y la altura sobre el elipsoide WGS84
 * mediante un algoritmo iterativo basado en la altitud y la excentricidad del elipsoide.
 *
 * @param r Vector de posición (3x1) en ECEF [km].
 * @return Tupla con:
 *         - Longitud (rad)
 *         - Latitud geodésica (rad)
 *         - Altura sobre el elipsoide (km)
 *
 * @note Si el vector `r` es el nulo, retorna lat = lon = 0, h = -R_Tierra.
 */
//------------------------------------------------------------------------------
tuple<double, double, double> Geodetic(Matrix& r){
    double R_equ = SAT_Const::R_Earth;         // radio ecuatorial [km]
    double f = SAT_Const::f_Earth;             // achatamiento del elipsoide

    double epsRequ = SAT_Const::eps * R_equ;   // tolerancia para iteración
    double e2 = f * (2.0 - f);                 // excentricidad cuadrada

    double X = r(1);
    double Y = r(2);
    double Z = r(3);

    double rho2 = X*X + Y*Y;

    // Validación de entrada
    if (r.norm() == 0){
        cout << "ERROR: invalid input in Geodetic constructor\n";
        double lon = 0.0;
        double lat = 0.0;
        double h   = -R_equ;
        return tie(lon, lat, h);
    }

    double dZ = e2 * Z;
    double ZdZ, Nh, SinPhi, N, dZ_new;

    // Iteración para ajustar la latitud (lat) y altura (h)
    while (true) {
        ZdZ = Z + dZ;
        Nh = sqrt(rho2 + ZdZ * ZdZ);
        SinPhi = ZdZ / Nh;
        N = R_equ / sqrt(1.0 - e2 * SinPhi * SinPhi);
        dZ_new = N * e2 * SinPhi;

        if (fabs(dZ - dZ_new) < epsRequ) {
            break;
        }

        dZ = dZ_new;
    }

    double lon = atan2(Y, X);                  // Longitud geodésica
    double lat = atan2(ZdZ, sqrt(rho2));       // Latitud geodésica
    double h = Nh - N;                         // Altura sobre el elipsoide

    return tie(lon, lat, h);
}
