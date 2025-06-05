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
 * @file Geodetic.hpp
 * @brief Declaración de la función que convierte coordenadas cartesianas
 *        en coordenadas geodésicas (latitud, longitud, altura).
 *
 * Esta conversión asume una Tierra elipsoidal definida por el modelo WGS.
 */
//------------------------------------------------------------------------------

#ifndef _GEODETIC_
#define _GEODETIC_

#include "SAT_Const.hpp"   ///< Constantes físicas y del elipsoide terrestre
#include <tuple>           ///< Tupla para retornar múltiples valores
#include <iostream>        ///< Entrada/salida (opcional para debugging)
#include "matrix.hpp"      ///< Clase Matrix para manejo de vectores

using namespace std;

/**
 * @brief Convierte un vector de posición cartesiano a coordenadas geodésicas.
 *
 * @param r Vector de posición (3x1) en coordenadas cartesianas (ECEF) [km].
 * @return Tupla con:
 *         - Latitud geodésica (rad)
 *         - Longitud (rad)
 *         - Altura sobre el elipsoide (km)
 */
tuple<double, double, double> Geodetic(Matrix& r);

#endif  // _GEODETIC_
