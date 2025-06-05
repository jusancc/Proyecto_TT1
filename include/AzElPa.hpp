//$Header$
//------------------------------------------------------------------------------
//                                AzElPa
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file AzElPa.hpp
 * @brief Declaración de la función que calcula el azimut, la elevación y los
 *        vectores de la base local (topocéntrica) a partir de un vector posición.
 *
 * Esta función devuelve los ángulos de observación (azimut y elevación) así como
 * los vectores unitarios del sistema local (sur, este, cenit) para observaciones terrestres.
 */
//------------------------------------------------------------------------------

#ifndef _AZ_EL_PA_
#define _AZ_EL_PA_

#include "matrix.hpp"        ///< Clase Matrix personalizada
#include "SAT_Const.hpp"     ///< Constantes astronómicas
#include <iostream>          ///< Entrada/salida estándar
#include <tuple>             ///< Uso de tuplas como valor de retorno
#include <cmath>             ///< Funciones matemáticas

using namespace std;

/**
 * @brief Calcula el azimut, la elevación y los vectores base (topocéntricos) locales.
 *
 * @param s Vector de posición en coordenadas cartesianas (3x1) [km].
 * @return Tupla con:
 *         - Azimut [rad]
 *         - Elevación [rad]
 *         - Vector sur (3x1)
 *         - Vector este (3x1)
 */
tuple<double, double, Matrix&, Matrix&> AzElPa(Matrix &s);

#endif  // _AZ_EL_PA_
