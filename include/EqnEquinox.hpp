//$Header$
//------------------------------------------------------------------------------
//                              EqnEquinox
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file EqnEquinox.hpp
 * @brief Declaración de la función que calcula la ecuación de los equinoccios.
 *
 * La ecuación de los equinoccios es la diferencia entre el ángulo de Greenwich
 * aparente y medio, y se utiliza en cálculos de tiempo sidéreo.
 */
//------------------------------------------------------------------------------

#ifndef _EQN_EQUINOX_
#define _EQN_EQUINOX_

#include "NutAngles.hpp"       ///< Ángulos de nutación Δψ y Δε
#include "MeanObliquity.hpp"   ///< Oblicuidad media de la eclíptica
#include <cmath>               ///< Funciones matemáticas estándar
#include <tuple>               ///< Uso de tuplas para resultados compuestos
#include <iostream>            ///< Entrada/salida estándar

using namespace std;

/**
 * @brief Calcula la ecuación de los equinoccios para una fecha dada.
 *
 * @param Mjd_TT Fecha en tiempo terrestre (Modified Julian Date TT).
 * @return Ecuación de los equinoccios (rad).
 */
double EqnEquinox(double Mjd_TT);

#endif  // _EQN_EQUINOX_
