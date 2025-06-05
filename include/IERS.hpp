//$Header$
//------------------------------------------------------------------------------
//                                  IERS
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file IERS.hpp
 * @brief Declaración de la función que interpola parámetros EOP del IERS.
 *
 * Esta función proporciona valores interpolados o discretos de parámetros de orientación
 * de la Tierra (EOP) necesarios para modelos de tiempo y transformaciones de marcos.
 */
//------------------------------------------------------------------------------

#ifndef _IERS_
#define _IERS_

#include "SAT_Const.hpp"   ///< Constantes astronómicas
#include "matrix.hpp"      ///< Clase Matrix personalizada
#include "iostream"        ///< Entrada/salida estándar
#include <tuple>           ///< Para retornar múltiples valores

using namespace std;

/**
 * @brief Recupera o interpola los parámetros IERS para una fecha UTC dada.
 *
 * Los parámetros incluyen desplazamientos del polo, diferencias de tiempo
 * y correcciones de nutación para una fecha en tiempo UTC.
 *
 * @param eop Matriz con los datos EOP cargados previamente.
 * @param Mjd_UTC Fecha UTC en formato Modified Julian Date.
 * @param interp Método de interpolación ('l' = lineal, 'n' = ninguno).
 * @return Tupla con:
 *         - x_pole [rad]
 *         - y_pole [rad]
 *         - UT1 - UTC [s]
 *         - LOD (Length of Day) [s]
 *         - dpsi (nutación en longitud) [rad]
 *         - deps (nutación en oblicuidad) [rad]
 *         - dx_pole [rad]
 *         - dy_pole [rad]
 *         - TAI - UTC [s]
 */
tuple<double, double, double, double, double, double, double, double, double>
IERS(Matrix &eop, double Mjd_UTC, char interp = 'n');

#endif  // _IERS_
