//$Header$
//------------------------------------------------------------------------------
//                           JPL_Eph_DE430
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file JPL_Eph_DE430.hpp
 * @brief Declaración de la función que evalúa posiciones planetarias usando
 *        coeficientes Chebyshev de las efemérides JPL DE430.
 *
 * Esta función devuelve las posiciones de los principales cuerpos del sistema solar
 * interpoladas para un instante dado en tiempo TDB.
 */
//------------------------------------------------------------------------------

#ifndef _JPL_EPH_DE430_
#define _JPL_EPH_DE430_

#include <tuple>             ///< Para retornar múltiples vectores
#include "global.hpp"        ///< Parámetros globales y datos PC
#include "matrix.hpp"        ///< Clase Matrix para vectores de posición
#include "Cheb3D.hpp"        ///< Evaluación Chebyshev en 3D

/**
 * @brief Interpola las posiciones de cuerpos del sistema solar usando efemérides DE430.
 *
 * @param Mjd_TDB Tiempo dinámico baricéntrico (Modified Julian Date TDB).
 * @return Tupla con las posiciones (3x1) de:
 *         - Mercurio
 *         - Venus
 *         - Tierra
 *         - Marte
 *         - Júpiter
 *         - Saturno
 *         - Urano
 *         - Neptuno
 *         - Plutón
 *         - Luna
 *         - Sol
 */
tuple<Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&> JPL_Eph_DE430(double Mjd_TDB);

#endif  // _JPL_EPH_DE430_
