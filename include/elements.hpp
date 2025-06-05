//$Header$
//------------------------------------------------------------------------------
//                                elements
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file elements.hpp
 * @brief Declaración de la función que calcula elementos orbitales keplerianos
 *        a partir del vector de estado (posición y velocidad).
 *
 * Esta cabecera define la interfaz de la función `elements`, que recibe un 
 * vector de estado y devuelve una tupla con los elementos keplerianos:
 * p, a, e, i, Omega, omega, M.
 */
//------------------------------------------------------------------------------

#ifndef _ELEMENTS_
#define _ELEMENTS_

#include "../include/matrix.hpp"       ///< Clase Matrix personalizada para álgebra lineal
#include "../include/SAT_Const.hpp"    ///< Constantes del sistema satelital
#include <cmath>                       ///< Funciones matemáticas estándar
#include <tuple>                       ///< Para devolver múltiples valores

using namespace std;

/**
 * @brief Calcula los elementos orbitales keplerianos a partir del vector de estado.
 * 
 * @param y Vector de estado (6x1): posición (x,y,z) y velocidad (vx,vy,vz).
 * @return tuple<double, double, double, double, double, double, double> con:
 *         p (semilatus rectum), a (semi-eje mayor), e (excentricidad),
 *         i (inclinación), Omega (longitud nodo ascendente), 
 *         omega (argumento del periastro), M (anomalía media).
 */
tuple<double, double, double, double, double, double, double> elements(Matrix& y);

#endif  // _ELEMENTS_
