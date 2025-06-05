//$Header$
//------------------------------------------------------------------------------
//                          AccelPointMass
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file AccelPointMass.hpp
 * @brief Declaración de la función que calcula la aceleración gravitacional
 *        debida a un cuerpo puntual como el Sol, la Luna o un planeta.
 *
 * Esta función es utilizada en dinámica orbital para sumar la perturbación
 * ejercida por un cuerpo externo sobre un satélite artificial.
 */
//------------------------------------------------------------------------------

#ifndef _ACCEL_POINT_MASS_
#define _ACCEL_POINT_MASS_

#include "matrix.hpp"    ///< Clase Matrix personalizada para álgebra lineal
#include <cmath>         ///< Funciones matemáticas estándar

using namespace std;

/**
 * @brief Calcula la aceleración gravitatoria causada por un cuerpo puntual.
 *
 * @param r Vector de posición (3x1) del satélite en el sistema inercial.
 * @param s Vector de posición (3x1) del cuerpo puntual en el mismo sistema.
 * @param GM Producto de la constante gravitacional por la masa del cuerpo [km³/s²].
 * @return Referencia a un vector (3x1) con la aceleración resultante [km/s²].
 */
Matrix& AccelPointMass(Matrix &r, Matrix &s, double GM);

#endif  // _ACCEL_POINT_MASS_
