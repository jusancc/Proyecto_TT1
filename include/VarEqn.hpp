//$Header$
//------------------------------------------------------------------------------
//                                  VarEqn.hpp
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file VarEqn.hpp
 * @brief Prototipo de la función de ecuaciones variacionales para propagación orbital.
 *
 * Esta función calcula la derivada del vector de estado extendido que incluye la posición,
 * velocidad y matriz de variación del estado Phi (6x6). Se usa en el contexto de
 * propagación numérica con perturbaciones gravitatorias armónicas, típicamente para
 * aplicaciones de dinámica orbital y estimación de estado.
 */
//------------------------------------------------------------------------------

#ifndef _VAR_EQN_
#define _VAR_EQN_

#include "matrix.hpp"
#include "IERS.hpp"
#include "global.hpp"
#include "SAT_Const.hpp"
#include "Timediff.hpp"
#include "PrecMatrix.hpp"
#include "NutMatrix.hpp"
#include "PoleMatrix.hpp"
#include "GHAMatrix.hpp"
#include "AccelHarmonic.hpp"
#include "G_AccelHarmonic.hpp"

/**
 * @brief Calcula la derivada del vector extendido y matriz de variación (ecuaciones variacionales).
 *
 * @param x Tiempo (en segundos desde la época inicial Mjd_0).
 * @param yPhi Vector columna de dimensión 42x1, compuesto por:
 *             - Posición (elementos 1-3)
 *             - Velocidad (elementos 4-6)
 *             - Matriz de variación del estado Phi (elementos 7-42, almacenada por columnas).
 * @return Referencia al vector derivada de yPhi con respecto al tiempo.
 */
Matrix& varEqn(double x, Matrix &yPhi);

#endif
