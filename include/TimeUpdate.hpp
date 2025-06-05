//$Header$
//------------------------------------------------------------------------------
//                                TimeUpdate
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file TimeUpdate.hpp
 * @brief Declaración de la función TimeUpdate para el paso de tiempo del filtro de Kalman.
 *
 * Esta función realiza la actualización del error de covarianza en el modelo de predicción
 * del filtro de Kalman extendido.
 */
//------------------------------------------------------------------------------

#ifndef _TIME_UPDATE_
#define _TIME_UPDATE_

#include "matrix.hpp"

/**
 * @brief Actualiza la matriz de covarianza P mediante el modelo de transición de estados Phi.
 *
 * @param P Matriz de covarianza (n x n)
 * @param Phi Matriz de transición del modelo (n x n)
 * @param Qdt Varianza del proceso multiplicada por el paso temporal (escalar)
 * @return Referencia a la nueva matriz de covarianza P actualizada
 */
Matrix& TimeUpdate(Matrix &P, Matrix &Phi, double Qdt);

#endif
