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
 * @file TimeUpdate.cpp
 * @brief Implementación de la función TimeUpdate para el filtro de Kalman.
 *
 * Esta función actualiza la matriz de covarianza del error P en el paso de predicción
 * del filtro de Kalman extendido, aplicando la matriz de transición Phi y añadiendo
 * el término de ruido de proceso Qdt.
 */
//------------------------------------------------------------------------------

#include "../include/TimeUpdate.hpp"

/**
 * @brief Realiza la actualización de la matriz de covarianza P.
 *
 * Aplica la fórmula P = Φ·P·Φᵗ + Qdt, donde Qdt es el ruido del proceso.
 *
 * @param P Matriz de covarianza del estado.
 * @param Phi Matriz de transición de estados.
 * @param Qdt Valor escalar que representa el efecto del ruido del proceso.
 * @return Referencia a la matriz de covarianza P actualizada.
 */
Matrix& TimeUpdate(Matrix &P, Matrix &Phi, double Qdt) {
    Matrix PTrans = transpose(Phi);
    P = Phi * P * PTrans + Qdt;
    return P;
}
