//$Header$
//------------------------------------------------------------------------------
//                              GHAMatrix
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file GHAMatrix.cpp
 * @brief Implementación de la función que genera la matriz de rotación
 *        correspondiente al Ángulo Horario de Greenwich (GHA).
 */
//------------------------------------------------------------------------------

#include "../include/GHAMatrix.hpp"

//------------------------------------------------------------------------------
//  Matrix& GHAMatrix(double Mjd_UT1)
//------------------------------------------------------------------------------
/**
 * @brief Devuelve la matriz de rotación que aplica el ángulo horario de Greenwich (GHA).
 *
 * Calcula el GAST (Greenwich Apparent Sidereal Time) a partir de la fecha UT1
 * y retorna la matriz de rotación correspondiente mediante una rotación alrededor
 * del eje Z. Esta transformación es esencial para convertir coordenadas del sistema
 * inercial al sistema de la Tierra rotante.
 *
 * @param Mjd_UT1 Fecha en tiempo universal UT1 (Modified Julian Date).
 * @return Referencia a la matriz de rotación (3x3) que representa el GHA.
 */
Matrix& GHAMatrix(double Mjd_UT1){
    return R_z(gast(Mjd_UT1));
}
