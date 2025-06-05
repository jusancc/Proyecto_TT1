//$Header$
//------------------------------------------------------------------------------
//                               PoleMatrix
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file PoleMatrix.cpp
 * @brief Implementación de la función que calcula la matriz de transformación por efecto del movimiento polar.
 *
 * Esta matriz es usada para corregir la orientación del sistema de referencia terrestre
 * debido al desplazamiento de los polos geográficos respecto al eje de rotación.
 */
//------------------------------------------------------------------------------

#include "../include/PoleMatrix.hpp"

/**
 * @brief Calcula la matriz de transformación por efecto del movimiento polar.
 *
 * El desplazamiento de los polos (xp, yp) afecta a la orientación del sistema de coordenadas 
 * terrestres. Esta función genera una matriz de rotación que corrige este efecto aplicando
 * dos rotaciones consecutivas: una en torno al eje Y y otra en torno al eje X.
 *
 * @param xp Desplazamiento del polo en dirección X [rad].
 * @param yp Desplazamiento del polo en dirección Y [rad].
 * @return Referencia a una matriz 3x3 resultante de la composición de ambas rotaciones.
 */
Matrix& PoleMatrix(double xp, double yp){
    Matrix &ry = R_y(-xp);   // Rotación alrededor del eje Y
    Matrix &rx = R_x(-yp);   // Rotación alrededor del eje X
    Matrix &PoleMat = ry * rx;
    return PoleMat;
}
