//$Header$
//------------------------------------------------------------------------------
//                                  LTC
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file LTC.cpp
 * @brief Implementación de la matriz de transformación LTC (Local Tangent Coordinate).
 *
 * Esta función construye la matriz que transforma coordenadas cartesianas ECEF
 * al sistema de coordenadas locales topocéntricas tangentes a la superficie terrestre
 * en una ubicación específica dada por latitud y longitud.
 */
//------------------------------------------------------------------------------

#include "../include/LTC.hpp"

/**
 * @brief Genera la matriz de transformación del sistema ECEF al sistema local topocéntrico.
 *
 * @param lon Longitud geodésica [rad]
 * @param lat Latitud geodésica [rad]
 * @return Referencia a una matriz 3x3 con la transformación LTC.
 */
Matrix& LTC(double lon, double lat){
    Matrix Ry = R_y(-lat);  ///< Rotación alrededor del eje Y negativo (ajuste de latitud)
    Matrix Rz = R_z(lon);   ///< Rotación alrededor del eje Z (ajuste de longitud)
    Matrix &M = Ry * Rz;    ///< Producto de rotaciones para formar la matriz de transformación

    // Reordenamiento de filas para formar el sistema topocéntrico: Este, Norte, Cenit
    Matrix& row1 = M.extract_row(1);
    Matrix& row2 = M.extract_row(2);
    Matrix& row3 = M.extract_row(3);

    M.assign_row(1, row2);  ///< Nueva fila 1 = fila 2 (Este)
    M.assign_row(2, row3);  ///< Nueva fila 2 = fila 3 (Norte)
    M.assign_row(3, row1);  ///< Nueva fila 3 = fila 1 (Cenit)

    return M;
}
