//$Header$
//------------------------------------------------------------------------------
//                                PrecMatrix
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file PrecMatrix.cpp
 * @brief Implementación de la función que calcula la matriz de precesión.
 *
 * Esta función genera la matriz de precesión entre dos fechas dadas, usando
 * el modelo IAU 1976. La precesión es una rotación lenta del eje de la Tierra
 * que cambia las coordenadas astronómicas a lo largo del tiempo.
 */
//------------------------------------------------------------------------------

#include "../include/PrecMatrix.hpp"

/**
 * @brief Calcula la matriz de precesión entre dos fechas dadas.
 *
 * @param Mjd_1 Fecha inicial en Tiempo Juliano Modificado (MJD).
 * @param Mjd_2 Fecha final en Tiempo Juliano Modificado (MJD).
 * @return Matriz de rotación 3x3 que representa la transformación debida a la precesión.
 *
 * La matriz se calcula como Rz(-z) * Ry(theta) * Rz(-zeta), donde los ángulos se
 * obtienen de las expresiones del modelo IAU 1976.
 */
Matrix& PrecMatrix(double Mjd_1, double Mjd_2){
    // Tiempo en siglos julianos desde J2000
    double T = (Mjd_1 - SAT_Const::MJD_J2000) / 36525;
    double dT = (Mjd_2 - Mjd_1) / 36525;

    // Ángulos de precesión en radianes
    double zeta = ( (2306.2181 + (1.39656 - 0.000139 * T) * T)
                  + ((0.30188 - 0.000344 * T) + 0.017998 * dT) * dT ) * dT / SAT_Const::Arcs;

    double z = zeta + ( (0.79280 + 0.000411 * T) + 0.000205 * dT ) * dT * dT / SAT_Const::Arcs;

    double theta = ( (2004.3109 - (0.85330 + 0.000217 * T) * T)
                   - ((0.42665 + 0.000217 * T) + 0.041833 * dT) * dT ) * dT / SAT_Const::Arcs;

    // Construcción de la matriz de precesión: Rz(-z) * Ry(theta) * Rz(-zeta)
    Matrix &rz1 = R_z(-z);
    Matrix &ry = R_y(theta);
    Matrix &rz2 = R_z(-zeta);

    Matrix &PrecMat = rz1 * ry * rz2;

    return PrecMat;
}
