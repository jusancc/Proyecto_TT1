//$Header$
//------------------------------------------------------------------------------
//                                 Mjday_TDB
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file Mjday_TDB.cpp
 * @brief Implementación de la conversión de MJD_TT a MJD_TDB.
 *
 * Esta función calcula la conversión desde el Tiempo Terrestre (TT) al
 * Tiempo Dinámico Baricéntrico (TDB) utilizando una expansión trigonométrica
 * más precisa basada en términos del JPL.
 * 
 * La expresión empleada es:
 * \f[
 * \text{TDB} = \text{TT} + \frac{1}{86400} \cdot \sum_{i=1}^{n} A_i \cdot \sin(B_i \cdot T + C_i)
 * \f]
 * 
 * donde \( T = \frac{MJD_{TT} - 51544.5}{36525} \) es el tiempo juliano en siglos.
 */
//------------------------------------------------------------------------------

#include "../include/Mjday_TDB.hpp"
#include <cmath>

double Mjday_TDB(double Mjd_TT) {
    double T_TT = (Mjd_TT - 51544.5) / 36525.0;

    Mjd_TT += (
          0.001658 * sin(628.3076 * T_TT + 6.2401)
        + 0.000022 * sin(575.3385 * T_TT + 4.2970)
        + 0.000014 * sin(1256.6152 * T_TT + 6.1969)
        + 0.000005 * sin(606.9777 * T_TT + 4.0212)
        + 0.000005 * sin(52.9691 * T_TT + 0.4444)
        + 0.000002 * sin(21.3299 * T_TT + 5.5431)
        + 0.000010 * sin(628.3076 * T_TT + 4.2490)
    ) / 86400.0;

    return Mjd_TT;
}
