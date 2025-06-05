//$Header$
//------------------------------------------------------------------------------
//                                  IERS
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file IERS.cpp
 * @brief Implementación de la función que devuelve parámetros interpolados
 *        de orientación terrestre a partir de los datos EOP.
 */
//------------------------------------------------------------------------------

#include "../include/IERS.hpp"

//------------------------------------------------------------------------------
//  tuple<double, double, double, double, double, double, double, double, double>
//  IERS(Matrix &eop, double Mjd_UTC, char interp)
//------------------------------------------------------------------------------
/**
 * @brief Recupera los parámetros EOP del IERS para una fecha UTC dada.
 *
 * Dependiendo del parámetro `interp`, esta función realiza interpolación lineal
 * o devuelve los valores directos para el día UTC dado. La información extraída
 * incluye desplazamientos polares, UT1-UTC, LOD, correcciones de nutación y
 * diferencia TAI-UTC.
 *
 * @param eop Matriz de datos EOP ya cargada (13 x N columnas).
 * @param Mjd_UTC Fecha en tiempo UTC (Modified Julian Date).
 * @param interp Modo de interpolación: 'l' = lineal, 'n' = sin interpolación.
 * @return Tupla con:
 *         - x_pole [rad]
 *         - y_pole [rad]
 *         - UT1 - UTC [s]
 *         - LOD [s]
 *         - Δψ [rad] (nutación en longitud)
 *         - Δε [rad] (nutación en oblicuidad)
 *         - dx_pole [rad]
 *         - dy_pole [rad]
 *         - TAI - UTC [s]
 *
 * @exception Termina el programa si la fecha `Mjd_UTC` no se encuentra en los datos.
 */
//------------------------------------------------------------------------------
tuple<double, double, double, double, double, double, double, double, double>
IERS(Matrix &eop, double Mjd_UTC, char interp) {
    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;

    if (interp == 'l') {
        // Interpolación lineal
        double mjd = floor(Mjd_UTC);
        int i = -1;
        for (int col = 1; col <= eop.n_column; col++) {
            if (mjd == eop(4, col)) {
                i = col;
                break;
            }
        }

        if (i == -1) {
            cout << "IERS: MJD not found in eop data.\n";
            exit(EXIT_FAILURE);
        }

        Matrix& preeop = eop.extract_column(i);
        Matrix& nexteop = eop.extract_column(i + 1);

        double mfme = 1440 * (Mjd_UTC - mjd);
        double fixf = mfme / 1440;

        x_pole   = preeop(5)  + (nexteop(5)  - preeop(5))  * fixf;
        y_pole   = preeop(6)  + (nexteop(6)  - preeop(6))  * fixf;
        UT1_UTC  = preeop(7)  + (nexteop(7)  - preeop(7))  * fixf;
        LOD      = preeop(8)  + (nexteop(8)  - preeop(8))  * fixf;
        dpsi     = preeop(9)  + (nexteop(9)  - preeop(9))  * fixf;
        deps     = preeop(10) + (nexteop(10) - preeop(10)) * fixf;
        dx_pole  = preeop(11) + (nexteop(11) - preeop(11)) * fixf;
        dy_pole  = preeop(12) + (nexteop(12) - preeop(12)) * fixf;
        TAI_UTC  = preeop(13);

        // Conversión de arcosegundos a radianes
        x_pole  /= SAT_Const::Arcs;
        y_pole  /= SAT_Const::Arcs;
        dpsi    /= SAT_Const::Arcs;
        deps    /= SAT_Const::Arcs;
        dx_pole /= SAT_Const::Arcs;
        dy_pole /= SAT_Const::Arcs;
    } 
    else if (interp == 'n') {
        // Sin interpolación, usar valores directos
        double mjd = floor(Mjd_UTC);
        int i = -1;
        for (int col = 1; col <= eop.n_column; col++) {
            if (mjd == eop(4, col)) {
                i = col;
                break;
            }
        }

        if (i == -1) {
            cout << "IERS: MJD not found in eop data.\n";
            exit(EXIT_FAILURE);
        }

        eop = eop.extract_column(i);

        x_pole   = eop(5)  / SAT_Const::Arcs;
        y_pole   = eop(6)  / SAT_Const::Arcs;
        UT1_UTC  = eop(7);
        LOD      = eop(8);
        dpsi     = eop(9)  / SAT_Const::Arcs;
        deps     = eop(10) / SAT_Const::Arcs;
        dx_pole  = eop(11) / SAT_Const::Arcs;
        dy_pole  = eop(12) / SAT_Const::Arcs;
        TAI_UTC  = eop(13);
    }

    return tie(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
}
