#include "../include/IERS.hpp"

tuple<double, double, double, double, double, double, double, double, double> IERS(Matrix &eop, double Mjd_UTC, char interp) {

    double x_pole;
    double y_pole;
    double UT1_UTC;
    double LOD;
    double dpsi;
    double deps;
    double dx_pole;
    double dy_pole;
    double TAI_UTC;

    if (interp == 'l') {  // linear interpolation
        double mjd = floor(Mjd_UTC);
        int i = -1;
        for (int col = 1; col <= eop.n_column; col++) {
            if (mjd == eop(4, col)) {  // Buscar el MJD en la columna
                i = col;
                break;
            }
        }

        if (i == -1) {
            cout << "IERS: MJD not found in eop data.\n";
            exit(EXIT_FAILURE);  // Salir si no se encuentra el MJD
        }

        Matrix& preeop = eop.extract_column(i);
        Matrix& nexteop = eop.extract_column(i + 1);

        double mfme = 1440 * (Mjd_UTC - floor(Mjd_UTC));
        double fixf = mfme / 1440;

        // Interpolación lineal de los parámetros de rotación de la Tierra
        x_pole = preeop(5) + (nexteop(5) - preeop(5)) * fixf;
        y_pole = preeop(6) + (nexteop(6) - preeop(6)) * fixf;
        UT1_UTC = preeop(7) + (nexteop(7) - preeop(7)) * fixf;
        LOD = preeop(8) + (nexteop(8) - preeop(8)) * fixf;
        dpsi = preeop(9) + (nexteop(9) - preeop(9)) * fixf;
        deps = preeop(10) + (nexteop(10) - preeop(10)) * fixf;
        dx_pole = preeop(11) + (nexteop(11) - preeop(11)) * fixf;
        dy_pole = preeop(12) + (nexteop(12) - preeop(12)) * fixf;
        TAI_UTC = preeop(13);

        // Convertir de grados a radianes usando la constante Arcs
        x_pole = x_pole / SAT_Const::Arcs;  // Pole coordinate [rad]
        y_pole = y_pole / SAT_Const::Arcs;  // Pole coordinate [rad]
        dpsi = dpsi / SAT_Const::Arcs;
        deps = deps / SAT_Const::Arcs;
        dx_pole = dx_pole / SAT_Const::Arcs;  // Pole coordinate [rad]
        dy_pole = dy_pole / SAT_Const::Arcs;  // Pole coordinate [rad]
    } 
    else if (interp == 'n') {  // No interpolation, just extract the data for the given MJD
        double mjd = floor(Mjd_UTC);
        int i = -1;
        for (int col = 1; col <= eop.n_column; col++) {
            if (mjd == eop(4, col)) {  // Buscar el MJD en la columna
                i = col;
                break;
            }
        }

        if (i == -1) {
            cout << "IERS: MJD not found in eop data.\n";
            exit(EXIT_FAILURE);  // Salir si no se encuentra el MJD
        }

        eop = eop.extract_column(i);

        // Extraer directamente los valores para el MJD encontrado
        x_pole = eop(5) / SAT_Const::Arcs;  // Pole coordinate [rad]
        y_pole = eop(6) / SAT_Const::Arcs;  // Pole coordinate [rad]
        UT1_UTC = eop(7);  // UT1-UTC time difference [s]
        LOD = eop(8);  // Length of day [s]
        dpsi = eop(9) / SAT_Const::Arcs;
        deps = eop(10) / SAT_Const::Arcs;
        dx_pole = eop(11) / SAT_Const::Arcs;  // Pole coordinate [rad]
        dy_pole = eop(12) / SAT_Const::Arcs;  // Pole coordinate [rad]
        TAI_UTC = eop(13);  // TAI-UTC time difference [s]
    }

    return tie(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
}
