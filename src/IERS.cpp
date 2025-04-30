#include "../include/IERS.hpp"

tuple<double, double, double, double, double, double, double, double, double> IERS(Matrix eop, double Mjd_UTC, char interp)
{

    double x_pole;
    double y_pole;
    double UT1_UTC;
    double LOD;
    double dpsi;
    double deps;
    double dx_pole;
    double dy_pole;
    double TAI_UTC;

    if (interp == '1')
    {
        double mjd = floor(Mjd_UTC);
        Matrix f4 = eop.extract_row(4);
        int i;
        for (i = 0; i < 4; i++)
        {
            if (f4(i, 1) == mjd)
            {
                i = f4(i, 1);
                break;
            }
        }

        Matrix preeop = eop.extract_column(i);
        Matrix nexteop = eop.extract_column(i + 1);

        double mfme = 1440 * (Mjd_UTC - floor(Mjd_UTC));
        double fixf = mfme / 1440;

        x_pole = preeop(1, 5) + (nexteop(1, 5) - preeop(1, 5)) * fixf;
        y_pole = preeop(1, 6) + (nexteop(1, 6) - preeop(1, 6)) * fixf;
        UT1_UTC = preeop(1, 7) + (nexteop(1, 7) - preeop(1, 7)) * fixf;
        LOD = preeop(1, 8) + (nexteop(1, 8) - preeop(1, 8)) * fixf;
        dpsi = preeop(1, 9) + (nexteop(1, 9) - preeop(1, 9)) * fixf;
        deps = preeop(1, 10) + (nexteop(1, 10) - preeop(1, 10)) * fixf;
        dx_pole = preeop(1, 11) + (nexteop(1, 11) - preeop(1, 11)) * fixf;
        dy_pole = preeop(1, 12) + (nexteop(1, 12) - preeop(1, 12)) * fixf;
        TAI_UTC = preeop(1, 13);

        x_pole = x_pole / SAT_Const::Arcs;
        y_pole = y_pole / SAT_Const::Arcs;
        dpsi = dpsi / SAT_Const::Arcs;
        deps = deps / SAT_Const::Arcs;
        dx_pole = SAT_Const::Arcs;
        dy_pole = dy_pole / SAT_Const::Arcs;
    }
    else if (interp == 'n')
    {
        double mjd = floor(Mjd_UTC);
        Matrix f4 = eop.extract_row(4);
        int i;
        for (i = 0; i < 4; i++)
        {
            if (f4(i, 1) == mjd)
            {
                i = f4(i, 1);
                break;
            }
        }
        eop = eop.extract_column(i);
        x_pole = eop(1, 5) / SAT_Const::Arcs;
        y_pole = eop(1, 6) / SAT_Const::Arcs;
        UT1_UTC = eop(1, 7);
        LOD = eop(1, 8);
        dpsi = eop(1, 9) / SAT_Const::Arcs;
        deps = eop(1, 10) / SAT_Const::Arcs;
        dx_pole = eop(1, 11) / SAT_Const::Arcs;
        dy_pole = eop(1, 12) / SAT_Const::Arcs;
        TAI_UTC = eop(1, 13);
    }

    return tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);
}