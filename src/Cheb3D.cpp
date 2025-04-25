#include "../include/Cheb3D.hpp"

Matrix Cheb3D(double t, int N, double Ta, double Tb, double Cx[], double Cy[], double Cz[]){
    if (t<Ta || Tb<t)
    {
        std:cout << "Time out of range in Cheb3D";
        exit(EXIT_FAILURE);
    }

    double tau;
    Matrix f1(1,3);
    Matrix f2(1,3);

    tau = (2*t-Ta-Tb)/(Tb-Ta);

    Matrix old_f1(1,3);

    for (int i = N - 1; i >= 1; --i) {
        old_f1 = f1;

        f1 = f1 * (2.0 * tau) - f2;

        f1(1, 1) += Cx[i];
        f1(1, 2) += Cy[i];
        f1(1, 3) += Cz[i];

        f2 = old_f1;
    }

    Matrix ChebApp(1,3);

    ChebApp = f1 * tau - f2;
    ChebApp(1,1) += Cx[0];
    ChebApp(1,2) += Cy[0];
    ChebApp(1,3) += Cx[0];

    return ChebApp;
    
}