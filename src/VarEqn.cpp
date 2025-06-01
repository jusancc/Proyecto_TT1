#include "../include/VarEqn.hpp"
#include <iostream>
#include <iomanip>
using namespace std;

Matrix& varEqn(double x, Matrix &yPhi) {
    cout << "[INFO] Entrando en varEqn..." << endl;

    // Tiempo y matrices
    auto [x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC] =
        IERS(eopdata, AuxParam.Mjd_UTC, 'l');
    auto [UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);
    double Mjd_UT1 = AuxParam.Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;

    Matrix &P = PrecMatrix(SAT_Const::MJD_J2000, AuxParam.Mjd_TT + x / 86400.0);
    Matrix &N = NutMatrix(AuxParam.Mjd_TT + x / 86400.0);
    Matrix &T = N * P;
    Matrix &E = (PoleMatrix(x_pole, y_pole)) * (GHAMatrix(Mjd_UT1)) * T;

    cout << "[INFO] Preparando r (3x1) y v (3x1)..." << endl;
    Matrix &r = zeros(3, 1);
    Matrix &v = zeros(3, 1);
    for (int i = 1; i <= 3; i++) {
        (r)(i, 1) = yPhi(i, 1);
        (v)(i, 1) = yPhi(i + 3, 1);
    }

    cout << "[DEBUG] r = \n" << r;
    cout << "[DEBUG] v = \n" << v;
    cout << "[DEBUG] E = \n" << E;
    cout << "[DEBUG] AuxParam.n = " << AuxParam.n << ", m = " << AuxParam.m << endl;
    cout << "[DEBUG] Cnm(3,2) = " << Cnm(3,2) << ", Snm(3,2) = " << Snm(3,2) << endl;

    // Matriz Phi
    Matrix &Phi = zeros(6, 6);
    for (int j = 1; j <= 6; j++) {
        Matrix &col = yPhi.extract_vector(6*j + 1, 6*j + 6); // fila 1x6
        Phi.assign_column(j, transpose(col));
    }

    cout << "[INFO] Calculando aceleración a..." << endl;
    Matrix &a = AccelHarmonic(r, E, AuxParam.n, AuxParam.m);

    cout << std::fixed << std::setprecision(12);
    cout << "[DEBUG] a = \n" << a;
    cout << "[TEST] a(1) = " << a(1,1) << " (esperado ≈ -5.1348367854085)" << endl;
    cout << "[TEST] a(2) = " << a(2,1) << " (esperado ≈ -2.97717622353621)" << endl;
    cout << "[TEST] a(3) = " << a(3,1) << " (esperado ≈ -3.70591776714204)" << endl;

    cout << "[INFO] Calculando gradiente G..." << endl;
    Matrix &G = G_AccelHarmonic(r, E, AuxParam.n, AuxParam.m);
    cout << "[DEBUG] G = \n" << G;

    // Jacobiano
    Matrix &dfdy = zeros(6, 6);
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            dfdy(i, j) = 0.0;
            dfdy(i+3, j) = G(i, j);
            dfdy(i, j+3) = (i == j) ? 1.0 : 0.0;
            dfdy(i+3, j+3) = 0.0;
        }
    }

    Matrix &Phip = dfdy * Phi;
    cout << "[DEBUG] dfdy = \n" << dfdy;
    cout << "[DEBUG] Phi = \n" << Phi;
    cout << "[DEBUG] Phip = \n" << Phip;

    // Derivada del vector completo
    Matrix &yPhip = zeros(42, 1);
    for (int i = 1; i <= 3; i++) {
        yPhip(i, 1) = (v)(i, 1);
        yPhip(i+3, 1) = a(i, 1);
    }

    for (int j = 1; j <= 6; j++) {
        for (int i = 1; i <= 6; i++) {
            yPhip(6*j + i, 1) = Phip(i, j);
        }
    }

    cout << "[DEBUG] yPhip (resultado final):" << endl;
    for (int i = 1; i <= 42; i++) {
        cout << "yPhip(" << i << ") = " << yPhip(i, 1) << endl;
    }

    cout << "[INFO] varEqn finalizado correctamente." << endl;
    return yPhip;
}
