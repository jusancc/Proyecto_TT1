#include "../include/AccelHarmonic.hpp"

#include "../include/AccelHarmonic.hpp"
#include <iostream>
using namespace std;

Matrix& AccelHarmonic(Matrix &r, Matrix &E, int n_max, int m_max) {
    double r_ref = 6378.1363e3;
    double gm = 398600.4415e9;

    Matrix &r_bf = E * r;

    double d = transpose(r_bf).norm();

    double latgc = asin(r_bf(3, 1) / d);
    double lon = atan2(r_bf(2, 1), r_bf(1, 1));
    auto [pnm, dpnm] = Legendre(n_max, m_max, latgc);

    double dUdr = 0.0;
    double dUdlatgc = 0.0;
    double dUdlon = 0.0;

    for (int n = 0; n <= n_max; n++) {
        double b1 = (-gm / (d * d)) * pow(r_ref / d, n) * (n + 1);
        double b2 = (gm / d) * pow(r_ref / d, n);
        double b3 = (gm / d) * pow(r_ref / d, n);

        double q1 = 0.0;
        double q2 = 0.0;
        double q3 = 0.0;

        for (int m = 0; m <= m_max; m++) {
            double cos_ml = cos(m * lon);
            double sin_ml = sin(m * lon);

            double term1 = pnm(n + 1, m + 1) * (Cnm(n + 1, m + 1) * cos_ml + Snm(n + 1, m + 1) * sin_ml);
            double term2 = dpnm(n + 1, m + 1) * (Cnm(n + 1, m + 1) * cos_ml + Snm(n + 1, m + 1) * sin_ml);
            double term3 = m * pnm(n + 1, m + 1) * (Snm(n + 1, m + 1) * cos_ml - Cnm(n + 1, m + 1) * sin_ml);

            q1 += term1;
            q2 += term2;
            q3 += term3;
        }

        dUdr += q1 * b1;
        dUdlatgc += q2 * b2;
        dUdlon += q3 * b3;

    }

    double r2xy = pow(r_bf(1, 1), 2) + pow(r_bf(2, 1), 2);

    double ax = (1.0 / d * dUdr - r_bf(3, 1) / (d * d * sqrt(r2xy)) * dUdlatgc) * r_bf(1, 1) - (1.0 / r2xy * dUdlon) * r_bf(2, 1);
    double ay = (1.0 / d * dUdr - r_bf(3, 1) / (d * d * sqrt(r2xy)) * dUdlatgc) * r_bf(2, 1) + (1.0 / r2xy * dUdlon) * r_bf(1, 1);
    double az = 1.0 / d * dUdr * r_bf(3, 1) + sqrt(r2xy) / (d * d) * dUdlatgc;


    Matrix a_bf(3, 1);
    a_bf(1, 1) = ax;
    a_bf(2, 1) = ay;
    a_bf(3, 1) = az;

    Matrix &a = transpose(E) * a_bf;

    return a;
}
