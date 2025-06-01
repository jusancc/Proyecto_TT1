#include "../include/G_AccelHarmonic.hpp"
#include <iostream>

Matrix& G_AccelHarmonic(Matrix &r, Matrix &U, double n_max, double m_max){
    double d = 1.0;

    Matrix &G = zeros(3,3);
    Matrix &dr = zeros(3,1);

    for (int i=1; i<=3; i++){
        dr = zeros(3,1);
        dr(i) = d;

        Matrix &a_plus = AccelHarmonic(r + dr/2, U, n_max, m_max);

        Matrix &a_minus = AccelHarmonic(r - dr/2, U, n_max, m_max);

        Matrix &da = a_plus - a_minus;

        G.assign_column(i, da/d);

    }

    return G;
}
