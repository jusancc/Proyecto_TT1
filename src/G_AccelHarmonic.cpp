#include "../include/G_AccelHarmonic.hpp"

Matrix& G_AccelHarmonic(Matrix &r, Matrix &U, double n_max, double m_max){
    double d = 1.0;

    Matrix &G = zeros(3);
    Matrix &dr = zeros(3,1);

    for (int i=1;i<=3;i++){
        dr = zeros(3);
        dr(i,1) = d;

        Matrix &da = AccelHarmonic(r+dr/2,U,n_max,m_max)-AccelHarmonic(r-dr/2,U,n_max,m_max);

        G.assign_column(i,da/d);
    }

    return G;
}