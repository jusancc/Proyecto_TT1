#include "../include/AccelPointMass.hpp"

Matrix& AccelPointMass(Matrix &r, Matrix &s, double GM){
    if (r.n_row != 1 || r.n_column != 3 || s.n_row != 1 || s.n_column != 3) {
        std::cout << "AccelPointMass: Los vectores r y s deben ser de tamaño 1,3\n";
        exit(EXIT_FAILURE);
    }

    Matrix& d = r - s;
    double norm_d = d.norm();
    double norm_s = s.norm();

    Matrix& term1 = d / pow(norm_d, 3);
    Matrix& term2 = s / pow(norm_s, 3);
    Matrix& acc = term1 + term2;

    Matrix& result = acc * (-GM);

    return result;
}