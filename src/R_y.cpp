#include "../include/R_y.hpp"

Matrix& R_y(double angle){
    Matrix &rotmat = zeros(3,3);
    double C,S;

    C = cos(angle);
    S = sin(angle);

    rotmat(1,1) = C;  rotmat(1,2) =    0.0;  rotmat(1,3) = -1.0 * S;
    rotmat(2,1) = 0.0;  rotmat(2,2) =  1.0;  rotmat(2,3) =   0.0;
    rotmat(3,1) = S;  rotmat(3,2) = 0.0;  rotmat(3,3) =   C;

    return rotmat;
}