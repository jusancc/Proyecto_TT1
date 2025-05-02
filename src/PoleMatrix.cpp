#include "../include/PoleMatrix.hpp"

Matrix& PoleMatrix(double xp, double yp){
    Matrix ry = R_y(-xp);
    Matrix rx = R_x(-yp);
    Matrix PoleMat = ry*rx;
    return PoleMat;
}