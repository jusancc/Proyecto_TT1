#include "../include/LTC.hpp"

Matrix& LTC(double lon, double lat){
    Matrix Ry = R_y(-lat);
    Matrix Rz = R_z(lon);
    Matrix &M = Ry * Rz;

    Matrix& row1 = M.extract_row(1);
    Matrix& row2 = M.extract_row(2);
    Matrix& row3 = M.extract_row(3);

    M.assign_row(1, row2);
    M.assign_row(2, row3);
    M.assign_row(3, row1);
    
    return M;
}