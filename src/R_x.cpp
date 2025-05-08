#include <iostream>
#include <iomanip>
#include "../include/Matrix.hpp"
#include "../include/R_x.hpp"

//-------------------------------------------------------------
//  input:
//    angle       - angle of rotation [rad]
//
//  output:
//    rotmat      - vector result
//-------------------------------------------------------------


Matrix& R_x(double alpha)
{
    Matrix &rotmat = zeros(3,3);
    double C, S;
    
    C = cos(alpha);
    S = sin(alpha);
    // rotmat = zeros(3,3);
    
    rotmat(1,1) = 1.0;  rotmat(1,2) =    0.0;  rotmat(1,3) = 0.0;
    rotmat(2,1) = 0.0;  rotmat(2,2) =      C;  rotmat(2,3) =   S;
    rotmat(3,1) = 0.0;  rotmat(3,2) = -1.0*S;  rotmat(3,3) =   C;
    
    return rotmat;
}
