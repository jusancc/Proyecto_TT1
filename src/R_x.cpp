//$Header$
//------------------------------------------------------------------------------
//                                   R_x
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file R_x.cpp
 * @brief Implementación de la matriz de rotación respecto al eje X.
 *
 * Esta función genera una matriz de rotación 3x3 para un ángulo dado,
 * que rota un vector en sentido antihorario alrededor del eje X.
 * Se utiliza comúnmente en transformaciones espaciales tridimensionales.
 */
//------------------------------------------------------------------------------

#include "../include/R_x.hpp"

/**
 * @brief Genera la matriz de rotación respecto al eje X.
 * 
 * @param alpha Ángulo de rotación en radianes.
 * @return Referencia a una matriz 3x3 con la rotación aplicada.
 * 
 * La forma de la matriz de rotación es:
 * 
 *     [ 1     0        0     ]
 *     [ 0   cos(α)   sin(α) ]
 *     [ 0  -sin(α)   cos(α) ]
 */

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
