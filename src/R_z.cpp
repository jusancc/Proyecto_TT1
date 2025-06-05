//$Header$
//------------------------------------------------------------------------------
//                                   R_z
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file R_z.cpp
 * @brief Implementación de la función R_z para generar una matriz de rotación sobre el eje Z.
 *
 * Esta función genera una matriz de rotación de 3x3 que representa una rotación en torno
 * al eje Z por un ángulo especificado en radianes. Es útil en transformaciones
 * de coordenadas en el espacio tridimensional, especialmente en mecánica orbital.
 */
//------------------------------------------------------------------------------

#include "../include/R_z.hpp"

/**
 * @brief Genera una matriz de rotación en torno al eje Z.
 * 
 * @param angle Ángulo de rotación (en radianes).
 * @return Referencia a una matriz 3x3 que representa la rotación.
 */
Matrix& R_z(double angle){
    Matrix &rotmat = zeros(3,3);
    double C = cos(angle);
    double S = sin(angle);

    rotmat(1,1) =  C;   rotmat(1,2) =  S;   rotmat(1,3) = 0.0;
    rotmat(2,1) = -S;   rotmat(2,2) =  C;   rotmat(2,3) = 0.0;
    rotmat(3,1) = 0.0;  rotmat(3,2) = 0.0;  rotmat(3,3) = 1.0;

    return rotmat;
}
