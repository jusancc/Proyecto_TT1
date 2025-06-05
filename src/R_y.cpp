//$Header$
//------------------------------------------------------------------------------
//                                   R_y
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file R_y.cpp
 * @brief Implementación de la función R_y para generar una matriz de rotación sobre el eje Y.
 *
 * Esta función construye una matriz 3x3 de rotación alrededor del eje Y para un
 * ángulo especificado en radianes. Se utiliza comúnmente en transformaciones
 * espaciales dentro de la mecánica orbital y sistemas de referencia tridimensionales.
 */
//------------------------------------------------------------------------------

#include "../include/R_y.hpp"

/**
 * @brief Genera la matriz de rotación sobre el eje Y.
 * 
 * @param angle Ángulo de rotación en radianes.
 * @return Referencia a una matriz 3x3 correspondiente a la rotación.
 */
Matrix& R_y(double angle){
    Matrix &rotmat = zeros(3,3);
    double C = cos(angle);
    double S = sin(angle);

    rotmat(1,1) =  C;  rotmat(1,2) = 0.0;  rotmat(1,3) = -S;
    rotmat(2,1) = 0.0; rotmat(2,2) = 1.0;  rotmat(2,3) = 0.0;
    rotmat(3,1) =  S;  rotmat(3,2) = 0.0;  rotmat(3,3) =  C;

    return rotmat;
}
