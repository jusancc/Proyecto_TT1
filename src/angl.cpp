//$Header$
//------------------------------------------------------------------------------
//                                  angl
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file angl.cpp
 * @brief Implementación de la función que calcula el ángulo entre dos vectores.
 */
//------------------------------------------------------------------------------

#include "../include/angl.hpp"

//------------------------------------------------------------------------------
//  double angl(Matrix& vec1, Matrix& vec2)
//------------------------------------------------------------------------------
/**
 * @brief Calcula el ángulo entre dos vectores tridimensionales.
 *
 * Utiliza el producto escalar para determinar el ángulo entre dos vectores dados.
 * En caso de que uno de los vectores tenga magnitud muy pequeña, se devuelve un
 * valor indefinido (999999.1) como marcador de error.
 * 
 * @param vec1 Primer vector (3x1 o 1x3).
 * @param vec2 Segundo vector (3x1 o 1x3).
 * @return Ángulo entre `vec1` y `vec2` en radianes. Si no se puede calcular, retorna 999999.1.
 */
//------------------------------------------------------------------------------
double angl(Matrix& vec1, Matrix& vec2){
    double small = 0.00000001;      // Umbral para considerar vectores nulos
    double undefined = 999999.1;    // Valor de ángulo indefinido

    double magv1 = vec1.norm();     // Magnitud del primer vector
    double magv2 = vec2.norm();     // Magnitud del segundo vector

    double theta;
    if (magv1 * magv2 > pow(small, 2)) {
        double temp = vec1.dot(vec2) / (magv1 * magv2);

        // Corrección numérica para asegurar dominio de acos
        if (fabs(temp) > 1.0) {
            if (temp > 0)
                temp = 1.0;
            else if (temp < 0)
                temp = -1.0;
            else
                temp = 0.0;
        }

        theta = acos(temp);  // Ángulo entre los vectores
    } else {
        theta = undefined;   // Ángulo no definido si alguno es nulo
    }

    return theta;
}
