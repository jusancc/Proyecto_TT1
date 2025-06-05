//$Header$
//------------------------------------------------------------------------------
//                          AccelPointMass
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file AccelPointMass.cpp
 * @brief Implementación de la función que calcula la aceleración producida
 *        por un cuerpo puntual sobre un satélite.
 */
//------------------------------------------------------------------------------

#include "../include/AccelPointMass.hpp"

//------------------------------------------------------------------------------
//  Matrix& AccelPointMass(Matrix &r, Matrix &s, double GM)
//------------------------------------------------------------------------------
/**
 * @brief Calcula la aceleración gravitatoria ejercida por un cuerpo puntual.
 *
 * Esta función evalúa la perturbación sobre un satélite debida a la atracción
 * gravitacional de un cuerpo puntual (como el Sol, la Luna o planetas). Se utiliza
 * el modelo de atracción diferencial, considerando la aceleración relativa entre
 * el satélite y el centro de masas del sistema.
 * 
 * @param r Vector (1x3) de posición del satélite en coordenadas cartesianas [km].
 * @param s Vector (1x3) de posición del cuerpo perturbador en el mismo sistema [km].
 * @param GM Producto de la constante gravitacional por la masa del cuerpo [km³/s²].
 * @return Referencia a una matriz (1x3) con la aceleración resultante [km/s²].
 */
//------------------------------------------------------------------------------
Matrix& AccelPointMass(Matrix &r, Matrix &s, double GM){
    // Validación de dimensiones: ambos vectores deben ser 1x3
    if (r.n_row != 1 || r.n_column != 3 || s.n_row != 1 || s.n_column != 3) {
        exit(EXIT_FAILURE);  // Salir si dimensiones incorrectas
    }

    // Diferencia de posición (vector relativo)
    Matrix& d = r - s;

    double norm_d = d.norm(); // Magnitud del vector relativo
    double norm_s = s.norm(); // Magnitud de la posición del cuerpo

    // Componentes del modelo de atracción diferencial
    Matrix& term1 = d / pow(norm_d, 3);  // componente satélite-cuerpo
    Matrix& term2 = s / pow(norm_s, 3);  // componente cuerpo-centro de masas

    Matrix& acc = term1 + term2;

    Matrix& result = acc * (-GM);  // aceleración total

    return result;
}
