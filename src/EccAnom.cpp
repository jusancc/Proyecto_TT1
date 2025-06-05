//$Header$
//------------------------------------------------------------------------------
//                                EccAnom
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file EccAnom.cpp
 * @brief Implementación de la función que resuelve la ecuación de Kepler
 *        para obtener la anomalía excéntrica a partir de la media.
 */
//------------------------------------------------------------------------------

#include "../include/EccAnom.hpp"

//------------------------------------------------------------------------------
//  double EccAnom(double M, double e)
//------------------------------------------------------------------------------
/**
 * @brief Resuelve la ecuación de Kepler: M = E - e·sin(E) usando Newton-Raphson.
 *
 * Esta función calcula la anomalía excéntrica `E` en radianes a partir de la
 * anomalía media `M` y la excentricidad `e`, mediante un proceso iterativo.
 * 
 * Si la excentricidad es menor que 0.8 se toma E₀ = M como aproximación inicial;
 * en caso contrario, se toma E₀ = π. La iteración se detiene al alcanzar
 * una tolerancia basada en la precisión de máquina.
 *
 * @param M Anomalía media (rad).
 * @param e Excentricidad (0 ≤ e < 1).
 * @return Anomalía excéntrica E (rad).
 *
 * @exception Termina con exit(EXIT_FAILURE) si no converge en 15 iteraciones.
 */
//------------------------------------------------------------------------------
double EccAnom(double M, double e){
    int maxit = 15;
    int i = 1;
    double E, f;
    const double epsilon = 1e2 * std::numeric_limits<double>::epsilon();

    // Normalizar M al intervalo [0, 2π]
    M = fmod(M, 2.0 * SAT_Const::pi);

    // Aproximación inicial según excentricidad
    if (e < 0.8)
        E = M;
    else
        E = SAT_Const::pi;

    // Primer paso del método de Newton-Raphson
    f = E - e * sin(E) - M;
    E = E - f / (1.0 - e * cos(E));

    // Iteración hasta convergencia
    while (abs(f) > 1e2 * epsilon) {
        f = E - e * sin(E) - M;
        E = E - f / (1.0 - e * cos(E));
        i++;
        if (i == maxit) {
            std::cout << "Convergence problems in EccAnom";
            exit(EXIT_FAILURE);
        }
    }

    return E;
}
