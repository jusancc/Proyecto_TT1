//$Header$
//------------------------------------------------------------------------------
//                                Cheb3D
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file Cheb3D.cpp
 * @brief Implementación de la función que evalúa una serie de Chebyshev en 3D.
 */
//------------------------------------------------------------------------------

#include "../include/Cheb3D.hpp"

//------------------------------------------------------------------------------
//  Matrix& Cheb3D(double t, int N, double Ta, double Tb,
//                 Matrix Cx, Matrix Cy, Matrix Cz)
//------------------------------------------------------------------------------
/**
 * @brief Evalúa una serie de interpolación de Chebyshev en 3D.
 *
 * Esta función calcula la posición 3D en el tiempo `t` mediante una serie de 
 * Chebyshev de segundo orden. Usa el método de Clenshaw para evaluar los polinomios
 * y aplica coeficientes para las tres componentes espaciales (X, Y, Z).
 *
 * @param t Tiempo de evaluación.
 * @param N Número de coeficientes (orden del polinomio + 1).
 * @param Ta Inicio del intervalo de interpolación.
 * @param Tb Fin del intervalo de interpolación.
 * @param Cx Coeficientes Chebyshev para la coordenada X.
 * @param Cy Coeficientes Chebyshev para la coordenada Y.
 * @param Cz Coeficientes Chebyshev para la coordenada Z.
 * @return Referencia a un vector (1x3) con la posición interpolada [x, y, z].
 */
//------------------------------------------------------------------------------
Matrix& Cheb3D(double t, int N, double Ta, double Tb, Matrix Cx, Matrix Cy, Matrix Cz){
    if (t < Ta || Tb < t) {
        std::cout << "Time out of range in Cheb3D";
        exit(EXIT_FAILURE);
    }

    double tau = (2 * t - Ta - Tb) / (Tb - Ta);  // Transformación a [-1, 1]

    Matrix f1 = zeros(1, 3);       // término T(i)
    Matrix f2 = zeros(1, 3);       // término T(i-1)
    Matrix old_f1 = zeros(1, 3);   // buffer temporal

    // Algoritmo de Clenshaw inverso
    for (int i = N; i >= 2; i--) {
        old_f1 = f1;
        f1 = f1 * 2 * tau - f2;

        f1(1,1) += Cx(1, i);
        f1(1,2) += Cy(1, i);
        f1(1,3) += Cz(1, i);

        f2 = old_f1;
    }

    Matrix &ChebApp = zeros(1, 3);  // Resultado final

    ChebApp = f1 * tau - f2;
    ChebApp(1,1) += Cx(1, 1);
    ChebApp(1,2) += Cy(1, 1);
    ChebApp(1,3) += Cz(1, 1);

    return ChebApp;
}
