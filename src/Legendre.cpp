//$Header$
//------------------------------------------------------------------------------
//                                Legendre
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file Legendre.cpp
 * @brief Implementación de los polinomios asociados de Legendre y sus derivadas.
 *
 * Esta función calcula los polinomios asociados de Legendre (normalizados) y sus
 * derivadas con respecto al ángulo `fi`, necesarios para el cálculo del potencial
 * gravitacional armónico de la Tierra.
 */
//------------------------------------------------------------------------------

#include "../include/Legendre.hpp"
#include <cmath>
#include <tuple>
using namespace std;

/**
 * @brief Calcula los polinomios de Legendre y sus derivadas con normalización de Schmidt.
 *
 * @param n Grado máximo de la expansión.
 * @param m Orden máximo de la expansión.
 * @param fi Ángulo de latitud geocéntrica (en radianes).
 * @return Una tupla con:
 *         - `pnm`: Matriz (n+1)x(m+1) con los valores de los polinomios \( P_n^m(\sin(\phi)) \)
 *         - `dpnm`: Matriz (n+1)x(m+1) con las derivadas respecto a `fi` de los polinomios.
 */
tuple<Matrix&, Matrix&> Legendre(int n, int m, double fi) {
    Matrix& pnm = zeros(n + 1, m + 1);    // Polinomios de Legendre
    Matrix& dpnm = zeros(n + 1, m + 1);   // Derivadas de los polinomios

    // Inicialización
    pnm(1, 1) = 1;
    dpnm(1, 1) = 0;
    pnm(2, 2) = sqrt(3) * cos(fi);
    dpnm(2, 2) = -sqrt(3) * sin(fi);

    // Componentes diagonales
    for (int i = 2; i <= n; i++) {
        pnm(i + 1, i + 1) = sqrt((2.0 * i + 1) / (2.0 * i)) * cos(fi) * pnm(i, i);
        dpnm(i + 1, i + 1) = sqrt((2.0 * i + 1) / (2.0 * i)) *
                             (cos(fi) * dpnm(i, i) - sin(fi) * pnm(i, i));
    }

    // Primera fila subdiagonal
    for (int i = 1; i <= n; i++) {
        pnm(i + 1, i) = sqrt(2.0 * i + 1) * sin(fi) * pnm(i, i);
        dpnm(i + 1, i) = sqrt(2.0 * i + 1) *
                         (cos(fi) * pnm(i, i) + sin(fi) * dpnm(i, i));
    }

    // Recurrencia general
    int j = 0, k = 2;
    while (j <= m) {
        for (int i = k; i <= n; i++) {
            pnm(i + 1, j + 1) = sqrt((2.0 * i + 1) / ((i - j) * (i + j))) *
                (sqrt(2.0 * i - 1) * sin(fi) * pnm(i, j + 1)
                 - sqrt(((i + j - 1) * (i - j - 1)) / (2.0 * i - 3)) * pnm(i - 1, j + 1));
        }
        j++; k++;
    }

    // Recurrencia de derivadas
    j = 0; k = 2;
    while (j <= m) {
        for (int i = k; i <= n; i++) {
            dpnm(i + 1, j + 1) = sqrt((2.0 * i + 1) / ((i - j) * (i + j))) *
                (sqrt(2.0 * i - 1) * sin(fi) * dpnm(i, j + 1)
                 + sqrt(2.0 * i - 1) * cos(fi) * pnm(i, j + 1)
                 - sqrt(((i + j - 1) * (i - j - 1)) / (2.0 * i - 3)) * dpnm(i - 1, j + 1));
        }
        j++; k++;
    }

    return tie(pnm, dpnm);
}
