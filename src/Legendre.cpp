#include "../include/Legendre.hpp"
#include <cmath>
#include <tuple>
using namespace std;

tuple<Matrix&, Matrix&> Legendre(int n, int m, double fi) {
    // Inicializar matrices con el tamaño adecuado
    Matrix& pnm = zeros(n + 1, m + 1);
    Matrix& dpnm = zeros(n + 1, m + 1);

    // Asignaciones iniciales
    pnm(1, 1) = 1;
    dpnm(1, 1) = 0;
    pnm(2, 2) = sqrt(3) * cos(fi);
    dpnm(2, 2) = -sqrt(3) * sin(fi);

    // Coeficientes diagonales
    for (int i = 2; i <= n; i++) {  
        pnm(i + 1, i + 1) = sqrt((2.0 * i + 1) / (2.0 * i)) * cos(fi) * pnm(i, i);
    }
    for (int i = 2; i <= n; i++) {
        dpnm(i + 1, i + 1) = sqrt((2.0 * i + 1) / (2.0 * i)) * ((cos(fi) * dpnm(i, i)) - (sin(fi) * pnm(i, i)));
    }

    // Coeficientes horizontales (primer paso)
    for (int i = 1; i <= n; i++) {
        pnm(i + 1, i) = sqrt(2.0 * i + 1) * sin(fi) * pnm(i, i);
    }
    for (int i = 1; i <= n; i++) {
        dpnm(i + 1, i) = sqrt(2.0 * i + 1) * ((cos(fi) * pnm(i, i)) + (sin(fi) * dpnm(i, i)));
    }

    // Coeficientes horizontales (segundo paso)
    int j = 0;
    int k = 2;
    while (true) {
        for (int i = k; i <= n; i++) {
            pnm(i + 1, j + 1) = sqrt((2.0 * i + 1) / ((i - j) * (i + j))) * (
                (sqrt(2.0 * i - 1) * sin(fi) * pnm(i, j + 1)) 
                - (sqrt(((i + j - 1) * (i - j - 1)) / (2.0 * i - 3)) * pnm(i - 1, j + 1))
            );
        }
        j = j + 1;
        k = k + 1;
        if (j > m) break;
    }

    // Coeficientes de la derivada horizontal (segundo paso)
    j = 0;
    k = 2;
    while (true) {
        for (int i = k; i <= n; i++) {
            dpnm(i + 1, j + 1) = sqrt((2.0 * i + 1) / ((i - j) * (i + j))) * (
                (sqrt(2.0 * i - 1) * sin(fi) * dpnm(i, j + 1)) 
                + (sqrt(2.0 * i - 1) * cos(fi) * pnm(i, j + 1)) 
                - (sqrt(((i + j - 1) * (i - j - 1)) / (2.0 * i - 3)) * dpnm(i - 1, j + 1))
            );
        }
        j = j + 1;
        k = k + 1;
        if (j > m) break;
    }

    return tie(pnm, dpnm);
}
