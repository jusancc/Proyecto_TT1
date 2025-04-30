#include "../include/Legendre.hpp"

tuple<Matrix &, Matrix &> Legendre(int n, int m, double fi)
{
    Matrix pnm = zeros(n + 1, m + 1);
    Matrix dpnm = zeros(n + 1, m + 1);

    pnm(1, 1) = 1;
    dpnm(1, 1) = 0;
    pnm(2, 2) = sqrt(3) * cos(fi);
    dpnm(2, 2) = -sqrt(3) * sin(fi);

    for (int i = 2; i < n; i++)
    {
        pnm(i + 1, i + 1) = sqrt((2 * i + 1) / (2 * i)) * cos(fi) * pnm(i, i);
    }
    for (int i = 2; i < n; i++)
    {
        dpnm(i + 1, i + 1) = sqrt((2 * i + 1) / (2 * i)) * ((cos(fi) * dpnm(i, i)) - (sin(fi) * pnm(i, i)));
    }
    for (int i = 1; i < n; i++)
    {
        pnm(i + 1, i) = sqrt(2 * i + 1) * sin(fi) * pnm(i, i);
    }
    for (int i = 1; i < n; i++)
    {
        dpnm(i + 1, i) = sqrt(2 * i + 1) * ((cos(fi) * pnm(i, i)) + (sin(fi) * dpnm(i, i)));
    }

    int j = 0;
    int k = 2;
    while (true)
    {
        for (int i = k; i <= n; i++)
        {
            double factor = sqrt((2 * i + 1.0) / ((i - j) * (i + j)));
            double term1 = sqrt(2 * i - 1.0) * sin(fi) * pnm(i, j + 1);
            double term2 = sqrt(((i + j - 1.0) * (i - j - 1.0)) / (2 * i - 3.0)) * pnm(i - 1, j + 1);

            pnm(i + 1, j + 1) = factor * (term1 - term2);
        }
        j++;
        k++;
        if (j > m)
            break;
    }

    j = 0;
    k = 2;

    while (true)
    {
        for (int i = k; i <= n; i++)
        {
            double factor = sqrt((2 * i + 1.0) / ((i - j) * (i + j)));
            double term1 = sqrt(2 * i - 1.0) * sin(fi) * dpnm(i, j + 1);
            double term2 = sqrt(2 * i - 1.0) * cos(fi) * pnm(i, j + 1);
            double term3 = sqrt(((i + j - 1.0) * (i - j - 1.0)) / (2 * i - 3.0)) * dpnm(i - 1, j + 1);

            dpnm(i + 1, j + 1) = factor * (term1 + term2 - term3);
        }
        j++;
        k++;
        if (j > m)
            break;
    }

    return tie(pnm,dpnm);
}