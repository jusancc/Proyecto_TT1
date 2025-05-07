#include "../include/EccAnom.hpp"

double EccAnom(double M, double e){
    int maxit = 15;
    int i = 1;
    double E, f;
    const double epsilon = 1e2 * std::numeric_limits<double>::epsilon();
    M = fmod(M, 2.0*SAT_Const::pi);

    if (e < 0.8)
    {
       E = M; 
    } else{
        E = SAT_Const::pi;
    }

    f = E - e*sin(E) - M;
    E = E - f / (1.0 - e*cos(E));

    while (abs(f) > 1e2*epsilon)
    {
        f = E - e*sin(E) - M;
        E = E - f / (1.0 - e*cos(E));
        i = i+1;
        if (i == maxit)
        {
            std:cout << "Convergence problems in EccAnom";
            exit(EXIT_FAILURE);
        }
    }
    return E;
}