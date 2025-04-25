#include "../include/Frac.hpp"

double Frac(double x){
    double res;
    res = x-floor(x);
    return res;
}