#ifndef _GLOBAL_
#define _GLOBAL_

#include "..\include\matrix.hpp"
#include <cmath>

typedef struct{
    double Mjd_UTC, Mjd_TT;
    int n, m, sun, moon, planets;
} Param;

extern Matrix eopdata;
extern Matrix Cnm, Snm;
extern Matrix PC;
extern Param AuxParam;

void eop19620101(int c);
void GGM03S(int c);
void DE430Coeff(int f, int c);
void initializeAuxParam();

#endif