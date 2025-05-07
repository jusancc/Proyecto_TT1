#ifndef _GLOBAL_
#define _GLOBAL_

#include "..\include\matrix.hpp"
#include <cmath>

extern Matrix eopdata;
extern Matrix Cnm, Snm;

void eop19620101(int c);
void GGM03S(int c);

#endif