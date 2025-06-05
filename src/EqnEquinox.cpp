//$Header$
//------------------------------------------------------------------------------
//                              EqnEquinox
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file EqnEquinox.cpp
 * @brief Implementación de la función que calcula la ecuación de los equinoccios.
 */
//------------------------------------------------------------------------------

#include "../include/EqnEquinox.hpp"

//------------------------------------------------------------------------------
//  double EqnEquinox(double Mjd_TT)
//------------------------------------------------------------------------------
/**
 * @brief Calcula la ecuación de los equinoccios para una fecha dada.
 *
 * Esta función devuelve el valor de la ecuación de los equinoccios (EqE) en radianes,
 * que representa la diferencia entre el ángulo horario aparente y medio de Greenwich.
 * Se calcula como:
 * 
 *     EqE = Δψ · cos(ε)
 *
 * donde Δψ es la nutación en longitud y ε es la oblicuidad media de la eclíptica.
 * 
 * @param Mjd_TT Fecha en Tiempo Terrestre (Modified Julian Date TT).
 * @return Ecuación de los equinoccios en radianes.
 */
//------------------------------------------------------------------------------
double EqnEquinox(double Mjd_TT){
    double EqE, dpsi, deps;
    tie(dpsi, deps) = NutAngles(Mjd_TT);  // Obtener ángulos de nutación
    EqE = dpsi * cos(MeanObliquity(Mjd_TT));  // Aplicar fórmula
    return EqE;
}
