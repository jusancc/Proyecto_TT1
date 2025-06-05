//$Header$
//------------------------------------------------------------------------------
//                                 SAT_Const
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file SAT_Const.hpp
 * @brief Definición de constantes astronómicas, físicas y matemáticas usadas en simulaciones orbitales.
 *
 * Esta clase contiene un conjunto de constantes utilizadas en mecánica celeste,
 * transformaciones de coordenadas y propagación orbital.
 */
//------------------------------------------------------------------------------

#ifndef _SAT_CONST_
#define _SAT_CONST_

#include <cmath>
#include <iostream>

/**
 * @class SAT_Const
 * @brief Clase estática que agrupa constantes fundamentales astronómicas y físicas.
 */
class SAT_Const {
public:
    //--------------------------------------------
    // Constantes matemáticas
    //--------------------------------------------
    static const double pi;     ///< π (3.1415...)
    static const double pi2;    ///< 2π
    static const double Rad;    ///< Conversión de grados a radianes
    static const double Deg;    ///< Conversión de radianes a grados
    static const double Arcs;   ///< Conversión de segundos de arco a radianes
    static const double eps;    ///< Tolerancia numérica

    //--------------------------------------------
    // Constantes generales
    //--------------------------------------------
    static const double MJD_J2000;   ///< Fecha juliana modificada del 1 de enero de 2000
    static const double T_B1950;     ///< Tiempo Besseliano 1950
    static const double c_light;     ///< Velocidad de la luz [km/s]
    static const double AU;          ///< Unidad astronómica [km]

    //--------------------------------------------
    // Parámetros físicos
    //--------------------------------------------
    static const double R_Earth;     ///< Radio medio de la Tierra [km]
    static const double f_Earth;     ///< Aplanamiento de la Tierra
    static const double R_Sun;       ///< Radio del Sol [km]
    static const double R_Moon;      ///< Radio de la Luna [km]
    static const double omega_Earth; ///< Velocidad angular de rotación de la Tierra [rad/s]

    //--------------------------------------------
    // Coeficientes gravitacionales [km^3/s^2]
    //--------------------------------------------
    static const double GM_Earth;    ///< Tierra
    static const double GM_Sun;      ///< Sol
    static const double GM_Moon;     ///< Luna
    static const double GM_Mercury;  ///< Mercurio
    static const double GM_Venus;    ///< Venus
    static const double GM_Mars;     ///< Marte
    static const double GM_Jupiter;  ///< Júpiter
    static const double GM_Saturn;   ///< Saturno
    static const double GM_Uranus;   ///< Urano
    static const double GM_Neptune;  ///< Neptuno
    static const double GM_Pluto;    ///< Plutón

    //--------------------------------------------
    // Radiación solar
    //--------------------------------------------
    static const double P_Sol;       ///< Presión de radiación solar a 1 UA [N/m^2]
};

#endif // _SAT_CONST_
