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
 * @file SAT_Const.cpp
 * @brief Implementación de las constantes definidas en SAT_Const.hpp.
 *
 * Este fichero define el valor de todas las constantes físicas, astronómicas
 * y matemáticas necesarias para la propagación orbital y transformaciones
 * de sistemas de referencia en el proyecto.
 */
//------------------------------------------------------------------------------

#include "../include/SAT_Const.hpp"
#include <limits> 

//--------------------------------------------
// Constantes matemáticas
//--------------------------------------------

const double SAT_Const::pi       = 3.14159265358979324;                     ///< π
const double SAT_Const::pi2      = 2 * SAT_Const::pi;                       ///< 2π
const double SAT_Const::Rad      = SAT_Const::pi / 180;                     ///< grados → radianes
const double SAT_Const::Deg      = 180 / SAT_Const::pi;                     ///< radianes → grados
const double SAT_Const::Arcs     = 3600 * 180 / SAT_Const::pi;              ///< segundos de arco → radianes
const double SAT_Const::eps      = std::numeric_limits<double>::epsilon(); ///< tolerancia numérica

//--------------------------------------------
// Constantes generales
//--------------------------------------------

const double SAT_Const::MJD_J2000 = 51544.5;                ///< Fecha juliana modificada (J2000.0)
const double SAT_Const::T_B1950   = -0.500002108;           ///< Tiempo Besseliano 1950
const double SAT_Const::c_light   = 299792458.0;            ///< Velocidad de la luz [m/s]
const double SAT_Const::AU        = 149597870700.0;         ///< Unidad astronómica [m]

//--------------------------------------------
// Parámetros físicos
//--------------------------------------------

const double SAT_Const::R_Earth   = 6378.1363e3;            ///< Radio medio terrestre [m]
const double SAT_Const::f_Earth   = 1.0 / 298.257223563;    ///< Aplanamiento de la Tierra
const double SAT_Const::R_Sun     = 696000e3;               ///< Radio del Sol [m]
const double SAT_Const::R_Moon    = 1738e3;                 ///< Radio de la Luna [m]
const double SAT_Const::omega_Earth = 15.04106717866910 / 3600 * SAT_Const::Rad; ///< Velocidad angular de la Tierra [rad/s]

//--------------------------------------------
// Coeficientes gravitacionales [m³/s²]
//--------------------------------------------

const double SAT_Const::GM_Earth    = 398600.435436e9;       ///< GM de la Tierra
const double SAT_Const::GM_Sun      = 132712440041.939400e9; ///< GM del Sol
const double SAT_Const::GM_Moon     = SAT_Const::GM_Earth / 81.30056907419062; ///< GM de la Luna
const double SAT_Const::GM_Mercury  = 22031.780000e9;        ///< GM de Mercurio
const double SAT_Const::GM_Venus    = 324858.592000e9;       ///< GM de Venus
const double SAT_Const::GM_Mars     = 42828.375214e9;        ///< GM de Marte
const double SAT_Const::GM_Jupiter  = 126712764.800000e9;    ///< GM de Júpiter
const double SAT_Const::GM_Saturn   = 37940585.200000e9;     ///< GM de Saturno
const double SAT_Const::GM_Uranus   = 5794548.600000e9;      ///< GM de Urano
const double SAT_Const::GM_Neptune  = 6836527.100580e9;      ///< GM de Neptuno
const double SAT_Const::GM_Pluto    = 977.0000000000009e9;   ///< GM de Plutón

//--------------------------------------------
// Radiación solar
//--------------------------------------------

const double SAT_Const::P_Sol = 1367 / SAT_Const::c_light; ///< Presión de radiación solar a 1 UA [N/m²]
