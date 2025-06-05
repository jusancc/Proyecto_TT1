//$Header$
//------------------------------------------------------------------------------
//                                 Mjday_TDB
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file Mjday_TDB.hpp
 * @brief Declaración de la función que convierte de MJD_TT a MJD_TDB.
 *
 * Esta función aplica una corrección periódica para transformar el Tiempo Terrestre (TT)
 * en Tiempo Dinámico Baricéntrico (TDB) en formato Modified Julian Date.
 */
//------------------------------------------------------------------------------

#ifndef _MJDAY_TDB_
#define _MJDAY_TDB_

#include <cmath>

/**
 * @brief Convierte Tiempo Terrestre (TT) a Tiempo Dinámico Baricéntrico (TDB).
 * 
 * @param Mjd_TT Tiempo Terrestre en MJD.
 * @return Tiempo TDB en MJD.
 */
double Mjday_TDB(double Mjd_TT);

#endif
