//$Header$
//------------------------------------------------------------------------------
//                                  global
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file global.hpp
 * @brief Declaración de variables y funciones globales usadas en dinámica orbital.
 *
 * Este archivo define estructuras de parámetros, datos globales y funciones auxiliares
 * necesarias para inicializar y cargar coeficientes y parámetros astronómicos.
 */
//------------------------------------------------------------------------------

#ifndef _GLOBAL_
#define _GLOBAL_

#include "..\include\matrix.hpp"  ///< Clase Matrix personalizada
#include <cmath>                  ///< Funciones matemáticas estándar

/**
 * @struct Param
 * @brief Estructura con parámetros astronómicos auxiliares.
 */
typedef struct {
    double Mjd_UTC;     ///< Fecha en tiempo universal (UTC) en MJD
    double Mjd_TT;      ///< Fecha en tiempo terrestre (TT) en MJD
    int n;              ///< Orden máximo del modelo armónico
    int m;              ///< Grado máximo del modelo armónico
    int sun;            ///< 1 si se considera la perturbación solar
    int moon;           ///< 1 si se considera la perturbación lunar
    int planets;        ///< 1 si se consideran los planetas
} Param;

//------------------------------------------------------------------------------
// Variables globales
//------------------------------------------------------------------------------

extern Matrix eopdata;  ///< Tabla de parámetros IERS (EOP)
extern Matrix Cnm;      ///< Coeficientes armónicos Cnm
extern Matrix Snm;      ///< Coeficientes armónicos Snm
extern Matrix PC;       ///< Coeficientes de precesión/nutación
extern Param AuxParam;  ///< Parámetros auxiliares de configuración

//------------------------------------------------------------------------------
// Funciones auxiliares de inicialización
//------------------------------------------------------------------------------

/**
 * @brief Carga datos de parámetros EOP desde 1962.
 * @param c Cantidad de registros a cargar.
 */
void eop19620101(int c);

/**
 * @brief Carga el modelo de campo gravitacional GGM03S.
 * @param c Número de coeficientes a cargar.
 */
void GGM03S(int c);

/**
 * @brief Carga coeficientes de efemérides JPL DE430.
 * @param f Fila inicial.
 * @param c Cantidad de columnas a cargar.
 */
void DE430Coeff(int f, int c);

/**
 * @brief Inicializa los parámetros globales definidos en AuxParam.
 */
void initializeAuxParam();

#endif  // _GLOBAL_
