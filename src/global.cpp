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
 * @file global.cpp
 * @brief Implementación de las funciones que inicializan los parámetros globales
 *        para el modelo dinámico: datos EOP, coeficientes armónicos y efemérides.
 */
//------------------------------------------------------------------------------

#include "..\include\global.hpp"

Matrix eopdata;

//------------------------------------------------------------------------------
//  void eop19620101(int c)
//------------------------------------------------------------------------------
/**
 * @brief Carga los parámetros EOP (Earth Orientation Parameters) desde 1962.
 *
 * Lee el archivo `eop19620101.txt` y almacena los datos en la matriz `eopdata`.
 * 
 * @param c Número de columnas (registros) a cargar.
 *
 * @exception Termina el programa si no se puede abrir el archivo.
 */
//------------------------------------------------------------------------------
void eop19620101(int c){
    eopdata = zeros(13, c);
    FILE *fid = fopen("../data/eop19620101.txt", "r");
    if (fid == NULL) {
        printf("Fail open eop19620101.txt file\n");
        exit(EXIT_FAILURE);
    }

    for (int j = 1; j <= c; j++) {
        fscanf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &(eopdata(1, j)), &(eopdata(2, j)), &(eopdata(3, j)),
               &(eopdata(4, j)), &(eopdata(5, j)), &(eopdata(6, j)),
               &(eopdata(7, j)), &(eopdata(8, j)), &(eopdata(9, j)),
               &(eopdata(10, j)), &(eopdata(11, j)), &(eopdata(12, j)),
               &(eopdata(13, j)));
    }

    fclose(fid);
}

Matrix Cnm;
Matrix Snm;

//------------------------------------------------------------------------------
//  void GGM03S(int n)
//------------------------------------------------------------------------------
/**
 * @brief Carga los coeficientes armónicos Cnm y Snm desde el archivo GGM03S.txt.
 *
 * @param n Orden y grado máximo del desarrollo armónico (matriz n x n).
 *
 * @exception Termina el programa si no se puede abrir el archivo.
 */
//------------------------------------------------------------------------------
void GGM03S(int n){
    Cnm = zeros(n, n);
    Snm = zeros(n, n);
    FILE *fid = fopen("../data/GGM03S.txt", "r");
    if (fid == NULL) {
        printf("Fail open GGM03S.txt file\n");
        exit(EXIT_FAILURE);
    }

    double aux;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= i; j++) {
            fscanf(fid, "%lf %lf %lf %lf %lf %lf", &aux, &aux, &(Cnm(i, j)), &(Snm(i, j)), &aux, &aux);
        }
    }

    fclose(fid);
}

Matrix PC;

//------------------------------------------------------------------------------
//  void DE430Coeff(int f, int c)
//------------------------------------------------------------------------------
/**
 * @brief Carga los coeficientes de interpolación de las efemérides JPL DE430.
 *
 * @param f Número de filas a cargar.
 * @param c Número de columnas a cargar.
 *
 * @exception Termina el programa si no se puede abrir el archivo.
 */
//------------------------------------------------------------------------------
void DE430Coeff(int f, int c){
    PC = zeros(f, c);
    FILE *fid = fopen("../data/DE430Coeff.txt", "r");
    if (fid == NULL) {
        printf("Fail open DE430Coeff.txt file\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 1; i <= f; i++) {
        for (int j = 1; j <= c; j++) {
            fscanf(fid, "%lf", &(PC(i, j)));
        }
    }

    fclose(fid);
}

Param AuxParam;

//------------------------------------------------------------------------------
//  void initializeAuxParam()
//------------------------------------------------------------------------------
/**
 * @brief Inicializa los parámetros globales auxiliares (`AuxParam`) por defecto.
 *
 * Define las fechas iniciales, el orden de armónicos y qué cuerpos considerar
 * como perturbadores (Sol, Luna, Planetas).
 */
//------------------------------------------------------------------------------
void initializeAuxParam(){
    AuxParam.Mjd_UTC = 49746.1163541665;
    AuxParam.Mjd_TT = 49746.1170623147;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;
}
