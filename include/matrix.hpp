//$Header$
//------------------------------------------------------------------------------
//                                   Matrix
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file matrix.hpp
 * @brief Declaración de la clase Matrix y sus operaciones auxiliares.
 *
 * Esta clase implementa una matriz dinámica con operadores y métodos para álgebra
 * lineal básica, diseñada para sustituir a librerías externas en entornos cerrados.
 */
//------------------------------------------------------------------------------

#ifndef _MATRIX_
#define _MATRIX_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

/**
 * @class Matrix
 * @brief Implementación básica de matrices dinámicas con operaciones elementales.
 */
class Matrix {
private:
    double **data; ///< Puntero a los datos de la matriz

public:
    int n_row;     ///< Número de filas
    int n_column;  ///< Número de columnas

    // Constructores
    Matrix(); ///< Constructor por defecto
    Matrix(const int v_size); ///< Constructor para vectores columna
    Matrix(const int n_row, const int n_column); ///< Constructor con dimensiones

    // Operadores miembro
    double& operator()(const int n); ///< Acceso para vectores
    double& operator()(const int row, const int column); ///< Acceso para matrices

    Matrix& operator+(Matrix &m); ///< Suma de matrices
    Matrix& operator+(double s); ///< Suma escalar
    Matrix& operator-(Matrix &m); ///< Resta de matrices
    Matrix& operator-(double s); ///< Resta escalar
    Matrix& operator*(Matrix &m); ///< Multiplicación de matrices
    Matrix& operator*(double s); ///< Multiplicación escalar
    Matrix& operator/(Matrix &m); ///< División de matrices (elemento a elemento)
    Matrix& operator/(double s); ///< División escalar
    Matrix& operator=(Matrix &m); ///< Asignación

    // Métodos
    double norm(); ///< Norma euclídea
    double dot(Matrix &m); ///< Producto escalar

    Matrix& extract_row(int row); ///< Extrae una fila como vector fila
    Matrix& extract_column(int column); ///< Extrae una columna como vector columna
    void assign_row(int row, Matrix &v); ///< Asigna una fila
    void assign_column(int column, Matrix &v); ///< Asigna una columna

    Matrix& extract_vector(int from, int to); ///< Extrae un subconjunto como vector

    Matrix& inv(); ///< Inversa de matriz (solo cuadradas)

    // Sobrecarga de operador de salida
    friend ostream& operator<<(ostream &o, Matrix &m);
};

/**
 * @brief Sobrecarga de << para imprimir matrices.
 * @param o Flujo de salida
 * @param m Matriz a imprimir
 * @return Referencia al flujo de salida
 */
ostream& operator<<(ostream &o, Matrix &m);

// Funciones auxiliares

/**
 * @brief Genera una matriz de ceros.
 * @param n_row Filas
 * @param n_column Columnas
 * @return Matriz llena de ceros
 */
Matrix& zeros(const int n_row, const int n_column);

/**
 * @brief Genera un vector columna de ceros.
 * @param n Tamaño del vector
 * @return Vector de ceros
 */
Matrix& zeros(const int n);

/**
 * @brief Genera una matriz identidad.
 * @param n Tamaño de la matriz cuadrada
 * @return Matriz identidad
 */
Matrix& eye(const int n);

/**
 * @brief Calcula la transpuesta de una matriz.
 * @param m Matriz a transponer
 * @return Matriz transpuesta
 */
Matrix& transpose(Matrix &m);

/**
 * @brief Une dos vectores columna en uno solo.
 * @param v1 Primer vector
 * @param v2 Segundo vector
 * @return Vector concatenado
 */
Matrix& union_vector(Matrix &v1, Matrix &v2);

/**
 * @brief Calcula el producto vectorial entre dos vectores 3D.
 * @param v Primer vector
 * @param w Segundo vector
 * @return Producto vectorial
 */
Matrix& v_cross(Matrix &v, Matrix &w);

#endif  // _MATRIX_
