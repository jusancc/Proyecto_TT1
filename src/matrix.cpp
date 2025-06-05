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
 * @file matrix.cpp
 * @brief Implementación de la clase Matrix y sus operaciones básicas.
 *
 * Esta clase maneja matrices dinámicas con operaciones comunes de álgebra lineal,
 * útiles en aplicaciones científicas donde no se desea depender de librerías externas.
 */
//------------------------------------------------------------------------------

#include "..\include\matrix.hpp"

//--------------------------------------------
// Constructores
//--------------------------------------------

/**
 * @brief Constructor por defecto.
 */

Matrix::Matrix(){
	this -> n_row = 0;
    this -> n_column = 0;
    this -> data = nullptr;
}

/**
 * @brief Constructor de matriz cuadrada (v_size x v_size).
 * @param v_size Tamaño del vector/cuadrado.
 */
Matrix::Matrix(const int v_size){
	if (v_size <= 0) {
        cout << "Matrix create: error in v_size\n";
        exit(EXIT_FAILURE);
    }
    
    this->n_row = v_size;
    this->n_column = v_size;
    this->data = (double **)malloc(v_size * sizeof(double *));
    
    if (this->data == NULL) {
        cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
    }
    
    for (int i = 0; i < v_size; i++) {
        this->data[i] = (double *)malloc(v_size * sizeof(double));
        if (this->data[i] == NULL) {
            cout << "Matrix create: error in row allocation\n";
            exit(EXIT_FAILURE);
        }
    }
}

/**
 * @brief Constructor con número de filas y columnas especificado.
 * @param n_row Número de filas.
 * @param n_column Número de columnas.
 */
Matrix::Matrix(const int n_row, const int n_column) {
    if (n_row <= 0 || n_column <= 0) {
		cout << "Matrix create: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	this->n_row = n_row;
	this->n_column = n_column;
	this->data = (double **) malloc(n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < n_row; i++) {
		this->data[i] = (double *) malloc(n_column*sizeof(double));
	}
}

//--------------------------------------------
// Accesores
//--------------------------------------------

/**
 * @brief Accede a un elemento de la matriz usando un único índice.
 * @param n Índice plano (por filas).
 * @return Referencia al valor.
 */
double& Matrix::operator () (const int n) {
    if (n <= 0 || n > this->n_row*this->n_column) {
		cout << "Matrix: error in row/column\n";
        exit(EXIT_FAILURE);
	}
	
	return this->data[(n - 1)/this->n_column][(n - 1)%this->n_column];
}

/**
 * @brief Accede a un elemento usando fila y columna.
 * @param row Fila (1-based)
 * @param column Columna (1-based)
 * @return Referencia al valor.
 */
double& Matrix::operator () (const int row, const int column) {
	if (row <= 0 || row > this->n_row || column <= 0 || column > this->n_column) {
		cout << "Matrix get: error in row/column\n";
        exit(EXIT_FAILURE);
	}
	
	return this->data[row - 1][column - 1];
}

//--------------------------------------------
// Operaciones básicas
//--------------------------------------------

/**
 * @brief Suma matricial.
 */
Matrix& Matrix::operator + (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sum: error in n_row+n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) + m(i,j);
		}
	}
	
	return *m_aux;
}

/**
 * @brief Resta matricial.
 */
Matrix& Matrix::operator - (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sub: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) - m(i,j);
		}
	}
	
	return *m_aux;
}

ostream& operator << (ostream &o, Matrix &m) {
	for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++)
			printf("%5.20lf ", m(i,j));
        o << "\n";
    }
	
    return o;
}

/**
 * @brief Crea una matriz de ceros.
 */
Matrix& zeros(const int n_row, const int n_column) {
	Matrix *m_aux = new Matrix(n_row, n_column);
	
	for(int i = 1; i <= n_row; i++) {
		for(int j = 1; j <= n_column; j++) {
			(*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
}

/**
 * @brief Crea una matriz cuadrada de ceros.
 */
Matrix& zeros(const int n){
	if (n <= 0) {
        cout << "zeros: error in size n\n";
        exit(EXIT_FAILURE);
    }
    
    Matrix *m_aux = new Matrix(n, n);
    
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            (*m_aux)(i, j) = 0;
        }
    }
    
    return *m_aux;
}

/**
 * @brief Producto matricial.
 */
Matrix& Matrix::operator * (Matrix &m) {
    if (this->n_column != m.n_row) {
        cout << "Matrix product: error in n_row*m_column\n";
        exit(EXIT_FAILURE);
    }

    Matrix *m_aux = new Matrix(this->n_row, m.n_column);

    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= m.n_column; j++) {
            double sum = 0;
            for(int k = 1; k <= this->n_column; k++) {
                sum += (*this)(i, k) * m(k, j);
            }
            (*m_aux)(i, j) = sum;
        }
    }

    return *m_aux;
}

/**
 * @brief División matricial (multiplicación por inversa).
 */
Matrix& Matrix::operator / (Matrix &m){
    // Verificar que las dimensiones sean compatibles
    if (this->n_column != m.n_row || m.n_row != m.n_column) {
        cout << "Matrix division: Incompatible dimensions or divisor is not square\n";
        exit(EXIT_FAILURE);
    }

    // Calcular la inversa de m
    Matrix& inv_m = m.inv();

    // Multiplicar (*this) por inv_m (A * B^{-1})
    Matrix *result = new Matrix(this->n_row, m.n_column);
    *result = (*this) * inv_m;

    
    return *result;
}

/**
 * @brief Asignación por copia profunda.
 */
Matrix& Matrix::operator = (Matrix &m){
	if (this == &m) {
        return *this;
    }
    
    this->n_row = m.n_row;
    this->n_column = m.n_column;
    
    this->data = (double **) malloc(this->n_row * sizeof(double *));
    if (this->data == NULL) {
        cout << "Matrix assignment: error in data allocation\n";
        exit(EXIT_FAILURE);
    }
    
    for (int i = 0; i < this->n_row; i++) {
        this->data[i] = (double *) malloc(this->n_column * sizeof(double));
        if (this->data[i] == NULL) {
            cout << "Matrix assignment: error in row allocation\n";
            exit(EXIT_FAILURE);
        }
        for (int j = 0; j < this->n_column; j++) {
            this->data[i][j] = m.data[i][j];
        }
    }
    
    return *this;
}

/**
 * @brief Crea una matriz identidad.
 */
Matrix& eye(const int n) {
    if (n <= 0) {
        cout << "eye: error in size n\n";
        exit(EXIT_FAILURE);
    }
    
    Matrix *m_aux = new Matrix(n, n);
    
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i == j) {
                (*m_aux)(i, j) = 1.0;
            } else {
                (*m_aux)(i, j) = 0.0;
            }
        }
    }
    
    return *m_aux;
}

/**
 * @brief Devuelve la transpuesta de una matriz.
 */
Matrix& transpose(Matrix &m) {
    Matrix *m_aux = new Matrix(m.n_column, m.n_row);
    
    for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++) {
            (*m_aux)(j, i) = m(i, j);
        }
    }
    
    return *m_aux;
}

/**
 * @brief Suma escalar.
 */
Matrix& Matrix::operator + (double s){
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
    
    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i, j) = (*this)(i, j) + s;
        }
    }
    
    return *m_aux;
}

/**
 * @brief Resta escalar.
 */
Matrix& Matrix::operator - (double s){
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
    
    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i, j) = (*this)(i, j) - s;
        }
    }
    
    return *m_aux;
}

/**
 * @brief Producto escalar por matriz.
 */
Matrix& Matrix::operator * (double s){
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
    
    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i, j) = (*this)(i, j) * s;
        }
    }
    
    return *m_aux;
}

/**
 * @brief División escalar.
 */
Matrix& Matrix::operator / (double s){
	if (s == 0) {
        cout << "Matrix division: division by zero detected\n";
        exit(EXIT_FAILURE);
    }
    
    Matrix *m_aux = new Matrix(this->n_row, this->n_column);
    
    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i, j) = (*this)(i, j) / s;
        }
    }
    
    return *m_aux;
}

//--------------------------------------------
// Métodos auxiliares
//--------------------------------------------

/**
 * @brief Calcula la norma euclídea.
 */
double Matrix::norm(){
    double r = 0.0;
    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            r += pow((*this)(i, j), 2);
        }
    }
    return sqrt(r);
}

/**
 * @brief Producto escalar (dot product).
 */
double Matrix::dot(Matrix &m){
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
        cout << "Matrix dot: error in dimensions\n";
        exit(EXIT_FAILURE);
    }
    
    double r = 0.0;
    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            r += (*this)(i, j) * m(i, j);
        }
    }
    return r;
}

/**
 * @brief Producto vectorial entre dos vectores 3D.
 */
Matrix& v_cross(Matrix &v, Matrix &w){
	if (v.n_row != 3 || v.n_column != 1 || w.n_row != 3 || w.n_column != 1) {
        cout << "v_cross: vectors must be 3x1\n";
        exit(EXIT_FAILURE);
    }
    Matrix *result = new Matrix(3, 1);
    (*result)(1, 1) = v(2, 1) * w(3, 1) - v(3, 1) * w(2, 1);
    (*result)(2, 1) = v(3, 1) * w(1, 1) - v(1, 1) * w(3, 1);
    (*result)(3, 1) = v(1, 1) * w(2, 1) - v(2, 1) * w(1, 1);
    return *result;
}

/**
 * @brief Extrae una fila como matriz (1 x n).
 */
Matrix& Matrix::extract_row(int row) {
    if (row <= 0 || row > this->n_row) {
        cout << "extract_row: error in row index\n";
        exit(EXIT_FAILURE);
    }
    
    Matrix *result = new Matrix(1, this->n_column);
    
    for (int j = 1; j <= this->n_column; j++) {
        (*result)(1, j) = (*this)(row, j);
    }
    
    return *result;
}

/**
 * @brief Extrae una columna como matriz (n x 1).
 */
Matrix& Matrix::extract_column(int column) {
    if (column <= 0 || column > this->n_column) {
        cout << "extract_column: error in column index\n";
        exit(EXIT_FAILURE);
    }
    
    Matrix *result = new Matrix(this->n_row, 1);
    
    for (int i = 1; i <= this->n_row; i++) {
        (*result)(i, 1) = (*this)(i, column);
    }
    
    return *result;
}

/**
 * @brief Asigna valores a una fila desde un vector.
 */
void Matrix::assign_row(int row, Matrix &v) {
    if (row <= 0 || row > this->n_row) {
        cout << "assign_row: error in row index\n";
        exit(EXIT_FAILURE);
    }
    if (v.n_row != 1 || v.n_column != this->n_column) {
        cout << "assign_row: vector dimensions mismatch\n";
        exit(EXIT_FAILURE);
    }
    
    for (int j = 1; j <= this->n_column; j++) {
        (*this)(row, j) = v(1, j);
    }
}

/**
 * @brief Asigna valores a una columna desde un vector.
 */
void Matrix::assign_column(int column, Matrix &v) {
    if (column <= 0 || column > this->n_column) {
        cout << "assign_column: error in column index\n";
        exit(EXIT_FAILURE);
    }
    if (v.n_row != this->n_row || v.n_column != 1) {
        cout << "assign_column: vector dimensions mismatch\n";
        exit(EXIT_FAILURE);
    }
    
    for (int i = 1; i <= this->n_row; i++) {
        (*this)(i, column) = v(i, 1);
    }
}

/**
 * @brief Concatena dos vectores.
 */
Matrix& union_vector(Matrix &v1, Matrix &v2) {
    if (v1.n_row != 1 || v2.n_row != 1){
		cout << "Matrix union_vector: error in n_row\n";
		exit(EXIT_FAILURE);
	}

	int tamano = v1.n_column + v2.n_column;
	Matrix *m_aux = new Matrix(1, tamano);

	int aux = 1;
	for (int i = 1; i <= v1.n_column; i++){
		(*m_aux)(1, aux) = (v1)(1, i);
		aux++;
	}
	for (int i = 1; i <= v2.n_column; i++){
		(*m_aux)(1, aux) = v2(1, i);
		aux++;
	}

	return *m_aux;
}

/**
 * @brief Extrae un vector desde un valor hasta otro.
 */
Matrix& Matrix::extract_vector(int from, int to) {
    int total_size = this->n_row * this->n_column;
    if (from <= 0 || from > total_size || to < from || to > total_size) {
        cout << "extract_vector: invalid range from " << from << " to " << to << "\n";
        exit(EXIT_FAILURE);
    }
    
    int length = to - from + 1;
    Matrix *result = new Matrix(1, length);
    
    for (int n = from; n <= to; n++) {
        (*result)(1, n - from + 1) = (*this)(n);
    }
    
    return *result;
}

/**
 * @brief Calcula la inversa de una matriz.
 */
Matrix& Matrix::inv() {
    // Verificar que la matriz sea cuadrada
    if (this->n_row != this->n_column) {
        cout << "Matrix inverse: Matrix must be square\n";
        exit(EXIT_FAILURE);
    }

    int n = this->n_row;
    Matrix *inv = new Matrix(n, n); // Matriz inversa (a calcular)
    Matrix temp(*this);             // Copia de la matriz original

    // Inicializar inv como matriz identidad
    
    *inv = eye(n);

    // Eliminación Gauss-Jordan con pivoteo parcial
    for (int col = 1; col <= n; col++) {
        // Pivoteo parcial: buscar el máximo en la columna actual
        int max_row = col;
        double max_val = fabs(temp(col, col));

        for (int row = col + 1; row <= n; row++) {
            if (fabs(temp(row, col)) > max_val) {
                max_val = fabs(temp(row, col));
                max_row = row;
            }
        }

        // Si el pivote es cero, la matriz es singular
        if (max_val < 1e-10) {
            cout << "Matrix inverse: Matrix is singular (non-invertible)\n";
            exit(EXIT_FAILURE);
        }

        // Intercambiar filas si es necesario
        if (max_row != col) {
            for (int j = 1; j <= n; j++) {
                std::swap(temp(col, j), temp(max_row, j));
                std::swap((*inv)(col, j), (*inv)(max_row, j));
            }
        }

        // Normalizar la fila del pivote
        double pivot = temp(col, col);
        for (int j = 1; j <= n; j++) {
            temp(col, j) /= pivot;
            (*inv)(col, j) /= pivot;
        }

        // Eliminación hacia adelante y atrás
        for (int row = 1; row <= n; row++) {
            if (row != col) {
                double factor = temp(row, col);
                for (int j = 1; j <= n; j++) {
                    temp(row, j) -= factor * temp(col, j);
                    (*inv)(row, j) -= factor * (*inv)(col, j);
                }
            }
        }
    }

    return *inv;
}