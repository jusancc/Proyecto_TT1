#include "..\include\matrix.hpp"

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

double& Matrix::operator () (const int n) {
    int total_size = this->n_row * this->n_column;
    if (n <= 0 || n > total_size) {
        cout << "Matrix get: error in index n\n";
        exit(EXIT_FAILURE);
    }
    
    int row = (n - 1) / this->n_column + 1;
    int column = (n - 1) % this->n_column + 1;
    
    return this->data[row - 1][column - 1];
}

double& Matrix::operator () (const int row, const int column) {
	if (row <= 0 || row > this->n_row || column <= 0 || column > this->n_column) {
		cout << "Matrix get: error in row/column\n";
        exit(EXIT_FAILURE);
	}
	
	return this->data[row - 1][column - 1];
}

Matrix& Matrix::operator + (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sum: error in n_row/n_column\n";
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

Matrix& zeros(const int n_row, const int n_column) {
	Matrix *m_aux = new Matrix(n_row, n_column);
	
	for(int i = 1; i <= n_row; i++) {
		for(int j = 1; j <= n_column; j++) {
			(*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
}

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

Matrix& Matrix::operator * (Matrix &m) {
	if (this->n_column != m.n_row) {
		cout << "Matrix sub: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			double sum = 0;
			for(int k = 1; k <= this->n_column; k++) {
                sum += (*this)(i,k) * m(k,j);
            }
            (*m_aux)(i,j) = sum;
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator / (Matrix &m){

}

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

Matrix& transpose(Matrix &m) {
    Matrix *m_aux = new Matrix(m.n_column, m.n_row);
    
    for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++) {
            (*m_aux)(j, i) = m(i, j);
        }
    }
    
    return *m_aux;
}

Matrix& Matrix::operator + (double s){
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
    
    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i, j) = (*this)(i, j) + s;
        }
    }
    
    return *m_aux;
}

Matrix& Matrix::operator - (double s){
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
    
    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i, j) = (*this)(i, j) - s;
        }
    }
    
    return *m_aux;
}

Matrix& Matrix::operator * (double s){
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
    
    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i, j) = (*this)(i, j) * s;
        }
    }
    
    return *m_aux;
}

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

double& Matrix::norm(){
    double r = 0.0;
    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            r += pow((*this)(i, j), 2);
        }
    }
    return sqrt(r);
}

double& Matrix::dot(){
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

double& Matrix::v_cross(Matrix &v, Matrix &w){
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

Matrix& Matrix::union_vector(Matrix &v1, Matrix &v2) {
    bool is_row_vector = (v1.n_row == 1 && v2.n_row == 1);
    bool is_column_vector = (v1.n_column == 1 && v2.n_column == 1);
    
    if (!is_row_vector && !is_column_vector) {
        cout << "union_vector: inputs must be row or column vectors\n";
        exit(EXIT_FAILURE);
    }
    
    if (is_row_vector) {
        if (v1.n_column != v2.n_column) {
            cout << "union_vector: row vectors must have equal columns\n";
            exit(EXIT_FAILURE);
        }
        
        Matrix *result = new Matrix(2, v1.n_column);
        
        for (int j = 1; j <= v1.n_column; j++) {
            (*result)(1, j) = v1(1, j);
        }
        
        for (int j = 1; j <= v2.n_column; j++) {
            (*result)(2, j) = v2(1, j);
        }
        
        return *result;
    } else {
        if (v1.n_row != v2.n_row) {
            cout << "union_vector: column vectors must have equal rows\n";
            exit(EXIT_FAILURE);
        }
        
        Matrix *result = new Matrix(v1.n_row, 2);
        
        for (int i = 1; i <= v1.n_row; i++) {
            (*result)(i, 1) = v1(i, 1);
        }
        
        for (int i = 1; i <= v2.n_row; i++) {
            (*result)(i, 2) = v2(i, 1);
        }
        
        return *result;
    }
}

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