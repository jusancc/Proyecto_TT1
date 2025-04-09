#ifndef _MATRIX_
#define _MATRIX_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

class Matrix {
private:	
	double **data;
public:
    int n_row, n_column;

    // Parameterized constructor
	Matrix(const int v_size);
    Matrix(const int n_row, const int n_column);
	
	// Member operators
	double& operator () (const int n);
	double& operator () (const int row, const int column);
	Matrix& operator + (Matrix &m);
	Matrix& operator + (double s);
	Matrix& operator - (Matrix &m);
	Matrix& operator - (double s);
	Matrix& operator * (Matrix &m);
	Matrix& operator * (double s);
	Matrix& operator / (Matrix &m);
	Matrix& operator / (double s);
	Matrix& operator = (Matrix &m);
	
	double& norm();
	double& dot();
	double& v_cross(Matrix &v, Matrix &w);
	
	// Non-member operators
	friend ostream& operator << (ostream &o, Matrix &m);
};

// Operator overloading
ostream& operator << (ostream &o, Matrix &m);

// Methods
Matrix& zeros(const int n_row, const int n_column);
Matrix& zeros(const int n);
Matrix& eye(const int n);
Matrix& transpose(Matrix &m);

#endif