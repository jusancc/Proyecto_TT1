#include "..\include\matrix.hpp"
#include <cstdio>
#include <cmath>

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

int m_equals(Matrix A, Matrix B, double p) {
	if (A.n_row != B.n_row || A.n_column != B.n_column)
		return 0;
	else
		for(int i = 1; i <= A.n_row; i++)
			for(int j = 1; j <= A.n_column; j++)
				if(fabs(A(i,j)-B(i,j)) > p) {
					printf("%2.20lf %2.20lf\n",A(i,j),B(i,j));
					return 0;
				}
	
	return 1;
}

int m_sum_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = 2; C(1,2) =  2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = 8; C(2,2) = -3; C(2,3) = 1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = -2; C(3,3) = 0; C(3,4) = 7;
	
	Matrix R = A + B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_sum_02() {
	int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 10; B(1,2) =  12; B(1,3) = 18; B(1,4) = 10;
	B(2,1) = 11; B(2,2) = 9; B(2,3) = 10; B(2,4) = 10;
	B(3,1) = 10; B(3,2) = 11; B(3,3) = 10; B(3,4) = 15;

	Matrix R = A + 10;
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_sub_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) = 0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = -2; C(1,2) = 2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = -6; C(2,2) = 1; C(2,3) = -1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = 4; C(3,3) = 0; C(3,4) = 3;
	
	Matrix R = A - B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_sub_02() {
	int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = -10; B(1,2) =  -8; B(1,3) = -2; B(1,4) = -10;
	B(2,1) = -9; B(2,2) = -11; B(2,3) = -10; B(2,4) = -10;
	B(3,1) = -10; B(3,2) = -9; B(3,3) = -10; B(3,4) = -5;

	Matrix R = A - 10;
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_zeros_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 0; A(1,3) = 0; A(1,4) = 0;
	A(2,1) = 0; A(2,2) = 0; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 0; A(3,4) = 0;
	
	Matrix B = zeros(3, 4);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

// Test para operator()(int n)
int m_access_1d_01() {
    Matrix A(2, 2);
    A(1, 1) = 1; A(1, 2) = 2;
    A(2, 1) = 3; A(2, 2) = 4;

    _assert(fabs(A(1) - 1) < 1e-10);
    _assert(fabs(A(2) - 2) < 1e-10);
    _assert(fabs(A(3) - 3) < 1e-10);
    _assert(fabs(A(4) - 4) < 1e-10);
    return 0;
}

// Test para eye(n)
int m_eye_01() {
    Matrix A = eye(3);
    Matrix B(3, 3);
    B(1, 1) = 1; B(1, 2) = 0; B(1, 3) = 0;
    B(2, 1) = 0; B(2, 2) = 1; B(2, 3) = 0;
    B(3, 1) = 0; B(3, 2) = 0; B(3, 3) = 1;

    _assert(m_equals(A, B, 1e-10));
    return 0;
}

// Test para transpose
int m_transpose_01() {
    Matrix A(2, 3);
    A(1, 1) = 1; A(1, 2) = 2; A(1, 3) = 3;
    A(2, 1) = 4; A(2, 2) = 5; A(2, 3) = 6;

    Matrix B = transpose(A);
    Matrix C(3, 2);
    C(1, 1) = 1; C(1, 2) = 4;
    C(2, 1) = 2; C(2, 2) = 5;
    C(3, 1) = 3; C(3, 2) = 6;

    _assert(m_equals(B, C, 1e-10));
    return 0;
}

// Test para operator+(double s)
int m_scalar_sum_01() {
    Matrix A(2, 2);
    A(1, 1) = 1; A(1, 2) = 2;
    A(2, 1) = 3; A(2, 2) = 4;

    Matrix B = A + 5;
    Matrix C(2, 2);
    C(1, 1) = 6; C(1, 2) = 7;
    C(2, 1) = 8; C(2, 2) = 9;

    _assert(m_equals(B, C, 1e-10));
    return 0;
}

// Test para operator-(double s)
int m_scalar_sub_01() {
    Matrix A(2, 2);
    A(1, 1) = 1; A(1, 2) = 2;
    A(2, 1) = 3; A(2, 2) = 4;

    Matrix B = A - 2;
    Matrix C(2, 2);
    C(1, 1) = -1; C(1, 2) = 0;
    C(2, 1) = 1;  C(2, 2) = 2;

    _assert(m_equals(B, C, 1e-10));
    return 0;
}

// Test para operator*(double s)
int m_scalar_mul_01() {
    Matrix A(2, 2);
    A(1, 1) = 1; A(1, 2) = 2;
    A(2, 1) = 3; A(2, 2) = 4;

    Matrix B = A * 3;
    Matrix C(2, 2);
    C(1, 1) = 3; C(1, 2) = 6;
    C(2, 1) = 9; C(2, 2) = 12;

    _assert(m_equals(B, C, 1e-10));
    return 0;
}

// Test para operator/(double s)
int m_scalar_div_01() {
    Matrix A(2, 2);
    A(1, 1) = 2; A(1, 2) = 4;
    A(2, 1) = 6; A(2, 2) = 8;

    Matrix B = A / 2;
    Matrix C(2, 2);
    C(1, 1) = 1; C(1, 2) = 2;
    C(2, 1) = 3; C(2, 2) = 4;

    _assert(m_equals(B, C, 1e-10));
    return 0;
}

// Test para norm()
int m_norm_01() {
    Matrix A(2, 2);
    A(1, 1) = 1; A(1, 2) = 2;
    A(2, 1) = 3; A(2, 2) = 4;

    double norm = A.norm();
    double expected = sqrt(1*1 + 2*2 + 3*3 + 4*4);  // sqrt(30) ≈ 5.477
    _assert(fabs(norm - expected) < 1e-10);
    return 0;
}

// Test para dot()
int m_dot_01() {
    Matrix A(2, 2);
    A(1, 1) = 1; A(1, 2) = 2;
    A(2, 1) = 3; A(2, 2) = 4;

    Matrix B(2, 2);
    B(1, 1) = 5; B(1, 2) = 6;
    B(2, 1) = 7; B(2, 2) = 8;

    double dot = A.dot(B);
    double expected = 1*5 + 2*6 + 3*7 + 4*8;  // 70
    _assert(fabs(dot - expected) < 1e-10);
    return 0;
}

// Test para v_cross()
int m_v_cross_01() {
    Matrix v(3, 1);
    v(1, 1) = 1; v(2, 1) = 2; v(3, 1) = 3;

    Matrix w(3, 1);
    w(1, 1) = 4; w(2, 1) = 5; w(3, 1) = 6;

    Matrix result = v.v_cross(v, w);
    Matrix expected(3, 1);
    expected(1, 1) = -3;  // 2*6 - 3*5
    expected(2, 1) = 6;   // 3*4 - 1*6
    expected(3, 1) = -3;  // 1*5 - 2*4

    _assert(m_equals(result, expected, 1e-10));
    return 0;
}

// Test para zeros(n)
int m_zeros_square_01() {
    Matrix A = zeros(3);
    Matrix B(3, 3);
    for (int i = 1; i <= 3; i++)
        for (int j = 1; j <= 3; j++)
            B(i, j) = 0;

    _assert(m_equals(A, B, 1e-10));
    return 0;
}

int all_tests()
{
    _verify(m_sum_01);
	_verify(m_sum_02);
    _verify(m_sub_01);
	_verify(m_sub_02);
    _verify(m_zeros_01);
	_verify(m_sum_01);
    _verify(m_sum_02);
    _verify(m_sub_01);
    _verify(m_sub_02);
    _verify(m_zeros_01);
    _verify(m_access_1d_01);
    _verify(m_eye_01);
    _verify(m_transpose_01);
    _verify(m_scalar_sum_01);
    _verify(m_scalar_sub_01);
    _verify(m_scalar_mul_01);
    _verify(m_scalar_div_01);
    _verify(m_norm_01);
    _verify(m_dot_01);
    _verify(m_v_cross_01);
    _verify(m_zeros_square_01);

    return 0;
}


int main()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
