#include "..\include\matrix.hpp"
#include "..\include\AccelPointMass.hpp"
#include "..\include\Cheb3D.hpp"
#include "..\include\EccAnom.hpp"
#include "..\include\Frac.hpp"
#include "..\include\MeanObliquity.hpp"
#include "..\include\Mjday.hpp"
#include "..\include\Mjday_TDB.hpp"
#include "..\include\Position.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include "..\include\Sign.hpp"
#include "..\include\Timediff.hpp"
#include "..\include\AzElPa.hpp"
#include "..\include\IERS.hpp"
#include "..\include\Legendre.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\TimeUpdate.hpp"
#include "..\include\EqnEquinox.hpp"
#include "..\include\LTC.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\PrecMatrix.hpp"
#include "..\include\GMST.hpp"
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
    double expected = 1*5 + 2*6 + 3*7 + 4*8;
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

int m_inverse_02() {
    Matrix C(3, 3);
    C(1,1)=1; C(1,2)=2; C(1,3)=3;
    C(2,1)=0; C(2,2)=1; C(2,3)=4;
    C(3,1)=5; C(3,2)=6; C(3,3)=0;


    Matrix invC = C.inv();

    Matrix I(3,3);
    I(1,1)=-24; I(1,2)=18; I(1,3)=5;
    I(2,1)=20; I(2,2)=-15; I(2,3)=-4;
    I(3,1)=-5; I(3,2)=4; I(3,3)=1;
    _assert(m_equals(invC, I, 1e-10));
    return 0;
}

int m_accel_point_mass_01(){
    Matrix A(1,3);
    A(1,1)=1;A(1,2)=2;A(1,3)=3;
    Matrix B(1,3);
    B(1,1)=4;B(1,2)=5;B(1,3)=6;
    int GM = 5;

    Matrix C = AccelPointMass(A,B,GM);

    Matrix result(1,3);
    result(1,1)=0.0773165667868213;result(1,2)=0.0699165293543773;result(1,3)=0.0625164919219332;

    _assert(m_equals(result,C,1e-10));
    return 0;
}

int m_cheb_3d_01(){
    Matrix Cx(1,3);
    Cx(1,1)=1;Cx(1,2)=2;Cx(1,3)=3;
    Matrix Cy(1,3);
    Cy(1,1)=4;Cy(1,2)=5;Cy(1,3)=6;
    Matrix Cz(1,3);
    Cz(1,1)=7;Cz(1,2)=8;Cz(1,3)=9;

    Matrix A = Cheb3D(5,3,0,10,Cx,Cy,Cz);
    cout << A << endl;
    Matrix result(1,3);
    result(1,1)=-2;result(1,2)=-2;result(1,3)=-2;

    _assert(m_equals(A,result,1e-10));
    return 0;
}

int m_ecc_anom_01(){
    _assert(m_equals(2.9504, EccAnom(2,5),1e-10));
    return 0;
}

int m_frac_01(){
    _assert(fabs(0.7-Frac(10.7))<1e-10);
    return 0;
}

int m_mean_obliquity_(){
    _assert(fabs(0.4090928042-MeanObliquity(51544.5))<1e-10);
    return 0;
}

int m_mjday_01(){
    _assert(fabs(51544.0 - Mjday(2000, 1, 1, 0, 0, 0)) < 1e-10);
    return 0;
}

int m_mjday_tdb_01(){
    _assert(fabs(2000.00000001537-Mjday_TDB(2000))<1e-10);
    return 0;
}

int m_position_01(){
    Matrix C(1,3);
    C(1,1)=6378136.3;C(1,2)=0.0;C(1,3)=0.0;

    Matrix res(3);
    res = Position(0.0,0.0,0.0);

    _assert(m_equals(C,res,1e-10));
    return 0;
}

int m_rx_01(){
    Matrix A(3,3);
    A(1,1)=1;A(1,2)=0;A(1,3)=0;
    A(2,1)=0;A(2,2)=-0.416146836547142;A(2,3)=0.909297426825682;
    A(3,1)=0;A(3,2)=-0.909297426825682;A(3,3)=-0.416146836547142;

    Matrix res(3,3);
    res = R_x(2);

    _assert(m_equals(A,res,1e-10));
    return 0;
}

int m_ry_01(){
    Matrix A(3,3);
    A(1,1)=-0.416146836547142 ;A(1,2)=0;A(1,3)=-0.909297426825682;
    A(2,1)=0;A(2,2)=1;A(2,3)=0;
    A(3,1)=0.909297426825682;A(3,2)=0;A(3,3)=-0.416146836547142;

    Matrix res(3,3);
    res = R_y(2);

    _assert(m_equals(A,res,1e-10));
    return 0;
}

int m_rz_01(){
    Matrix A(3,3);
    A(1,1)=-0.416146836547142 ;A(1,2)=0.909297426825682;A(1,3)=0;
    A(2,1)=-0.909297426825682;A(2,2)=-0.416146836547142;A(2,3)=0;
    A(3,1)=0;A(3,2)=0;A(3,3)=1;

    Matrix res(3,3);
    res = R_z(2);

    _assert(m_equals(A,res,1e-10));
    return 0; 
}

int m_sign_01(){
    _assert(fabs(-2-sign_(2,-3))<1e-10);
    return 0;
}

int m_timediff_01(){

    auto [UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(1, 2);

    double ans_UT1_TAI = -1;
    double ans_UTC_GPS = 17;
    double ans_UT1_GPS = 18;
    double ans_TT_UTC = 34.184;
    double ans_GPS_UTC = -17;

    _assert(fabs(UT1_TAI - ans_UT1_TAI)< 1e-10);
    _assert(fabs(UTC_GPS - ans_UTC_GPS)< 1e-10);
    _assert(fabs(UT1_GPS - ans_UT1_GPS)< 1e-10);
    _assert(fabs(TT_UTC - ans_TT_UTC)< 1e-10);
    _assert(fabs(GPS_UTC - ans_GPS_UTC)< 1e-10);
	
    return 0;
}

int m_azElPa_01(){
    Matrix A(1,3);
    A(1,1)=1;A(1,2)=2;A(1,3)=3;

    auto [Az, El, dAds, dEds] = AzElPa(A);
    double ans_Az = 0.463647609000806;
    double ans_El = 0.930274014115472;
    Matrix ans_dAds(1,3);
    ans_dAds(1,1)=0.4;ans_dAds(1,2)=-0.2;ans_dAds(1,3)=0;
    Matrix ans_dEds(1,3);
    ans_dEds(1,1)=-0.095831484749991;ans_dEds(1,2)=-0.191662969499982;ans_dEds(1,3)=0.159719141249985;

    _assert(fabs(ans_Az-Az)<1e-10);
    _assert(fabs(ans_El-El)<1e-10);
    _assert(m_equals(ans_dAds,dAds,1e-10));
    _assert(m_equals(ans_dEds,dEds,1e-10));

    return 0;
}

int m_iers_01(){
    Matrix eop(13,2);

    eop(1,1)=0;         eop(1,2)=0;
    eop(2,1)=0;         eop(2,2)=0;
    eop(3,1)=0;         eop(3,2)=0;
    eop(4,1)=58000;     eop(4,2)=58001;
    eop(5,1)=0.055;     eop(5,2)=0.056;
    eop(6,1)=0.325;     eop(6,2)=0.326;
    eop(7,1)=-0.1;      eop(7,2)=-0.09;
    eop(8,1)=0.00015;   eop(8,2)=0.00014;
    eop(9,1)=-0.054;    eop(9,2)=-0.053;
    eop(10,1)=0.004;    eop(10,2)=0.005;
    eop(11,1)=0.003;    eop(11,2)=0.0031;
    eop(12,1)=-0.001;   eop(12,2)=-0.0011;
    eop(13,1)=37;       eop(13,2)=37;
    
    double Mjd_UTC = 58000.5;
    char interp = 'l';

    auto [x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC] = IERS(eop, Mjd_UTC, interp);

    double expected_x_pole = 2.69071593015792e-07;
    double expected_y_pole = 1.57806853201154e-06;
    double expected_UT1_UTC = -0.095;
    double expected_LOD = 0.000145;
    double expected_dpsi = -2.59375319393602e-07;
    double expected_deps = 2.18166156499291e-08;
    double expected_dx_pole = 1.47868172738408e-08;
    double expected_dy_pole = -5.09054365165013e-09;
    double expected_TAI_UTC = 37;

    _assert(fabs(expected_x_pole - x_pole)< 1e-10);
    _assert(fabs(expected_y_pole - y_pole)< 1e-10);
    _assert(fabs(expected_UT1_UTC - UT1_UTC)< 1e-10);
    _assert(fabs(expected_LOD - LOD)< 1e-10);
    _assert(fabs(expected_dpsi - dpsi)< 1e-10);
    _assert(fabs(expected_deps - deps)< 1e-10);
    _assert(fabs(expected_dx_pole - dx_pole)< 1e-10);
    _assert(fabs(expected_dy_pole - dy_pole)< 1e-10);
    _assert(fabs(expected_TAI_UTC - TAI_UTC)< 1e-10);
    return 0;
}

int m_legendre_01(){

    auto [pnm, dpnm] = Legendre(1,2,3);

    Matrix ans_pnm(2,3);
    ans_pnm(1,1)=1;ans_pnm(1,2)=0;ans_pnm(1,3)=0;
    ans_pnm(2,1)=0.244427023924219;ans_pnm(2,2)=-1.71471730322393;ans_pnm(2,3)=0;

    Matrix ans_dpnm(2,3);
    ans_dpnm(1,1)=0;ans_dpnm(1,2)=0;ans_dpnm(1,3)=0;
    ans_dpnm(2,1)=-1.71471730322393;ans_dpnm(2,2)=-0.244427023924219;ans_dpnm(2,3)=0;

    _assert(m_equals(ans_pnm,pnm,1e-9));
    _assert(m_equals(ans_dpnm,dpnm,1e-9));
    return 0;
}

int m_nutangles_01(){
    auto [dpsi, deps] = NutAngles(2);

    double ans_dpsi = 2.7179807523643e-05;
    double ans_deps = 3.91872311875582e-05;

    _assert(fabs(ans_dpsi-dpsi)<1e-10);
    _assert(fabs(ans_deps-deps)<1e-10);

    return 0;
}

int m_timeupdate_01(){
    Matrix P(2,2);
    P(1,1) = 1; P(1,2) = 2;
    P(2,1) = 3; P(2,2) = 4;

    Matrix Phi(2,2);
    Phi(1,1) = 5; Phi(1,2) = 6;
    Phi(2,1) = 7; Phi(2,2) = 8;

    double Qdt = 0.5;
    
    Matrix R = TimeUpdate(P, Phi, Qdt);

    Matrix expected(2,2);
    expected(1,1) = 319.5; expected(1,2) = 433.5;
    expected(2,1) = 431.5; expected(2,2) = 585.5;


    _assert(m_equals(expected, R, 1e-10));
    return 0;
}

int m_eqnEquinox_01(){
    _assert(fabs(2.49335221515174e-05 - EqnEquinox(2))<1e-10);
    return 0;
}

int m_ltc_01(){
    Matrix A(3,3);
    A(1,1)=-0.909297426825682;A(1,2)=-0.416146836547142;A(1,3)=0;
    A(2,1)=0.058726644927621;A(2,2)=-0.128320060202457;A(2,3)=-0.989992496600445;
    A(3,1)=0.411982245665683;A(3,2)=-0.900197629735517;A(3,3)=0.141120008059867;

    Matrix res(3,3);
    res = LTC(2,3);

    _assert(m_equals(A,res,1e-10));
    return 0;
}

int m_nutMatrix_01(){
    Matrix A(3,3);
    A(1,1)=0.999999999630629;A(1,2)=-2.49335221484475e-05;A(1,3)=-1.08194921374918e-05;
    A(2,1)=2.49330981433634e-05;A(2,2)=0.999999998921346;A(2,3)=-3.91873660592346e-05;
    A(3,1)=1.08204692048809e-05;A(3,2)=3.91870962813123e-05;A(3,3)=0.999999999173645;

    Matrix res(3,3);
    res = NutMatrix(2);

    _assert(m_equals(A,res,1e-10));
    return 0;
}

int m_poleMatrix_01(){
    Matrix A(3,3);
    A(1,1)=-0.416146836547142;A(1,2)=0.128320060202457;A(1,3)=-0.900197629735517;
    A(2,1)=0;A(2,2)=-0.989992496600445;A(2,3)=-0.141120008059867;
    A(3,1)=-0.909297426825682;A(3,2)=-0.058726644927621;A(3,3)=0.411982245665683;

    Matrix res(3,3);
    res = PoleMatrix(2,3);

    _assert(m_equals(A,res,1e-10));
    return 0;
}

int m_precMatrix_01(){
    Matrix A(3,3);
    A(1,1)=0.999999999999778;A(1,2)=-6.11707327974946e-07;A(1,3)=-2.66201482252295e-07;
    A(2,1)=6.11707327974946e-07;A(2,2)=0.999999999999813;A(2,3)=-8.14186990889656e-14;
    A(3,1)=2.66201482252295e-07;A(3,2)=-8.1418698322574e-14;A(3,3)=0.999999999999965;

    Matrix res(3,3);
    res = PrecMatrix(2,3);

    _assert(m_equals(A,res,1e-10));
    return 0;
}

int m_gmst_01(){
    _assert(fabs(1.00761373073367-gmst(2))<1e-10);
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
    _verify(m_inverse_02);

    
    _verify(m_accel_point_mass_01);
    _verify(m_cheb_3d_01);
    _verify(m_ecc_anom_01);
    _verify(m_frac_01);
    _verify(m_mjday_01);
    _verify(m_mjday_tdb_01);
    _verify(m_position_01);
    _verify(m_rx_01);
    _verify(m_ry_01);
    _verify(m_rz_01);
    _verify(m_sign_01);
    _verify(m_timediff_01);
    _verify(m_azElPa_01);
    _verify(m_iers_01);
    _verify(m_legendre_01);
    _verify(m_nutangles_01);
    _verify(m_timeupdate_01);
    _verify(m_eqnEquinox_01);
    _verify(m_ltc_01);
    _verify(m_nutMatrix_01);
    _verify(m_poleMatrix_01);
    _verify(m_precMatrix_01);
    _verify(m_gmst_01);
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
