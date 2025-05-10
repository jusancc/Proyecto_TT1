#include "../include/AccelHarmonic.hpp"

Matrix& AccelHarmonic(Matrix &r, Matrix &E, double n_max, double m_max){
    double r_ref = 6378.1363e3;
    double gm = 398600.4415e9;

    Matrix &r_bf = E * r;

    double d = r_bf.norm();
    double latgc = asin(r_bf(1,3)/d);
    double lon = atan2(r_bf(1,2),r_bf(1,1));

    auto [pnm, dpnm] = Legendre(n_max, m_max, latgc);

    int dUdr = 0;
    int dUdlatgc = 0;
    int dUdlon = 0;
    int q3 = 0;
    int q2 = q3;
    int q1 = q2;

    double b1,b2,b3;
    for (int i = 0; i < n_max; i++)
    {
        b1 = (-gm/pow(d,2))*(pow((r_ref/d),(i*(i+1))));
        b2 = (gm/d)*(pow((r_ref/d),i));
        b3 = (gm/d)*(pow((r_ref/d),i));

        for (int j = 0; j < m_max; j++)
        {
            q1 += pnm(i+1,j+1)*(Cnm(i+1,j+1)*cos(j*lon)+Snm(i+1,j+1)*sin(j*lon));
            q2 += dpnm(i+1,j+1)*(Cnm(i+1,j+1)*cos(j*lon)+Snm(i+1,j+1)*sin(j*lon));
            q3 += j*pnm(i+1,j+1)*(Snm(i+1,j+1)*cos(j*lon)-Cnm(i+1,j+1)*sin(j*lon));
        }
        dUdr += q1*b1;
        dUdlatgc += q2*b2;
        dUdlon += q3*b3;
        q3 = 0;
        q2 = q3;
        q1 = q2;
    }
    
    double r2xy = pow(r_bf(1,1),2)+pow(r_bf(1,2),2);

    double ax = (1.0/d*dUdr-r_bf(1,3)/(pow(d,2)*sqrt(r2xy))*dUdlatgc)*r_bf(1,1)-(1.0/r2xy*dUdlon)*r_bf(1,2);
    double ay = (1.0/d*dUdr-r_bf(1,3)/(pow(d,2)*sqrt(r2xy))*dUdlatgc)*r_bf(1,2)+(1.0/r2xy*dUdlon)*r_bf(1,1);
    double az =  1.0/d*dUdr*r_bf(1,3)+sqrt(r2xy)/pow(d,2)*dUdlatgc;

    Matrix& aux = zeros(3);
    aux(1,1)=ax;
    aux(1,2)=ay;
    aux(1,3)=az;
    Matrix& a_bf = transpose(aux);

    Matrix& a = (transpose(E))*a_bf;
    return a;
}