#include "../include/AzElPa.hpp"

tuple<double,double,Matrix&,Matrix&> AzElPa(Matrix &s){
    double pi2 = SAT_Const::pi2;
    double rho = sqrt(s(1,1)*s(1,1)+s(1,2)*s(1,2));

    double Az = atan2(s(1,1),s(1,2));
    if (Az < 0.0)
    {
        Az = Az + pi2;
    }

    double El = atan(s(1,3)/rho);

    Matrix &dAds = zeros(1,3);
    Matrix &dEds = zeros(1,3);

    dAds(1,1) = s(1,2)/(rho*rho);
    dAds(1,2) = -s(1,1)/(rho*rho);
    dAds(1,3) = 0.0;

    dEds(1,1) = -s(1,1)*s(1,3)/rho;
    dEds(1,2) = -s(1,2)*s(1,3)/rho;
    dEds(1,3) = rho;

    dEds = dEds/s.dot(s);

    return tie(Az,El,dAds,dEds);
}