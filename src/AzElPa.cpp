#include "../include/AzElPa.hpp"

tuple<double,double,Matrix&,Matrix&> AzElPa(Matrix &s){
    double pi2 = SAT_Const::pi2;
    double rho = sqrt(s(1,1)*s(1,1)+s(2,1)*s(2,1));

    double Az = atan2(s(1,1),s(2,1));
    if (Az < 0.0)
    {
        Az = Az + pi2;
    }

    double El = atan(s(3,1)/rho);

    Matrix dAds(3,1);
    Matrix dEds(3,1);

    dAds(1,1) = s(2,1)/(rho*rho);
    dAds(2,1) = -s(1,1)/(rho*rho);
    dAds(3,1) = 0.0;

    dEds(1,1) = -s(1,1)*s(3,1)/rho;
    dEds(2,1) = -s(2,1)*s(3,1)/rho;
    dEds(3,1) = rho;

    dEds = dEds/s.dot(s);

    return tie(Az,El,dAds,dEds);
}