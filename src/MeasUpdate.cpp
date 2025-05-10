#include "../include/MeasUpdate.hpp"

tuple<Matrix&,Matrix&,Matrix&> MeasUpdate(Matrix &x, double z, double g, double s, Matrix &G, Matrix &P, int n){
    double Inv_W = s*s;

    Matrix& K = transpose((P*(transpose(G))*(G*P*(transpose(G))+Inv_W).inv()));

    x = x + (K*(z-g));

    P = (eye(n)-transpose(K)*G)*P;

    return tie(K,x,P);
}