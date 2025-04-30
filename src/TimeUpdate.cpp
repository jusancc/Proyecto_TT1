#include "../include/TimeUpdate.hpp"

Matrix& TimeUpdate(Matrix &P, Matrix &Phi, double Qdt=0.0){
    Matrix PTrans = transpose(Phi);
    P = Phi*P*PTrans + Qdt;

    return P;
}