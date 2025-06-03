#include "../include/angl.hpp"

double angl(Matrix& vec1, Matrix& vec2){
    double small = 0.00000001;
    double undefined = 999999.1;

    double magv1 = vec1.norm();
    double magv2 = vec2.norm();

    double theta;
    if (magv1*magv2 > pow(small,2))
    {
        double temp = vec1.dot(vec2) / (magv1*magv2);
        if (fabs(temp) > 1.0){
            if (temp > 0)
                temp = 1.0;
            else if (temp < 0)
                temp = -1.0;
            else
                temp = 0.0;
        }
        theta = acos(temp);
    } else{
        theta = undefined;
    }

    return theta;
    
}