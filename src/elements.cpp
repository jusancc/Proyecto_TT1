#include "../include/elements.hpp"

/**
 * @brief Calcula los elementos orbitales keplerianos a partir del vector de estado.
 * 
 * @param y Vector de estado (6x1): posición (x,y,z) y velocidad (vx,vy,vz)
 * @return tuple<double, double, double, double, double, double, double> con:
 *         p (semilatus rectum), a (semi-eje mayor), e (excentricidad),
 *         i (inclinación), Omega (longitud nodo ascendente), 
 *         omega (argumento del periastro), M (anomalía media)
 */
tuple<double, double, double, double, double, double, double> elements(Matrix& y) {
    const double pi2 = 2.0 * SAT_Const::pi;

    Matrix r = transpose(y.extract_vector(1, 3)); // posición
    Matrix v = transpose(y.extract_vector(4, 6)); // velocidad

    Matrix h = v_cross(r, v);                        // momento angular específico
    double magh = h.norm();
    double H = magh;

    double p = (magh * magh) / SAT_Const::GM_Earth;

    double Omega = atan2(h(1), -h(2));             // longitud nodo ascendente
    Omega = fmod(Omega + pi2, pi2);

    double i = atan2(sqrt(h(1)*h(1) + h(2)*h(2)), h(3));  // inclinación

    double u = atan2(r(3) * H, -r(1) * h(2) + r(2) * h(1)); // argumento de la latitud

    double R = r.norm();                           // distancia

    double a = 1.0 / (2.0 / R - v.dot(v) / SAT_Const::GM_Earth); // semi-eje mayor

    double eCosE = 1.0 - R / a;
    double eSinE = r.dot(v) / sqrt(SAT_Const::GM_Earth * a);

    double e2 = eCosE * eCosE + eSinE * eSinE;
    double e = sqrt(e2);                           // excentricidad

    double E = atan2(eSinE, eCosE);                // anomalía excéntrica

    double M = fmod(E - eSinE + pi2, pi2);          // anomalía media

    double nu = atan2(sqrt(1.0 - e2) * eSinE, eCosE - e2); // anomalía verdadera

    double omega = fmod(u - nu + pi2, pi2);         // argumento del periastro

    return make_tuple(p, a, e, i, Omega, omega, M);
}
