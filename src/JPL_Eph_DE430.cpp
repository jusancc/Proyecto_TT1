//$Header$
//------------------------------------------------------------------------------
//                              JPL_Eph_DE430
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file JPL_Eph_DE430.cpp
 * @brief Interpola las posiciones de los cuerpos del sistema solar usando la efeméride DE430.
 *
 * Utiliza coeficientes de Chebyshev almacenados en la tabla `PC` para calcular las posiciones
 * de los planetas, la Luna y el Sol en un instante dado expresado en tiempo TDB.
 * Las posiciones se devuelven en kilómetros en un sistema de referencia centrado en la Tierra.
 */
//------------------------------------------------------------------------------

#include "../include/JPL_Eph_DE430.hpp"

/**
 * @brief Genera un vector de índices para extracción secuencial.
 * @param a Primer índice
 * @param b Paso
 * @param c Último índice
 * @return Vector con los índices generados
 */
Matrix &_temp(int a, int b, int c)
{
    int len = (c - a) / b + 1;
    Matrix &temp = zeros(len);
    for (int i = 1; i <= len; i++)
    {
        temp(i) = (i - 1) * b + a;
    }
    return temp;
}

/**
 * @brief Evalúa las posiciones de los planetas, Sol y Luna usando DE430 para un tiempo TDB.
 *
 * @param Mjd_TDB Tiempo dinámico baricéntrico (MJD).
 * @return Tupla con las posiciones (en km) de:
 *         - Mercurio
 *         - Venus
 *         - Tierra
 *         - Marte
 *         - Júpiter
 *         - Saturno
 *         - Urano
 *         - Neptuno
 *         - Plutón
 *         - Luna
 *         - Sol
 */
tuple<Matrix &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &> JPL_Eph_DE430(double Mjd_TDB)
{
    double JD = Mjd_TDB + 2400000.5;

    int i = 1;
    bool found = false;
    while (i <= PC.n_row && !found) {
        if (PC(i, 1) <= JD && JD <= PC(i, 2)) {
            found = true;
        } else {
            i++;
        }
    }
    Matrix& PCtemp = PC.extract_row(i);

    double t1 = PCtemp(1, 1) - 2400000.5;
    double dt = Mjd_TDB - t1;

    // Earth
    Matrix& temp = _temp(231, 13, 270);
    Matrix& Cx_Earth = PCtemp.extract_vector(temp(1,1), temp(1,2)-1);
    Matrix& Cy_Earth = PCtemp.extract_vector(temp(1,2), temp(1,3)-1);
    Matrix& Cz_Earth = PCtemp.extract_vector(temp(1,3), temp(1,4)-1);

    temp = temp + 39;
    Matrix& Cx = PCtemp.extract_vector(temp(1,1), temp(1,2)-1);
    Matrix& Cy = PCtemp.extract_vector(temp(1,2), temp(1,3)-1);
    Matrix& Cz = PCtemp.extract_vector(temp(1,3), temp(1,4)-1);

    Cx_Earth = union_vector(Cx_Earth, Cx);
    Cy_Earth = union_vector(Cy_Earth, Cy);
    Cz_Earth = union_vector(Cz_Earth, Cz);

    int j;
    double Mjd0;
    if (0 <= dt && dt <= 16) {
        j = 0;
        Mjd0 = t1;
    } else {
        j = 1;
        Mjd0 = t1 + 16 * j;
    }

    Matrix &r_Earth = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+16, 
                             Cx_Earth.extract_vector(13*j+1, 13*j+13),
                             Cy_Earth.extract_vector(13*j+1, 13*j+13),
                             Cz_Earth.extract_vector(13*j+1, 13*j+13)) * 1e3;

    // Moon
    temp = _temp(441, 13, 480);
    Matrix& Cx_Moon = PCtemp.extract_vector(temp(1,1), temp(1,2)-1);
    Matrix& Cy_Moon = PCtemp.extract_vector(temp(1,2), temp(1,3)-1);
    Matrix& Cz_Moon = PCtemp.extract_vector(temp(1,3), temp(1,4)-1);
    for (int i = 1; i <= 7; i++) {
        temp = temp + 39;
        Cx = PCtemp.extract_vector(temp(1,1), temp(1,2)-1);
        Cy = PCtemp.extract_vector(temp(1,2), temp(1,3)-1);
        Cz = PCtemp.extract_vector(temp(1,3), temp(1,4)-1);

        Cx_Moon = union_vector(Cx_Moon, Cx);
        Cy_Moon = union_vector(Cy_Moon, Cy);
        Cz_Moon = union_vector(Cz_Moon, Cz);
    }
    if (0 <= dt && dt <= 4) {
        j = 0;
        Mjd0 = t1;
    } else if (4 < dt && dt <= 8) {
        j = 1;
        Mjd0 = t1 + 4 * j;
    } else if (8 < dt && dt <= 12) {
        j = 2;
        Mjd0 = t1 + 4 * j;
    } else if (12 < dt && dt <= 16) {
        j = 3;
        Mjd0 = t1 + 4 * j;
    } else if (16 < dt && dt <= 20) {
        j = 4;
        Mjd0 = t1 + 4 * j;
    } else if (20 < dt && dt <= 24) {
        j = 5;
        Mjd0 = t1 + 4 * j;
    } else if (24 < dt && dt <= 28) {
        j = 6;
        Mjd0 = t1 + 4 * j;
    } else {
        j = 7;
        Mjd0 = t1 + 4 * j;
    }

    Matrix &r_Moon = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+4,
                            Cx_Moon.extract_vector(13*j+1, 13*j+13),
                            Cy_Moon.extract_vector(13*j+1, 13*j+13),
                            Cz_Moon.extract_vector(13*j+1, 13*j+13)) * 1e3;

    // Sun
    temp = _temp(753, 11, 786);
    Matrix& Cx_Sun = PCtemp.extract_vector(temp(1,1), temp(1,2)-1);
    Matrix& Cy_Sun = PCtemp.extract_vector(temp(1,2), temp(1,3)-1);
    Matrix& Cz_Sun = PCtemp.extract_vector(temp(1,3), temp(1,4)-1);
    temp = temp + 33;
    Cx = PCtemp.extract_vector(temp(1,1), temp(1,2)-1);
    Cy = PCtemp.extract_vector(temp(1,2), temp(1,3)-1);
    Cz = PCtemp.extract_vector(temp(1,3), temp(1,4)-1);
    Cx_Sun = union_vector(Cx_Sun, Cx);
    Cy_Sun = union_vector(Cy_Sun, Cy);
    Cz_Sun = union_vector(Cz_Sun, Cz);

    if (0 <= dt && dt <= 16) {
        j = 0;
        Mjd0 = t1;
    } else {
        j = 1;
        Mjd0 = t1 + 16 * j;
    }

    Matrix &r_Sun = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+16,
                           Cx_Sun.extract_vector(11*j+1, 11*j+11),
                           Cy_Sun.extract_vector(11*j+1, 11*j+11),
                           Cz_Sun.extract_vector(11*j+1, 11*j+11)) * 1e3;

    // Mercury
    temp = _temp(3, 14, 45);
    Matrix& Cx_Mercury = PCtemp.extract_vector(temp(1,1), temp(1,2)-1);
    Matrix& Cy_Mercury = PCtemp.extract_vector(temp(1,2), temp(1,3)-1);
    Matrix& Cz_Mercury = PCtemp.extract_vector(temp(1,3), temp(1,4)-1);
    for (int i = 1; i <= 3; i++) {
        temp = temp + 42;
        Cx = PCtemp.extract_vector(temp(1,1), temp(1,2)-1);
        Cy = PCtemp.extract_vector(temp(1,2), temp(1,3)-1);
        Cz = PCtemp.extract_vector(temp(1,3), temp(1,4)-1);
        Cx_Mercury = union_vector(Cx_Mercury, Cx);
        Cy_Mercury = union_vector(Cy_Mercury, Cy);
        Cz_Mercury = union_vector(Cz_Mercury, Cz);
    }

    if (0 <= dt && dt <= 8) {
        j = 0;
        Mjd0 = t1;
    } else if (8 < dt && dt <= 16) {
        j = 1;
        Mjd0 = t1 + 8 * j;
    } else if (16 < dt && dt <= 24) {
        j = 2;
        Mjd0 = t1 + 8 * j;
    } else {
        j = 3;
        Mjd0 = t1 + 8 * j;
    }

    Matrix &r_Mercury = Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0+8,
                                Cx_Mercury.extract_vector(14*j+1, 14*j+14),
                                Cy_Mercury.extract_vector(14*j+1, 14*j+14),
                                Cz_Mercury.extract_vector(14*j+1, 14*j+14)) * 1e3;

    // Venus
    temp = _temp(171, 10, 201);
    Matrix& Cx_Venus = PCtemp.extract_vector(temp(1,1), temp(1,2)-1);
    Matrix& Cy_Venus = PCtemp.extract_vector(temp(1,2), temp(1,3)-1);
    Matrix& Cz_Venus = PCtemp.extract_vector(temp(1,3), temp(1,4)-1);
    temp = temp + 30;
    Cx = PCtemp.extract_vector(temp(1,1), temp(1,2)-1);
    Cy = PCtemp.extract_vector(temp(1,2), temp(1,3)-1);
    Cz = PCtemp.extract_vector(temp(1,3), temp(1,4)-1);
    Cx_Venus = union_vector(Cx_Venus, Cx);
    Cy_Venus = union_vector(Cy_Venus, Cy);
    Cz_Venus = union_vector(Cz_Venus, Cz);

    if (0 <= dt && dt <= 16) {
        j = 0;
        Mjd0 = t1;
    } else {
        j = 1;
        Mjd0 = t1 + 16 * j;
    }

    Matrix &r_Venus = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+16,
                              Cx_Venus.extract_vector(10*j+1, 10*j+10),
                              Cy_Venus.extract_vector(10*j+1, 10*j+10),
                              Cz_Venus.extract_vector(10*j+1, 10*j+10)) * 1e3;

    // Mars
    temp = _temp(309, 11, 342);
    Matrix& Cx_Mars = PCtemp.extract_vector(temp(1,1), temp(1,2)-1);
    Matrix& Cy_Mars = PCtemp.extract_vector(temp(1,2), temp(1,3)-1);
    Matrix& Cz_Mars = PCtemp.extract_vector(temp(1,3), temp(1,4)-1);
    j = 0;
    Mjd0 = t1;
    Matrix &r_Mars = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+32,
                            Cx_Mars.extract_vector(11*j+1, 11*j+11),
                            Cy_Mars.extract_vector(11*j+1, 11*j+11),
                            Cz_Mars.extract_vector(11*j+1, 11*j+11)) * 1e3;

    // Jupiter
    temp = _temp(342, 8, 366);
    Matrix& Cx_Jupiter = PCtemp.extract_vector(temp(1,1), temp(1,2)-1);
    Matrix& Cy_Jupiter = PCtemp.extract_vector(temp(1,2), temp(1,3)-1);
    Matrix& Cz_Jupiter = PCtemp.extract_vector(temp(1,3), temp(1,4)-1);
    j = 0;
    Mjd0 = t1;
    Matrix &r_Jupiter = Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0+32,
                               Cx_Jupiter.extract_vector(8*j+1, 8*j+8),
                               Cy_Jupiter.extract_vector(8*j+1, 8*j+8),
                               Cz_Jupiter.extract_vector(8*j+1, 8*j+8)) * 1e3;

    // Saturn
    temp = _temp(366, 7, 387);
    Matrix& Cx_Saturn = PCtemp.extract_vector(temp(1,1), temp(1,2)-1);
    Matrix& Cy_Saturn = PCtemp.extract_vector(temp(1,2), temp(1,3)-1);
    Matrix& Cz_Saturn = PCtemp.extract_vector(temp(1,3), temp(1,4)-1);
    j = 0;
    Mjd0 = t1;
    Matrix &r_Saturn = Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0+32,
                              Cx_Saturn.extract_vector(7*j+1, 7*j+7),
                              Cy_Saturn.extract_vector(7*j+1, 7*j+7),
                              Cz_Saturn.extract_vector(7*j+1, 7*j+7)) * 1e3;

    // Uranus
    temp = _temp(387, 6, 405);
    Matrix& Cx_Uranus = PCtemp.extract_vector(temp(1,1), temp(1,2)-1);
    Matrix& Cy_Uranus = PCtemp.extract_vector(temp(1,2), temp(1,3)-1);
    Matrix& Cz_Uranus = PCtemp.extract_vector(temp(1,3), temp(1,4)-1);
    j = 0;
    Mjd0 = t1;
    Matrix &r_Uranus = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32,
                              Cx_Uranus.extract_vector(6*j+1, 6*j+6),
                              Cy_Uranus.extract_vector(6*j+1, 6*j+6),
                              Cz_Uranus.extract_vector(6*j+1, 6*j+6)) * 1e3;

    // Neptune
    temp = _temp(405, 6, 423);
    Matrix& Cx_Neptune = PCtemp.extract_vector(temp(1,1), temp(1,2)-1);
    Matrix& Cy_Neptune = PCtemp.extract_vector(temp(1,2), temp(1,3)-1);
    Matrix& Cz_Neptune = PCtemp.extract_vector(temp(1,3), temp(1,4)-1);
    j = 0;
    Mjd0 = t1;
    Matrix &r_Neptune = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32,
                                Cx_Neptune.extract_vector(6*j+1, 6*j+6),
                                Cy_Neptune.extract_vector(6*j+1, 6*j+6),
                                Cz_Neptune.extract_vector(6*j+1, 6*j+6)) * 1e3;

    // Pluto
    temp = _temp(423, 6, 441);
    Matrix& Cx_Pluto = PCtemp.extract_vector(temp(1,1), temp(1,2)-1);
    Matrix& Cy_Pluto = PCtemp.extract_vector(temp(1,2), temp(1,3)-1);
    Matrix& Cz_Pluto = PCtemp.extract_vector(temp(1,3), temp(1,4)-1);
    j = 0;
    Mjd0 = t1;
    Matrix &r_Pluto = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32,
                             Cx_Pluto.extract_vector(6*j+1, 6*j+6),
                             Cy_Pluto.extract_vector(6*j+1, 6*j+6),
                             Cz_Pluto.extract_vector(6*j+1, 6*j+6)) * 1e3;

    // Ajuste posiciones respecto a Tierra
    const double EMRAT = 81.30056907419062;
    const double EMRAT1 = 1.0 / (1.0 + EMRAT);

    for (int k = 1; k <= 3; k++) {
        r_Earth(k) -= EMRAT1 * r_Moon(k);
        r_Mercury(k) = -r_Earth(k) + r_Mercury(k);
        r_Venus(k) = -r_Earth(k) + r_Venus(k);
        r_Mars(k) = -r_Earth(k) + r_Mars(k);
        r_Jupiter(k) = -r_Earth(k) + r_Jupiter(k);
        r_Saturn(k) = -r_Earth(k) + r_Saturn(k);
        r_Uranus(k) = -r_Earth(k) + r_Uranus(k);
        r_Neptune(k) = -r_Earth(k) + r_Neptune(k);
        r_Pluto(k) = -r_Earth(k) + r_Pluto(k);
        r_Sun(k) = -r_Earth(k) + r_Sun(k);
    }

    return tie(r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun);
}
