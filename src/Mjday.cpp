//$Header$
//------------------------------------------------------------------------------
//                                   Mjday
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file Mjday.cpp
 * @brief Implementación de la función que calcula el Modified Julian Date (MJD)
 *        a partir de una fecha y hora del calendario gregoriano.
 *
 * El MJD se calcula a partir de la fecha y hora en formato calendario gregoriano.
 * La fórmula utilizada transforma esta fecha en el Julian Date (JD) y después 
 * aplica el desplazamiento estándar para obtener el MJD:
 * 
 *     MJD = JD - 2400000.5
 * 
 * Esta función es fundamental en aplicaciones astronómicas y de dinámica orbital.
 */
//------------------------------------------------------------------------------

#include "../include/Mjday.hpp"
#include <math.h>

/**
 * @brief Calcula el Modified Julian Date (MJD) a partir de una fecha y hora.
 *
 * Utiliza una fórmula estándar basada en el algoritmo de Fliegel y Van Flandern
 * para calcular el Julian Date, al que luego se le resta 2400000.5 para obtener
 * el MJD. Es compatible con fechas del calendario gregoriano.
 *
 * @param yr  Año (ejemplo: 2025)
 * @param mon Mes del año (1 a 12)
 * @param day Día del mes
 * @param hr  Hora del día (0 a 23)
 * @param min Minutos (0 a 59)
 * @param sec Segundos (0.0 a 59.999...)
 * @return Modified Julian Date (días desde el 17 de noviembre de 1858 a las 00:00 UTC)
 */
double Mjday(int yr, int mon, int day, int hr, int min, double sec)
{
    double jd, Mjd;
	
    jd = 367.0 * yr
        - floor((7 * (yr + floor((mon + 9) / 12.0))) * 0.25)
        + floor(275 * mon / 9.0)
        + day + 1721013.5
        + (((sec / 60.0) + min) / 60.0 + hr) / 24.0;

    Mjd = jd - 2400000.5;

    return Mjd;
}
