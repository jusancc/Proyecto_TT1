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
 * @file Mjday.hpp
 * @brief Declaración de la función que calcula el Modified Julian Date (MJD)
 *        a partir de una fecha y hora del calendario gregoriano.
 *
 * La fecha se proporciona en formato año, mes, día, hora, minuto y segundo,
 * y se transforma en el número correspondiente de días julianos modificados,
 * con el día cero en 1858-11-17 a las 00:00:00 UTC.
 */
//------------------------------------------------------------------------------

#ifndef _MJDAY_
#define _MJDAY_

/**
 * @brief Calcula el Modified Julian Date (MJD) a partir de una fecha y hora.
 * 
 * @param yr  Año (ej. 2025)
 * @param mon Mes (1–12)
 * @param day Día del mes
 * @param hr  Hora (opcional, por defecto 0)
 * @param min Minutos (opcional, por defecto 0)
 * @param sec Segundos (opcional, por defecto 0.0)
 * @return Número de días julianos modificados (MJD)
 */
double Mjday(int yr, int mon, int day, int hr = 0, int min = 0, double sec = 0.0);

#endif
