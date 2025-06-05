//$Header$
//------------------------------------------------------------------------------
//                                DEInteg
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// **Legal**
//
// Author: Juan Sánchez de Corta
//
/**
 * @file DEInteg.cpp
 * @brief Implementación del integrador multistep de Shampine y Gordon.
 *
 * Este integrador numérico resuelve sistemas de EDOs usando pasos adaptativos
 * y una estrategia de control de errores basada en estimaciones embebidas.
 * La implementación gestiona casos de rigidez, control de paso y condiciones
 * de error relativo/absoluto.
 */
//------------------------------------------------------------------------------

#include "../include/DEInteg.hpp"

//------------------------------------------------------------------------------
//  Matrix& DEInteg(Matrix &func(double, Matrix &), double &t, double tout,
//                  double &relerr, double &abserr, int n_eqn, Matrix &y)
//------------------------------------------------------------------------------
/**
 * @brief Resuelve sistemas de EDOs utilizando un método multistep adaptativo
 *        de Shampine y Gordon.
 *
 * Esta función resuelve un sistema de ecuaciones diferenciales ordinarias
 * de la forma y' = f(t, y), utilizando control de error y paso variable.
 * Está basado en el algoritmo de integración con predicción-corrección
 * y control de orden. El integrador se detiene en `tout`.
 *
 * @param func Función del sistema EDO a integrar. Debe recibir (t, y) y devolver dy/dt.
 * @param t Tiempo inicial (modificado al terminar la integración).
 * @param tout Tiempo final de integración.
 * @param relerr Tolerancia relativa de error.
 * @param abserr Tolerancia absoluta de error.
 * @param n_eqn Número de ecuaciones del sistema.
 * @param y Vector de estado inicial (modificado con el valor final).
 * @return Referencia a un `Matrix` con el estado integrado en `tout`.
 *
 * @note Si `t == tout`, se devuelve el estado original. En caso de errores de configuración
 *       (tolerancias negativas, valores inválidos), se devuelve sin integración.
 *
 * @warning Este integrador es sensible a errores numéricos y requiere vectores columna.
 */
//------------------------------------------------------------------------------

Matrix &DEInteg(Matrix &func(double, Matrix &), double t, double tout, double relerr, double abserr, int n_eqn, Matrix &y)
{
    y = transpose(y);

    double twou = 2 * SAT_Const::eps;
    double fouru = 4 * SAT_Const::eps;

    struct
    {
        int DE_INIT = 1;
        int DE_DONE = 2;
        int DE_BADACC = 3;
        int DE_NUMSTEPS = 4;
        int DE_STIFF = 5;
        int DE_INVPARAM = 6;
    } DE_STATE;

    int State_ = DE_STATE.DE_INIT;
    bool PermitTOUT = true;
    double told = 0;

    Matrix &two = zeros(14);

    two(1) = 1.0;
    two(2) = 2.0;
    two(3) = 4.0;
    two(4) = 8.0;
    two(5) = 16.0;
    two(6) = 32.0;
    two(7) = 64.0;
    two(8) = 128.0;
    two(9) = 256.0;
    two(10) = 512.0;
    two(11) = 1024.0;
    two(12) = 2048.0;
    two(13) = 4096.0;
    two(14) = 8192.0;

    Matrix &gstr = zeros(14);

    gstr(1) = 1.0;
    gstr(2) = 0.5;
    gstr(3) = 0.0833;
    gstr(4) = 0.0417;
    gstr(5) = 0.0264;
    gstr(6) = 0.0188;
    gstr(7) = 0.0143;
    gstr(8) = 0.0114;
    gstr(9) = 0.00936;
    gstr(10) = 0.00789;
    gstr(11) = 0.00679;
    gstr(12) = 0.00592;
    gstr(13) = 0.00524;
    gstr(14) = 0.00468;

    Matrix &yy = zeros(n_eqn, 1);
    Matrix &wt = zeros(n_eqn, 1);
    Matrix &p = zeros(n_eqn, 1);
    Matrix &yp = zeros(n_eqn, 1);
    Matrix &phi = zeros(n_eqn, 17);
    Matrix &g = zeros(14, 1);
    Matrix &sig = zeros(14, 1);
    Matrix &rho = zeros(14, 1);
    Matrix &w = zeros(13, 1);
    Matrix &alpha = zeros(13, 1);
    Matrix &beta = zeros(13, 1);
    Matrix &v = zeros(13, 1);
    Matrix &psi_ = zeros(13, 1);

    if (t == tout)
    {
        y = transpose(y);
        return y;
    }

    double epsilon = fmax(relerr, abserr);

    if ((relerr < 0.0) ||
        (abserr < 0.0) ||
        (epsilon <= 0.0) ||
        (State_ > DE_STATE.DE_INVPARAM) ||
        ((State_ != DE_STATE.DE_INIT) &&
         (t != told)))
    {
        State_ = DE_STATE.DE_INVPARAM;
        y = transpose(y);
        return y;
    }
    double del = tout - t;
    double absdel = fabs(del);

    double tend = t + 100.0 * del;
    if (!PermitTOUT)
    {
        tend = tout;
    }

    int nostep = 0;
    int kle4 = 0;
    bool stiff = false;
    double releps = relerr / epsilon;
    double abseps = abserr / epsilon;

    bool OldPermit = PermitTOUT;
    double delsgn = 0.0;
    bool start = false;
    double x = 0.0;
    double h;

    if ((State_ == DE_STATE.DE_INIT) || (!OldPermit) || (delsgn * del <= 0.0))
    {
        start = true;
        x = t;
        yy = y;
        delsgn = sign_(1.0, del);
        h = sign_(fmax(fouru * fabs(x), fabs(tout - x)), tout - x);
    }

    double hi;
    double ki;
    int kold = 0;
    double temp1;
    double term;
    double psijm1;
    double gamma;
    double eta;
    int aux;
    bool crash;
    double p5eps;
    int ifail;
    double round;
    double absh;
    double hnew;
    double sum;
    double hold;
    double knew;
    int k;
    bool phase1;
    bool nornd;
    int kp1;
    int kp2;
    int km1;
    int km2;
    int ns;
    int nsp1;
    int realns;
    int im1;
    double temp2;
    int reali;
    double temp4;
    int nsm2;
    int limit1;
    double temp5;
    double temp3;
    int nsp2;
    int limit2;
    double temp6;
    int ip1;
    double tau;
    double xold;
    double erkm2;
    double erkm1;
    double erk;
    double err;
    bool success;
    double rho_;
    double erkp1;
    double r;
    while (true)
    {
        if (fabs(x - t) >= absdel)
        {

            Matrix &yout = zeros(n_eqn, 1);
            Matrix &ypout = zeros(n_eqn, 1);

            g(2) = 1.0;
            rho(2) = 1.0;

            hi = tout - x;
            ki = kold + 1;

            for (int i = 1; i <= ki; i++)
            {
                temp1 = i;
                w(i + 1) = 1.0 / temp1;
            }

            term = 0.0;

            for (int j = 2; j <= ki; j++)
            {
                psijm1 = psi_(j);
                gamma = (hi + term) / psijm1;
                eta = hi / psijm1;

                for (int i = 1; i <= ki + 1 - j; i++)
                {
                    w(i + 1) = gamma * w(i + 1) - eta * w(i + 2);
                }
                g(j + 1) = w(2);
                rho(j + 1) = gamma * rho(j);
                term = psijm1;
            }
            // Interpolate for the solution yout and for
            // the derivative of the solution ypout
            for (int j = 1; j <= ki; j++)
            {
                int i = ki + 1 - j;

                yout = yout + transpose(transpose(phi.extract_column(i + 1)) * g(i + 1));
                ypout = ypout + transpose(transpose(phi.extract_column(i + 1)) * rho(i + 1));
            }

            yout = y + (yout * hi);
            y = yout;
            State_ = DE_STATE.DE_DONE; // Set return code
            t = tout;                  // Set independent variable
            told = t;                  // Store independent variable
            OldPermit = PermitTOUT;
            y = transpose(y);
            return y;
        }

        if (!PermitTOUT && (fabs(tout - x) < fouru * fabs(x)))
        {
            h = tout - x;
            yp = transpose(func(x, transpose(yy)));
            y = yy + yp * h;
            State_ = DE_STATE.DE_DONE;
            t = tout;
            told = t;
            OldPermit = PermitTOUT;
            y = transpose(y);
            return y;
        }

        h = sign_(fmin(fabs(h), fabs(tend - x)), h);
        for (int l = 1; l <= n_eqn; l++)
        {
            wt(l) = releps * fabs(yy(l)) + abseps;
        }

        bool crash;

        if (fabs(h) < fouru * fabs(x))
        {
            h = sign_(fouru * fabs(x), h);
            crash = true;
            y = transpose(y);
            return y;
        }

        p5eps = 0.5 * epsilon;
        crash = false;
        g(2) = 1.0;
        g(3) = 0.5;
        sig(2) = 1.0;

        ifail = 0;

        double round = 0.0;
        for (int l = 1; l <= n_eqn; l++)
        {
            round = round + (y(l) * y(l)) / (wt(l) * wt(l));
        }
        round = twou * sqrt(round);
        if (p5eps < round)
        {
            epsilon = 2.0 * round * (1.0 + fouru);
            crash = true;
            y = transpose(y);
            return y;
        }
        double absh;

        if (start)
        {
            yp = transpose(func(x, transpose(y)));
            double sum = 0.0;
            for (int l = 1; l < n_eqn; l++)
            {
                phi(l, 2) = yp(l);
                phi(l, 3) = 0.0;
                sum = sum + (yp(l) * yp(l)) / (wt(l) * wt(l));
            }

            sum = sqrt(sum);
            absh = fabs(h);
            if (epsilon < 16.0 * sum * h * h)
            {
                absh = 0.25 * sqrt(epsilon / sum);
            }

            h = sign_(fmax(absh, fouru * fabs(x)), h);
            hold = 0.0;
            hnew = 0.0;
            k = 1;
            kold = 0;
            start = false;
            phase1 = true;
            nornd = true;
            if (p5eps <= 100.0 * round)
            {
                nornd = false;
                for (int l = 1; l <= n_eqn; l++)
                {
                    phi(l, 16) = 0.0;
                }
            }
        }

        while (true)
        {
            kp1 = k + 1;
            kp2 = k + 2;
            km1 = k - 1;
            km2 = k - 2;

            if (h != hold)
            {
                ns = 0;
            }

            if (ns <= kold)
            {
                ns = ns + 1;
            }

            nsp1 = ns + 1;

            if (k >= ns)
            {

                beta(ns + 1) = 1.0;
                realns = ns;
                alpha(ns + 1) = 1.0 / realns;
                temp1 = h * realns;
                sig(nsp1 + 1) = 1.0;
                if (k >= nsp1)
                {
                    for (int i = nsp1; i <= k; i++)
                    {
                        im1 = i - 1;
                        temp2 = psi_(im1 + 1);
                        psi_(im1 + 1) = temp1;
                        beta(i + 1) = beta(im1 + 1) * psi_(im1 + 1) / temp2;
                        temp1 = temp2 + h;
                        alpha(i + 1) = h / temp1;
                        reali = i;
                        sig(i + 2) = reali * alpha(i + 1) * sig(i + 1);
                    }
                }
                psi_(k + 1) = temp1;

                int i;
                if (ns > 1)
                {
                    if (k > kold)
                    {
                        temp4 = k * kp1;
                        v(k + 1) = 1.0 / temp4;
                        int nsm2 = ns - 2;
                        for (int j = 1; j <= nsm2; j++)
                        {
                            i = k - j;
                            v(i + 1) = v(i + 1) - alpha(j + 2) * v(i + 2);
                        }
                    }

                    limit1 = kp1 - ns;
                    temp5 = alpha(ns + 1);
                    for (int iq = 1; iq <= limit1; iq++)
                    {
                        v(iq + 1) = v(iq + 1) - temp5 * v(iq + 2);
                        w(iq + 1) = v(iq + 1);
                    }
                    g(nsp1 + 1) = w(2);
                }
                else
                {
                    for (int iq = 1; iq <= k; iq++)
                    {
                        temp3 = iq * (iq + 1);
                        v(iq + 1) = 1.0 / temp3;
                        w(iq + 1) = v(iq + 1);
                    }
                    nsp2 = ns + 2;
                    if (kp1 >= nsp2)
                    {
                        for (int i = nsp2; i <= kp1; i++)
                        {
                            limit2 = kp2 - i;
                            temp6 = alpha(i);
                            for (int iq = 1; iq <= limit2; iq++)
                            {
                                w(iq + 1) = w(iq + 1) - temp6 * w(iq + 2);
                            }

                            g(i + 1) = w(2);
                        }
                    }
                }
                if (k >= nsp1)
                {
                    for (int i = nsp1; i <= k; i++)
                    {
                        temp1 = beta(i + 1);
                        for (int l = 1; l <= n_eqn; l++)
                        {
                            phi(l, i + 1) = temp1 * phi(l, i + 1);
                        }
                    }
                }
                for (int l = 1; l <= n_eqn; l++)
                {
                    phi(l, kp2 + 1) = phi(l, kp1 + 1);
                    phi(l, kp1 + 1) = 0.0;
                    p(l) = 0.0;
                }
                for (int j = 1; j <= k; j++)
                {
                    int i = kp1 - j;
                    int ip1 = i + 1;
                    temp2 = g(i + 1);
                    for (int l = 1; l <= n_eqn; l++)
                    {
                        p(l) = p(l) + temp2 * phi(l, i + 1);
                        phi(l, i + 1) = phi(l, i + 1) + phi(l, ip1 + 1);
                    }
                }
                if (nornd)
                {
                    p = y + (p * h);
                }
                else
                {
                    for (int l = 1; l <= n_eqn; l++)
                    {
                        tau = h * p(l) - phi(l, 16);
                        p(l) = y(l) + tau;
                        phi(l, 17) = (p(l) - y(l)) - tau;
                    }
                }
                xold = x;
                x = x + h;
                absh = fabs(h);
                yp = transpose(func(x, transpose(p)));

                erkm2 = 0.0;
                erkm1 = 0.0;
                erk = 0.0;
                for (int l = 1; l <= n_eqn; l++)
                {
                    temp3 = 1.0 / wt(l);
                    temp4 = yp(l) - phi(l, 1 + 1);
                    if (km2 > 0)
                    {
                        erkm2 = erkm2 + ((phi(l, km1 + 1) + temp4) * temp3) * ((phi(l, km1 + 1) + temp4) * temp3);
                    }

                    if (km2 >= 0)
                    {
                        erkm1 = erkm1 + ((phi(l, k + 1) + temp4) * temp3) * ((phi(l, k + 1) + temp4) * temp3);
                    }

                    erk = erk + (temp4 * temp3) * (temp4 * temp3);
                }
                if (km2 > 0)
                {
                    erkm2 = absh * sig(km1 + 1) * gstr(km2 + 1) * sqrt(erkm2);
                }

                if (km2 >= 0)
                {
                    erkm1 = absh * sig(k + 1) * gstr(km1 + 1) * sqrt(erkm1);
                }

                temp5 = absh * sqrt(erk);
                err = temp5 * (g(k + 1) - g(kp1 + 1));
                erk = temp5 * sig(kp1 + 1) * gstr(k + 1);
                knew = k;
                if (km2 > 0)
                {
                    if (fmax(erkm1, erkm2) <= erk)
                    {
                        knew = km1;
                    }
                }

                if (km2 == 0)
                {
                    if (erkm1 <= 0.5 * erk)
                    {
                        knew = km1;
                    }
                }

                success = (err <= epsilon);

                if (!success)
                {

                    phase1 = false;
                    x = xold;
                    for (int i = 1; i <= k; i++)
                    {
                        temp1 = 1.0 / beta(i + 1);
                        ip1 = i + 1;
                        for (int l = 1; l <= n_eqn; l++)
                        {
                            phi(l, i + 1) = temp1 * (phi(l, i + 1) - phi(l, ip1 + 1));
                        }
                    }

                    if (k >= 2)
                    {
                        for (int i = 2; i <= k; i++)
                        {
                            psi_(i) = psi_(i + 1) - h;
                        }
                    }

                    ifail = ifail + 1;
                    temp2 = 0.5;
                    if (ifail > 3)
                    {
                        if (p5eps < 0.25 * erk)
                        {
                            temp2 = sqrt(p5eps / erk);
                        }
                    }

                    if (ifail >= 3)
                    {
                        knew = 1;
                    }

                    h = temp2 * h;
                    k = knew;
                    if (fabs(h) < fouru * fabs(x))
                    {
                        crash = true;
                        h = sign_(fouru * fabs(x), h);
                        epsilon = epsilon * 2.0;
                        y = transpose(y);
                        return y;
                    }
                }
                if (success)
                {
                    break;
                }
            }

            kold = k;
            hold = h;

            temp1 = h * g(kp1 + 1);
            if (nornd)
            {
                for (int l = 1; l <= n_eqn; l++)
                {
                    y(l) = p(l) + temp1 * (yp(l) - phi(l, 2));
                }
            }
            else
            {
                for (int l = 1; l <= n_eqn; l++)
                {
                    rho_ = temp1 * (yp(l) - phi(l, 2)) - phi(l, 17);
                    y(l) = p(l) + rho_;
                    phi(l, 16) = (y(l) - p(l)) - rho_;
                }
            }

            yp = transpose(func(x, transpose(y)));

            for (int l = 1; l <= n_eqn; l++)
            {
                phi(l, kp1 + 1) = yp(l) - phi(l, 2);
                phi(l, kp2 + 1) = phi(l, kp1 + 1) - phi(l, kp2 + 1);
            }

            for (int i = 1; i <= k; i++)
            {
                for (int l = 1; l <= n_eqn; l++)
                {
                    phi(l, i + 1) = phi(l, i + 1) + phi(l, kp1 + 1);
                }
            }

            erkp1 = 0.0;
            if ((knew == km1) || (k == 12))
            {
                phase1 = false;
            }

            if (phase1)
            {
                k = kp1;
                erk = erkp1;
            }
            else
            {
                if (knew == km1)
                {
                    k = km1;
                    erk = erkm1;
                }
                else
                {
                    if (kp1 <= ns)
                        for (int l = 1; l <= n_eqn; l++)
                        {
                            erkp1 = erkp1 + (phi(l, kp2 + 1) / wt(l)) * (phi(l, kp2 + 1) / wt(l));
                        }

                    erkp1 = absh * gstr(kp1 + 1) * sqrt(erkp1);

                    if (k > 1)
                    {
                        if (erkm1 <= fmin(erk, erkp1))
                        {
                            k = km1;
                            erk = erkm1;
                        }
                        else
                        {
                            if ((erkp1 < erk) && (k != 12))
                            {
                                k = kp1;
                                erk = erkp1;
                            }
                        }
                    }
                    else if (erkp1 < 0.5 * erk)
                    {
                        k = kp1;
                        erk = erkp1;
                    }
                }
            }
            double hnew;
            if (phase1 || (p5eps >= erk * two(k + 2)))
            {
                hnew = 2.0 * h;
            }
            else
            {
                if (p5eps < erk)
                {
                    temp2 = k + 1;
                    r = p5eps / pow(erk, 1.0 / temp2);
                    hnew = absh * fmax(0.5, fmin(0.9, r));
                    hnew = sign_(fmax(hnew, fouru * fabs(x)), h);
                }
                else
                {
                    hnew = h;
                }
            }

            h = hnew;

            if (crash)
            {
                State_ = DE_STATE.DE_BADACC;
                relerr = epsilon * releps;
                abserr = epsilon * abseps;
                y = yy;
                t = x;
                told = t;
                OldPermit = true;
                y = transpose(y);
                return y;
            }

            nostep = nostep + 1;

            kle4 = kle4 + 1;
            if (kold > 4)
            {
                kle4 = 0;
            }

            if (kle4 >= 50)
            {
                stiff = true;
            }
        }
        y = transpose(y);
        return y;
    }
}