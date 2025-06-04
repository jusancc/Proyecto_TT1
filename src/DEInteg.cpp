#include "../include/DEInteg.hpp"

using namespace std;

Matrix &DEInteg(Matrix &func(double, Matrix &), double &t, double tout, double &relerr, double &abserr, int n_eqn, Matrix &y)
{

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
    int told = 0;

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
        return transpose(y);
    }

    double epsilon = max(relerr, abserr);

    if ((relerr < 0.0) ||
        (abserr < 0.0) ||
        (epsilon <= 0.0) ||
        (State_ > DE_STATE.DE_INVPARAM) ||
        ((State_ != DE_STATE.DE_INIT) &&
         (t != told)))
    {
        State_ = DE_STATE.DE_INVPARAM;
        return transpose(y);
    }

    double del = tout - t;
    double absdel = abs(del);

    double tend = t + 100.0 * del;
    if (~PermitTOUT)
    {
        tend = tout;
    }

    int kold = 0;
    double nostep = 0;
    double kle4 = 0;
    bool stiff = false;
    double releps = relerr / epsilon;
    double abseps = abserr / epsilon;

    bool nornd = false;
    bool phase1 = false;
    int knew = 0;
    int km1 = 0;
    double erkp1, erk;
    int ns = 0;
    double hold;

    bool OldPermit = false;
    double delsgn = 0.0;
    double x = 0.0;
    double h;
    bool start = false;

    if ((State_ == DE_STATE.DE_INIT) || (~OldPermit) || (delsgn * del <= 0.0))
    {
        start = true;
        x = t;
        yy = y;
        delsgn = sign_(1.0, del);
        h = sign_(max(fouru * abs(x), abs(tout - x)), tout - x);
    }

    int k=0;
    while (true)
    {
        double temp1;

        if (fabs(x - t) >= absdel)
        {
            Matrix &yout = zeros(n_eqn, 1);
            Matrix &ypout = zeros(n_eqn, 1);
            g(2) = 1.0;
            rho(2) = 1.0;
            double hi = tout - x;
            int ki = kold + 1;
            // Initialize w[*] for computing g[*]
            for (int i = 1; i <= ki; i++)
            {
                temp1 = i;
                w(i + 1) = 1.0 / temp1;
            }
            // Compute g[*]
            double term = 0.0;
            for (int j = 2; j <= ki; j++)
            {
                double psijm1 = psi_(j);
                double gamma = (hi + term) / psijm1;
                double eta = hi / psijm1;
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
                yout = yout + (transpose(phi.extract_column(i + 1)) * g(i + 1));
                ypout = ypout + (transpose(phi.extract_column(i + 1)) * rho(i + 1));
            }
            yout = y + (yout * hi);
            y = yout;
            State_ = DE_STATE.DE_DONE; // Set return code
            t = tout;                  // Set independent variable
            told = t;                  // Store independent variable
            OldPermit = PermitTOUT;
            return transpose(y); // Normal exit
        }

        if (~PermitTOUT && (abs(tout - x) < fouru * abs(x)))
        {
            h = tout - x;
            yp = func(x, yy);
            y = yy + yp * h;
            State_ = DE_STATE.DE_DONE;
            t = tout;
            told = t;
            OldPermit = PermitTOUT;
            return transpose(y);
        }

        h = sign_(min(abs(h), abs(tend - x)), h);
        for (int l = 1; l < n_eqn; l++)
        {
            wt(l) = releps * abs(yy(l)) + abseps;
        }

        bool crash;

        if (abs(h) < fouru * abs(x))
        {
            h = sign_(fouru * abs(x), h);
            crash = true;
            return transpose(y);
        }

        double p5eps = 0.5 * epsilon;
        crash = false;
        g(2) = 1.0;
        g(3) = 0.5;
        sig(2) = 1.0;

        int ifail = 0;

        double round = 0.0;
        for (int l = 1; l < n_eqn; l++)
        {
            round = round + (y(l) * y(l)) / (wt(l) * wt(l));
        }

        round = twou * sqrt(round);
        if (p5eps < round)
        {
            epsilon = 2.0 * round * (1.0 + fouru);
            crash = true;
            return transpose(y);
        }
        double absh;

        if (start)
        {

            yp = func(x, y);
            double sum = 0.0;
            for (int l = 1; l < n_eqn; l++)
            {
                phi(l, 2) = yp(l);
                phi(l, 3) = 0.0;
                sum = sum + (yp(l) * yp(l)) / (wt(l) * wt(l));
            }

            sum = sqrt(sum);
            absh = abs(h);
            if (epsilon < 16.0 * sum * h * h)
            {
                absh = 0.25 * sqrt(epsilon / sum);
            }

            h = sign_(max(absh, fouru * abs(x)), h);
            double hold = 0.0;
            double hnew = 0.0;
            k = 1;
            kold = 0;
            start = false;
            bool phase1 = true;
            bool nornd = true;
            if (p5eps <= 100.0 * round)
            {
                nornd = false;
                for (int l = 1; l < n_eqn; l++)
                {
                    phi(l, 16) = 0.0;
                }
            }
        }

        double km1;
        double km2;
        int kp1 = 0;
        int kp2 = 0;
        double erkm1 = 0.0;
        double erk = 0.0;

        double temp2, temp3, temp4, temp5, temp6;
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

            double nsp1 = ns + 1;

            if (k >= ns)
            {

                beta(ns + 1) = 1.0;
                int realns = ns;
                alpha(ns + 1) = 1.0 / realns;
                temp1 = h * realns;
                sig(nsp1 + 1) = 1.0;
                if (k >= nsp1)
                {
                    for (int i = nsp1; i < k; i++)
                    {
                        int im1 = i - 1;
                        double temp2 = psi_(im1 + 1);
                        psi_(im1 + 1) = temp1;
                        beta(i + 1) = beta(im1 + 1) * psi_(im1 + 1) / temp2;
                        temp1 = temp2 + h;
                        alpha(i + 1) = h / temp1;
                        int reali = i;
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
                        for (int j = 1; j < nsm2; j++)
                        {
                            i = k - j;
                            v(i + 1) = v(i + 1) - alpha(j + 2) * v(i + 2);
                        }
                    }

                    double limit1 = kp1 - ns;
                    temp5 = alpha(ns + 1);
                    for (int iq = 1; iq < limit1; iq++)
                    {
                        v(iq + 1) = v(iq + 1) - temp5 * v(iq + 2);
                        w(iq + 1) = v(iq + 1);
                    }

                    g(nsp1 + 1) = w(2);
                }
                else
                {
                    for (int iq = 1; iq < k; iq++)
                    {
                        temp3 = iq * (iq + 1);
                        v(iq + 1) = 1.0 / temp3;
                        w(iq + 1) = v(iq + 1);
                    }
                }

                int nsp2 = ns + 2;
                if (kp1 >= nsp2)
                {
                    for (int i = nsp2; i < kp1; i++)
                    {
                        double limit2 = kp2 - i;
                        temp6 = alpha(i);
                        for (int iq = 1; iq < limit2; iq++)
                        {
                            w(iq + 1) = w(iq + 1) - temp6 * w(iq + 2);
                        }

                        g(i + 1) = w(2);
                    }
                }
            }

            if (k >= nsp1)
            {
                for (int i = nsp1; i < k; i++)
                {
                    temp1 = beta(i + 1);
                    for (int l = 1; l < n_eqn; l++)
                    {
                        phi(l, i + 1) = temp1 * phi(l, i + 1);
                    }
                }
            }

            for (int l = 1; l < n_eqn; l++)
            {
                phi(l, kp2 + 1) = phi(l, kp1 + 1);
                phi(l, kp1 + 1) = 0.0;
                p(l) = 0.0;
            }

            for (int j = 1; j < k; j++)
            {
                int i = kp1 - j;
                int ip1 = i + 1;
                double temp2 = g(i + 1);
                for (int l = 1; l < n_eqn; l++)
                {
                    p(l) = p(l) + temp2 * phi(l, i + 1);
                    phi(l, i + 1) = phi(l, i + 1) + phi(l, ip1 + 1);
                }
            }

            if (nornd)
            {
                p = y + p * h;
            }
            else
            {
                for (int l = 1; l < n_eqn; l++)
                {
                    double tau = h * p(l) - phi(l, 16);
                    p(l) = y(l) + tau;
                    phi(l, 17) = (p(l) - y(l)) - tau;
                }
            }

            double xold = x;
            x = x + h;
            double absh = abs(h);
            yp = func(x, p);

            double erkm2 = 0.0;
            double erkm1 = 0.0;
            double erk = 0.0;

            for (int l = 1; l < n_eqn; l++)
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
            double err = temp5 * (g(k + 1) - g(kp1 + 1));
            erk = temp5 * sig(kp1 + 1) * gstr(k + 1);
            knew = k;

            if (km2 > 0)
            {
                if (max(erkm1, erkm2) <= erk)
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

            bool success = (err <= epsilon);

            if (~success)

                phase1 = false;
            x = xold;
            for (int i = 1; i < k; i++)
            {
                temp1 = 1.0 / beta(i + 1);
                int ip1 = i + 1;
                for (int l = 1; l < n_eqn; l++)
                {
                    phi(l, i + 1) = temp1 * (phi(l, i + 1) - phi(l, ip1 + 1));
                }
            }

            if (k >= 2)
            {
                for (int i = 2; i < k; i++)
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
            if (abs(h) < fouru * abs(x))
            {
                crash = true;
                h = sign_(fouru * abs(x), h);
                epsilon = epsilon * 2.0;
                return transpose(y);
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
            for (int l = 1; l < n_eqn; l++)
            {
                y(l) = p(l) + temp1 * (yp(l) - phi(l, 2));
            }
        }
        else
        {
            for (int l = 1; l <= n_eqn; l++)
            {
                double rho_ = temp1 * (yp(l) - phi(l, 2)) - phi(l, 17);
                y(l) = p(l) + rho_;
                phi(l, 16) = (y(l) - p(l)) - rho_;
            }
        }

        yp = func(x, y);

        for (int l = 1; l < n_eqn; l++)
        {
            phi(l, kp1 + 1) = yp(l) - phi(l, 2);
            phi(l, kp2 + 1) = phi(l, kp1 + 1) - phi(l, kp2 + 1);
        }

        for (int i = 1; i < k; i++)
        {
            for (int l = 1; l < n_eqn; l++)
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
                    for (int l = 1; l < n_eqn; l++)
                    {
                        erkp1 = erkp1 + (phi(l, kp2 + 1) / wt(l)) * (phi(l, kp2 + 1) / wt(l));
                    }

                erkp1 = absh * gstr(kp1 + 1) * sqrt(erkp1);

                if (k > 1)
                {
                    if (erkm1 <= min(erk, erkp1))
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
                double r = p5eps / pow(erk, 1.0 / temp2);
                hnew = absh * max(0.5, min(0.9, r));
                hnew = sign_(max(hnew, fouru * abs(x)), h);
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
            return transpose(y);
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
}