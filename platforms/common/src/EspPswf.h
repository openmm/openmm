/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2026 Stanford University and the Authors.           *
 * Authors: Libin Lu                                                          *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

// Prolate spheroidal wave function coefficient builder for ESP.  This runs at
// Context initialization; GPU kernels only use the generated polynomial tables.

#ifndef OPENMM_ESP_PSWF_H_
#define OPENMM_ESP_PSWF_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pswf {

// ============================================================================
// Horner polynomial evaluation (runtime)
// ============================================================================

static inline double horner_eval(const double* coeffs, int n, double x) {
    if (n <= 0) return 0.0;
    double y = coeffs[0];
    for (int i = 1; i < n; ++i)
        y = y * x + coeffs[i];
    return y;
}

static inline float horner_eval_f(const float* coeffs, int n, float x) {
    if (n <= 0) return 0.0f;
    float y = coeffs[0];
    for (int i = 1; i < n; ++i)
        y = y * x + coeffs[i];
    return y;
}


static inline float clenshaw_eval_f(const float* coeffs, int n, float x) {
    if (n <= 0) return 0.0f;
    const float u = 2.0f*x - 1.0f;
    float b_kplus1 = 0.0f;
    float b_kplus2 = 0.0f;
    for (int k = n-1; k >= 1; --k) {
        const float b_k = 2.0f*u*b_kplus1 - b_kplus2 + coeffs[k];
        b_kplus2 = b_kplus1;
        b_kplus1 = b_k;
    }
    return u*b_kplus1 - b_kplus2 + coeffs[0];
}

// ============================================================================
// Legendre polynomial helpers (literal from GROMACS pswf.cpp)
// ============================================================================

static inline void legepol(double x, int n, double& pol, double& der) {
    double pkm1 = 1.0;
    double pk = x;
    double pkp1;

    if (n == 0) {
        pol = 1.0;
        der = 0.0;
        return;
    }

    if (n == 1) {
        pol = x;
        der = 1.0;
        return;
    }

    pk = 1.0;
    pkp1 = x;

    for (int k = 1; k < n; ++k) {
        pkm1 = pk;
        pk = pkp1;
        pkp1 = ((2 * k + 1) * x * pk - k * pkm1) / (k + 1);
    }

    pol = pkp1;
    der = n * (x * pkp1 - pk) / (x * x - 1);
}

static inline void legetayl(double pol, double der, double x, double h,
                             int n, int k, double& sum, double& sumder) {
    double done = 1.0;
    double q0 = pol;
    double q1 = der * h;
    double q2 = (2 * x * der - n * (n + done) * pol) / (1 - x * x);
    q2 = q2 * h * h / 2;

    sum = q0 + q1 + q2;
    sumder = q1 / h + q2 * 2 / h;

    if (k <= 2) return;

    double qi = q1;
    double qip1 = q2;

    for (int i = 1; i <= k - 2; ++i) {
        double d = 2 * x * (i + 1) * (i + 1) / h * qip1 - (n * (n + done) - i * (i + 1)) * qi;
        d = d / (i + 1) / (i + 2) * h * h / (1 - x * x);
        double qip2 = d;

        sum += qip2;
        sumder += d * (i + 2) / h;

        qi = qip1;
        qip1 = qip2;
    }
}

// Gauss-Legendre quadrature: itype=1 computes both roots and weights.
static inline void legerts(int itype, int n, double* ts, double* whts) {
    int k = 30;
    double d = 1.0;
    double d2 = d + 1.0e-24;
    if (d2 != d) {
        k = 54;
    }

    int half = n / 2;
    int ifodd = n - 2 * half;
    double pi_val = atan(1.0) * 4.0;
    double h = pi_val / (2.0 * n);

    int ii = 0;
    for (int i = 1; i <= n; i++) {
        if (i < (n / 2 + 1)) {
            continue;
        }
        ii++;
        double t = (2.0 * i - 1.0) * h;
        ts[ii - 1] = -cos(t);
    }

    double pol = 1.0, der = 0.0;
    double x0 = 0.0;
    legepol(x0, n, pol, der);
    double x1 = ts[0];

    int n2 = (n + 1) / 2;
    double pol3 = pol, der3 = der;

    for (int kk = 1; kk <= n2; kk++) {
        if ((ifodd == 1) && (kk == 1)) {
            ts[kk - 1] = x0;
            if (itype > 0) {
                whts[kk - 1] = der;
            }
            x0 = x1;
            x1 = ts[kk];
            pol3 = pol;
            der3 = der;
            continue;
        }

        int ifstop = 0;
        for (int i = 1; i <= 10; i++) {
            double hh = x1 - x0;

            legetayl(pol3, der3, x0, hh, n, k, pol, der);
            x1 = x1 - pol / der;

            if (fabs(pol) < 1.0e-12) {
                ifstop++;
            }
            if (ifstop == 3) {
                break;
            }
        }

        ts[kk - 1] = x1;
        if (itype > 0) {
            whts[kk - 1] = der;
        }

        x0 = x1;
        x1 = ts[kk];
        pol3 = pol;
        der3 = der;
    }

    for (int i = n2; i >= 1; i--) {
        ts[i - 1 + half] = ts[i - 1];
    }
    for (int i = 1; i <= half; i++) {
        ts[i - 1] = -ts[n - i];
    }
    if (itype <= 0) {
        return;
    }

    for (int i = n2; i >= 1; i--) {
        whts[i - 1 + half] = whts[i - 1];
    }
    for (int i = 1; i <= half; i++) {
        whts[i - 1] = whts[n - i];
    }

    for (int i = 0; i < n; i++) {
        double tmp = 1.0 - ts[i] * ts[i];
        whts[i] = 2.0 / tmp / (whts[i] * whts[i]);
    }
}

// Convenience: compute Gauss-Legendre nodes and weights.
static inline void legerts(int n, double* ts, double* whts) {
    legerts(1, n, ts, whts);
}

static inline void legepols(double x, int n, double* pols) {
    double pkm1 = 1.0;
    double pk = x;

    if (n == 0) {
        pols[0] = 1.0;
        return;
    }

    if (n == 1) {
        pols[0] = 1.0;
        pols[1] = x;
        return;
    }

    pols[0] = 1.0;
    pols[1] = x;

    for (int k = 1; k < n; ++k) {
        double pkp1 = ((2 * k + 1) * x * pk - k * pkm1) / (k + 1);
        pols[k + 1] = pkp1;
        pkm1 = pk;
        pk = pkp1;
    }
}

static inline void legeexps(int itype, int n, double* x,
                             std::vector<std::vector<double>>& u,
                             std::vector<std::vector<double>>& v, double* whts) {
    int itype_rts = (itype > 0) ? 1 : 0;
    legerts(itype_rts, n, x, whts);

    if (itype != 2) return;

    u.resize(n, std::vector<double>(n));
    v.resize(n, std::vector<double>(n));

    for (int i = 0; i < n; ++i) {
        std::vector<double> pols(n);
        legepols(x[i], n - 1, pols.data());
        for (int j = 0; j < n; ++j) {
            u[j][i] = pols[j];
        }
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            v[i][j] = u[j][i];
        }
    }

    for (int i = 0; i < n; ++i) {
        double dd = 1.0 * (2 * (i + 1) - 1) / 2;
        for (int j = 0; j < n; ++j) {
            u[i][j] = v[j][i] * whts[j] * dd;
        }
    }
}

static inline void legeexev(double x, double& val, const double* pexp, int n) {
    double pjm2 = 1.0;
    double pjm1 = x;

    val = pexp[0] * pjm2 + pexp[1] * pjm1;

    for (int j = 2; j <= n; ++j) {
        double pj = ((2 * j - 1) * x * pjm1 - (j - 1) * pjm2) / j;
        val += pexp[j] * pj;
        pjm2 = pjm1;
        pjm1 = pj;
    }
}

static inline void legeFDER(double x, double& val, double& der, const double* pexp, int n) {
    double pjm2 = 1.0;
    double pjm1 = x;
    double derjm2 = 0.0;
    double derjm1 = 1.0;

    val = pexp[0] * pjm2 + pexp[1] * pjm1;
    der = pexp[1];

    for (int j = 2; j <= n; ++j) {
        double pj = ((2 * j - 1) * x * pjm1 - (j - 1) * pjm2) / j;
        val += pexp[j] * pj;

        double derj = (2 * j - 1) * (pjm1 + x * derjm1) - (j - 1) * derjm2;
        derj /= j;
        der += pexp[j] * derj;

        pjm2 = pjm1;
        pjm1 = pj;
        derjm2 = derjm1;
        derjm1 = derj;
    }
}

// ============================================================================
// Prolate eigensolver (literal from GROMACS pswf.cpp)
// ============================================================================

static inline void legeFDER2(double x, double& val, double& der, double& der2, const double* pexp, int n) {
    double pjm2 = 1.0;
    double pjm1 = x;
    double derjm2 = 0.0;
    double derjm1 = 1.0;
    double der2jm2 = 0.0;
    double der2jm1 = 0.0;

    val = pexp[0] * pjm2 + pexp[1] * pjm1;
    der = pexp[1];
    der2 = 0.0;

    for (int j = 2; j <= n; ++j) {
        double pj = ((2 * j - 1) * x * pjm1 - (j - 1) * pjm2) / j;
        val += pexp[j] * pj;

        double derj = (2 * j - 1) * (pjm1 + x * derjm1) - (j - 1) * derjm2;
        derj /= j;
        der += pexp[j] * derj;

        double der2j = (2 * j - 1) * (2.0 * derjm1 + x * der2jm1) - (j - 1) * der2jm2;
        der2j /= j;
        der2 += pexp[j] * der2j;

        pjm2 = pjm1;
        pjm1 = pj;
        derjm2 = derjm1;
        derjm1 = derj;
        der2jm2 = der2jm1;
        der2jm1 = der2j;
    }
}

static inline void legeFDER3(double x, double& val, double& der, double& der2, double& der3, const double* pexp, int n) {
    double pjm2 = 1.0;
    double pjm1 = x;
    double derjm2 = 0.0;
    double derjm1 = 1.0;
    double der2jm2 = 0.0;
    double der2jm1 = 0.0;
    double der3jm2 = 0.0;
    double der3jm1 = 0.0;

    val = pexp[0] * pjm2 + pexp[1] * pjm1;
    der = pexp[1];
    der2 = 0.0;
    der3 = 0.0;

    for (int j = 2; j <= n; ++j) {
        double pj = ((2 * j - 1) * x * pjm1 - (j - 1) * pjm2) / j;
        val += pexp[j] * pj;

        double derj = (2 * j - 1) * (pjm1 + x * derjm1) - (j - 1) * derjm2;
        derj /= j;
        der += pexp[j] * derj;

        double der2j = (2 * j - 1) * (2.0 * derjm1 + x * der2jm1) - (j - 1) * der2jm2;
        der2j /= j;
        der2 += pexp[j] * der2j;

        double der3j = (2 * j - 1) * (3.0 * der2jm1 + x * der3jm1) - (j - 1) * der3jm2;
        der3j /= j;
        der3 += pexp[j] * der3j;

        pjm2 = pjm1;
        pjm1 = pj;
        derjm2 = derjm1;
        derjm1 = derj;
        der2jm2 = der2jm1;
        der2jm1 = der2j;
        der3jm2 = der3jm1;
        der3jm1 = der3j;
    }
}

static inline void legeFDER4(double x, double& val, double& der, double& der2, double& der3, double& der4, const double* pexp, int n) {
    double pjm2 = 1.0;
    double pjm1 = x;
    double derjm2 = 0.0;
    double derjm1 = 1.0;
    double der2jm2 = 0.0;
    double der2jm1 = 0.0;
    double der3jm2 = 0.0;
    double der3jm1 = 0.0;
    double der4jm2 = 0.0;
    double der4jm1 = 0.0;

    val = pexp[0] * pjm2 + pexp[1] * pjm1;
    der = pexp[1];
    der2 = 0.0;
    der3 = 0.0;
    der4 = 0.0;

    for (int j = 2; j <= n; ++j) {
        double pj = ((2 * j - 1) * x * pjm1 - (j - 1) * pjm2) / j;
        val += pexp[j] * pj;

        double derj = (2 * j - 1) * (pjm1 + x * derjm1) - (j - 1) * derjm2;
        derj /= j;
        der += pexp[j] * derj;

        double der2j = (2 * j - 1) * (2.0 * derjm1 + x * der2jm1) - (j - 1) * der2jm2;
        der2j /= j;
        der2 += pexp[j] * der2j;

        double der3j = (2 * j - 1) * (3.0 * der2jm1 + x * der3jm1) - (j - 1) * der3jm2;
        der3j /= j;
        der3 += pexp[j] * der3j;

        double der4j = (2 * j - 1) * (4.0 * der3jm1 + x * der4jm1) - (j - 1) * der4jm2;
        der4j /= j;
        der4 += pexp[j] * der4j;

        pjm2 = pjm1;
        pjm1 = pj;
        derjm2 = derjm1;
        derjm1 = derj;
        der2jm2 = der2jm1;
        der2jm1 = der2j;
        der3jm2 = der3jm1;
        der3jm1 = der3j;
        der4jm2 = der4jm1;
        der4jm1 = der4j;
    }
}

static inline void prosinin(double c, const double* ts, const double* whts, const double* fs,
                            double x, int n, double& rint, double& derrint) {
    rint = 0.0;
    derrint = 0.0;

    for (int i = 0; i < n; ++i) {
        double diff = x - ts[i];
        double sin_term = sin(c * diff);
        double cos_term = cos(c * diff);

        rint += whts[i] * fs[i] * sin_term / diff;
        derrint += whts[i] * fs[i] / (diff * diff) * (c * diff * cos_term - sin_term);
    }
}

static inline void prosinin2(double c, const double* ts, const double* whts, const double* fs,
                             double x, int n, double& rint, double& derrint, double& der2rint) {
    rint = 0.0;
    derrint = 0.0;
    der2rint = 0.0;

    for (int i = 0; i < n; ++i) {
        double diff = x - ts[i];
        double sin_term = sin(c * diff);
        double cos_term = cos(c * diff);

        rint += whts[i] * fs[i] * sin_term / diff;
        derrint += whts[i] * fs[i] / (diff * diff) * (c * diff * cos_term - sin_term);
        der2rint += whts[i] * fs[i] * ((2.0 - c * c * diff * diff) * sin_term - 2.0 * c * diff * cos_term) / (diff * diff * diff);
    }
}

static inline void prosinin3(double c, const double* ts, const double* whts, const double* fs,
                             double x, int n, double& rint, double& derrint, double& der2rint, double& der3rint) {
    rint = 0.0;
    derrint = 0.0;
    der2rint = 0.0;
    der3rint = 0.0;

    for (int i = 0; i < n; ++i) {
        double diff = x - ts[i];
        double sin_term = sin(c * diff);
        double cos_term = cos(c * diff);

        rint += whts[i] * fs[i] * sin_term / diff;
        derrint += whts[i] * fs[i] / (diff * diff) * (c * diff * cos_term - sin_term);
        der2rint += whts[i] * fs[i] * ((2.0 - c * c * diff * diff) * sin_term - 2.0 * c * diff * cos_term) / (diff * diff * diff);
        der3rint += whts[i] * fs[i] * ((3.0 * c * c * diff * diff - 6.0) * sin_term + c * diff * (6.0 - c * c * diff * diff) * cos_term) / (diff * diff * diff * diff);
    }
}

static inline void prosinin4(double c, const double* ts, const double* whts, const double* fs,
                             double x, int n, double& rint, double& derrint, double& der2rint,
                             double& der3rint, double& der4rint) {
    rint = 0.0;
    derrint = 0.0;
    der2rint = 0.0;
    der3rint = 0.0;
    der4rint = 0.0;

    for (int i = 0; i < n; ++i) {
        double diff = x - ts[i];
        double sin_term = sin(c * diff);
        double cos_term = cos(c * diff);
        double diff2 = diff * diff;
        double diff3 = diff2 * diff;
        double diff4 = diff3 * diff;
        double diff5 = diff4 * diff;
        double c2 = c * c;
        double c4 = c2 * c2;

        rint += whts[i] * fs[i] * sin_term / diff;
        derrint += whts[i] * fs[i] / diff2 * (c * diff * cos_term - sin_term);
        der2rint += whts[i] * fs[i] * ((2.0 - c2 * diff2) * sin_term - 2.0 * c * diff * cos_term) / diff3;
        der3rint += whts[i] * fs[i] * ((3.0 * c2 * diff2 - 6.0) * sin_term + c * diff * (6.0 - c2 * diff2) * cos_term) / diff4;
        der4rint += whts[i] * fs[i] *
                    ((c4 * diff4 - 12.0 * c2 * diff2 + 24.0) * sin_term +
                     4.0 * c * diff * (c2 * diff2 - 6.0) * cos_term) / diff5;
    }
}
static inline void prolcoef(double rlam, int k, double c, double& alpha0, double& beta0,
                            double& gamma0, double& alpha, double& beta, double& gamma) {
    double d = k * (k - 1);
    d = d / (2 * k + 1) / (2 * k - 1);
    double uk = d;

    d = (k + 1) * (k + 1);
    d = d / (2 * k + 3);
    double d2 = k * k;
    d2 = d2 / (2 * k - 1);
    double vk = (d + d2) / (2 * k + 1);

    d = (k + 1) * (k + 2);
    d = d / (2 * k + 1) / (2 * k + 3);
    double wk = d;

    alpha = -c * c * uk;
    beta = rlam - k * (k + 1) - c * c * vk;
    gamma = -c * c * wk;

    alpha0 = uk;
    beta0 = vk;
    gamma0 = wk;
}

static inline void prolmatr(double* as, double* bs, double* cs, int n, double c, double rlam,
                            int ifsymm, int ifodd) {
    double done = 1.0;
    double half = done / 2.0;
    int k = 0;

    if (ifodd > 0) {
        for (int k0 = 1; k0 <= n + 2; k0 += 2) {
            k++;
            double alpha0, beta0, gamma0, alpha, beta, gamma;
            prolcoef(rlam, k0, c, alpha0, beta0, gamma0, alpha, beta, gamma);

            as[k - 1] = alpha;
            bs[k - 1] = beta;
            cs[k - 1] = gamma;

            if (ifsymm != 0) {
                if (k0 > 1) {
                    as[k - 1] = as[k - 1] / std::sqrt(k0 - 2 + half) * std::sqrt(k0 + half);
                }
                cs[k - 1] = cs[k - 1] * std::sqrt(k0 + half) / std::sqrt(k0 + half + 2);
            }
        }
    } else {
        for (int k0 = 0; k0 <= n + 2; k0 += 2) {
            k++;
            double alpha0, beta0, gamma0, alpha, beta, gamma;
            prolcoef(rlam, k0, c, alpha0, beta0, gamma0, alpha, beta, gamma);

            as[k - 1] = alpha;
            bs[k - 1] = beta;
            cs[k - 1] = gamma;

            if (ifsymm != 0) {
                if (k0 != 0) {
                    as[k - 1] = as[k - 1] / std::sqrt(k0 - 2 + half) * std::sqrt(k0 + half);
                }
                cs[k - 1] = cs[k - 1] * std::sqrt(k0 + half) / std::sqrt(k0 + half + 2);
            }
        }
    }
}

static inline void prolql1(int n, double* d, double* e, int& ierr) {
    ierr = 0;
    if (n == 1) return;

    for (int i = 1; i < n; ++i) {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;

    for (int l = 0; l < n; ++l) {
        int j = 0;
        while (true) {
            int m;
            for (m = l; m < n - 1; ++m) {
                double tst1 = std::abs(d[m]) + std::abs(d[m + 1]);
                double tst2 = tst1 + std::abs(e[m]);
                if (tst2 == tst1) break;
            }

            if (m == l) break;
            if (j == 30) {
                ierr = l + 1;
                return;
            }
            ++j;

            double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            double r = std::sqrt(g * g + 1.0);
            g = d[m] - d[l] + e[l] / (g + std::copysign(r, g));
            double s = 1.0;
            double c = 1.0;
            double p = 0.0;

            for (int i = m - 1; i >= l; --i) {
                double f = s * e[i];
                double b = c * e[i];
                r = std::sqrt(f * f + g * g);
                e[i + 1] = r;
                if (r == 0.0) {
                    d[i + 1] -= p;
                    e[m] = 0.0;
                    break;
                }
                s = f / r;
                c = g / r;
                g = d[i + 1] - p;
                r = (d[i] - g) * s + 2.0 * c * b;
                p = s * r;
                d[i + 1] = g + p;
                g = c * r - b;
            }

            if (r == 0.0) break;
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
        }

        if (l == 0) continue;
        for (int i = l; i > 0; --i) {
            if (d[i] >= d[i - 1]) break;
            std::swap(d[i], d[i - 1]);
        }
    }
}

static inline void prolfact(double* a, double* b, double* c, int n, double* u, double* v,
                            double* w) {
    for (int i = 0; i < n - 1; ++i) {
        double d = c[i + 1] / a[i];
        a[i + 1] -= b[i] * d;
        u[i] = d;
    }

    for (int i = n - 1; i > 0; --i) {
        double d = b[i - 1] / a[i];
        v[i] = d;
    }

    double done = 1.0;
    for (int i = 0; i < n; ++i) {
        w[i] = done / a[i];
    }
}

static inline void prolsolv(const double* u, const double* v, const double* w, int n, double* rhs) {
    for (int i = 0; i < n - 1; ++i) {
        rhs[i + 1] -= u[i] * rhs[i];
    }

    for (int i = n - 1; i > 0; --i) {
        rhs[i - 1] -= rhs[i] * v[i];
    }

    for (int i = 0; i < n; ++i) {
        rhs[i] *= w[i];
    }
}

static inline void prolfun0(int& ier, int n, double c, double* as, double* bs, double* cs,
                            double* xk, double* u, double* v, double* w, double eps, int& nterms,
                            double& rkhi) {
    ier = 0;
    double delta = 1.0e-8;
    int ifsymm = 1;
    int numit = 4;
    double rlam = 0;
    int ifodd = -1;

    prolmatr(as, bs, cs, n, c, rlam, ifsymm, ifodd);

    prolql1(n / 2, bs, as, ier);
    if (ier != 0) {
        ier = 2048;
        return;
    }

    rkhi = -bs[n / 2 - 1];
    rlam = -bs[n / 2 - 1] + delta;

    std::fill(xk, xk + n, 1.0);

    prolmatr(as, bs, cs, n, c, rlam, ifsymm, ifodd);

    prolfact(bs, cs, as, n / 2, u, v, w);

    for (int ijk = 0; ijk < numit; ++ijk) {
        prolsolv(u, v, w, n / 2, xk);

        double d = 0;
        for (int j = 0; j < n / 2; ++j) {
            d += xk[j] * xk[j];
        }

        d = std::sqrt(d);
        for (int j = 0; j < n / 2; ++j) {
            xk[j] /= d;
        }

        double err = 0;
        for (int j = 0; j < n / 2; ++j) {
            err += (as[j] - xk[j]) * (as[j] - xk[j]);
            as[j] = xk[j];
        }
        err = std::sqrt(err);
    }

    double half = 0.5;
    for (int i = 0; i < n / 2; ++i) {
        if (std::abs(xk[i]) > eps) nterms = i + 1;
        xk[i] *= std::sqrt(i * 2 + half);
        cs[i] = xk[i];
    }

    int j = 0;
    for (int i = 0; i <= nterms; ++i) {
        xk[j++] = cs[i];
        xk[j++] = 0;
    }

    nterms *= 2;
}

static inline void prolps0i(int& ier, double c, double* w, int lenw, int& nterms, int& ltot,
                            double& rkhi) {
    static const int ns[] = {48,  64,  80,  92,  106, 120, 130, 144, 156, 168,
                             178, 190, 202, 214, 224, 236, 248, 258, 268, 280};

    double eps = 1.0e-16;
    int n = static_cast<int>(c * 3);
    n = n / 2;

    int i = static_cast<int>(c / 10);
    if (i <= 19) n = ns[i];

    ier = 0;
    int ixk = 1;
    int lxk = n + 2;

    int ias = ixk + lxk;
    int las = n + 2;

    int ibs = ias + las;
    int lbs = n + 2;

    int ics = ibs + lbs;
    int lcs = n + 2;

    int iu = ics + lcs;
    int lu = n + 2;

    int iv = iu + lu;
    int lv = n + 2;

    int iw = iv + lv;
    int lw = n + 2;

    ltot = iw + lw;

    if (ltot >= lenw) {
        ier = 512;
        return;
    }

    prolfun0(ier, n, c, w + ias - 1, w + ibs - 1, w + ics - 1, w + ixk - 1, w + iu - 1, w + iv - 1,
             w + iw - 1, eps, nterms, rkhi);

    if (ier != 0) return;
}

// ============================================================================
// Prolate initialization and evaluation (literal from GROMACS pswf.cpp)
// ============================================================================

static inline void prol0ini(int& ier, double c, double* w, double& rlam20, double& rkhi, int lenw,
                            int& keep, int& ltot) {
    ier = 0;
    double thresh = 45;
    int iw = 11;
    w[0] = iw + 0.1;
    w[8] = thresh;

    int nterms = 0;
    prolps0i(ier, c, w + iw - 1, lenw, nterms, ltot, rkhi);

    if (ier != 0) return;

    if (c >= thresh) {
        w[7] = c;
        w[4] = nterms + 0.1;
        keep = nterms + 3;
        return;
    }

    int ngauss = nterms * 2;
    int lw = nterms + 2;
    int its = iw + lw;
    int lts = ngauss + 2;
    int iwhts = its + lts;
    int lwhts = ngauss + 2;
    int ifs = iwhts + lwhts;
    int lfs = ngauss + 2;

    keep = ifs + lfs;
    if (keep > ltot) ltot = keep;
    if (keep >= lenw) {
        ier = 1024;
        return;
    }

    w[1] = its + 0.1;
    w[2] = iwhts + 0.1;
    w[3] = ifs + 0.1;

    int itype = 1;
    std::vector<std::vector<double>> u, v;
    legeexps(itype, ngauss, w + its - 1, u, v, w + iwhts - 1);

    for (int i = 0; i < ngauss; ++i) {
        legeexev(w[its + i - 1], w[ifs + i - 1], w + iw - 1, nterms - 1);
    }

    double rlam = 0;
    double x0 = 0;
    double f0;
    legeexev(x0, f0, w + iw - 1, nterms - 1);
    double der;
    prosinin(c, w + its - 1, w + iwhts - 1, w + ifs - 1, x0, ngauss, rlam, der);

    rlam = rlam / f0;
    rlam20 = rlam;

    w[4] = nterms + 0.1;
    w[5] = ngauss + 0.1;
    w[6] = rlam;
    w[7] = c;
}

static inline void prol0eva(double x, const double* w, double& psi0, double& derpsi0) {
    int iw = static_cast<int>(w[0]);
    int its = static_cast<int>(w[1]);
    int iwhts = static_cast<int>(w[2]);
    int ifs = static_cast<int>(w[3]);

    int nterms = static_cast<int>(w[4]);
    int ngauss = static_cast<int>(w[5]);
    double rlam = w[6];
    double c = w[7];
    double thresh = w[8];

    if (std::abs(x) > 1) {
        if (c >= thresh - 1.0e-10) {
            psi0 = 0;
            derpsi0 = 0;
            return;
        }

        prosinin(c, &w[its - 1], &w[iwhts - 1], &w[ifs - 1], x, ngauss, psi0, derpsi0);
        psi0 /= rlam;
        derpsi0 /= rlam;
        return;
    }

    legeFDER(x, psi0, derpsi0, &w[iw - 1], nterms - 2);
}

static inline void prol0eva2(double x, const double* w, double& psi0, double& derpsi0, double& der2psi0) {
    int iw = static_cast<int>(w[0]);
    int its = static_cast<int>(w[1]);
    int iwhts = static_cast<int>(w[2]);
    int ifs = static_cast<int>(w[3]);

    int nterms = static_cast<int>(w[4]);
    int ngauss = static_cast<int>(w[5]);
    double rlam = w[6];
    double c = w[7];
    double thresh = w[8];

    if (std::abs(x) > 1) {
        if (c >= thresh - 1.0e-10) {
            psi0 = 0.0;
            derpsi0 = 0.0;
            der2psi0 = 0.0;
            return;
        }

        prosinin2(c, &w[its - 1], &w[iwhts - 1], &w[ifs - 1], x, ngauss, psi0, derpsi0, der2psi0);
        psi0 /= rlam;
        derpsi0 /= rlam;
        der2psi0 /= rlam;
        return;
    }

    legeFDER2(x, psi0, derpsi0, der2psi0, &w[iw - 1], nterms - 2);
}

static inline void prol0eva3(double x, const double* w, double& psi0, double& derpsi0, double& der2psi0, double& der3psi0) {
    int iw = static_cast<int>(w[0]);
    int its = static_cast<int>(w[1]);
    int iwhts = static_cast<int>(w[2]);
    int ifs = static_cast<int>(w[3]);

    int nterms = static_cast<int>(w[4]);
    int ngauss = static_cast<int>(w[5]);
    double rlam = w[6];
    double c = w[7];
    double thresh = w[8];

    if (std::abs(x) > 1) {
        if (c >= thresh - 1.0e-10) {
            psi0 = 0.0;
            derpsi0 = 0.0;
            der2psi0 = 0.0;
            der3psi0 = 0.0;
            return;
        }

        prosinin3(c, &w[its - 1], &w[iwhts - 1], &w[ifs - 1], x, ngauss, psi0, derpsi0, der2psi0, der3psi0);
        psi0 /= rlam;
        derpsi0 /= rlam;
        der2psi0 /= rlam;
        der3psi0 /= rlam;
        return;
    }

    legeFDER3(x, psi0, derpsi0, der2psi0, der3psi0, &w[iw - 1], nterms - 2);
}

static inline void prol0eva4(double x, const double* w, double& psi0, double& derpsi0, double& der2psi0, double& der3psi0, double& der4psi0) {
    int iw = static_cast<int>(w[0]);
    int its = static_cast<int>(w[1]);
    int iwhts = static_cast<int>(w[2]);
    int ifs = static_cast<int>(w[3]);

    int nterms = static_cast<int>(w[4]);
    int ngauss = static_cast<int>(w[5]);
    double rlam = w[6];
    double c = w[7];
    double thresh = w[8];

    if (std::abs(x) > 1) {
        if (c >= thresh - 1.0e-10) {
            psi0 = 0.0;
            derpsi0 = 0.0;
            der2psi0 = 0.0;
            der3psi0 = 0.0;
            der4psi0 = 0.0;
            return;
        }

        prosinin4(c, &w[its - 1], &w[iwhts - 1], &w[ifs - 1], x, ngauss, psi0, derpsi0, der2psi0, der3psi0, der4psi0);
        psi0 /= rlam;
        derpsi0 /= rlam;
        der2psi0 /= rlam;
        der3psi0 /= rlam;
        der4psi0 /= rlam;
        return;
    }

    legeFDER4(x, psi0, derpsi0, der2psi0, der3psi0, der4psi0, &w[iw - 1], nterms - 2);
}

static inline void prol0int0r(const double* w, double r, double& val) {
    int npts = 200;
    std::vector<double> xs(npts), ws(npts), fvals(npts);

    int itype = 1;
    std::vector<std::vector<double>> u, v;
    legeexps(itype, npts, xs.data(), u, v, ws.data());

    double derpsi0;
    for (int i = 0; i < npts; ++i) {
        double xs_r = (xs[i] + 1) * r / 2;
        prol0eva(xs_r, w, fvals[i], derpsi0);
    }

    val = 0;
    for (int i = 0; i < npts; ++i) {
        val += ws[i] * r / 2 * fvals[i];
    }
}

// Prolate0Fun: wraps the workarray-based prolate evaluation (literal from GROMACS)
struct Prolate0Fun {
    Prolate0Fun() : c(0), lenw(0), keep(0), ltot(0), rlam20(0), rkhi(0) {}

    explicit Prolate0Fun(double c_, int lenw_ = 10000) : c(c_), lenw(lenw_) {
        int ier;
        workarray.resize(lenw);
        prol0ini(ier, c, workarray.data(), rlam20, rkhi, lenw, keep, ltot);
        if (ier) throw std::runtime_error("Unable to init Prolate0Fun");
    }

    std::pair<double, double> eval_val_derivative(double x) const {
        double psi0, derpsi0;
        prol0eva(x, workarray.data(), psi0, derpsi0);
        return {psi0, derpsi0};
    }

    double eval_val(double x) const {
        std::pair<double, double> valDer = eval_val_derivative(x);
        return valDer.first;
    }

    double eval_derivative(double x) const {
        std::pair<double, double> valDer = eval_val_derivative(x);
        return valDer.second;
    }

    double eval_second_derivative(double x) const {
        double psi0, derpsi0, der2psi0;
        prol0eva2(x, workarray.data(), psi0, derpsi0, der2psi0);
        return der2psi0;
    }

    double eval_third_derivative(double x) const {
        double psi0, derpsi0, der2psi0, der3psi0;
        prol0eva3(x, workarray.data(), psi0, derpsi0, der2psi0, der3psi0);
        return der3psi0;
    }

    double eval_fourth_derivative(double x) const {
        double psi0, derpsi0, der2psi0, der3psi0, der4psi0;
        prol0eva4(x, workarray.data(), psi0, derpsi0, der2psi0, der3psi0, der4psi0);
        return der4psi0;
    }

    double int_eval(double r) const {
        double val;
        prol0int0r(workarray.data(), r, val);
        return val;
    }

    double c;
    int lenw, keep, ltot;
    std::vector<double> workarray;
    double rlam20, rkhi;
};

// ============================================================================
// Tolerance-to-c lookup table (literal from GROMACS prolc180)
// ============================================================================

static inline void prolc180(double eps, double& c) {
    static const double cs[] = {
        0.43368E-16, 0.10048E+01, 0.17298E+01, 0.22271E+01, 0.26382E+01, 0.30035E+01, 0.33409E+01,
        0.36598E+01, 0.39658E+01, 0.42621E+01, 0.45513E+01, 0.48347E+01, 0.51136E+01, 0.53887E+01,
        0.56606E+01, 0.59299E+01, 0.61968E+01, 0.64616E+01, 0.67247E+01, 0.69862E+01, 0.72462E+01,
        0.75049E+01, 0.77625E+01, 0.80189E+01, 0.82744E+01, 0.85289E+01, 0.87826E+01, 0.90355E+01,
        0.92877E+01, 0.95392E+01, 0.97900E+01, 0.10040E+02, 0.10290E+02, 0.10539E+02, 0.10788E+02,
        0.11036E+02, 0.11284E+02, 0.11531E+02, 0.11778E+02, 0.12024E+02, 0.12270E+02, 0.12516E+02,
        0.12762E+02, 0.13007E+02, 0.13251E+02, 0.13496E+02, 0.13740E+02, 0.13984E+02, 0.14228E+02,
        0.14471E+02, 0.14714E+02, 0.14957E+02, 0.15200E+02, 0.15443E+02, 0.15685E+02, 0.15927E+02,
        0.16169E+02, 0.16411E+02, 0.16652E+02, 0.16894E+02, 0.17135E+02, 0.17376E+02, 0.17617E+02,
        0.17858E+02, 0.18098E+02, 0.18339E+02, 0.18579E+02, 0.18819E+02, 0.19059E+02, 0.19299E+02,
        0.19539E+02, 0.19778E+02, 0.20018E+02, 0.20257E+02, 0.20496E+02, 0.20736E+02, 0.20975E+02,
        0.21214E+02, 0.21452E+02, 0.21691E+02, 0.21930E+02, 0.22168E+02, 0.22407E+02, 0.22645E+02,
        0.22884E+02, 0.23122E+02, 0.23360E+02, 0.23598E+02, 0.23836E+02, 0.24074E+02, 0.24311E+02,
        0.24549E+02, 0.24787E+02, 0.25024E+02, 0.25262E+02, 0.25499E+02, 0.25737E+02, 0.25974E+02,
        0.26211E+02, 0.26448E+02, 0.26685E+02, 0.26922E+02, 0.27159E+02, 0.27396E+02, 0.27633E+02,
        0.27870E+02, 0.28106E+02, 0.28343E+02, 0.28580E+02, 0.28816E+02, 0.29053E+02, 0.29289E+02,
        0.29526E+02, 0.29762E+02, 0.29998E+02, 0.30234E+02, 0.30471E+02, 0.30707E+02, 0.30943E+02,
        0.31179E+02, 0.31415E+02, 0.31651E+02, 0.31887E+02, 0.32123E+02, 0.32358E+02, 0.32594E+02,
        0.32830E+02, 0.33066E+02, 0.33301E+02, 0.33537E+02, 0.33773E+02, 0.34008E+02, 0.34244E+02,
        0.34479E+02, 0.34714E+02, 0.34950E+02, 0.35185E+02, 0.35421E+02, 0.35656E+02, 0.35891E+02,
        0.36126E+02, 0.36362E+02, 0.36597E+02, 0.36832E+02, 0.37067E+02, 0.37302E+02, 0.37537E+02,
        0.37772E+02, 0.38007E+02, 0.38242E+02, 0.38477E+02, 0.38712E+02, 0.38947E+02, 0.39181E+02,
        0.39416E+02, 0.39651E+02, 0.39886E+02, 0.40120E+02, 0.40355E+02, 0.40590E+02, 0.40824E+02,
        0.41059E+02, 0.41294E+02, 0.41528E+02, 0.41763E+02, 0.41997E+02, 0.42232E+02, 0.42466E+02,
        0.42700E+02, 0.42935E+02, 0.43169E+02, 0.43404E+02, 0.43638E+02, 0.43872E+02, 0.44107E+02,
        0.44341E+02, 0.44575E+02, 0.44809E+02, 0.45044E+02, 0.45278E+02};

    double e = eps;
    if (e < 1.0e-18) e = 1e-18;
    double d = -log10(e);
    int i = static_cast<int>(d * 10 + 0.1);
    c = cs[i - 1];
}

static inline double get_prolate_c(double tol) {
    double c;
    prolc180(tol, c);
    return c;
}

// ============================================================================
// Chebyshev nodes and interpolation (literal from GROMACS pswf.cpp)
// ============================================================================

static const int MAX_CHEB_ORDER = 30;

static inline void cheb_nodes_1d(int order, std::vector<double>& nodes, double a = 0, double b = 1) {
    const double pi = 3.1415926535897932384626433832795028841;
    nodes.resize(order);
    for (int i = 0; i < order; i++) {
        nodes[i] = -cos((i + 0.5) * pi / order) * 0.5 + 0.5;
        nodes[i] = nodes[i] * (b - a) + a;
    }
}

static inline void cheb_basis_1d(int order, const std::vector<double>& x, std::vector<double>& y,
                                  double a = 0, double b = 1) {
    int n = (int)x.size();
    y.resize(order * n);

    if (order > 0) {
        for (int i = 0; i < n; i++) {
            y[i] = 1.0;
        }
    }
    if (order > 1) {
        for (int i = 0; i < n; i++) {
            y[i + n] = x[i] * 2 / (b - a) - 2 * a / (b - a) - 1;
        }
    }
    for (int i = 2; i < order; i++) {
        for (int j = 0; j < n; j++) {
            y[i * n + j] = 2 * y[n + j] * y[i * n - n + j] - y[i * n - 2 * n + j];
        }
    }
}

static inline void cheb_interp_1d(int order, std::vector<double>& fn_v,
                                    std::vector<double>& coeff) {
    std::vector<double> x, p;
    cheb_nodes_1d(order, x);
    cheb_basis_1d(order, x, p);

    const size_t dof = fn_v.size() / order;
    assert(fn_v.size() == dof * (size_t)order);
    coeff.resize(dof * (size_t)order);

    const double invOrder = 1.0 / (double)order;
    const double twoInvOrder = 2.0 * invOrder;
    for (size_t idof = 0; idof < dof; ++idof) {
        const size_t offset = idof * (size_t)order;
        for (int k = 0; k < order; ++k) {
            double sum = 0.0;
            const double* pk = &p[(size_t)k * (size_t)order];
            for (int j = 0; j < order; ++j) {
                sum += fn_v[offset + j] * pk[j];
            }
            coeff[offset + (size_t)k] = (k == 0) ? (sum * invOrder) : (sum * twoInvOrder);
        }
    }
}

// ============================================================================
// Monomial interpolation via Newton divided differences (literal from GROMACS)
// ============================================================================

static inline void monomial_interp_1d(int order, int nnodes, std::vector<double>& fn_v,
                                       std::vector<double>& coeff,
                                       double a = 0, double b = 1) {
    assert(order == nnodes);

    std::vector<double> x;
    cheb_nodes_1d(nnodes, x, a, b);

    auto multiply_x = [](const std::vector<double>& p, double x_in) {
        std::vector<double> r(p.size() + 1, 0.0);
        for (size_t i = 0; i < p.size(); ++i) {
            r[i] += -x_in * p[i];
            r[i + 1] += p[i];
        }
        return r;
    };

    const size_t dof = fn_v.size() / nnodes;
    assert(fn_v.size() == dof * (size_t)nnodes);

    std::vector<double> newton_coeffs = fn_v;
    coeff.assign(dof * (size_t)order, 0.0);

    for (size_t idof = 0; idof < dof; ++idof) {
        const size_t fn_offset = idof * (size_t)nnodes;
        const size_t coeff_offset = idof * (size_t)order;

        for (int j = 1; j < nnodes; ++j) {
            for (int i = nnodes - 1; i >= j; --i) {
                newton_coeffs[fn_offset + i] =
                    (newton_coeffs[fn_offset + i] - newton_coeffs[fn_offset + i - 1]) / (x[i] - x[i - j]);
            }
        }

        coeff[coeff_offset + (size_t)(order - 1)] = newton_coeffs[fn_offset];
        std::vector<double> basis{1.0};
        for (int j = 1; j < nnodes; ++j) {
            basis = multiply_x(basis, x[j - 1]);
            const double newton_coeff = newton_coeffs[fn_offset + j];
            for (size_t m = 0; m < basis.size(); ++m) {
                coeff[coeff_offset + (size_t)(order - 1) - m] += newton_coeff * basis[m];
            }
        }
    }
}

// ============================================================================
// Prolate0Data: public API wrapping Prolate0Fun
// ============================================================================

struct Prolate0Data {
    Prolate0Fun pfun;

    double eval_val(double x) const { return pfun.eval_val(x); }
    double eval_derivative(double x) const { return pfun.eval_derivative(x); }
    double eval_second_derivative(double x) const { return pfun.eval_second_derivative(x); }
    double eval_third_derivative(double x) const { return pfun.eval_third_derivative(x); }
    double eval_fourth_derivative(double x) const { return pfun.eval_fourth_derivative(x); }

    void eval(double x, double& val, double& der) const {
        std::pair<double, double> valDer = pfun.eval_val_derivative(x);
        val = valDer.first;
        der = valDer.second;
    }
    double integral(double r) const { return pfun.int_eval(r); }
};

static inline Prolate0Data prolate0_prepare(double c) {
    Prolate0Data pdata;
    pdata.pfun = Prolate0Fun(c);
    return pdata;
}

// ============================================================================
// EspCoefficients and builder
// ============================================================================

enum SpreadDerivativeFitMode {
    DifferentiateSpreadPolynomial = 0,
    DirectFiniteDifferenceDerivativePolynomial = 1,
    DirectAnalyticDerivativePolynomial = 2,
    DirectAnalyticDerivativePolynomialSeparateOrder = 3
};

struct EspCoefficients {
    double c_s;         // Splitting bandwidth (larger, controls Fourier truncation)
    double c_w;         // Window bandwidth (smaller, controls aliasing; c_w = pi/2 * P)
    double c;           // Alias for c_s (backwards compatibility)
    double c0;          // Integral: int_0^1 psi_s(x) dx (using c_s)
    double psi0;        // psi_s(0) (using c_s)
    double lambda;      // Fourier eigenvalue (using c_s)
    int P;              // Stencil size

    // Spreading uses c_w (monomial, descending Horner order)
    int spread_poly_order;
    std::vector<float> spread_coeffs;       // [P * spread_poly_order]
    std::vector<float> spread_der_coeffs;   // [P * spread_der_poly_order]
    int spread_der_poly_order;
    std::vector<float> spread_der2_coeffs;  // [P * spread_der2_poly_order]
    int spread_der2_poly_order;
    std::vector<float> spread_der3_coeffs;  // [P * spread_der3_poly_order]
    int spread_der3_poly_order;

    // Fourier splitting uses c_s (monomial, descending Horner order)
    int split_fourier_poly_order;
    std::vector<float> split_fourier_coeffs;     // [split_fourier_poly_order]
    std::vector<float> split_fourier_der_coeffs; // [split_fourier_der_poly_order]
    int split_fourier_der_poly_order;

    // Fourier splitting Chebyshev coefficients (ascending order, on [0,1])
    int split_fourier_cheb_order;
    std::vector<double> split_fourier_cheb_coeffs;  // [split_fourier_cheb_order]

    // Window function properties (using c_w)
    double c0_w;        // Integral: int_0^1 psi_w(x) dx
    double psi0_w;      // psi_w(0)
    double lambda_w;    // Fourier eigenvalue for window

    // Self-energy: psi_s(0) / c0_s (uses splitting PSWF)
    double self_energy_factor() const { return psi0 / c0; }
    double dipole_far_field_A0 = 0.0;
    double quadrupole_far_field_B0 = 0.0;

    double dipole_self_field_scale() const { return -dipole_far_field_A0; }
    double dipole_self_energy_scale() const { return 0.5*dipole_far_field_A0; }
    double quadrupole_self_energy_scale() const { return -quadrupole_far_field_B0; }
};

// Helper: compute lambda (Fourier eigenvalue) for a given prolate function
static inline double compute_lambda(const Prolate0Fun& pfun, double c) {
    int quad_npts = 200;
    std::vector<double> xs(quad_npts), ws(quad_npts);
    legerts(quad_npts, xs.data(), ws.data());
    double lam = 0.0;
    for (int i = 0; i < quad_npts; i++) {
        lam += ws[i] * pfun.eval_val(xs[i]) * std::cos(c * xs[i] * 0.5);
    }
    lam /= pfun.eval_val(0.5);
    return lam;
}

static inline double finite_difference_second_derivative(const Prolate0Fun& pfun, double x, double h = 1.0e-5) {
    const double fp = pfun.eval_val(x+h);
    const double f0 = pfun.eval_val(x);
    const double fm = pfun.eval_val(x-h);
    return (fp - 2.0*f0 + fm)/(h*h);
}

static inline double finite_difference_third_derivative(const Prolate0Fun& pfun, double x, double h = 1.0e-5) {
    const double fpp = pfun.eval_val(x+2.0*h);
    const double fp = pfun.eval_val(x+h);
    const double fm = pfun.eval_val(x-h);
    const double fmm = pfun.eval_val(x-2.0*h);
    return (fpp - 2.0*fp + 2.0*fm - fmm)/(2.0*h*h*h);
}

static inline double finite_difference_fourth_derivative_even0(const Prolate0Fun& pfun, double h = 1.0e-3) {
    return (pfun.eval_val(-2.0*h)-4.0*pfun.eval_val(-h)+6.0*pfun.eval_val(0.0)-4.0*pfun.eval_val(h)+pfun.eval_val(2.0*h))/(h*h*h*h);
}

// Build ESP polynomial coefficients with separate splitting (c_s) and window (c_w) parameters.
// Per Bostrom, Tornberg, af Klinteberg (arXiv:2602.16591):
//   c_s controls the Ewald split (Fourier truncation error ~ e^{-c_s})
//   c_w controls the window/spreading (aliasing error ~ e^{-c_w}), c_w = pi/2 * P
// Splitting coefficients (split_fourier, self-energy) use c_s.
// Spreading coefficients (spread_coeffs, deconvolution) use c_w.
static inline EspCoefficients build_esp_coefficients_split(double c_s, double c_w, int P,
                                                            double tol,
                                                            int max_poly_order = 16,
                                                            int min_spread_poly_order = -1,
                                                            int min_spread_der_poly_order = -1,
                                                            int min_split_fourier_cheb_order = -1,
                                                            SpreadDerivativeFitMode spread_derivative_fit_mode = DifferentiateSpreadPolynomial,
                                                            int max_spread_poly_order = -1,
                                                            int max_spread_der_poly_order = -1) {
    EspCoefficients esp;
    esp.c_s = c_s;
    esp.c_w = c_w;
    esp.c = c_s;  // backwards compatibility
    esp.P = P;

    // --- Splitting PSWF (c_s) for Fourier kernel and self-energy ---
    Prolate0Fun pfun_s(c_s);
    esp.c0 = pfun_s.int_eval(1.0);
    esp.psi0 = pfun_s.eval_val(0.0);
    esp.lambda = compute_lambda(pfun_s, c_s);
    {
        const double psi2 = pfun_s.eval_second_derivative(0.0);
        const double psi4 = pfun_s.eval_fourth_derivative(0.0);
        esp.dipole_far_field_A0 = psi2/(3.0*esp.c0);
        esp.quadrupole_far_field_B0 = psi4/(15.0*esp.c0);
    }

    // --- Window PSWF (c_w) for spreading/interpolation ---
    Prolate0Fun pfun_w(c_w);
    esp.c0_w = pfun_w.int_eval(1.0);
    esp.psi0_w = pfun_w.eval_val(0.0);
    esp.lambda_w = compute_lambda(pfun_w, c_w);

    // Polynomial order selection: Chebyshev coefficient decay filtering.
    // The filter tolerance (tol_coeff) determines how many polynomial terms
    // are kept: coefficients below tol_coeff * max_coeff are truncated.
    //
    // Spread polynomials use a 10x looser tolerance than the ESP tolerance.
    // Benchmarks across tol=1e-3..1e-6 show this saves 1 polynomial term
    // while keeping force L2 relative error well within the ESP tolerance:
    //   tol=1e-3: spread 6->5, force_rel=8.1e-4 < 1e-3
    //   tol=1e-4: spread 7->6, force_rel=6.1e-5 < 1e-4
    //   tol=1e-5: spread 8->7, force_rel=3.9e-6 < 1e-5
    //   tol=1e-6: spread 9->8, force_rel=3.1e-7 < 1e-6
    //
    // Fourier polynomials are evaluated once per box change (0.3% of runtime),
    // so they keep the full tolerance — no speed benefit from loosening.
    const double tol_coeff_spread  = tol;          // spread window & derivative
    const double tol_coeff_fourier = tol;          // fourier splitting & derivative

    // --- Spreading coefficients use c_w (window PSWF) ---
    {
        int order = MAX_CHEB_ORDER;
        std::vector<double> nodes;
        cheb_nodes_1d(order, nodes);

        int dof = P;
        std::vector<double> fn_v(dof * order);
        for (int idof = 0; idof < dof; idof++) {
            for (int i = 0; i < order; i++) {
                double arg = nodes[i] - P / 2.0 + (dof - idof - 1);
                arg /= P / 2.0;
                fn_v[idof * order + i] = pfun_w.eval_val(arg);
            }
        }

        // Chebyshev interpolation for order estimation
        std::vector<double> cheb_coeff;
        cheb_interp_1d(order, fn_v, cheb_coeff);

        int est_order = -1;
        for (int idof = 0; idof < dof; idof++) {
            double max_c = 0.0;
            for (int i = 0; i < order; i++)
                max_c = std::max(max_c, std::abs(cheb_coeff[idof * order + i]));
            for (int i = 0; i < order; i++) {
                if (std::abs(cheb_coeff[idof * order + i]) > tol_coeff_spread * max_c)
                    est_order = std::max(est_order, i + 1);
            }
        }
        if (min_spread_poly_order > 0)
            est_order = std::max(est_order, min_spread_poly_order);
        if (max_spread_poly_order > 0)
            est_order = std::min(est_order, max_spread_poly_order);
        if (est_order > max_poly_order) est_order = max_poly_order;

        // Build monomial coefficients at estimated order
        int nnodes = est_order;
        cheb_nodes_1d(nnodes, nodes, 0, 1);
        fn_v.resize(dof * nnodes);
        for (int idof = 0; idof < dof; idof++) {
            for (int i = 0; i < nnodes; i++) {
                double arg = nodes[i] - P / 2.0 + (dof - idof - 1);
                arg /= P / 2.0;
                fn_v[idof * nnodes + i] = pfun_w.eval_val(arg);
            }
        }

        std::vector<double> coeffs_tmp;
        monomial_interp_1d(est_order, nnodes, fn_v, coeffs_tmp);

        // Store in dense layout [P * est_order], direct order (no storage reversal).
        // The fn_v building already uses (dof-idof-1) to reverse the stencil order
        // (matching GROMACS convention). An additional reversal here would cancel it.
        // scoeffs[ip] for grid point gridIndex+ip must evaluate
        // spread_window_ref(c, P, P-1-ip, dr) = psi((dr + P/2 - 1 - ip)/(P/2)).
        esp.spread_poly_order = est_order;
        esp.spread_coeffs.resize(P * est_order);
        for (int ip = 0; ip < P; ++ip) {
            int idof = ip;
            for (int j = 0; j < est_order; ++j) {
                esp.spread_coeffs[ip * est_order + j] = (float)coeffs_tmp[idof * est_order + j];
            }
        }
    }

    // --- Spreading derivative coefficients ---
    // Default path: analytically differentiate the same monomial spread
    // polynomials to preserve discrete adjoint consistency between source and
    // gather operators.  Optional A/B path: independently fit each derivative
    // level while keeping the same polynomial orders, to isolate whether the
    // residual error is dominated by derivative approximation or operator
    // mismatch.
    if (spread_derivative_fit_mode == DifferentiateSpreadPolynomial) {
        const int derOrder = std::max(1, esp.spread_poly_order-1);
        esp.spread_der_poly_order = derOrder;
        if (max_spread_der_poly_order > 0)
            esp.spread_der_poly_order = std::min(esp.spread_der_poly_order, max_spread_der_poly_order);
        esp.spread_der_coeffs.assign(P * esp.spread_der_poly_order, 0.0f);
        for (int ip = 0; ip < P; ++ip) {
            const float* src = &esp.spread_coeffs[ip * esp.spread_poly_order];
            float* dst = &esp.spread_der_coeffs[ip * esp.spread_der_poly_order];
            const int degree = esp.spread_poly_order-1;
            for (int j = 0; j < esp.spread_der_poly_order; ++j)
                dst[j] = (float) ((degree-j) * src[j]);
        }

        const int der2Order = std::max(1, esp.spread_der_poly_order-1);
        esp.spread_der2_poly_order = der2Order;
        esp.spread_der2_coeffs.assign(P * der2Order, 0.0f);
        for (int ip = 0; ip < P; ++ip) {
            const float* src = &esp.spread_der_coeffs[ip * esp.spread_der_poly_order];
            float* dst = &esp.spread_der2_coeffs[ip * der2Order];
            const int degree = esp.spread_der_poly_order-1;
            for (int j = 0; j < der2Order; ++j)
                dst[j] = (float) ((degree-j) * src[j]);
        }

        const int der3Order = std::max(1, esp.spread_der2_poly_order-1);
        esp.spread_der3_poly_order = der3Order;
        esp.spread_der3_coeffs.assign(P * der3Order, 0.0f);
        for (int ip = 0; ip < P; ++ip) {
            const float* src = &esp.spread_der2_coeffs[ip * esp.spread_der2_poly_order];
            float* dst = &esp.spread_der3_coeffs[ip * der3Order];
            const int degree = esp.spread_der2_poly_order-1;
            for (int j = 0; j < der3Order; ++j)
                dst[j] = (float) ((degree-j) * src[j]);
        }
    }
    else {
        const bool useAnalyticDerivatives = (spread_derivative_fit_mode == DirectAnalyticDerivativePolynomial ||
                                             spread_derivative_fit_mode == DirectAnalyticDerivativePolynomialSeparateOrder);
        const bool useSeparateDerivativeOrders = (spread_derivative_fit_mode == DirectAnalyticDerivativePolynomialSeparateOrder);

        auto evalSpreadDerivativeValue = [&](int derivativeOrder, double arg) {
            const double scale = 2.0/P;
            const double scalePow = std::pow(scale, derivativeOrder);
            if (derivativeOrder == 1)
                return scalePow*pfun_w.eval_derivative(arg);
            if (derivativeOrder == 2)
                return scalePow*(useAnalyticDerivatives ? pfun_w.eval_second_derivative(arg) : finite_difference_second_derivative(pfun_w, arg));
            return scalePow*(useAnalyticDerivatives ? pfun_w.eval_third_derivative(arg) : finite_difference_third_derivative(pfun_w, arg));
        };

        auto estimateSpreadDerivativePolyOrder = [&](int derivativeOrder, int minOrder) {
            const int order = MAX_CHEB_ORDER;
            std::vector<double> nodes;
            cheb_nodes_1d(order, nodes);
            std::vector<double> fn_v(P * order);
            for (int idof = 0; idof < P; ++idof) {
                for (int i = 0; i < order; ++i) {
                    double arg = nodes[i] - P / 2.0 + (P - idof - 1);
                    arg /= P / 2.0;
                    fn_v[idof * order + i] = evalSpreadDerivativeValue(derivativeOrder, arg);
                }
            }

            std::vector<double> cheb_coeff;
            cheb_interp_1d(order, fn_v, cheb_coeff);

            int est_order = 1;
            for (int idof = 0; idof < P; ++idof) {
                double max_c = 0.0;
                for (int i = 0; i < order; ++i)
                    max_c = std::max(max_c, std::abs(cheb_coeff[idof * order + i]));
                for (int i = 0; i < order; ++i) {
                    if (std::abs(cheb_coeff[idof * order + i]) > tol_coeff_spread * max_c)
                        est_order = std::max(est_order, i + 1);
                }
            }
            if (minOrder > 0)
                est_order = std::max(est_order, minOrder);
            if (max_spread_der_poly_order > 0)
                est_order = std::min(est_order, max_spread_der_poly_order);
            if (est_order > max_poly_order)
                est_order = max_poly_order;
            return std::max(est_order, 1);
        };

        auto fitSpreadDerivative = [&](int derivativeOrder, int polyOrder, std::vector<float>& out) {
            std::vector<double> nodes;
            cheb_nodes_1d(polyOrder, nodes, 0, 1);
            std::vector<double> fn_v(P * polyOrder);
            for (int idof = 0; idof < P; ++idof) {
                for (int i = 0; i < polyOrder; ++i) {
                    double arg = nodes[i] - P / 2.0 + (P - idof - 1);
                    arg /= P / 2.0;
                    fn_v[idof * polyOrder + i] = evalSpreadDerivativeValue(derivativeOrder, arg);
                }
            }
            std::vector<double> coeffs_tmp;
            monomial_interp_1d(polyOrder, polyOrder, fn_v, coeffs_tmp);
            out.resize(P * polyOrder);
            for (int ip = 0; ip < P; ++ip)
                for (int j = 0; j < polyOrder; ++j)
                    out[ip * polyOrder + j] = (float) coeffs_tmp[ip * polyOrder + j];
        };

        if (useSeparateDerivativeOrders) {
            esp.spread_der_poly_order = estimateSpreadDerivativePolyOrder(1, min_spread_der_poly_order > 0 ? min_spread_der_poly_order : 1);
            esp.spread_der2_poly_order = estimateSpreadDerivativePolyOrder(2, 1);
            esp.spread_der3_poly_order = estimateSpreadDerivativePolyOrder(3, 1);
        }
        else {
            esp.spread_der_poly_order = std::max(1, esp.spread_poly_order-1);
            if (min_spread_der_poly_order > 0)
                esp.spread_der_poly_order = std::max(esp.spread_der_poly_order, min_spread_der_poly_order);
            esp.spread_der2_poly_order = std::max(1, esp.spread_der_poly_order-1);
            esp.spread_der3_poly_order = std::max(1, esp.spread_der2_poly_order-1);
        }

        fitSpreadDerivative(1, esp.spread_der_poly_order, esp.spread_der_coeffs);
        fitSpreadDerivative(2, esp.spread_der2_poly_order, esp.spread_der2_coeffs);
        fitSpreadDerivative(3, esp.spread_der3_poly_order, esp.spread_der3_coeffs);
    }

    // --- Fourier splitting coefficients use c_s (splitting PSWF) ---
    // split_fourier(x) = lambda_s * psi_s(x) / c0_s for x in [0, 1]
    {
        int order = MAX_CHEB_ORDER;
        std::vector<double> nodes;
        cheb_nodes_1d(order, nodes);

        std::vector<double> fn_v(order);
        for (int i = 0; i < order; i++) {
            fn_v[i] = esp.lambda * pfun_s.eval_val(nodes[i]) / esp.c0;
        }

        std::vector<double> cheb_coeff;
        cheb_interp_1d(order, fn_v, cheb_coeff);

        int est_order = -1;
        double max_c = 0.0;
        for (int i = 0; i < order; i++)
            max_c = std::max(max_c, std::abs(cheb_coeff[i]));
        for (int i = 0; i < order; i++) {
            if (std::abs(cheb_coeff[i]) > tol_coeff_fourier * max_c)
                est_order = std::max(est_order, i + 1);
        }
        if (min_split_fourier_cheb_order > 0)
            est_order = std::max(est_order, min_split_fourier_cheb_order);
        if (est_order > max_poly_order) est_order = max_poly_order;

        int nnodes = est_order;
        cheb_nodes_1d(nnodes, nodes, 0, 1);
        fn_v.resize(nnodes);
        for (int i = 0; i < nnodes; i++) {
            fn_v[i] = esp.lambda * pfun_s.eval_val(nodes[i]) / esp.c0;
        }

        std::vector<double> coeffs_tmp;
        monomial_interp_1d(est_order, nnodes, fn_v, coeffs_tmp);

        esp.split_fourier_poly_order = est_order;
        esp.split_fourier_coeffs.resize(est_order);
        for (int i = 0; i < est_order; ++i) {
            esp.split_fourier_coeffs[i] = (float)coeffs_tmp[i];
        }

        // Also store Chebyshev coefficients at the same order for Clenshaw evaluation
        // Re-evaluate at est_order Chebyshev nodes and compute Chebyshev expansion
        cheb_nodes_1d(est_order, nodes);
        fn_v.resize(est_order);
        for (int i = 0; i < est_order; i++) {
            fn_v[i] = esp.lambda * pfun_s.eval_val(nodes[i]) / esp.c0;
        }
        cheb_interp_1d(est_order, fn_v, cheb_coeff);
        esp.split_fourier_cheb_order = est_order;
        esp.split_fourier_cheb_coeffs.resize(est_order);
        for (int i = 0; i < est_order; ++i) {
            esp.split_fourier_cheb_coeffs[i] = cheb_coeff[i];
        }
    }

    // --- Fourier splitting derivative coefficients ---
    {
        int order = MAX_CHEB_ORDER;
        std::vector<double> nodes;
        cheb_nodes_1d(order, nodes);

        std::vector<double> fn_v(order);
        for (int i = 0; i < order; i++) {
            fn_v[i] = esp.lambda * nodes[i] * pfun_s.eval_derivative(nodes[i]) / esp.c0;
        }

        std::vector<double> cheb_coeff;
        cheb_interp_1d(order, fn_v, cheb_coeff);

        int est_order = -1;
        double max_c = 0.0;
        for (int i = 0; i < order; i++)
            max_c = std::max(max_c, std::abs(cheb_coeff[i]));
        for (int i = 0; i < order; i++) {
            if (std::abs(cheb_coeff[i]) > tol_coeff_fourier * max_c)
                est_order = std::max(est_order, i + 1);
        }
        if (est_order > max_poly_order) est_order = max_poly_order;

        int nnodes = est_order;
        cheb_nodes_1d(nnodes, nodes, 0, 1);
        fn_v.resize(nnodes);
        for (int i = 0; i < nnodes; i++) {
            fn_v[i] = esp.lambda * nodes[i] * pfun_s.eval_derivative(nodes[i]) / esp.c0;
        }

        std::vector<double> coeffs_tmp;
        monomial_interp_1d(est_order, nnodes, fn_v, coeffs_tmp);

        esp.split_fourier_der_poly_order = est_order;
        esp.split_fourier_der_coeffs.resize(est_order);
        for (int i = 0; i < est_order; ++i) {
            esp.split_fourier_der_coeffs[i] = (float)coeffs_tmp[i];
        }
    }

    return esp;
}

// Tolerance-based with explicit P override.
static inline EspCoefficients build_esp_coefficients_with_P(double tol, int P,
                                                             int max_poly_order = 16,
                                                             int min_spread_poly_order = -1,
                                                             int min_spread_der_poly_order = -1,
                                                             int min_split_fourier_cheb_order = -1,
                                                             SpreadDerivativeFitMode spread_derivative_fit_mode = DifferentiateSpreadPolynomial,
                                                             int max_spread_poly_order = -1,
                                                             int max_spread_der_poly_order = -1) {
    double c = get_prolate_c(tol);
    return build_esp_coefficients_split(c, c, P, tol, max_poly_order,
                                        min_spread_poly_order, min_spread_der_poly_order,
                                        min_split_fourier_cheb_order,
                                        spread_derivative_fit_mode,
                                        max_spread_poly_order, max_spread_der_poly_order);
}

// ============================================================================
// PSWF moduli for Fourier deconvolution
// ============================================================================

// Inverse PSWF moduli for Fourier deconvolution (multiply instead of divide).
// Modulus = (h * lambda * psi(zarg*k))^2 where h = P/2.
// NOTE: This uses the continuous FT value (single alias n=0 only).
// For typical ESP parameters (P>=16), scale*gridSize = pi*P/c > 2,
// meaning no aliases fall within PSWF support [-1,1], so the continuous
// approximation is exact. For small P or large c where scale*gridSize < 2,
// an aliased sum would improve accuracy.
static inline std::vector<float> compute_pswf_inv_moduli(int n, int order,
                                                          const Prolate0Data& pdata,
                                                          double c, double lambda) {
    std::vector<float> inv_mod(n);
    double scale = M_PI * order / ((double)n * c);
    double h = order / 2.0;  // P/2 factor: DFT of spread window = h * lambda * psi(zarg*k)
    int maxk = (n + 1) / 2;

    // Threshold: if f² is below this, treat as zero to avoid huge inv_mod values.
    // The k=0 value (h * lambda * psi(0))² is the largest; use it as reference.
    double f0 = h * lambda * pdata.eval_val(0.0);
    double f0_sq = f0 * f0;
    double threshold = f0_sq * 1e-20;  // 10 orders of magnitude below peak

    // k = 0
    inv_mod[0] = (f0_sq > threshold) ? (float)(1.0 / f0_sq) : 0.0f;

    // Positive k
    for (int k = 1; k < maxk; ++k) {
        double arg = scale * k;
        if (arg > 1.0) {
            inv_mod[k] = 0.0f;  // Out of support: multiply by 0
        } else {
            double f = h * lambda * pdata.eval_val(arg);
            double f2 = f * f;
            inv_mod[k] = (f2 > threshold) ? (float)(1.0 / f2) : 0.0f;
        }
    }

    // Negative k (mapped to n+k for k < 0)
    for (int k = -maxk; k < 0; ++k) {
        double arg = -scale * k;
        if (arg > 1.0) {
            inv_mod[n + k] = 0.0f;
        } else {
            double f = h * lambda * pdata.eval_val(arg);
            double f2 = f * f;
            inv_mod[n + k] = (f2 > threshold) ? (float)(1.0 / f2) : 0.0f;
        }
    }

    return inv_mod;
}
} // namespace pswf

#endif // OPENMM_ESP_PSWF_H_
