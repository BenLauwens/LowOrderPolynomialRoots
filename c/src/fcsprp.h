#include <stdlib.h>
#include <math.h>
#include <float.h>

static inline double quadratic(double const coeffs[3]) {
    double a1 = 1.0 / coeffs[0];
    double b = -0.5 * coeffs[1] * a1;
    double c = coeffs[2] * a1;
    double delta = b * b - c;
    double res = DBL_MAX;
    if (delta < -4 * DBL_EPSILON) {
    } else if (delta > 4 * DBL_EPSILON) {
        if (b > DBL_EPSILON) {
            if (c > DBL_EPSILON) {
                res = c / (b + sqrt(delta));
            } else {
                res = b + sqrt(delta);
            } 
        } else if (b < -DBL_EPSILON) {
            if (c < DBL_EPSILON) {
                res = c / (b - sqrt(delta));
            }
        } else {
                res = sqrt(-c);
        }
    } else {
        if (b > -DBL_EPSILON) {
            res = b;
        }
    }
    return res;
}

static inline void horner(size_t const N, double const coeffs[N], double x, double * const v) {
    double vv = coeffs[0];
    for (size_t i = 1; i < N; ++i) {
        vv = fma(vv, x, coeffs[i]); //vv * x + coeffs[i];
    }
    *v = vv;
}

static inline void horner2(size_t const N, double const coeffs[N], double const x, double const y, double * const v, double * const w) {
    double vv = coeffs[0];
    double ww = coeffs[0];
    for (size_t i = 1; i < N; ++i) {
        vv = fma(vv, x, coeffs[i]); //vv * x + coeffs[i];
        ww = fma(ww, y, coeffs[i]); //ww * y + coeffs[i];
    }
    *v = vv;
    *w = ww;
}

static inline void hornerd(size_t const N, double const coeffs[N], double const x, double * const v, double * const d) {
    double vv = coeffs[0];
    double dd = 0.0;
    for (size_t i = 1; i < N; ++i) {
        dd = fma(dd, x, vv); //dd * x + vv;
        vv = fma(vv, x, coeffs[i]); //vv * x + coeffs[i];
    }
    *v = vv;
    *d = dd;
}

static inline void posintervalhorner(size_t const N, double const coeffs[N], double const low, double const high, double * const colow, double * const cohigh) {
    double col = coeffs[0];
    double coh = coeffs[0];
    double tmp;
    for (size_t i = 1; i < N; ++i) {
        if (col > -DBL_EPSILON) {
            tmp = fma(col, low, coeffs[i]); //col * low + coeffs[i];
            coh = fma(coh, high, coeffs[i]); //coh * high + coeffs[i];
        } else if (coh < DBL_EPSILON) {
            tmp = fma(col, high, coeffs[i]); //col * high + coeffs[i];
            coh = fma(coh, low, coeffs[i]); //coh * low + coeffs[i];
        } else {
            tmp = fma(col, high, coeffs[i]); //col * high + coeffs[i];
            coh = fma(coh, high, coeffs[i]); //coh * high + coeffs[i];
        }
        col = tmp;
    }
    *colow = col;
    *cohigh = coh;
}

static inline double min(double a, double b) {
    return a < b ? a : b;
}

static inline double max(double a, double b) {
    return b > a ? b : a;
}

static inline void intervalhorner(size_t const N, double const coeffs[N], double const low, double const high, double * const colow, double * const cohigh) {
    double col = coeffs[0];
    double coh = coeffs[0];
    double tmp;
    for (size_t i = 1; i < N; ++i) {
        if (col > DBL_EPSILON) {
            if (low > -DBL_EPSILON) {
                tmp = fma(col, low, coeffs[i]);
                coh = fma(coh, high, coeffs[i]);
            } else if (high < DBL_EPSILON) {
                tmp = fma(coh, low, coeffs[i]); 
                coh = fma(col, high, coeffs[i]);
            } else {
                tmp = fma(coh, low, coeffs[i]); 
                coh = fma(coh, high, coeffs[i]);
            }
        } else if (coh < DBL_EPSILON) {
            if (low > -DBL_EPSILON) {
                tmp = fma(col, high, coeffs[i]);
                coh = fma(coh, low, coeffs[i]);
            } else if (high < DBL_EPSILON) {
                tmp = fma(coh, high, coeffs[i]);
                coh = fma(col, low, coeffs[i]);
            } else {
                tmp = fma(col, high, coeffs[i]);
                coh = fma(col, low, coeffs[i]);
            }
        } else {
            if (low > -DBL_EPSILON) {
                tmp = fma(col, high, coeffs[i]);
                coh = fma(coh, high, coeffs[i]);
            } else if (high < DBL_EPSILON) {
                tmp = fma(coh, low, coeffs[i]);
                coh = fma(col, low, coeffs[i]);
            } else {
                tmp = min(fma(coh, low, coeffs[i]), fma(col, high, coeffs[i]));
                coh = max(fma(col, low, coeffs[i]), fma(coh, high, coeffs[i]));
            }    
        }
        col = tmp;
    }
    *colow = col;
    *cohigh = coh;
}

double smallestpositiveroot(size_t const N, double const coeffs[N]);
size_t realroots(size_t const N, double const coeffs[N], double (* const roots)[N-1]);