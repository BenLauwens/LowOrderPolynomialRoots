#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <stdio.h>

#include "fcsprp.h"

double smallestpositiveroot(size_t const N, double const coeffs[N]) {
    if (N == 1) {
        return DBL_MAX;
    } else if (N == 2) {
        double res = - coeffs[1] / coeffs[0];
        if (res < 0.0) {
            res = DBL_MAX;
        }
        return res;
    }
    double const coeffs1 = 1.0 / coeffs[0];
    double poly[N];
    poly[0] = 1.0;
    size_t const N1 = N - 1;
    double dpoly[N1];
    double low = 0.0;
    double high = 1.0;
    for (size_t i = 1; i < N; ++i) {
        poly[i] = coeffs[i] * coeffs1;
        dpoly[i-1] = (N-i) * poly[i-1];
        if (poly[i] < high) high = poly[i];
    }
    if (high > 0.0) {
        return DBL_MAX;
    }
    
    high = 1.0 - high;
    int index = -1;
    double list[2*N+1][2];
    double mid, comid;
    double colow, cohigh, dcolow, dcohigh;
    double leftlow, lefthigh, rightlow, righthigh;
    while (true) {
        mid = 0.5 * (low + high);
        horner(N, poly, mid, &comid);
        posintervalhorner(N1, dpoly, low, high, &dcolow, &dcohigh);
        if (dcolow < 0.0 && dcohigh > 0.0) {
            if (comid > 0.0) {
                leftlow = - DBL_MAX;
                lefthigh = mid - comid / dcohigh;
                rightlow = mid - comid / dcolow;
                righthigh = DBL_MAX;
            } else if (comid < 0.0) {
                leftlow = - DBL_MAX;
                lefthigh = mid - comid / dcolow;
                rightlow = mid - comid / dcohigh;
                righthigh = DBL_MAX;
            } else {
                return mid;
            }
            if (!(high < leftlow || lefthigh < low)) {
                if (!(high < rightlow || righthigh < low)) {
                    index += 1;
                    if (low < rightlow) {
                        list[index][0] = rightlow;
                    } else {
                        list[index][0] = low;
                    }
                    if (high < righthigh) {
                        list[index][1] = high; 
                    } else {
                        list[index][1] = righthigh;
                    }
                }
                if (low < leftlow) low = leftlow;
                if (!(high < lefthigh)) high = lefthigh;
            } else if (!(high < rightlow || righthigh < low)) {
                if (low < rightlow) low = rightlow;
                if (!(high < righthigh)) high = righthigh;
            } else if (index == -1) {
                return DBL_MAX;
            } else {
                low = list[index][0];
                high = list[index][1];
                index -= 1;
            }
        } else {
            horner2(N, poly, low, high, &colow, &cohigh);
            if (colow * cohigh < 0.0) {
                double dcomid, delta, res;
                while (true) {
                    hornerd(N, poly, mid, &comid, &dcomid);
                    delta = comid / dcomid;
                    res = mid - delta;
                    if (fabs(delta) < 1.0e-8 * mid) {
                        return res;
                    } else if (low < res && res < high) {
                        mid = res;
                    } else {
                        if (comid * colow > 0.0) {
                            low = mid;
                            colow = comid;
                        } else if (comid * cohigh > 0.0) {
                            high = mid;
                            cohigh = comid;
                        } else {
                            return mid;
                        }
                        mid = (low * cohigh - high * colow) / (cohigh - colow);
                        if (!(low < mid && mid < high)) mid = 0.5 * (low + high);
                    }
                }
            } else if (index == -1) {
                return DBL_MAX;
            } else {
                low = list[index][0];
                high = list[index][1];
                index -= 1;
            }
        }
    }
}

size_t realroots(size_t const N, double const coeffs[N], double (* const roots)[N-1]) {
    if (N == 1) {
        return 0;
    } else if (N == 2) {
        (*roots)[0] = - coeffs[1] / coeffs[0];
        return 1;
    }
    double const coeffs1 = 1.0 / coeffs[0];
    double poly[N];
    poly[0] = 1.0;
    size_t const N1 = N - 1;
    double dpoly[N1];
    double low = 1.0;
    double high = 1.0;
    double s = -1.0;
    double poly1;
    double poly2;
    for (size_t i = 1; i < N; ++i) {
        dpoly[i-1] = (N-i) * poly[i-1];
        poly1 = coeffs[i] * coeffs1;
        poly[i] = poly1;
        if (poly1 < high) high = poly[i];
        poly2 = s * poly1;
        s = -s;
        if (poly2 < low) low = poly2;
    }
    if (low > 0.0) low = 1.0;
    if (high > 0.0) high = 1.0;
    if (low == 1.0 && high == 1.0) return 0;
    low = low - 1.0;
    high = 1.0 - high;
    int index = -1;
    double list[2*N+1][2];
    size_t counter = 0;
    double mid, comid;
    double colow, cohigh, dcolow, dcohigh;
    double leftlow, lefthigh, rightlow, righthigh;
    while (true) {
        //printf("low: %+.16f high: %+.16f\n", low, high);
        mid = 0.5 * (low + high);
        horner(N, poly, mid, &comid);
        intervalhorner(N1, dpoly, low, high, &dcolow, &dcohigh);
        //printf("%+.16f %+.16f %+.16f %+.16f\n", mid, comid, dcolow, dcohigh);
        if (dcolow < 0.0 && dcohigh > 0.0) {
            if (comid > 0.0) {
                leftlow = - DBL_MAX;
                lefthigh = mid - comid / dcohigh;
                rightlow = mid - comid / dcolow;
                righthigh = DBL_MAX;
            } else if (comid < 0.0) {
                leftlow = - DBL_MAX;
                lefthigh = mid - comid / dcolow;
                rightlow = mid - comid / dcohigh;
                righthigh = DBL_MAX;
            } else {
                (*roots)[counter] = mid;
                counter += 1;
                if (index == -1) return counter;
                low = list[index][0];
                high = list[index][1];
                index -= 1;
                continue;
            }
            if (!(high < leftlow || lefthigh < low)) {
                if (!(high < rightlow || righthigh < low)) {
                    index += 1;
                    if (low < rightlow) {
                        list[index][0] = rightlow;
                    } else {
                        list[index][0] = low;
                    }
                    if (high < righthigh) {
                        list[index][1] = high; 
                    } else {
                        list[index][1] = righthigh;
                    }
                }
                if (low < leftlow) low = leftlow;
                if (!(high < lefthigh)) high = lefthigh;
            } else if (!(high < rightlow || righthigh < low)) {
                if (low < rightlow) low = rightlow;
                if (!(high < righthigh)) high = righthigh;
            } else if (index == -1) {
                return counter;
            } else {
                low = list[index][0];
                high = list[index][1];
                index -= 1;
            }
        } else {
            horner2(N, poly, low, high, &colow, &cohigh);
            //printf("%+.16f %+.16f\n", colow, cohigh);
            if (colow * cohigh < 0.0) {
                double dcomid, delta, newmid;
                while (true) {
                    hornerd(N, poly, mid, &comid, &dcomid);
                    delta = comid / dcomid;
                    newmid = mid - delta;
                    if (fabs(delta) < 1.0e-8 * fabs(mid)) {
                        (*roots)[counter] = newmid;
                        counter += 1;
                        if (index == -1) return counter;
                        low = list[index][0];
                        high = list[index][1];
                        index -= 1;
                        break;
                    } else if (low < newmid && newmid < high) {
                        mid = newmid;
                    } else {
                        if (comid * colow > 0.0) {
                            low = mid;
                            colow = comid;
                        } else if (comid * cohigh > 0.0) {
                            high = mid;
                            cohigh = comid;
                        } else {
                            (*roots)[counter] = mid;
                            counter += 1;
                            if (index == -1) return counter;
                            low = list[index][0];
                            high = list[index][1];
                            index -= 1;
                            break;
                        }
                        mid = (low * cohigh - high * colow) / (cohigh - colow);
                        if (!(low < mid && mid < high)) mid = 0.5 * (low + high);
                    }
                }
            } else if (index == -1) {
                return counter;
            } else {
                low = list[index][0];
                high = list[index][1];
                index -= 1;
            }
        }
    }
}
