#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_poly.h>

#include "fcsprp.h"

double average(size_t const N, double const values[N]) {
    double res = 0.0;
    for (size_t i = 0; i < N; i++) {
        res += values[i];
    }
    return res / N;
}

double minval(size_t const N, double const values[N]) {
    double res = DBL_MAX;
    for (size_t i = 0; i < N; i++) {
        if (values[i] < res) {
            res = values[i];
        }   
    }
    return res;
}

double maxval(size_t const N, double const values[N]) {
    double res = -DBL_MAX;
    for (size_t i = 0; i < N; i++) {
        if (values[i] > res) {
            res = values[i];
        }   
    }
    return res;
}

double stddev(size_t const N, double const values[N]) {
    double res = 0.0;
    for (size_t i = 0; i < N; i++) {
        res += pow(values[i], 2);
    }
    return sqrt((res - pow(average(N, values), 2) / N) / (N-1));
}

int main(void) {
    int const N = 10000;
    int const S = 100;

    double times[N];
    double res;

    double bench[S];
    /*
    double benchlib[100];
    for (size_t k = 0; k < S; ++k) {
        double const coeffs[6]= {
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
        };
        printf("%f %f %f %f %f %f\n", coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4], coeffs[5]);
        double const coeffs_1[6] = {
            coeffs[5], coeffs[4], coeffs[3], coeffs[2], coeffs[1], coeffs[0]
        };
        double z[10];
        gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(6);
        for (size_t i = 0; i < N; ++i) {
            clock_t start = clock();
            for (size_t j = 0; j < 1000; ++j) {
                gsl_poly_complex_solve (coeffs_1, 6, w, z);
            }
            clock_t finish = clock();
            times[i] = 1000000 * (finish - start) / CLOCKS_PER_SEC;
        }
        for (size_t i = 0; i < 5; ++i) {
            printf("z%zu = %+.16f %+.16f\n", i, z[2*i], z[2*i+1]);
        }
        benchlib[k] = minval(N, times);
        double r[5];
        size_t n;
        for (size_t i = 0; i < N; ++i) {
            clock_t start = clock();
            for (size_t j = 0; j < 1000; ++j) {
                n = realroots(6, coeffs, &r);
            }
            clock_t finish = clock();
            times[i] = 1000000 * (finish - start) / CLOCKS_PER_SEC;
        }
        for (size_t i = 0; i < n; ++i) {
            printf("r%zu = %+.16f\n", i, r[i]);
        }
        bench[k] = minval(N, times);
    }
    printf("%i %16f %16f %16f %16f\n", 6, minval(100, benchlib), average(100, benchlib), maxval(100, benchlib), stddev(100, benchlib));
    printf("%i %16f %16f %16f %16f\n", 6, minval(100, bench), average(100, bench), maxval(100, bench), stddev(100, bench));*/
    
    for (size_t k = 0; k < S; ++k) {
        double const coeffs[3]= {//
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
        };
        double r[3];
        for (size_t i = 0; i < N; ++i) {
            clock_t start = clock();
            for (size_t j = 0; j < 1000; ++j) {
                realroots(3, coeffs, &r);
            }
            clock_t finish = clock();
            times[i] = 1000000 * (finish - start) / CLOCKS_PER_SEC;
        }
        bench[k] = minval(N, times);
    }
    printf("%i %16f %16f %16f %16f\n", 3, minval(100, bench), average(100, bench), maxval(100, bench), stddev(100, bench));

    for (size_t k = 0; k < S; ++k) {
        double const coeffs[4]= {//
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
        };
        double r[4];
        for (size_t i = 0; i < N; ++i) {
            clock_t start = clock();
            for (size_t j = 0; j < 1000; ++j) {
                realroots(4, coeffs, &r);
            }
            clock_t finish = clock();
            times[i] = 1000000 * (finish - start) / CLOCKS_PER_SEC;
        }
        bench[k] = minval(N, times);
    }
    printf("%i %16f %16f %16f %16f\n", 4, minval(100, bench), average(100, bench), maxval(100, bench), stddev(100, bench));

    for (size_t k = 0; k < S; ++k) {
        double const coeffs[5]= {//
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
        };
        double r[5];
        for (size_t i = 0; i < N; ++i) {
            clock_t start = clock();
            for (size_t j = 0; j < 1000; ++j) {
                realroots(5, coeffs, &r);
            }
            clock_t finish = clock();
            times[i] = 1000000 * (finish - start) / CLOCKS_PER_SEC;
        }
        bench[k] = minval(N, times);
    }
    printf("%i %16f %16f %16f %16f\n", 5, minval(100, bench), average(100, bench), maxval(100, bench), stddev(100, bench));

    for (size_t k = 0; k < S; ++k) {
        double const coeffs[6]= {//
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
        };
        double r[6];
        for (size_t i = 0; i < N; ++i) {
            clock_t start = clock();
            for (size_t j = 0; j < 1000; ++j) {
                realroots(6, coeffs, &r);
            }
            clock_t finish = clock();
            times[i] = 1000000 * (finish - start) / CLOCKS_PER_SEC;
        }
        bench[k] = minval(N, times);
    }
    printf("%i %16f %16f %16f %16f\n", 6, minval(100, bench), average(100, bench), maxval(100, bench), stddev(100, bench));

    for (size_t k = 0; k < S; ++k) {
        double const coeffs[7]= {//
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
        };
        double r[7];
        for (size_t i = 0; i < N; ++i) {
            clock_t start = clock();
            for (size_t j = 0; j < 1000; ++j) {
                realroots(7, coeffs, &r);
            }
            clock_t finish = clock();
            times[i] = 1000000 * (finish - start) / CLOCKS_PER_SEC;
        }
        bench[k] = minval(N, times);
    }
    printf("%i %16f %16f %16f %16f\n", 7, minval(100, bench), average(100, bench), maxval(100, bench), stddev(100, bench));

    for (size_t k = 0; k < S; ++k) {
        double const coeffs[8]= {//
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
        };
        double r[8];
        for (size_t i = 0; i < N; ++i) {
            clock_t start = clock();
            for (size_t j = 0; j < 1000; ++j) {
                realroots(8, coeffs, &r);
            }
            clock_t finish = clock();
            times[i] = 1000000 * (finish - start) / CLOCKS_PER_SEC;
        }
        bench[k] = minval(N, times);
    }
    printf("%i %16f %16f %16f %16f\n", 8, minval(100, bench), average(100, bench), maxval(100, bench), stddev(100, bench));

    for (size_t k = 0; k < S; ++k) {
        double const coeffs[9]= {//
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
        };
        double r[9];
        for (size_t i = 0; i < N; ++i) {
            clock_t start = clock();
            for (size_t j = 0; j < 1000; ++j) {
                realroots(9, coeffs, &r);
            }
            clock_t finish = clock();
            times[i] = 1000000 * (finish - start) / CLOCKS_PER_SEC;
        }
        bench[k] = minval(N, times);
    }
    printf("%i %16f %16f %16f %16f\n", 9, minval(100, bench), average(100, bench), maxval(100, bench), stddev(100, bench));

    for (size_t k = 0; k < S; ++k) {
        double const coeffs[3]= {
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
        };
        for (size_t i = 0; i < N; ++i) {
            clock_t start = clock();
            for (size_t j = 0; j < 1000; ++j) {
                res = smallestpositiveroot(3, coeffs);
            }
            clock_t finish = clock();
            times[i] = 1000000 * (finish - start) / CLOCKS_PER_SEC;
        }
        bench[k] = minval(N, times);
    }
    printf("%i %16f %16f %16f %16f\n", 3, minval(100, bench), average(100, bench), maxval(100, bench), stddev(100, bench));

    for (size_t k = 0; k < S; ++k) {
        double const coeffs[4]= {
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
        };
        for (size_t i = 0; i < N; ++i) {
            clock_t start = clock();
            for (size_t j = 0; j < 1000; ++j) {
                res = smallestpositiveroot(4, coeffs);
            }
            clock_t finish = clock();
            times[i] = 1000000 * (finish - start) / CLOCKS_PER_SEC;
        }
        bench[k] = minval(N, times);
    }
    printf("%i %16f %16f %16f %16f\n", 4, minval(100, bench), average(100, bench), maxval(100, bench), stddev(100, bench));

    for (size_t k = 0; k < S; ++k) {
        double const coeffs[5]= {
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
        };
        for (size_t i = 0; i < N; ++i) {
            clock_t start = clock();
            for (size_t j = 0; j < 1000; ++j) {
                res = smallestpositiveroot(5, coeffs);
            }
            clock_t finish = clock();
            times[i] = 1000000 * (finish - start) / CLOCKS_PER_SEC;
        }
        bench[k] = minval(N, times);
    }
    printf("%i %16f %16f %16f %16f\n", 5, minval(100, bench), average(100, bench), maxval(100, bench), stddev(100, bench));

    for (size_t k = 0; k < S; ++k) {
        double const coeffs[6]= {
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
        };
        for (size_t i = 0; i < N; ++i) {
            clock_t start = clock();
            for (size_t j = 0; j < 1000; ++j) {
                res = smallestpositiveroot(6, coeffs);
            }
            clock_t finish = clock();
            times[i] = 1000000 * (finish - start) / CLOCKS_PER_SEC;
        }
        bench[k] = minval(N, times);
    }
    printf("%i %16f %16f %16f %16f\n", 6, minval(100, bench), average(100, bench), maxval(100, bench), stddev(100, bench));

    for (size_t k = 0; k < S; ++k) {
        double const coeffs[7]= {
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
        };
        for (size_t i = 0; i < N; ++i) {
            clock_t start = clock();
            for (size_t j = 0; j < 1000; ++j) {
                res = smallestpositiveroot(7, coeffs);
            }
            clock_t finish = clock();
            times[i] = 1000000 * (finish - start) / CLOCKS_PER_SEC;
        }
        bench[k] = minval(N, times);
    }
    printf("%i %16f %16f %16f %16f\n", 7, minval(100, bench), average(100, bench), maxval(100, bench), stddev(100, bench));

    for (size_t k = 0; k < S; ++k) {
        double const coeffs[8]= {
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
        };
        for (size_t i = 0; i < N; ++i) {
            clock_t start = clock();
            for (size_t j = 0; j < 1000; ++j) {
                res = smallestpositiveroot(8, coeffs);
            }
            clock_t finish = clock();
            times[i] = 1000000 * (finish - start) / CLOCKS_PER_SEC;
        }
        bench[k] = minval(N, times);
    }
    printf("%i %16f %16f %16f %16f\n", 8, minval(100, bench), average(100, bench), maxval(100, bench), stddev(100, bench));

    for (size_t k = 0; k < S; ++k) {
        double const coeffs[9]= {
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
            2.0 * (double) rand() / (double) RAND_MAX - 1.0,
        };
        for (size_t i = 0; i < N; ++i) {
            clock_t start = clock();
            for (size_t j = 0; j < 1000; ++j) {
                res = smallestpositiveroot(9, coeffs);
            }
            clock_t finish = clock();
            times[i] = 1000000 * (finish - start) / CLOCKS_PER_SEC;
        }
        bench[k] = minval(N, times);
    }
    printf("%i %16f %16f %16f %16f\n", 9, minval(100, bench), average(100, bench), maxval(100, bench), stddev(100, bench));


    return EXIT_SUCCESS;
}