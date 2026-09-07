/* test.h — minimal assertion harness. Each test file defines ls_test_main(). */
#ifndef LIGHTSIM_TEST_H
#define LIGHTSIM_TEST_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern int ls_test_failures;
extern int ls_test_count;

#define CHECK(cond)                                                            \
    do {                                                                       \
        ls_test_count++;                                                       \
        if (!(cond)) {                                                         \
            ls_test_failures++;                                                \
            printf("  FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);           \
        }                                                                      \
    } while (0)

/* Relative comparison, falling back to absolute near zero. */
#define CHECK_NEAR(got, want, tol)                                             \
    do {                                                                       \
        ls_test_count++;                                                       \
        double g_ = (double)(got), w_ = (double)(want), t_ = (double)(tol);     \
        double d_ = fabs(g_ - w_);                                             \
        double r_ = (fabs(w_) > 1e-12) ? d_ / fabs(w_) : d_;                   \
        if (!(r_ <= t_)) {                                                     \
            ls_test_failures++;                                                \
            printf("  FAIL %s:%d  %s\n        got  %.10g\n        want %.10g"  \
                   "  (rel err %.3g > tol %.3g)\n",                            \
                   __FILE__, __LINE__, #got, g_, w_, r_, t_);                  \
        }                                                                      \
    } while (0)

#define SECTION(name) printf("\n[%s]\n", name)
#define NOTE(...)     do { printf("  note: "); printf(__VA_ARGS__); printf("\n"); } while (0)

#endif /* LIGHTSIM_TEST_H */
