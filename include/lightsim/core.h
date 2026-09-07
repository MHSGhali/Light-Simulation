/* core.h — fundamental types, config and assertions for the lightsim engine. */
#ifndef LIGHTSIM_CORE_H
#define LIGHTSIM_CORE_H

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/* All internal computation is double; only bulk spectral storage is float. */
typedef double ls_real;

#define LS_PI       3.14159265358979323846
#define LS_INV_PI   0.31830988618379067154
#define LS_TWO_PI   6.28318530717958647692
#define LS_INV_4PI  0.07957747154594766788

/* Photometric constant: maximum luminous efficacy of radiation at 555 nm. */
#define LS_KM_LM_PER_W 683.0

/* Physical constants (SI, CODATA / SI-2019 exact values). */
#define LS_PLANCK_H    6.62607015e-34   /* J s   */
#define LS_LIGHT_C     2.99792458e8     /* m/s   */
#define LS_BOLTZMANN_K 1.380649e-23     /* J/K   */
#define LS_STEFAN_BOLTZMANN 5.670374419e-8 /* W/(m^2 K^4) */

#define LS_UNUSED(x) ((void)(x))

/* Debug-only invariant. Compiled out with NDEBUG. */
#define LS_ASSERT(cond) assert(cond)

static inline ls_real ls_clamp(ls_real x, ls_real lo, ls_real hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}
static inline ls_real ls_min(ls_real a, ls_real b) { return a < b ? a : b; }
static inline ls_real ls_max(ls_real a, ls_real b) { return a > b ? a : b; }
static inline ls_real ls_lerp(ls_real t, ls_real a, ls_real b) { return a + t * (b - a); }
static inline ls_real ls_sqr(ls_real x) { return x * x; }

#endif /* LIGHTSIM_CORE_H */
