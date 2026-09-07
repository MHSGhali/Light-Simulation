/* rng.h — PCG32. Small, fast, and streamable: each thread/tile gets its own
 * stream from a distinct sequence constant, so results are reproducible
 * regardless of how work is scheduled. */
#ifndef LIGHTSIM_RNG_H
#define LIGHTSIM_RNG_H

#include "core.h"

typedef struct { uint64_t state, inc; } Rng;

static inline uint32_t ls_rng_u32(Rng *r) {
    uint64_t old = r->state;
    r->state = old * 6364136223846793005ULL + r->inc;
    uint32_t xorshifted = (uint32_t)(((old >> 18u) ^ old) >> 27u);
    uint32_t rot = (uint32_t)(old >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((32u - rot) & 31u));
}

/* Uniform in [0,1). The 2^-32 scaling never returns exactly 1.0. */
static inline ls_real ls_rng_f(Rng *r) {
    return (ls_real)ls_rng_u32(r) * 2.3283064365386963e-10;
}

static inline Rng ls_rng_seed(uint64_t seed, uint64_t stream) {
    Rng r;
    r.state = 0u;
    r.inc = (stream << 1u) | 1u;
    (void)ls_rng_u32(&r);
    r.state += seed;
    (void)ls_rng_u32(&r);
    return r;
}

#endif /* LIGHTSIM_RNG_H */
