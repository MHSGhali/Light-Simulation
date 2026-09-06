/* integrator.h — the transport estimators.
 *
 * M2 provides direct lighting only. Multi-bounce path tracing with MIS lands in
 * M3, after the direct estimator has been validated against closed forms.
 */
#ifndef LIGHTSIM_INTEGRATOR_H
#define LIGHTSIM_INTEGRATOR_H

#include "scene.h"
#include "rng.h"

/* Direct spectral irradiance at `p` on a surface with normal `n`:
 *
 *     E = integral over the hemisphere of L_i(p,w) cos(theta) dw
 *
 * estimated by next-event estimation against every light:
 *
 *     E ~= sum over lights of  (1/N) sum over samples of
 *              li_over_pdf * cos(theta_p) * V(p,y)
 *
 * where li_over_pdf is L/pdf_w for area lights and I/r^2 for delta lights, so
 * both share one code path (see LightSample).
 *
 * M2 loops over ALL lights rather than selecting one stochastically. With a
 * handful of lights that is lower variance, and it keeps the light-selection
 * probability q_i out of the validation entirely -- q_i arrives with MIS in M3.
 *
 * Delta lights are sampled exactly once: their sample is deterministic, so
 * averaging N copies of it would only waste work.
 *
 * If `a_row` is non-NULL it must have room for sc->nlights doubles, and each
 * light's band-integrated contribution is deposited into its own column. This
 * is the per-source attribution that makes the contribution matrix fall out of
 * a single traversal.
 */
void ls_estimate_irradiance(const Scene *sc, vec3 p, vec3 n, int nsamples,
                            Rng *rng, SpectrumAcc *out, ls_real *a_row);

/* Which sampling strategies the path tracer uses to find emitted light.
 *
 * NEE and BSDF are each independently unbiased, so rendering the same scene
 * three ways and requiring agreement is what catches MIS bookkeeping errors --
 * a wrong weight makes MIS disagree with two strategies that cannot both be
 * wrong in the same direction. MIS should also show the lowest variance. */
typedef enum {
    LS_STRAT_BSDF = 1,   /* find emitters only by scattering into them   */
    LS_STRAT_NEE  = 2,   /* find emitters only by explicit connection    */
    LS_STRAT_MIS  = 3    /* both, combined with the power heuristic      */
} LsStrategy;

/* Power heuristic (beta = 2). Exposed so the tests can assert its partition
 * of unity: mis_power2(a,b) + mis_power2(b,a) == 1. */
static inline ls_real ls_mis_power2(ls_real pa, ls_real pb) {
    ls_real a = pa * pa, b = pb * pb, d = a + b;
    return d > 0.0 ? a / d : 0.0;
}

/* Estimate spectral radiance arriving along `ray`, in W/(m^2 sr nm).
 *
 * `a_row`, if non-NULL, receives per-source attribution: every contribution
 * terminates at exactly one emitter, whether found by NEE or by scattering and
 * whether after zero or twelve bounces, so one traversal fills every column. */
void ls_trace_radiance(const Scene *sc, Ray ray, Rng *rng, int max_depth,
                       LsStrategy strat, SpectrumAcc *out, ls_real *a_row);

/* As above, but emitted radiance is ignored at path vertices shallower than
 * `skip_emission_before`. Used by the probe estimator: direct light is already
 * counted by explicit light sampling at the probe, so collecting emission at
 * the first hit as well would double count it. */
void ls_trace_radiance_ex(const Scene *sc, Ray ray, Rng *rng, int max_depth,
                          LsStrategy strat, int skip_emission_before,
                          SpectrumAcc *out, ls_real *a_row);

/* Full spectral irradiance at `p` on a surface with normal `n`, INCLUDING
 * interreflection -- the quantity the analytic cos(theta)/r^2 model cannot
 * produce.
 *
 *   E = E_direct + E_indirect
 *
 * E_direct comes from explicit light sampling with visibility (so it accounts
 * for shadowing). E_indirect cosine-samples the hemisphere and traces:
 * with pdf = cos/pi the cosine cancels and the estimator is (pi/N) sum L_i,
 * with first-hit emission suppressed so the two halves do not overlap.
 *
 * `a_row`, if non-NULL, receives each source's own contribution. */
void ls_estimate_irradiance_full(const Scene *sc, vec3 p, vec3 n,
                                 int direct_samples, int indirect_samples,
                                 int max_depth, Rng *rng,
                                 SpectrumAcc *out, ls_real *a_row);

#endif /* LIGHTSIM_INTEGRATOR_H */
