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

#endif /* LIGHTSIM_INTEGRATOR_H */
