#include "lightsim/integrator.h"

static bool light_is_delta(const Light *l) {
    return l->kind == LS_LIGHT_POINT
        || l->kind == LS_LIGHT_SPOT
        || l->kind == LS_LIGHT_DIRECTIONAL;
}

void ls_estimate_irradiance(const Scene *sc, vec3 p, vec3 n, int nsamples,
                            Rng *rng, SpectrumAcc *out, ls_real *a_row) {
    *out = ls_acc_zero();
    if (a_row)
        for (int i = 0; i < sc->nlights; ++i) a_row[i] = 0.0;

    for (int li = 0; li < sc->nlights; ++li) {
        const Light *l = &sc->lights[li];
        const bool delta = light_is_delta(l);
        const int  ns = delta ? 1 : nsamples;

        SpectrumAcc contrib = ls_acc_zero();
        for (int k = 0; k < ns; ++k) {
            LightSample s;
            ls_real u1 = delta ? 0.0 : ls_rng_f(rng);
            ls_real u2 = delta ? 0.0 : ls_rng_f(rng);
            if (!ls_light_sample(l, p, u1, u2, &s)) continue;

            ls_real cos_p = v3dot(n, s.wi);
            if (cos_p <= 0.0) continue;              /* below the horizon */
            if (ls_scene_occluded(sc, p, n, s.wi, s.dist)) continue;

            ls_acc_add_scaled(&contrib, &s.li_over_pdf, cos_p);
        }

        ls_real inv = 1.0 / (ls_real)ns;
        for (int i = 0; i < LS_NBINS; ++i) out->v[i] += contrib.v[i] * inv;

        if (a_row) {
            /* Band integral of this light's contribution: its column entry. */
            Spectrum m = ls_acc_mean(&contrib, (uint64_t)ns);
            a_row[li] = ls_spectrum_integrate(&m);
        }
    }
}
