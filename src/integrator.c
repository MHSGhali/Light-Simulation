#include "lightsim/integrator.h"
#include <stdlib.h>

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

/* ---- path tracing ---- */

static const Material *mat_of(const Scene *sc, int mat_id) {
    return (mat_id >= 0 && mat_id < sc->nmats) ? &sc->mats[mat_id] : NULL;
}

/* Deposit a contribution into the radiance accumulator and, if requested, into
 * the emitting source's own column. */
static void deposit(SpectrumAcc *out, ls_real *a_row, const Spectrum *beta,
                    const Spectrum *value, ls_real w, int light_id) {
    Spectrum c = ls_spectrum_mul(*beta, *value);
    ls_acc_add_scaled(out, &c, w);
    if (a_row && light_id >= 0)
        a_row[light_id] += ls_spectrum_integrate(&c) * w;
}

void ls_trace_radiance(const Scene *sc, Ray ray, Rng *rng, int max_depth,
                       LsStrategy strat, SpectrumAcc *out, ls_real *a_row) {
    ls_trace_radiance_ex(sc, ray, rng, max_depth, strat, 0, out, a_row);
}

void ls_trace_radiance_ex(const Scene *sc, Ray ray, Rng *rng, int max_depth,
                          LsStrategy strat, int skip_emission_before,
                          SpectrumAcc *out, ls_real *a_row) {
    *out = ls_acc_zero();
    if (a_row)
        for (int i = 0; i < sc->nlights; ++i) a_row[i] = 0.0;

    Spectrum beta = ls_spectrum_const(1.0);   /* path throughput, dimensionless */
    bool    prev_delta = true;                /* camera vertex: emission counts */
    ls_real pdf_bsdf_prev = 0.0;
    const bool use_nee  = (strat & LS_STRAT_NEE)  != 0;
    const bool use_bsdf = (strat & LS_STRAT_BSDF) != 0;

    for (int depth = 0; ; ++depth) {
        Hit h;
        if (!ls_scene_intersect(sc, &ray, &h)) break;

        const Material *m = mat_of(sc, h.mat_id);
        if (!m) break;

        /* Orient the shading frame against the incoming ray so that both sides
         * of a surface shade correctly. ng itself is never mutated. */
        vec3 ns = h.backface ? v3neg(h.ng) : h.ng;

        /* ---- emitted radiance found by scattering ---- */
        if (m->emissive && use_bsdf && depth >= skip_emission_before) {
            ls_real w = 1.0;
            if (!prev_delta && h.light_id >= 0 && (strat & LS_STRAT_NEE)) {
                /* h.ng, NOT ns. ls_light_pdf_w returns 0 when the emitter is
                 * seen from behind; handing it the ray-facing normal makes that
                 * cosine positive on a backface hit, so BSDF sampling would
                 * claim a density for a connection NEE refuses to make, and the
                 * MIS weights would stop summing to one. */
                ls_real pl = ls_light_pdf_w(&sc->lights[h.light_id], ray.o, h.p, h.ng);
                w = ls_mis_power2(pdf_bsdf_prev, pl);
            }
            Spectrum le = m->le;
            /* Also h.ng. `ns` has already been flipped to face the incoming ray,
             * so dot(ns, -ray.d) is >= 0 for every hit and this gate could never
             * fire -- emitters found by scattering were effectively two-sided
             * while ls_light_radiance and ls_light_pdf_w treated them as
             * one-sided, so NEE and BSDF sampling disagreed about the back of
             * every area light. */
            if (v3dot(h.ng, v3neg(ray.d)) <= 0.0) le = ls_spectrum_zero();
            deposit(out, a_row, &beta, &le, w, h.light_id);
        }

        if (depth >= max_depth) break;

        Basis fr = ls_basis(ns);
        vec3  wo = ls_basis_to_local(fr, v3neg(ray.d));
        if (wo.z <= 0.0) break;

        /* ---- next-event estimation ---- */
        if (use_nee && !ls_bsdf_is_delta(&m->bsdf)) {
            for (int li = 0; li < sc->nlights; ++li) {
                const Light *l = &sc->lights[li];
                LightSample s;
                if (!ls_light_sample(l, h.p, ls_rng_f(rng), ls_rng_f(rng), &s)) continue;
                vec3 wi = ls_basis_to_local(fr, s.wi);
                if (wi.z <= 0.0) continue;
                if (ls_scene_occluded(sc, h.p, ns, s.wi, s.dist)) continue;

                Spectrum f;
                ls_bsdf_eval(&m->bsdf, wo, wi, &f);
                if (ls_spectrum_is_black(&f)) continue;

                ls_real w = 1.0;
                if (s.pdf_w > 0.0 && use_bsdf) {
                    /* Not a delta light, and BSDF sampling could also have
                     * found it: combine. Delta lights keep w = 1 because BSDF
                     * sampling can never hit them, so there is no double count. */
                    w = ls_mis_power2(s.pdf_w, ls_bsdf_pdf(&m->bsdf, wo, wi));
                }
                Spectrum c = ls_spectrum_mul(f, s.li_over_pdf);
                deposit(out, a_row, &beta, &c, w * wi.z, l->index);
            }
        }

        /* ---- continue the path ---- */
        vec3     wi_l;
        Spectrum f;
        ls_real  pdf;
        if (!ls_bsdf_sample(&m->bsdf, wo, ls_rng_f(rng), ls_rng_f(rng), &wi_l, &f, &pdf))
            break;
        if (pdf <= 0.0) break;

        beta = ls_spectrum_mul(beta, ls_spectrum_scale(f, wi_l.z / pdf));
        prev_delta = ls_bsdf_is_delta(&m->bsdf);
        pdf_bsdf_prev = prev_delta ? 0.0 : pdf;

        vec3 wi_w = ls_basis_to_world(fr, wi_l);
        ray.o = ls_offset_origin(h.p, ns, wi_w);
        ray.d = wi_w;
        ray.tmin = 0.0;
        ray.tmax = HUGE_VAL;

        /* ---- Russian roulette ----
         * The survival probability is a SCALAR function of the throughput.
         * A per-bin decision would decorrelate the wavelengths and destroy the
         * meaning of the spectral throughput. Applied after both deposits. */
        if (depth >= 3) {
            ls_real q = ls_min(0.95, ls_spectrum_max(&beta));
            if (q <= 1e-6 || ls_rng_f(rng) >= q) break;
            beta = ls_spectrum_scale(beta, 1.0 / q);
        }
    }
}


void ls_estimate_irradiance_full(const Scene *sc, vec3 p, vec3 n,
                                 int direct_samples, int indirect_samples,
                                 int max_depth, Rng *rng,
                                 SpectrumAcc *out, ls_real *a_row) {
    /* Direct: explicit light sampling with visibility. */
    ls_estimate_irradiance(sc, p, n, direct_samples, rng, out, a_row);
    if (indirect_samples <= 0 || max_depth <= 0) return;

    /* Indirect: cosine-sample the hemisphere. With pdf = cos/pi the cosine in
     * the irradiance integral cancels exactly, leaving E = (pi/N) sum L_i. */
    Basis fr = ls_basis(n);
    SpectrumAcc ind = ls_acc_zero();
    ls_real *row = NULL;
    if (a_row) {
        row = calloc((size_t)sc->nlights, sizeof *row);
        if (!row) return;
    }

    for (int k = 0; k < indirect_samples; ++k) {
        vec3 wl = ls_sample_hemisphere_cosine(ls_rng_f(rng), ls_rng_f(rng));
        if (wl.z <= 0.0) continue;
        vec3 wi = ls_basis_to_world(fr, wl);

        Ray r;
        r.o = ls_offset_origin(p, n, wi);
        r.d = wi;
        r.tmin = 0.0;
        r.tmax = HUGE_VAL;

        SpectrumAcc L;
        /* skip_emission_before = 1: the first hit's emission is direct light,
         * already counted above. Everything deeper is genuine indirect. */
        ls_trace_radiance_ex(sc, r, rng, max_depth, LS_STRAT_MIS, 1, &L, row);
        for (int b = 0; b < LS_NBINS; ++b) ind.v[b] += L.v[b];
        if (a_row && row)
            for (int i = 0; i < sc->nlights; ++i) a_row[i] += row[i] * LS_PI / (ls_real)indirect_samples;
    }

    ls_real w = LS_PI / (ls_real)indirect_samples;
    for (int b = 0; b < LS_NBINS; ++b) out->v[b] += ind.v[b] * w;
    free(row);
}
