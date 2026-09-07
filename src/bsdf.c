#include "lightsim/bsdf.h"
#include <math.h>

ls_real ls_fresnel_dielectric(ls_real cos_theta, ls_real eta) {
    cos_theta = ls_clamp(cos_theta, -1.0, 1.0);
    if (cos_theta < 0.0) { eta = 1.0 / eta; cos_theta = -cos_theta; }
    ls_real sin2_t = (1.0 - cos_theta * cos_theta) / (eta * eta);
    if (sin2_t >= 1.0) return 1.0;                    /* total internal reflection */
    ls_real cos_t = sqrt(1.0 - sin2_t);
    ls_real rs = (cos_theta - eta * cos_t) / (cos_theta + eta * cos_t);
    ls_real rp = (eta * cos_theta - cos_t) / (eta * cos_theta + cos_t);
    return 0.5 * (rs * rs + rp * rp);
}

ls_real ls_fresnel_conductor(ls_real cos_theta, ls_real n, ls_real k) {
    /* Unpolarised Fresnel for an absorbing medium. Standard formulation in
     * terms of a^2+b^2; numerically stable across the whole angle range. */
    cos_theta = ls_clamp(cos_theta, 0.0, 1.0);
    ls_real c2 = cos_theta * cos_theta;
    ls_real s2 = 1.0 - c2;
    ls_real n2 = n * n, k2 = k * k;

    ls_real t0 = n2 - k2 - s2;
    ls_real a2b2 = sqrt(ls_max(0.0, t0 * t0 + 4.0 * n2 * k2));
    ls_real t1 = a2b2 + c2;
    ls_real a = sqrt(ls_max(0.0, 0.5 * (a2b2 + t0)));
    ls_real t2 = 2.0 * a * cos_theta;
    ls_real rs = (t1 - t2) / (t1 + t2);

    ls_real t3 = c2 * a2b2 + s2 * s2;
    ls_real t4 = t2 * s2;
    ls_real rp = rs * (t3 - t4) / (t3 + t4);

    return 0.5 * (rs + rp);
}

ls_real ls_ggx_d(ls_real cos_theta_m, ls_real alpha) {
    if (cos_theta_m <= 0.0) return 0.0;
    ls_real a2 = alpha * alpha;
    ls_real c2 = cos_theta_m * cos_theta_m;
    ls_real t = c2 * (a2 - 1.0) + 1.0;
    return a2 / (LS_PI * t * t);
}

/* Smith Lambda term for GGX. */
static ls_real ggx_lambda(vec3 w, ls_real alpha) {
    ls_real c = fabs(w.z);
    if (c >= 1.0) return 0.0;
    ls_real tan2 = (1.0 - c * c) / (c * c);
    return 0.5 * (-1.0 + sqrt(1.0 + alpha * alpha * tan2));
}

ls_real ls_ggx_g1(vec3 w, ls_real alpha) {
    return 1.0 / (1.0 + ggx_lambda(w, alpha));
}

ls_real ls_ggx_g2(vec3 wo, vec3 wi, ls_real alpha) {
    /* Height-correlated Smith: less energy loss than the separable form. */
    return 1.0 / (1.0 + ggx_lambda(wo, alpha) + ggx_lambda(wi, alpha));
}

/* Sample the GGX distribution of VISIBLE normals (Heitz 2018). Sampling the
 * visible normals rather than D itself removes the samples that would be
 * masked, which is where most of the variance in a naive GGX sampler lives. */
static vec3 sample_ggx_vndf(vec3 wo, ls_real alpha, ls_real u1, ls_real u2) {
    vec3 vh = v3norm(v3(alpha * wo.x, alpha * wo.y, wo.z));
    ls_real len2 = vh.x * vh.x + vh.y * vh.y;
    vec3 t1 = len2 > 0.0 ? v3scale(v3(-vh.y, vh.x, 0.0), 1.0 / sqrt(len2))
                         : v3(1.0, 0.0, 0.0);
    vec3 t2 = v3cross(vh, t1);

    ls_real r = sqrt(u1);
    ls_real phi = LS_TWO_PI * u2;
    ls_real p1 = r * cos(phi);
    ls_real p2 = r * sin(phi);
    ls_real s = 0.5 * (1.0 + vh.z);
    p2 = (1.0 - s) * sqrt(ls_max(0.0, 1.0 - p1 * p1)) + s * p2;

    vec3 nh = v3add(v3add(v3scale(t1, p1), v3scale(t2, p2)),
                    v3scale(vh, sqrt(ls_max(0.0, 1.0 - p1 * p1 - p2 * p2))));
    return v3norm(v3(alpha * nh.x, alpha * nh.y, ls_max(1e-12, nh.z)));
}

bool ls_bsdf_is_delta(const Bsdf *b) {
    return b->kind == LS_BSDF_CONDUCTOR && b->alpha <= 0.0;
}

static void conductor_fresnel_spectrum(const Bsdf *b, ls_real cos_theta, Spectrum *out) {
    for (int i = 0; i < LS_NBINS; ++i)
        out->v[i] = (float)ls_fresnel_conductor(cos_theta,
                                                (ls_real)b->eta.v[i],
                                                (ls_real)b->kappa.v[i]);
}

void ls_bsdf_eval(const Bsdf *b, vec3 wo, vec3 wi, Spectrum *f_out) {
    *f_out = ls_spectrum_zero();
    if (wo.z <= 0.0 || wi.z <= 0.0) return;          /* same hemisphere only */

    if (b->kind == LS_BSDF_LAMBERT) {
        *f_out = ls_spectrum_scale(b->rho, LS_INV_PI);
        return;
    }
    if (b->alpha <= 0.0) return;                     /* delta: no density */

    vec3 m = v3norm(v3add(wo, wi));                  /* half vector */
    if (m.z <= 0.0) return;
    ls_real d  = ls_ggx_d(m.z, b->alpha);
    ls_real g2 = ls_ggx_g2(wo, wi, b->alpha);
    ls_real denom = 4.0 * wo.z * wi.z;
    if (denom <= 0.0) return;

    Spectrum fr;
    conductor_fresnel_spectrum(b, v3dot(wi, m), &fr);
    *f_out = ls_spectrum_scale(fr, d * g2 / denom);
}

ls_real ls_bsdf_pdf(const Bsdf *b, vec3 wo, vec3 wi) {
    if (wo.z <= 0.0 || wi.z <= 0.0) return 0.0;
    if (b->kind == LS_BSDF_LAMBERT) return ls_pdf_hemisphere_cosine(wi.z);
    if (b->alpha <= 0.0) return 0.0;                 /* delta */

    vec3 m = v3norm(v3add(wo, wi));
    if (m.z <= 0.0) return 0.0;
    /* pdf(wi) = D_visible(m) / (4 |wo . m|) */
    ls_real dot_om = v3dot(wo, m);
    if (dot_om <= 0.0) return 0.0;
    ls_real dv = ls_ggx_g1(wo, b->alpha) * dot_om * ls_ggx_d(m.z, b->alpha) / wo.z;
    return dv / (4.0 * dot_om);
}

bool ls_bsdf_sample(const Bsdf *b, vec3 wo, ls_real u1, ls_real u2,
                    vec3 *wi_out, Spectrum *f_out, ls_real *pdf_out) {
    if (wo.z <= 0.0) return false;

    if (b->kind == LS_BSDF_LAMBERT) {
        vec3 wi = ls_sample_hemisphere_cosine(u1, u2);
        *wi_out = wi;
        *pdf_out = ls_pdf_hemisphere_cosine(wi.z);
        if (*pdf_out <= 0.0) return false;
        *f_out = ls_spectrum_scale(b->rho, LS_INV_PI);
        return true;
    }

    if (b->alpha <= 0.0) {
        /* Perfect mirror. f is a delta; return f such that f*cos/pdf gives the
         * Fresnel term exactly, with pdf reported as 1 by convention. */
        vec3 wi = v3(-wo.x, -wo.y, wo.z);
        *wi_out = wi;
        *pdf_out = 1.0;
        Spectrum fr;
        conductor_fresnel_spectrum(b, wi.z, &fr);
        *f_out = ls_spectrum_scale(fr, 1.0 / ls_max(1e-12, wi.z));
        return true;
    }

    vec3 m = sample_ggx_vndf(wo, b->alpha, u1, u2);
    ls_real dot_om = v3dot(wo, m);
    if (dot_om <= 0.0) return false;
    vec3 wi = v3sub(v3scale(m, 2.0 * dot_om), wo);   /* reflect wo about m */
    if (wi.z <= 0.0) return false;

    *wi_out = wi;
    *pdf_out = ls_bsdf_pdf(b, wo, wi);
    if (*pdf_out <= 0.0) return false;
    ls_bsdf_eval(b, wo, wi, f_out);
    return true;
}

/* ---- spectral complex IOR presets (Rakic et al., abridged) ---- */
static void metal_from_table(const ls_real *lam, const ls_real *n, const ls_real *k,
                             int count, Spectrum *eta, Spectrum *kappa) {
    *eta   = ls_spectrum_from_samples(lam, n, count);
    *kappa = ls_spectrum_from_samples(lam, k, count);
}

void ls_metal_aluminium(Spectrum *eta, Spectrum *kappa) {
    static const ls_real lam[] = {400,450,500,550,600,650,700,750,800};
    static const ls_real n[]   = {0.490,0.618,0.769,0.958,1.200,1.470,1.830,2.400,2.750};
    static const ls_real k[]   = {4.860,5.470,6.080,6.690,7.260,7.790,8.310,8.620,8.310};
    metal_from_table(lam, n, k, 9, eta, kappa);
}

void ls_metal_copper(Spectrum *eta, Spectrum *kappa) {
    static const ls_real lam[] = {400,450,500,550,600,650,700,750,800};
    static const ls_real n[]   = {1.180,1.180,1.120,0.826,0.468,0.243,0.214,0.223,0.260};
    static const ls_real k[]   = {2.210,2.210,2.600,2.600,2.810,3.310,3.750,4.140,4.550};
    metal_from_table(lam, n, k, 9, eta, kappa);
}

void ls_metal_gold(Spectrum *eta, Spectrum *kappa) {
    static const ls_real lam[] = {400,450,500,550,600,650,700,750,800};
    static const ls_real n[]   = {1.658,1.578,0.849,0.371,0.242,0.192,0.166,0.157,0.153};
    static const ls_real k[]   = {1.956,1.867,1.871,2.399,2.966,3.470,3.929,4.357,4.759};
    metal_from_table(lam, n, k, 9, eta, kappa);
}
