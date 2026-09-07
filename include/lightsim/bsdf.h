/* bsdf.h — surface scattering.
 *
 * CONVENTIONS (stated once, relied on everywhere):
 *
 *   - All directions are in the LOCAL shading frame, +z along the normal, and
 *     point AWAY from the surface.
 *   - ls_bsdf_eval returns f, the BRDF, in units of 1/sr. It does NOT include
 *     the cosine factor. The cosine belongs to the estimator, and appears there
 *     exactly twice (once in NEE, once in the path continuation).
 *   - ls_bsdf_pdf is ALWAYS in solid-angle measure, 1/sr. Never projected
 *     solid angle.
 *   - `alpha` is the GGX width parameter, stored directly. The alpha =
 *     roughness^2 remapping is an artist convention and belongs in a scene
 *     parser, not in the physics.
 */
#ifndef LIGHTSIM_BSDF_H
#define LIGHTSIM_BSDF_H

#include "spectrum.h"
#include "vec.h"

typedef enum {
    LS_BSDF_LAMBERT,     /* ideal diffuse                                   */
    LS_BSDF_CONDUCTOR    /* GGX microfacet metal, spectral complex IOR n+ik */
} BsdfKind;

typedef struct {
    BsdfKind kind;
    Spectrum rho;        /* LAMBERT: diffuse reflectance, dimensionless [0,1] */
    Spectrum eta;        /* CONDUCTOR: real part of the IOR, n(lambda)        */
    Spectrum kappa;      /* CONDUCTOR: extinction coefficient, k(lambda)      */
    ls_real  alpha;      /* GGX width. 0 => perfect specular (delta)          */
} Bsdf;

/* Fresnel reflectance of a conductor at normal-to-surface angle cos_theta,
 * for a single wavelength with complex IOR n + ik. Unpolarised average. */
ls_real ls_fresnel_conductor(ls_real cos_theta, ls_real n, ls_real k);

/* Fresnel reflectance of a dielectric interface, real IOR ratio eta. */
ls_real ls_fresnel_dielectric(ls_real cos_theta, ls_real eta);

/* GGX normal distribution and height-correlated Smith masking-shadowing. */
ls_real ls_ggx_d(ls_real cos_theta_m, ls_real alpha);
ls_real ls_ggx_g2(vec3 wo, vec3 wi, ls_real alpha);
ls_real ls_ggx_g1(vec3 w, ls_real alpha);

bool    ls_bsdf_is_delta(const Bsdf *b);
/* f, in 1/sr, WITHOUT the cosine. */
void    ls_bsdf_eval(const Bsdf *b, vec3 wo, vec3 wi, Spectrum *f_out);
/* Solid-angle pdf of sampling wi given wo. Zero for delta lobes. */
ls_real ls_bsdf_pdf(const Bsdf *b, vec3 wo, vec3 wi);
/* Draw a direction. Returns false if the sample is degenerate. */
bool    ls_bsdf_sample(const Bsdf *b, vec3 wo, ls_real u1, ls_real u2,
                       vec3 *wi_out, Spectrum *f_out, ls_real *pdf_out);

/* Spectral complex IOR presets for common metals, resampled onto the bin grid. */
void ls_metal_aluminium(Spectrum *eta, Spectrum *kappa);
void ls_metal_copper(Spectrum *eta, Spectrum *kappa);
void ls_metal_gold(Spectrum *eta, Spectrum *kappa);

#endif /* LIGHTSIM_BSDF_H */
