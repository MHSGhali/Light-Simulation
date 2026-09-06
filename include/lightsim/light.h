/* light.h — emitters.
 *
 * FLUX/SHAPE FACTORIZATION (the central anti-foot-gun):
 *   Every light stores its total RADIANT flux `phi_e` in watts, separately from
 *   a normalised spectral shape `s_hat` whose band integral is exactly 1.
 *   Spectral radiant intensity is therefore always phi_e * s_hat(lambda) * d(omega),
 *   and there is no other way to spell it.
 *
 *   A light specified in lumens is converted to watts by ls_watts_from_lumens()
 *   in the units layer BEFORE it reaches a constructor, so lumens are never
 *   stored and this file stays free of photometric constants entirely.
 *
 *   ls_light_finalize() re-derives the emitted flux from the light's geometry
 *   and asserts it matches phi_e. That catches the classic normalisation bugs
 *   (I0 = phi/4pi used for a spot, a two-sided area light emitting double)
 *   at scene-build time rather than as a plausible-looking wrong number.
 */
#ifndef LIGHTSIM_LIGHT_H
#define LIGHTSIM_LIGHT_H

#include "geom.h"
#include "spectrum.h"

typedef enum {
    LS_LIGHT_POINT,        /* isotropic delta source                        */
    LS_LIGHT_DIRECTIONAL,  /* delta direction, infinitely far (sun)         */
    LS_LIGHT_SPOT,         /* delta position, cone with smoothstep falloff  */
    LS_LIGHT_SPHERE,       /* uniform-radiance sphere                       */
    LS_LIGHT_DISK,         /* one-sided Lambertian disk                     */
    LS_LIGHT_RECT          /* one-sided Lambertian parallelogram            */
} LightKind;

typedef struct {
    LightKind kind;
    Spectrum  s_hat;       /* INVARIANT: ls_spectrum_integrate(&s_hat) == 1 */
    ls_real   phi_e;       /* total radiant flux, W (unused by DIRECTIONAL) */
    ls_real   e_perp;      /* DIRECTIONAL only: irradiance on a perpendicular
                            * surface, W/m^2 */
    vec3      p;           /* position / centre                             */
    vec3      n;           /* axis (spot, directional) or normal (disk/rect)*/
    vec3      ex, ey;      /* rect half-edge vectors                        */
    ls_real   radius;      /* sphere / disk                                 */
    ls_real   cos_total;   /* spot outer cone                               */
    ls_real   cos_falloff; /* spot inner cone; == cos_total for a hard edge  */
    ls_real   omega_eff;   /* spot: integral of the falloff over the sphere  */
    ls_real   area;        /* emitting area, m^2                            */
    ls_real   radiance;    /* area lights: uniform radiance scale, W/(m^2 sr)*/
    int       index;       /* column in the contribution matrix             */
} Light;

/* A sample of a light, taken from a shading point.
 *
 * `li_over_pdf` is the quantity such that a sample's contribution to
 * irradiance is EXACTLY  li_over_pdf * cos(theta_at_receiver).  For area
 * lights it is L / pdf_w; for delta lights it is I / r^2. Using one convention
 * for both means the estimator has a single code path.
 *
 * `pdf_w == 0` marks a delta light. Encoding it that way rather than as a bool
 * means a caller who forgets to check divides by zero (loud) instead of
 * computing a silently wrong MIS weight (quiet). */
typedef struct {
    vec3     wi;           /* unit, from the shading point toward the light  */
    ls_real  dist;         /* to the sampled point; INFINITY for directional */
    Spectrum li_over_pdf;  /* W/(m^2 sr nm)                                  */
    ls_real  pdf_w;        /* solid-angle pdf; EXACTLY 0 => delta light      */
    int      light_index;
} LightSample;

/* Constructors. `spd` is any non-negative spectrum; it is normalised to unit
 * band integral internally, so only its shape matters. */
Light ls_light_point(vec3 p, ls_real phi_e_w, Spectrum spd);
Light ls_light_directional(vec3 dir, ls_real e_perp, Spectrum spd);
Light ls_light_spot(vec3 p, vec3 dir, ls_real cone_total_rad,
                    ls_real cone_falloff_rad, ls_real phi_e_w, Spectrum spd);
Light ls_light_sphere(vec3 c, ls_real radius, ls_real phi_e_w, Spectrum spd);
Light ls_light_disk(vec3 c, vec3 n, ls_real radius, ls_real phi_e_w, Spectrum spd);
Light ls_light_rect(vec3 c, vec3 ex, vec3 ey, ls_real phi_e_w, Spectrum spd);

/* Normalise s_hat, derive area/radiance/omega_eff, and assert that the flux
 * implied by the geometry equals phi_e. Call once per light after construction. */
void ls_light_finalize(Light *l, int index);

/* Total radiant flux implied by the light's geometry and radiance. Used by
 * finalize's self-check and by the tests. */
ls_real ls_light_emitted_flux(const Light *l);

/* Radiant intensity in direction `w` (unit, pointing away from the light).
 * Defined for the delta kinds; area lights return their on-axis equivalent. */
ls_real ls_light_intensity(const Light *l, vec3 w);

/* Sample the light as seen from `p`. Returns false if the sample cannot
 * contribute (back face, degenerate geometry, outside a spot cone). */
bool ls_light_sample(const Light *l, vec3 p, ls_real u1, ls_real u2, LightSample *s);

/* Solid-angle pdf of having sampled the point `y` (with normal `ny`) on this
 * light from `ref`. Needed for the MIS weight applied to emission found by BSDF
 * sampling. Returns 0 for delta lights, which BSDF sampling can never hit. */
ls_real ls_light_pdf_w(const Light *l, vec3 ref, vec3 y, vec3 ny);

/* Emitted spectral radiance leaving this light in direction `w` (unit, away
 * from the surface with normal `ny`). Zero behind a one-sided emitter. */
Spectrum ls_light_radiance(const Light *l, vec3 ny, vec3 w);

#endif /* LIGHTSIM_LIGHT_H */
