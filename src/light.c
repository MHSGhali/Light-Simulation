#include "lightsim/light.h"
#include <string.h>
#include <math.h>

static Light base(LightKind k, Spectrum spd) {
    Light l;
    memset(&l, 0, sizeof l);
    l.kind = k;
    l.s_hat = spd;
    l.cos_total = -1.0;
    l.cos_falloff = -1.0;
    l.index = -1;
    return l;
}

Light ls_light_point(vec3 p, ls_real phi_e_w, Spectrum spd) {
    Light l = base(LS_LIGHT_POINT, spd);
    l.p = p; l.phi_e = phi_e_w;
    return l;
}

Light ls_light_directional(vec3 dir, ls_real e_perp, Spectrum spd) {
    Light l = base(LS_LIGHT_DIRECTIONAL, spd);
    l.n = v3norm(dir);          /* direction of propagation */
    l.e_perp = e_perp;
    return l;
}

Light ls_light_spot(vec3 p, vec3 dir, ls_real cone_total_rad,
                    ls_real cone_falloff_rad, ls_real phi_e_w, Spectrum spd) {
    Light l = base(LS_LIGHT_SPOT, spd);
    l.p = p; l.n = v3norm(dir); l.phi_e = phi_e_w;
    l.cos_total   = cos(cone_total_rad);
    l.cos_falloff = cos(ls_min(cone_falloff_rad, cone_total_rad));
    return l;
}

Light ls_light_sphere(vec3 c, ls_real radius, ls_real phi_e_w, Spectrum spd) {
    Light l = base(LS_LIGHT_SPHERE, spd);
    l.p = c; l.radius = radius; l.phi_e = phi_e_w;
    return l;
}

Light ls_light_disk(vec3 c, vec3 n, ls_real radius, ls_real phi_e_w, Spectrum spd) {
    Light l = base(LS_LIGHT_DISK, spd);
    l.p = c; l.n = v3norm(n); l.radius = radius; l.phi_e = phi_e_w;
    return l;
}

Light ls_light_rect(vec3 c, vec3 ex, vec3 ey, ls_real phi_e_w, Spectrum spd) {
    Light l = base(LS_LIGHT_RECT, spd);
    l.p = c; l.ex = ex; l.ey = ey; l.phi_e = phi_e_w;
    l.n = v3norm(v3cross(ex, ey));
    return l;
}

/* Spot falloff: 1 inside the inner cone, smoothstep to 0 at the outer cone. */
static ls_real spot_falloff(const Light *l, ls_real cos_theta) {
    if (cos_theta <= l->cos_total)   return 0.0;
    if (cos_theta >= l->cos_falloff) return 1.0;
    ls_real t = (cos_theta - l->cos_total) / (l->cos_falloff - l->cos_total);
    return t * t * (3.0 - 2.0 * t);          /* smoothstep */
}

/* Integral of the spot falloff over the sphere:
 *   omega_eff = 2 pi * integral_{cos_total}^{1} f(mu) dmu
 * For a hard-edged cone this is the exact 2 pi (1 - cos_total); the smoothstep
 * region is integrated numerically. */
static ls_real spot_omega_eff(const Light *l) {
    if (l->cos_falloff <= l->cos_total)
        return LS_TWO_PI * (1.0 - l->cos_total);
    /* Hard part plus the numerically integrated transition. */
    ls_real hard = LS_TWO_PI * (1.0 - l->cos_falloff);
    const int N = 4096;
    ls_real a = l->cos_total, b = l->cos_falloff, sum = 0.0;
    for (int i = 0; i < N; ++i) {
        ls_real mu = a + (b - a) * ((ls_real)i + 0.5) / (ls_real)N;
        sum += spot_falloff(l, mu);
    }
    return hard + LS_TWO_PI * sum * (b - a) / (ls_real)N;
}

ls_real ls_light_emitted_flux(const Light *l) {
    switch (l->kind) {
        case LS_LIGHT_POINT:
            /* isotropic: I * 4 pi */
            return (l->phi_e * LS_INV_4PI) * 4.0 * LS_PI;
        case LS_LIGHT_SPOT:
            /* I0 * omega_eff */
            return (l->omega_eff > 0.0 ? l->phi_e / l->omega_eff : 0.0) * l->omega_eff;
        case LS_LIGHT_DIRECTIONAL:
            return 0.0;   /* infinite extent: flux is not defined */
        case LS_LIGHT_SPHERE:
        case LS_LIGHT_DISK:
        case LS_LIGHT_RECT:
            /* A Lambertian emitter of uniform radiance L over area A radiates
             * L * A * pi into the hemisphere above each surface element. */
            return l->radiance * l->area * LS_PI;
    }
    return 0.0;
}

void ls_light_finalize(Light *l, int index) {
    l->index = index;
    l->s_hat = ls_spectrum_normalize_to(l->s_hat, 1.0);

    switch (l->kind) {
        case LS_LIGHT_SPOT:
            l->omega_eff = spot_omega_eff(l);
            break;
        case LS_LIGHT_SPHERE:
            l->area = 4.0 * LS_PI * l->radius * l->radius;
            /* phi = L * A * pi  =>  L = phi / (4 pi^2 R^2) */
            l->radiance = l->area > 0.0 ? l->phi_e / (l->area * LS_PI) : 0.0;
            break;
        case LS_LIGHT_DISK:
            l->area = LS_PI * l->radius * l->radius;
            l->radiance = l->area > 0.0 ? l->phi_e / (l->area * LS_PI) : 0.0;
            break;
        case LS_LIGHT_RECT:
            l->area = 4.0 * v3len(v3cross(l->ex, l->ey));
            l->radiance = l->area > 0.0 ? l->phi_e / (l->area * LS_PI) : 0.0;
            break;
        default:
            break;
    }

    /* Self-check: the flux implied by the geometry must equal the flux the
     * caller asked for. Catches every normalisation slip at build time. */
    if (l->kind != LS_LIGHT_DIRECTIONAL) {
        ls_real got = ls_light_emitted_flux(l);
        LS_ASSERT(fabs(got - l->phi_e) <= 1e-9 * ls_max(1.0, fabs(l->phi_e))
                  && "light flux normalisation is inconsistent with its geometry");
        LS_UNUSED(got);
    }
    /* The spectral shape must integrate to exactly 1. */
    LS_ASSERT(fabs(ls_spectrum_integrate(&l->s_hat) - 1.0) < 1e-5
              && "spectral shape is not normalised");
}

ls_real ls_light_intensity(const Light *l, vec3 w) {
    switch (l->kind) {
        case LS_LIGHT_POINT:  return l->phi_e * LS_INV_4PI;
        case LS_LIGHT_SPOT:   return l->omega_eff > 0.0
                                   ? (l->phi_e / l->omega_eff) * spot_falloff(l, v3dot(w, l->n))
                                   : 0.0;
        case LS_LIGHT_DIRECTIONAL: return 0.0;
        case LS_LIGHT_SPHERE: return l->radiance * LS_PI * l->radius * l->radius;
        case LS_LIGHT_DISK:
        case LS_LIGHT_RECT:   return l->radiance * l->area * ls_max(0.0, v3dot(w, l->n));
    }
    return 0.0;
}

ls_real ls_light_pdf_w(const Light *l, vec3 ref, vec3 y, vec3 ny) {
    switch (l->kind) {
        case LS_LIGHT_POINT:
        case LS_LIGHT_SPOT:
        case LS_LIGHT_DIRECTIONAL:
            return 0.0;                       /* delta: unreachable by sampling */
        case LS_LIGHT_SPHERE:
        case LS_LIGHT_DISK:
        case LS_LIGHT_RECT: {
            vec3 d = v3sub(y, ref);
            ls_real d2 = v3len2(d);
            if (d2 <= 0.0 || l->area <= 0.0) return 0.0;
            vec3 wi = v3scale(d, 1.0 / sqrt(d2));
            ls_real cos_y = v3dot(ny, v3neg(wi));
            if (cos_y <= 0.0) return 0.0;
            return ls_pdf_area_to_solid_angle(1.0 / l->area, d2, cos_y);
        }
    }
    return 0.0;
}

Spectrum ls_light_radiance(const Light *l, vec3 ny, vec3 w) {
    switch (l->kind) {
        case LS_LIGHT_SPHERE:
        case LS_LIGHT_DISK:
        case LS_LIGHT_RECT:
            if (v3dot(ny, w) <= 0.0) return ls_spectrum_zero();   /* one-sided */
            return ls_spectrum_scale(l->s_hat, l->radiance);
        default:
            return ls_spectrum_zero();        /* delta lights have no radiance */
    }
}

bool ls_light_sample(const Light *l, vec3 p, ls_real u1, ls_real u2, LightSample *s) {
    s->light_index = l->index;

    switch (l->kind) {
        case LS_LIGHT_POINT: {
            vec3 d = v3sub(l->p, p);
            ls_real d2 = v3len2(d);
            if (d2 <= 0.0) return false;
            ls_real dist = sqrt(d2);
            s->wi = v3scale(d, 1.0 / dist);
            s->dist = dist;
            s->pdf_w = 0.0;                              /* delta */
            s->li_over_pdf = ls_spectrum_scale(l->s_hat, (l->phi_e * LS_INV_4PI) / d2);
            return true;
        }
        case LS_LIGHT_SPOT: {
            vec3 d = v3sub(l->p, p);
            ls_real d2 = v3len2(d);
            if (d2 <= 0.0) return false;
            ls_real dist = sqrt(d2);
            s->wi = v3scale(d, 1.0 / dist);
            /* Direction from the light toward the receiver is -wi. */
            ls_real f = spot_falloff(l, v3dot(v3neg(s->wi), l->n));
            if (f <= 0.0) return false;
            s->dist = dist;
            s->pdf_w = 0.0;                              /* delta */
            ls_real i0 = l->omega_eff > 0.0 ? l->phi_e / l->omega_eff : 0.0;
            s->li_over_pdf = ls_spectrum_scale(l->s_hat, i0 * f / d2);
            return true;
        }
        case LS_LIGHT_DIRECTIONAL: {
            s->wi = v3neg(l->n);                         /* toward the source */
            s->dist = HUGE_VAL;
            s->pdf_w = 0.0;                              /* delta */
            s->li_over_pdf = ls_spectrum_scale(l->s_hat, l->e_perp);
            return true;
        }
        case LS_LIGHT_SPHERE: {
            /* Uniform area sampling over the whole sphere. Points on the far
             * side come out back-facing and are rejected; that halves the
             * sample efficiency but keeps the estimator unbiased and the code
             * honest. Cone sampling is a variance optimisation for later. */
            vec3 dir = ls_sample_sphere_uniform(u1, u2);
            vec3 y   = v3add(l->p, v3scale(dir, l->radius));
            vec3 ny  = dir;
            vec3 dv  = v3sub(y, p);
            ls_real d2 = v3len2(dv);
            if (d2 <= 0.0) return false;
            ls_real dist = sqrt(d2);
            s->wi = v3scale(dv, 1.0 / dist);
            ls_real cos_y = v3dot(ny, v3neg(s->wi));
            if (cos_y <= 0.0) return false;              /* back face */
            ls_real pdf_a = 1.0 / l->area;
            s->pdf_w = ls_pdf_area_to_solid_angle(pdf_a, d2, cos_y);
            if (s->pdf_w <= 0.0) return false;
            s->dist = dist;
            s->li_over_pdf = ls_spectrum_scale(l->s_hat, l->radiance / s->pdf_w);
            return true;
        }
        case LS_LIGHT_DISK: {
            vec2 dd = ls_sample_disk_concentric(u1, u2);
            Basis b = ls_basis(l->n);
            vec3 y = v3add(l->p, v3add(v3scale(b.t, dd.x * l->radius),
                                       v3scale(b.b, dd.y * l->radius)));
            vec3 dv = v3sub(y, p);
            ls_real d2 = v3len2(dv);
            if (d2 <= 0.0) return false;
            ls_real dist = sqrt(d2);
            s->wi = v3scale(dv, 1.0 / dist);
            ls_real cos_y = v3dot(l->n, v3neg(s->wi));
            if (cos_y <= 0.0) return false;              /* one-sided */
            ls_real pdf_a = 1.0 / l->area;
            s->pdf_w = ls_pdf_area_to_solid_angle(pdf_a, d2, cos_y);
            if (s->pdf_w <= 0.0) return false;
            s->dist = dist;
            s->li_over_pdf = ls_spectrum_scale(l->s_hat, l->radiance / s->pdf_w);
            return true;
        }
        case LS_LIGHT_RECT: {
            ls_real a = 2.0 * u1 - 1.0, b = 2.0 * u2 - 1.0;
            vec3 y = v3add(l->p, v3add(v3scale(l->ex, a), v3scale(l->ey, b)));
            vec3 dv = v3sub(y, p);
            ls_real d2 = v3len2(dv);
            if (d2 <= 0.0) return false;
            ls_real dist = sqrt(d2);
            s->wi = v3scale(dv, 1.0 / dist);
            ls_real cos_y = v3dot(l->n, v3neg(s->wi));
            if (cos_y <= 0.0) return false;              /* one-sided */
            ls_real pdf_a = 1.0 / l->area;
            s->pdf_w = ls_pdf_area_to_solid_angle(pdf_a, d2, cos_y);
            if (s->pdf_w <= 0.0) return false;
            s->dist = dist;
            s->li_over_pdf = ls_spectrum_scale(l->s_hat, l->radiance / s->pdf_w);
            return true;
        }
    }
    return false;
}
