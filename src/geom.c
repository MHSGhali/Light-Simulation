#include "lightsim/geom.h"

ls_real ls_prim_area(const Prim *p) {
    switch (p->kind) {
        case LS_PRIM_SPHERE: return 4.0 * LS_PI * p->r * p->r;
        case LS_PRIM_DISK:   return LS_PI * p->r * p->r;
        case LS_PRIM_QUAD:   return 4.0 * v3len(v3cross(p->ex, p->ey));
        case LS_PRIM_PLANE:  return 0.0;   /* infinite */
    }
    return 0.0;
}

static bool hit_sphere(const Prim *p, const Ray *ray, ls_real *t_out, vec3 *n_out) {
    vec3 oc = v3sub(ray->o, p->c);
    ls_real b = v3dot(oc, ray->d);
    ls_real c = v3dot(oc, oc) - p->r * p->r;
    ls_real disc = b * b - c;
    if (disc < 0.0) return false;
    ls_real sq = sqrt(disc);
    ls_real t = -b - sq;
    if (t <= ray->tmin || t >= ray->tmax) {
        t = -b + sq;
        if (t <= ray->tmin || t >= ray->tmax) return false;
    }
    *t_out = t;
    *n_out = v3scale(v3sub(v3add(ray->o, v3scale(ray->d, t)), p->c), 1.0 / p->r);
    return true;
}

static bool hit_plane_t(const Prim *p, const Ray *ray, ls_real *t_out) {
    ls_real denom = v3dot(ray->d, p->n);
    if (fabs(denom) < 1e-12) return false;          /* parallel */
    ls_real t = v3dot(v3sub(p->c, ray->o), p->n) / denom;
    if (t <= ray->tmin || t >= ray->tmax) return false;
    *t_out = t;
    return true;
}

static bool hit_disk(const Prim *p, const Ray *ray, ls_real *t_out) {
    ls_real t;
    if (!hit_plane_t(p, ray, &t)) return false;
    vec3 q = v3sub(v3add(ray->o, v3scale(ray->d, t)), p->c);
    if (v3len2(q) > p->r * p->r) return false;
    *t_out = t;
    return true;
}

static bool hit_quad(const Prim *p, const Ray *ray, ls_real *t_out) {
    ls_real t;
    if (!hit_plane_t(p, ray, &t)) return false;
    vec3 q = v3sub(v3add(ray->o, v3scale(ray->d, t)), p->c);
    /* Project onto the half-edge vectors; |coord| <= 1 is inside. */
    ls_real ex2 = v3len2(p->ex), ey2 = v3len2(p->ey);
    if (ex2 <= 0.0 || ey2 <= 0.0) return false;
    ls_real a = v3dot(q, p->ex) / ex2;
    ls_real b = v3dot(q, p->ey) / ey2;
    if (fabs(a) > 1.0 || fabs(b) > 1.0) return false;
    *t_out = t;
    return true;
}

bool ls_prim_intersect(const Prim *p, int prim_id, const Ray *ray, Hit *hit) {
    ls_real t = 0.0;
    vec3 n = p->n;
    switch (p->kind) {
        case LS_PRIM_SPHERE: if (!hit_sphere(p, ray, &t, &n)) return false; break;
        case LS_PRIM_PLANE:  if (!hit_plane_t(p, ray, &t))    return false; break;
        case LS_PRIM_DISK:   if (!hit_disk(p, ray, &t))       return false; break;
        case LS_PRIM_QUAD:   if (!hit_quad(p, ray, &t))       return false; break;
        default: return false;
    }
    hit->t        = t;
    hit->p        = v3add(ray->o, v3scale(ray->d, t));
    hit->ng       = n;
    hit->prim_id  = prim_id;
    hit->mat_id   = p->mat_id;
    hit->light_id = p->light_id;
    hit->backface = v3dot(ray->d, n) > 0.0;
    return true;
}

bool ls_prim_occludes(const Prim *p, const Ray *ray) {
    ls_real t;
    vec3 n;
    switch (p->kind) {
        case LS_PRIM_SPHERE: return hit_sphere(p, ray, &t, &n);
        case LS_PRIM_PLANE:  return hit_plane_t(p, ray, &t);
        case LS_PRIM_DISK:   return hit_disk(p, ray, &t);
        case LS_PRIM_QUAD:   return hit_quad(p, ray, &t);
    }
    return false;
}
