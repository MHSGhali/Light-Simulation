/* vec.h — 3D vector / matrix maths and the sampling warps used by the
 * estimators. All warps document the measure their PDF is expressed in;
 * confusing solid angle with projected solid angle is the classic way to get a
 * renderer that looks right and integrates wrong. */
#ifndef LIGHTSIM_VEC_H
#define LIGHTSIM_VEC_H

#include "core.h"
#include <math.h>

typedef struct { ls_real x, y, z; } vec3;
typedef struct { ls_real x, y; } vec2;

static inline vec3 v3(ls_real x, ls_real y, ls_real z) { vec3 r = {x,y,z}; return r; }
static inline vec3 v3add(vec3 a, vec3 b) { return v3(a.x+b.x, a.y+b.y, a.z+b.z); }
static inline vec3 v3sub(vec3 a, vec3 b) { return v3(a.x-b.x, a.y-b.y, a.z-b.z); }
static inline vec3 v3mul(vec3 a, vec3 b) { return v3(a.x*b.x, a.y*b.y, a.z*b.z); }
static inline vec3 v3scale(vec3 a, ls_real s) { return v3(a.x*s, a.y*s, a.z*s); }
static inline vec3 v3neg(vec3 a) { return v3(-a.x, -a.y, -a.z); }
static inline ls_real v3dot(vec3 a, vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline vec3 v3cross(vec3 a, vec3 b) {
    return v3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}
static inline ls_real v3len2(vec3 a) { return v3dot(a, a); }
static inline ls_real v3len(vec3 a) { return sqrt(v3dot(a, a)); }
static inline vec3 v3norm(vec3 a) {
    ls_real l = v3len(a);
    return l > 0.0 ? v3scale(a, 1.0 / l) : a;
}
static inline ls_real v3dist(vec3 a, vec3 b) { return v3len(v3sub(a, b)); }
static inline vec3 v3lerp(ls_real t, vec3 a, vec3 b) { return v3add(a, v3scale(v3sub(b,a), t)); }
static inline ls_real v3maxc(vec3 a) { return ls_max(a.x, ls_max(a.y, a.z)); }
static inline ls_real v3minc(vec3 a) { return ls_min(a.x, ls_min(a.y, a.z)); }

/* Rotate `v` about a unit `axis` by `angle` radians (Rodrigues). Right-handed:
 * positive angle turns counter-clockwise looking down the axis toward the
 * origin. */
static inline vec3 v3rotate(vec3 v, vec3 axis, ls_real angle) {
    ls_real c = cos(angle), s = sin(angle);
    return v3add(v3add(v3scale(v, c), v3scale(v3cross(axis, v), s)),
                 v3scale(axis, v3dot(axis, v) * (1.0 - c)));
}

/* Closest point on the line (p, dir) to the line (o, d), returned as the
 * parameter t along the first. Returns false when the lines are near-parallel,
 * where the closest point is not well defined. Used to drag along an axis
 * handle: the axis is one line, the pick ray the other. */
static inline bool ls_line_closest_t(vec3 p, vec3 dir, vec3 o, vec3 d,
                                     ls_real *t_out) {
    ls_real b = v3dot(dir, d);
    ls_real denom = 1.0 - b * b;
    if (fabs(denom) < 1e-7) return false;
    vec3 w0 = v3sub(p, o);
    ls_real dd = v3dot(dir, w0), e = v3dot(d, w0);
    *t_out = (b * e - dd) / denom;
    return true;
}

/* Intersect the ray (o, d) with the plane through `p` with unit normal `n`.
 * Returns false when the ray is near-parallel to the plane, where the hit is
 * numerically useless. Used to drag a rotation ring. */
static inline bool ls_ray_plane(vec3 o, vec3 d, vec3 p, vec3 n, vec3 *hit) {
    ls_real dn = v3dot(d, n);
    if (fabs(dn) < 1e-6) return false;
    ls_real t = v3dot(v3sub(p, o), n) / dn;
    if (t <= 0.0) return false;
    *hit = v3add(o, v3scale(d, t));
    return true;
}

/* Orthonormal basis around a unit normal (Duff et al., branchless, stable). */
typedef struct { vec3 t, b, n; } Basis;
static inline Basis ls_basis(vec3 n) {
    Basis o;
    o.n = n;
    ls_real sign = copysign(1.0, n.z);
    ls_real a = -1.0 / (sign + n.z);
    ls_real b = n.x * n.y * a;
    o.t = v3(1.0 + sign * n.x * n.x * a, sign * b, -sign * n.x);
    o.b = v3(b, sign + n.y * n.y * a, -n.y);
    return o;
}
static inline vec3 ls_basis_to_world(Basis o, vec3 v) {
    return v3add(v3add(v3scale(o.t, v.x), v3scale(o.b, v.y)), v3scale(o.n, v.z));
}
static inline vec3 ls_basis_to_local(Basis o, vec3 v) {
    return v3(v3dot(v, o.t), v3dot(v, o.b), v3dot(v, o.n));
}

/* ---- sampling warps ---- */

/* Concentric (Shirley-Chiu) unit-disk map: low distortion, area-uniform. */
static inline vec2 ls_sample_disk_concentric(ls_real u1, ls_real u2) {
    ls_real ox = 2.0 * u1 - 1.0, oy = 2.0 * u2 - 1.0;
    vec2 r = { 0.0, 0.0 };
    if (ox == 0.0 && oy == 0.0) return r;
    ls_real rad, theta;
    if (fabs(ox) > fabs(oy)) { rad = ox; theta = (LS_PI / 4.0) * (oy / ox); }
    else                     { rad = oy; theta = (LS_PI / 2.0) - (LS_PI / 4.0) * (ox / oy); }
    r.x = rad * cos(theta);
    r.y = rad * sin(theta);
    return r;
}

/* Cosine-weighted hemisphere about +z. PDF is cos(theta)/pi in SOLID ANGLE. */
static inline vec3 ls_sample_hemisphere_cosine(ls_real u1, ls_real u2) {
    vec2 d = ls_sample_disk_concentric(u1, u2);
    ls_real z = sqrt(ls_max(0.0, 1.0 - d.x * d.x - d.y * d.y));
    return v3(d.x, d.y, z);
}
static inline ls_real ls_pdf_hemisphere_cosine(ls_real cos_theta) {
    return cos_theta > 0.0 ? cos_theta * LS_INV_PI : 0.0;
}

/* Uniform on the unit sphere. PDF 1/(4 pi) in solid angle. */
static inline vec3 ls_sample_sphere_uniform(ls_real u1, ls_real u2) {
    ls_real z = 1.0 - 2.0 * u1;
    ls_real r = sqrt(ls_max(0.0, 1.0 - z * z));
    ls_real phi = LS_TWO_PI * u2;
    return v3(r * cos(phi), r * sin(phi), z);
}

/* Uniform inside a cone of half-angle acos(cos_max) about +z.
 * PDF 1/(2 pi (1 - cos_max)) in solid angle. */
static inline vec3 ls_sample_cone_uniform(ls_real u1, ls_real u2, ls_real cos_max) {
    ls_real cos_theta = (1.0 - u1) + u1 * cos_max;
    ls_real sin_theta = sqrt(ls_max(0.0, 1.0 - cos_theta * cos_theta));
    ls_real phi = LS_TWO_PI * u2;
    return v3(sin_theta * cos(phi), sin_theta * sin(phi), cos_theta);
}
static inline ls_real ls_pdf_cone_uniform(ls_real cos_max) {
    return 1.0 / (LS_TWO_PI * (1.0 - cos_max));
}

/* Uniform barycentric coordinates on a triangle. */
static inline vec2 ls_sample_triangle_uniform(ls_real u1, ls_real u2) {
    ls_real su = sqrt(u1);
    vec2 r = { 1.0 - su, u2 * su };
    return r;
}

/* Solid angle of a cone of half-angle alpha: 2 pi (1 - cos alpha).
 * This is the relation used to convert a source's angular span into steradians. */
static inline ls_real ls_cone_solid_angle(ls_real half_angle_rad) {
    return LS_TWO_PI * (1.0 - cos(half_angle_rad));
}

/* Convert a PDF from area measure on a light to solid angle measure at the
 * shading point. THE single place this Jacobian is applied.
 *   pdf_omega = pdf_area * d^2 / |cos(theta_light)|
 * Returns 0 when the light element is edge-on (the conversion is singular). */
static inline ls_real ls_pdf_area_to_solid_angle(ls_real pdf_area, ls_real dist2,
                                                 ls_real cos_theta_light) {
    ls_real c = fabs(cos_theta_light);
    if (c <= 1e-9) return 0.0;
    return pdf_area * dist2 / c;
}

#endif /* LIGHTSIM_VEC_H */
