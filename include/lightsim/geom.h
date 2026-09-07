/* geom.h — rays, intersections, and analytic primitives.
 *
 * M2 uses a linear scan over primitives. A BVH lands in M4, deliberately after
 * the transport physics is validated: an energy-conservation failure with no
 * acceleration structure present can only be a BSDF, PDF, or estimator bug.
 */
#ifndef LIGHTSIM_GEOM_H
#define LIGHTSIM_GEOM_H

#include "vec.h"

typedef struct {
    vec3   o, d;          /* d is unit length for analytic primitives */
    ls_real tmin, tmax;
} Ray;

typedef struct {
    ls_real t;
    vec3    p;
    vec3    ng;           /* GEOMETRIC normal, unit, in raw winding order.
                           * Never flipped by the intersector: flipping it
                           * silently is a whole family of energy-leak bugs.
                           * `backface` records the orientation instead. */
    int     prim_id;
    int     mat_id;
    int     light_id;     /* >= 0 iff this surface is an emitter */
    bool    backface;     /* dot(ray.d, ng) > 0 */
} Hit;

typedef enum {
    LS_PRIM_SPHERE,
    LS_PRIM_PLANE,        /* infinite plane through c with normal n */
    LS_PRIM_DISK,         /* radius r, in the plane through c with normal n */
    LS_PRIM_QUAD,         /* parallelogram: c +/- ex +/- ey (ex,ey half-edges) */
    LS_PRIM_MESH          /* triangle soup; the geometry lives in Scene.meshes,
                           * named by mesh_id. See mesh.h for why it cannot live
                           * in this struct. */
} PrimKind;

typedef struct {
    PrimKind kind;
    vec3     c;           /* centre */
    vec3     n;           /* normal (plane/disk/quad), unit */
    vec3     ex, ey;      /* quad half-edge vectors, perpendicular to n */
    ls_real  r;           /* sphere / disk radius */
    int      mat_id;
    int      light_id;
    /* MESH only: index into Scene.meshes, and -1 for every other kind. For a
     * mesh, c/n/ex/ey are reused as the placement -- c translates and
     * (ex, ey, n) is an ORTHONORMAL object-to-world basis -- so a mesh rotates
     * through exactly the code a quad already rotates through. */
    int      mesh_id;
} Prim;

/* Surface area of the primitive; 0 for an infinite plane. */
ls_real ls_prim_area(const Prim *p);

/* Nearest intersection in (ray->tmin, ray->tmax). Fills `hit` and returns true
 * on a hit. Does not modify the ray. */
bool ls_prim_intersect(const Prim *p, int prim_id, const Ray *ray, Hit *hit);

/* Any-hit test for shadow rays: returns on the first hit and ignores ordering. */
bool ls_prim_occludes(const Prim *p, const Ray *ray);

/* Offset a ray origin off a surface to avoid self-intersection. The offset is
 * along the GEOMETRIC normal and scaled by the point magnitude, so it holds up
 * from millimetre to kilometre scene scales. */
static inline vec3 ls_offset_origin(vec3 p, vec3 ng, vec3 dir) {
    ls_real s = (v3dot(ng, dir) < 0.0) ? -1.0 : 1.0;
    ls_real mag = ls_max(1.0, ls_max(fabs(p.x), ls_max(fabs(p.y), fabs(p.z))));
    return v3add(p, v3scale(ng, s * 1e-9 * mag));
}

#endif /* LIGHTSIM_GEOM_H */
