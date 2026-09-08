/* scene.h — primitive and light aggregation, plus visibility queries.
 * M2 intersects by linear scan; the BVH replaces ls_scene_intersect's body in M4
 * without changing this interface. */
#ifndef LIGHTSIM_SCENE_H
#define LIGHTSIM_SCENE_H

#include "geom.h"
#include "light.h"
#include "bsdf.h"
#include "mesh.h"

typedef struct {
    Bsdf     bsdf;
    Spectrum le;        /* emitted radiance W/(m^2 sr nm); zero if not emissive */
    bool     emissive;
    char     metal[8];  /* conductor preset ("al"/"cu"/"au"), for round-tripping */
} Material;

typedef struct {
    Prim     *prims;
    int       nprims;
    Material *mats;
    int       nmats;
    Light    *lights;
    int       nlights;
    /* Triangle geometry, named by Prim.mesh_id. An entry may be NULL: deleting
     * a mesh prim leaves a hole rather than compacting, so that mesh_id stays
     * stable for every undo snapshot still holding one. */
    Mesh    **meshes;
    int       nmeshes;
} Scene;

/* Nearest hit, optionally ignoring one primitive (`skip_prim < 0` ignores
 * nothing). THE traversal entry point -- the BVH replaces this body in M4, and
 * anything that walks the prim array itself has to be re-plumbed then, so
 * callers that need an exclusion take it here rather than rolling their own.
 * Returns false if the ray escapes. */
bool ls_scene_intersect_ex(const Scene *sc, const Ray *ray, int skip_prim, Hit *hit);

static inline bool ls_scene_intersect(const Scene *sc, const Ray *ray, Hit *hit) {
    return ls_scene_intersect_ex(sc, ray, -1, hit);
}

/* Is the segment from `p` (offset off the surface along `ng`) toward `wi`,
 * of length `dist`, blocked? `dist` may be HUGE_VAL for directional lights.
 * The segment is shortened slightly at the far end so that a light's own
 * geometry does not shadow it. */
bool ls_scene_occluded(const Scene *sc, vec3 p, vec3 ng, vec3 wi, ls_real dist);

#endif /* LIGHTSIM_SCENE_H */
