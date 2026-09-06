/* scene.h — primitive and light aggregation, plus visibility queries.
 * M2 intersects by linear scan; the BVH replaces ls_scene_intersect's body in M4
 * without changing this interface. */
#ifndef LIGHTSIM_SCENE_H
#define LIGHTSIM_SCENE_H

#include "geom.h"
#include "light.h"

typedef struct {
    Prim   *prims;
    int     nprims;
    Light  *lights;
    int     nlights;
} Scene;

/* Nearest hit. Returns false if the ray escapes. */
bool ls_scene_intersect(const Scene *sc, const Ray *ray, Hit *hit);

/* Is the segment from `p` (offset off the surface along `ng`) toward `wi`,
 * of length `dist`, blocked? `dist` may be HUGE_VAL for directional lights.
 * The segment is shortened slightly at the far end so that a light's own
 * geometry does not shadow it. */
bool ls_scene_occluded(const Scene *sc, vec3 p, vec3 ng, vec3 wi, ls_real dist);

#endif /* LIGHTSIM_SCENE_H */
