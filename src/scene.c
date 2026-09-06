#include "lightsim/scene.h"

bool ls_scene_intersect(const Scene *sc, const Ray *ray, Hit *hit) {
    Ray r = *ray;
    Hit best;
    bool found = false;
    for (int i = 0; i < sc->nprims; ++i) {
        Hit h;
        if (ls_prim_intersect(&sc->prims[i], i, &r, &h)) {
            r.tmax = h.t;         /* shrink so later prims must beat it */
            best = h;
            found = true;
        }
    }
    if (found) *hit = best;
    return found;
}

bool ls_scene_occluded(const Scene *sc, vec3 p, vec3 ng, vec3 wi, ls_real dist) {
    Ray r;
    r.o = ls_offset_origin(p, ng, wi);
    r.d = wi;
    r.tmin = 0.0;
    /* Stop just short of the light so its own surface does not occlude it. */
    r.tmax = (dist == HUGE_VAL) ? HUGE_VAL : dist * (1.0 - 1e-6);
    for (int i = 0; i < sc->nprims; ++i)
        if (ls_prim_occludes(&sc->prims[i], &r)) return true;
    return false;
}
