#include "lightsim/scene.h"

/* The mesh a prim names, or NULL if it names none or names a hole. */
static const Mesh *mesh_of(const Scene *sc, const Prim *p) {
    if (p->kind != LS_PRIM_MESH) return NULL;
    if (p->mesh_id < 0 || p->mesh_id >= sc->nmeshes) return NULL;
    return sc->meshes[p->mesh_id];
}

bool ls_scene_intersect_ex(const Scene *sc, const Ray *ray, int skip_prim, Hit *hit) {
    Ray r = *ray;
    Hit best;
    bool found = false;
    for (int i = 0; i < sc->nprims; ++i) {
        if (i == skip_prim) continue;
        const Prim *p = &sc->prims[i];
        Hit h;
        bool got;
        if (p->kind == LS_PRIM_MESH) {
            const Mesh *m = mesh_of(sc, p);
            got = m && ls_mesh_intersect(m, p, &r, i, &h);
        } else {
            got = ls_prim_intersect(p, i, &r, &h);
        }
        if (got) {
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
    for (int i = 0; i < sc->nprims; ++i) {
        const Prim *pr = &sc->prims[i];        /* not `p`: that is the shading point */
        if (pr->kind == LS_PRIM_MESH) {
            const Mesh *m = mesh_of(sc, pr);
            if (m && ls_mesh_occludes(m, pr, &r)) return true;
        } else if (ls_prim_occludes(pr, &r)) return true;
    }
    return false;
}
