#include "lightsim/sceneedit.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

void ls_scene_rebuild(SceneDesc *d) {
    d->scene.prims  = d->prims;   d->scene.nprims  = d->nprims;
    d->scene.mats   = d->mats;    d->scene.nmats   = d->nmats;
    d->scene.lights = d->lights;  d->scene.nlights = d->nlights;
}

/* ------------------------------------------------------------- growth ---- */

static bool grow_prims(SceneDesc *d) {
    if (d->nprims < d->cap_prims) return true;
    int cap = d->cap_prims ? d->cap_prims * 2 : 8;
    Prim *p = realloc(d->prims, (size_t)cap * sizeof *p);
    if (!p) return false;
    d->prims = p; d->cap_prims = cap;
    return true;
}

static bool grow_lights(SceneDesc *d) {
    if (d->nlights < d->cap_lights) return true;
    int cap = d->cap_lights ? d->cap_lights * 2 : 8;
    Light *l = realloc(d->lights, (size_t)cap * sizeof *l);
    if (!l) return false;
    d->lights = l; d->cap_lights = cap;
    return true;
}

/* Materials and their names are parallel arrays sharing one capacity, so they
 * must grow together -- growing them separately is how the parser once left the
 * name array unallocated. */
static bool grow_mats(SceneDesc *d) {
    if (d->nmats < d->cap_mats) return true;
    int cap = d->cap_mats ? d->cap_mats * 2 : 8;
    Material *m = realloc(d->mats, (size_t)cap * sizeof *m);
    if (!m) return false;
    d->mats = m;
    void *n = realloc(d->names, (size_t)cap * sizeof *d->names);
    if (!n) return false;
    d->names = n;
    d->cap_mats = cap;
    return true;
}

/* ------------------------------------------------------------- queries --- */

int ls_scene_light_prim(const SceneDesc *d, int index) {
    for (int i = 0; i < d->nprims; ++i)
        if (d->prims[i].light_id == index) return i;
    return -1;
}

static bool light_has_geometry(const Light *l) {
    return l->kind == LS_LIGHT_RECT || l->kind == LS_LIGHT_DISK
        || l->kind == LS_LIGHT_SPHERE;
}

/* ------------------------------------------------------------ mutators --- */

int ls_scene_add_material(SceneDesc *d, Material m, const char *name) {
    if (!grow_mats(d)) return -1;
    int id = d->nmats;
    d->mats[id] = m;
    snprintf(d->names[id], sizeof d->names[id], "%s", name ? name : "mat");
    d->nmats++;
    ls_scene_rebuild(d);
    return id;
}

int ls_scene_add_prim(SceneDesc *d, Prim p) {
    if (!grow_prims(d)) return -1;
    int id = d->nprims;
    d->prims[id] = p;
    d->nprims++;
    ls_scene_rebuild(d);
    return id;
}

void ls_scene_remove_prim(SceneDesc *d, int index) {
    if (index < 0 || index >= d->nprims) return;
    for (int i = index; i + 1 < d->nprims; ++i) d->prims[i] = d->prims[i + 1];
    d->nprims--;
    ls_scene_rebuild(d);
}

void ls_scene_sync_light_geom(SceneDesc *d, int index) {
    if (index < 0 || index >= d->nlights) return;
    Light *l = &d->lights[index];
    int pi = ls_scene_light_prim(d, index);

    if (!light_has_geometry(l)) {
        /* A light that used to be an area source and is now a point one must
         * lose its geometry, or a phantom emitter stays in the scene. */
        if (pi >= 0) ls_scene_remove_prim(d, pi);
        return;
    }

    /* The emissive material carries the light's own radiance, so the surface
     * seen by the camera and the source sampled by NEE agree by construction. */
    Material em;
    memset(&em, 0, sizeof em);
    em.bsdf.kind = LS_BSDF_LAMBERT;
    em.bsdf.rho  = ls_spectrum_zero();
    em.le        = ls_spectrum_scale(l->s_hat, l->radiance);
    em.emissive  = true;

    int mid;
    if (pi >= 0) {
        mid = d->prims[pi].mat_id;
        d->mats[mid] = em;                 /* refresh in place */
    } else {
        char nm[32];
        snprintf(nm, sizeof nm, "__emit%d", index);
        mid = ls_scene_add_material(d, em, nm);
        if (mid < 0) return;
    }

    Prim pr;
    memset(&pr, 0, sizeof pr);
    pr.mat_id = mid;
    pr.light_id = index;
    switch (l->kind) {
        case LS_LIGHT_RECT:
            pr.kind = LS_PRIM_QUAD; pr.c = l->p; pr.n = l->n;
            pr.ex = l->ex; pr.ey = l->ey; break;
        case LS_LIGHT_DISK:
            pr.kind = LS_PRIM_DISK; pr.c = l->p; pr.n = l->n; pr.r = l->radius; break;
        case LS_LIGHT_SPHERE:
            pr.kind = LS_PRIM_SPHERE; pr.c = l->p; pr.r = l->radius; break;
        default: return;
    }

    if (pi >= 0) d->prims[pi] = pr;
    else         ls_scene_add_prim(d, pr);
    ls_scene_rebuild(d);
}

int ls_scene_add_light(SceneDesc *d, Light l) {
    if (!grow_lights(d)) return -1;
    int id = d->nlights;
    ls_light_finalize(&l, id);
    d->lights[id] = l;
    d->nlights++;
    ls_scene_rebuild(d);
    ls_scene_sync_light_geom(d, id);
    return id;
}

void ls_scene_update_light(SceneDesc *d, int index) {
    if (index < 0 || index >= d->nlights) return;
    ls_light_finalize(&d->lights[index], index);
    ls_scene_sync_light_geom(d, index);
    ls_scene_rebuild(d);
}

void ls_scene_remove_light(SceneDesc *d, int index) {
    if (index < 0 || index >= d->nlights) return;

    /* Drop the paired geometry first, while `index` still identifies it. */
    for (int i = d->nprims - 1; i >= 0; --i)
        if (d->prims[i].light_id == index) ls_scene_remove_prim(d, i);

    for (int i = index; i + 1 < d->nlights; ++i) d->lights[i] = d->lights[i + 1];
    d->nlights--;

    /* Compacting shifts every later light down one, so both the light's own
     * index and every reference to it have to follow. Missing either leaves a
     * prim bound to the wrong emitter. */
    for (int i = index; i < d->nlights; ++i) d->lights[i].index = i;
    for (int i = 0; i < d->nprims; ++i)
        if (d->prims[i].light_id > index) d->prims[i].light_id--;

    ls_scene_rebuild(d);
}

bool ls_scene_rotate_light(SceneDesc *d, int index, vec3 axis, ls_real angle) {
    if (index < 0 || index >= d->nlights) return false;
    Light *l = &d->lights[index];
    vec3 a = v3norm(axis);
    switch (l->kind) {
        case LS_LIGHT_RECT:
            /* Turn the edge vectors and let the normal follow from them, so the
             * emitting face cannot end up disagreeing with the geometry. */
            l->ex = v3rotate(l->ex, a, angle);
            l->ey = v3rotate(l->ey, a, angle);
            l->n = v3norm(v3cross(l->ex, l->ey));
            break;
        case LS_LIGHT_DISK:
        case LS_LIGHT_SPOT:
        case LS_LIGHT_DIRECTIONAL:
            l->n = v3norm(v3rotate(l->n, a, angle));
            break;
        case LS_LIGHT_POINT:
        case LS_LIGHT_SPHERE:
            return false;                 /* isotropic: nothing to turn */
    }
    ls_scene_update_light(d, index);
    return true;
}

bool ls_scene_rotate_prim(SceneDesc *d, int index, vec3 axis, ls_real angle) {
    if (index < 0 || index >= d->nprims) return false;
    Prim *p = &d->prims[index];
    vec3 a = v3norm(axis);
    switch (p->kind) {
        /* A mesh turns by exactly the same three vectors: for it they are an
         * orthonormal object-to-world basis rather than half-edges plus a
         * normal, and a rotation preserves orthonormality either way. */
        case LS_PRIM_MESH:
        case LS_PRIM_QUAD:
            p->ex = v3rotate(p->ex, a, angle);
            p->ey = v3rotate(p->ey, a, angle);
            p->n = v3norm(v3rotate(p->n, a, angle));
            break;
        case LS_PRIM_PLANE:
        case LS_PRIM_DISK:
            p->n = v3norm(v3rotate(p->n, a, angle));
            break;
        case LS_PRIM_SPHERE:
            return false;
    }
    /* A prim bound to an area light is that light's own face; turn the light so
     * the two cannot drift apart. */
    if (p->light_id >= 0) return ls_scene_rotate_light(d, p->light_id, a, angle);
    ls_scene_rebuild(d);
    return true;
}

int ls_scene_duplicate_light(SceneDesc *d, int index, vec3 offset) {
    if (index < 0 || index >= d->nlights) return -1;
    Light copy = d->lights[index];
    copy.p = v3add(copy.p, offset);
    return ls_scene_add_light(d, copy);
}

/* ---------------------------------------------------------------- clone --- */

bool ls_scene_clone(const SceneDesc *src, SceneDesc *dst) {
    memset(dst, 0, sizeof *dst);

    /* Everything that is not a pointer copies wholesale; the arrays are then
     * replaced with private copies below. */
    *dst = *src;
    dst->prims = NULL; dst->mats = NULL; dst->names = NULL; dst->lights = NULL;
    dst->cap_prims = dst->nprims ? dst->nprims : 0;
    dst->cap_mats  = dst->nmats  ? dst->nmats  : 0;
    dst->cap_lights= dst->nlights? dst->nlights: 0;

    if (src->nprims > 0) {
        dst->prims = malloc((size_t)src->nprims * sizeof *dst->prims);
        if (!dst->prims) return false;
        memcpy(dst->prims, src->prims, (size_t)src->nprims * sizeof *dst->prims);
    }
    if (src->nmats > 0) {
        dst->mats = malloc((size_t)src->nmats * sizeof *dst->mats);
        dst->names = malloc((size_t)src->nmats * sizeof *dst->names);
        if (!dst->mats || !dst->names) return false;
        memcpy(dst->mats, src->mats, (size_t)src->nmats * sizeof *dst->mats);
        memcpy(dst->names, src->names, (size_t)src->nmats * sizeof *dst->names);
    }
    if (src->nlights > 0) {
        dst->lights = malloc((size_t)src->nlights * sizeof *dst->lights);
        if (!dst->lights) return false;
        memcpy(dst->lights, src->lights, (size_t)src->nlights * sizeof *dst->lights);
    }
    ls_scene_rebuild(dst);
    return true;
}
