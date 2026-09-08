#include "lightsim/mesh.h"
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------- triangle -- */

/* Moller-Trumbore, two-sided.
 *
 * NOT backface-culled: the transport core depends on hitting both sides. The
 * integrator builds its shading frame as `backface ? -ng : ng`, and the
 * emitter one-sidedness gate reads the raw ng -- culling here would make a
 * surface invisible from behind rather than merely dark, which is a different
 * and much worse thing.
 *
 * The inside test is INCLUSIVE on all three edges. Moller-Trumbore is not
 * watertight, and with strict comparisons a ray meeting the diagonal shared by
 * two triangles can miss BOTH. The measurement grid is a regular lattice, and
 * an axis-aligned imported floor split on an axis-aligned diagonal is exactly
 * the arrangement that lines those up -- the symptom would be scattered dead
 * cells in the field map that read as noise. Inclusive comparisons hit the
 * shared edge twice instead, and since both triangles return the same t the
 * nearest-hit answer is unchanged. */
static bool tri_hit(const Mesh *m, int t, const Ray *r, ls_real *t_out) {
    const Tri *tr = &m->tris[t];
    vec3 a = m->verts[tr->v[0]];
    vec3 e1 = v3sub(m->verts[tr->v[1]], a);
    vec3 e2 = v3sub(m->verts[tr->v[2]], a);

    vec3 pv = v3cross(r->d, e2);
    ls_real det = v3dot(e1, pv);
    if (fabs(det) < 1e-15) return false;          /* ray parallel to the plane */
    ls_real inv = 1.0 / det;

    vec3 tv = v3sub(r->o, a);
    ls_real u = v3dot(tv, pv) * inv;
    if (u < 0.0 || u > 1.0) return false;

    vec3 qv = v3cross(tv, e1);
    ls_real v = v3dot(r->d, qv) * inv;
    if (v < 0.0 || u + v > 1.0) return false;

    ls_real tt = v3dot(e2, qv) * inv;
    if (tt <= r->tmin || tt >= r->tmax) return false;
    *t_out = tt;
    return true;
}

bool ls_mesh_intersect_brute(const Mesh *m, const Ray *ray,
                             ls_real *t_out, vec3 *ng_out, int *tri_out) {
    Ray r = *ray;
    bool found = false;
    for (int i = 0; i < m->ntris; ++i) {
        ls_real t;
        if (tri_hit(m, i, &r, &t)) {
            r.tmax = t;
            *t_out = t;
            if (ng_out)  *ng_out = m->ng[i];
            if (tri_out) *tri_out = i;
            found = true;
        }
    }
    return found;
}

/* ------------------------------------------------------------------ bvh -- */

#define BVH_LEAF_MAX 4
#define BVH_MAX_DEPTH 63          /* what licenses the fixed traversal stack */

static void box_reset(vec3 *lo, vec3 *hi) {
    *lo = v3( HUGE_VAL,  HUGE_VAL,  HUGE_VAL);
    *hi = v3(-HUGE_VAL, -HUGE_VAL, -HUGE_VAL);
}
static void box_add(vec3 *lo, vec3 *hi, vec3 p) {
    if (p.x < lo->x) lo->x = p.x;  if (p.x > hi->x) hi->x = p.x;
    if (p.y < lo->y) lo->y = p.y;  if (p.y > hi->y) hi->y = p.y;
    if (p.z < lo->z) lo->z = p.z;  if (p.z > hi->z) hi->z = p.z;
}
static ls_real axis_of(vec3 v, int a) { return a == 0 ? v.x : (a == 1 ? v.y : v.z); }

typedef struct {
    Mesh *m;
    vec3 *cent;                   /* per-triangle centroid, build scratch */
} Build;

/* Bounds of the triangles in order[first .. first+count). */
static void node_bounds(const Build *b, int first, int count, vec3 *lo, vec3 *hi) {
    box_reset(lo, hi);
    for (int i = first; i < first + count; ++i) {
        const Tri *tr = &b->m->tris[b->m->order[i]];
        for (int k = 0; k < 3; ++k) box_add(lo, hi, b->m->verts[tr->v[k]]);
    }
}

/* Median split on the widest axis of the centroid bounds. Simple and adequate;
 * a binned SAH slots in behind this same signature without touching traversal
 * or the agreement test that guards it.
 *
 * Fills the node at `self`, which the CALLER has already reserved. The pair of
 * children is reserved together, before either subtree is built, so that they
 * really are adjacent -- allocating the right child after the left subtree had
 * finished would put an entire subtree between them, and `first + 1` would
 * address a grandchild. */
static void build_into(Build *b, int self, int first, int count, int depth) {
    Mesh *m = b->m;
    BvhNode *nd = &m->nodes[self];
    node_bounds(b, first, count, &nd->lo, &nd->hi);
    nd->first = first;
    nd->count = count;
    nd->axis = 0;

    if (count <= BVH_LEAF_MAX || depth >= BVH_MAX_DEPTH) return;

    vec3 clo, chi;
    box_reset(&clo, &chi);
    for (int i = first; i < first + count; ++i) box_add(&clo, &chi, b->cent[m->order[i]]);
    vec3 ext = v3sub(chi, clo);
    int axis = (ext.x > ext.y)
             ? ((ext.x > ext.z) ? 0 : 2)
             : ((ext.y > ext.z) ? 1 : 2);
    if (axis_of(ext, axis) <= 0.0) return;         /* all centroids coincide */

    /* Partition about the centroid median. Selection by repeated partition
     * would be faster; at these sizes the build is not the bottleneck. */
    int mid = first + count / 2;
    for (int i = first; i < first + count; ++i)
        for (int j = i + 1; j < first + count; ++j)
            if (axis_of(b->cent[m->order[j]], axis) <
                axis_of(b->cent[m->order[i]], axis)) {
                int tmp = m->order[i]; m->order[i] = m->order[j]; m->order[j] = tmp;
            }

    int left = m->nnodes;
    m->nnodes += 2;
    build_into(b, left,     first, mid - first,          depth + 1);
    build_into(b, left + 1, mid,   first + count - mid,  depth + 1);

    nd = &m->nodes[self];
    nd->first = left;
    nd->count = 0;                                 /* 0 marks an interior node */
    nd->axis = axis;
}

/* Slab test. inv_d may be +/-inf for an axis-aligned ray; the min/max ordering
 * resolves that correctly, and the multiply by (1 + 2*eps) on the exit distance
 * keeps a ray originating exactly on a slab plane from being rejected by
 * rounding. */
static bool box_hit(const BvhNode *n, vec3 o, vec3 inv, ls_real tmin, ls_real tmax) {
    ls_real t0, t1, lo, hi;
    t0 = (n->lo.x - o.x) * inv.x;  t1 = (n->hi.x - o.x) * inv.x;
    lo = t0 < t1 ? t0 : t1;        hi = t0 < t1 ? t1 : t0;
    if (lo > tmin) tmin = lo;      if (hi < tmax) tmax = hi;
    t0 = (n->lo.y - o.y) * inv.y;  t1 = (n->hi.y - o.y) * inv.y;
    lo = t0 < t1 ? t0 : t1;        hi = t0 < t1 ? t1 : t0;
    if (lo > tmin) tmin = lo;      if (hi < tmax) tmax = hi;
    t0 = (n->lo.z - o.z) * inv.z;  t1 = (n->hi.z - o.z) * inv.z;
    lo = t0 < t1 ? t0 : t1;        hi = t0 < t1 ? t1 : t0;
    if (lo > tmin) tmin = lo;      if (hi < tmax) tmax = hi;
    return tmin <= tmax * 1.0000000000000004;      /* 1 + 2*DBL_EPSILON */
}

static vec3 inv_dir(vec3 d) {
    return v3(1.0 / d.x, 1.0 / d.y, 1.0 / d.z);    /* +/-inf is intended */
}

bool ls_mesh_intersect_local(const Mesh *m, const Ray *ray,
                             ls_real *t_out, vec3 *ng_out, int *tri_out) {
    if (!m || m->nnodes <= 0) return false;
    vec3 inv = inv_dir(ray->d);
    ls_real tmax = ray->tmax;
    bool found = false;

    int stack[BVH_MAX_DEPTH + 2], sp = 0;
    stack[sp++] = 0;
    while (sp > 0) {
        const BvhNode *n = &m->nodes[stack[--sp]];
        if (!box_hit(n, ray->o, inv, ray->tmin, tmax)) continue;
        if (n->count == 0) {
            /* Near child first, so the running tmax culls the far one on pop. */
            int a = n->first, bch = n->first + 1;
            if (axis_of(ray->d, n->axis) < 0.0) { int t = a; a = bch; bch = t; }
            stack[sp++] = bch;
            stack[sp++] = a;
            continue;
        }
        for (int i = n->first; i < n->first + n->count; ++i) {
            int tri = m->order[i];
            Ray q = *ray;
            q.tmax = tmax;
            ls_real t;
            if (tri_hit(m, tri, &q, &t)) {
                tmax = t;
                *t_out = t;
                if (ng_out)  *ng_out = m->ng[tri];
                if (tri_out) *tri_out = tri;
                found = true;
            }
        }
    }
    return found;
}

/* Any-hit. Differs from nearest-hit in ways that compound: it returns on the
 * first triangle in range, so there is no running tmax to shrink, so ordering
 * buys nothing and both children are pushed unconditionally, and no normal or
 * hit point is ever computed. Shadow rays are where this engine spends its
 * time -- integrator.c loops every light on every bounce -- so the specialised
 * form is worth having rather than aliasing to nearest-hit-and-discard. */
bool ls_mesh_occludes_local(const Mesh *m, const Ray *ray) {
    if (!m || m->nnodes <= 0) return false;
    vec3 inv = inv_dir(ray->d);

    int stack[BVH_MAX_DEPTH + 2], sp = 0;
    stack[sp++] = 0;
    while (sp > 0) {
        const BvhNode *n = &m->nodes[stack[--sp]];
        if (!box_hit(n, ray->o, inv, ray->tmin, ray->tmax)) continue;
        if (n->count == 0) {
            stack[sp++] = n->first;
            stack[sp++] = n->first + 1;
            continue;
        }
        for (int i = n->first; i < n->first + n->count; ++i) {
            ls_real t;
            if (tri_hit(m, m->order[i], ray, &t)) return true;
        }
    }
    return false;
}

/* ---------------------------------------------------------------- build -- */

Mesh *ls_mesh_build(const vec3 *verts, int nv, const int *idx, int ntri) {
    if (!verts || !idx || nv <= 0 || ntri <= 0) return NULL;

    Mesh *m = calloc(1, sizeof *m);
    if (!m) return NULL;
    m->refs = 1;

    m->verts = malloc((size_t)nv * sizeof *m->verts);
    m->tris  = malloc((size_t)ntri * sizeof *m->tris);
    m->ng    = malloc((size_t)ntri * sizeof *m->ng);
    if (!m->verts || !m->tris || !m->ng) { ls_mesh_release(m); return NULL; }
    memcpy(m->verts, verts, (size_t)nv * sizeof *m->verts);
    m->nverts = nv;

    /* Copy triangles, dropping degenerates and refusing bad indices. A stray
     * index is the commonest way a malformed OBJ crashes a reader, so it is
     * rejected here rather than trusted into the traversal. */
    m->ntris = 0;
    m->area = 0.0;
    for (int i = 0; i < ntri; ++i) {
        int a = idx[i*3+0], b = idx[i*3+1], c = idx[i*3+2];
        if (a < 0 || a >= nv || b < 0 || b >= nv || c < 0 || c >= nv) {
            ls_mesh_release(m);
            return NULL;
        }
        Tri t = { { a, b, c } };
        m->tris[m->ntris] = t;
        vec3 n = v3cross(v3sub(verts[b], verts[a]), v3sub(verts[c], verts[a]));
        ls_real len = v3len(n);
        if (len <= 0.0) continue;                  /* degenerate; drop it */
        m->ng[m->ntris] = v3scale(n, 1.0 / len);
        m->area += 0.5 * len;
        m->ntris++;
    }
    if (m->ntris == 0) { ls_mesh_release(m); return NULL; }

    m->order = malloc((size_t)m->ntris * sizeof *m->order);
    /* A median split makes at most 2n-1 nodes; +1 keeps the arithmetic honest
     * for a single-triangle mesh. */
    m->nodes = malloc((size_t)(2 * m->ntris + 1) * sizeof *m->nodes);
    vec3 *cent = malloc((size_t)m->ntris * sizeof *cent);
    if (!m->order || !m->nodes || !cent) {
        free(cent); ls_mesh_release(m); return NULL;
    }
    for (int i = 0; i < m->ntris; ++i) {
        m->order[i] = i;
        const Tri *tr = &m->tris[i];
        cent[i] = v3scale(v3add(v3add(m->verts[tr->v[0]], m->verts[tr->v[1]]),
                                m->verts[tr->v[2]]), 1.0 / 3.0);
    }

    Build b = { m, cent };
    m->nnodes = 1;                                 /* the root occupies slot 0 */
    build_into(&b, 0, 0, m->ntris, 0);
    free(cent);

    m->lo = m->nodes[0].lo;
    m->hi = m->nodes[0].hi;
    return m;
}

Mesh *ls_mesh_retain(Mesh *m) {
    if (m) m->refs++;
    return m;
}

void ls_mesh_release(Mesh *m) {
    if (!m || --m->refs > 0) return;
    free(m->verts); free(m->tris); free(m->ng);
    free(m->order); free(m->nodes);
    free(m);
}

/* ------------------------------------------------------- placed in world -- */

/* World -> object. The basis is orthonormal, so its inverse is its transpose,
 * and `r` carries a uniform scale.
 *
 * The object-space direction is deliberately NOT renormalised. Dividing both
 * the origin offset and the direction by the scale leaves the ray parameter
 * meaning the same number in both spaces:
 *
 *     Q(t) = R^T(o + t*d - c)/s  =  [R^T(o-c)/s] + t*[R^T d / s]
 *
 * so `t` comes back directly comparable with every other primitive's, and
 * ray.tmin/tmax pass through untouched. Renormalising d' would scale t by s and
 * quietly break both -- including ls_scene_occluded's `dist * (1 - 1e-6)`,
 * which is a length in world units. */
static Ray to_object(const Prim *p, const Ray *r) {
    ls_real s = (p->r > 0.0) ? p->r : 1.0;
    ls_real inv = 1.0 / s;
    vec3 rel = v3sub(r->o, p->c);
    Ray q;
    q.o = v3scale(v3(v3dot(rel, p->ex), v3dot(rel, p->ey), v3dot(rel, p->n)), inv);
    q.d = v3scale(v3(v3dot(r->d, p->ex), v3dot(r->d, p->ey), v3dot(r->d, p->n)), inv);
    q.tmin = r->tmin;
    q.tmax = r->tmax;
    return q;
}

/* World-space bounds of a placed mesh, for the editor's outline and for
 * deciding where a dropped object sits. */
void ls_mesh_world_bounds(const Mesh *m, const Prim *p, vec3 *lo, vec3 *hi) {
    ls_real s = (p->r > 0.0) ? p->r : 1.0;
    box_reset(lo, hi);
    for (int k = 0; k < 8; ++k) {
        vec3 o = v3((k & 1) ? m->hi.x : m->lo.x,
                    (k & 2) ? m->hi.y : m->lo.y,
                    (k & 4) ? m->hi.z : m->lo.z);
        o = v3scale(o, s);
        box_add(lo, hi, v3add(p->c, v3add(v3scale(p->ex, o.x),
                              v3add(v3scale(p->ey, o.y), v3scale(p->n, o.z)))));
    }
}

bool ls_mesh_intersect(const Mesh *m, const Prim *p, const Ray *r,
                       int prim_id, Hit *hit) {
    Ray q = to_object(p, r);
    ls_real t;
    vec3 ng;
    if (!ls_mesh_intersect_local(m, &q, &t, &ng, NULL)) return false;

    /* Object -> world for the normal is the same orthonormal basis applied
     * forward. A UNIFORM scale leaves the inverse-transpose parallel to R, so
     * it changes the normal's length and not its direction, and ng is already
     * unit. A non-uniform scale would not, which is why only one is offered. */
    vec3 nw = v3add(v3scale(p->ex, ng.x),
                    v3add(v3scale(p->ey, ng.y), v3scale(p->n, ng.z)));
    hit->t        = t;
    hit->p        = v3add(r->o, v3scale(r->d, t));
    hit->ng       = nw;
    hit->prim_id  = prim_id;
    hit->mat_id   = p->mat_id;
    hit->light_id = p->light_id;
    hit->backface = v3dot(r->d, nw) > 0.0;
    return true;
}

bool ls_mesh_occludes(const Mesh *m, const Prim *p, const Ray *r) {
    Ray q = to_object(p, r);
    return ls_mesh_occludes_local(m, &q);
}
