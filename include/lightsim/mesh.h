/* mesh.h — triangle meshes and the BVH that makes them affordable.
 *
 * WHY A MESH IS NOT A Prim
 *   Prim is a fixed-size POD copied by value into flat arrays and duplicated
 *   wholesale by ls_scene_clone for the editor's undo stack. There is nowhere
 *   in it to put a vertex buffer, and a raw pointer inside it would be shared
 *   by every snapshot and freed more than once. So a mesh lives out of line: the
 *   scene owns an array of Mesh pointers, and a Prim of kind LS_PRIM_MESH names
 *   one by index. The Prim still carries the material, the selection identity
 *   and the placement, so editing, picking and undo work unchanged.
 *
 * WHY THE GEOMETRY IS IMMUTABLE
 *   Vertices never change after ls_mesh_build. Placement lives on the Prim
 *   instead -- `c` translates and (ex, ey, n) is an ORTHONORMAL object-to-world
 *   basis -- so rays are transformed into object space rather than the mesh
 *   being transformed into world space. That keeps undo snapshots cheap (a
 *   100k-triangle mesh is refcounted, not copied 32 times) and lets the editor's
 *   existing rotation code turn a mesh by rotating those three vectors, exactly
 *   as it already turns a quad. Non-uniform scale is baked in at import, so the
 *   runtime transform stays rigid and `t` is directly comparable with the `t` of
 *   every other primitive.
 *
 * ACCELERATION
 *   Each mesh owns a BVH over its OWN triangles. The scene-level loop stays a
 *   linear scan over the (small) prim array, so this does not pre-empt the
 *   scene-wide BVH that scene.h reserves for M4 -- it just stops one imported
 *   model from costing 50k prim tests per ray.
 */
#ifndef LIGHTSIM_MESH_H
#define LIGHTSIM_MESH_H

#include "geom.h"

/* `count == 0` marks an interior node, whose children are at `first` and
 * `first + 1`. A leaf owns order[first .. first+count).
 *
 * Bounds are ls_real, not float. Float bounds are the usual optimisation but
 * they need outward rounding on every conversion, and a missed rounding is a
 * rare, ray-dependent MISS -- exactly the class of bug this codebase's
 * -ffp-contract=off and 1e-12 test bar exist to exclude. Take the memory. */
typedef struct {
    vec3 lo, hi;
    int  first, count, axis;
} BvhNode;

typedef struct { int v[3]; } Tri;

typedef struct {
    int      refs;                 /* shared by undo snapshots */
    vec3    *verts;   int nverts;  /* object space */
    Tri     *tris;    int ntris;
    vec3    *ng;                   /* per-triangle unit normal, raw winding */
    int     *order;                /* triangle permutation the leaves index */
    BvhNode *nodes;   int nnodes;
    vec3     lo, hi;               /* object-space AABB (== nodes[0]) */
    ls_real  area;                 /* sum of triangle areas, cached */
    /* Where it came from, so the scene writer can name the file rather than
     * serialise the vertices. Set by the importer; empty for a mesh built in
     * memory, which the writer then has to skip. */
    char     src_path[512];
    char     group[64];
    ls_real  scale;                /* what src_path's units were multiplied by */
} Mesh;

/* Copy `verts` and `idx` (3 indices per triangle) into a new mesh and build its
 * BVH. Returns NULL on allocation failure, empty input, or an out-of-range
 * index. The result has one reference.
 *
 * Degenerate triangles are dropped rather than kept with a garbage normal: they
 * cannot be hit and would only pollute the bounds. */
Mesh *ls_mesh_build(const vec3 *verts, int nv, const int *idx, int ntri);

/* Brute-force nearest hit, ignoring the BVH. Exists so the accelerated path has
 * something to be checked against; the test asserts bit-identical results. */
bool ls_mesh_intersect_brute(const Mesh *m, const Ray *r,
                             ls_real *t_out, vec3 *ng_out, int *tri_out);

Mesh *ls_mesh_retain(Mesh *m);
void  ls_mesh_release(Mesh *m);

/* Object-space queries. These are the tested core; the Prim-aware wrappers
 * below are a transform around them. `ng` is the raw geometric normal in
 * winding order and is never flipped, matching the Hit contract in geom.h. */
bool ls_mesh_intersect_local(const Mesh *m, const Ray *r,
                             ls_real *t_out, vec3 *ng_out, int *tri_out);
bool ls_mesh_occludes_local(const Mesh *m, const Ray *r);

/* World-space queries for a mesh placed by `p` (kind LS_PRIM_MESH, with an
 * orthonormal ex/ey/n). Fills `hit` exactly as ls_prim_intersect does. */
bool ls_mesh_intersect(const Mesh *m, const Prim *p, const Ray *r,
                       int prim_id, Hit *hit);
bool ls_mesh_occludes(const Mesh *m, const Prim *p, const Ray *r);

#endif /* LIGHTSIM_MESH_H */
