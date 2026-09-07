/* Triangle meshes and their BVH.
 *
 * The BVH is pure optimisation: it must not change a single answer. So the
 * gate here is not a tolerance but exact agreement -- the accelerated
 * traversal and a brute-force scan run the same triangle test, so any
 * difference at all is a traversal bug, and asserting `==` catches it on the
 * first ray rather than as a faint bias later. */
#include "test.h"
#include "tests.h"
#include "lightsim/mesh.h"
#include "lightsim/rng.h"
#include <stdlib.h>
#include <string.h>

/* An axis-aligned box of 12 triangles, wound outward. */
static Mesh *build_box(double h) {
    vec3 v[8] = {
        v3(-h,-h,-h), v3( h,-h,-h), v3( h, h,-h), v3(-h, h,-h),
        v3(-h,-h, h), v3( h,-h, h), v3( h, h, h), v3(-h, h, h)
    };
    static const int idx[36] = {
        0,3,2, 0,2,1,   4,5,6, 4,6,7,      /* -z, +z */
        0,1,5, 0,5,4,   3,7,6, 3,6,2,      /* -y, +y */
        0,4,7, 0,7,3,   1,2,6, 1,6,5       /* -x, +x */
    };
    return ls_mesh_build(v, 8, idx, 12);
}

/* A UV sphere, deliberately big enough that the BVH is more than one leaf. */
static Mesh *build_sphere(double r, int nu, int nv) {
    int nverts = (nu + 1) * (nv + 1);
    vec3 *v = malloc((size_t)nverts * sizeof *v);
    int *idx = malloc((size_t)(nu * nv * 6) * sizeof *idx);
    int n = 0;
    for (int j = 0; j <= nv; ++j)
        for (int i = 0; i <= nu; ++i) {
            double phi = LS_TWO_PI * i / nu, th = LS_PI * j / nv;
            v[j * (nu + 1) + i] = v3(r * sin(th) * cos(phi),
                                     r * sin(th) * sin(phi),
                                     r * cos(th));
        }
    for (int j = 0; j < nv; ++j)
        for (int i = 0; i < nu; ++i) {
            int a = j * (nu + 1) + i, b = a + 1, c = a + nu + 1, d = c + 1;
            idx[n++] = a; idx[n++] = c; idx[n++] = b;
            idx[n++] = b; idx[n++] = c; idx[n++] = d;
        }
    Mesh *m = ls_mesh_build(v, nverts, idx, n / 3);
    free(v); free(idx);
    return m;
}

void test_mesh(void) {
    SECTION("a mesh knows its own size");
    {
        Mesh *box = build_box(0.5);
        CHECK(box != NULL);
        if (!box) return;
        CHECK(box->ntris == 12);
        /* Six faces of a unit cube. Exact: the areas are sums of halves. */
        CHECK_NEAR(box->area, 6.0, 1e-12);
        CHECK_NEAR(box->lo.x, -0.5, 1e-12);
        CHECK_NEAR(box->hi.z,  0.5, 1e-12);
        ls_mesh_release(box);
    }

    SECTION("malformed input is refused, not trusted");
    {
        vec3 v[3] = { v3(0,0,0), v3(1,0,0), v3(0,1,0) };
        int  ok[3] = { 0, 1, 2 };
        int  oob[3] = { 0, 1, 7 };            /* index past the vertex array */
        int  neg[3] = { 0, -1, 2 };

        Mesh *good = ls_mesh_build(v, 3, ok, 1);
        CHECK(good != NULL);
        if (good) { CHECK(good->ntris == 1); ls_mesh_release(good); }

        /* An out-of-range index is the commonest way a malformed OBJ crashes a
         * reader. It must be a NULL, not a traversal into someone else's heap. */
        CHECK(ls_mesh_build(v, 3, oob, 1) == NULL);
        CHECK(ls_mesh_build(v, 3, neg, 1) == NULL);
        CHECK(ls_mesh_build(v, 3, ok, 0) == NULL);
        CHECK(ls_mesh_build(NULL, 3, ok, 1) == NULL);

        /* A degenerate triangle is dropped, not kept with a garbage normal. */
        vec3 flat[3] = { v3(0,0,0), v3(1,0,0), v3(2,0,0) };
        CHECK(ls_mesh_build(flat, 3, ok, 1) == NULL);   /* nothing left */
    }

    SECTION("the BVH changes nothing");
    {
        Mesh *m = build_sphere(1.0, 40, 20);
        CHECK(m != NULL);
        if (!m) return;
        /* 40 of the 1600 are dropped: at the north pole sin(0) is exactly 0 so
         * the fan's triangles are exactly degenerate, while at the south pole
         * sin(pi) is 1.22e-16 and the same triangles survive as slivers. That
         * asymmetry is a fact about doubles, not a bug -- and it is why the
         * drop test is `len <= 0` rather than a tolerance: a sliver is still a
         * real triangle with a well-defined normal, and discarding it would
         * punch a hole in a closed surface. */
        CHECK(m->ntris < 1600 && m->ntris >= 1520);
        /* A tessellated sphere approaches 4 pi r^2 from below. */
        NOTE("sphere mesh: %d of 1600 triangles kept, %d BVH nodes, "
             "area %.4f vs 4 pi = %.4f",
             m->ntris, m->nnodes, m->area, 4.0 * LS_PI);
        CHECK(m->area < 4.0 * LS_PI);
        CHECK(m->area > 3.9 * LS_PI);
        CHECK(m->nnodes > 1);            /* it really did subdivide */

        Rng rng = ls_rng_seed(0x9E3779B97F4A7C15ull, 42);
        int hits = 0, misses = 0;
        for (int k = 0; k < 20000; ++k) {
            /* Origins outside and inside, directions all over, so that both
             * the empty-space rejections and the deep descents are exercised. */
            Ray r;
            r.o = v3(4.0 * ls_rng_f(&rng) - 2.0,
                     4.0 * ls_rng_f(&rng) - 2.0,
                     4.0 * ls_rng_f(&rng) - 2.0);
            vec3 d = v3(2.0 * ls_rng_f(&rng) - 1.0,
                        2.0 * ls_rng_f(&rng) - 1.0,
                        2.0 * ls_rng_f(&rng) - 1.0);
            if (v3len2(d) < 1e-12) continue;
            r.d = v3norm(d);
            r.tmin = 0.0;
            r.tmax = HUGE_VAL;

            ls_real tb = 0, tf = 0;
            vec3 nb = v3(0,0,0), nf = v3(0,0,0);
            bool hb = ls_mesh_intersect_local(m, &r, &tb, &nb, NULL);
            bool hf = ls_mesh_intersect_brute(m, &r, &tf, &nf, NULL);

            CHECK(hb == hf);
            if (hb && hf) {
                CHECK(tb == tf);                 /* bit-identical, not near */
                CHECK(nb.x == nf.x && nb.y == nf.y && nb.z == nf.z);
                hits++;
            } else misses++;

            /* Any-hit must agree with nearest-hit on whether anything is there.
             * They are separate traversals, so this is a real check. */
            CHECK(ls_mesh_occludes_local(m, &r) == hf);
        }
        NOTE("%d rays hit, %d missed; BVH and brute force agree exactly", hits, misses);
        CHECK(hits > 1000);                      /* the test actually hit things */
        CHECK(misses > 1000);                    /* and actually missed things */
        ls_mesh_release(m);
    }

    SECTION("a placed mesh is the same mesh, moved");
    {
        /* The Prim carries an orthonormal object->world basis, so a world ray
         * must hit at exactly the t its object-space equivalent does -- no
         * scaling means t means the same number in both frames, which is what
         * lets it be compared against every other primitive's t. */
        Mesh *box = build_box(0.5);
        CHECK(box != NULL);
        if (!box) return;

        Prim p;
        memset(&p, 0, sizeof p);
        p.kind = LS_PRIM_MESH;
        p.mat_id = 0;
        p.light_id = -1;
        /* A quarter turn about z, then translated. */
        p.ex = v3(0, 1, 0);
        p.ey = v3(-1, 0, 0);
        p.n  = v3(0, 0, 1);
        p.c  = v3(2.0, 3.0, 0.0);

        Rng rng = ls_rng_seed(0x2545F4914F6CDD1Dull, 7);
        int hits = 0;
        for (int k = 0; k < 4000; ++k) {
            vec3 o = v3(2.0 + 3.0 * ls_rng_f(&rng) - 1.5,
                        3.0 + 3.0 * ls_rng_f(&rng) - 1.5,
                        3.0 * ls_rng_f(&rng) - 1.5);
            vec3 target = v3(2.0 + ls_rng_f(&rng) - 0.5,
                             3.0 + ls_rng_f(&rng) - 0.5,
                             ls_rng_f(&rng) - 0.5);
            vec3 d = v3sub(target, o);
            if (v3len2(d) < 1e-12) continue;
            Ray r = { o, v3norm(d), 0.0, HUGE_VAL };

            Hit h;
            if (!ls_mesh_intersect(box, &p, &r, 3, &h)) continue;
            hits++;

            CHECK(h.prim_id == 3);
            CHECK(h.mat_id == 0);
            CHECK(h.light_id == -1);
            /* The hit point is on the ray at t, and the normal is a unit
             * vector consistent with the recorded orientation. */
            vec3 want = v3add(r.o, v3scale(r.d, h.t));
            CHECK_NEAR(v3dist(h.p, want), 0.0, 1e-12);
            CHECK_NEAR(v3len(h.ng), 1.0, 1e-12);
            CHECK(h.backface == (v3dot(r.d, h.ng) > 0.0));
            /* The box is a metre wide at (2,3,0); every hit is on its surface. */
            CHECK(fabs(h.p.x - 2.0) <= 0.5 + 1e-9);
            CHECK(fabs(h.p.y - 3.0) <= 0.5 + 1e-9);
            CHECK(fabs(h.p.z)       <= 0.5 + 1e-9);

            /* And occlusion agrees with intersection through the same placement. */
            CHECK(ls_mesh_occludes(box, &p, &r));
        }
        NOTE("placed box: %d of 4000 rays hit", hits);
        CHECK(hits > 500);
        ls_mesh_release(box);
    }

    SECTION("references keep geometry alive");
    {
        /* Undo snapshots share one payload; the last release frees it. Run
         * under ASan (make test-asan) this is the test that would catch a
         * double free or a use-after-free in the snapshot path. */
        Mesh *m = build_box(0.25);
        CHECK(m != NULL);
        if (!m) return;
        CHECK(m->refs == 1);
        Mesh *snap = ls_mesh_retain(m);
        CHECK(snap == m && m->refs == 2);

        ls_mesh_release(m);                       /* the live scene drops it */
        CHECK(snap->refs == 1);

        Ray r = { v3(0, 0, -2), v3(0, 0, 1), 0.0, HUGE_VAL };
        ls_real t;
        CHECK(ls_mesh_intersect_local(snap, &r, &t, NULL, NULL));
        CHECK_NEAR(t, 1.75, 1e-12);               /* still usable from the snapshot */
        ls_mesh_release(snap);
    }
}
