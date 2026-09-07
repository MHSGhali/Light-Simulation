/* E1 validation: a scene can be mutated, saved and reloaded without drift.
 *
 * The invariant these guard is the one that is invisible when it breaks: an
 * area light and its paired emissive geometry must stay in step, or the render
 * shows a source the field map does not sample. Neither picture looks wrong. */
#include "test.h"
#include "tests.h"
#include "lightsim/sceneedit.h"
#include "lightsim/integrator.h"
#include <string.h>
#include "lightsim/units.h"

/* Set a light's luminous flux the way the inspector does: the authored value is
 * what the user asked for, and the watts follow from the spectrum. */
static bool ls_inspect_set_shim(SceneDesc *d, int i, double lumens) {
    if (i < 0 || i >= d->nlights) return false;
    d->lights[i].flux_in_lumens = true;
    d->lights[i].flux_authored = lumens;
    d->lights[i].phi_e = ls_watts_from_lumens(lumens, &d->lights[i].s_hat);
    ls_scene_update_light(d, i);
    return true;
}

#define SCENE_IN  "scenes/workcell.scene"
#define SCENE_OUT "/tmp/lightsim_roundtrip.scene"

static double probe(const SceneDesc *d, vec3 p) {
    SpectrumAcc acc;
    Rng rng = ls_rng_seed(0xA5A5A5A5u, 11);
    ls_estimate_irradiance_full(&d->scene, p, v3(0, 0, 1), 64, 64, 4, &rng, &acc, NULL);
    Spectrum E = ls_acc_mean(&acc, 1);
    return ls_spectrum_integrate(&E);
}

static void check_light_equal(const Light *x, const Light *y, const char *what) {
    LS_UNUSED(what);
    CHECK(x->kind == y->kind);
    CHECK_NEAR(x->phi_e, y->phi_e, 1e-9);
    CHECK_NEAR(x->area, y->area, 1e-9);
    CHECK_NEAR(x->radiance, y->radiance, 1e-9);
    CHECK_NEAR(x->p.x, y->p.x, 1e-9);
    CHECK_NEAR(x->p.y, y->p.y, 1e-9);
    CHECK_NEAR(x->p.z, y->p.z, 1e-9);
    CHECK(x->spd_kind == y->spd_kind);
    CHECK_NEAR(x->spd_a, y->spd_a, 1e-9);
}

/* Every area light has exactly one prim bound to it, and that prim's geometry
 * and emissive radiance match the light. */
static void check_pairing(const SceneDesc *d) {
    for (int i = 0; i < d->nlights; ++i) {
        const Light *l = &d->lights[i];
        CHECK(l->index == i);
        bool area = (l->kind == LS_LIGHT_RECT || l->kind == LS_LIGHT_DISK
                     || l->kind == LS_LIGHT_SPHERE);
        int count = 0, pi = -1;
        for (int k = 0; k < d->nprims; ++k)
            if (d->prims[k].light_id == i) { count++; pi = k; }
        CHECK(count == (area ? 1 : 0));
        if (!area || pi < 0) continue;

        const Prim *pr = &d->prims[pi];
        CHECK_NEAR(pr->c.x, l->p.x, 1e-9);
        CHECK_NEAR(pr->c.y, l->p.y, 1e-9);
        CHECK_NEAR(pr->c.z, l->p.z, 1e-9);
        if (l->kind == LS_LIGHT_RECT) {
            CHECK(pr->kind == LS_PRIM_QUAD);
            CHECK_NEAR(v3len(pr->ex), v3len(l->ex), 1e-9);
            CHECK_NEAR(v3len(pr->ey), v3len(l->ey), 1e-9);
        } else if (l->kind == LS_LIGHT_DISK) {
            CHECK(pr->kind == LS_PRIM_DISK);
            CHECK_NEAR(pr->r, l->radius, 1e-9);
        } else {
            CHECK(pr->kind == LS_PRIM_SPHERE);
            CHECK_NEAR(pr->r, l->radius, 1e-9);
        }
        /* The emissive material must carry this light's radiance. */
        const Material *m = &d->mats[pr->mat_id];
        CHECK(m->emissive);
        Spectrum want = ls_spectrum_scale(l->s_hat, l->radiance);
        CHECK_NEAR(ls_spectrum_integrate(&m->le), ls_spectrum_integrate(&want), 1e-6);
    }
}

void test_scene(void) {
    SECTION("rotation maths and gizmo ray queries");
    {
        /* Rodrigues: a quarter turn about +z takes +x to +y. */
        vec3 r = v3rotate(v3(1, 0, 0), v3(0, 0, 1), LS_PI / 2.0);
        CHECK_NEAR(r.x, 0.0, 1e-12);
        CHECK_NEAR(r.y, 1.0, 1e-12);
        CHECK_NEAR(r.z, 0.0, 1e-12);

        /* A vector on the axis is unmoved, and length is preserved. */
        vec3 onaxis = v3rotate(v3(0, 0, 3), v3(0, 0, 1), 1.234);
        CHECK_NEAR(onaxis.z, 3.0, 1e-12);
        Rng rng = ls_rng_seed(0xB0Bu, 3);
        for (int i = 0; i < 500; ++i) {
            vec3 v = v3(ls_rng_f(&rng) - 0.5, ls_rng_f(&rng) - 0.5, ls_rng_f(&rng) - 0.5);
            vec3 ax = v3norm(v3(ls_rng_f(&rng) - 0.5, ls_rng_f(&rng) - 0.5,
                                ls_rng_f(&rng) - 0.5));
            double ang = (ls_rng_f(&rng) - 0.5) * 6.0;
            vec3 w = v3rotate(v, ax, ang);
            CHECK_NEAR(v3len(w), v3len(v), 1e-12);
            /* Rotating back returns the original. */
            vec3 back = v3rotate(w, ax, -ang);
            CHECK_NEAR(v3dist(back, v), 0.0, 1e-9);
        }

        /* Closest point on an axis to a ray: a ray crossing the x axis at
         * x = 2 must give t = 2. */
        double t = 0.0;
        CHECK(ls_line_closest_t(v3(0, 0, 0), v3(1, 0, 0),
                                v3(2, -5, 0), v3(0, 1, 0), &t));
        CHECK_NEAR(t, 2.0, 1e-12);
        /* Parallel lines have no well-defined closest point. */
        CHECK(!ls_line_closest_t(v3(0, 0, 0), v3(1, 0, 0),
                                 v3(0, 1, 0), v3(1, 0, 0), &t));

        /* Ray-plane, and the parallel rejection. */
        vec3 hit;
        CHECK(ls_ray_plane(v3(0, 0, 5), v3(0, 0, -1), v3(0, 0, 1), v3(0, 0, 1), &hit));
        CHECK_NEAR(hit.z, 1.0, 1e-12);
        CHECK(!ls_ray_plane(v3(0, 0, 5), v3(1, 0, 0), v3(0, 0, 1), v3(0, 0, 1), &hit));
    }

    SECTION("rotating an object keeps it physically consistent");
    {
        SceneDesc d;
        if (!ls_scene_load(&d, SCENE_IN)) return;

        const Light *l0 = &d.lights[0];
        CHECK(l0->kind == LS_LIGHT_RECT);
        double area0 = l0->area, phi0 = l0->phi_e, rad0 = l0->radiance;
        vec3 n0 = l0->n;

        /* A quarter turn about +x should tip a downward-facing panel to face
         * along -y, and must change nothing about how much light it emits. */
        CHECK(ls_scene_rotate_light(&d, 0, v3(1, 0, 0), LS_PI / 2.0));
        CHECK_NEAR(d.lights[0].area, area0, 1e-12);
        CHECK_NEAR(d.lights[0].phi_e, phi0, 1e-12);
        CHECK_NEAR(d.lights[0].radiance, rad0, 1e-12);
        CHECK_NEAR(v3len(d.lights[0].n), 1.0, 1e-12);
        CHECK_NEAR(v3dot(d.lights[0].n, n0), 0.0, 1e-9);   /* turned 90 degrees */
        /* The flux self-check inside finalize still holds. */
        CHECK_NEAR(ls_light_emitted_flux(&d.lights[0]), phi0, 1e-9);
        /* And the paired geometry turned with it. */
        check_pairing(&d);

        /* Turning back restores the original orientation exactly. */
        CHECK(ls_scene_rotate_light(&d, 0, v3(1, 0, 0), -LS_PI / 2.0));
        CHECK_NEAR(v3dist(d.lights[0].n, n0), 0.0, 1e-9);

        /* An isotropic source has no orientation to turn. */
        Light sph = ls_light_sphere(v3(0, 0, 0.2), 0.02, 0.5, ls_spectrum_const(1.0));
        int si = ls_scene_add_light(&d, sph);
        CHECK(!ls_scene_rotate_light(&d, si, v3(0, 0, 1), 1.0));

        /* Rotating an area light's own face turns the LIGHT, so the two cannot
         * drift apart -- the failure that is invisible in the field map. */
        int pi = ls_scene_light_prim(&d, 0);
        CHECK(pi >= 0);
        vec3 before = d.lights[0].n;
        CHECK(ls_scene_rotate_prim(&d, pi, v3(0, 1, 0), 0.3));
        CHECK(v3dist(d.lights[0].n, before) > 1e-6);
        check_pairing(&d);

        ls_scene_desc_free(&d);
    }

    SECTION("camera projection is the exact inverse of the pick ray");
    {
        /* Picking a light that has no geometry, and drawing every overlay
         * gizmo, both go through ls_camera_project. If it disagrees with
         * ls_camera_pick_ray by even a pixel the markers drift off the objects
         * they label, and picking misses at the edges of the frame -- so this
         * is asserted at double precision rather than "close enough". */
        Rng rng = ls_rng_seed(0x4242u, 1);
        double worst = 0.0;
        int n = 0;
        for (int trial = 0; trial < 120; ++trial) {
            double az = ls_rng_f(&rng) * LS_TWO_PI;
            double el = (ls_rng_f(&rng) - 0.5) * 2.6;
            double dist = 0.5 + 3.0 * ls_rng_f(&rng);
            vec3 tgt = v3(ls_rng_f(&rng) - 0.5, ls_rng_f(&rng) - 0.5, ls_rng_f(&rng) - 0.5);
            vec3 eye = v3add(tgt, v3(dist * cos(el) * cos(az),
                                     dist * cos(el) * sin(az),
                                     dist * sin(el)));
            Camera c = ls_camera_look_at(eye, tgt, v3(0, 0, 1),
                                         20.0 + 60.0 * ls_rng_f(&rng), 640, 480);
            for (int k = 0; k < 20; ++k) {
                double px = ls_rng_f(&rng) * 640.0, py = ls_rng_f(&rng) * 480.0;
                Ray r = ls_camera_pick_ray(&c, px, py);
                vec3 p = v3add(r.o, v3scale(r.d, 0.05 + 5.0 * ls_rng_f(&rng)));
                double qx, qy;
                CHECK(ls_camera_project(&c, p, &qx, &qy));
                double e = ls_max(fabs(qx - px), fabs(qy - py));
                if (e > worst) worst = e;
                n++;
            }
        }
        CHECK(worst < 1e-9);
        NOTE("projection round-trip over %d samples: worst error %.2e px", n, worst);

        /* A point behind the eye has no pixel, and must say so rather than
         * returning a plausible coordinate from a negative depth. */
        Camera c = ls_camera_look_at(v3(0, 0, 0), v3(0, 1, 0), v3(0, 0, 1), 45, 100, 100);
        double qx, qy;
        CHECK(!ls_camera_project(&c, v3(0, -1, 0), &qx, &qy));
        CHECK(ls_camera_project(&c, v3(0, 1, 0), &qx, &qy));
        CHECK_NEAR(qx, 50.0, 1e-9);
        CHECK_NEAR(qy, 50.0, 1e-9);
    }

    SECTION("scene round-trip through the writer");
    {
        SceneDesc a, b;
        if (!ls_scene_load(&a, SCENE_IN)) {
            printf("  SKIP: cannot load %s (%s)\n", SCENE_IN, a.err);
            return;
        }
        CHECK(ls_scene_save(&a, SCENE_OUT));
        CHECK(ls_scene_load(&b, SCENE_OUT));

        CHECK(a.nlights == b.nlights);
        CHECK(a.nprims  == b.nprims);
        CHECK(a.nmats   == b.nmats);
        CHECK(a.has_grid == b.has_grid);
        CHECK(a.grid_nu == b.grid_nu && a.grid_nv == b.grid_nv);
        for (int i = 0; i < a.nlights && i < b.nlights; ++i)
            check_light_equal(&a.lights[i], &b.lights[i], "round-trip");

        /* The real test: the reloaded scene must simulate identically. Same
         * seed, same estimator -- any drift in geometry or flux shows here. */
        double ea = probe(&a, v3(0.1, 0.05, 0.001));
        double eb = probe(&b, v3(0.1, 0.05, 0.001));
        CHECK_NEAR(eb, ea, 1e-12);
        NOTE("round-trip irradiance %.9g vs %.9g W/m^2", ea, eb);

        check_pairing(&a);
        check_pairing(&b);
        ls_scene_desc_free(&a);
        ls_scene_desc_free(&b);
    }

    SECTION("light and geometry stay paired through edits");
    {
        SceneDesc d;
        if (!ls_scene_load(&d, SCENE_IN)) return;
        int n0 = d.nlights, p0 = d.nprims;
        check_pairing(&d);

        /* Add a rect light: gains one light, one prim, one material. */
        Light l = ls_light_rect(v3(0, 0, 0.3), v3(0.04, 0, 0), v3(0, -0.04, 0),
                                0.5, ls_spectrum_const(1.0));
        int id = ls_scene_add_light(&d, l);
        CHECK(id == n0);
        CHECK(d.nlights == n0 + 1);
        CHECK(d.nprims == p0 + 1);
        check_pairing(&d);

        /* Move it: the paired prim must follow. */
        d.lights[id].p = v3(0.1, -0.05, 0.25);
        ls_scene_update_light(&d, id);
        check_pairing(&d);

        /* Resize it: extent and radiance both change. */
        d.lights[id].ex = v3(0.08, 0, 0);
        ls_scene_update_light(&d, id);
        check_pairing(&d);
        CHECK_NEAR(d.lights[id].area, 4.0 * 0.08 * 0.04, 1e-12);

        /* Retype it to a point light: the geometry must go away entirely, or a
         * phantom emitter stays in the scene radiating on its own. */
        d.lights[id].kind = LS_LIGHT_POINT;
        ls_scene_update_light(&d, id);
        CHECK(d.nprims == p0);
        CHECK(ls_scene_light_prim(&d, id) == -1);
        check_pairing(&d);

        ls_scene_remove_light(&d, id);
        CHECK(d.nlights == n0);
        check_pairing(&d);
        ls_scene_desc_free(&d);
    }

    SECTION("removal compacts indices and remaps references");
    {
        SceneDesc d;
        if (!ls_scene_load(&d, SCENE_IN)) return;
        CHECK(d.nlights >= 3);
        vec3 kept = d.lights[2].p;          /* remember a light past the hole */
        int n0 = d.nlights, p0 = d.nprims;

        ls_scene_remove_light(&d, 1);
        CHECK(d.nlights == n0 - 1);
        CHECK(d.nprims == p0 - 1);
        /* What was light 2 is now light 1, still with its own position, and
         * every prim must point at the light it actually belongs to. */
        CHECK_NEAR(d.lights[1].p.x, kept.x, 1e-12);
        CHECK_NEAR(d.lights[1].p.y, kept.y, 1e-12);
        check_pairing(&d);
        for (int i = 0; i < d.nprims; ++i)
            CHECK(d.prims[i].light_id < d.nlights);
        ls_scene_desc_free(&d);
    }

    SECTION("build, edit, save, reload, re-solve");
    {
        /* The loop the editor exists to support, exercised end to end through
         * exactly the calls the UI gestures make. */
        SceneDesc d;
        if (!ls_scene_load(&d, SCENE_IN)) return;

        /* Place a light, as the ADD LIGHT tool does. */
        Spectrum spd = ls_spectrum_daylight(4000.0);
        Light l = ls_light_rect(v3(0.0, 0.0, 0.35), v3(0.05, 0, 0), v3(0, -0.05, 0),
                                ls_watts_from_lumens(400.0, &spd), spd);
        l.spd_kind = LS_SPD_DAYLIGHT;
        l.spd_a = 4000.0;
        l.flux_in_lumens = true;
        l.flux_authored = 400.0;
        int id = ls_scene_add_light(&d, l);
        CHECK(id >= 0);

        /* Retune it through the inspector, as scrubbing or typing does. */
        CHECK(ls_inspect_set_shim(&d, id, 850.0));
        check_pairing(&d);

        double before = probe(&d, v3(0.0, 0.0, 0.001));
        CHECK(before > 0.0);

        /* Save, reload, and the reloaded scene must simulate identically --
         * which is what makes a session recoverable rather than merely
         * screenshot-able. */
        CHECK(ls_scene_save(&d, SCENE_OUT));
        SceneDesc r;
        CHECK(ls_scene_load(&r, SCENE_OUT));
        CHECK(r.nlights == d.nlights);
        double after = probe(&r, v3(0.0, 0.0, 0.001));
        CHECK_NEAR(after, before, 1e-12);
        NOTE("edited scene re-solves to %.9g W/m^2 after a save/reload", after);
        check_pairing(&r);

        /* And the light we placed came back with the flux we asked for. */
        CHECK_NEAR(r.lights[id].flux_authored, 850.0, 1e-9);
        CHECK(r.lights[id].flux_in_lumens);

        ls_scene_desc_free(&d);
        ls_scene_desc_free(&r);
    }

    SECTION("clone is a deep copy");
    {
        SceneDesc a, b;
        if (!ls_scene_load(&a, SCENE_IN)) return;
        CHECK(ls_scene_clone(&a, &b));
        CHECK(b.lights != a.lights);
        CHECK(b.prims  != a.prims);
        CHECK(b.mats   != a.mats);
        CHECK(b.nlights == a.nlights);

        double before = b.lights[0].phi_e;
        /* Mutating the original must not touch the snapshot -- this is what
         * undo depends on. */
        a.lights[0].phi_e *= 3.0;
        ls_scene_update_light(&a, 0);
        CHECK_NEAR(b.lights[0].phi_e, before, 1e-12);
        CHECK(b.lights[0].phi_e != a.lights[0].phi_e);
        check_pairing(&b);

        ls_scene_desc_free(&a);
        ls_scene_desc_free(&b);
    }
}
