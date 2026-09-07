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
