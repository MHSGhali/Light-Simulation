/* E4 validation: the property inspector's model.
 *
 * Kept free of SDL so the field lists and the edit semantics can be checked
 * without a window. The thing most worth guarding is that one set of stored
 * values really does serve all three tiers -- the tier must change what is
 * SHOWN, never what is simulated. */
#include "test.h"
#include "tests.h"
#include "../viewer/inspect.h"
#include "lightsim/units.h"
#include <string.h>

static int find(const Field *f, int n, FieldId id) {
    for (int i = 0; i < n; ++i) if (f[i].id == id) return i;
    return -1;
}

static double get(const SceneDesc *d, int li, LsTier t, FieldId id) {
    Field f[LS_INSPECT_MAX];
    int n = ls_inspect_fields(d, li, -1, t, f, LS_INSPECT_MAX);
    int k = find(f, n, id);
    return k >= 0 ? f[k].value : -1e30;
}

void test_inspect(void) {
    SceneDesc d;
    if (!ls_scene_load(&d, "scenes/workcell.scene")) {
        printf("  SKIP: cannot load the demo scene\n");
        return;
    }

    SECTION("tiers reveal, never change");
    {
        Field f[LS_INSPECT_MAX];
        int ns = ls_inspect_fields(&d, 0, -1, LS_TIER_SIMPLE, f, LS_INSPECT_MAX);
        int na = ls_inspect_fields(&d, 0, -1, LS_TIER_ADVANCED, f, LS_INSPECT_MAX);
        int nc = ls_inspect_fields(&d, 0, -1, LS_TIER_SCIENTIFIC, f, LS_INSPECT_MAX);
        CHECK(ns > 0);
        CHECK(na > ns);
        CHECK(nc > na);
        NOTE("fields per tier: simple %d, advanced %d, scientific %d", ns, na, nc);

        /* Every SIMPLE field must still be present, with the same value, at the
         * deeper tiers -- a deeper tier adds rows, it does not restate them. */
        Field s[LS_INSPECT_MAX];
        int nsimple = ls_inspect_fields(&d, 0, -1, LS_TIER_SIMPLE, s, LS_INSPECT_MAX);
        nc = ls_inspect_fields(&d, 0, -1, LS_TIER_SCIENTIFIC, f, LS_INSPECT_MAX);
        for (int i = 0; i < nsimple; ++i) {
            if (s[i].heading) continue;
            int k = find(f, nc, s[i].id);
            CHECK(k >= 0);
            if (k >= 0) CHECK_NEAR(f[k].value, s[i].value, 1e-12);
        }

        /* The simple tier must NOT show watts; the scientific tier must. */
        CHECK(find(s, nsimple, FLD_L_W) < 0);
        CHECK(find(f, nc, FLD_L_W) >= 0);
        CHECK(find(s, nsimple, FLD_L_LM) >= 0);   /* lumens lead everywhere */
    }

    SECTION("lumens are what the user asked for");
    {
        /* Set 850 lm and the light must actually deliver 850 lm. */
        CHECK(ls_inspect_set(&d, 0, -1, FLD_L_LM, 850.0));
        CHECK_NEAR(get(&d, 0, LS_TIER_SIMPLE, FLD_L_LM), 850.0, 1e-6);

        double w_before = get(&d, 0, LS_TIER_SCIENTIFIC, FLD_L_W);
        double eff_before = get(&d, 0, LS_TIER_SCIENTIFIC, FLD_L_EFFICACY);

        /* Change the colour temperature. The lumens the user asked for must be
         * preserved; the WATTS required change, because efficacy changed. That
         * is the whole point of storing spectra rather than a colour. */
        CHECK(ls_inspect_set(&d, 0, -1, FLD_L_CCT, 2700.0));
        CHECK_NEAR(get(&d, 0, LS_TIER_SIMPLE, FLD_L_LM), 850.0, 1e-6);
        double w_after = get(&d, 0, LS_TIER_SCIENTIFIC, FLD_L_W);
        double eff_after = get(&d, 0, LS_TIER_SCIENTIFIC, FLD_L_EFFICACY);
        CHECK(w_after > w_before);          /* warmer source, worse efficacy */
        CHECK(eff_after < eff_before);
        NOTE("850 lm at 5000K needs %.3f W (%.1f lm/W); at 2700K, %.3f W (%.1f lm/W)",
             w_before, eff_before, w_after, eff_after);

        /* And the efficacy reported is exactly lumens over watts. Asserted at
         * the float-Spectrum storage floor, not machine precision -- the lumens
         * are recomputed by integrating a float spectrum, so they come back as
         * 850 to about 1e-7, not to 1e-15. */
        CHECK_NEAR(eff_after, 850.0 / w_after, 1e-6);

        /* Cross-check against the blackbody efficacy computed from first
         * principles elsewhere in this suite: a 2700 K Planckian radiator has a
         * band luminous efficacy of 112.66 lm/W. The inspector's colour-temperature
         * switch is meant to pick a Planckian radiator below 4000 K, and this is
         * what proves it actually did. */
        CHECK_NEAR(eff_after, 112.663, 2e-3);
    }

    SECTION("editing a light keeps its geometry in step");
    {
        CHECK(ls_inspect_set(&d, 0, -1, FLD_L_SIZEU, 0.2));
        CHECK_NEAR(get(&d, 0, LS_TIER_ADVANCED, FLD_L_SIZEU), 0.2, 1e-9);
        int pi = ls_scene_light_prim(&d, 0);
        CHECK(pi >= 0);
        CHECK_NEAR(2.0 * v3len(d.prims[pi].ex), 0.2, 1e-9);

        CHECK(ls_inspect_set(&d, 0, -1, FLD_L_Z, 0.42));
        CHECK_NEAR(d.lights[0].p.z, 0.42, 1e-12);
        pi = ls_scene_light_prim(&d, 0);
        CHECK_NEAR(d.prims[pi].c.z, 0.42, 1e-12);
    }

    SECTION("beam angle round-trips through the inspector");
    {
        CHECK(ls_inspect_set(&d, 0, -1, FLD_L_KIND, (double)LS_LIGHT_SPOT));
        for (double beam = 10.0; beam <= 150.0; beam += 20.0) {
            CHECK(ls_inspect_set(&d, 0, -1, FLD_L_BEAM, beam));
            CHECK_NEAR(get(&d, 0, LS_TIER_SIMPLE, FLD_L_BEAM), beam, 1e-9);
            /* And the flux self-check inside ls_light_finalize still holds,
             * which is what stops a beam edit quietly changing the output. */
            CHECK_NEAR(ls_light_emitted_flux(&d.lights[0]), d.lights[0].phi_e, 1e-9);
        }
        double field = get(&d, 0, LS_TIER_SIMPLE, FLD_L_FIELD);
        CHECK(field > get(&d, 0, LS_TIER_SIMPLE, FLD_L_BEAM));
    }

    SECTION("part and material editing");
    {
        /* Find a non-emissive primitive to edit. */
        int pi = -1;
        for (int i = 0; i < d.nprims; ++i)
            if (d.prims[i].light_id < 0) { pi = i; break; }
        CHECK(pi >= 0);

        Field f[LS_INSPECT_MAX];
        int n = ls_inspect_fields(&d, -1, pi, LS_TIER_ADVANCED, f, LS_INSPECT_MAX);
        CHECK(n > 0);
        CHECK(find(f, n, FLD_M_ALBEDO) >= 0);

        CHECK(ls_inspect_set(&d, -1, pi, FLD_M_ALBEDO, 0.62));
        CHECK_NEAR(ls_spectrum_mean(&d.mats[d.prims[pi].mat_id].bsdf.rho), 0.62, 1e-6);

        /* Convert it to a metal and back. */
        CHECK(ls_inspect_set(&d, -1, pi, FLD_M_KIND, 1.0));
        CHECK(d.mats[d.prims[pi].mat_id].bsdf.kind == LS_BSDF_CONDUCTOR);
        CHECK(ls_inspect_set(&d, -1, pi, FLD_M_ROUGH, 0.25));
        CHECK_NEAR(d.mats[d.prims[pi].mat_id].bsdf.alpha, 0.25, 1e-12);
        CHECK(ls_inspect_set(&d, -1, pi, FLD_M_KIND, 0.0));
        CHECK(d.mats[d.prims[pi].mat_id].bsdf.kind == LS_BSDF_LAMBERT);

        /* Moving a part works. */
        CHECK(ls_inspect_set(&d, -1, pi, FLD_P_X, 0.11));
        CHECK_NEAR(d.prims[pi].c.x, 0.11, 1e-12);
    }

    SECTION("formatting");
    {
        Field f = { 0 };
        f.value = 1234.5;
        char buf[32];
        ls_inspect_format(&f, buf, sizeof buf);
        CHECK(strcmp(buf, "1235") == 0 || strcmp(buf, "1234") == 0);
        f.value = 0.25;
        ls_inspect_format(&f, buf, sizeof buf);
        CHECK(strcmp(buf, "0.250") == 0);
    }

    ls_scene_desc_free(&d);
}
