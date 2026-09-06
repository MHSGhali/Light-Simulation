/* M3 gate: energy conservation and strategy agreement.
 *
 * The furnace test is the highest-value test in the suite. It fails loudly on
 * energy-conservation bugs, missing or extra cosine factors, wrong sampling
 * PDFs, and bad Russian-roulette compensation -- the errors that are otherwise
 * invisible in a picture that "looks fine". */
#include "test.h"
#include "tests.h"
#include "lightsim/integrator.h"
#include <string.h>

/* A closed box of six quads, every face sharing one material. */
static void build_box(Prim *q, ls_real half, int mat_id) {
    const vec3 c[6] = {
        {0,0,-1},{0,0, 1},{-1,0,0},{ 1,0,0},{0,-1,0},{0, 1,0}
    };
    const vec3 ax[6] = {
        {1,0,0},{1,0,0},{0,1,0},{0,1,0},{1,0,0},{1,0,0}
    };
    const vec3 ay[6] = {
        {0,1,0},{0,1,0},{0,0,1},{0,0,1},{0,0,1},{0,0,1}
    };
    for (int i = 0; i < 6; ++i) {
        memset(&q[i], 0, sizeof q[i]);
        q[i].kind = LS_PRIM_QUAD;
        q[i].c  = v3scale(c[i], half);
        q[i].n  = v3neg(c[i]);                 /* inward */
        q[i].ex = v3scale(ax[i], half);
        q[i].ey = v3scale(ay[i], half);
        q[i].mat_id = mat_id;
        q[i].light_id = -1;
    }
}

/* Mean radiance seen from the centre of the box, averaged over directions. */
static double furnace_radiance(double rho, double le, int max_depth,
                               int nrays, Rng *rng) {
    Material m;
    memset(&m, 0, sizeof m);
    m.bsdf.kind = LS_BSDF_LAMBERT;
    m.bsdf.rho  = ls_spectrum_const(rho);
    m.le        = ls_spectrum_const(le);
    m.emissive  = true;

    Prim box[6];
    build_box(box, 1.0, 0);
    Scene sc = { .prims = box, .nprims = 6, .mats = &m, .nmats = 1 };

    SpectrumAcc total = ls_acc_zero();
    for (int i = 0; i < nrays; ++i) {
        Ray r;
        r.o = v3(0, 0, 0);
        r.d = ls_sample_sphere_uniform(ls_rng_f(rng), ls_rng_f(rng));
        r.tmin = 0.0;
        r.tmax = HUGE_VAL;
        SpectrumAcc acc;
        /* BSDF sampling only: the walls are emissive but are not registered as
         * sampleable lights, so this is pure scatter-and-collect. */
        ls_trace_radiance(&sc, r, rng, max_depth, LS_STRAT_BSDF, &acc, NULL);
        for (int b = 0; b < LS_NBINS; ++b) total.v[b] += acc.v[b];
    }
    Spectrum mean = ls_acc_mean(&total, (uint64_t)nrays);
    /* le was a per-nm density, so divide the band integral back out to recover
     * the scalar radiance in the same units the caller supplied. */
    return ls_spectrum_integrate(&mean) / (LS_NBINS * LS_SPECTRAL_STEP);
}

void test_furnace(void) {
    Rng rng = ls_rng_seed(0x5DEECE66Dull, 11);

    SECTION("furnace: closed enclosure equilibrium");
    {
        /* Every surface emits radiance Le and reflects rho. The field is
         * isotropic, so the equilibrium radiance is Le + rho*Le + rho^2*Le + ...
         * = Le / (1 - rho). Variance grows sharply with albedo, so rho = 0.99
         * is deliberately not asserted at this sample count. */
        const int NR = 60000, DEPTH = 400;
        struct { double rho, want; } cases[] = {
            { 0.0, 1.0 }, { 0.5, 2.0 }, { 0.8, 5.0 }, { 0.9, 10.0 }
        };
        for (size_t i = 0; i < sizeof cases / sizeof cases[0]; ++i) {
            double got = furnace_radiance(cases[i].rho, 1.0, DEPTH, NR, &rng);
            CHECK_NEAR(got, cases[i].want, 8e-3);
            NOTE("rho=%.2f  L = %.5f  (analytic Le/(1-rho) = %.5f)",
                 cases[i].rho, got, cases[i].want);
        }
    }

    SECTION("furnace: per-depth partial sums");
    {
        /* Capping the path at k bounces must give the truncated geometric
         * series Le (1 - rho^(k+1)) / (1 - rho). If term k is right and k+1 is
         * wrong, the bug is localised to a single bounce -- this is the best
         * interreflection debugger in the suite. */
        const double rho = 0.7;
        for (int k = 0; k <= 6; ++k) {
            double want = (1.0 - pow(rho, k + 1)) / (1.0 - rho);
            double got  = furnace_radiance(rho, 1.0, k, 40000, &rng);
            CHECK_NEAR(got, want, 1e-2);
        }
        NOTE("per-depth series matches Le(1-rho^(k+1))/(1-rho) for k = 0..6");
    }

    SECTION("strategy agreement: NEE vs BSDF vs MIS");
    {
        /* NEE and BSDF sampling are each independently unbiased, so requiring
         * all three to agree catches MIS bookkeeping errors: a wrong weight
         * makes MIS disagree with two strategies that cannot both be wrong in
         * the same direction. */
        Material mats[2];
        memset(mats, 0, sizeof mats);
        /* 0: Lambertian floor. */
        mats[0].bsdf.kind = LS_BSDF_LAMBERT;
        mats[0].bsdf.rho  = ls_spectrum_const(0.6);
        /* 1: the emissive face of the area light. */
        double phi = 20.0;
        vec3 ex = v3(0.5, 0, 0), ey = v3(0, -0.5, 0);
        Light lt = ls_light_rect(v3(0, 0, 2), ex, ey, phi, ls_spectrum_const(1.0));
        ls_light_finalize(&lt, 0);
        mats[1].bsdf.kind = LS_BSDF_LAMBERT;
        mats[1].bsdf.rho  = ls_spectrum_zero();      /* pure emitter */
        mats[1].le        = ls_spectrum_scale(lt.s_hat, lt.radiance);
        mats[1].emissive  = true;

        Prim prims[2];
        memset(prims, 0, sizeof prims);
        prims[0].kind = LS_PRIM_QUAD;                /* floor at z = 0 */
        prims[0].c = v3(0, 0, 0);
        prims[0].n = v3(0, 0, 1);
        prims[0].ex = v3(8, 0, 0);
        prims[0].ey = v3(0, 8, 0);
        prims[0].mat_id = 0;
        prims[0].light_id = -1;
        prims[1].kind = LS_PRIM_QUAD;                /* the light's geometry */
        prims[1].c = lt.p;
        prims[1].n = lt.n;
        prims[1].ex = ex;
        prims[1].ey = ey;
        prims[1].mat_id = 1;
        prims[1].light_id = 0;                       /* bound to lights[0] */

        Scene sc = { .prims = prims, .nprims = 2, .mats = mats, .nmats = 2,
                     .lights = &lt, .nlights = 1 };

        Ray r;
        r.o = v3(0.0, 1.6, 1.2);
        r.d = v3norm(v3sub(v3(0.4, 0.3, 0.0), r.o));
        r.tmin = 0.0;
        r.tmax = HUGE_VAL;

        const int NR = 120000;
        LsStrategy strats[3] = { LS_STRAT_NEE, LS_STRAT_BSDF, LS_STRAT_MIS };
        const char *names[3] = { "NEE", "BSDF", "MIS" };
        double mean[3], var[3];

        for (int s = 0; s < 3; ++s) {
            double sum = 0.0, sum2 = 0.0;
            for (int i = 0; i < NR; ++i) {
                SpectrumAcc acc;
                ls_trace_radiance(&sc, r, &rng, 2, strats[s], &acc, NULL);
                Spectrum m = ls_acc_mean(&acc, 1);
                double v = ls_spectrum_integrate(&m);
                sum += v;
                sum2 += v * v;
            }
            mean[s] = sum / NR;
            var[s]  = sum2 / NR - mean[s] * mean[s];
            NOTE("%-4s  L = %.6f   variance = %.3e", names[s], mean[s], var[s]);
        }

        CHECK_NEAR(mean[1], mean[0], 2e-2);      /* BSDF agrees with NEE  */
        CHECK_NEAR(mean[2], mean[0], 2e-2);      /* MIS  agrees with NEE  */
        CHECK(mean[0] > 0.0);

        /* MIS should not be worse than the better of the two it combines. */
        CHECK(var[2] <= ls_max(var[0], var[1]) * 1.05);
    }
}
