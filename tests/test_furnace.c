/* M3 gate: energy conservation and strategy agreement.
 *
 * The furnace test is the highest-value test in the suite. It fails loudly on
 * energy-conservation bugs, missing or extra cosine factors, wrong sampling
 * PDFs, and bad Russian-roulette compensation -- the errors that are otherwise
 * invisible in a picture that "looks fine". */
#include "test.h"
#include "tests.h"
#include "lightsim/integrator.h"
#include "lightsim/mesh.h"
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

/* The SAME box, as 12 triangles in one mesh.
 *
 * `inward` reverses each triangle's winding. It matters: emission is gated on
 * the raw geometric normal (integrator.c), never on the ray-facing one, so a
 * box wound outward is a box whose walls do not emit toward the inside. That
 * makes this a two-sided check on the mesh intersector's normals rather than a
 * one-sided check that a plausible number came out. */
static Mesh *build_box_mesh(ls_real h, bool inward) {
    vec3 v[8] = {
        v3(-h,-h,-h), v3( h,-h,-h), v3( h, h,-h), v3(-h, h,-h),
        v3(-h,-h, h), v3( h,-h, h), v3( h, h, h), v3(-h, h, h)
    };
    static const int out[36] = {
        0,3,2, 0,2,1,   4,5,6, 4,6,7,      /* -z, +z */
        0,1,5, 0,5,4,   3,7,6, 3,6,2,      /* -y, +y */
        0,4,7, 0,7,3,   1,2,6, 1,6,5       /* -x, +x */
    };
    int idx[36];
    for (int i = 0; i < 12; ++i) {
        idx[i*3+0] = out[i*3+0];
        idx[i*3+1] = inward ? out[i*3+2] : out[i*3+1];
        idx[i*3+2] = inward ? out[i*3+1] : out[i*3+2];
    }
    return ls_mesh_build(v, 8, idx, 12);
}

/* Equilibrium radiance of the mesh box, by the same measurement. One Prim of
 * kind LS_PRIM_MESH replaces the six quads; everything downstream is identical,
 * which is the point -- the answer must not move. */
static double furnace_radiance_mesh(double rho, double le, int max_depth,
                                    int nrays, Rng *rng, bool inward) {
    Material m;
    memset(&m, 0, sizeof m);
    m.bsdf.kind = LS_BSDF_LAMBERT;
    m.bsdf.rho  = ls_spectrum_const(rho);
    m.le        = ls_spectrum_const(le);
    m.emissive  = true;

    Mesh *mesh = build_box_mesh(1.0, inward);
    if (!mesh) return -1.0;

    Prim p;
    memset(&p, 0, sizeof p);
    p.kind = LS_PRIM_MESH;
    p.c  = v3(0, 0, 0);
    p.ex = v3(1, 0, 0);
    p.ey = v3(0, 1, 0);
    p.n  = v3(0, 0, 1);
    p.mat_id = 0;
    p.light_id = -1;
    p.mesh_id = 0;

    Scene sc = { .prims = &p, .nprims = 1, .mats = &m, .nmats = 1,
                 .meshes = &mesh, .nmeshes = 1 };

    SpectrumAcc total = ls_acc_zero();
    for (int i = 0; i < nrays; ++i) {
        Ray r;
        r.o = v3(0, 0, 0);
        r.d = ls_sample_sphere_uniform(ls_rng_f(rng), ls_rng_f(rng));
        r.tmin = 0.0;
        r.tmax = HUGE_VAL;
        SpectrumAcc acc;
        ls_trace_radiance(&sc, r, rng, max_depth, LS_STRAT_BSDF, &acc, NULL);
        for (int b = 0; b < LS_NBINS; ++b) total.v[b] += acc.v[b];
    }
    ls_mesh_release(mesh);
    Spectrum mean = ls_acc_mean(&total, (uint64_t)nrays);
    return ls_spectrum_integrate(&mean) / (LS_NBINS * LS_SPECTRAL_STEP);
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

    SECTION("furnace: the same box as a triangle mesh");
    {
        /* The gate for the mesh intersector. Six analytic quads and twelve
         * triangles describe the same enclosure, so they must reach the same
         * equilibrium -- and reaching it exercises winding, geometric normals,
         * backface handling, self-intersection offsets and t ordering at once,
         * over hundreds of bounces, with no new physics to get wrong.
         *
         * Same seed sequence and sample count as the quad case above, so the
         * two numbers are comparable directly and not merely both plausible. */
        const int NR = 60000, DEPTH = 400;
        struct { double rho, want; } cases[] = {
            { 0.0, 1.0 }, { 0.5, 2.0 }, { 0.8, 5.0 }
        };
        for (size_t i = 0; i < sizeof cases / sizeof cases[0]; ++i) {
            Rng qr = ls_rng_seed(0x5DEECE66Dull, 900 + (uint64_t)i);
            Rng mr = ls_rng_seed(0x5DEECE66Dull, 900 + (uint64_t)i);
            double quads = furnace_radiance(cases[i].rho, 1.0, DEPTH, NR, &qr);
            double mesh  = furnace_radiance_mesh(cases[i].rho, 1.0, DEPTH, NR,
                                                 &mr, true);
            CHECK_NEAR(mesh, cases[i].want, 8e-3);
            CHECK_NEAR(mesh, quads, 2e-2);
            NOTE("rho=%.2f  quads %.5f  mesh %.5f  (analytic %.5f)",
                 cases[i].rho, quads, mesh, cases[i].want);
        }

        /* And the normals are not being quietly flipped. Wound the other way
         * the walls face away from the interior, emission is gated off, and the
         * enclosure goes black. An intersector that returned a ray-facing
         * normal would light this up and pass the test above regardless. */
        Rng br = ls_rng_seed(0x5DEECE66Dull, 77);
        double backwards = furnace_radiance_mesh(0.5, 1.0, 8, 20000, &br, false);
        CHECK_NEAR(backwards, 0.0, 1e-12);
        NOTE("outward-wound box collects %.3g -- emission is one-sided on ng",
             backwards);
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

    SECTION("one-sided emitters agree between strategies");
    {
        /* An area light emits from one face only. NEE enforces that through
         * ls_light_pdf_w's cos_y test; BSDF sampling has to agree, or the two
         * strategies disagree about the back of every light and MIS silently
         * loses energy. The regression this guards: the emitter gate used the
         * ray-facing shading normal, which is flipped on a backface hit and so
         * could never fire -- making emitters two-sided to scattering only. */
        Light lt = ls_light_rect(v3(0, 0, 0), v3(0.5, 0, 0), v3(0, -0.5, 0),
                                 100.0, ls_spectrum_const(1.0));
        ls_light_finalize(&lt, 0);
        CHECK_NEAR(lt.n.z, -1.0, 1e-12);              /* faces -z */

        Material m;
        memset(&m, 0, sizeof m);
        m.bsdf.kind = LS_BSDF_LAMBERT;
        m.bsdf.rho  = ls_spectrum_zero();
        m.le        = ls_spectrum_scale(lt.s_hat, lt.radiance);
        m.emissive  = true;

        Prim pr;
        memset(&pr, 0, sizeof pr);
        pr.kind = LS_PRIM_QUAD; pr.c = lt.p; pr.n = lt.n;
        pr.ex = lt.ex; pr.ey = lt.ey; pr.mat_id = 0; pr.light_id = 0;

        Scene sc = { .prims = &pr, .nprims = 1, .mats = &m, .nmats = 1,
                     .lights = &lt, .nlights = 1 };

        LsStrategy all[3] = { LS_STRAT_NEE, LS_STRAT_BSDF, LS_STRAT_MIS };
        for (int s = 0; s < 3; ++s) {
            /* Looking up at the emitting face: all three see the same radiance. */
            Ray front = { v3(0, 0, -1), v3(0, 0, 1), 0.0, HUGE_VAL };
            /* Looking down at the dark back: all three must see exactly zero. */
            Ray back  = { v3(0, 0,  1), v3(0, 0, -1), 0.0, HUGE_VAL };
            SpectrumAcc a;
            ls_trace_radiance(&sc, back, &rng, 1, all[s], &a, NULL);
            Spectrum L = ls_acc_mean(&a, 1);
            CHECK_NEAR(ls_spectrum_integrate(&L), 0.0, 1e-12);
            if (all[s] & LS_STRAT_BSDF) {
                ls_trace_radiance(&sc, front, &rng, 1, all[s], &a, NULL);
                L = ls_acc_mean(&a, 1);
                CHECK(ls_spectrum_integrate(&L) > 0.0);
            }
        }
        /* And the irradiance behind the panel is zero, which is what NEE says. */
        SpectrumAcc acc;
        ls_estimate_irradiance(&sc, v3(0, 0, 0.5), v3(0, 0, 1), 512, &rng, &acc, NULL);
        Spectrum E = ls_acc_mean(&acc, 1);
        CHECK_NEAR(ls_spectrum_integrate(&E), 0.0, 1e-12);
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
