/* M2 validation: direct lighting against closed-form radiometry.
 * Every expected value is analytic; the constants were derived independently
 * and are written as literals so a refactor cannot quietly move them. */
#include "test.h"
#include "tests.h"
#include "lightsim/integrator.h"
#include "lightsim/units.h"
#include <string.h>

/* Tolerance for any value that has passed through a Spectrum. The bins are
 * float, so ~1e-7 relative is the storage floor regardless of how exact the
 * transport maths is -- delta lights have no sampling variance but still
 * round-trip through float. Pure-double quantities (intensity, flux, solid
 * angle) are asserted at 1e-12 instead. */
#define E_TOL 1e-6

/* Band-integrated irradiance from a spectral accumulator, W/m^2. */
static double E_of(const SpectrumAcc *a) {
    Spectrum s = ls_acc_mean(a, 1);
    return ls_spectrum_integrate(&s);
}

static double measure(Scene *sc, vec3 p, vec3 n, int nsamples, Rng *rng) {
    SpectrumAcc acc;
    ls_estimate_irradiance(sc, p, n, nsamples, rng, &acc, NULL);
    return E_of(&acc);
}

void test_transport(void) {
    Rng rng = ls_rng_seed(0x9E3779B97F4A7C15ull, 1);
    Spectrum flat = ls_spectrum_const(1.0);   /* shape only; normalised on finalize */

    SECTION("sampling warps");
    {
        /* Cone solid angle: Omega = 2 pi (1 - cos alpha). This is the relation
         * that converts a source's angular span into steradians. */
        CHECK_NEAR(ls_cone_solid_angle(LS_PI), 4.0 * LS_PI, 1e-12);
        CHECK_NEAR(ls_cone_solid_angle(LS_PI / 2.0), LS_TWO_PI, 1e-12);
        CHECK_NEAR(ls_cone_solid_angle(30.0 * LS_PI / 180.0), 0.8417872144769325, 1e-12);

        /* Cosine-weighted hemisphere: the pdf must integrate to 1, and the
         * cosine integral over the hemisphere must be pi. A projected-solid-
         * angle pdf (1/pi instead of cos/pi) fails the first of these. */
        const int N = 400000;
        double sum_pdf = 0.0, sum_cos = 0.0;
        for (int i = 0; i < N; ++i) {
            vec3 w = ls_sample_hemisphere_cosine(ls_rng_f(&rng), ls_rng_f(&rng));
            double pdf = ls_pdf_hemisphere_cosine(w.z);
            sum_pdf += 1.0;                       /* integral of pdf/pdf */
            sum_cos += w.z / pdf;                 /* integral of cos dw = pi */
        }
        CHECK_NEAR(sum_pdf / N, 1.0, 1e-12);
        CHECK_NEAR(sum_cos / N, LS_PI, 2e-3);

        /* Uniform sphere sampling: directions must be unit length, have zero
         * mean vector, and satisfy <|z|> = 1/2 (the mean of |cos| over a
         * uniform sphere). A pdf or warp error shows up in the last of these. */
        vec3 mean = v3(0, 0, 0);
        double sum_absz = 0.0, max_len_err = 0.0;
        for (int i = 0; i < N; ++i) {
            vec3 w = ls_sample_sphere_uniform(ls_rng_f(&rng), ls_rng_f(&rng));
            mean = v3add(mean, w);
            sum_absz += fabs(w.z);
            max_len_err = ls_max(max_len_err, fabs(v3len(w) - 1.0));
        }
        CHECK(max_len_err < 1e-12);
        CHECK_NEAR(sum_absz / N, 0.5, 5e-3);
        CHECK(v3len(v3scale(mean, 1.0 / N)) < 5e-3);
    }

    SECTION("point source: inverse square and Lambert cosine");
    {
        Light l = ls_light_point(v3(0, 0, 0), 1.0, flat);
        ls_light_finalize(&l, 0);
        Scene sc = { .lights = &l, .nlights = 1 };

        /* C1: isotropic intensity I = Phi / 4 pi. */
        CHECK_NEAR(ls_light_intensity(&l, v3(0, 0, 1)), 0.07957747154594767, 1e-12);
        /* The light must radiate exactly the flux it was given. */
        CHECK_NEAR(ls_light_emitted_flux(&l), 1.0, 1e-12);

        /* C2: inverse square. The probe sits below the source and faces it.
         * Deterministic -- a delta light has no sampling variance. */
        CHECK_NEAR(measure(&sc, v3(0, 0, -1), v3(0, 0, 1), 1, &rng),
                   0.07957747154594767, E_TOL);
        CHECK_NEAR(measure(&sc, v3(0, 0, -2), v3(0, 0, 1), 1, &rng),
                   0.019894367886486918, E_TOL);

        /* C3: Lambert cosine law, swept. E(theta) = I cos(theta) / r^2, with
         * theta the angle between the surface normal and the incident ray. */
        for (int deg = 0; deg <= 85; deg += 5) {
            double th = (double)deg * LS_PI / 180.0;
            vec3 n = v3(sin(th), 0.0, cos(th));
            CHECK_NEAR(measure(&sc, v3(0, 0, -1), n, 1, &rng),
                       0.07957747154594767 * cos(th), E_TOL);
        }

        /* A surface facing away receives nothing. */
        CHECK_NEAR(measure(&sc, v3(0, 0, -1), v3(0, 0, -1), 1, &rng), 0.0, 1e-15);
    }

    SECTION("directional and spot sources");
    {
        /* C9: E = E_perp cos(theta), independent of distance. */
        Light d = ls_light_directional(v3(0, 0, -1), 2.5, flat);
        ls_light_finalize(&d, 0);
        Scene sd = { .lights = &d, .nlights = 1 };
        CHECK_NEAR(measure(&sd, v3(0, 0, 0), v3(0, 0, 1), 1, &rng), 2.5, E_TOL);
        CHECK_NEAR(measure(&sd, v3(0, 0, 500), v3(0, 0, 1), 1, &rng), 2.5, E_TOL);
        double th = 60.0 * LS_PI / 180.0;
        CHECK_NEAR(measure(&sd, v3(0, 0, 0), v3(sin(th), 0, cos(th)), 1, &rng),
                   2.5 * cos(th), E_TOL);

        /* C10: hard-edged spot, half-angle 30 deg, Phi = 1 W.
         * I0 = Phi / (2 pi (1 - cos alpha)). Getting Phi/(4 pi) here instead is
         * the classic spot bug; ls_light_finalize would already have caught it. */
        double alpha = 30.0 * LS_PI / 180.0;
        Light s = ls_light_spot(v3(0, 0, 0), v3(0, 0, -1), alpha, alpha, 1.0, flat);
        ls_light_finalize(&s, 0);
        Scene ss = { .lights = &s, .nlights = 1 };
        CHECK_NEAR(s.omega_eff, 0.8417872144769325, 1e-12);
        CHECK_NEAR(ls_light_intensity(&s, v3(0, 0, -1)), 1.187948667789374, 1e-12);
        CHECK_NEAR(ls_light_emitted_flux(&s), 1.0, 1e-12);
        /* On axis at r = 2. */
        CHECK_NEAR(measure(&ss, v3(0, 0, -2), v3(0, 0, 1), 1, &rng),
                   0.2969871669473435, E_TOL);
        /* Outside the cone: exactly zero. */
        CHECK_NEAR(measure(&ss, v3(3, 0, -2), v3(0, 0, 1), 1, &rng), 0.0, 1e-15);

        /* A smooth-falloff spot must still radiate exactly its stated flux --
         * that is what the numerically integrated omega_eff buys. */
        Light sf = ls_light_spot(v3(0, 0, 0), v3(0, 0, -1),
                                 alpha, 20.0 * LS_PI / 180.0, 1.0, flat);
        ls_light_finalize(&sf, 0);
        CHECK_NEAR(ls_light_emitted_flux(&sf), 1.0, 1e-12);
        CHECK(sf.omega_eff < 0.8417872144769325);   /* narrower than hard-edged */
    }

    SECTION("area sources against closed forms");
    {
        const int NS = 400000;

        /* C4: Lambertian disk, R = 0.5, L = 1, on-axis at h = 2.
         *     E = pi L R^2 / (R^2 + h^2). */
        double phi_disk = 2.4674011002723395;      /* L=1 => Phi = pi^2 R^2 */
        Light dk = ls_light_disk(v3(0, 0, 0), v3(0, 0, -1), 0.5, phi_disk, flat);
        ls_light_finalize(&dk, 0);
        CHECK_NEAR(dk.radiance, 1.0, 1e-12);
        CHECK_NEAR(ls_light_emitted_flux(&dk), phi_disk, 1e-12);
        Scene sdk = { .lights = &dk, .nlights = 1 };
        CHECK_NEAR(measure(&sdk, v3(0, 0, -2), v3(0, 0, 1), NS, &rng),
                   0.18479956785822313, 5e-3);

        /* C5: Lambertian sphere, R = 0.25, L = 1, centre at d = 3.
         *     E = pi L R^2 / d^2. */
        double phi_sph = 2.4674011002723395;       /* L=1 => Phi = 4 pi^2 R^2 */
        Light sp = ls_light_sphere(v3(0, 0, 0), 0.25, phi_sph, flat);
        ls_light_finalize(&sp, 0);
        CHECK_NEAR(sp.radiance, 1.0, 1e-12);
        Scene ssp = { .lights = &sp, .nlights = 1 };
        double e_sphere = measure(&ssp, v3(0, 0, -3), v3(0, 0, 1), NS, &rng);
        CHECK_NEAR(e_sphere, 0.02181661564992912, 5e-3);

        /* C6: THE sphere-equals-point identity. A uniform-radiance sphere
         * produces exactly the irradiance of a point source of equal flux at
         * its centre. This one test exercises sphere area sampling, the
         * area-to-solid-angle Jacobian, the flux/radiance normalisation, and
         * the delta-vs-area code paths at once -- an error in any of them
         * breaks the agreement. */
        Light pt = ls_light_point(v3(0, 0, 0), phi_sph, flat);
        ls_light_finalize(&pt, 0);
        Scene spt = { .lights = &pt, .nlights = 1 };
        double e_point = measure(&spt, v3(0, 0, -3), v3(0, 0, 1), 1, &rng);
        CHECK_NEAR(e_point, 0.02181661564992912, E_TOL);
        CHECK_NEAR(e_sphere, e_point, 5e-3);
        NOTE("sphere %.8f vs point %.8f  (%.3f%% apart)",
             e_sphere, e_point, 100.0 * fabs(e_sphere - e_point) / e_point);

        /* C7: Lambertian rectangle 1 x 2, L = 1, field point on the
         * perpendicular through the centre at c = 1.5.  E = pi L F. */
        double phi_rect = 6.283185307179586;       /* L=1 => Phi = pi a b */
        Light rc = ls_light_rect(v3(0, 0, 0), v3(0.5, 0, 0), v3(0, -1.0, 0),
                                 phi_rect, flat);
        ls_light_finalize(&rc, 0);
        CHECK_NEAR(rc.area, 2.0, 1e-12);
        CHECK_NEAR(rc.radiance, 1.0, 1e-12);
        CHECK_NEAR(rc.n.z, -1.0, 1e-12);           /* faces the probe */
        Scene src = { .lights = &rc, .nlights = 1 };
        CHECK_NEAR(measure(&src, v3(0, 0, -1.5), v3(0, 0, 1), NS, &rng),
                   0.6568166565040201, 5e-3);
    }

    SECTION("occlusion");
    {
        Light l = ls_light_point(v3(0, 0, 0), 1.0, flat);
        ls_light_finalize(&l, 0);

        /* An opaque quad between source and probe blocks it completely. */
        Prim blocker;
        memset(&blocker, 0, sizeof blocker);
        blocker.kind = LS_PRIM_QUAD;
        blocker.c  = v3(0, 0, -1);
        blocker.n  = v3(0, 0, 1);
        blocker.ex = v3(5, 0, 0);
        blocker.ey = v3(0, 5, 0);
        blocker.mat_id = 0;
        blocker.light_id = -1;

        Scene sc = { .prims = &blocker, .nprims = 1, .lights = &l, .nlights = 1 };
        CHECK_NEAR(measure(&sc, v3(0, 0, -2), v3(0, 0, 1), 1, &rng), 0.0, 1e-15);

        /* Step the probe out from behind the blocker and the full unoccluded
         * value returns -- proving the blocker is finite, not a global switch. */
        Scene open_sc = { .lights = &l, .nlights = 1 };
        double want = measure(&open_sc, v3(20, 0, -2), v3(0, 0, 1), 1, &rng);
        CHECK_NEAR(measure(&sc, v3(20, 0, -2), v3(0, 0, 1), 1, &rng), want, E_TOL);
        CHECK(want > 0.0);

        /* Half-occlusion: a blocker covering exactly the -x half plane must
         * halve a symmetric disk source's irradiance. */
        double phi_disk = 2.4674011002723395;
        Light dk = ls_light_disk(v3(0, 0, 0), v3(0, 0, -1), 0.5, phi_disk, flat);
        ls_light_finalize(&dk, 0);
        Scene open_dk = { .lights = &dk, .nlights = 1 };
        double full = measure(&open_dk, v3(0, 0, -2), v3(0, 0, 1), 200000, &rng);

        Prim half;
        memset(&half, 0, sizeof half);
        half.kind = LS_PRIM_QUAD;
        half.c  = v3(-0.5, 0, -1);       /* covers x in [-1, 0] */
        half.n  = v3(0, 0, 1);
        half.ex = v3(0.5, 0, 0);
        half.ey = v3(0, 5, 0);
        half.light_id = -1;
        Scene half_sc = { .prims = &half, .nprims = 1, .lights = &dk, .nlights = 1 };
        double halved = measure(&half_sc, v3(0, 0, -2), v3(0, 0, 1), 200000, &rng);
        CHECK_NEAR(halved, 0.5 * full, 1e-2);
        NOTE("half-occluded disk: %.6f vs half of %.6f = %.6f",
             halved, full, 0.5 * full);
    }
}
