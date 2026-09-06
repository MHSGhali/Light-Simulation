/* M3 validation: surface scattering. */
#include "test.h"
#include "tests.h"
#include "lightsim/bsdf.h"
#include "lightsim/integrator.h"
#include <string.h>

void test_bsdf(void) {
    Rng rng = ls_rng_seed(0xDEADBEEFCAFEull, 7);

    SECTION("Fresnel");
    {
        /* Dielectric normal incidence: R0 = ((n-1)/(n+1))^2 = 0.04 for n=1.5. */
        CHECK_NEAR(ls_fresnel_dielectric(1.0, 1.5), 0.04, 1e-12);
        /* Grazing incidence goes to unity for any interface. */
        CHECK_NEAR(ls_fresnel_dielectric(1e-9, 1.5), 1.0, 1e-6);

        /* At Brewster's angle the p-polarised term vanishes, leaving Rs/2. */
        double n = 1.5, thB = atan(n);
        CHECK_NEAR(ls_fresnel_dielectric(cos(thB), n), 0.07396449704142013, 1e-9);

        /* Conductor: aluminium at 550 nm, n = 0.958, k = 6.69.
         * R = ((n-1)^2 + k^2) / ((n+1)^2 + k^2) at normal incidence. */
        double nn = 0.958, kk = 6.69;
        double r0 = ((nn - 1) * (nn - 1) + kk * kk) / ((nn + 1) * (nn + 1) + kk * kk);
        CHECK_NEAR(ls_fresnel_conductor(1.0, nn, kk), r0, 1e-9);
        CHECK_NEAR(r0, 0.9211358, 1e-6);
        NOTE("aluminium normal-incidence reflectance at 550 nm = %.4f", r0);

        /* Reflectance at true grazing exceeds normal incidence. But it is NOT
         * monotone in between: an absorbing medium has a pseudo-Brewster
         * minimum, and aluminium's sits near cos(theta) = 0.2, where R dips to
         * ~0.870 below its normal-incidence 0.921 before climbing to ~0.98.
         * Asserting monotonicity here would be asserting the wrong physics. */
        CHECK(ls_fresnel_conductor(0.01, nn, kk) > r0);
        CHECK(ls_fresnel_conductor(0.20, nn, kk) < r0);
        CHECK(ls_fresnel_conductor(0.01, nn, kk) > ls_fresnel_conductor(0.20, nn, kk));
        NOTE("Al pseudo-Brewster dip: R(cos=0.2)=%.4f < R(cos=1)=%.4f, R(cos=0.01)=%.4f",
             ls_fresnel_conductor(0.20, nn, kk), r0, ls_fresnel_conductor(0.01, nn, kk));
        for (int i = 0; i <= 10; ++i) {
            double c = (double)i / 10.0;
            double r = ls_fresnel_conductor(c, nn, kk);
            CHECK(r >= 0.0 && r <= 1.0);
        }
    }

    SECTION("GGX distribution");
    {
        /* The normal distribution must satisfy  integral D(m) (n.m) dw = 1. */
        for (double alpha = 0.1; alpha <= 0.81; alpha += 0.35) {
            const int N = 300000;
            double sum = 0.0;
            for (int i = 0; i < N; ++i) {
                vec3 m = ls_sample_hemisphere_cosine(ls_rng_f(&rng), ls_rng_f(&rng));
                double pdf = ls_pdf_hemisphere_cosine(m.z);
                if (pdf <= 0.0) continue;
                sum += ls_ggx_d(m.z, alpha) * m.z / pdf;
            }
            CHECK_NEAR(sum / N, 1.0, 1e-2);
        }

        /* Masking is bounded and monotone in the expected direction. */
        vec3 straight = v3(0, 0, 1);
        CHECK_NEAR(ls_ggx_g1(straight, 0.3), 1.0, 1e-9);
        double th = 70.0 * LS_PI / 180.0;
        vec3 oblique = v3(sin(th), 0, cos(th));
        CHECK(ls_ggx_g1(oblique, 0.3) < 1.0);
        CHECK(ls_ggx_g1(oblique, 0.3) > 0.0);
    }

    SECTION("BSDF contracts");
    {
        Bsdf lam;
        memset(&lam, 0, sizeof lam);
        lam.kind = LS_BSDF_LAMBERT;
        lam.rho = ls_spectrum_const(1.0);   /* white */

        /* Lambert white furnace: integral of f cos(theta) dw over the
         * hemisphere must be exactly the albedo. Anything other than rho/pi in
         * ls_bsdf_eval, or a projected-solid-angle pdf, breaks this. */
        const int N = 400000;
        vec3 wo = v3norm(v3(0.3, 0.1, 0.8));
        double sum = 0.0;
        for (int i = 0; i < N; ++i) {
            vec3 wi = ls_sample_hemisphere_cosine(ls_rng_f(&rng), ls_rng_f(&rng));
            double pdf = ls_pdf_hemisphere_cosine(wi.z);
            if (pdf <= 0.0) continue;
            Spectrum f;
            ls_bsdf_eval(&lam, wo, wi, &f);
            sum += (double)f.v[0] * wi.z / pdf;
        }
        CHECK_NEAR(sum / N, 1.0, 3e-3);

        /* Albedo recovery at a few reflectances. */
        for (double rho = 0.18; rho < 0.95; rho += 0.36) {
            lam.rho = ls_spectrum_const(rho);
            Spectrum f;
            ls_bsdf_eval(&lam, wo, v3(0, 0, 1), &f);
            CHECK_NEAR((double)f.v[0], rho / LS_PI, 1e-6);
        }

        /* Helmholtz reciprocity: f(wo,wi) == f(wi,wo) for every BSDF. */
        Bsdf metal;
        memset(&metal, 0, sizeof metal);
        metal.kind = LS_BSDF_CONDUCTOR;
        metal.alpha = 0.25;
        ls_metal_aluminium(&metal.eta, &metal.kappa);
        for (int i = 0; i < 2000; ++i) {
            vec3 a = ls_sample_hemisphere_cosine(ls_rng_f(&rng), ls_rng_f(&rng));
            vec3 b = ls_sample_hemisphere_cosine(ls_rng_f(&rng), ls_rng_f(&rng));
            Spectrum fab, fba;
            ls_bsdf_eval(&metal, a, b, &fab);
            ls_bsdf_eval(&metal, b, a, &fba);
            CHECK_NEAR((double)fab.v[40], (double)fba.v[40], 1e-9);
        }

        /* sample() must report the same pdf that pdf() computes. */
        for (int i = 0; i < 2000; ++i) {
            vec3 wi;
            Spectrum f;
            ls_real pdf;
            if (!ls_bsdf_sample(&metal, wo, ls_rng_f(&rng), ls_rng_f(&rng),
                                &wi, &f, &pdf)) continue;
            CHECK_NEAR(pdf, ls_bsdf_pdf(&metal, wo, wi), 1e-9);
        }

        /* GGX directional albedo must never exceed 1.
         * NOTE: single-scattering GGX+Smith is NOT energy conserving at high
         * roughness -- the multiple-scattering deficit is real physics of the
         * model, not a bug. It is measured and reported here so that a genuine
         * GGX error cannot hide behind an expected shortfall. */
        for (double alpha = 0.1; alpha <= 0.91; alpha += 0.4) {
            metal.alpha = alpha;
            double e = 0.0;
            int n = 0;
            for (int i = 0; i < 200000; ++i) {
                vec3 wi;
                Spectrum f;
                ls_real pdf;
                if (!ls_bsdf_sample(&metal, wo, ls_rng_f(&rng), ls_rng_f(&rng),
                                    &wi, &f, &pdf)) { n++; continue; }
                e += (double)f.v[40] * wi.z / pdf;
                n++;
            }
            e /= n;
            CHECK(e <= 1.0 + 1e-6);
            NOTE("GGX alpha=%.2f directional albedo = %.4f (Fresnel-weighted)",
                 alpha, e);
        }
    }

    SECTION("MIS weights");
    {
        /* Partition of unity: the two weights for a pair of strategies must
         * sum to exactly 1, or energy is created or lost at every vertex. */
        for (int i = 0; i < 1000; ++i) {
            double a = ls_rng_f(&rng) * 10.0 + 1e-3;
            double b = ls_rng_f(&rng) * 10.0 + 1e-3;
            CHECK_NEAR(ls_mis_power2(a, b) + ls_mis_power2(b, a), 1.0, 1e-12);
        }
        /* A strategy that cannot generate the sample takes all the weight. */
        CHECK_NEAR(ls_mis_power2(1.0, 0.0), 1.0, 1e-12);
        CHECK_NEAR(ls_mis_power2(0.0, 1.0), 0.0, 1e-12);
    }
}
