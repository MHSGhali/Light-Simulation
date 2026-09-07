/* M1 validation: spectral integration, photometry, colorimetry.
 * Every expected value here is a closed form or a published CIE constant. */
#include "test.h"
#include "lightsim/spectrum.h"
#include "lightsim/units.h"
#include "lightsim/color.h"
#include "tests.h"

/* Peak of the blackbody LER curve, found by golden-section over T. */
static double peak_ler(double *t_at_peak) {
    double lo = 4000.0, hi = 9000.0;
    for (int it = 0; it < 200; ++it) {
        double m1 = lo + (hi - lo) / 3.0, m2 = hi - (hi - lo) / 3.0;
        Spectrum s1 = ls_spectrum_blackbody(m1), s2 = ls_spectrum_blackbody(m2);
        double e1 = ls_luminous_efficacy_total(&s1, ls_blackbody_total_radiance(m1));
        double e2 = ls_luminous_efficacy_total(&s2, ls_blackbody_total_radiance(m2));
        if (e1 < e2) lo = m1; else hi = m2;
    }
    *t_at_peak = 0.5 * (lo + hi);
    Spectrum s = ls_spectrum_blackbody(*t_at_peak);
    return ls_luminous_efficacy_total(&s, ls_blackbody_total_radiance(*t_at_peak));
}

void test_spectral(void) {
    SECTION("spectral integration");
    {
        /* A flat spectrum of 1.0 per nm over the band integrates to the band width. */
        Spectrum flat = ls_spectrum_const(1.0);
        double width = LS_LAMBDA_MAX - LS_LAMBDA_MIN + LS_SPECTRAL_STEP;
        CHECK_NEAR(ls_spectrum_integrate(&flat), width, 1e-12);

        /* normalize_to must make the band integral exact. */
        Spectrum n = ls_spectrum_normalize_to(flat, 12.5);
        CHECK_NEAR(ls_spectrum_integrate(&n), 12.5, 1e-6);

        /* A monochromatic spike carries exactly its stated power.
         * Tolerance floor is float32 storage (~1.2e-7 eps), not the maths:
         * power/step = 0.2 is not exactly representable. 1e-6 is the tightest
         * bound Spectrum's float backing can honour -- see the note in
         * test_spectral's photometry section. */
        Spectrum m = ls_spectrum_monochromatic(555.0, 1.0);
        CHECK_NEAR(ls_spectrum_integrate(&m), 1.0, 1e-6);
    }

    SECTION("photometry");
    {
        /* 555 nm is the definition point of the 683 lm/W constant, so the exact
         * identity is only assertable when 555 lands on the bin grid. It does at
         * 5 nm; it does NOT at 10 nm (the grid runs 550, 560), where the spike
         * snaps to 550 nm and V = 0.99495. That is inherent to binning, not a
         * bug, so the bound relaxes rather than the test disappearing. */
        const float *V = ls_cmf_ybar();
        int i555 = (int)(((555 - LS_LAMBDA_MIN_NM) + LS_SPECTRAL_STEP_NM / 2)
                         / LS_SPECTRAL_STEP_NM);
#if ((555 - LS_LAMBDA_MIN_NM) % LS_SPECTRAL_STEP_NM) == 0
        const double mono_tol = 1e-6;
        CHECK_NEAR(V[i555], 1.0, 1e-9);
#else
        const double mono_tol = 0.01;
        CHECK_NEAR(V[i555], 1.0, 0.01);
        NOTE("555 nm is off-grid at %d nm bins; V(peak bin) = %.5f",
             LS_SPECTRAL_STEP_NM, (double)V[i555]);
#endif

        /* THE defining photometric check: 1 W at 555 nm is exactly 683 lm.
         * Asserted at 1e-6, which is the float32 storage floor rather than a
         * physics tolerance -- the integral itself is done in double. Storing
         * Spectrum as double would buy another 9 digits of a precision that no
         * radiometric result needs, at 2x the inner-loop footprint. */
        Spectrum mono = ls_spectrum_monochromatic(555.0, 1.0);
        CHECK_NEAR(ls_photometric(&mono), 683.0, mono_tol);
        CHECK_NEAR(ls_radiometric(&mono), 1.0, 1e-6);
        CHECK_NEAR(ls_luminous_efficacy_band(&mono), 683.0, mono_tol);

        /* Integral of V(lambda) over the band. The CIE 1931 tables are built so
         * that the three CMF integrals are near-equal at ~106.857 nm; this pins
         * the table transcription and the integration rule together. */
        Spectrum unity = ls_spectrum_const(1.0);
        double ybar_int = ls_spectrum_integrate_weighted(&unity, V);
        CHECK_NEAR(ybar_int, 106.857, 1e-4);
        CHECK_NEAR(ls_cmf_ybar_integral(), 106.857, 1e-4);
        NOTE("integral of V(lambda) d(lambda) = %.4f nm", ybar_int);

        /* Equal-energy illuminant E: LER = Km * int(ybar) / band width.
         * Asserted as an identity (which checks that the efficacy function
         * composes correctly at any resolution) and, at the reference 5 nm
         * grid, against the literal value. The band width is NBINS*step, so
         * the literal shifts to 152.06 at 10 nm -- it is a property of the
         * quadrature grid, not a physical constant. */
        double band_nm = (double)LS_NBINS * LS_SPECTRAL_STEP;
        CHECK_NEAR(ls_luminous_efficacy_band(&unity),
                   LS_KM_LM_PER_W * ybar_int / band_nm, 1e-6);
        if (LS_SPECTRAL_STEP_NM == 5)
            CHECK_NEAR(ls_luminous_efficacy_band(&unity), 153.649, 1e-4);
        NOTE("illuminant E LER = %.3f lm/W over a %.0f nm band",
             ls_luminous_efficacy_band(&unity), band_nm);

        /* Far red carries power but almost no luminous flux. */
        Spectrum ir = ls_spectrum_monochromatic(800.0, 1.0);
        CHECK(ls_radiometric(&ir) > 0.99);
        CHECK(ls_photometric(&ir) < 0.01);
    }

    SECTION("blackbody luminous efficacy");
    {
        Spectrum a = ls_spectrum_blackbody(2856.0);   /* CIE illuminant A */
        double total = ls_blackbody_total_radiance(2856.0);

        /* Against TOTAL radiant power (all wavelengths). */
        double ler_total = ls_luminous_efficacy_total(&a, total);
        CHECK_NEAR(ler_total, 16.45, 0.01);

        /* Against visible-band power only. These differ by 7.4x; conflating them
         * is a live failure mode, so both are pinned. */
        /* Band-relative LER is a quadrature over the sharply peaked V(lambda),
         * so its accuracy tracks bin width directly: ~0.1% at 5 nm, ~1.1% at
         * 10 nm. 5 nm is the reference configuration. */
        double ler_band = ls_luminous_efficacy_band(&a);
        CHECK_NEAR(ler_band, 121.6, LS_SPECTRAL_STEP_NM <= 5 ? 0.01 : 0.02);
        NOTE("2856 K: LER_total = %.2f lm/W, LER_band = %.1f lm/W, ratio %.2f",
             ler_total, ler_band, ler_band / ler_total);

        /* The blackbody LER curve peaks near 6630 K at about 95.4 lm/W. */
        double t_peak = 0.0, ler_peak = peak_ler(&t_peak);
        CHECK_NEAR(ler_peak, 95.4, 0.01);
        CHECK_NEAR(t_peak, 6630.0, 0.02);
        NOTE("blackbody LER peaks at %.0f K -> %.2f lm/W", t_peak, ler_peak);
    }

    SECTION("colorimetry");
    {
        /* Equal-energy illuminant E sits at x = y = 1/3, but only to the extent
         * that the tabulated CMF integrals are equal -- do not tighten this. */
        Spectrum e = ls_spectrum_const(1.0);
        XYZ xyz_e = ls_spectrum_to_xyz(&e);
        double x, y;
        ls_xyz_chromaticity(xyz_e, &x, &y);
        CHECK_NEAR(x, 1.0 / 3.0, 1e-3);
        CHECK_NEAR(y, 1.0 / 3.0, 1e-3);
        NOTE("illuminant E chromaticity = (%.5f, %.5f)", x, y);

        /* D65 must land on its defining chromaticity. */
        Spectrum d65 = ls_spectrum_daylight(6504.0);
        XYZ xyz_d = ls_spectrum_to_xyz(&d65);
        ls_xyz_chromaticity(xyz_d, &x, &y);
        CHECK_NEAR(x, 0.31272, 1e-3);
        CHECK_NEAR(y, 0.32903, 1e-3);
        NOTE("D65 chromaticity = (%.5f, %.5f)", x, y);

        /* CIE illuminant A (2856 K Planckian) has a published chromaticity of
         * (0.44758, 0.40745) -- an independent check on the Planck evaluation,
         * the CMF tables and the integration all at once. */
        Spectrum ill_a = ls_spectrum_blackbody(2856.0);
        ls_xyz_chromaticity(ls_spectrum_to_xyz(&ill_a), &x, &y);
        CHECK_NEAR(x, 0.44758, 1e-3);
        CHECK_NEAR(y, 0.40745, 1e-3);
        NOTE("illuminant A chromaticity = (%.5f, %.5f)", x, y);

        /* D65 normalised to Y=1 is sRGB white: linear (1,1,1). */
        XYZ w = { xyz_d.x / xyz_d.y, 1.0, xyz_d.z / xyz_d.y };
        RGB rgb = ls_xyz_to_linear_srgb(w);
        CHECK_NEAR(rgb.r, 1.0, 3e-3);
        CHECK_NEAR(rgb.g, 1.0, 3e-3);
        CHECK_NEAR(rgb.b, 1.0, 3e-3);

        /* sRGB transfer function endpoints. */
        CHECK_NEAR(ls_srgb_encode(0.0), 0.0, 1e-12);
        CHECK_NEAR(ls_srgb_encode(1.0), 1.0, 1e-12);
    }

    SECTION("source models");
    {
        /* A Gaussian LED lobe carries exactly the flux it was given. */
        Spectrum led = ls_spectrum_gaussian(520.0, 30.0, 2.5);
        CHECK_NEAR(ls_spectrum_integrate(&led), 2.5, 1e-6);

        /* A green LED is far more luminous per watt than a blue one. */
        Spectrum green = ls_spectrum_gaussian(530.0, 30.0, 1.0);
        Spectrum blue  = ls_spectrum_gaussian(465.0, 25.0, 1.0);
        CHECK(ls_luminous_efficacy_band(&green) > 3.0 * ls_luminous_efficacy_band(&blue));
        NOTE("LED efficacy: green(530nm) %.0f lm/W, blue(465nm) %.0f lm/W",
             ls_luminous_efficacy_band(&green), ls_luminous_efficacy_band(&blue));

        /* Hotter blackbody puts a larger fraction of its power in the visible. */
        Spectrum b3000 = ls_spectrum_blackbody(3000.0);
        Spectrum b6000 = ls_spectrum_blackbody(6000.0);
        double f3000 = ls_radiometric(&b3000) / ls_blackbody_total_radiance(3000.0);
        double f6000 = ls_radiometric(&b6000) / ls_blackbody_total_radiance(6000.0);
        CHECK(f6000 > f3000);
        NOTE("visible fraction: 3000 K %.1f%%, 6000 K %.1f%%", f3000 * 100, f6000 * 100);
    }
}
