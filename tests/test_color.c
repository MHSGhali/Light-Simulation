/* Colour reduction, and its approximate inverse.
 *
 * Imported geometry (Blender, OBJ) carries only RGB, while every material in
 * this engine is spectral. ls_spectrum_from_rgb_reflectance bridges that gap,
 * and the thing worth guarding is that it is a real inverse of the reduction
 * the exporter already uses -- not merely "some spectrum of about that colour"
 * -- and that it can never hand the transport core an albedo above 1. */
#include "test.h"
#include "tests.h"
#include "lightsim/color.h"
#include <string.h>

/* The forward direction, as the exporter and the viewer use it. */
static RGB round_trip(RGB in) {
    Spectrum s = ls_spectrum_from_rgb_reflectance(in);
    return ls_rgb_from_spectrum_reflectance(&s);
}

static double bin_max(const Spectrum *s) {
    double m = 0.0;
    for (int i = 0; i < LS_NBINS; ++i) if (s->v[i] > m) m = s->v[i];
    return m;
}
static double bin_min(const Spectrum *s) {
    double m = s->v[0];
    for (int i = 1; i < LS_NBINS; ++i) if (s->v[i] < m) m = s->v[i];
    return m;
}

/* Bins are float; the cap is applied in double. The nearest float to 0.99 sits
 * just above it, so the bound is asserted with one float ulp of slack. */
#define RHO_CEIL (LS_RHO_MAX + 1e-6)

static double luminance(RGB c) { return 0.2126*c.r + 0.7152*c.g + 0.0722*c.b; }

void test_color(void) {
    SECTION("a grey uplifts to a flat spectrum of the same albedo");
    {
        /* The one case with an exact answer: a neutral colour has no hue to
         * reconstruct, so the smoothest spectrum reproducing it is the constant
         * one, and it must come back through the reduction unchanged. */
        const double grey[] = { 0.0, 0.18, 0.5, 0.8, 0.9 };
        for (unsigned k = 0; k < sizeof grey / sizeof *grey; ++k) {
            RGB in = { grey[k], grey[k], grey[k] };
            Spectrum s = ls_spectrum_from_rgb_reflectance(in);
            CHECK_NEAR(bin_max(&s) - bin_min(&s), 0.0, 1e-3);   /* flat */
            RGB out = round_trip(in);
            /* Luminance is exact to the precision of the representation --
             * that is what the scalar correction in the uplift buys, and it is
             * the quantity this simulator reports. The floor is Spectrum's
             * float bins at ~1.2e-7 relative, not the method. */
            CHECK_NEAR(luminance(out), grey[k], 1e-6);
            /* The channels carry a little more: ls_spectrum_daylight(6504) is a
             * very good D65 but not exactly the sRGB white point, so a neutral
             * does not land on a perfectly equal triple. */
            CHECK_NEAR(out.r, grey[k], 2e-3);
            CHECK_NEAR(out.g, grey[k], 2e-3);
            CHECK_NEAR(out.b, grey[k], 2e-3);
        }

        /* A perfect reflector is deliberately not honoured: it is capped, so a
         * closed imported room cannot have a divergent equilibrium. */
        RGB white = round_trip((RGB){ 1.0, 1.0, 1.0 });
        CHECK_NEAR(luminance(white), LS_RHO_MAX, 2e-3);
        NOTE("Kd 1 1 1 imports as reflectance %.4f, not 1.0", luminance(white));
    }

    SECTION("no imported colour can create energy");
    {
        /* The furnace test's guarantee reaches only as far as the albedos it is
         * given. This is what keeps an imported material inside it. */
        double worst = 0.0;
        for (int r = 0; r <= 4; ++r)
        for (int g = 0; g <= 4; ++g)
        for (int b = 0; b <= 4; ++b) {
            RGB in = { r / 4.0, g / 4.0, b / 4.0 };
            Spectrum s = ls_spectrum_from_rgb_reflectance(in);
            CHECK(bin_min(&s) >= 0.0);
            CHECK(bin_max(&s) <= RHO_CEIL);
            if (bin_max(&s) > worst) worst = bin_max(&s);
        }
        NOTE("highest reflectance bin produced over the RGB cube: %.6f", worst);

        /* Out-of-gamut and over-unity inputs are clamped, not trusted. */
        RGB hot = { 4.0, -1.0, 0.5 };
        Spectrum s = ls_spectrum_from_rgb_reflectance(hot);
        CHECK(bin_min(&s) >= 0.0);
        CHECK(bin_max(&s) <= RHO_CEIL);
    }

    SECTION("colour survives the round trip");
    {
        /* Smits' basis is an approximation, so this pins the achieved error
         * rather than asserting an exact inverse. A regression that made the
         * uplift meaningfully worse would show up here as a tolerance failure. */
        double worst = 0.0;
        const char *worst_at = "";
        static char buf[64];
        for (int r = 0; r <= 4; ++r)
        for (int g = 0; g <= 4; ++g)
        for (int b = 0; b <= 4; ++b) {
            RGB in = { r / 4.0, g / 4.0, b / 4.0 };
            RGB out = round_trip(in);
            double e = fabs(out.r - in.r);
            if (fabs(out.g - in.g) > e) e = fabs(out.g - in.g);
            if (fabs(out.b - in.b) > e) e = fabs(out.b - in.b);
            if (e > worst) {
                worst = e;
                snprintf(buf, sizeof buf, "(%.2f %.2f %.2f)", in.r, in.g, in.b);
                worst_at = buf;
            }
        }
        NOTE("worst round-trip error %.4f at %s", worst, worst_at);
        CHECK(worst < 0.10);
    }

    SECTION("hue is preserved, not just luminance");
    {
        /* The failure mode a luminance-only check would miss: an uplift that
         * returned grey for everything would pass an energy test and be
         * useless. Red must read as red. */
        RGB red   = round_trip((RGB){ 0.9, 0.1, 0.1 });
        RGB green = round_trip((RGB){ 0.1, 0.9, 0.1 });
        RGB blue  = round_trip((RGB){ 0.1, 0.1, 0.9 });
        CHECK(red.r   > red.g   && red.r   > red.b);
        CHECK(green.g > green.r && green.g > green.b);
        CHECK(blue.b  > blue.r  && blue.b  > blue.g);
    }
}
