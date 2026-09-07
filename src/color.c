#include "lightsim/color.h"
#include "lightsim/cie_data.h"
#include <math.h>
#include <pthread.h>

static float g_xbar[LS_NBINS], g_ybar[LS_NBINS], g_zbar[LS_NBINS];
static ls_real g_ybar_integral;
static pthread_once_t g_once = PTHREAD_ONCE_INIT;

static ls_real interp(const double *x, const double *y, int n, ls_real q) {
    if (q <= x[0]) return y[0];
    if (q >= x[n - 1]) return y[n - 1];
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (x[mid] <= q) lo = mid; else hi = mid;
    }
    ls_real t = (q - x[lo]) / (x[hi] - x[lo]);
    return y[lo] + t * (y[hi] - y[lo]);
}

static void build_tables(void) {
    ls_real sum_y = 0.0;
    for (int i = 0; i < LS_NBINS; ++i) {
        ls_real l = ls_bin_lambda(i);
        g_xbar[i] = (float)interp(ls_cie_lambda, ls_cie_xbar, ls_cie_count, l);
        g_ybar[i] = (float)interp(ls_cie_lambda, ls_cie_ybar, ls_cie_count, l);
        g_zbar[i] = (float)interp(ls_cie_lambda, ls_cie_zbar, ls_cie_count, l);
        sum_y += (ls_real)g_ybar[i];
    }
    g_ybar_integral = sum_y * LS_SPECTRAL_STEP;
}

const float *ls_cmf_xbar(void) { pthread_once(&g_once, build_tables); return g_xbar; }
const float *ls_cmf_ybar(void) { pthread_once(&g_once, build_tables); return g_ybar; }
const float *ls_cmf_zbar(void) { pthread_once(&g_once, build_tables); return g_zbar; }
ls_real ls_cmf_ybar_integral(void) { pthread_once(&g_once, build_tables); return g_ybar_integral; }

XYZ ls_spectrum_to_xyz(const Spectrum *s) {
    XYZ c;
    c.x = ls_spectrum_integrate_weighted(s, ls_cmf_xbar());
    c.y = ls_spectrum_integrate_weighted(s, ls_cmf_ybar());
    c.z = ls_spectrum_integrate_weighted(s, ls_cmf_zbar());
    return c;
}

void ls_xyz_chromaticity(XYZ c, ls_real *x, ls_real *y) {
    ls_real sum = c.x + c.y + c.z;
    if (sum == 0.0) { *x = 0.0; *y = 0.0; return; }
    *x = c.x / sum;
    *y = c.y / sum;
}

RGB ls_xyz_to_linear_srgb(XYZ c) {
    RGB o;
    o.r =  3.2404542 * c.x - 1.5371385 * c.y - 0.4985314 * c.z;
    o.g = -0.9692660 * c.x + 1.8760108 * c.y + 0.0415560 * c.z;
    o.b =  0.0556434 * c.x - 0.2040259 * c.y + 1.0572252 * c.z;
    return o;
}

ls_real ls_srgb_encode(ls_real u) {
    if (u <= 0.0031308) return 12.92 * u;
    return 1.055 * pow(u, 1.0 / 2.4) - 0.055;
}

RGB ls_rgb_gamma_encode(RGB c) {
    RGB o = { ls_srgb_encode(c.r), ls_srgb_encode(c.g), ls_srgb_encode(c.b) };
    return o;
}

RGB ls_rgb_from_spectrum_reflectance(const Spectrum *rho) {
    Spectrum d65 = ls_spectrum_daylight(6504.0);
    Spectrum seen = ls_spectrum_mul(*rho, d65);
    XYZ c = ls_spectrum_to_xyz(&seen);
    XYZ w = ls_spectrum_to_xyz(&d65);
    RGB zero = { 0.0, 0.0, 0.0 };
    if (w.y <= 0.0) return zero;
    XYZ norm = { c.x / w.y, c.y / w.y, c.z / w.y };
    return ls_xyz_to_linear_srgb(norm);
}

/* ---------------------------------------------- RGB -> reflectance ------- */

/* Smits, "An RGB-to-Spectrum Conversion for Reflectances" (1999).
 *
 * Seven basis curves tabulated at ten evenly spaced wavelengths over 380-720
 * nm. The case analysis below spends the colour from the least-saturated
 * component outward -- white first, then a secondary (cyan/magenta/yellow),
 * then a primary -- which is what keeps the result smooth and inside [0,1]
 * instead of reaching for spikes at the primaries.
 *
 * The band this engine samples (360-830 nm) is wider than Smits' table on both
 * sides; ls_spectrum_from_samples holds the endpoint value out to the edges,
 * which is the conventional and physically harmless choice for a reflectance
 * that is by definition flat-ish outside the visible range. */
#define SMITS_N 10

static const ls_real SMITS_L[SMITS_N] = {
    380.000000, 417.777778, 455.555556, 493.333333, 531.111111,
    568.888889, 606.666667, 644.444444, 682.222222, 720.000000
};
static const ls_real SMITS_WHITE[SMITS_N] = {
    1.0000, 1.0000, 0.9999, 0.9993, 0.9992, 0.9998, 1.0000, 1.0000, 1.0000, 1.0000
};
static const ls_real SMITS_CYAN[SMITS_N] = {
    0.9710, 0.9426, 1.0007, 1.0007, 1.0007, 1.0007, 0.1564, 0.0000, 0.0000, 0.0000
};
static const ls_real SMITS_MAGENTA[SMITS_N] = {
    1.0000, 1.0000, 0.9685, 0.2229, 0.0000, 0.0458, 0.8369, 1.0000, 1.0000, 0.9959
};
static const ls_real SMITS_YELLOW[SMITS_N] = {
    0.0001, 0.0000, 0.1088, 0.6651, 1.0000, 1.0000, 0.9996, 0.9586, 0.9685, 0.9840
};
static const ls_real SMITS_RED[SMITS_N] = {
    0.1012, 0.0515, 0.0000, 0.0000, 0.0000, 0.0000, 0.8325, 1.0149, 1.0149, 1.0149
};
static const ls_real SMITS_GREEN[SMITS_N] = {
    0.0000, 0.0000, 0.0273, 0.7937, 1.0000, 0.9418, 0.1719, 0.0000, 0.0000, 0.0025
};
static const ls_real SMITS_BLUE[SMITS_N] = {
    1.0000, 1.0000, 0.8916, 0.3323, 0.0000, 0.0000, 0.0003, 0.0369, 0.0483, 0.0496
};

static void smits_add(ls_real *acc, const ls_real *basis, ls_real w) {
    for (int i = 0; i < SMITS_N; ++i) acc[i] += w * basis[i];
}

Spectrum ls_spectrum_from_rgb_reflectance(RGB c) {
    /* Clamp on the way in: an out-of-gamut or over-unity colour must not be
     * able to produce a reflectance above 1. */
    ls_real r = ls_clamp(c.r, 0.0, 1.0);
    ls_real g = ls_clamp(c.g, 0.0, 1.0);
    ls_real b = ls_clamp(c.b, 0.0, 1.0);

    ls_real acc[SMITS_N] = { 0 };
    if (r <= g && r <= b) {
        smits_add(acc, SMITS_WHITE, r);
        if (g <= b) { smits_add(acc, SMITS_CYAN, g - r); smits_add(acc, SMITS_BLUE,  b - g); }
        else        { smits_add(acc, SMITS_CYAN, b - r); smits_add(acc, SMITS_GREEN, g - b); }
    } else if (g <= r && g <= b) {
        smits_add(acc, SMITS_WHITE, g);
        if (r <= b) { smits_add(acc, SMITS_MAGENTA, r - g); smits_add(acc, SMITS_BLUE, b - r); }
        else        { smits_add(acc, SMITS_MAGENTA, b - g); smits_add(acc, SMITS_RED,  r - b); }
    } else {
        smits_add(acc, SMITS_WHITE, b);
        if (r <= g) { smits_add(acc, SMITS_YELLOW, r - b); smits_add(acc, SMITS_GREEN, g - r); }
        else        { smits_add(acc, SMITS_YELLOW, g - b); smits_add(acc, SMITS_RED,   r - g); }
    }

    Spectrum s = ls_spectrum_from_samples(SMITS_L, acc, SMITS_N);
    /* The basis curves themselves exceed 1 slightly (cyan peaks at 1.0007, red
     * at 1.0149), so clamp per bin rather than trusting the construction. */
    for (int i = 0; i < LS_NBINS; ++i)
        s.v[i] = (float)ls_clamp((ls_real)s.v[i], 0.0, 1.0);

    /* Correct the LUMINANCE exactly, leaving only chromaticity error.
     *
     * Smits' basis reproduces a colour approximately, and the approximation is
     * worst in the saturated corners. But the quantity this simulator reports
     * is photometric, so of the two errors it is luminance that must not drift:
     * a wall imported at 60% grey has to reflect 60% of the light, whatever its
     * hue does. One scalar fixes it. */
    Spectrum d65 = ls_spectrum_daylight(6504.0);
    Spectrum seen = ls_spectrum_mul(s, d65);
    ls_real y_white = ls_spectrum_to_xyz(&d65).y;
    ls_real y_got = (y_white > 0.0) ? ls_spectrum_to_xyz(&seen).y / y_white : 0.0;
    ls_real y_want = 0.2126 * r + 0.7152 * g + 0.0722 * b;   /* sRGB luminance row */
    if (y_got > 1e-9 && y_want > 0.0) s = ls_spectrum_scale(s, y_want / y_got);

    /* Final cap, which is what makes the guarantee in color.h true.
     *
     * The bound is just below 1, not at it. A perfect reflector is an
     * idealisation no real surface reaches, and `Kd 1 1 1` in an MTL file is a
     * modelling default rather than a measurement -- but taken literally it
     * makes the equilibrium radiance of a closed room, Le/(1-rho), diverge. The
     * path tracer is depth-bounded so this cannot hang, it just returns a large
     * number that grows with max_depth. Capping here keeps an imported scene
     * inside the energy budget the furnace test verifies. */
    for (int i = 0; i < LS_NBINS; ++i)
        s.v[i] = (float)ls_clamp((ls_real)s.v[i], 0.0, LS_RHO_MAX);
    return s;
}
