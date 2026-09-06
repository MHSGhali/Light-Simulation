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
