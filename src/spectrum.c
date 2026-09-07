#include "lightsim/spectrum.h"
#include "lightsim/cie_data.h"
#include <math.h>
#include <string.h>

Spectrum ls_spectrum_zero(void) {
    Spectrum s;
    memset(s.v, 0, sizeof s.v);
    return s;
}

Spectrum ls_spectrum_const(ls_real value) {
    Spectrum s;
    for (int i = 0; i < LS_NBINS; ++i) s.v[i] = (float)value;
    return s;
}

Spectrum ls_spectrum_monochromatic(ls_real lambda_nm, ls_real power) {
    /* Stored as power/step so that the band integral is exactly `power`:
     *   integral = (power/step) * step = power. */
    Spectrum s = ls_spectrum_zero();
    if (lambda_nm < LS_LAMBDA_MIN || lambda_nm > LS_LAMBDA_MAX) return s;
    int i = (int)floor((lambda_nm - LS_LAMBDA_MIN) / LS_SPECTRAL_STEP + 0.5);
    if (i < 0) i = 0;
    if (i >= LS_NBINS) i = LS_NBINS - 1;
    s.v[i] = (float)(power / LS_SPECTRAL_STEP);
    return s;
}

Spectrum ls_spectrum_blackbody(ls_real temperature_k) {
    /* Planck spectral radiance, W/(m^2 sr m). The 1e-9 converts per-metre to
     * per-nanometre -- the ONE place this conversion is applied. */
    Spectrum s = ls_spectrum_zero();
    if (temperature_k <= 0.0) return s;
    const ls_real c1 = 2.0 * LS_PLANCK_H * LS_LIGHT_C * LS_LIGHT_C;
    const ls_real c2 = LS_PLANCK_H * LS_LIGHT_C / LS_BOLTZMANN_K;
    for (int i = 0; i < LS_NBINS; ++i) {
        ls_real l = ls_bin_lambda(i) * 1e-9;               /* metres */
        ls_real l5 = l * l * l * l * l;
        ls_real e = exp(c2 / (l * temperature_k)) - 1.0;
        s.v[i] = (float)((c1 / (l5 * e)) * 1e-9);          /* -> per nm */
    }
    return s;
}

ls_real ls_blackbody_total_radiance(ls_real temperature_k) {
    ls_real t2 = temperature_k * temperature_k;
    return LS_STEFAN_BOLTZMANN * t2 * t2 / LS_PI;
}

static ls_real interp_table(const double *x, const double *y, int n, ls_real q) {
    /* Linear interpolation, clamped to the endpoint values outside the table. */
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

Spectrum ls_spectrum_daylight(ls_real cct_k) {
    /* CIE method: chromaticity from CCT, then S = S0 + M1*S1 + M2*S2. */
    ls_real t = cct_k, xd;
    if (t < 4000.0) t = 4000.0;
    if (t > 25000.0) t = 25000.0;
    if (t <= 7000.0)
        xd = -4.6070e9 / (t*t*t) + 2.9678e6 / (t*t) + 0.09911e3 / t + 0.244063;
    else
        xd = -2.0064e9 / (t*t*t) + 1.9018e6 / (t*t) + 0.24748e3 / t + 0.237040;
    ls_real yd = -3.000 * xd * xd + 2.870 * xd - 0.275;

    ls_real denom = 0.0241 + 0.2562 * xd - 0.7341 * yd;
    ls_real m1 = (-1.3515 -  1.7703 * xd +  5.9114 * yd) / denom;
    ls_real m2 = ( 0.0300 - 31.4424 * xd + 30.0717 * yd) / denom;

    Spectrum s = ls_spectrum_zero();
    for (int i = 0; i < LS_NBINS; ++i) {
        ls_real l = ls_bin_lambda(i);
        ls_real s0 = interp_table(ls_cie_daylight_lambda, ls_cie_s0, ls_cie_daylight_count, l);
        ls_real s1 = interp_table(ls_cie_daylight_lambda, ls_cie_s1, ls_cie_daylight_count, l);
        ls_real s2 = interp_table(ls_cie_daylight_lambda, ls_cie_s2, ls_cie_daylight_count, l);
        s.v[i] = (float)(s0 + m1 * s1 + m2 * s2);
    }
    return s;
}

Spectrum ls_spectrum_gaussian(ls_real center_nm, ls_real fwhm_nm, ls_real power) {
    Spectrum s = ls_spectrum_zero();
    if (fwhm_nm <= 0.0) return s;
    ls_real sigma = fwhm_nm / 2.354820045030949; /* 2*sqrt(2*ln2) */
    for (int i = 0; i < LS_NBINS; ++i) {
        ls_real d = (ls_bin_lambda(i) - center_nm) / sigma;
        s.v[i] = (float)exp(-0.5 * d * d);
    }
    return ls_spectrum_normalize_to(s, power);
}

Spectrum ls_spectrum_from_samples(const ls_real *lambda_nm, const ls_real *value, int n) {
    Spectrum s = ls_spectrum_zero();
    if (n <= 0) return s;
    for (int i = 0; i < LS_NBINS; ++i)
        s.v[i] = (float)interp_table(lambda_nm, value, n, ls_bin_lambda(i));
    return s;
}

Spectrum ls_spectrum_add(Spectrum a, Spectrum b) {
    for (int i = 0; i < LS_NBINS; ++i) a.v[i] += b.v[i];
    return a;
}
Spectrum ls_spectrum_sub(Spectrum a, Spectrum b) {
    for (int i = 0; i < LS_NBINS; ++i) a.v[i] -= b.v[i];
    return a;
}
Spectrum ls_spectrum_mul(Spectrum a, Spectrum b) {
    for (int i = 0; i < LS_NBINS; ++i) a.v[i] *= b.v[i];
    return a;
}
Spectrum ls_spectrum_scale(Spectrum a, ls_real s) {
    for (int i = 0; i < LS_NBINS; ++i) a.v[i] = (float)((ls_real)a.v[i] * s);
    return a;
}
void ls_spectrum_add_inplace(Spectrum *a, const Spectrum *b) {
    for (int i = 0; i < LS_NBINS; ++i) a->v[i] += b->v[i];
}
bool ls_spectrum_is_black(const Spectrum *a) {
    for (int i = 0; i < LS_NBINS; ++i) if (a->v[i] != 0.0f) return false;
    return true;
}
ls_real ls_spectrum_max(const Spectrum *a) {
    /* Reductions accumulate in double deliberately: the bins are float, but a
     * 95-term sum in float would lose several digits before integration. */
    ls_real m = (ls_real)a->v[0];
    for (int i = 1; i < LS_NBINS; ++i) {
        ls_real x = (ls_real)a->v[i];
        if (x > m) m = x;
    }
    return m;
}
ls_real ls_spectrum_mean(const Spectrum *a) {
    ls_real sum = 0.0;
    for (int i = 0; i < LS_NBINS; ++i) sum += (ls_real)a->v[i];
    return sum / (ls_real)LS_NBINS;
}

SpectrumAcc ls_acc_zero(void) {
    SpectrumAcc a;
    memset(a.v, 0, sizeof a.v);
    return a;
}
void ls_acc_add_scaled(SpectrumAcc *acc, const Spectrum *s, ls_real w) {
    for (int i = 0; i < LS_NBINS; ++i) acc->v[i] += (ls_real)s->v[i] * w;
}
Spectrum ls_acc_mean(const SpectrumAcc *acc, uint64_t n) {
    Spectrum s = ls_spectrum_zero();
    if (n == 0) return s;
    ls_real inv = 1.0 / (ls_real)n;
    for (int i = 0; i < LS_NBINS; ++i) s.v[i] = (float)(acc->v[i] * inv);
    return s;
}

/* Rectangle (bin-centred) rule: each sample represents a bin of width
 * LS_SPECTRAL_STEP centred on its wavelength.
 *
 * A trapezoid rule was considered and REJECTED. It weights the first and last
 * samples by 1/2, so a monochromatic spike placed in the edge bin would
 * integrate to half its power -- energy silently vanishing at the band edges,
 * with no way for a caller to notice. The rectangle rule makes
 * ls_spectrum_monochromatic() carry exactly its stated power in EVERY bin,
 * which is the invariant the 683 lm/W test depends on. */
ls_real ls_spectrum_integrate(const Spectrum *a) {
    ls_real sum = 0.0;
    for (int i = 0; i < LS_NBINS; ++i) sum += (ls_real)a->v[i];
    return sum * LS_SPECTRAL_STEP;   /* bin width in NANOMETRES */
}

ls_real ls_spectrum_integrate_weighted(const Spectrum *a, const float *weight) {
    ls_real sum = 0.0;
    for (int i = 0; i < LS_NBINS; ++i) sum += (ls_real)a->v[i] * (ls_real)weight[i];
    return sum * LS_SPECTRAL_STEP;
}

Spectrum ls_spectrum_normalize_to(Spectrum a, ls_real target) {
    ls_real cur = ls_spectrum_integrate(&a);
    if (cur == 0.0) return a;
    return ls_spectrum_scale(a, target / cur);
}
