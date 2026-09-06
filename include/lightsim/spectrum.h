/* spectrum.h — fixed-size sampled spectral distribution.
 *
 * UNIT CONVENTION (read this before touching any spectral code):
 *
 *   Every Spectrum holds a SPECTRAL DENSITY PER NANOMETRE.
 *   So an irradiance spectrum is W/(m^2 nm), a radiance spectrum W/(m^2 sr nm),
 *   a flux spectrum W/nm, and so on.  The band-integrated quantity is recovered
 *   by ls_spectrum_integrate(), which multiplies by the bin width in NANOMETRES.
 *
 *   Choosing per-nm (rather than per-metre) keeps values human-scaled and avoids
 *   the 1e9 factor that silently wrecks radiometric integrals.  Planck's law is
 *   expressed per-metre by physics convention, so ls_spectrum_blackbody()
 *   converts on the way in -- that is the ONE place the conversion happens.
 *
 * REPRESENTATION
 *   Fixed-size POD (float v[LS_NBINS]).  Never a pointer+count: this keeps
 *   Spectrum stack-allocatable, forbids allocation in the path-tracing inner
 *   loop, and lets the compiler vectorise the elementwise ops.
 *
 * DISPERSION
 *   Carrying all bins along one geometric path is only valid because every
 *   material in this engine is NON-DISPERSIVE: wavelength changes the value of
 *   the BSDF but never the sampled direction.  If a dispersive dielectric
 *   (refraction through glass, thin films) is ever added, this whole approach
 *   must be replaced by hero-wavelength sampling.  Do not add one casually.
 */
#ifndef LIGHTSIM_SPECTRUM_H
#define LIGHTSIM_SPECTRUM_H

#include "core.h"

/* The band edges and step are INTEGER nanometres so that LS_NBINS is a true
 * integer constant expression -- float arithmetic here makes `float v[LS_NBINS]`
 * a folded VLA rather than a fixed-size array, which is a compiler extension.
 * The step must divide the band evenly; 5 and 10 nm both do. */
#define LS_LAMBDA_MIN_NM 360
#define LS_LAMBDA_MAX_NM 830

#ifndef LS_SPECTRAL_STEP_NM
#define LS_SPECTRAL_STEP_NM 5   /* override with -DLS_SPECTRAL_STEP_NM=10 */
#endif

#define LS_NBINS ((LS_LAMBDA_MAX_NM - LS_LAMBDA_MIN_NM) / LS_SPECTRAL_STEP_NM + 1)

/* Real-valued aliases for use inside the maths. */
#define LS_LAMBDA_MIN    ((ls_real)LS_LAMBDA_MIN_NM)
#define LS_LAMBDA_MAX    ((ls_real)LS_LAMBDA_MAX_NM)
#define LS_SPECTRAL_STEP ((ls_real)LS_SPECTRAL_STEP_NM)

typedef struct { float v[LS_NBINS]; } Spectrum;

/* Centre wavelength of bin i, in nanometres. */
static inline ls_real ls_bin_lambda(int i) {
    return LS_LAMBDA_MIN + (ls_real)i * LS_SPECTRAL_STEP;
}

/* ---- construction ---- */
Spectrum ls_spectrum_zero(void);
Spectrum ls_spectrum_const(ls_real value);
/* Single-bin spike whose band integral equals `power`. Used for monochromatic
 * sources; the value stored is power/step so the integral is exactly `power`. */
Spectrum ls_spectrum_monochromatic(ls_real lambda_nm, ls_real power);
/* Planck spectral radiance of a blackbody at T, converted to per-nm. */
Spectrum ls_spectrum_blackbody(ls_real temperature_k);
/* Total radiance of a blackbody integrated over ALL wavelengths, W/(m^2 sr).
 * Stefan-Boltzmann: sigma T^4 / pi. Needed because a Spectrum is truncated to
 * the visible band and therefore cannot supply the out-of-band power. */
ls_real ls_blackbody_total_radiance(ls_real temperature_k);
/* CIE D-series daylight illuminant from the S0/S1/S2 basis (relative units). */
Spectrum ls_spectrum_daylight(ls_real cct_k);
/* Gaussian lobe, a serviceable model of a single-die LED. `fwhm_nm` is the
 * full width at half maximum; the band integral equals `power`. */
Spectrum ls_spectrum_gaussian(ls_real center_nm, ls_real fwhm_nm, ls_real power);
/* Resample an arbitrary tabulated SPD (sorted by wavelength) onto the bin grid. */
Spectrum ls_spectrum_from_samples(const ls_real *lambda_nm, const ls_real *value, int n);

/* ---- elementwise ops ---- */
Spectrum ls_spectrum_add(Spectrum a, Spectrum b);
Spectrum ls_spectrum_sub(Spectrum a, Spectrum b);
Spectrum ls_spectrum_mul(Spectrum a, Spectrum b);
Spectrum ls_spectrum_scale(Spectrum a, ls_real s);
void     ls_spectrum_add_inplace(Spectrum *a, const Spectrum *b);
bool     ls_spectrum_is_black(const Spectrum *a);
ls_real  ls_spectrum_max(const Spectrum *a);
ls_real  ls_spectrum_mean(const Spectrum *a);

/* ---- integration ----
 * THE single place bin width is applied. Every spectral integral in the engine
 * must go through one of these two functions; never hand-roll the loop. */
ls_real ls_spectrum_integrate(const Spectrum *a);
/* Integral of a against a per-bin weight table (e.g. the CIE ybar / V(lambda)). */
ls_real ls_spectrum_integrate_weighted(const Spectrum *a, const float *weight);

/* Rescale so that the band integral equals `target`. No-op on a black spectrum. */
Spectrum ls_spectrum_normalize_to(Spectrum a, ls_real target);

#endif /* LIGHTSIM_SPECTRUM_H */
