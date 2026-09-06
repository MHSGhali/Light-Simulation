/* color.h — CIE colorimetry: spectrum -> XYZ -> sRGB.
 *
 * The colour matching functions are resampled onto the active spectral bin grid
 * once, on first use, and exposed as per-bin weight tables so that every
 * spectral integral in the engine goes through ls_spectrum_integrate_weighted().
 */
#ifndef LIGHTSIM_COLOR_H
#define LIGHTSIM_COLOR_H

#include "spectrum.h"

typedef struct { ls_real x, y, z; } XYZ;
typedef struct { ls_real r, g, b; } RGB;

/* Per-bin CMF weight tables. ybar IS V(lambda). */
const float *ls_cmf_xbar(void);
const float *ls_cmf_ybar(void);
const float *ls_cmf_zbar(void);

/* Integral of ybar over the band, in nm. Needed to normalise relative SPDs. */
ls_real ls_cmf_ybar_integral(void);

XYZ ls_spectrum_to_xyz(const Spectrum *s);
/* Chromaticity coordinates; independent of overall scale. */
void ls_xyz_chromaticity(XYZ c, ls_real *x, ls_real *y);

RGB ls_xyz_to_linear_srgb(XYZ c);
/* Apply the sRGB transfer function (gamma encode) to a linear value. */
ls_real ls_srgb_encode(ls_real linear);
RGB ls_rgb_gamma_encode(RGB c);

#endif /* LIGHTSIM_COLOR_H */
