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

/* A reflectance shown as a display colour.
 *
 * A reflectance is not an emitter: it has no colour until something shines on
 * it. Pushing rho straight through XYZ implicitly views it under an
 * equal-energy illuminant, and since sRGB's white point is D65 that renders a
 * neutral grey as noticeably warm. Viewing it under D65 and normalising by that
 * illuminant's own luminance is the correct reduction, and it sends a flat rho
 * of 0.8 to exactly (0.8, 0.8, 0.8). Not clamped -- a saturated reflectance can
 * fall outside the sRGB gamut and the caller decides what to do about it. */
RGB ls_rgb_from_spectrum_reflectance(const Spectrum *rho);

/* Ceiling on an imported reflectance. Just below 1 rather than at it: a perfect
 * reflector is an idealisation, `Kd 1 1 1` in an MTL file is a modelling
 * default rather than a measurement, and taken literally it makes a closed
 * room's equilibrium radiance Le/(1-rho) diverge. */
#define LS_RHO_MAX 0.99

/* The other direction, for REFLECTANCE only: a plausible spectrum for a linear
 * sRGB colour. Needed because imported geometry (Blender, OBJ) carries nothing
 * but RGB, while every material in this engine is spectral.
 *
 * This is not an inverse -- infinitely many spectra share a colour. It returns
 * the SMOOTHEST one, which is the right choice here: a spectrum with narrow
 * features would beat against a narrow-band source and make the result depend
 * on the bin grid rather than on the physics.
 *
 * The result is bounded to [0,1] per bin, so an imported albedo can never
 * create energy and break the furnace test. Nothing authored in a .scene goes
 * through this path; those materials are spectral from the start. */
Spectrum ls_spectrum_from_rgb_reflectance(RGB c);

#endif /* LIGHTSIM_COLOR_H */
