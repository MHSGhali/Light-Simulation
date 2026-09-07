/* units.h — the ONLY radiometric -> photometric conversion in the engine.
 *
 * Design rule (do not violate): every internal quantity is stored as a
 * radiometric Spectrum. No struct anywhere holds a photometric value, and no
 * accumulator carries a unit tag. Photometric numbers exist only as the return
 * value of the functions below, which are called at report/export time.
 *
 * Tagging scalars with a unit enum was considered and rejected: it lets a
 * photometric value enter an accumulator that is later scaled by a radiometric
 * BSDF, and nothing catches it. Making the photometric value unreachable
 * outside this layer prevents that error class structurally.
 */
#ifndef LIGHTSIM_UNITS_H
#define LIGHTSIM_UNITS_H

#include "spectrum.h"

typedef enum { LS_UNITS_RADIOMETRIC, LS_UNITS_PHOTOMETRIC } LsUnitSystem;

/* The kind of quantity a Spectrum represents. Used only for labelling output;
 * it is never stored alongside the data. */
typedef enum {
    LS_Q_FLUX,        /* Phi : W        <-> lm       */
    LS_Q_IRRADIANCE,  /* E   : W/m^2    <-> lx       */
    LS_Q_INTENSITY,   /* I   : W/sr     <-> cd       */
    LS_Q_RADIANCE,    /* L   : W/m^2/sr <-> cd/m^2   */
    LS_Q_EXITANCE     /* M   : W/m^2    <-> lm/m^2   */
} LsQuantity;

/* Band integral of the spectrum: the radiometric scalar. */
ls_real ls_radiometric(const Spectrum *s);

/* Km * integral( S(lambda) * V(lambda) d lambda ), with V == CIE ybar. */
ls_real ls_photometric(const Spectrum *s);

/* Luminous efficacy of radiation, lm/W, relative to the power inside the
 * simulated band [LS_LAMBDA_MIN, LS_LAMBDA_MAX] only.
 *
 * CAUTION: the literature usually quotes LER against TOTAL radiant power across
 * all wavelengths. For a source with significant out-of-band emission (any
 * incandescent/blackbody source) the two differ enormously -- 121.6 vs 16.4
 * lm/W for a 2856 K blackbody. Use ls_luminous_efficacy_total() when the
 * out-of-band power is known. */
ls_real ls_luminous_efficacy_band(const Spectrum *s);

/* LER against a caller-supplied total radiant power (same units as the
 * spectrum's band integral). */
ls_real ls_luminous_efficacy_total(const Spectrum *s, ls_real total_radiant_power);

/* Radiant flux (W) equivalent to a luminous flux (lm) for a given spectral
 * shape. This is the one place a lumen enters the system; callers convert here
 * and hand watts to the light constructors, so no lumen is ever stored. */
ls_real ls_watts_from_lumens(ls_real lumens, const Spectrum *spd);

/* Report a quantity in the requested system. */
ls_real     ls_quantity_value(const Spectrum *s, LsUnitSystem sys);
const char *ls_quantity_unit(LsQuantity q, LsUnitSystem sys);
const char *ls_quantity_symbol(LsQuantity q, LsUnitSystem sys);

#endif /* LIGHTSIM_UNITS_H */
