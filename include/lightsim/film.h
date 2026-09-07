/* film.h — image buffers and writers.
 *
 * Film contents are ABSOLUTE spectral radiance, W/(m^2 sr nm). PFM dumps them
 * unmodified in physical units; only the PPM path applies exposure and the sRGB
 * transfer function. Keeping the buffer calibrated is what makes downstream
 * luminance work (glare, contrast) trustworthy.
 */
#ifndef LIGHTSIM_FILM_H
#define LIGHTSIM_FILM_H

#include "spectrum.h"
#include "color.h"

typedef struct {
    int          width, height;
    SpectrumAcc *pix;      /* running sum of spectral radiance */
    uint64_t    *n;        /* samples per pixel                */
} Film;

bool ls_film_init(Film *f, int width, int height);
void ls_film_free(Film *f);
void ls_film_add(Film *f, int x, int y, const Spectrum *L);
Spectrum ls_film_mean(const Film *f, int x, int y);

/* Tone-mapped sRGB. `exposure` scales linear luminance before encoding;
 * pass 0 to auto-expose so the 99th percentile maps to white. */
bool ls_film_write_ppm(const Film *f, const char *path, ls_real exposure);
/* Raw 32-bit float RGB in physical units, no tone mapping. */
bool ls_film_write_pfm(const Film *f, const char *path);

/* ---- false-colour mapping for scalar fields ---- */

/* Perceptually uniform viridis ramp; t is clamped to [0,1]. */
RGB ls_colormap_viridis(ls_real t);

/* Write a scalar field as a false-colour PPM, with `lo`..`hi` spanning the ramp.
 * `upscale` repeats each cell to make small grids legible. */
bool ls_write_falsecolor_ppm(const char *path, const ls_real *v, int nu, int nv,
                             ls_real lo, ls_real hi, int upscale);

#endif /* LIGHTSIM_FILM_H */
