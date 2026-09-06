/* lightsim — CLI driver.
 *
 * M1 scope: the `source` subcommand, which characterises a light source's
 * spectrum and reports it in either unit system. Scene rendering arrives with
 * the geometry and transport milestones.
 */
#include "lightsim/spectrum.h"
#include "lightsim/units.h"
#include "lightsim/color.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void usage(const char *prog) {
    printf("usage: %s source <model> [args] [--units radiometric|photometric]\n\n", prog);
    printf("models:\n");
    printf("  blackbody <T_kelvin>              Planckian radiator\n");
    printf("  daylight  <CCT_kelvin>            CIE D-series illuminant\n");
    printf("  led       <center_nm> <fwhm_nm>   Gaussian LED lobe\n");
    printf("  mono      <lambda_nm>             monochromatic spike\n\n");
    printf("Sources are normalised to 1 W of radiant flux within the %d-%d nm band.\n",
           LS_LAMBDA_MIN_NM, LS_LAMBDA_MAX_NM);
}

static int cmd_source(int argc, char **argv) {
    LsUnitSystem sys = LS_UNITS_RADIOMETRIC;
    const char *model = NULL;
    double a[2] = { 0.0, 0.0 };
    int na = 0;

    for (int i = 0; i < argc; ++i) {
        if (strcmp(argv[i], "--units") == 0 && i + 1 < argc) {
            ++i;
            if (strcmp(argv[i], "photometric") == 0) sys = LS_UNITS_PHOTOMETRIC;
            else if (strcmp(argv[i], "radiometric") == 0) sys = LS_UNITS_RADIOMETRIC;
            else { fprintf(stderr, "unknown unit system: %s\n", argv[i]); return 1; }
        } else if (!model) {
            model = argv[i];
        } else if (na < 2) {
            a[na++] = atof(argv[i]);
        }
    }
    if (!model) { fprintf(stderr, "no source model given\n"); return 1; }

    Spectrum s;
    double total_radiant = 0.0;   /* including out-of-band power, when known */
    if (strcmp(model, "blackbody") == 0) {
        if (na < 1) { fprintf(stderr, "blackbody needs a temperature\n"); return 1; }
        s = ls_spectrum_blackbody(a[0]);
        double band = ls_radiometric(&s);
        total_radiant = ls_blackbody_total_radiance(a[0]);
        /* Renormalise the band to 1 W, carrying the out-of-band ratio with it. */
        if (band > 0.0) { total_radiant /= band; s = ls_spectrum_normalize_to(s, 1.0); }
    } else if (strcmp(model, "daylight") == 0) {
        if (na < 1) { fprintf(stderr, "daylight needs a CCT\n"); return 1; }
        s = ls_spectrum_normalize_to(ls_spectrum_daylight(a[0]), 1.0);
    } else if (strcmp(model, "led") == 0) {
        if (na < 2) { fprintf(stderr, "led needs center_nm and fwhm_nm\n"); return 1; }
        s = ls_spectrum_gaussian(a[0], a[1], 1.0);
    } else if (strcmp(model, "mono") == 0) {
        if (na < 1) { fprintf(stderr, "mono needs a wavelength\n"); return 1; }
        s = ls_spectrum_monochromatic(a[0], 1.0);
    } else {
        fprintf(stderr, "unknown model: %s\n", model);
        return 1;
    }

    XYZ xyz = ls_spectrum_to_xyz(&s);
    double cx, cy;
    ls_xyz_chromaticity(xyz, &cx, &cy);

    printf("source: %s", model);
    for (int i = 0; i < na; ++i) printf(" %g", a[i]);
    printf("\n  band            %d-%d nm, %d bins at %d nm\n",
           LS_LAMBDA_MIN_NM, LS_LAMBDA_MAX_NM, LS_NBINS, LS_SPECTRAL_STEP_NM);
    printf("  %-6s %-8s %12.4f %s\n",
           ls_quantity_symbol(LS_Q_FLUX, sys), "(flux)",
           ls_quantity_value(&s, sys), ls_quantity_unit(LS_Q_FLUX, sys));
    printf("  efficacy (band) %12.2f lm/W\n", ls_luminous_efficacy_band(&s));
    if (total_radiant > 0.0)
        printf("  efficacy (total)%12.2f lm/W   [out-of-band power included]\n",
               ls_luminous_efficacy_total(&s, total_radiant));
    printf("  chromaticity     (%.5f, %.5f)\n", cx, cy);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) { usage(argv[0]); return 1; }
    if (strcmp(argv[1], "source") == 0) return cmd_source(argc - 2, argv + 2);
    usage(argv[0]);
    return 1;
}
