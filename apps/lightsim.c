/* lightsim — CLI driver.
 *
 *   lightsim source <model> [args]        characterise a source spectrum
 *   lightsim grid   <scene.txt>           irradiance/illuminance over a grid
 *   lightsim render <scene.txt>           camera image
 */
#include "lightsim/scenefile.h"
#include "lightsim/integrator.h"
#include "lightsim/analysis.h"
#include "lightsim/film.h"
#include "lightsim/units.h"
#include "lightsim/thread.h"
#include "lightsim/export.h"
#include "lightsim/import.h"
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static void usage(const char *prog) {
    printf("usage:\n");
    printf("  %s source <model> [args] [--units radiometric|photometric]\n", prog);
    printf("  %s grid   <scene.txt> [--out PREFIX] [--units U] [--spp N]\n", prog);
    printf("           [--indirect N] [--depth N] [--threads N] [--blender FILE.py]\n");
    printf("           [--import FILE.obj]   add a mesh from Blender/OBJ\n");
    printf("  %s render <scene.txt> [--out FILE.ppm] [--spp N] [--depth N]\n\n", prog);
    printf("source models:\n");
    printf("  blackbody <T_K>                 Planckian radiator\n");
    printf("  daylight  <CCT_K>               CIE D-series illuminant\n");
    printf("  led       <center_nm> <fwhm_nm> Gaussian LED lobe\n");
    printf("  mono      <lambda_nm>           monochromatic spike\n");
}

/* ---------------------------------------------------------------- options */
typedef struct {
    LsUnitSystem units;
    const char  *out;
    int spp, indirect, depth, threads;
    const char *blender;
    const char *import;
} Opts;

static Opts parse_opts(int argc, char **argv, int from, const char *def_out) {
    Opts o = { LS_UNITS_PHOTOMETRIC, def_out, 256, 64, 4, 0, NULL, NULL };
    for (int i = from; i < argc; ++i) {
        if (!strcmp(argv[i], "--units") && i + 1 < argc) {
            ++i;
            o.units = !strcmp(argv[i], "radiometric") ? LS_UNITS_RADIOMETRIC
                                                      : LS_UNITS_PHOTOMETRIC;
        } else if (!strcmp(argv[i], "--out")      && i + 1 < argc) o.out = argv[++i];
        else if   (!strcmp(argv[i], "--spp")      && i + 1 < argc) o.spp = atoi(argv[++i]);
        else if   (!strcmp(argv[i], "--indirect") && i + 1 < argc) o.indirect = atoi(argv[++i]);
        else if   (!strcmp(argv[i], "--depth")    && i + 1 < argc) o.depth = atoi(argv[++i]);
        else if   (!strcmp(argv[i], "--threads")  && i + 1 < argc) o.threads = atoi(argv[++i]);
        else if   (!strcmp(argv[i], "--blender")  && i + 1 < argc) o.blender = argv[++i];
        else if   (!strcmp(argv[i], "--import")   && i + 1 < argc) o.import = argv[++i];
    }
    return o;
}

/* ----------------------------------------------------------------- source */
static int cmd_source(int argc, char **argv) {
    Opts o = parse_opts(argc, argv, 0, NULL);
    const char *model = NULL;
    double a[2] = { 0, 0 };
    int na = 0;
    for (int i = 0; i < argc; ++i) {
        if (!strcmp(argv[i], "--units")) { ++i; continue; }
        if (!model) model = argv[i];
        else if (na < 2) a[na++] = atof(argv[i]);
    }
    if (!model) { fprintf(stderr, "no source model given\n"); return 1; }

    Spectrum s;
    double total = 0.0;
    if (!strcmp(model, "blackbody")) {
        s = ls_spectrum_blackbody(a[0]);
        double band = ls_radiometric(&s);
        total = ls_blackbody_total_radiance(a[0]);
        if (band > 0) { total /= band; s = ls_spectrum_normalize_to(s, 1.0); }
    } else if (!strcmp(model, "daylight")) {
        s = ls_spectrum_normalize_to(ls_spectrum_daylight(a[0]), 1.0);
    } else if (!strcmp(model, "led")) {
        s = ls_spectrum_gaussian(a[0], a[1], 1.0);
    } else if (!strcmp(model, "mono")) {
        s = ls_spectrum_monochromatic(a[0], 1.0);
    } else { fprintf(stderr, "unknown model: %s\n", model); return 1; }

    XYZ xyz = ls_spectrum_to_xyz(&s);
    double cx, cy;
    ls_xyz_chromaticity(xyz, &cx, &cy);
    printf("source: %s", model);
    for (int i = 0; i < na; ++i) printf(" %g", a[i]);
    printf("\n  band            %d-%d nm, %d bins at %d nm\n",
           LS_LAMBDA_MIN_NM, LS_LAMBDA_MAX_NM, LS_NBINS, LS_SPECTRAL_STEP_NM);
    printf("  %-6s %-8s %12.4f %s\n", ls_quantity_symbol(LS_Q_FLUX, o.units),
           "(flux)", ls_quantity_value(&s, o.units),
           ls_quantity_unit(LS_Q_FLUX, o.units));
    printf("  efficacy (band) %12.2f lm/W\n", ls_luminous_efficacy_band(&s));
    if (total > 0)
        printf("  efficacy (total)%12.2f lm/W   [out-of-band power included]\n",
               ls_luminous_efficacy_total(&s, total));
    printf("  chromaticity     (%.5f, %.5f)\n", cx, cy);
    return 0;
}

/* ------------------------------------------------------------------- grid */
typedef struct {
    SceneDesc *d;
    Opts      *o;
    vec3       nrm;
    ls_real   *val, *rad;
    int        nu, nv;
    atomic_int done;
} GridJob;

static void grid_row(int j, void *user) {
    GridJob *g = user;
    for (int i = 0; i < g->nu; ++i) {
        ls_real fu = ((ls_real)i + 0.5) / (ls_real)g->nu;
        ls_real fv = ((ls_real)j + 0.5) / (ls_real)g->nv;
        vec3 p = v3add(g->d->grid_o, v3add(v3scale(g->d->grid_u, fu),
                                           v3scale(g->d->grid_v, fv)));
        /* Seeded from the POINT index, never a thread id, so the result is
         * bit-identical no matter how the rows are scheduled. */
        Rng rng = ls_rng_seed(0x2545F4914F6CDD1Dull, (uint64_t)(j * g->nu + i) + 1);
        SpectrumAcc acc;
        ls_estimate_irradiance_full(&g->d->scene, p, g->nrm, g->o->spp,
                                    g->o->indirect, g->o->depth, &rng, &acc, NULL);
        Spectrum E = ls_acc_mean(&acc, 1);
        int idx = j * g->nu + i;
        g->val[idx] = ls_quantity_value(&E, g->o->units);
        g->rad[idx] = ls_radiometric(&E);
    }
    int n = atomic_fetch_add(&g->done, 1) + 1;
    printf("\r  %3d%%", n * 100 / g->nv);
    fflush(stdout);
}

static int cmd_grid(int argc, char **argv) {
    if (argc < 1) { fprintf(stderr, "grid needs a scene file\n"); return 1; }
    Opts o = parse_opts(argc, argv, 1, "out/grid");

    SceneDesc d;
    if (!ls_scene_load(&d, argv[0])) { fprintf(stderr, "scene error: %s\n", d.err); return 1; }
    if (o.import) {
        int added = ls_scene_import_obj(&d, o.import);
        if (added < 0) { fprintf(stderr, "import error: %s\n", d.err); return 1; }
        printf("import  %s -> %d mesh%s\n", o.import, added, added == 1 ? "" : "es");
    }
    if (!d.has_grid) { fprintf(stderr, "scene has no `grid` directive\n"); return 1; }

    int nu = d.grid_nu, nv = d.grid_nv, n = nu * nv;
    vec3 nrm = v3norm(v3cross(d.grid_u, d.grid_v));
    ls_real *val = malloc((size_t)n * sizeof *val);
    ls_real *rad = malloc((size_t)n * sizeof *rad);
    if (!val || !rad) return 1;

    const char *uname = ls_quantity_unit(LS_Q_IRRADIANCE, o.units);
    printf("scene   %s\n", argv[0]);
    printf("        %d prims, %d materials, %d lights\n", d.nprims, d.nmats, d.nlights);
    printf("grid    %d x %d = %d points, normal (%.3f %.3f %.3f)\n", nu, nv, n,
           nrm.x, nrm.y, nrm.z);
    printf("samples %d direct, %d indirect, depth %d\n", o.spp, o.indirect, o.depth);
    fflush(stdout);

    GridJob job = { .d = &d, .o = &o, .nrm = nrm, .val = val,
                    .rad = rad, .nu = nu, .nv = nv };
    atomic_init(&job.done, 0);
    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    int used = ls_parallel_for(nv, o.threads, grid_row, &job);
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    double secs = (double)(ts1.tv_sec - ts0.tv_sec)
                + (double)(ts1.tv_nsec - ts0.tv_nsec) * 1e-9;
    printf("\r        done in %.1f s on %d threads\n", secs, used);

    LsStats st = ls_stats(val, n);
    int nzero = 0;
    for (int k = 0; k < n; ++k) if (val[k] <= 0.0) nzero++;
    printf("\n%-14s %s\n", "quantity", o.units == LS_UNITS_PHOTOMETRIC
                                       ? "illuminance E_v" : "irradiance E_e");
    printf("%-14s %.4g %s\n", "min",     st.min,    uname);
    printf("%-14s %.4g %s\n", "max",     st.max,    uname);
    printf("%-14s %.4g %s\n", "mean",    st.mean,   uname);
    printf("%-14s %.4g %s\n", "std dev", st.stddev, uname);
    printf("%-14s %.4f\n",    "U0 min/mean", st.u0);
    printf("%-14s %.4f\n",    "Ud min/max",  st.ud);
    printf("%-14s %.4f\n",    "Michelson",   st.michelson);
    if (nzero > 0)
        printf("%-14s %d of %d points receive no light at all, which is why the\n"
               "%-14s uniformity ratios are zero -- they are fully occluded.\n",
               "note", nzero, n, "");

    char path[512];
    snprintf(path, sizeof path, "%s.csv", o.out);
    FILE *f = fopen(path, "w");
    if (f) {
        fprintf(f, "# lightsim grid  units=%s  nu=%d nv=%d\n", uname, nu, nv);
        fprintf(f, "# min=%.6g max=%.6g mean=%.6g stddev=%.6g U0=%.6g\n",
                st.min, st.max, st.mean, st.stddev, st.u0);
        fprintf(f, "i,j,x,y,z,value_%s,irradiance_W_m2\n",
                o.units == LS_UNITS_PHOTOMETRIC ? "lx" : "W_m2");
        for (int j = 0; j < nv; ++j)
            for (int i = 0; i < nu; ++i) {
                ls_real fu = ((ls_real)i + 0.5) / nu, fv = ((ls_real)j + 0.5) / nv;
                vec3 p = v3add(d.grid_o, v3add(v3scale(d.grid_u, fu), v3scale(d.grid_v, fv)));
                int k = j * nu + i;
                fprintf(f, "%d,%d,%.6g,%.6g,%.6g,%.8g,%.8g\n",
                        i, j, p.x, p.y, p.z, val[k], rad[k]);
            }
        fclose(f);
        printf("\nwrote %s\n", path);
    }

    snprintf(path, sizeof path, "%s.ppm", o.out);
    int up = nu < 64 ? (64 + nu - 1) / nu : 1;
    if (ls_write_falsecolor_ppm(path, val, nu, nv, st.min, st.max, up))
        printf("wrote %s  (false colour, %.4g..%.4g %s)\n", path, st.min, st.max, uname);

    snprintf(path, sizeof path, "%s.json", o.out);
    f = fopen(path, "w");
    if (f) {
        fprintf(f, "{\n  \"units\": \"%s\",\n  \"quantity\": \"%s\",\n",
                uname, o.units == LS_UNITS_PHOTOMETRIC ? "illuminance" : "irradiance");
        fprintf(f, "  \"nu\": %d,\n  \"nv\": %d,\n", nu, nv);
        fprintf(f, "  \"scene\": \"%s\",\n", argv[0]);
        fprintf(f, "  \"nlights\": %d,\n  \"nprims\": %d,\n", d.nlights, d.nprims);
        fprintf(f, "  \"spp\": %d,\n  \"indirect\": %d,\n  \"depth\": %d,\n",
                o.spp, o.indirect, o.depth);
        fprintf(f, "  \"seconds\": %.2f,\n", secs);
        fprintf(f, "  \"stats\": {\"min\": %.8g, \"max\": %.8g, \"mean\": %.8g,"
                   " \"stddev\": %.8g, \"u0\": %.6g, \"ud\": %.6g, \"michelson\": %.6g, \"occluded\": %d},\n",
                st.min, st.max, st.mean, st.stddev, st.u0, st.ud, st.michelson, nzero);
        fprintf(f, "  \"extent\": {\"u\": %.6g, \"v\": %.6g},\n",
                v3len(d.grid_u), v3len(d.grid_v));
        fprintf(f, "  \"values\": [");
        for (int k = 0; k < n; ++k)
            fprintf(f, "%s%.6g", k ? "," : "", val[k]);
        fprintf(f, "]\n}\n");
        fclose(f);
        printf("wrote %s\n", path);
    }

    if (o.blender) {
        if (ls_export_blender(&d, val, nu, nv, st.min, st.max, uname,
                              o.units == LS_UNITS_PHOTOMETRIC ? "illuminance" : "irradiance",
                              o.blender))
            printf("wrote %s  (open in Blender: Scripting > Open > Run)\n", o.blender);
        else
            fprintf(stderr, "could not write %s\n", o.blender);
    }

    free(val); free(rad);
    ls_scene_desc_free(&d);
    return 0;
}

/* ----------------------------------------------------------------- render */
typedef struct {
    SceneDesc *d;
    Opts      *o;
    Film      *film;
    atomic_int done;
} RenderJob;

static void render_row(int y, void *user) {
    RenderJob *r = user;
    int W = r->d->camera.width;
    for (int x = 0; x < W; ++x) {
        /* Seeded per pixel, so the image is identical at any thread count. */
        Rng rng = ls_rng_seed(0x853C49E6748FEA9Bull, (uint64_t)(y * W + x) + 1);
        for (int s = 0; s < r->o->spp; ++s) {
            Ray ray = ls_camera_ray(&r->d->camera, x, y, ls_rng_f(&rng), ls_rng_f(&rng));
            SpectrumAcc acc;
            ls_trace_radiance(&r->d->scene, ray, &rng, r->o->depth,
                              LS_STRAT_MIS, &acc, NULL);
            Spectrum L = ls_acc_mean(&acc, 1);
            ls_film_add(r->film, x, y, &L);
        }
    }
    int n = atomic_fetch_add(&r->done, 1) + 1;
    printf("\r  %3d%%", n * 100 / r->d->camera.height);
    fflush(stdout);
}

static int cmd_render(int argc, char **argv) {
    if (argc < 1) { fprintf(stderr, "render needs a scene file\n"); return 1; }
    Opts o = parse_opts(argc, argv, 1, "out/render.ppm");

    SceneDesc d;
    if (!ls_scene_load(&d, argv[0])) { fprintf(stderr, "scene error: %s\n", d.err); return 1; }
    if (o.import) {
        int added = ls_scene_import_obj(&d, o.import);
        if (added < 0) { fprintf(stderr, "import error: %s\n", d.err); return 1; }
        printf("import  %s -> %d mesh%s\n", o.import, added, added == 1 ? "" : "es");
    }
    if (!d.has_camera) { fprintf(stderr, "scene has no `camera` directive\n"); return 1; }

    Film film;
    if (!ls_film_init(&film, d.camera.width, d.camera.height)) return 1;
    printf("scene   %s\n", argv[0]);
    printf("render  %d x %d, %d spp, depth %d\n",
           d.camera.width, d.camera.height, o.spp, o.depth);
    fflush(stdout);

    RenderJob job = { .d = &d, .o = &o, .film = &film };
    atomic_init(&job.done, 0);
    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    int used = ls_parallel_for(d.camera.height, o.threads, render_row, &job);
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    printf("\r        done in %.1f s on %d threads\n",
           (double)(ts1.tv_sec - ts0.tv_sec)
         + (double)(ts1.tv_nsec - ts0.tv_nsec) * 1e-9, used);

    if (ls_film_write_ppm(&film, o.out, 0.0)) printf("wrote %s\n", o.out);
    char pfm[512];
    snprintf(pfm, sizeof pfm, "%s.pfm", o.out);
    if (ls_film_write_pfm(&film, pfm))
        printf("wrote %s  (raw radiance, physical units)\n", pfm);

    ls_film_free(&film);
    ls_scene_desc_free(&d);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) { usage(argv[0]); return 1; }
    if (!strcmp(argv[1], "source")) return cmd_source(argc - 2, argv + 2);
    if (!strcmp(argv[1], "grid"))   return cmd_grid(argc - 2, argv + 2);
    if (!strcmp(argv[1], "render")) return cmd_render(argc - 2, argv + 2);
    usage(argv[0]);
    return 1;
}
