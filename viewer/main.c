/* lightsim-view — an interactive window onto the simulation.
 *
 *   FIELD MAP : the measurement grid, false-coloured, with a live cross-section
 *               under the cursor and the distribution statistics beside it.
 *   RENDER    : a progressive path-traced view you can orbit; every pass adds
 *               samples to the same film, so the image refines while you watch.
 *
 * SDL is used for the window only. Everything physical comes from the same
 * library the CLI and the test suite use, so what is on screen is the simulated
 * quantity and not a preview approximation of it.
 */
#include <SDL2/SDL.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

#include "lightsim/scenefile.h"
#include "lightsim/sceneedit.h"
#include "lightsim/integrator.h"
#include "lightsim/analysis.h"
#include "lightsim/units.h"
#include "lightsim/film.h"
#include "lightsim/thread.h"
#include "lightsim/export.h"
#include "lightsim/import.h"
#include "ui.h"
#include "inspect.h"
#include "draw.h"
#include "font.h"

#define WIN_W 1360
#define WIN_H 860
#define PLOT_H 168
#define STATS_W 268
#define CBAR_H 14
#define LS_UNDO_MAX 32
#define ROW_H 17
#define GIZMO_PX 72.0            /* on-screen size of the gizmo, window pixels */
#define GIZMO_HIT 11             /* handle pick tolerance, window pixels */

/* Handles, in pick priority order: the translate axes sit inside the rotation
 * rings, so testing them first makes the inner handle win where they cross. */
typedef enum {
    GZ_NONE = 0,
    GZ_MOVE_X, GZ_MOVE_Y, GZ_MOVE_Z,
    GZ_ROT_X,  GZ_ROT_Y,  GZ_ROT_Z
} GizmoHandle;

typedef struct {
    SceneDesc d;
    Toolbar   t;
    UiState   st;

    /* ---- field map ---- */
    int       nu, nv;
    Spectrum *g_full, *g_direct;      /* one spectrum per measurement point */
    double   *disp, *other;           /* current view, and the other transport */
    double    lo, hi;
    LsStats   stats;
    int       hover_i, hover_j;
    bool      have_hover;

    /* ---- progressive render ---- */
    Film            film;
    unsigned char  *rgb;
    pthread_mutex_t lock;
    atomic_int      passes;
    atomic_int      restart;
    atomic_bool     quit;
    double          exposure;

    /* ---- orbit camera ---- */
    vec3   target;
    double dist, az, el;

    /* ---- selection ----
     * Held as indices rather than flags on the objects, because a light's index
     * IS its identity here (it is the contribution-matrix column), and removal
     * compacts the array. -1 means nothing selected. */
    int    sel_light;
    int    sel_prim;

    /* Canvas rect of the 3D view, in window coordinates, so picking and the
     * overlay share one mapping with the blit. */
    int    view_x, view_y, view_w, view_h;
    /* Canvas rect of the field map, likewise shared by its hover mapping. */
    int    fmap_x, fmap_y, fmap_side;

    /* Cached vertex buffer for the field drape, so a 64x64 grid does not
     * allocate half a megabyte every frame. */
    SDL_Vertex *drape;
    int         drape_cap;
    bool        show_drape;

    int    view;                 /* 0 = 3D perspective, 1 = top plan */
    double top_zoom;             /* plan-view framing, 1 = fit the grid */

    /* Heat shading. Not a separate view but a different thing to SHOW through
     * whichever camera is active: instead of the radiance leaving each visible
     * point, the illuminance arriving at it, false-coloured. That covers every
     * surface in the scene rather than only the measurement plane. */
    bool     shade_heat;
    double  *heat_sum;           /* accumulated per pixel */
    uint32_t *heat_n;
    double   heat_lo, heat_hi;   /* current colour range, for the bar */

    /* ---- transform gizmo ---- */
    int    gz_active;            /* GizmoHandle, or 0 for none */
    int    gz_hover;
    vec3   gz_origin;            /* object position when the drag began */
    double gz_grab;              /* axis parameter, or ring angle, at grab */

    /* ---- editing ---- */
    LsTier tier;
    UiTool tool;                 /* armed placement tool */
    LightKind new_kind;          /* what ADD LIGHT places next */
    Field  fields[LS_INSPECT_MAX];
    int    nfields;
    int    focus;                /* index into fields[], or -1 */
    char   entry[24];            /* in-progress typed value */
    int    entry_len;
    bool   typing;

    /* Whole-scene snapshots, as in the Linkage project: one per gesture,
     * pushed before the edit. Simple and correct at this scene size. */
    SceneDesc undo[LS_UNDO_MAX];
    int       undo_n;
    SceneDesc redo[LS_UNDO_MAX];
    int       redo_n;

    /* The render thread reads the scene continuously. Mutating it from the main
     * thread would be a data race, so the thread parks between passes on
     * `paused` and confirms with `acked`. Nothing takes a lock per ray. */
    atomic_bool paused, acked;

    const char *scene_path;

    /* A transient message under the canvas. Everything else this viewer says
     * ("wrote out/...", "selected light 3") goes to stdout, which a user who
     * launched from the Finder never sees. An import that silently failed
     * would be the worst instance of that. */
    char   status[128];
    Uint32 status_at;
} App;

static void say(App *a, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(a->status, sizeof a->status, fmt, ap);
    va_end(ap);
    a->status_at = SDL_GetTicks();
    printf("%s\n", a->status);
}

/* ------------------------------------------------------------------ grid */

typedef struct { App *a; vec3 n; int spp, ind, dep; } GridJob;

/* One grid row. Both transport variants are solved in the same sweep so the
 * FULL/DIRECT toggle is instant, and each point is seeded from its own index
 * so the result does not depend on how rows are scheduled. */
static void solve_grid_row(int j, void *user) {
    GridJob *g = user;
    App *ap = g->a;
    for (int i = 0; i < ap->nu; ++i) {
        double fu = (i + 0.5) / ap->nu, fv = (j + 0.5) / ap->nv;
        vec3 p = v3add(ap->d.grid_o, v3add(v3scale(ap->d.grid_u, fu),
                                           v3scale(ap->d.grid_v, fv)));
        SpectrumAcc acc;
        Rng rng = ls_rng_seed(0x2545F4914F6CDD1Dull, (uint64_t)(j * ap->nu + i) + 1);
        ls_estimate_irradiance_full(&ap->d.scene, p, g->n, g->spp, g->ind,
                                    g->dep, &rng, &acc, NULL);
        ap->g_full[j * ap->nu + i] = ls_acc_mean(&acc, 1);

        Rng r2 = ls_rng_seed(0x9E3779B97F4A7C15ull, (uint64_t)(j * ap->nu + i) + 1);
        ls_estimate_irradiance(&ap->d.scene, p, g->n, g->spp, &r2, &acc, NULL);
        ap->g_direct[j * ap->nu + i] = ls_acc_mean(&acc, 1);
    }
}

static void solve_grid(App *a) {
    GridJob job;
    job.a = a;
    job.n = v3norm(v3cross(a->d.grid_u, a->d.grid_v));
    job.spp = a->st.high_quality ? 512 : 96;
    job.ind = a->st.high_quality ? 512 : 64;
    job.dep = a->st.high_quality ? 6 : 3;
    ls_parallel_for(a->nv, 0, solve_grid_row, &job);
}

/* Project the stored spectra into the currently selected quantity. */
static void refresh_display(App *a) {
    int n = a->nu * a->nv;
    LsUnitSystem u = a->st.photometric ? LS_UNITS_PHOTOMETRIC : LS_UNITS_RADIOMETRIC;
    const Spectrum *src = a->st.direct_only ? a->g_direct : a->g_full;
    const Spectrum *oth = a->st.direct_only ? a->g_full : a->g_direct;
    for (int i = 0; i < n; ++i) {
        a->disp[i]  = ls_quantity_value(&src[i], u);
        a->other[i] = ls_quantity_value(&oth[i], u);
    }
    a->stats = ls_stats(a->disp, n);
    /* One scale across both transport modes, so toggling shows the real
     * difference instead of renormalising it away. */
    LsStats so = ls_stats(a->other, n);
    a->lo = fmin(a->stats.min, so.min);
    a->hi = fmax(a->stats.max, so.max);
    if (a->hi <= a->lo) a->hi = a->lo + 1.0;
}

/* --------------------------------------------------------------- render */

/* An orthographic plan view framed on the measurement grid, with +y up and +x
 * right -- the same orientation the heat map is drawn in, so switching between
 * the two compares geometry against result without the picture moving. */
static Camera top_camera(const App *a, int w, int h) {
    vec3 c, up_hint = v3(0, 1, 0);
    double extent;
    if (a->d.has_grid) {
        c = v3add(a->d.grid_o, v3scale(v3add(a->d.grid_u, a->d.grid_v), 0.5));
        double eu = v3len(a->d.grid_u), ev = v3len(a->d.grid_v);
        extent = eu > ev ? eu : ev;
    } else {
        c = v3(0, 0, 0);
        extent = 2.0;
    }
    if (extent < 1e-6) extent = 1.0;

    /* Height of the eye. Parking it far above the scene is the obvious choice
     * and the wrong one: in an enclosure it looks straight at the OUTSIDE of the
     * ceiling, which is unlit, so the plan view comes out black. Physically
     * right, practically useless.
     *
     * Instead sit just above the highest light -- under any ceiling, above
     * everything being lit -- so the plan shows the work plane, the parts and
     * the luminaires. */
    double top = c.z + 0.1 * extent;
    for (int i = 0; i < a->d.nlights; ++i)
        if (a->d.lights[i].p.z > top) top = a->d.lights[i].p.z;
    top += 0.04 * extent;

    vec3 eye = v3(c.x, c.y, top);
    return ls_camera_ortho(eye, c, up_hint, extent * a->top_zoom, w, h);
}

static Camera orbit_camera(const App *a, int w, int h) {
    double ce = cos(a->el), se = sin(a->el);
    vec3 eye = v3add(a->target, v3(a->dist * ce * cos(a->az),
                                   a->dist * ce * sin(a->az),
                                   a->dist * se));
    return ls_camera_look_at(eye, a->target, v3(0, 0, 1),
                             a->d.has_camera ? a->d.camera.fov_y * 180.0 / LS_PI : 42.0,
                             w, h);
}

/* The camera the active view is looking through. Everything -- rendering,
 * picking, the overlay and the gizmo -- goes through this, so a view change
 * cannot leave one of them looking somewhere else. */
static Camera active_camera(const App *a, int w, int h) {
    return (a->view == 1) ? top_camera(a, w, h) : orbit_camera(a, w, h);
}

typedef struct { App *a; Camera cam; int spp; } RJob;

static void render_row(int y, void *user) {
    RJob *r = user;
    App *a = r->a;
    int W = r->cam.width;
    int depth = a->st.high_quality ? 8 : 4;
    for (int x = 0; x < W; ++x) {
        int pass = atomic_load(&a->passes);
        Rng rng = ls_rng_seed(0x853C49E6748FEA9Bull,
                              (uint64_t)((y * W + x) * 977 + pass) + 1);
        for (int s = 0; s < r->spp; ++s) {
            Ray ray = ls_camera_ray(&r->cam, x, y, ls_rng_f(&rng), ls_rng_f(&rng));
            SpectrumAcc acc;
            ls_trace_radiance(&a->d.scene, ray, &rng, depth, LS_STRAT_MIS, &acc, NULL);
            Spectrum L = ls_acc_mean(&acc, 1);
            ls_film_add(&a->film, x, y, &L);
        }
    }
}

/* One row of the illuminance pass. For each pixel: find the surface the camera
 * sees, then measure the illuminance arriving there. The expensive part is the
 * irradiance estimate, so the per-pass sample counts are small and the result
 * accumulates over passes exactly as the radiance render does. */
static void heat_row(int y, void *user) {
    RJob *r = user;
    App *a = r->a;
    int W = r->cam.width;
    LsUnitSystem u = a->st.photometric ? LS_UNITS_PHOTOMETRIC : LS_UNITS_RADIOMETRIC;
    int nd = a->st.high_quality ? 48 : 12;
    int ni = a->st.high_quality ? 32 : 6;
    int dep = a->st.high_quality ? 4 : 2;

    for (int x = 0; x < W; ++x) {
        size_t idx = (size_t)y * (size_t)W + (size_t)x;
        int pass = atomic_load(&a->passes);
        Rng rng = ls_rng_seed(0x27BB2EE687B0B0FDull,
                              (uint64_t)((y * W + x) * 131 + pass) + 1);
        Ray ray = ls_camera_ray(&r->cam, x, y, ls_rng_f(&rng), ls_rng_f(&rng));
        Hit h;
        if (!ls_scene_intersect(&a->d.scene, &ray, &h)) continue;   /* background */

        /* Measure on the side the camera can see. */
        vec3 n = h.backface ? v3neg(h.ng) : h.ng;
        SpectrumAcc acc;
        ls_estimate_irradiance_full(&a->d.scene, h.p, n, nd, ni, dep, &rng, &acc, NULL);
        Spectrum E = ls_acc_mean(&acc, 1);
        a->heat_sum[idx] += ls_quantity_value(&E, u);
        a->heat_n[idx]++;
    }
}

/* Colour the accumulated illuminance through the same viridis ramp the field
 * map uses. The range comes from robust percentiles of the pixels that actually
 * hit something, so a single blown-out highlight cannot flatten the rest. */
static void heat_colour(App *a) {
    int W = a->film.width, H = a->film.height;
    size_t np = (size_t)W * (size_t)H;

    double *vals = malloc(np * sizeof *vals);
    if (!vals) return;
    size_t nv = 0;
    for (size_t i = 0; i < np; ++i)
        if (a->heat_n[i] > 0) vals[nv++] = a->heat_sum[i] / a->heat_n[i];
    if (nv == 0) { free(vals); return; }

    double lo = vals[0], hi = vals[0];
    for (size_t i = 1; i < nv; ++i) {
        if (vals[i] < lo) lo = vals[i];
        if (vals[i] > hi) hi = vals[i];
    }
    if (hi > lo) {
        int hist[512] = { 0 };
        for (size_t i = 0; i < nv; ++i) {
            int b = (int)((vals[i] - lo) / (hi - lo) * 511.0);
            hist[b < 0 ? 0 : (b > 511 ? 511 : b)]++;
        }
        size_t want_lo = (size_t)((double)nv * 0.02), want_hi = (size_t)((double)nv * 0.98);
        size_t run = 0;
        int b0 = 0, b1 = 511;
        for (int b = 0; b < 512; ++b) { run += (size_t)hist[b]; if (run >= want_lo) { b0 = b; break; } }
        run = 0;
        for (int b = 0; b < 512; ++b) { run += (size_t)hist[b]; if (run >= want_hi) { b1 = b; break; } }
        double nlo = lo + (hi - lo) * b0 / 511.0;
        double nhi = lo + (hi - lo) * b1 / 511.0;
        if (nhi > nlo) { lo = nlo; hi = nhi; }
    }
    free(vals);
    if (hi <= lo) hi = lo + 1.0;
    a->heat_lo = lo;
    a->heat_hi = hi;

    pthread_mutex_lock(&a->lock);
    for (size_t i = 0; i < np; ++i) {
        if (a->heat_n[i] == 0) {
            a->rgb[i*3+0] = 12; a->rgb[i*3+1] = 16; a->rgb[i*3+2] = 18;
            continue;
        }
        Uint8 r, g, b;
        draw_viridis((a->heat_sum[i] / a->heat_n[i] - lo) / (hi - lo), &r, &g, &b);
        a->rgb[i*3+0] = r; a->rgb[i*3+1] = g; a->rgb[i*3+2] = b;
    }
    pthread_mutex_unlock(&a->lock);
}

static void tonemap(App *a) {
    int W = a->film.width, H = a->film.height;
    size_t np = (size_t)W * (size_t)H;
    RGB *lin = malloc(np * sizeof *lin);
    if (!lin) return;
    double *lum = malloc(np * sizeof *lum);
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
            Spectrum s = ls_film_mean(&a->film, x, y);
            size_t i = (size_t)y * (size_t)W + (size_t)x;
            lin[i] = ls_xyz_to_linear_srgb(ls_spectrum_to_xyz(&s));
            if (lum) lum[i] = 0.2126 * lin[i].r + 0.7152 * lin[i].g + 0.0722 * lin[i].b;
        }
    /* Fix the exposure on the first pass so later passes refine the image
     * rather than making it pulse. */
    if (a->exposure <= 0.0 && lum) {
        double p = 0.0;
        for (size_t i = 0; i < np; ++i) if (lum[i] > p) p = lum[i];
        /* 92nd percentile via a coarse histogram; robust to a few hot pixels. */
        if (p > 0.0) {
            int hist[256] = { 0 };
            for (size_t i = 0; i < np; ++i) {
                int b = (int)(lum[i] / p * 255.0);
                hist[b < 0 ? 0 : (b > 255 ? 255 : b)]++;
            }
            size_t want = (size_t)((double)np * 0.92), run = 0;
            int b = 0;
            for (; b < 256; ++b) { run += (size_t)hist[b]; if (run >= want) break; }
            double v = (b + 0.5) / 256.0 * p;
            a->exposure = v > 0.0 ? 1.0 / v : 1.0;
        } else a->exposure = 1.0;
    }
    free(lum);

    pthread_mutex_lock(&a->lock);
    for (size_t i = 0; i < np; ++i) {
        RGB c = { lin[i].r * a->exposure, lin[i].g * a->exposure, lin[i].b * a->exposure };
        c = ls_rgb_gamma_encode(c);
        a->rgb[i*3+0] = (unsigned char)(ls_clamp(c.r, 0, 1) * 255.0 + 0.5);
        a->rgb[i*3+1] = (unsigned char)(ls_clamp(c.g, 0, 1) * 255.0 + 0.5);
        a->rgb[i*3+2] = (unsigned char)(ls_clamp(c.b, 0, 1) * 255.0 + 0.5);
    }
    pthread_mutex_unlock(&a->lock);
    free(lin);
}

static void *render_thread(void *arg) {
    App *a = arg;
    while (!atomic_load(&a->quit)) {
        if (atomic_load(&a->paused)) {
            atomic_store(&a->acked, true);
            SDL_Delay(2);
            continue;
        }
        atomic_store(&a->acked, false);
        if (atomic_exchange(&a->restart, 0)) {
            size_t np = (size_t)a->film.width * (size_t)a->film.height;
            memset(a->film.pix, 0, np * sizeof *a->film.pix);
            memset(a->film.n,   0, np * sizeof *a->film.n);
            if (a->heat_sum) memset(a->heat_sum, 0, np * sizeof *a->heat_sum);
            if (a->heat_n)   memset(a->heat_n,   0, np * sizeof *a->heat_n);
            atomic_store(&a->passes, 0);
            a->exposure = 0.0;
        }
        RJob job = { a, (a->view == 1) ? top_camera(a, a->film.width, a->film.height)
                                       : orbit_camera(a, a->film.width, a->film.height), 4 };
        if (a->shade_heat) {
            ls_parallel_for(a->film.height, 0, heat_row, &job);
            atomic_fetch_add(&a->passes, 1);
            heat_colour(a);
        } else {
            ls_parallel_for(a->film.height, 0, render_row, &job);
            atomic_fetch_add(&a->passes, 1);
            tonemap(a);
        }
    }
    return NULL;
}

/* A bare stage to build on: a floor, a measurement grid over it, a camera, and
 * one panel overhead so the first render is not black. */
static void make_default_scene(SceneDesc *d) {
    memset(d, 0, sizeof *d);

    Material floor_mat;
    memset(&floor_mat, 0, sizeof floor_mat);
    floor_mat.bsdf.kind = LS_BSDF_LAMBERT;
    floor_mat.bsdf.rho = ls_spectrum_const(0.5);
    int mid = ls_scene_add_material(d, floor_mat, "floor");

    Prim floor_prim;
    memset(&floor_prim, 0, sizeof floor_prim);
    floor_prim.kind = LS_PRIM_QUAD;
    floor_prim.c = v3(0, 0, 0);
    floor_prim.n = v3(0, 0, 1);
    floor_prim.ex = v3(0.5, 0, 0);
    floor_prim.ey = v3(0, 0.5, 0);
    floor_prim.mat_id = mid;
    floor_prim.light_id = -1;
    ls_scene_add_prim(d, floor_prim);

    Spectrum spd = ls_spectrum_daylight(4000.0);
    Light l = ls_light_rect(v3(0, 0, 0.7), v3(0.08, 0, 0), v3(0, -0.08, 0),
                            ls_watts_from_lumens(800.0, &spd), spd);
    l.spd_kind = LS_SPD_DAYLIGHT;
    l.spd_a = 4000.0;
    l.flux_in_lumens = true;
    l.flux_authored = 800.0;
    ls_scene_add_light(d, l);

    d->grid_o = v3(-0.5, -0.5, 0.001);
    d->grid_u = v3(1.0, 0, 0);
    d->grid_v = v3(0, 1.0, 0);
    d->grid_nu = d->grid_nv = 48;
    d->has_grid = true;

    d->cam_eye = v3(0, -1.6, 0.9);
    d->cam_target = v3(0, 0, 0.1);
    d->cam_fov_deg = 42.0;
    d->camera = ls_camera_look_at(d->cam_eye, d->cam_target, v3(0, 0, 1),
                                  d->cam_fov_deg, 480, 360);
    d->has_camera = true;
    ls_scene_rebuild(d);
}

/* ------------------------------------------------------------------ main */

/* ------------------------------------------------------- overlay + picking */

/* Project a world point to window coordinates inside the 3D view. The render
 * film is rw x rh and is blitted into (view_x, view_y, view_w, view_h), so the
 * film pixel has to be scaled the same way the blit scales it -- otherwise the
 * gizmos drift away from the pixels they annotate as the window is resized. */
static bool project_view(const App *a, const Camera *cam, vec3 p, int *sx, int *sy) {
    ls_real fx, fy;
    if (!ls_camera_project(cam, p, &fx, &fy)) return false;
    *sx = a->view_x + (int)(fx * (ls_real)a->view_w / (ls_real)cam->width);
    *sy = a->view_y + (int)(fy * (ls_real)a->view_h / (ls_real)cam->height);
    return true;
}

/* Window coordinates -> a pick ray through the 3D view. */
static bool view_pick_ray(const App *a, const Camera *cam, int mx, int my, Ray *out) {
    if (a->view_w <= 0 || a->view_h <= 0) return false;
    int lx = mx - a->view_x, ly = my - a->view_y;
    if (lx < 0 || ly < 0 || lx >= a->view_w || ly >= a->view_h) return false;
    *out = ls_camera_pick_ray(cam,
        (ls_real)lx * (ls_real)cam->width  / (ls_real)a->view_w,
        (ls_real)ly * (ls_real)cam->height / (ls_real)a->view_h);
    return true;
}

#define LIGHT_PICK_RADIUS 14

/* Select whatever is under the cursor.
 *
 * Delta lights have no geometry to intersect, so they are picked by projecting
 * their position and taking the nearest marker within a screen radius. That is
 * tested FIRST, so a marker drawn in front of a wall wins over the wall -- if
 * geometry were tested first, a point light would be unselectable whenever
 * anything lay behind it. */
static void pick_at(App *a, const Camera *cam, int mx, int my) {
    Ray r;
    if (!view_pick_ray(a, cam, mx, my, &r)) return;

    int best = -1;
    double best_d2 = (double)(LIGHT_PICK_RADIUS * LIGHT_PICK_RADIUS);
    for (int i = 0; i < a->d.nlights; ++i) {
        int sx, sy;
        if (!project_view(a, cam, a->d.lights[i].p, &sx, &sy)) continue;
        double dx = sx - mx, dy = sy - my, d2 = dx * dx + dy * dy;
        if (d2 <= best_d2) { best_d2 = d2; best = i; }
    }
    if (best >= 0) {
        a->sel_light = best;
        a->sel_prim = -1;
        printf("selected light %d\n", best);
        return;
    }

    Hit h;
    if (ls_scene_intersect(&a->d.scene, &r, &h)) {
        /* A hit on an area light's own geometry selects the LIGHT, not the
         * quad -- clicking a luminaire should give you the luminaire. */
        if (h.light_id >= 0) { a->sel_light = h.light_id; a->sel_prim = -1;
                               printf("selected light %d\n", h.light_id); }
        else                 { a->sel_light = -1; a->sel_prim = h.prim_id;
                               printf("selected surface %d\n", h.prim_id); }
        return;
    }
    a->sel_light = -1;
    a->sel_prim = -1;
}

/* The primitive that IS the current selection: a part, or an area light's own
 * emissive face. -1 when the selection has no geometry (a point or spot).
 *
 * Handed to ls_scene_intersect_ex as the prim to skip while dragging. A drag
 * slides the selection along the surface under the cursor -- but the
 * selection's OWN surface is under the cursor too, because that is what was
 * grabbed. Letting it answer the query makes the drag chase itself: the hit
 * lands on the near face, the result is offset along that face's normal, which
 * points back at the camera, and so the object advances toward the eye on every
 * motion event for as long as the button is held. Skipping it lands the object
 * on the surface BEHIND it, which is the surface the gesture is about. */
static int selection_prim(const App *a) {
    if (a->sel_light >= 0) return ls_scene_light_prim(&a->d, a->sel_light);
    if (a->sel_prim  >= 0) return a->sel_prim;
    return -1;
}

/* Clearance between the surface under the cursor and the dragged object's
 * CENTRE, so the object rests on that surface rather than sinking into it. A
 * sphere is offset by its own radius, which is what placement already does. */
static double drop_clearance(const App *a) {
    if (a->sel_light >= 0 && a->sel_light < a->d.nlights) {
        const Light *l = &a->d.lights[a->sel_light];
        if (l->kind == LS_LIGHT_SPHERE) return l->radius;
    } else if (a->sel_prim >= 0 && a->sel_prim < a->d.nprims) {
        const Prim *p = &a->d.prims[a->sel_prim];
        if (p->kind == LS_PRIM_SPHERE) return p->r;
    }
    return 0.02;
}

/* A closed polyline through world points. */
static void overlay_loop(SDL_Renderer *ren, const App *a, const Camera *cam,
                         const vec3 *pts, int n, Col c, Uint8 alpha) {
    int px = 0, py = 0, fx = 0, fy = 0;
    bool have_prev = false, have_first = false;
    for (int i = 0; i < n; ++i) {
        int sx, sy;
        if (!project_view(a, cam, pts[i], &sx, &sy)) { have_prev = false; continue; }
        if (!have_first) { fx = sx; fy = sy; have_first = true; }
        if (have_prev) draw_line(ren, px, py, sx, sy, c, alpha);
        px = sx; py = sy; have_prev = true;
    }
    if (have_first && have_prev) draw_line(ren, px, py, fx, fy, c, alpha);
}

static void draw_overlay(SDL_Renderer *ren, const App *a, const Camera *cam) {
    char buf[64];

    /* The measurement grid's footprint, so its relationship to the geometry is
     * visible in 3D rather than only in the field map. */
    if (a->d.has_grid) {
        vec3 o = a->d.grid_o, u = a->d.grid_u, v = a->d.grid_v;
        vec3 corner[4] = { o, v3add(o, u), v3add(v3add(o, u), v), v3add(o, v) };
        overlay_loop(ren, a, cam, corner, 4, COL_MUTED, 150);
    }

    for (int i = 0; i < a->d.nlights; ++i) {
        const Light *l = &a->d.lights[i];
        bool sel = (i == a->sel_light);
        Col c = sel ? COL_ACCENT : COL_WARN;
        Uint8 alpha = sel ? 255 : 170;

        /* Emitting extent. */
        if (l->kind == LS_LIGHT_RECT) {
            vec3 q[4] = { v3add(v3add(l->p, l->ex), l->ey),
                          v3sub(v3add(l->p, l->ex), l->ey),
                          v3sub(v3sub(l->p, l->ex), l->ey),
                          v3add(v3sub(l->p, l->ex), l->ey) };
            overlay_loop(ren, a, cam, q, 4, c, alpha);
        } else if (l->kind == LS_LIGHT_DISK || l->kind == LS_LIGHT_SPHERE) {
            vec3 ring[24];
            Basis b = ls_basis(l->kind == LS_LIGHT_DISK ? l->n : v3(0, 0, 1));
            for (int k = 0; k < 24; ++k) {
                double t = LS_TWO_PI * k / 24.0;
                ring[k] = v3add(l->p, v3add(v3scale(b.t, l->radius * cos(t)),
                                            v3scale(b.b, l->radius * sin(t))));
            }
            overlay_loop(ren, a, cam, ring, 24, c, alpha);
        }

        /* Aim, and the outer cone for a spot. */
        if (l->kind == LS_LIGHT_SPOT || l->kind == LS_LIGHT_DIRECTIONAL
            || l->kind == LS_LIGHT_RECT || l->kind == LS_LIGHT_DISK) {
            double len = (l->kind == LS_LIGHT_SPOT) ? 0.25 : 0.12;
            int x0, y0, x1, y1;
            if (project_view(a, cam, l->p, &x0, &y0) &&
                project_view(a, cam, v3add(l->p, v3scale(l->n, len)), &x1, &y1))
                draw_line(ren, x0, y0, x1, y1, c, alpha);
        }
        if (l->kind == LS_LIGHT_SPOT) {
            double half = acos(ls_clamp(l->cos_total, -1.0, 1.0));
            double len = 0.4, rad = len * tan(half);
            Basis b = ls_basis(l->n);
            vec3 base = v3add(l->p, v3scale(l->n, len));
            vec3 ring[20];
            for (int k = 0; k < 20; ++k) {
                double t = LS_TWO_PI * k / 20.0;
                ring[k] = v3add(base, v3add(v3scale(b.t, rad * cos(t)),
                                            v3scale(b.b, rad * sin(t))));
            }
            overlay_loop(ren, a, cam, ring, 20, c, alpha / 2);
            int x0, y0, x1, y1;
            if (project_view(a, cam, l->p, &x0, &y0))
                for (int k = 0; k < 20; k += 5)
                    if (project_view(a, cam, ring[k], &x1, &y1))
                        draw_line(ren, x0, y0, x1, y1, c, alpha / 2);
        }

        int sx, sy;
        if (project_view(a, cam, l->p, &sx, &sy)) {
            draw_marker(ren, sx, sy, sel ? 7 : 5, c, alpha);
            snprintf(buf, sizeof buf, "L%d", i);
            draw_text(ren, sx + 10, sy - 4, 1, buf, c, alpha);
        }
    }

    /* Outline the selected surface, so a material edit has a visible target. */
    if (a->sel_prim >= 0 && a->sel_prim < a->d.nprims) {
        const Prim *p = &a->d.prims[a->sel_prim];
        if (p->kind == LS_PRIM_QUAD) {
            vec3 q[4] = { v3add(v3add(p->c, p->ex), p->ey),
                          v3sub(v3add(p->c, p->ex), p->ey),
                          v3sub(v3sub(p->c, p->ex), p->ey),
                          v3add(v3sub(p->c, p->ex), p->ey) };
            overlay_loop(ren, a, cam, q, 4, COL_ACCENT, 255);
        } else if (p->kind == LS_PRIM_MESH && p->mesh_id >= 0 &&
                   p->mesh_id < a->d.nmeshes && a->d.meshes[p->mesh_id]) {
            /* An imported mesh has no single outline worth drawing, so show its
             * bounding box, carried through the same object-to-world basis the
             * intersector uses. It also makes the transform visible: a rotated
             * mesh has a visibly rotated box. */
            const Mesh *m = a->d.meshes[p->mesh_id];
            vec3 lo = m->lo, hi = m->hi, corner[8];
            for (int k = 0; k < 8; ++k) {
                vec3 o = v3((k & 1) ? hi.x : lo.x,
                            (k & 2) ? hi.y : lo.y,
                            (k & 4) ? hi.z : lo.z);
                corner[k] = v3add(p->c, v3add(v3scale(p->ex, o.x),
                                   v3add(v3scale(p->ey, o.y), v3scale(p->n, o.z))));
            }
            /* Four faces trace all twelve edges without repeating one. */
            const int face[3][4] = { {0,1,3,2}, {4,5,7,6}, {0,1,5,4} };
            for (int fi = 0; fi < 3; ++fi) {
                vec3 q[4] = { corner[face[fi][0]], corner[face[fi][1]],
                              corner[face[fi][2]], corner[face[fi][3]] };
                overlay_loop(ren, a, cam, q, 4, COL_ACCENT, 255);
            }
            const int rung[2][2] = { {2,6}, {3,7} };
            for (int ri = 0; ri < 2; ++ri) {
                int x0, y0, x1, y1;
                if (project_view(a, cam, corner[rung[ri][0]], &x0, &y0) &&
                    project_view(a, cam, corner[rung[ri][1]], &x1, &y1))
                    draw_line(ren, x0, y0, x1, y1, COL_ACCENT, 255);
            }
        } else {
            int sx, sy;
            if (project_view(a, cam, p->c, &sx, &sy))
                draw_marker(ren, sx, sy, 8, COL_ACCENT, 255);
        }
    }
}

/* Defined with the editing helpers below; the gizmo and inspector commit
 * through them. */
static void scene_pause(App *a);
static void scene_resume(App *a);
static void solve_and_refresh(App *a);
static void discard_last_undo(App *a);

/* Draw the measured field where it was actually measured: one coloured quad per
 * measurement point, laid on the grid plane in the 3D view. This is the same
 * data and the same viridis ramp as the field map beside it, so the two panes
 * are visibly one result rather than two pictures.
 *
 * Semi-transparent, so the path-traced geometry stays readable underneath and
 * it reads as a measurement laid over the scene rather than as paint. */
static void draw_field_drape(SDL_Renderer *ren, App *a, const Camera *cam) {
    if (!a->d.has_grid || a->nu <= 0 || a->nv <= 0 || !a->disp) return;

    /* Cap the drawn resolution: past a certain density the cells are smaller
     * than a pixel and the triangles are wasted work. */
    int step = 1;
    while ((a->nu / step) * (a->nv / step) > 4096) step++;
    int nu = a->nu / step, nv = a->nv / step;
    if (nu < 1 || nv < 1) return;

    int need = nu * nv * 6;
    if (need > a->drape_cap) {
        SDL_Vertex *nv2 = realloc(a->drape, (size_t)need * sizeof *nv2);
        if (!nv2) return;
        a->drape = nv2;
        a->drape_cap = need;
    }

    double span = (a->hi > a->lo) ? (a->hi - a->lo) : 1.0;
    int n = 0;
    for (int j = 0; j < nv; ++j) {
        for (int i = 0; i < nu; ++i) {
            double u0 = (double)(i * step) / a->nu, u1 = (double)((i + 1) * step) / a->nu;
            double v0 = (double)(j * step) / a->nv, v1 = (double)((j + 1) * step) / a->nv;
            vec3 c[4] = {
                v3add(a->d.grid_o, v3add(v3scale(a->d.grid_u, u0), v3scale(a->d.grid_v, v0))),
                v3add(a->d.grid_o, v3add(v3scale(a->d.grid_u, u1), v3scale(a->d.grid_v, v0))),
                v3add(a->d.grid_o, v3add(v3scale(a->d.grid_u, u1), v3scale(a->d.grid_v, v1))),
                v3add(a->d.grid_o, v3add(v3scale(a->d.grid_u, u0), v3scale(a->d.grid_v, v1)))
            };
            SDL_FPoint p[4];
            bool ok = true;
            for (int k = 0; k < 4; ++k) {
                int sx, sy;
                if (!project_view(a, cam, c[k], &sx, &sy)) { ok = false; break; }
                p[k].x = (float)sx; p[k].y = (float)sy;
            }
            if (!ok) continue;

            /* Sample the cell's value at its centre. */
            int si = i * step + step / 2, sj = j * step + step / 2;
            if (si >= a->nu) si = a->nu - 1;
            if (sj >= a->nv) sj = a->nv - 1;
            Uint8 r, g, b;
            draw_viridis((a->disp[sj * a->nu + si] - a->lo) / span, &r, &g, &b);
            SDL_Color col = { r, g, b, 205 };

            const int tri[6] = { 0, 1, 2, 0, 2, 3 };
            for (int k = 0; k < 6; ++k) {
                a->drape[n].position = p[tri[k]];
                a->drape[n].color = col;
                a->drape[n].tex_coord.x = 0.0f;
                a->drape[n].tex_coord.y = 0.0f;
                n++;
            }
        }
    }
    if (n > 0) SDL_RenderGeometry(ren, NULL, a->drape, n, NULL, 0);
}

/* -------------------------------------------------------------- gizmo ---- */

static const vec3 GZ_AXIS[3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

/* Centre of whatever is selected, or false if nothing is. */
static bool gizmo_center(const App *a, vec3 *c) {
    if (a->sel_light >= 0 && a->sel_light < a->d.nlights) { *c = a->d.lights[a->sel_light].p; return true; }
    if (a->sel_prim  >= 0 && a->sel_prim  < a->d.nprims)  { *c = a->d.prims[a->sel_prim].c;   return true; }
    return false;
}

/* World length that projects to GIZMO_PX window pixels at the gizmo's depth, so
 * the handles stay the same size on screen however far away the object is. */
static double gizmo_scale(const App *a, const Camera *cam, vec3 c) {
    double px_film = GIZMO_PX * (double)cam->width / (double)(a->view_w > 0 ? a->view_w : 1);
    return px_film * ls_camera_world_per_pixel(cam, c);
}

/* Distance in window pixels from (mx,my) to the projected segment ab. */
static double seg_dist_px(const App *a, const Camera *cam, vec3 p0, vec3 p1,
                          int mx, int my) {
    int x0, y0, x1, y1;
    if (!project_view(a, cam, p0, &x0, &y0)) return 1e9;
    if (!project_view(a, cam, p1, &x1, &y1)) return 1e9;
    double vx = x1 - x0, vy = y1 - y0;
    double L2 = vx * vx + vy * vy;
    double t = (L2 > 0.0) ? ((mx - x0) * vx + (my - y0) * vy) / L2 : 0.0;
    t = ls_clamp(t, 0.0, 1.0);
    double dx = mx - (x0 + t * vx), dy = my - (y0 + t * vy);
    return sqrt(dx * dx + dy * dy);
}

static void ring_points(vec3 c, vec3 axis, double r, vec3 *out, int n) {
    Basis b = ls_basis(axis);
    for (int i = 0; i < n; ++i) {
        double t = LS_TWO_PI * i / (double)n;
        out[i] = v3add(c, v3add(v3scale(b.t, r * cos(t)), v3scale(b.b, r * sin(t))));
    }
}

static int gizmo_pick(const App *a, const Camera *cam, int mx, int my) {
    vec3 c;
    if (!gizmo_center(a, &c)) return GZ_NONE;
    double L = gizmo_scale(a, cam, c);
    int best = GZ_NONE;
    double best_d = (double)GIZMO_HIT;

    for (int i = 0; i < 3; ++i) {
        double d = seg_dist_px(a, cam, c, v3add(c, v3scale(GZ_AXIS[i], L)), mx, my);
        if (d < best_d) { best_d = d; best = GZ_MOVE_X + i; }
    }
    /* Rings are only offered where rotating means something. */
    bool rotatable = (a->sel_light >= 0)
                   ? (a->d.lights[a->sel_light].kind != LS_LIGHT_POINT &&
                      a->d.lights[a->sel_light].kind != LS_LIGHT_SPHERE)
                   : (a->sel_prim >= 0 && a->d.prims[a->sel_prim].kind != LS_PRIM_SPHERE);
    if (rotatable) {
        vec3 pts[28];
        for (int i = 0; i < 3; ++i) {
            ring_points(c, GZ_AXIS[i], L * 0.92, pts, 28);
            for (int k = 0; k < 28; ++k) {
                double d = seg_dist_px(a, cam, pts[k], pts[(k + 1) % 28], mx, my);
                if (d < best_d) { best_d = d; best = GZ_ROT_X + i; }
            }
        }
    }
    return best;
}

static void gizmo_draw(SDL_Renderer *ren, const App *a, const Camera *cam) {
    vec3 c;
    if (!gizmo_center(a, &c)) return;
    double L = gizmo_scale(a, cam, c);
    const Col axis_col[3] = { COL_AXIS_X, COL_AXIS_Y, COL_AXIS_Z };
    const char *axis_name[3] = { "X", "Y", "Z" };

    bool rotatable = (a->sel_light >= 0)
                   ? (a->d.lights[a->sel_light].kind != LS_LIGHT_POINT &&
                      a->d.lights[a->sel_light].kind != LS_LIGHT_SPHERE)
                   : (a->sel_prim >= 0 && a->d.prims[a->sel_prim].kind != LS_PRIM_SPHERE);

    if (rotatable) {
        vec3 pts[28];
        for (int i = 0; i < 3; ++i) {
            bool hot = (a->gz_active == GZ_ROT_X + i) || (a->gz_hover == GZ_ROT_X + i);
            ring_points(c, GZ_AXIS[i], L * 0.92, pts, 28);
            overlay_loop(ren, a, cam, pts, 28, axis_col[i], hot ? 255 : 90);
        }
    }

    int cx0, cy0;
    if (!project_view(a, cam, c, &cx0, &cy0)) return;
    for (int i = 0; i < 3; ++i) {
        bool hot = (a->gz_active == GZ_MOVE_X + i) || (a->gz_hover == GZ_MOVE_X + i);
        vec3 tip = v3add(c, v3scale(GZ_AXIS[i], L));
        int tx, ty;
        if (!project_view(a, cam, tip, &tx, &ty)) continue;
        Uint8 al = hot ? 255 : 200;
        draw_line(ren, cx0, cy0, tx, ty, axis_col[i], al);
        /* A second line one pixel over reads as a thicker stroke without
         * needing a polygon. */
        draw_line(ren, cx0, cy0 + 1, tx, ty + 1, axis_col[i], al);
        draw_marker(ren, tx, ty, hot ? 5 : 3, axis_col[i], al);
        draw_text(ren, tx + 7, ty - 3, 1, axis_name[i], axis_col[i], al);
    }
}

/* Move whatever is selected to `np`. A part that is an area light's own face
 * moves the light, so the pairing cannot come apart. */
static void move_selection(App *a, vec3 np) {
    scene_pause(a);
    if (a->sel_light >= 0) {
        a->d.lights[a->sel_light].p = np;
        ls_scene_update_light(&a->d, a->sel_light);
    } else if (a->sel_prim >= 0) {
        int lid = a->d.prims[a->sel_prim].light_id;
        if (lid >= 0) {
            a->d.lights[lid].p = np;
            ls_scene_update_light(&a->d, lid);
        } else {
            a->d.prims[a->sel_prim].c = np;
            ls_scene_rebuild(&a->d);
        }
    }
    scene_resume(a);
}

static void gizmo_begin(App *a, const Camera *cam, int mx, int my, int handle) {
    Ray r;
    vec3 c;
    if (!gizmo_center(a, &c) || !view_pick_ray(a, cam, mx, my, &r)) return;
    a->gz_active = handle;
    a->gz_origin = c;
    a->gz_grab = 0.0;
    if (handle <= GZ_MOVE_Z) {
        vec3 ax = GZ_AXIS[handle - GZ_MOVE_X];
        double t;
        if (ls_line_closest_t(c, ax, r.o, r.d, &t)) a->gz_grab = t;
    } else {
        vec3 ax = GZ_AXIS[handle - GZ_ROT_X];
        vec3 q;
        if (ls_ray_plane(r.o, r.d, c, ax, &q)) {
            Basis b = ls_basis(ax);
            vec3 w = v3sub(q, c);
            a->gz_grab = atan2(v3dot(w, b.b), v3dot(w, b.t));
        }
    }
}

static void gizmo_drag(App *a, const Camera *cam, int mx, int my) {
    Ray r;
    if (a->gz_active == GZ_NONE || !view_pick_ray(a, cam, mx, my, &r)) return;

    if (a->gz_active <= GZ_MOVE_Z) {
        /* Slide along the axis: the closest point on the axis line to the pick
         * ray, minus where it was when the handle was grabbed. */
        vec3 ax = GZ_AXIS[a->gz_active - GZ_MOVE_X];
        double t;
        if (!ls_line_closest_t(a->gz_origin, ax, r.o, r.d, &t)) return;
        move_selection(a, v3add(a->gz_origin, v3scale(ax, t - a->gz_grab)));
        return;
    }

    vec3 c;
    if (!gizmo_center(a, &c)) return;
    vec3 ax = GZ_AXIS[a->gz_active - GZ_ROT_X];
    vec3 q;
    if (!ls_ray_plane(r.o, r.d, c, ax, &q)) return;
    Basis b = ls_basis(ax);
    vec3 w = v3sub(q, c);
    double ang = atan2(v3dot(w, b.b), v3dot(w, b.t));
    double delta = ang - a->gz_grab;
    /* Applied incrementally, so dragging past the wrap point keeps turning the
     * same way instead of snapping back a full revolution. */
    while (delta >  LS_PI) delta -= LS_TWO_PI;
    while (delta < -LS_PI) delta += LS_TWO_PI;
    scene_pause(a);
    if (a->sel_light >= 0)     ls_scene_rotate_light(&a->d, a->sel_light, ax, delta);
    else if (a->sel_prim >= 0) ls_scene_rotate_prim(&a->d, a->sel_prim, ax, delta);
    scene_resume(a);
    a->gz_grab = ang;
}

/* ---------------------------------------------------------- inspector ---- */

#define INSP_TOP 212

static void insp_rect(const App *a, int *x, int *y, int *w, int *h) {
    *x = WIN_W - STATS_W - 12;
    *y = 44 + INSP_TOP + ((a->view == 2) ? 100 : 0);
    *w = STATS_W;
    *h = (WIN_H - PLOT_H - 26) - *y;
}

static int insp_row_at(const App *a, int mx, int my) {
    int x, y, w, h;
    insp_rect(a, &x, &y, &w, &h);
    if (mx < x || mx >= x + w || my < y + 24 || my >= y + h) return -1;
    int i = (my - (y + 24)) / ROW_H;
    return (i >= 0 && i < a->nfields) ? i : -1;
}

/* Move a value by a horizontal drag. A wide positive range (a flux in lumens)
 * scrubs multiplicatively so the whole range is reachable; everything else
 * scrubs linearly across its own span. */
static double scrub_value(const Field *f, double v, int dx) {
    if (f->is_enum) return v;
    if (f->lo >= 0.0 && f->hi > 1000.0) {
        double nv = (v > 1e-9 ? v : 1e-3) * pow(1.012, (double)dx);
        return ls_clamp(nv, f->lo, f->hi);
    }
    return ls_clamp(v + (double)dx * (f->hi - f->lo) * 0.0025, f->lo, f->hi);
}

static void insp_rebuild(App *a) {
    a->nfields = ls_inspect_fields(&a->d, a->sel_light, a->sel_prim,
                                   a->st.photometric ? LS_UNITS_PHOTOMETRIC
                                                     : LS_UNITS_RADIOMETRIC,
                                   a->tier, a->fields, LS_INSPECT_MAX);
    if (a->focus >= a->nfields) a->focus = -1;
}

static void insp_commit(App *a, double v) {
    if (a->focus < 0 || a->focus >= a->nfields) return;
    Field *f = &a->fields[a->focus];
    if (f->readonly || f->heading) return;
    /* A typed value is held to the same range a scrub is. Without this, a
     * colour temperature of 2 K reaches ls_spectrum_blackbody, every visible
     * bin underflows to zero, and the shape that cannot be normalised aborts
     * the next ls_light_finalize -- an assert three calls away from the typing
     * that caused it. inspect.h documents lo/hi as "clamped on set"; the scrub
     * path kept that promise and this one did not. */
    v = ls_clamp(v, f->lo, f->hi);
    scene_pause(a);
    bool changed = ls_inspect_set(&a->d, a->sel_light, a->sel_prim, f->id, v);
    scene_resume(a);
    if (changed) { insp_rebuild(a); solve_and_refresh(a); }
    else discard_last_undo(a);
}

static void insp_draw(SDL_Renderer *ren, const App *a) {
    int x, y, w, h;
    insp_rect(a, &x, &y, &w, &h);
    draw_rect_fill(ren, x, y, w, h, COL_PANEL, 255);
    draw_rect_line(ren, x, y, w, h, COL_RULE, 255);

    const char *title = a->sel_light >= 0 ? "LIGHT" : a->sel_prim >= 0 ? "PART" : "INSPECTOR";
    draw_text(ren, x + 14, y + 9, 1, title, COL_ACCENT, 255);
    const char *tier = a->tier == LS_TIER_SIMPLE ? "SIMPLE"
                     : a->tier == LS_TIER_ADVANCED ? "ADVANCED" : "SCIENTIFIC";
    draw_text_right(ren, x + w - 14, y + 9, 1, tier, COL_MUTED, 255);

    if (a->nfields == 0) {
        draw_text(ren, x + 14, y + 40, 1,
                  (a->view == 2) ? "SWITCH TO RENDER TO EDIT" : "CLICK A LIGHT OR PART",
                  COL_MUTED, 160);
        return;
    }

    char buf[48];
    for (int i = 0; i < a->nfields; ++i) {
        const Field *f = &a->fields[i];
        int ry = y + 24 + i * ROW_H;
        if (ry + ROW_H > y + h) break;

        if (f->heading) {
            draw_line(ren, x + 12, ry + ROW_H - 4, x + w - 12, ry + ROW_H - 4,
                      COL_RULE, 255);
            draw_text(ren, x + 12, ry + 3, 1, f->label, COL_MUTED, 200);
            continue;
        }
        bool focused = (i == a->focus);
        if (focused)
            draw_rect_fill(ren, x + 6, ry - 1, w - 12, ROW_H - 1, COL_RULE, 160);

        draw_text(ren, x + 12, ry + 3, 1, f->label,
                  f->readonly ? COL_MUTED : COL_INK, f->readonly ? 190 : 255);

        if (focused && a->typing) {
            snprintf(buf, sizeof buf, "%s_", a->entry);
            draw_text_right(ren, x + w - 12, ry + 3, 1, buf, COL_ACCENT, 255);
        } else {
            ls_inspect_format(f, buf, sizeof buf);
            int vx = x + w - 12;
            if (f->unit && f->unit[0]) {
                draw_text_right(ren, vx, ry + 3, 1, f->unit, COL_MUTED, 170);
                vx -= font_text_width(f->unit, 1) + 6;
            }
            draw_text_right(ren, vx, ry + 3, 1, buf,
                            f->readonly ? COL_MUTED : COL_INK, 255);
        }
    }
    draw_text(ren, x + 12, y + h - 14, 1, "DRAG A ROW TO SCRUB   TYPE TO SET",
              COL_MUTED, 130);
}

/* ------------------------------------------------------------ editing ---- */

/* Park the render thread so the scene can be mutated safely. Checked between
 * passes, never inside one, so the tracing loop stays lock-free. */
static void scene_pause(App *a) {
    atomic_store(&a->paused, true);
    for (int i = 0; i < 500 && !atomic_load(&a->acked); ++i) SDL_Delay(1);
}
static void scene_resume(App *a) {
    atomic_store(&a->paused, false);
    atomic_store(&a->restart, 1);
}

static void solve_and_refresh(App *a) {
    if (!a->d.has_grid) return;
    solve_grid(a);
    refresh_display(a);
}

static void clear_stack(SceneDesc *st, int *n) {
    for (int i = 0; i < *n; ++i) ls_scene_desc_free(&st[i]);
    *n = 0;
}

static void push_snapshot(SceneDesc *st, int *n, const SceneDesc *src) {
    if (*n >= LS_UNDO_MAX) {
        ls_scene_desc_free(&st[0]);
        for (int i = 1; i < LS_UNDO_MAX; ++i) st[i - 1] = st[i];
        *n = LS_UNDO_MAX - 1;
    }
    if (ls_scene_clone(src, &st[*n])) (*n)++;
}

static void push_undo(App *a) {
    push_snapshot(a->undo, &a->undo_n, &a->d);
    clear_stack(a->redo, &a->redo_n);
}

/* Drop the most recent snapshot without restoring it, for the case where a
 * gesture optimistically snapshotted and then turned out to change nothing. */
static void discard_last_undo(App *a) {
    if (a->undo_n <= 0) return;
    a->undo_n--;
    ls_scene_desc_free(&a->undo[a->undo_n]);
}

static void adopt(App *a, SceneDesc *restored) {
    ls_scene_desc_free(&a->d);
    a->d = *restored;
    ls_scene_rebuild(&a->d);
    if (a->sel_light >= a->d.nlights) a->sel_light = -1;
    if (a->sel_prim  >= a->d.nprims)  a->sel_prim  = -1;
}

static void app_undo(App *a) {
    if (a->undo_n <= 0) { printf("nothing to undo\n"); return; }
    scene_pause(a);
    push_snapshot(a->redo, &a->redo_n, &a->d);
    a->undo_n--;
    adopt(a, &a->undo[a->undo_n]);
    scene_resume(a);
    solve_and_refresh(a);
    printf("undo\n");
}

static void app_redo(App *a) {
    if (a->redo_n <= 0) { printf("nothing to redo\n"); return; }
    scene_pause(a);
    /* Deliberately not push_undo: redoing must not clear the redo stack. */
    push_snapshot(a->undo, &a->undo_n, &a->d);
    a->redo_n--;
    adopt(a, &a->redo[a->redo_n]);
    scene_resume(a);
    solve_and_refresh(a);
    printf("redo\n");
}

/* A light of the currently selected kind, sized and aimed for the surface it
 * was dropped on. A click always produces a working luminaire rather than
 * nothing, in the spirit of the Linkage tools. */
static Light make_light(const App *a, vec3 p, vec3 aim) {
    Spectrum spd = ls_spectrum_daylight(4000.0);
    double watts = ls_watts_from_lumens(400.0, &spd);
    Light l;
    switch (a->new_kind) {
        case LS_LIGHT_SPOT:
            l = ls_light_beam(p, aim, 40.0, watts, spd); break;
        case LS_LIGHT_POINT:
            l = ls_light_point(p, watts, spd); break;
        case LS_LIGHT_SPHERE:
            l = ls_light_sphere(p, 0.03, watts, spd); break;
        case LS_LIGHT_DISK:
            l = ls_light_disk(p, aim, 0.05, watts, spd); break;
        case LS_LIGHT_DIRECTIONAL:
            l = ls_light_directional(aim, 5.0, spd); break;
        case LS_LIGHT_RECT:
        default: {
            /* Build the rect's edges in the plane perpendicular to the aim, so
             * a panel dropped on a wall lies flat against it. */
            Basis b = ls_basis(v3neg(v3norm(aim)));
            l = ls_light_rect(p, v3scale(b.t, 0.05), v3scale(b.b, -0.05), watts, spd);
            break;
        }
    }
    l.spd_kind = LS_SPD_DAYLIGHT;
    l.spd_a = 4000.0;
    l.flux_in_lumens = true;
    l.flux_authored = 400.0;
    return l;
}

static void app_place(App *a, const Camera *cam, int mx, int my) {
    Ray r;
    if (!view_pick_ray(a, cam, mx, my, &r)) return;
    Hit h;
    vec3 p, ng;
    if (ls_scene_intersect(&a->d.scene, &r, &h)) {
        ng = h.backface ? v3neg(h.ng) : h.ng;
        p = v3add(h.p, v3scale(ng, 0.02));
    } else {
        /* Nothing under the cursor: drop it a little way down the ray rather
         * than doing nothing, so a click always produces something. */
        p = v3add(r.o, v3scale(r.d, a->dist));
        ng = v3(0, 0, 1);
    }
    scene_pause(a);
    push_undo(a);
    if (a->tool == UI_TOOL_LIGHT) {
        int id = ls_scene_add_light(&a->d, make_light(a, p, v3neg(ng)));
        a->sel_light = id; a->sel_prim = -1;
        printf("placed light %d\n", id);
    } else {
        Material m;
        memset(&m, 0, sizeof m);
        m.bsdf.kind = LS_BSDF_LAMBERT;
        m.bsdf.rho = ls_spectrum_const(0.5);
        char nm[32];
        snprintf(nm, sizeof nm, "part%d", a->d.nprims);
        int mid = ls_scene_add_material(&a->d, m, nm);
        Prim pr;
        memset(&pr, 0, sizeof pr);
        pr.kind = LS_PRIM_SPHERE;
        pr.c = v3add(p, v3scale(ng, 0.03));
        pr.r = 0.05;
        pr.mat_id = mid;
        pr.light_id = -1;
        int id = ls_scene_add_prim(&a->d, pr);
        a->sel_prim = id; a->sel_light = -1;
        printf("placed part %d\n", id);
    }
    a->tool = UI_TOOL_NONE;
    scene_resume(a);
    solve_and_refresh(a);
}

static void app_delete(App *a) {
    if (a->sel_light < 0 && a->sel_prim < 0) return;
    scene_pause(a);
    push_undo(a);
    if (a->sel_light >= 0) { ls_scene_remove_light(&a->d, a->sel_light); a->sel_light = -1; }
    else                   { ls_scene_remove_prim(&a->d, a->sel_prim);   a->sel_prim = -1; }
    scene_resume(a);
    solve_and_refresh(a);
}

static void app_duplicate(App *a) {
    if (a->sel_light < 0) return;
    scene_pause(a);
    push_undo(a);
    int id = ls_scene_duplicate_light(&a->d, a->sel_light, v3(0.05, 0.05, 0.0));
    if (id >= 0) a->sel_light = id;
    scene_resume(a);
    solve_and_refresh(a);
}

static void save_outputs(App *a) {
    if (a->view == 2) {
        int up = a->nu < 64 ? (64 + a->nu - 1) / a->nu : 1;
        if (ls_write_falsecolor_ppm("out/view_field.ppm", a->disp, a->nu, a->nv,
                                    a->lo, a->hi, up))
            printf("wrote out/view_field.ppm\n");
    } else {
        if (ls_film_write_ppm(&a->film, "out/view_render.ppm", a->exposure))
            printf("wrote out/view_render.ppm\n");
    }
}

/* Write the edited scene back beside the one it was loaded from, so a session
 * can be reopened in this tool and also re-run headlessly from the CLI. The
 * original is never overwritten -- an edit should not silently rewrite the file
 * the user authored. */
static void save_scene(App *a) {
    char path[512];
    snprintf(path, sizeof path, "%s", a->scene_path);
    size_t n = strlen(path);
    const char *suffix = ".edited.scene";
    if (n > 6 && strcmp(path + n - 6, ".scene") == 0) path[n - 6] = '\0';
    strncat(path, suffix, sizeof path - strlen(path) - 1);
    if (ls_scene_save(&a->d, path)) printf("wrote %s\n", path);
    else fprintf(stderr, "could not write %s\n", path);
}

/* Bring in an OBJ, or open a whole .scene. Goes through the same
 * pause/snapshot/resume discipline as every other edit, so ^Z undoes an import
 * in one step and the render thread is never reading a half-built scene. */
static void app_import(App *a, const char *path) {
    size_t n = strlen(path);
    if (n > 6 && strcmp(path + n - 6, ".scene") == 0) {
        SceneDesc nd;
        if (!ls_scene_load(&nd, path)) { say(a, "%s", nd.err); return; }
        scene_pause(a);
        push_undo(a);
        ls_scene_desc_free(&a->d);
        a->d = nd;
        ls_scene_rebuild(&a->d);
        a->sel_light = a->d.nlights > 0 ? 0 : -1;
        a->sel_prim = -1;
        scene_resume(a);
        solve_and_refresh(a);
        say(a, "OPENED %s", path);
        return;
    }
    scene_pause(a);
    push_undo(a);
    int added = ls_scene_import_obj(&a->d, path);
    scene_resume(a);
    if (added < 0) {
        app_undo(a);              /* nothing was added; drop the snapshot */
        say(a, "%s", a->d.err);
        return;
    }
    /* Select the first thing that arrived, so the inspector has something to
     * show and the gizmo lands on it. */
    for (int i = a->d.nprims - 1; i >= 0; --i)
        if (a->d.prims[i].kind == LS_PRIM_MESH) { a->sel_prim = i; a->sel_light = -1; }
    solve_and_refresh(a);
    say(a, "IMPORTED %d MESH%s FROM %s", added, added == 1 ? "" : "ES", path);
}

/* A file chooser without a dependency. SDL2 has none, and this project has no
 * third-party code beyond SDL, so on macOS ask the system for one. Anywhere
 * else, say so and point at drag and drop, which always works. */
static void app_import_dialog(App *a) {
#ifdef __APPLE__
    FILE *p = popen("osascript -e 'POSIX path of (choose file with prompt "
                    "\"Import geometry\" of type {\"obj\", \"scene\"})' "
                    "2>/dev/null", "r");
    if (!p) { say(a, "COULD NOT OPEN A FILE CHOOSER -- DRAG A FILE IN"); return; }
    char path[1024];
    if (fgets(path, sizeof path, p)) {
        size_t n = strlen(path);
        while (n > 0 && (path[n-1] == '\n' || path[n-1] == '\r')) path[--n] = '\0';
        pclose(p);
        if (n > 0) app_import(a, path);
        return;
    }
    pclose(p);
    say(a, "IMPORT CANCELLED");
#else
    (void)a;
    say(a, "DRAG AN .OBJ ONTO THE WINDOW TO IMPORT");
#endif
}

static void export_blender(App *a) {
    const char *unit = ls_quantity_unit(LS_Q_IRRADIANCE,
        a->st.photometric ? LS_UNITS_PHOTOMETRIC : LS_UNITS_RADIOMETRIC);
    if (ls_export_blender(&a->d, a->disp, a->nu, a->nv, a->lo, a->hi, unit,
                          a->st.photometric ? "illuminance" : "irradiance",
                          "out/view_scene.py"))
        printf("wrote out/view_scene.py  (Blender: Scripting > Open > Run)\n");
}

int main(int argc, char **argv) {
    const char *scene_arg = NULL;
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
            printf("usage: %s [scene.scene] [--3d|--top|--heat] [--fine]\n\n", argv[0]);
            printf("  With no scene file, opens an empty stage to build on.\n\n");
            printf("  VIEW    1 3D perspective   2 orthographic plan\n");
    printf("  SHOW    3 illuminance on every surface (toggle)\n");
            printf("  EDIT    A add light   P add part   I import OBJ   D duplicate\n");
            printf("          DEL delete   (or drag an .obj/.scene onto the window)\n");
            printf("          TAB cycle selection   ^Z undo   ^Y redo\n");
            printf("          U lux/watt   T full/direct   Q draft/fine\n");
            printf("          V tier: simple / advanced / scientific\n");
            printf("          F drape the measured field on the geometry\n");
            printf("  FILE    R re-solve   S save PPM   W save scene   B export Blender\n");
            printf("  ESC     cancel, then deselect, then quit\n\n");
            printf("  Click to select. Drag the selection to move it, or use the axis\n");
            printf("  handles and rotation rings on the gizmo. Drag elsewhere to orbit\n");
            printf("  (3D view only). In the inspector, drag a row to scrub or type.\n");
            return 0;
        }
        if (argv[i][0] != '-' && !scene_arg) scene_arg = argv[i];
    }
    int want_view = -1;
    bool want_fine = false;
    for (int i = 1; i < argc; ++i) {
        if      (!strcmp(argv[i], "--3d")   || !strcmp(argv[i], "--render")) want_view = 0;
        else if (!strcmp(argv[i], "--top"))  want_view = 1;
        else if (!strcmp(argv[i], "--heat") || !strcmp(argv[i], "--field"))  want_view = 2;
        else if (!strcmp(argv[i], "--fine")) want_fine = true;
    }
    App a;
    memset(&a, 0, sizeof a);
    a.scene_path = scene_arg ? scene_arg : "untitled.scene";
    if (scene_arg) {
        if (!ls_scene_load(&a.d, scene_arg)) {
            fprintf(stderr, "scene error: %s\n", a.d.err);
            return 1;
        }
    } else {
        /* Launched with no file: start on an empty stage rather than refusing,
         * so a scene can be built from nothing in the same tool. */
        make_default_scene(&a.d);
        printf("no scene given -- starting an empty stage. "
               "ADD LIGHT / ADD PART to build one, SAVE SCENE to keep it.\n");
    }
    a.st.has_grid = a.d.has_grid;
    a.st.has_camera = a.d.has_camera;
    a.view = (want_view == 1) ? 1 : 0;
    a.shade_heat = (want_view == 2);
    a.st.photometric = true;
    a.st.direct_only = false;
    a.st.high_quality = want_fine;
    a.hover_i = a.hover_j = -1;
    a.sel_light = a.sel_prim = -1;
    a.focus = -1;
    a.tier = LS_TIER_SIMPLE;
    a.tool = UI_TOOL_NONE;
    a.new_kind = LS_LIGHT_RECT;
    /* Off by default: the three views should each show a different thing --
     * geometry in perspective, geometry in plan, and the measured field. F
     * lays the field over either geometry view when the combination is what
     * you want. */
    a.show_drape = false;
    a.top_zoom = 1.0;

    if (a.d.has_grid) {
        a.nu = a.d.grid_nu; a.nv = a.d.grid_nv;
        int n = a.nu * a.nv;
        a.g_full   = calloc((size_t)n, sizeof *a.g_full);
        a.g_direct = calloc((size_t)n, sizeof *a.g_direct);
        a.disp     = calloc((size_t)n, sizeof *a.disp);
        a.other    = calloc((size_t)n, sizeof *a.other);
        printf("solving %dx%d field...\n", a.nu, a.nv);
        solve_grid(&a);
        refresh_display(&a);
        printf("  %.4g .. %.4g lx, mean %.4g\n", a.stats.min, a.stats.max, a.stats.mean);
    }

    /* Open with the first light selected, so the inspector shows something to
     * read and edit rather than an empty panel the user has to discover. */
    if (a.d.nlights > 0) a.sel_light = 0;

    /* Orbit starts wherever the scene's own camera is pointing. */
    if (a.d.has_camera) {
        vec3 e = a.d.camera.eye;
        a.target = v3add(e, v3scale(a.d.camera.fwd, 1.0));
        vec3 off = v3sub(e, a.target);
        a.dist = v3len(off);
        a.az = atan2(off.y, off.x);
        a.el = asin(ls_clamp(off.z / (a.dist > 0 ? a.dist : 1.0), -1.0, 1.0));
    } else {
        a.target = v3(0, 0, 0); a.dist = 2.0; a.az = -LS_PI / 2; a.el = 0.3;
    }

    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        fprintf(stderr, "SDL init failed: %s\n", SDL_GetError());
        return 1;
    }
    SDL_Window *win = SDL_CreateWindow("lightsim", SDL_WINDOWPOS_CENTERED,
        SDL_WINDOWPOS_CENTERED, WIN_W, WIN_H, SDL_WINDOW_SHOWN | SDL_WINDOW_ALLOW_HIGHDPI);
    SDL_Renderer *ren = SDL_CreateRenderer(win, -1,
        SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    SDL_SetRenderDrawBlendMode(ren, SDL_BLENDMODE_BLEND);
    /* With ALLOW_HIGHDPI the renderer's backing store is larger than the
     * window (2x on a Retina display), so drawing in window coordinates would
     * fill only the top-left quarter. A logical size keeps every coordinate
     * below in window units and lets SDL scale up, which also keeps the result
     * crisp rather than upscaled. */
    SDL_RenderSetLogicalSize(ren, WIN_W, WIN_H);
    SDL_EventState(SDL_DROPFILE, SDL_ENABLE);   /* drag an .obj onto the window */

    ui_init(&a.t);
    SDL_StartTextInput();

    /* Canvas geometry. */
    int cx = UI_TOOLBAR_W + 18, cy = 44;

    int rw = 512, rh = 512;
    ls_film_init(&a.film, rw, rh);
    a.rgb = calloc((size_t)rw * (size_t)rh * 3u, 1);
    a.heat_sum = calloc((size_t)rw * (size_t)rh, sizeof *a.heat_sum);
    a.heat_n   = calloc((size_t)rw * (size_t)rh, sizeof *a.heat_n);
    pthread_mutex_init(&a.lock, NULL);
    atomic_store(&a.restart, 1);
    atomic_store(&a.paused, false);
    atomic_store(&a.acked, false);
    pthread_t rt;
    pthread_create(&rt, NULL, render_thread, &a);

    /* The canvas rect is recomputed each frame with the layout; this just
     * gives picking something sane before the first one. */
    a.view_x = a.fmap_x = cx;
    a.view_y = a.fmap_y = cy;
    a.view_w = a.view_h = a.fmap_side = 400;

    SDL_Texture *rend_tex = SDL_CreateTexture(ren, SDL_PIXELFORMAT_RGB24,
        SDL_TEXTUREACCESS_STREAMING, rw, rh);
    SDL_SetTextureScaleMode(rend_tex, SDL_ScaleModeLinear);

    bool dragging = false, pending = false, moving = false;
    bool scrubbing = false, scrub_moved = false, scrub_snapped = false;
    int last_x = 0, last_y = 0, press_x = 0, press_y = 0;
    int scrub_start = 0;
    double scrub_base = 0.0;
    const int DRAG_THRESHOLD = 4;
    bool running = true;
    char buf[128];

    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) running = false;
            else if (e.type == SDL_DROPFILE) {
                /* The gesture people actually reach for. SDL hands us a string
                 * it allocated; we own it. */
                app_import(&a, e.drop.file);
                SDL_free(e.drop.file);
            }
            else if (e.type == SDL_KEYDOWN) {
                switch (e.key.keysym.sym) {
                    case SDLK_ESCAPE:
                        /* A cancel ladder, most transient state first, so one
                         * key backs out of whatever is in progress before it
                         * closes the window. */
                        if (a.typing)                { a.typing = false; a.entry_len = 0; }
                        else if (a.tool != UI_TOOL_NONE) a.tool = UI_TOOL_NONE;
                        else if (a.focus >= 0)       a.focus = -1;
                        else if (a.sel_light >= 0 || a.sel_prim >= 0)
                                                     { a.sel_light = a.sel_prim = -1; }
                        else running = false;
                        break;
                    case SDLK_1: a.view = 0; atomic_store(&a.restart, 1); break;
                    case SDLK_2: a.view = 1; atomic_store(&a.restart, 1); break;
                    case SDLK_u: a.st.photometric = !a.st.photometric; refresh_display(&a); break;
                    case SDLK_t: a.st.direct_only = !a.st.direct_only; refresh_display(&a); break;
                    case SDLK_q: a.st.high_quality = !a.st.high_quality;
                                 if (a.view == 2) { solve_grid(&a); refresh_display(&a); }
                                 else atomic_store(&a.restart, 1);
                                 break;
                    case SDLK_r: if (a.view == 2) { solve_grid(&a); refresh_display(&a); } break;
                    case SDLK_s: save_outputs(&a); break;
                    case SDLK_w: save_scene(&a); break;
                    case SDLK_3: a.shade_heat = !a.shade_heat;
                                 atomic_store(&a.restart, 1); break;
                    case SDLK_v: a.tier = (LsTier)((a.tier + 1) % 3); break;
                    case SDLK_f: a.show_drape = !a.show_drape; break;
                    case SDLK_TAB: {
                        /* Step through the lights, then the parts, then back to
                         * nothing -- so everything is reachable without hunting
                         * for it with the cursor. */
                        int nl = a.d.nlights;
                        if (a.sel_light >= 0) {
                            a.sel_light++;
                            if (a.sel_light >= nl) { a.sel_light = -1; a.sel_prim = 0; }
                        } else if (a.sel_prim >= 0) {
                            do { a.sel_prim++; }
                            while (a.sel_prim < a.d.nprims &&
                                   a.d.prims[a.sel_prim].light_id >= 0);
                            if (a.sel_prim >= a.d.nprims) a.sel_prim = -1;
                        } else if (nl > 0) {
                            a.sel_light = 0;
                        }
                        if (a.sel_light >= 0) a.sel_prim = -1;
                        a.focus = -1;
                        break;
                    }
                    case SDLK_a: if ((a.view != 2))
                                     a.tool = (a.tool == UI_TOOL_LIGHT)
                                            ? UI_TOOL_NONE : UI_TOOL_LIGHT;
                                 break;
                    case SDLK_p: if ((a.view != 2))
                                     a.tool = (a.tool == UI_TOOL_PART)
                                            ? UI_TOOL_NONE : UI_TOOL_PART;
                                 break;
                    case SDLK_d: app_duplicate(&a); break;
                    case SDLK_DELETE:
                    case SDLK_BACKSPACE:
                        if (a.typing) {
                            if (a.entry_len > 0) a.entry[--a.entry_len] = '\0';
                        } else app_delete(&a);
                        break;
                    case SDLK_RETURN:
                    case SDLK_KP_ENTER:
                        if (a.typing && a.entry_len > 0) {
                            scene_pause(&a); push_undo(&a); scene_resume(&a);
                            insp_commit(&a, atof(a.entry));
                        }
                        a.typing = false; a.entry_len = 0; a.entry[0] = '\0';
                        break;
                    case SDLK_z:
                        if (SDL_GetModState() & (KMOD_CTRL | KMOD_GUI)) {
                            if (SDL_GetModState() & KMOD_SHIFT) app_redo(&a);
                            else app_undo(&a);
                        }
                        break;
                    case SDLK_y:
                        if (SDL_GetModState() & (KMOD_CTRL | KMOD_GUI)) app_redo(&a);
                        break;
                    case SDLK_b: export_blender(&a); break;
                    case SDLK_i: app_import_dialog(&a); break;
                    default: break;
                }
            }
            else if (e.type == SDL_TEXTINPUT) {
                /* Only meaningful while a numeric row has focus. */
                if (a.focus >= 0 && a.focus < a.nfields &&
                    !a.fields[a.focus].readonly && !a.fields[a.focus].is_enum) {
                    for (const char *c = e.text.text; *c; ++c) {
                        if (!((*c >= '0' && *c <= '9') || *c == '.' || *c == '-')) continue;
                        if (a.entry_len < (int)sizeof a.entry - 1) {
                            a.entry[a.entry_len++] = *c;
                            a.entry[a.entry_len] = '\0';
                            a.typing = true;
                        }
                    }
                }
            }
            else if (e.type == SDL_MOUSEMOTION) {
                a.t.hover = ui_hit_test(&a.t, e.motion.x, e.motion.y);
                if (a.gz_active != GZ_NONE) {
                    Camera gcam = active_camera(&a, rw, rh);
                    gizmo_drag(&a, &gcam, e.motion.x, e.motion.y);
                    last_x = e.motion.x; last_y = e.motion.y;
                    break;
                }
                if ((a.view != 2) && !dragging && !moving && !scrubbing) {
                    Camera gcam = active_camera(&a, rw, rh);
                    a.gz_hover = gizmo_pick(&a, &gcam, e.motion.x, e.motion.y);
                } else a.gz_hover = GZ_NONE;
                if (scrubbing && a.focus >= 0 && a.focus < a.nfields) {
                    int dx = e.motion.x - scrub_start;
                    if (dx != 0) scrub_moved = true;
                    double nv = scrub_value(&a.fields[a.focus], scrub_base, dx);
                    /* Snapshot once, on the first movement of the gesture. */
                    if (scrub_moved && !scrub_snapped) {
                        scene_pause(&a); push_undo(&a); scene_resume(&a);
                        scrub_snapped = true;
                    }
                    scene_pause(&a);
                    ls_inspect_set(&a.d, a.sel_light, a.sel_prim,
                                   a.fields[a.focus].id, nv);
                    scene_resume(&a);
                    insp_rebuild(&a);
                }
                if (moving && (a.view != 2)) {
                    /* Slide along whatever surface is under the cursor: the
                     * surface IS the constraint, which keeps a 3D drag
                     * unambiguous without axis gizmos. The object's own face is
                     * excluded from the query, or the drag walks it into the
                     * camera -- see selection_prim. */
                    Camera cam = active_camera(&a, rw, rh);
                    Ray pr;
                    Hit hh;
                    if (view_pick_ray(&a, &cam, e.motion.x, e.motion.y, &pr) &&
                        ls_scene_intersect_ex(&a.d.scene, &pr,
                                              selection_prim(&a), &hh)) {
                        vec3 ng = hh.backface ? v3neg(hh.ng) : hh.ng;
                        move_selection(&a, v3add(hh.p, v3scale(ng, drop_clearance(&a))));
                    }
                    pending = false;
                }
                if (pending && !moving &&
                    (abs(e.motion.x - press_x) > DRAG_THRESHOLD ||
                     abs(e.motion.y - press_y) > DRAG_THRESHOLD)) {
                    pending = false;
                    dragging = true;
                }
                if (dragging && a.view == 0) {
                    a.az -= (e.motion.x - last_x) * 0.008;
                    a.el = ls_clamp(a.el + (e.motion.y - last_y) * 0.008, -1.45, 1.45);
                    atomic_store(&a.restart, 1);
                }
                last_x = e.motion.x; last_y = e.motion.y;
                a.have_hover = false;
                if (a.shade_heat && a.view_w > 0) {
                    /* Read back the value the heat pass already measured for
                     * this pixel, so probing any surface is free. */
                    int gx = e.motion.x - a.view_x, gy = e.motion.y - a.view_y;
                    if (gx >= 0 && gx < a.view_w && gy >= 0 && gy < a.view_h) {
                        a.hover_i = gx * rw / a.view_w;
                        a.hover_j = gy * rh / a.view_h;
                        a.have_hover = true;
                    }
                }
            }
            else if (e.type == SDL_MOUSEBUTTONDOWN) {
                int hit = ui_hit_test(&a.t, e.button.x, e.button.y);
                if (hit >= 0) {
                    if (a.t.buttons[hit].enabled) switch (a.t.buttons[hit].action) {
                        case UI_VIEW_3D:     a.view = 0; atomic_store(&a.restart, 1); break;
                        case UI_VIEW_TOP:    a.view = 1; atomic_store(&a.restart, 1); break;
                        case UI_VIEW_HEAT:   a.shade_heat = !a.shade_heat;
                                             atomic_store(&a.restart, 1); break;
                        case UI_UNITS:       a.st.photometric = !a.st.photometric;
                                             refresh_display(&a); break;
                        case UI_TRANSPORT:   a.st.direct_only = !a.st.direct_only;
                                             refresh_display(&a); break;
                        case UI_QUALITY:     a.st.high_quality = !a.st.high_quality;
                                             if (a.view == 2) { solve_grid(&a); refresh_display(&a); }
                                             else atomic_store(&a.restart, 1);
                                             break;
                        case UI_SOLVE:       solve_grid(&a); refresh_display(&a); break;
                        case UI_TIER:        a.tier = (LsTier)((a.tier + 1) % 3); break;
                        case UI_ADD_LIGHT:   a.tool = (a.tool == UI_TOOL_LIGHT)
                                                    ? UI_TOOL_NONE : UI_TOOL_LIGHT; break;
                        case UI_ADD_PART:    a.tool = (a.tool == UI_TOOL_PART)
                                                    ? UI_TOOL_NONE : UI_TOOL_PART; break;
                        case UI_DUPLICATE:   app_duplicate(&a); break;
                        case UI_DELETE:      app_delete(&a); break;
                        case UI_UNDO:        app_undo(&a); break;
                        case UI_REDO:        app_redo(&a); break;
                        case UI_SAVE:        save_outputs(&a); break;
                        case UI_SAVE_SCENE:  save_scene(&a); break;
                        case UI_BLENDER:     export_blender(&a); break;
                        case UI_IMPORT:      app_import_dialog(&a); break;
                        default: break;
                    }
                } else {
                    int row = insp_row_at(&a, e.button.x, e.button.y);
                    if (row >= 0) {
                        const Field *f = &a.fields[row];
                        a.typing = false; a.entry_len = 0;
                        if (f->heading || f->readonly) { a.focus = -1; }
                        else if (f->is_enum) {
                            /* An enum has no continuum to scrub: clicking it
                             * steps to the next value. */
                            a.focus = row;
                            scene_pause(&a); push_undo(&a); scene_resume(&a);
                            double nv = f->value + 1.0;
                            if (nv > f->hi) nv = f->lo;
                            insp_commit(&a, nv);
                        } else {
                            a.focus = row;
                            scrubbing = true;
                            scrub_start = e.button.x;
                            scrub_base = f->value;
                            scrub_moved = false;
                        }
                    } else if ((a.view != 2)) {
                        Camera gcam = active_camera(&a, rw, rh);
                        int handle = (a.tool == UI_TOOL_NONE)
                                   ? gizmo_pick(&a, &gcam, e.button.x, e.button.y)
                                   : GZ_NONE;
                        if (handle != GZ_NONE) {
                            scene_pause(&a); push_undo(&a); scene_resume(&a);
                            gizmo_begin(&a, &gcam, e.button.x, e.button.y, handle);
                            break;
                        }
                        pending = true;
                        press_x = e.button.x; press_y = e.button.y;
                        /* Pressing on an already-selected object begins a move;
                         * pressing anywhere else begins an orbit. Click to
                         * select, then drag it -- so orbiting never requires
                         * finding empty space in a closed scene. */
                        Camera cam = active_camera(&a, rw, rh);
                        Ray pr;
                        if (a.sel_light >= 0 &&
                            view_pick_ray(&a, &cam, e.button.x, e.button.y, &pr)) {
                            int px, py;
                            if (project_view(&a, &cam, a.d.lights[a.sel_light].p, &px, &py) &&
                                abs(px - e.button.x) < 18 && abs(py - e.button.y) < 18)
                                moving = true;
                        }
                        if (!moving && a.sel_prim >= 0 &&
                            view_pick_ray(&a, &cam, e.button.x, e.button.y, &pr)) {
                            Hit hh;
                            if (ls_scene_intersect(&a.d.scene, &pr, &hh) &&
                                hh.prim_id == a.sel_prim)
                                moving = true;
                        }
                        if (moving) { scene_pause(&a); push_undo(&a); scene_resume(&a); }
                    }
                }
            }
            else if (e.type == SDL_MOUSEBUTTONUP) {
                if (a.gz_active != GZ_NONE) {
                    a.gz_active = GZ_NONE;
                    solve_and_refresh(&a);
                    pending = false;
                    dragging = false;
                    break;
                }
                if (scrubbing) {
                    if (scrub_moved) solve_and_refresh(&a);
                    else if (scrub_snapped) discard_last_undo(&a);
                    scrubbing = false; scrub_moved = false; scrub_snapped = false;
                }
                if (moving) { moving = false; solve_and_refresh(&a); }
                else if (pending && (a.view != 2)) {
                    Camera cam = active_camera(&a, rw, rh);
                    if (a.tool != UI_TOOL_NONE) app_place(&a, &cam, e.button.x, e.button.y);
                    else { pick_at(&a, &cam, e.button.x, e.button.y); a.focus = -1; }
                }
                pending = false;
                dragging = false;
            }
            else if (e.type == SDL_MOUSEWHEEL) {
                if (a.view == 0) {
                    a.dist *= (e.wheel.y > 0) ? 0.9 : 1.111;
                    atomic_store(&a.restart, 1);
                } else if (a.view == 1) {
                    a.top_zoom = ls_clamp(a.top_zoom * ((e.wheel.y > 0) ? 0.9 : 1.111),
                                          0.1, 20.0);
                    atomic_store(&a.restart, 1);
                }
            }
        }

        if (a.sel_light >= 0 && a.sel_light < a.d.nlights)
            a.new_kind = a.d.lights[a.sel_light].kind;
        a.st.view = a.view;
        a.st.shade_heat = a.shade_heat;
        a.st.tier = (int)a.tier;
        a.st.tool = a.tool;
        a.st.has_selection = (a.sel_light >= 0 || a.sel_prim >= 0);
        a.st.can_undo = (a.undo_n > 0);
        a.st.can_redo = (a.redo_n > 0);
        insp_rebuild(&a);
        ui_apply_state(&a.t, a.st);

        SDL_SetRenderDrawColor(ren, COL_BG.r, COL_BG.g, COL_BG.b, 255);
        SDL_RenderClear(ren);
        draw_toolbar(ren, &a.t, a.st, WIN_H);

        const char *unit = ls_quantity_unit(LS_Q_IRRADIANCE,
            a.st.photometric ? LS_UNITS_PHOTOMETRIC : LS_UNITS_RADIOMETRIC);

        /* ---- title ---- */
        snprintf(buf, sizeof buf, "%s", a.scene_path);
        draw_text(ren, cx, 16, 2,
                  a.view == 0 ? (a.shade_heat ? "3D  ILLUMINANCE" : "3D VIEW")
                              : (a.shade_heat ? "PLAN  ILLUMINANCE" : "TOP VIEW"),
                  COL_INK, 255);
        draw_text_right(ren, WIN_W - 18, 18, 1, buf, COL_MUTED, 255);

        /* One canvas, three views. The plan view and the heat map are framed
         * identically -- same square rect, same orientation, +y up and +x right
         * -- so switching between them compares geometry against result without
         * the picture moving under you. */
        {
            int CW = WIN_W - cx - STATS_W - 30;
            int CH = WIN_H - cy - PLOT_H - 46 - CBAR_H;
            int fit = CW < CH ? CW : CH;
            a.view_x = a.fmap_x = cx;
            a.view_y = a.fmap_y = cy;
            a.view_w = a.view_h = a.fmap_side = fit;
        }

        int side_c = a.view_w;
        {
            pthread_mutex_lock(&a.lock);
            SDL_UpdateTexture(rend_tex, NULL, a.rgb, rw * 3);
            pthread_mutex_unlock(&a.lock);
            SDL_Rect dst = { a.view_x, a.view_y, side_c, side_c };
            SDL_RenderCopy(ren, rend_tex, NULL, &dst);
            draw_rect_line(ren, a.view_x, a.view_y, side_c, side_c, COL_RULE, 255);

            /* Gizmos for objects outside the frustum still project to a
             * coordinate, so the overlay has to be clipped to its own canvas or
             * it draws over the title bar and the panels beside it. */
            SDL_Rect clip = { a.view_x, a.view_y, side_c, side_c };
            SDL_RenderSetClipRect(ren, &clip);
            Camera cam = active_camera(&a, rw, rh);
            /* The drape would be redundant on top of illuminance shading, which
             * already colours every surface by the same quantity. */
            if (a.show_drape && !a.shade_heat) draw_field_drape(ren, &a, &cam);
            draw_overlay(ren, &a, &cam);
            gizmo_draw(ren, &a, &cam);
            SDL_RenderSetClipRect(ren, NULL);

            if (a.shade_heat)
                draw_colorbar(ren, a.view_x, a.view_y + side_c + 10, side_c, CBAR_H,
                              a.heat_lo, a.heat_hi, unit);

            int hint_y = a.view_y + side_c + (a.shade_heat ? CBAR_H + 26 : 12);
            if (a.status[0] && SDL_GetTicks() - a.status_at < 6000) {
                draw_text(ren, a.view_x, hint_y, 1, a.status, COL_ACCENT, 255);
            } else if (a.tool != UI_TOOL_NONE) {
                snprintf(buf, sizeof buf, "CLICK TO PLACE %s   ESC CANCELS",
                         a.tool == UI_TOOL_LIGHT ? "LIGHT" : "PART");
                draw_text(ren, a.view_x, hint_y, 1, buf, COL_ACCENT, 255);
            } else {
                snprintf(buf, sizeof buf, "%d PASSES   %s   %s",
                         atomic_load(&a.passes),
                         a.shade_heat ? "ILLUMINANCE ON EVERY SURFACE" : "RADIANCE",
                         a.view == 1 ? "ORTHOGRAPHIC PLAN" : "DRAG ORBIT");
                draw_text(ren, a.view_x, hint_y, 1, buf, COL_MUTED, 255);
            }
        }

        /* ---- statistics ---- */
        int sx = WIN_W - STATS_W - 12, sy = cy;
        draw_rect_fill(ren, sx, sy, STATS_W, 250, COL_PANEL, 255);
        draw_rect_line(ren, sx, sy, STATS_W, 250, COL_RULE, 255);
        draw_text(ren, sx + 14, sy + 14, 1,
                  a.st.photometric ? "ILLUMINANCE E V" : "IRRADIANCE E E", COL_ACCENT, 255);
        draw_text(ren, sx + 14, sy + 30, 1,
                  a.st.direct_only ? "DIRECT ONLY" : "FULL TRANSPORT", COL_MUTED, 255);
        struct { const char *k; double v; int dec; } rows[] = {
            { "MIN",     a.stats.min,    2 }, { "MAX",   a.stats.max,    2 },
            { "MEAN",    a.stats.mean,   2 }, { "STDDEV",a.stats.stddev, 2 },
            { "U0",      a.stats.u0,     3 }, { "UD",    a.stats.ud,     3 },
            { "CONTRAST",a.stats.michelson, 3 }
        };
        for (int i = 0; i < 7; ++i) {
            int ry = sy + 56 + i * 22;
            draw_text(ren, sx + 14, ry, 1, rows[i].k, COL_MUTED, 255);
            if (rows[i].dec == 3) snprintf(buf, sizeof buf, "%.3f", rows[i].v);
            else if (rows[i].v >= 100.0) snprintf(buf, sizeof buf, "%.0f", rows[i].v);
            else snprintf(buf, sizeof buf, "%.3f", rows[i].v);
            draw_text_right(ren, sx + STATS_W - 14, ry, 2, buf, COL_INK, 255);
        }
        draw_text_right(ren, sx + STATS_W - 14, sy + 224, 1, unit, COL_MUTED, 255);

        insp_draw(ren, &a);

        /* ---- probe readout: only meaningful while illuminance is being shown */
        int py = sy + 208;
        if (!a.shade_heat) py = -1000;
        draw_rect_fill(ren, sx, py, STATS_W, 92, COL_PANEL, 255);
        draw_rect_line(ren, sx, py, STATS_W, 92, COL_RULE, 255);
        draw_text(ren, sx + 14, py + 14, 1, "PROBE", COL_ACCENT, 255);
        if (a.shade_heat && a.have_hover) {
            size_t k = (size_t)a.hover_j * (size_t)rw + (size_t)a.hover_i;
            if (a.hover_i >= 0 && a.hover_i < rw && a.hover_j >= 0 && a.hover_j < rh &&
                a.heat_n[k] > 0) {
                /* The surface point under the cursor, found the same way the
                 * heat pass found it. */
                Camera pcam = active_camera(&a, rw, rh);
                Ray pr = ls_camera_pick_ray(&pcam, a.hover_i + 0.5, a.hover_j + 0.5);
                Hit ph;
                if (ls_scene_intersect(&a.d.scene, &pr, &ph)) {
                    snprintf(buf, sizeof buf, "%.3f %.3f %.3f M", ph.p.x, ph.p.y, ph.p.z);
                    draw_text(ren, sx + 14, py + 34, 1, buf, COL_MUTED, 255);
                }
                double v = a.heat_sum[k] / a.heat_n[k];
                snprintf(buf, sizeof buf, v >= 100 ? "%.0f" : "%.3f", v);
                draw_text(ren, sx + 14, py + 52, 3, buf, COL_INK, 255);
                draw_text_right(ren, sx + STATS_W - 14, py + 58, 1, unit, COL_MUTED, 255);
            } else {
                draw_text(ren, sx + 14, py + 48, 2, "NO SURFACE", COL_MUTED, 160);
            }
        } else if (a.shade_heat) {
            draw_text(ren, sx + 14, py + 48, 2, "HOVER A SURFACE", COL_MUTED, 160);
        }

        /* ---- cross-section ---- */
        int plot_y = WIN_H - PLOT_H - 14;
        if ((a.view == 2) && a.nu > 0) {
            int j = a.have_hover ? a.hover_j : a.nv / 2;
            double *full = malloc((size_t)a.nu * sizeof *full);
            double *dir  = malloc((size_t)a.nu * sizeof *dir);
            for (int i = 0; i < a.nu; ++i) {
                int k = j * a.nu + i;
                full[i] = a.st.direct_only ? a.other[k] : a.disp[k];
                dir[i]  = a.st.direct_only ? a.disp[k]  : a.other[k];
            }
            double fv = (j + 0.5) / a.nv;
            vec3 p = v3add(a.d.grid_o, v3scale(a.d.grid_v, fv));
            snprintf(buf, sizeof buf, "CROSS SECTION AT Y %.3f M", p.y);
            draw_plot(ren, cx, plot_y, WIN_W - cx - 12, PLOT_H,
                      full, dir, a.nu, a.hi, unit, buf);
            free(full); free(dir);
        } else {
            draw_rect_fill(ren, cx, plot_y, WIN_W - cx - 12, PLOT_H, COL_PANEL, 255);
            draw_rect_line(ren, cx, plot_y, WIN_W - cx - 12, PLOT_H, COL_RULE, 255);
            draw_text(ren, cx + 14, plot_y + 14, 1,
                      "CROSS SECTION AVAILABLE IN FIELD MAP MODE", COL_MUTED, 255);
        }

        SDL_RenderPresent(ren);
    }

    atomic_store(&a.quit, true);
    pthread_join(rt, NULL);
    SDL_DestroyTexture(rend_tex);
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
    clear_stack(a.undo, &a.undo_n);
    clear_stack(a.redo, &a.redo_n);
    free(a.heat_sum); free(a.heat_n);
    free(a.drape);
    free(a.g_full); free(a.g_direct); free(a.disp); free(a.other); free(a.rgb);
    ls_film_free(&a.film);
    ls_scene_desc_free(&a.d);
    return 0;
}
