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

#include "lightsim/scenefile.h"
#include "lightsim/integrator.h"
#include "lightsim/analysis.h"
#include "lightsim/units.h"
#include "lightsim/film.h"
#include "lightsim/thread.h"
#include "lightsim/export.h"
#include "ui.h"
#include "draw.h"
#include "font.h"

#define WIN_W 1360
#define WIN_H 860
#define PLOT_H 168
#define STATS_W 268
#define CBAR_H 14

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

    const char *scene_path;
} App;

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

static Camera orbit_camera(const App *a, int w, int h) {
    double ce = cos(a->el), se = sin(a->el);
    vec3 eye = v3add(a->target, v3(a->dist * ce * cos(a->az),
                                   a->dist * ce * sin(a->az),
                                   a->dist * se));
    return ls_camera_look_at(eye, a->target, v3(0, 0, 1),
                             a->d.has_camera ? a->d.camera.fov_y * 180.0 / LS_PI : 42.0,
                             w, h);
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
        if (a->st.grid_mode) { SDL_Delay(40); continue; }
        if (atomic_exchange(&a->restart, 0)) {
            size_t np = (size_t)a->film.width * (size_t)a->film.height;
            memset(a->film.pix, 0, np * sizeof *a->film.pix);
            memset(a->film.n,   0, np * sizeof *a->film.n);
            atomic_store(&a->passes, 0);
            a->exposure = 0.0;
        }
        RJob job = { a, orbit_camera(a, a->film.width, a->film.height), 4 };
        ls_parallel_for(a->film.height, 0, render_row, &job);
        atomic_fetch_add(&a->passes, 1);
        tonemap(a);
    }
    return NULL;
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
        } else {
            int sx, sy;
            if (project_view(a, cam, p->c, &sx, &sy))
                draw_marker(ren, sx, sy, 8, COL_ACCENT, 255);
        }
    }
}

static void save_outputs(App *a) {
    if (a->st.grid_mode) {
        int up = a->nu < 64 ? (64 + a->nu - 1) / a->nu : 1;
        if (ls_write_falsecolor_ppm("out/view_field.ppm", a->disp, a->nu, a->nv,
                                    a->lo, a->hi, up))
            printf("wrote out/view_field.ppm\n");
    } else {
        if (ls_film_write_ppm(&a->film, "out/view_render.ppm", a->exposure))
            printf("wrote out/view_render.ppm\n");
    }
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
    if (argc < 2) {
        fprintf(stderr, "usage: %s <scene.txt> [--render|--field] [--fine]\n", argv[0]);
        return 1;
    }
    bool want_render = false, want_fine = false;
    for (int i = 2; i < argc; ++i) {
        if      (!strcmp(argv[i], "--render")) want_render = true;
        else if (!strcmp(argv[i], "--field"))  want_render = false;
        else if (!strcmp(argv[i], "--fine"))   want_fine = true;
    }
    App a;
    memset(&a, 0, sizeof a);
    a.scene_path = argv[1];
    if (!ls_scene_load(&a.d, argv[1])) {
        fprintf(stderr, "scene error: %s\n", a.d.err);
        return 1;
    }
    a.st.has_grid = a.d.has_grid;
    a.st.has_camera = a.d.has_camera;
    a.st.grid_mode = a.d.has_grid && !want_render;
    if (!a.d.has_camera) a.st.grid_mode = true;
    a.st.photometric = true;
    a.st.direct_only = false;
    a.st.high_quality = want_fine;
    a.hover_i = a.hover_j = -1;
    a.sel_light = a.sel_prim = -1;

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

    ui_init(&a.t);

    /* Canvas geometry. */
    int cx = UI_TOOLBAR_W + 18, cy = 44;
    int cw = WIN_W - cx - STATS_W - 30;
    int ch = WIN_H - cy - PLOT_H - 46 - CBAR_H;
    int side = cw < ch ? cw : ch;

    int rw = 560, rh = (int)(side > 0 ? side * 0.75 : 420);
    ls_film_init(&a.film, rw, rh);
    a.rgb = calloc((size_t)rw * (size_t)rh * 3u, 1);
    pthread_mutex_init(&a.lock, NULL);
    atomic_store(&a.restart, 1);
    pthread_t rt;
    pthread_create(&rt, NULL, render_thread, &a);

    a.view_x = cx; a.view_y = cy;
    a.view_w = side; a.view_h = side * rh / rw;

    SDL_Texture *grid_tex = SDL_CreateTexture(ren, SDL_PIXELFORMAT_RGB24,
        SDL_TEXTUREACCESS_STREAMING, a.nu > 0 ? a.nu : 1, a.nv > 0 ? a.nv : 1);
    SDL_Texture *rend_tex = SDL_CreateTexture(ren, SDL_PIXELFORMAT_RGB24,
        SDL_TEXTUREACCESS_STREAMING, rw, rh);
    SDL_SetTextureScaleMode(grid_tex, SDL_ScaleModeNearest);
    SDL_SetTextureScaleMode(rend_tex, SDL_ScaleModeLinear);

    bool dragging = false, pending = false;
    int last_x = 0, last_y = 0, press_x = 0, press_y = 0;
    const int DRAG_THRESHOLD = 4;
    bool running = true;
    char buf[128];

    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) running = false;
            else if (e.type == SDL_KEYDOWN) {
                switch (e.key.keysym.sym) {
                    case SDLK_ESCAPE: running = false; break;
                    case SDLK_1: if (a.st.has_grid) a.st.grid_mode = true; break;
                    case SDLK_2: if (a.st.has_camera) { a.st.grid_mode = false;
                                     atomic_store(&a.restart, 1); } break;
                    case SDLK_u: a.st.photometric = !a.st.photometric; refresh_display(&a); break;
                    case SDLK_t: a.st.direct_only = !a.st.direct_only; refresh_display(&a); break;
                    case SDLK_q: a.st.high_quality = !a.st.high_quality;
                                 if (a.st.grid_mode) { solve_grid(&a); refresh_display(&a); }
                                 else atomic_store(&a.restart, 1);
                                 break;
                    case SDLK_r: if (a.st.grid_mode) { solve_grid(&a); refresh_display(&a); } break;
                    case SDLK_s: save_outputs(&a); break;
                    case SDLK_b: export_blender(&a); break;
                    default: break;
                }
            }
            else if (e.type == SDL_MOUSEMOTION) {
                a.t.hover = ui_hit_test(&a.t, e.motion.x, e.motion.y);
                if (pending && (abs(e.motion.x - press_x) > DRAG_THRESHOLD ||
                                abs(e.motion.y - press_y) > DRAG_THRESHOLD)) {
                    pending = false;
                    dragging = true;
                }
                if (dragging && !a.st.grid_mode) {
                    a.az -= (e.motion.x - last_x) * 0.008;
                    a.el = ls_clamp(a.el + (e.motion.y - last_y) * 0.008, -1.45, 1.45);
                    atomic_store(&a.restart, 1);
                }
                last_x = e.motion.x; last_y = e.motion.y;
                a.have_hover = false;
                if (a.st.grid_mode && a.nu > 0 && side > 0) {
                    int gx = e.motion.x - cx, gy = e.motion.y - cy;
                    if (gx >= 0 && gx < side && gy >= 0 && gy < side) {
                        a.hover_i = gx * a.nu / side;
                        a.hover_j = a.nv - 1 - (gy * a.nv / side);
                        a.have_hover = true;
                    }
                }
            }
            else if (e.type == SDL_MOUSEBUTTONDOWN) {
                int hit = ui_hit_test(&a.t, e.button.x, e.button.y);
                if (hit >= 0) {
                    if (a.t.buttons[hit].enabled) switch (a.t.buttons[hit].action) {
                        case UI_MODE_GRID:   a.st.grid_mode = true; break;
                        case UI_MODE_RENDER: a.st.grid_mode = false;
                                             atomic_store(&a.restart, 1); break;
                        case UI_UNITS:       a.st.photometric = !a.st.photometric;
                                             refresh_display(&a); break;
                        case UI_TRANSPORT:   a.st.direct_only = !a.st.direct_only;
                                             refresh_display(&a); break;
                        case UI_QUALITY:     a.st.high_quality = !a.st.high_quality;
                                             if (a.st.grid_mode) { solve_grid(&a); refresh_display(&a); }
                                             else atomic_store(&a.restart, 1);
                                             break;
                        case UI_SOLVE:       solve_grid(&a); refresh_display(&a); break;
                        case UI_SAVE:        save_outputs(&a); break;
                        case UI_BLENDER:     export_blender(&a); break;
                        default: break;
                    }
                } else if (!a.st.grid_mode) {
                    pending = true;
                    press_x = e.button.x; press_y = e.button.y;
                }
            }
            else if (e.type == SDL_MOUSEBUTTONUP) {
                if (pending && !a.st.grid_mode) {
                    /* A press that never moved: select whatever is under it. */
                    Camera cam = orbit_camera(&a, rw, rh);
                    pick_at(&a, &cam, e.button.x, e.button.y);
                }
                pending = false;
                dragging = false;
            }
            else if (e.type == SDL_MOUSEWHEEL && !a.st.grid_mode) {
                a.dist *= (e.wheel.y > 0) ? 0.9 : 1.111;
                atomic_store(&a.restart, 1);
            }
        }

        ui_apply_state(&a.t, a.st);

        SDL_SetRenderDrawColor(ren, COL_BG.r, COL_BG.g, COL_BG.b, 255);
        SDL_RenderClear(ren);
        draw_toolbar(ren, &a.t, a.st, WIN_H);

        const char *unit = ls_quantity_unit(LS_Q_IRRADIANCE,
            a.st.photometric ? LS_UNITS_PHOTOMETRIC : LS_UNITS_RADIOMETRIC);

        /* ---- title ---- */
        snprintf(buf, sizeof buf, "%s", a.scene_path);
        draw_text(ren, cx, 16, 2, a.st.grid_mode ? "FIELD MAP" : "PROGRESSIVE RENDER",
                  COL_INK, 255);
        draw_text_right(ren, WIN_W - 18, 18, 1, buf, COL_MUTED, 255);

        if (a.st.grid_mode && a.nu > 0) {
            /* false-colour field */
            unsigned char *px = malloc((size_t)a.nu * (size_t)a.nv * 3u);
            double span = a.hi - a.lo;
            for (int j = 0; j < a.nv; ++j)
                for (int i = 0; i < a.nu; ++i) {
                    int src = (a.nv - 1 - j) * a.nu + i;   /* +y up on screen */
                    Uint8 r, g, b;
                    draw_viridis((a.disp[src] - a.lo) / (span > 0 ? span : 1), &r, &g, &b);
                    size_t o = ((size_t)j * (size_t)a.nu + (size_t)i) * 3u;
                    px[o] = r; px[o+1] = g; px[o+2] = b;
                }
            SDL_UpdateTexture(grid_tex, NULL, px, a.nu * 3);
            free(px);
            SDL_Rect dst = { cx, cy, side, side };
            SDL_RenderCopy(ren, grid_tex, NULL, &dst);
            draw_rect_line(ren, cx, cy, side, side, COL_RULE, 255);

            if (a.have_hover) {
                int hx = cx + a.hover_i * side / a.nu;
                int hy = cy + (a.nv - 1 - a.hover_j) * side / a.nv;
                draw_rect_line(ren, hx, cy, side / a.nu > 2 ? side / a.nu : 2, side,
                               COL_ACCENT, 70);
                draw_rect_line(ren, cx, hy, side, side / a.nv > 2 ? side / a.nv : 2,
                               COL_ACCENT, 140);
            }
            draw_colorbar(ren, cx, cy + side + 12, side, CBAR_H, a.lo, a.hi, unit);
        } else if (!a.st.grid_mode) {
            pthread_mutex_lock(&a.lock);
            SDL_UpdateTexture(rend_tex, NULL, a.rgb, rw * 3);
            pthread_mutex_unlock(&a.lock);
            int dw = a.view_w, dh = a.view_h;
            SDL_Rect dst = { a.view_x, a.view_y, dw, dh };
            SDL_RenderCopy(ren, rend_tex, NULL, &dst);
            draw_rect_line(ren, a.view_x, a.view_y, dw, dh, COL_RULE, 255);

            Camera cam = orbit_camera(&a, rw, rh);
            draw_overlay(ren, &a, &cam);

            snprintf(buf, sizeof buf, "%d PASSES  %d SPP  CLICK TO SELECT  DRAG TO ORBIT",
                     atomic_load(&a.passes), atomic_load(&a.passes) * 4);
            draw_text(ren, a.view_x, a.view_y + dh + 12, 1, buf, COL_MUTED, 255);
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

        /* ---- probe readout ---- */
        int py = sy + 262;
        draw_rect_fill(ren, sx, py, STATS_W, 92, COL_PANEL, 255);
        draw_rect_line(ren, sx, py, STATS_W, 92, COL_RULE, 255);
        draw_text(ren, sx + 14, py + 14, 1,
                  a.st.grid_mode ? "PROBE" : "SELECTION", COL_ACCENT, 255);
        if (a.have_hover) {
            int k = a.hover_j * a.nu + a.hover_i;
            double fu = (a.hover_i + 0.5) / a.nu, fv = (a.hover_j + 0.5) / a.nv;
            vec3 p = v3add(a.d.grid_o, v3add(v3scale(a.d.grid_u, fu),
                                             v3scale(a.d.grid_v, fv)));
            snprintf(buf, sizeof buf, "%.3f %.3f %.3f M", p.x, p.y, p.z);
            draw_text(ren, sx + 14, py + 34, 1, buf, COL_MUTED, 255);
            snprintf(buf, sizeof buf, a.disp[k] >= 100 ? "%.0f" : "%.3f", a.disp[k]);
            draw_text(ren, sx + 14, py + 52, 3, buf, COL_INK, 255);
            draw_text_right(ren, sx + STATS_W - 14, py + 58, 1, unit, COL_MUTED, 255);
        } else if (!a.st.grid_mode && a.sel_light >= 0) {
            /* A first cut at the inspector: what the selected source is, in the
             * terms the SIMPLE tier will use. */
            const Light *l = &a.d.lights[a.sel_light];
            static const char *kind_name[] = { "POINT", "SUN", "SPOT",
                                               "SPHERE", "DISK", "RECT" };
            snprintf(buf, sizeof buf, "LIGHT %d  %s", a.sel_light,
                     kind_name[l->kind]);
            draw_text(ren, sx + 14, py + 32, 1, buf, COL_MUTED, 255);
            Spectrum phi = ls_spectrum_scale(l->s_hat, l->phi_e);
            snprintf(buf, sizeof buf, "%.0f", ls_photometric(&phi));
            draw_text(ren, sx + 14, py + 50, 3, buf, COL_INK, 255);
            draw_text_right(ren, sx + STATS_W - 14, py + 58, 1, "LM", COL_MUTED, 255);
        } else if (!a.st.grid_mode && a.sel_prim >= 0) {
            snprintf(buf, sizeof buf, "SURFACE %d", a.sel_prim);
            draw_text(ren, sx + 14, py + 32, 1, buf, COL_MUTED, 255);
            const Prim *pr = &a.d.prims[a.sel_prim];
            snprintf(buf, sizeof buf, "%.2f",
                     ls_spectrum_mean(&a.d.mats[pr->mat_id].bsdf.rho));
            draw_text(ren, sx + 14, py + 50, 3, buf, COL_INK, 255);
            draw_text_right(ren, sx + STATS_W - 14, py + 58, 1, "ALBEDO", COL_MUTED, 255);
        } else {
            draw_text(ren, sx + 14, py + 48, 2,
                      a.st.grid_mode ? "HOVER THE MAP" : "CLICK A LIGHT", COL_MUTED, 160);
        }

        /* ---- cross-section ---- */
        int plot_y = WIN_H - PLOT_H - 14;
        if (a.st.grid_mode && a.nu > 0) {
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
    SDL_DestroyTexture(grid_tex);
    SDL_DestroyTexture(rend_tex);
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
    free(a.g_full); free(a.g_direct); free(a.disp); free(a.other); free(a.rgb);
    ls_film_free(&a.film);
    ls_scene_desc_free(&a.d);
    return 0;
}
