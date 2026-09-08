#include "lightsim/scenefile.h"
#include "lightsim/sceneedit.h"
#include "lightsim/units.h"
#include "lightsim/bsdf.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { char *cur; const char *end; SceneDesc *d; bool ok; } P;

static char *tok(P *p) {
    if (!p->ok) return NULL;
    while (p->cur < p->end && (*p->cur == ' ' || *p->cur == '\t' ||
                               *p->cur == '\n' || *p->cur == '\r')) p->cur++;
    if (p->cur >= p->end) return NULL;
    if (*p->cur == '#') {                       /* comment to end of line */
        while (p->cur < p->end && *p->cur != '\n') p->cur++;
        return tok(p);
    }
    char *s = p->cur;
    while (p->cur < p->end && *p->cur != ' ' && *p->cur != '\t' &&
           *p->cur != '\n' && *p->cur != '\r') p->cur++;
    if (p->cur < p->end) *p->cur++ = '\0';
    return s;
}

static ls_real num(P *p) {
    char *t = tok(p);
    if (!t) { p->ok = false; snprintf(p->d->err, sizeof p->d->err,
                                      "unexpected end of file, expected a number"); return 0; }
    return atof(t);
}
static vec3 vec(P *p) { ls_real x = num(p), y = num(p), z = num(p); return v3(x, y, z); }

#define GROW(arr, n, cap, type)                                            \
    do { if ((n) >= (cap)) { (cap) = (cap) ? (cap) * 2 : 8;                 \
         (arr) = realloc((arr), (size_t)(cap) * sizeof(type)); } } while (0)

/* Materials and their names are parallel arrays sharing one capacity, so they
 * must grow together. Growing them with two separate GROW calls does not work:
 * the first bumps the capacity, and the second then sees n < cap and skips the
 * allocation entirely, leaving the name array short. */
static void ensure_mat_capacity(SceneDesc *d) {
    if (d->nmats < d->cap_mats) return;
    d->cap_mats = d->cap_mats ? d->cap_mats * 2 : 8;
    d->mats  = realloc(d->mats,  (size_t)d->cap_mats * sizeof *d->mats);
    d->names = realloc(d->names, (size_t)d->cap_mats * sizeof *d->names);
}

static int find_mat(SceneDesc *d, const char *name) {
    for (int i = 0; i < d->nmats; ++i)
        if (strcmp(d->names[i], name) == 0) return i;
    return -1;
}

/* The spectrum, plus a record of how it was spelled, so ls_scene_save can
 * write the scene back in the form it was authored rather than as 95 bins. */
static Spectrum parse_spd(P *p, LsSpdKind *kind, ls_real *a, ls_real *b) {
    *kind = LS_SPD_FLAT; *a = 0.0; *b = 0.0;
    char *t = tok(p);
    if (!t) { p->ok = false; return ls_spectrum_const(1.0); }
    Spectrum s;
    if (strcmp(t, "flat") == 0) {
        s = ls_spectrum_const(1.0);
    } else if (strcmp(t, "blackbody") == 0) {
        *kind = LS_SPD_BLACKBODY; *a = num(p);
        s = ls_spectrum_blackbody(*a);
    } else if (strcmp(t, "daylight") == 0) {
        *kind = LS_SPD_DAYLIGHT; *a = num(p);
        s = ls_spectrum_daylight(*a);
    } else if (strcmp(t, "led") == 0) {
        *kind = LS_SPD_LED; *a = num(p); *b = num(p);
        s = ls_spectrum_gaussian(*a, *b, 1.0);
    } else {
        p->ok = false;
        snprintf(p->d->err, sizeof p->d->err, "unknown spectrum '%s'", t);
        return ls_spectrum_const(1.0);
    }

    /* A shape with no power in the sampled band cannot be normalised to unit
     * integral, and a light built from it trips the s_hat invariant inside
     * ls_light_finalize -- an abort, with the offending number nowhere in
     * sight. Refuse it here, where there is a directive to name. Reachable two
     * ways: a Planckian below about 130 K, whose every visible bin underflows
     * the float bins to zero, and an LED lobe with a non-positive width or a
     * centre so far outside the band that even its tail underflows. */
    if (!(ls_spectrum_integrate(&s) > 0.0)) {
        p->ok = false;
        if (*kind == LS_SPD_LED)
            snprintf(p->d->err, sizeof p->d->err,
                     "spectrum 'led %g %g' has no power in the sampled band", *a, *b);
        else
            snprintf(p->d->err, sizeof p->d->err,
                     "spectrum '%s %g' has no power in the sampled band", t, *a);
        return ls_spectrum_const(1.0);
    }
    return s;
}

/* Reads "<W|lm> <value> <spd>" and returns radiant flux in watts. Lumens are
 * converted here and never stored, so no photometric value reaches a light. */
typedef struct {
    LsSpdKind kind;
    ls_real   a, b;
    bool      in_lumens;
    ls_real   authored;      /* the number given, in the unit given */
} Authored;

static ls_real parse_flux(P *p, Spectrum *spd_out, Authored *rec) {
    char *unit = tok(p);
    if (!unit) { p->ok = false; return 0.0; }
    ls_real value = num(p);
    *spd_out = parse_spd(p, &rec->kind, &rec->a, &rec->b);
    rec->authored = value;
    if (strcmp(unit, "W") == 0)  { rec->in_lumens = false; return value; }
    if (strcmp(unit, "lm") == 0) { rec->in_lumens = true;
                                   return ls_watts_from_lumens(value, spd_out); }
    p->ok = false;
    snprintf(p->d->err, sizeof p->d->err, "flux unit must be W or lm, got '%s'", unit);
    return 0.0;
}

static void add_prim(SceneDesc *d, Prim pr) {
    GROW(d->prims, d->nprims, d->cap_prims, Prim);
    d->prims[d->nprims++] = pr;
}

bool ls_scene_load(SceneDesc *d, const char *path) {
    memset(d, 0, sizeof *d);
    FILE *fp = fopen(path, "rb");
    if (!fp) { snprintf(d->err, sizeof d->err, "cannot open '%s'", path); return false; }
    fseek(fp, 0, SEEK_END);
    long len = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    char *buf = malloc((size_t)len + 1);
    if (!buf) { fclose(fp); return false; }
    size_t got = fread(buf, 1, (size_t)len, fp);
    buf[got] = '\0';
    fclose(fp);

    P p = { buf, buf + got, d, true };
    char *t;
    while (p.ok && (t = tok(&p)) != NULL) {
        if (strcmp(t, "camera") == 0) {
            vec3 eye = vec(&p), target = vec(&p);
            ls_real fov = num(&p);
            int w = (int)num(&p), h = (int)num(&p);
            d->camera = ls_camera_look_at(eye, target, v3(0, 0, 1), fov, w, h);
            d->cam_eye = eye; d->cam_target = target; d->cam_fov_deg = fov;
            d->has_camera = true;
        } else if (strcmp(t, "grid") == 0) {
            d->grid_o = vec(&p); d->grid_u = vec(&p); d->grid_v = vec(&p);
            d->grid_nu = (int)num(&p); d->grid_nv = (int)num(&p);
            d->has_grid = true;
        } else if (strcmp(t, "material") == 0) {
            char *name = tok(&p);
            char *kind = tok(&p);
            if (!name || !kind) { p.ok = false; break; }
            Material m;
            memset(&m, 0, sizeof m);
            if (strcmp(kind, "lambert") == 0) {
                m.bsdf.kind = LS_BSDF_LAMBERT;
                m.bsdf.rho = ls_spectrum_const(num(&p));
            } else if (strcmp(kind, "metal") == 0) {
                char *which = tok(&p);
                m.bsdf.kind = LS_BSDF_CONDUCTOR;
                if      (which && strcmp(which, "cu") == 0) { ls_metal_copper(&m.bsdf.eta, &m.bsdf.kappa);
                                                              snprintf(m.metal, sizeof m.metal, "cu"); }
                else if (which && strcmp(which, "au") == 0) { ls_metal_gold(&m.bsdf.eta, &m.bsdf.kappa);
                                                              snprintf(m.metal, sizeof m.metal, "au"); }
                else                                        { ls_metal_aluminium(&m.bsdf.eta, &m.bsdf.kappa);
                                                              snprintf(m.metal, sizeof m.metal, "al"); }
                m.bsdf.alpha = num(&p);
            } else if (strcmp(kind, "emit") == 0) {
                m.bsdf.kind = LS_BSDF_LAMBERT;
                m.bsdf.rho = ls_spectrum_zero();
                m.le = ls_spectrum_const(num(&p));
                m.emissive = true;
            } else {
                p.ok = false;
                snprintf(d->err, sizeof d->err, "unknown material kind '%s'", kind);
                break;
            }
            ensure_mat_capacity(d);
            d->mats[d->nmats] = m;
            snprintf(d->names[d->nmats], 32, "%s", name);
            d->nmats++;
        } else if (strcmp(t, "plane") == 0 || strcmp(t, "quad") == 0 ||
                   strcmp(t, "sphere") == 0) {
            char *mname = tok(&p);
            int mid = mname ? find_mat(d, mname) : -1;
            if (mid < 0) { p.ok = false;
                snprintf(d->err, sizeof d->err, "unknown material '%s'",
                         mname ? mname : "(none)"); break; }
            Prim pr;
            memset(&pr, 0, sizeof pr);
            pr.mat_id = mid;
            pr.light_id = -1;
            if (strcmp(t, "plane") == 0) {
                pr.kind = LS_PRIM_PLANE; pr.c = vec(&p); pr.n = v3norm(vec(&p));
            } else if (strcmp(t, "quad") == 0) {
                pr.kind = LS_PRIM_QUAD; pr.c = vec(&p); pr.n = v3norm(vec(&p));
                pr.ex = vec(&p); pr.ey = vec(&p);
            } else {
                pr.kind = LS_PRIM_SPHERE; pr.c = vec(&p); pr.r = num(&p);
            }
            add_prim(d, pr);
        } else if (strcmp(t, "light") == 0) {
            char *kind = tok(&p);
            if (!kind) { p.ok = false; break; }
            Light l;
            Spectrum spd;
            Authored rec;
            memset(&rec, 0, sizeof rec);
            if (strcmp(kind, "point") == 0) {
                vec3 pos = vec(&p);
                ls_real w = parse_flux(&p, &spd, &rec);
                l = ls_light_point(pos, w, spd);
            } else if (strcmp(kind, "spot") == 0) {
                vec3 pos = vec(&p), dir = vec(&p);
                ls_real total = num(&p) * LS_PI / 180.0;
                ls_real fall  = num(&p) * LS_PI / 180.0;
                ls_real w = parse_flux(&p, &spd, &rec);
                l = ls_light_spot(pos, dir, total, fall, w, spd);
            } else if (strcmp(kind, "rect") == 0) {
                vec3 c = vec(&p), ex = vec(&p), ey = vec(&p);
                ls_real w = parse_flux(&p, &spd, &rec);
                l = ls_light_rect(c, ex, ey, w, spd);
            } else if (strcmp(kind, "disk") == 0) {
                vec3 c = vec(&p), n = vec(&p);
                ls_real r = num(&p);
                ls_real w = parse_flux(&p, &spd, &rec);
                l = ls_light_disk(c, n, r, w, spd);
            } else if (strcmp(kind, "sphere") == 0) {
                vec3 c = vec(&p);
                ls_real r = num(&p);
                ls_real w = parse_flux(&p, &spd, &rec);
                l = ls_light_sphere(c, r, w, spd);
            } else if (strcmp(kind, "sun") == 0) {
                vec3 dir = vec(&p);
                ls_real e = num(&p);
                spd = parse_spd(&p, &rec.kind, &rec.a, &rec.b);
                rec.authored = e; rec.in_lumens = false;
                l = ls_light_directional(dir, e, spd);
            } else {
                p.ok = false;
                snprintf(d->err, sizeof d->err, "unknown light kind '%s'", kind);
                break;
            }
            if (!p.ok) break;
            l.spd_kind = rec.kind;
            l.spd_a = rec.a;
            l.spd_b = rec.b;
            l.flux_in_lumens = rec.in_lumens;
            l.flux_authored = rec.authored;
            GROW(d->lights, d->nlights, d->cap_lights, Light);
            int idx = d->nlights;
            ls_light_finalize(&l, idx);
            d->lights[d->nlights++] = l;
            /* One implementation of the light/geometry pairing, shared with
             * the editor -- see ls_scene_sync_light_geom. */
            ls_scene_rebuild(d);
            ls_scene_sync_light_geom(d, idx);
        } else {
            p.ok = false;
            snprintf(d->err, sizeof d->err, "unknown directive '%s'", t);
        }
    }

    free(buf);
    if (!p.ok) return false;
    ls_scene_rebuild(d);
    return true;
}

void ls_scene_desc_free(SceneDesc *d) {
    for (int i = 0; i < d->nmeshes; ++i) ls_mesh_release(d->meshes[i]);
    free(d->meshes);
    free(d->prims); free(d->mats); free(d->names); free(d->lights);
    memset(d, 0, sizeof *d);
}

/* ------------------------------------------------------------------ writer */

static void wvec(FILE *f, vec3 v) {
    fprintf(f, " %.6g %.6g %.6g", v.x, v.y, v.z);
}

/* The spectrum clause, spelled the way it was authored. */
static void wspd(FILE *f, const Light *l) {
    switch (l->spd_kind) {
        case LS_SPD_BLACKBODY: fprintf(f, " blackbody %.6g", l->spd_a); break;
        case LS_SPD_DAYLIGHT:  fprintf(f, " daylight %.6g", l->spd_a);  break;
        case LS_SPD_LED:       fprintf(f, " led %.6g %.6g", l->spd_a, l->spd_b); break;
        case LS_SPD_FLAT:
        default:               fprintf(f, " flat"); break;
    }
}

/* The flux clause, in the unit it was given in. A light created in the editor
 * rather than parsed has flux_authored == 0, so fall back to its watts. */
static void wflux(FILE *f, const Light *l) {
    if (l->flux_authored > 0.0)
        fprintf(f, " %s %.6g", l->flux_in_lumens ? "lm" : "W", l->flux_authored);
    else
        fprintf(f, " W %.8g", l->kind == LS_LIGHT_DIRECTIONAL ? l->e_perp : l->phi_e);
}

bool ls_scene_save(const SceneDesc *d, const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) return false;

    fprintf(f, "# lightsim scene\n\n");

    for (int i = 0; i < d->nmats; ++i) {
        const Material *m = &d->mats[i];
        /* Emissive materials generated for an area light's geometry are
         * recreated on load; writing them would duplicate the light. */
        if (strncmp(d->names[i], "__emit", 6) == 0) continue;
        if (m->emissive)
            fprintf(f, "material %s emit %.6g\n", d->names[i], ls_spectrum_mean(&m->le));
        else if (m->bsdf.kind == LS_BSDF_CONDUCTOR)
            fprintf(f, "material %s metal %s %.6g\n", d->names[i],
                    m->metal[0] ? m->metal : "al", m->bsdf.alpha);
        else
            fprintf(f, "material %s lambert %.6g\n", d->names[i],
                    ls_spectrum_mean(&m->bsdf.rho));
    }

    fprintf(f, "\n");
    for (int i = 0; i < d->nprims; ++i) {
        const Prim *p = &d->prims[i];
        if (p->light_id >= 0) continue;          /* an area light's own geometry */
        if (p->mat_id < 0 || p->mat_id >= d->nmats) continue;
        const char *mat = d->names[p->mat_id];
        switch (p->kind) {
            case LS_PRIM_PLANE:
                fprintf(f, "plane %s", mat);  wvec(f, p->c); wvec(f, p->n);
                fprintf(f, "\n"); break;
            case LS_PRIM_QUAD:
                fprintf(f, "quad %s", mat);   wvec(f, p->c); wvec(f, p->n);
                wvec(f, p->ex); wvec(f, p->ey); fprintf(f, "\n"); break;
            case LS_PRIM_SPHERE:
                fprintf(f, "sphere %s", mat); wvec(f, p->c);
                fprintf(f, " %.6g\n", p->r); break;
            case LS_PRIM_DISK:
                /* Not a parseable primitive directive today; skip rather than
                 * write something the loader would reject. */
                break;
            case LS_PRIM_MESH:
                /* The `mesh` directive lands with the OBJ importer; until then
                 * a mesh has no textual form, so skip it rather than write
                 * something that would not load. */
                break;
        }
    }

    fprintf(f, "\n");
    for (int i = 0; i < d->nlights; ++i) {
        const Light *l = &d->lights[i];
        switch (l->kind) {
            case LS_LIGHT_POINT:
                fprintf(f, "light point");  wvec(f, l->p); wflux(f, l); break;
            case LS_LIGHT_SPOT:
                fprintf(f, "light spot");   wvec(f, l->p); wvec(f, l->n);
                fprintf(f, " %.6g %.6g",
                        acos(ls_clamp(l->cos_total, -1.0, 1.0)) * 180.0 / LS_PI,
                        acos(ls_clamp(l->cos_falloff, -1.0, 1.0)) * 180.0 / LS_PI);
                wflux(f, l); break;
            case LS_LIGHT_RECT:
                fprintf(f, "light rect");   wvec(f, l->p); wvec(f, l->ex); wvec(f, l->ey);
                wflux(f, l); break;
            case LS_LIGHT_DISK:
                fprintf(f, "light disk");   wvec(f, l->p); wvec(f, l->n);
                fprintf(f, " %.6g", l->radius); wflux(f, l); break;
            case LS_LIGHT_SPHERE:
                fprintf(f, "light sphere"); wvec(f, l->p);
                fprintf(f, " %.6g", l->radius); wflux(f, l); break;
            case LS_LIGHT_DIRECTIONAL:
                fprintf(f, "light sun");    wvec(f, l->n);
                fprintf(f, " %.6g", l->e_perp); break;
        }
        wspd(f, l);
        fprintf(f, "\n");
    }

    if (d->has_grid) {
        fprintf(f, "\ngrid");
        wvec(f, d->grid_o); wvec(f, d->grid_u); wvec(f, d->grid_v);
        fprintf(f, " %d %d\n", d->grid_nu, d->grid_nv);
    }
    if (d->has_camera) {
        fprintf(f, "camera");
        wvec(f, d->cam_eye); wvec(f, d->cam_target);
        fprintf(f, " %.6g %d %d\n", d->cam_fov_deg, d->camera.width, d->camera.height);
    }

    fclose(f);
    return true;
}
