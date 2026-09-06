#include "lightsim/scenefile.h"
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

static Spectrum parse_spd(P *p) {
    char *t = tok(p);
    if (!t) { p->ok = false; return ls_spectrum_const(1.0); }
    if (strcmp(t, "flat") == 0)      return ls_spectrum_const(1.0);
    if (strcmp(t, "blackbody") == 0) return ls_spectrum_blackbody(num(p));
    if (strcmp(t, "daylight") == 0)  return ls_spectrum_daylight(num(p));
    if (strcmp(t, "led") == 0)     { ls_real c = num(p), w = num(p);
                                     return ls_spectrum_gaussian(c, w, 1.0); }
    p->ok = false;
    snprintf(p->d->err, sizeof p->d->err, "unknown spectrum '%s'", t);
    return ls_spectrum_const(1.0);
}

/* Reads "<W|lm> <value> <spd>" and returns radiant flux in watts. Lumens are
 * converted here and never stored, so no photometric value reaches a light. */
static ls_real parse_flux(P *p, Spectrum *spd_out) {
    char *unit = tok(p);
    if (!unit) { p->ok = false; return 0.0; }
    ls_real value = num(p);
    *spd_out = parse_spd(p);
    if (strcmp(unit, "W") == 0)  return value;
    if (strcmp(unit, "lm") == 0) return ls_watts_from_lumens(value, spd_out);
    p->ok = false;
    snprintf(p->d->err, sizeof p->d->err, "flux unit must be W or lm, got '%s'", unit);
    return 0.0;
}

static void add_prim(SceneDesc *d, Prim pr) {
    GROW(d->prims, d->nprims, d->cap_prims, Prim);
    d->prims[d->nprims++] = pr;
}

/* Emissive geometry bound to an area light, so it is visible to the camera and
 * findable by BSDF sampling. */
static void add_emissive_geom(SceneDesc *d, const Light *l, int light_index) {
    Material m;
    memset(&m, 0, sizeof m);
    m.bsdf.kind = LS_BSDF_LAMBERT;
    m.bsdf.rho  = ls_spectrum_zero();
    m.le        = ls_spectrum_scale(l->s_hat, l->radiance);
    m.emissive  = true;
    ensure_mat_capacity(d);
    int mid = d->nmats;
    d->mats[d->nmats] = m;
    snprintf(d->names[d->nmats], 32, "__emit%d", light_index);
    d->nmats++;

    Prim pr;
    memset(&pr, 0, sizeof pr);
    pr.mat_id = mid;
    pr.light_id = light_index;
    switch (l->kind) {
        case LS_LIGHT_RECT:
            pr.kind = LS_PRIM_QUAD; pr.c = l->p; pr.n = l->n;
            pr.ex = l->ex; pr.ey = l->ey; break;
        case LS_LIGHT_DISK:
            pr.kind = LS_PRIM_DISK; pr.c = l->p; pr.n = l->n; pr.r = l->radius; break;
        case LS_LIGHT_SPHERE:
            pr.kind = LS_PRIM_SPHERE; pr.c = l->p; pr.r = l->radius; break;
        default: return;              /* delta lights have no geometry */
    }
    add_prim(d, pr);
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
                if      (which && strcmp(which, "cu") == 0) ls_metal_copper(&m.bsdf.eta, &m.bsdf.kappa);
                else if (which && strcmp(which, "au") == 0) ls_metal_gold(&m.bsdf.eta, &m.bsdf.kappa);
                else                                        ls_metal_aluminium(&m.bsdf.eta, &m.bsdf.kappa);
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
            if (strcmp(kind, "point") == 0) {
                vec3 pos = vec(&p);
                ls_real w = parse_flux(&p, &spd);
                l = ls_light_point(pos, w, spd);
            } else if (strcmp(kind, "spot") == 0) {
                vec3 pos = vec(&p), dir = vec(&p);
                ls_real total = num(&p) * LS_PI / 180.0;
                ls_real fall  = num(&p) * LS_PI / 180.0;
                ls_real w = parse_flux(&p, &spd);
                l = ls_light_spot(pos, dir, total, fall, w, spd);
            } else if (strcmp(kind, "rect") == 0) {
                vec3 c = vec(&p), ex = vec(&p), ey = vec(&p);
                ls_real w = parse_flux(&p, &spd);
                l = ls_light_rect(c, ex, ey, w, spd);
            } else if (strcmp(kind, "disk") == 0) {
                vec3 c = vec(&p), n = vec(&p);
                ls_real r = num(&p);
                ls_real w = parse_flux(&p, &spd);
                l = ls_light_disk(c, n, r, w, spd);
            } else if (strcmp(kind, "sphere") == 0) {
                vec3 c = vec(&p);
                ls_real r = num(&p);
                ls_real w = parse_flux(&p, &spd);
                l = ls_light_sphere(c, r, w, spd);
            } else if (strcmp(kind, "sun") == 0) {
                vec3 dir = vec(&p);
                ls_real e = num(&p);
                spd = parse_spd(&p);
                l = ls_light_directional(dir, e, spd);
            } else {
                p.ok = false;
                snprintf(d->err, sizeof d->err, "unknown light kind '%s'", kind);
                break;
            }
            if (!p.ok) break;
            GROW(d->lights, d->nlights, d->cap_lights, Light);
            int idx = d->nlights;
            ls_light_finalize(&l, idx);
            d->lights[d->nlights++] = l;
            add_emissive_geom(d, &d->lights[idx], idx);
        } else {
            p.ok = false;
            snprintf(d->err, sizeof d->err, "unknown directive '%s'", t);
        }
    }

    free(buf);
    if (!p.ok) return false;
    d->scene.prims = d->prims;   d->scene.nprims = d->nprims;
    d->scene.mats  = d->mats;    d->scene.nmats  = d->nmats;
    d->scene.lights = d->lights; d->scene.nlights = d->nlights;
    return true;
}

void ls_scene_desc_free(SceneDesc *d) {
    free(d->prims); free(d->mats); free(d->names); free(d->lights);
    memset(d, 0, sizeof *d);
}
