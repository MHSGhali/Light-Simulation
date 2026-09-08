#include "lightsim/import.h"
#include "lightsim/sceneedit.h"
#include "lightsim/color.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------- paths ---- */

void ls_dir_of(const char *path, char *out, size_t n) {
    const char *slash = strrchr(path, '/');
    if (!slash) { snprintf(out, n, "%s", ""); return; }
    size_t len = (size_t)(slash - path);
    if (len >= n) len = n - 1;
    memcpy(out, path, len);
    out[len] = '\0';
}

void ls_resolve_path(const char *base_dir, const char *rel, char *out, size_t n) {
    if (rel[0] == '/' || !base_dir || !base_dir[0]) {
        snprintf(out, n, "%s", rel);
        return;
    }
    snprintf(out, n, "%s/%s", base_dir, rel);
}

/* -------------------------------------------------------------- lexing -- */

/* The whole file in one buffer, tokenised in place, exactly as ls_scene_load
 * does it. No getline, no per-line allocation. */
typedef struct {
    char *buf;
    size_t len;
} Slurp;

static bool slurp(const char *path, Slurp *s, char err[256]) {
    FILE *f = fopen(path, "rb");
    if (!f) { snprintf(err, 256, "cannot open '%s'", path); return false; }
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (n < 0) { fclose(f); snprintf(err, 256, "cannot size '%s'", path); return false; }
    s->buf = malloc((size_t)n + 1);
    if (!s->buf) { fclose(f); snprintf(err, 256, "out of memory"); return false; }
    s->len = fread(s->buf, 1, (size_t)n, f);
    s->buf[s->len] = '\0';
    fclose(f);
    return true;
}

/* Line first, then tokens within it.
 *
 * OBJ is line-oriented -- `f 1 2 3` and `v 1 2 3` differ only by keyword, and a
 * face ends where its line ends. Tokenising the whole buffer as one stream and
 * hoping to notice newlines is how the first version of this read the `f` of
 * the NEXT line as another vertex of the current face. Splitting into lines
 * first makes that unrepresentable. */
static char *next_line(char **cur, char *end) {
    if (*cur >= end) return NULL;
    char *s = *cur, *p = s;
    while (p < end && *p != '\n') p++;
    if (p < end) { *p = '\0'; *cur = p + 1; } else { *cur = end; }
    size_t len = strlen(s);
    if (len && s[len - 1] == '\r') s[len - 1] = '\0';
    return s;
}

/* Next space-delimited token within one line, or NULL at its end. */
static char *word(char **p) {
    while (**p == ' ' || **p == '\t') (*p)++;
    if (**p == '\0') return NULL;
    char *s = *p;
    while (**p && **p != ' ' && **p != '\t') (*p)++;
    if (**p) { **p = '\0'; (*p)++; }
    return s;
}

/* ----------------------------------------------------------- obj faces -- */

typedef struct {
    vec3 *v;   int nv, cap_v;
    int  *idx; int ntri, cap_tri;
} Build;

static bool push_vert(Build *b, vec3 p) {
    if (b->nv >= b->cap_v) {
        int cap = b->cap_v ? b->cap_v * 2 : 256;
        vec3 *nv = realloc(b->v, (size_t)cap * sizeof *nv);
        if (!nv) return false;
        b->v = nv; b->cap_v = cap;
    }
    b->v[b->nv++] = p;
    return true;
}

static bool push_tri(Build *b, int a, int c, int d) {
    if (b->ntri >= b->cap_tri) {
        int cap = b->cap_tri ? b->cap_tri * 2 : 256;
        int *ni = realloc(b->idx, (size_t)cap * 3 * sizeof *ni);
        if (!ni) return false;
        b->idx = ni; b->cap_tri = cap;
    }
    b->idx[b->ntri*3+0] = a;
    b->idx[b->ntri*3+1] = c;
    b->idx[b->ntri*3+2] = d;
    b->ntri++;
    return true;
}

/* An OBJ face vertex is "v", "v/vt", "v//vn" or "v/vt/vn"; only the position
 * index is used. Indices are 1-based, and negative means relative to the end.
 * Returns a 0-based index, or -1 if it is out of range -- an unchecked index is
 * the commonest way a malformed OBJ walks a reader off its own heap. */
static int face_index(const char *tok, int nverts) {
    char *slash = strchr(tok, '/');
    long i = strtol(tok, NULL, 10);
    (void)slash;
    if (i > 0) i -= 1;
    else if (i < 0) i += nverts;
    else return -1;
    if (i < 0 || i >= nverts) return -1;
    return (int)i;
}

/* ------------------------------------------------------------- reading -- */

/* Read one group's triangles. `want` is the usemtl name to keep, or NULL for
 * everything. Positions are file-global (OBJ indices are), so all `v` records
 * are read whatever the group. */
static bool read_obj(const char *path, const char *want, Build *b,
                     char *mtllib, size_t mtllib_n, char err[256]) {
    Slurp s;
    if (!slurp(path, &s, err)) return false;
    char *cur = s.buf, *end = s.buf + s.len, *line;
    char active[64] = "";
    bool ok = true;
    int lineno = 0;

    while (ok && (line = next_line(&cur, end)) != NULL) {
        lineno++;
        char *lp = line;
        char *kw = word(&lp);
        if (!kw || kw[0] == '#') continue;

        if (!strcmp(kw, "v")) {
            char *x = word(&lp), *y = word(&lp), *z = word(&lp);
            if (!x || !y || !z) {
                snprintf(err, 256, "%s:%d: vertex needs three numbers", path, lineno);
                ok = false; break;
            }
            if (!push_vert(b, v3(atof(x), atof(y), atof(z)))) {
                snprintf(err, 256, "out of memory"); ok = false; break;
            }
        } else if (!strcmp(kw, "usemtl")) {
            char *n = word(&lp);
            snprintf(active, sizeof active, "%s", n ? n : "");
        } else if (!strcmp(kw, "mtllib")) {
            char *n = word(&lp);
            if (n && mtllib) snprintf(mtllib, mtllib_n, "%s", n);
        } else if (!strcmp(kw, "f")) {
            if (want && strcmp(active, want) != 0) continue;
            /* Fan-triangulate. Correct for a convex face; Blender's exporter
             * with "Triangulated Mesh" ticked makes the question moot. */
            int first = -1, prev = -1, made = 0;
            char *t;
            while ((t = word(&lp)) != NULL) {
                int vi = face_index(t, b->nv);
                if (vi < 0) {
                    snprintf(err, 256, "%s:%d: face index '%s' out of range "
                             "(%d vertices so far)", path, lineno, t, b->nv);
                    ok = false; break;
                }
                if (first < 0) first = vi;
                else if (prev < 0) prev = vi;
                else {
                    if (!push_tri(b, first, prev, vi)) {
                        snprintf(err, 256, "out of memory"); ok = false; break;
                    }
                    prev = vi;
                    made++;
                }
            }
            if (!ok) break;
            if (made == 0) {
                snprintf(err, 256, "%s:%d: face has fewer than three vertices",
                         path, lineno);
                ok = false; break;
            }
        }
        /* vn, vt, s, g, o and anything else are read and ignored; see
         * import.h for why vn in particular is not honoured. */
    }
    free(s.buf);
    return ok;
}

Mesh *ls_obj_load(const char *path, const char *group, char err[256]) {
    Build b;
    memset(&b, 0, sizeof b);
    char mtllib[256] = "";
    if (!read_obj(path, group, &b, mtllib, sizeof mtllib, err)) {
        free(b.v); free(b.idx);
        return NULL;
    }
    if (b.ntri == 0) {
        snprintf(err, 256, "'%s' has no faces%s%s", path,
                 group ? " in group " : "", group ? group : "");
        free(b.v); free(b.idx);
        return NULL;
    }
    Mesh *m = ls_mesh_build(b.v, b.nv, b.idx, b.ntri);
    free(b.v); free(b.idx);
    if (!m) { snprintf(err, 256, "'%s' built no usable geometry", path); return NULL; }
    /* Store the ABSOLUTE path. The scene writer emits whatever is here, and a
     * scene saved to a different directory than it was loaded from must still
     * find its geometry -- which a path relative to the old working directory
     * would not. A hand-authored scene may still use a relative path; it is
     * resolved against the scene's own directory on the way in, and comes back
     * out absolute. */
    char abs[512];
    if (realpath(path, abs)) snprintf(m->src_path, sizeof m->src_path, "%s", abs);
    else                     snprintf(m->src_path, sizeof m->src_path, "%s", path);
    snprintf(m->group, sizeof m->group, "%s", group ? group : "");
    return m;
}

/* --------------------------------------------------------------- mtl ---- */

/* Kd of each material in an .mtl, added to the scene. Names collide with the
 * scene's own materials often enough (everyone has a "floor") that a policy is
 * needed: an identical material is reused, a different one is renamed, and
 * nothing is ever silently overwritten. */
static int mtl_material(SceneDesc *d, const char *name, ls_real kd[3],
                        bool emissive_seen, char err[256]) {
    if (emissive_seen) {
        snprintf(err, 256, "material '%s' is emissive (Ke): a mesh emitter is "
                 "not sampled as a light and would read low. Use `light rect`.",
                 name);
        return -1;
    }
    Material m;
    memset(&m, 0, sizeof m);
    m.bsdf.kind = LS_BSDF_LAMBERT;
    m.rgb[0] = kd[0]; m.rgb[1] = kd[1]; m.rgb[2] = kd[2];
    m.from_rgb = true;
    RGB c = { kd[0], kd[1], kd[2] };
    m.bsdf.rho = ls_spectrum_from_rgb_reflectance(c);

    char want[32];
    snprintf(want, sizeof want, "%s", name);
    for (char *p = want; *p; ++p)
        if (*p == ' ' || *p == '\t' || *p == '#') *p = '_';

    for (int attempt = 0; attempt < 1000; ++attempt) {
        char tryname[32];
        if (attempt == 0) snprintf(tryname, sizeof tryname, "%s", want);
        else              snprintf(tryname, sizeof tryname, "%.26s.%03d", want, attempt);

        int existing = -1;
        for (int i = 0; i < d->nmats; ++i)
            if (!strcmp(d->names[i], tryname)) { existing = i; break; }
        if (existing < 0) return ls_scene_add_material(d, m, tryname);

        /* Same name, same material: legitimate dedup. Same name, different
         * material: rename rather than bind this mesh to a stranger's albedo
         * or overwrite a surface the user authored. */
        const Material *e = &d->mats[existing];
        if (e->from_rgb && e->rgb[0] == m.rgb[0] &&
            e->rgb[1] == m.rgb[1] && e->rgb[2] == m.rgb[2]) return existing;
    }
    snprintf(err, 256, "too many materials named '%s'", name);
    return -1;
}

/* ------------------------------------------------------------- scene ---- */

int ls_scene_import_obj(SceneDesc *d, const char *path) {
    /* Pass one: the group names and the mtllib, so each group can be loaded
     * into its own Mesh. */
    Slurp s;
    if (!slurp(path, &s, d->err)) return -1;

    char groups[64][64];
    int ngroups = 0;
    char mtllib[256] = "";
    char *cur = s.buf, *end = s.buf + s.len, *line;
    bool any_face_before_usemtl = false;
    while ((line = next_line(&cur, end)) != NULL) {
        char *lp = line;
        char *kw = word(&lp);
        if (!kw) continue;
        if (!strcmp(kw, "usemtl")) {
            char *n = word(&lp);
            if (n && ngroups < 64) {
                bool seen = false;
                for (int i = 0; i < ngroups; ++i)
                    if (!strcmp(groups[i], n)) { seen = true; break; }
                if (!seen) snprintf(groups[ngroups++], 64, "%s", n);
            }
        } else if (!strcmp(kw, "mtllib")) {
            char *n = word(&lp);
            if (n) snprintf(mtllib, sizeof mtllib, "%s", n);
        } else if (!strcmp(kw, "f") && ngroups == 0) {
            any_face_before_usemtl = true;
        }
    }
    free(s.buf);

    char dir[512];
    ls_dir_of(path, dir, sizeof dir);

    /* Kd per material name, from the sidecar .mtl if there is one. */
    char mtl_names[64][64];
    ls_real mtl_kd[64][3];
    bool mtl_emit[64];
    int nmtl = 0;
    if (mtllib[0]) {
        char mpath[1024];
        ls_resolve_path(dir, mtllib, mpath, sizeof mpath);
        Slurp ms;
        char ignored[256];
        if (slurp(mpath, &ms, ignored)) {
            char *mc = ms.buf, *me = ms.buf + ms.len, *mline;
            int curm = -1;
            while ((mline = next_line(&mc, me)) != NULL) {
                char *lp = mline;
                char *kw = word(&lp);
                if (!kw) continue;
                if (!strcmp(kw, "newmtl") && nmtl < 64) {
                    char *n = word(&lp);
                    curm = nmtl++;
                    snprintf(mtl_names[curm], 64, "%s", n ? n : "");
                    mtl_kd[curm][0] = mtl_kd[curm][1] = mtl_kd[curm][2] = 0.5;
                    mtl_emit[curm] = false;
                } else if (!strcmp(kw, "Kd") && curm >= 0) {
                    char *r = word(&lp), *g = word(&lp), *b = word(&lp);
                    if (r && g && b) {
                        mtl_kd[curm][0] = atof(r);
                        mtl_kd[curm][1] = atof(g);
                        mtl_kd[curm][2] = atof(b);
                    }
                } else if (!strcmp(kw, "Ke") && curm >= 0) {
                    char *r = word(&lp), *g = word(&lp), *b = word(&lp);
                    if (r && g && b && (atof(r) > 0 || atof(g) > 0 || atof(b) > 0))
                        mtl_emit[curm] = true;
                }
            }
            free(ms.buf);
        }
    }

    /* A file with no usemtl is one group named "-". */
    if (ngroups == 0 && any_face_before_usemtl) snprintf(groups[ngroups++], 64, "-");
    if (ngroups == 0) {
        snprintf(d->err, sizeof d->err, "'%s' has no faces", path);
        return -1;
    }

    int added = 0;
    for (int g = 0; g < ngroups; ++g) {
        const char *gname = strcmp(groups[g], "-") ? groups[g] : NULL;

        ls_real kd[3] = { 0.5, 0.5, 0.5 };
        bool emissive = false;
        for (int i = 0; i < nmtl; ++i)
            if (gname && !strcmp(mtl_names[i], gname)) {
                kd[0] = mtl_kd[i][0]; kd[1] = mtl_kd[i][1]; kd[2] = mtl_kd[i][2];
                emissive = mtl_emit[i];
                break;
            }

        int mat = mtl_material(d, gname ? gname : "imported", kd, emissive, d->err);
        if (mat < 0) return -1;

        Mesh *m = ls_obj_load(path, gname, d->err);
        if (!m) return -1;
        if (ls_scene_add_mesh_prim(d, m, mat) < 0) {
            ls_mesh_release(m);
            snprintf(d->err, sizeof d->err, "out of memory adding '%s'", path);
            return -1;
        }
        added++;
    }
    return added;
}
