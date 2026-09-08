#include "inspect.h"
#include "lightsim/units.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

static const char *const KIND_NAMES[] = { "POINT", "SUN", "SPOT", "SPHERE", "DISK", "RECT" };
static const char *const SPD_NAMES[]  = { "FLAT", "BLACKBODY", "DAYLIGHT", "LED" };
static const char *const MAT_NAMES[]  = { "LAMBERT", "METAL" };
static const char *const METAL_NAMES[]= { "AL", "CU", "AU" };
static const char *const PRIM_NAMES[] = { "SPHERE", "PLANE", "DISK", "QUAD", "MESH" };

static void push(Field *out, int max, int *n, Field f) {
    if (*n < max) out[(*n)++] = f;
}
static Field val(FieldId id, const char *label, const char *unit,
                 double v, double lo, double hi) {
    Field f; memset(&f, 0, sizeof f);
    f.id = id; f.label = label; f.unit = unit; f.value = v; f.lo = lo; f.hi = hi;
    return f;
}
static Field ro(FieldId id, const char *label, const char *unit, double v) {
    Field f = val(id, label, unit, v, 0, 0);
    f.readonly = true;
    return f;
}
static Field en(FieldId id, const char *label, double v,
                const char *const *names, int n) {
    Field f = val(id, label, "", v, 0, n - 1);
    f.is_enum = true; f.names = names; f.nnames = n;
    return f;
}
static Field head(const char *label) {
    Field f; memset(&f, 0, sizeof f);
    f.label = label; f.heading = true; f.readonly = true;
    return f;
}

/* Luminous flux of a light, from the one place that converts. */
static double light_lumens(const Light *l) {
    Spectrum phi = ls_spectrum_scale(l->s_hat, l->phi_e);
    return ls_photometric(&phi);
}

static bool light_is_area(const Light *l) {
    return l->kind == LS_LIGHT_RECT || l->kind == LS_LIGHT_DISK
        || l->kind == LS_LIGHT_SPHERE;
}

/* Illuminance a directional source delivers on a perpendicular plane, in lux.
 * The stored quantity is irradiance in W/m^2; this is the same conversion the
 * flux rows do, through the same efficacy. */
static double sun_lux(const Light *l) {
    Spectrum e = ls_spectrum_scale(l->s_hat, l->e_perp);
    return ls_photometric(&e);
}

int ls_inspect_fields(const SceneDesc *d, int sel_light, int sel_prim,
                      LsUnitSystem units, LsTier tier, Field *out, int max) {
    int n = 0;
    const bool photo = (units == LS_UNITS_PHOTOMETRIC);

    if (sel_light >= 0 && sel_light < d->nlights) {
        const Light *l = &d->lights[sel_light];
        const bool adv = (tier >= LS_TIER_ADVANCED);
        const bool sci = (tier >= LS_TIER_SCIENTIFIC);

        push(out, max, &n, head("SOURCE"));
        if (adv) push(out, max, &n, en(FLD_L_KIND, "TYPE", (double)l->kind,
                                       KIND_NAMES, 6));

        /* Flux, in whichever system the viewer is set to. Lumens is the number
         * on the box and watts is the number in the integral; both are the same
         * stored radiant flux, so either row edits it. The scientific tier gets
         * both at once with the efficacy that relates them. */
        Field f_lm = val(FLD_L_LM, "FLUX",        "LM", light_lumens(l), 0.0, 1e7);
        Field f_w  = val(FLD_L_W,  "RADIANT FLUX", "W", l->phi_e,        0.0, 1e5);
        Field f_lx = val(FLD_L_EV, "ILLUMINANCE", "LX", sun_lux(l),      0.0, 1e6);
        Field f_e  = val(FLD_L_W,  "IRRADIANCE", "W/M2", l->e_perp,      0.0, 1e5);

        if (l->kind != LS_LIGHT_DIRECTIONAL) {
            push(out, max, &n, photo ? f_lm : f_w);
            if (sci) {
                push(out, max, &n, photo ? f_w : f_lm);
                double lm = light_lumens(l);
                push(out, max, &n, ro(FLD_L_EFFICACY, "EFFICACY", "LM/W",
                                      l->phi_e > 0 ? lm / l->phi_e : 0.0));
            }
        } else {
            /* A directional source is defined by what it delivers on a plane
             * facing it, not by a total flux -- it has no finite extent. */
            push(out, max, &n, photo ? f_lx : f_e);
            if (sci) push(out, max, &n, photo ? f_e : f_lx);
        }

        if (l->kind == LS_LIGHT_SPOT) {
            double k = l->beam_k > 0.0 ? l->beam_k
                     : ls_light_k_from_beam_angle(
                           2.0 * acos(ls_clamp(l->cos_total, -1.0, 1.0)) * 180.0 / LS_PI);
            push(out, max, &n, val(FLD_L_BEAM, "BEAM ANGLE", "DEG",
                                   ls_light_beam_angle_from_k(k), 2.0, 178.0));
            push(out, max, &n, ro(FLD_L_FIELD, "FIELD AT 10%", "DEG",
                                  ls_light_field_angle_from_k(k)));
            if (sci) {
                push(out, max, &n, ro(FLD_L_K, "BEAM EXPONENT K", "", k));
                push(out, max, &n, ro(FLD_L_OMEGA, "SOLID ANGLE", "SR",
                                      l->omega_eff));
            }
        }

        /* Colour. SIMPLE gets one number; ADVANCED gets the model behind it. */
        if (adv) push(out, max, &n, en(FLD_L_SPD, "SPECTRUM", (double)l->spd_kind,
                                       SPD_NAMES, 4));
        if (l->spd_kind == LS_SPD_LED && adv) {
            push(out, max, &n, val(FLD_L_LED_C, "PEAK", "NM", l->spd_a, 380.0, 780.0));
            push(out, max, &n, val(FLD_L_LED_W, "WIDTH FWHM", "NM", l->spd_b, 1.0, 300.0));
        } else if (l->spd_kind != LS_SPD_LED) {
            push(out, max, &n, val(FLD_L_CCT, "COLOUR TEMP", "K",
                                   l->spd_a > 0 ? l->spd_a : 4000.0, 1200.0, 20000.0));
        }

        push(out, max, &n, head("PLACEMENT"));
        push(out, max, &n, val(FLD_L_X, "X", "M", l->p.x, -50, 50));
        push(out, max, &n, val(FLD_L_Y, "Y", "M", l->p.y, -50, 50));
        push(out, max, &n, val(FLD_L_Z, "Z", "M", l->p.z, -50, 50));
        if (adv && l->kind != LS_LIGHT_POINT && l->kind != LS_LIGHT_SPHERE) {
            push(out, max, &n, val(FLD_L_AIMX, "AIM X", "", l->n.x, -1, 1));
            push(out, max, &n, val(FLD_L_AIMY, "AIM Y", "", l->n.y, -1, 1));
            push(out, max, &n, val(FLD_L_AIMZ, "AIM Z", "", l->n.z, -1, 1));
        }

        if (adv && light_is_area(l)) {
            push(out, max, &n, head("EMITTER"));
            if (l->kind == LS_LIGHT_RECT) {
                push(out, max, &n, val(FLD_L_SIZEU, "WIDTH", "M",
                                       2.0 * v3len(l->ex), 1e-4, 20.0));
                push(out, max, &n, val(FLD_L_SIZEV, "HEIGHT", "M",
                                       2.0 * v3len(l->ey), 1e-4, 20.0));
            } else {
                push(out, max, &n, val(FLD_L_RADIUS, "RADIUS", "M",
                                       l->radius, 1e-4, 10.0));
            }
            if (sci) {
                push(out, max, &n, ro(FLD_L_RADIANCE, "RADIANCE", "W/M2/SR",
                                      l->radiance));
                Spectrum i0 = ls_spectrum_scale(l->s_hat,
                                                ls_light_intensity(l, l->n));
                push(out, max, &n, ro(FLD_L_I0_CD, "INTENSITY", "CD",
                                      ls_photometric(&i0)));
                push(out, max, &n, ro(FLD_L_I0_W, "INTENSITY", "W/SR",
                                      ls_light_intensity(l, l->n)));
            }
        } else if (sci && l->kind != LS_LIGHT_DIRECTIONAL) {
            Spectrum i0 = ls_spectrum_scale(l->s_hat, ls_light_intensity(l, l->n));
            push(out, max, &n, ro(FLD_L_I0_CD, "INTENSITY", "CD", ls_photometric(&i0)));
            push(out, max, &n, ro(FLD_L_I0_W, "INTENSITY", "W/SR",
                                  ls_light_intensity(l, l->n)));
        }
        return n;
    }

    if (sel_prim >= 0 && sel_prim < d->nprims) {
        const Prim *p = &d->prims[sel_prim];
        const Material *m = &d->mats[p->mat_id];
        const bool adv = (tier >= LS_TIER_ADVANCED);

        push(out, max, &n, head("PART"));
        push(out, max, &n, ro(FLD_P_KIND, "SHAPE", "", (double)p->kind));
        push(out, max, &n, val(FLD_P_X, "X", "M", p->c.x, -50, 50));
        push(out, max, &n, val(FLD_P_Y, "Y", "M", p->c.y, -50, 50));
        push(out, max, &n, val(FLD_P_Z, "Z", "M", p->c.z, -50, 50));
        if (p->kind == LS_PRIM_SPHERE || p->kind == LS_PRIM_DISK)
            push(out, max, &n, val(FLD_P_RADIUS, "RADIUS", "M", p->r, 1e-4, 20.0));
        if (p->kind == LS_PRIM_QUAD) {
            push(out, max, &n, val(FLD_P_SIZEU, "WIDTH", "M",
                                   2.0 * v3len(p->ex), 1e-4, 100.0));
            push(out, max, &n, val(FLD_P_SIZEV, "HEIGHT", "M",
                                   2.0 * v3len(p->ey), 1e-4, 100.0));
        }
        if (p->kind == LS_PRIM_MESH && p->mesh_id >= 0 &&
            p->mesh_id < d->nmeshes && d->meshes[p->mesh_id]) {
            /* Imported geometry is read-only here: its size and shape came from
             * the file, and the editable placement is the X/Y/Z above plus the
             * gizmo. Showing the count is what tells you the import worked. */
            const Mesh *msh = d->meshes[p->mesh_id];   /* `m` is the Material */
            push(out, max, &n, ro(FLD_P_TRIS, "TRIANGLES", "", (double)msh->ntris));
            if (adv) push(out, max, &n, ro(FLD_P_AREA, "SURFACE AREA", "M2", msh->area));
        }

        push(out, max, &n, head("MATERIAL"));
        int mk = (m->bsdf.kind == LS_BSDF_CONDUCTOR) ? 1 : 0;
        if (adv) push(out, max, &n, en(FLD_M_KIND, "TYPE", mk, MAT_NAMES, 2));
        if (mk == 0) {
            push(out, max, &n, val(FLD_M_ALBEDO, "REFLECTANCE", "",
                                   ls_spectrum_mean(&m->bsdf.rho), 0.0, 1.0));
        } else {
            int mi = (m->metal[0] == 'c') ? 1 : (m->metal[0] == 'a' && m->metal[1] == 'u') ? 2 : 0;
            if (adv) push(out, max, &n, en(FLD_M_METAL, "METAL", mi, METAL_NAMES, 3));
            push(out, max, &n, val(FLD_M_ROUGH, "ROUGHNESS", "",
                                   m->bsdf.alpha, 0.0, 1.0));
        }
        return n;
    }
    return 0;
}

void ls_inspect_format(const Field *f, char *buf, size_t n) {
    if (f->heading) { snprintf(buf, n, "%s", f->label); return; }
    if (f->is_enum) {
        int i = (int)(f->value + 0.5);
        if (i < 0) i = 0;
        if (i >= f->nnames) i = f->nnames - 1;
        snprintf(buf, n, "%s", f->names[i]);
        return;
    }
    if (f->id == FLD_P_KIND) {
        int i = (int)(f->value + 0.5);
        snprintf(buf, n, "%s", (i >= 0 && i < 5) ? PRIM_NAMES[i] : "?");
        return;
    }
    double v = fabs(f->value);
    if (v >= 1000.0)      snprintf(buf, n, "%.0f", f->value);
    else if (v >= 10.0)   snprintf(buf, n, "%.1f", f->value);
    else if (v >= 0.1)    snprintf(buf, n, "%.3f", f->value);
    else                  snprintf(buf, n, "%.4f", f->value);
}

/* Rebuild s_hat from the authoring record, and keep the AUTHORED flux fixed --
 * a user who asked for 850 lm still wants 850 lm after changing the colour
 * temperature, even though that changes the watts required. */
static void respec(Light *l) {
    Spectrum s;
    switch (l->spd_kind) {
        case LS_SPD_BLACKBODY: s = ls_spectrum_blackbody(l->spd_a); break;
        case LS_SPD_DAYLIGHT:  s = ls_spectrum_daylight(l->spd_a);  break;
        case LS_SPD_LED:       s = ls_spectrum_gaussian(l->spd_a, l->spd_b, 1.0); break;
        case LS_SPD_FLAT:
        default:               s = ls_spectrum_const(1.0); break;
    }
    /* Last line of defence for the s_hat invariant. A shape with no power in
     * the sampled band cannot be normalised to unit integral, and installing it
     * aborts inside ls_light_finalize -- an assert nowhere near the row that
     * was clicked. Every combination the inspector can author is re-homed by
     * the setters above; if one ever slips through, keep the spectrum the light
     * already had rather than taking the process down. */
    if (!(ls_spectrum_integrate(&s) > 0.0)) return;
    l->s_hat = ls_spectrum_normalize_to(s, 1.0);
    if (l->flux_in_lumens && l->flux_authored > 0.0)
        l->phi_e = ls_watts_from_lumens(l->flux_authored, &l->s_hat);
}

static void set_aim(Light *l, int axis, double v) {
    double a[3] = { l->n.x, l->n.y, l->n.z };
    a[axis] = v;
    vec3 d = v3(a[0], a[1], a[2]);
    if (v3len2(d) < 1e-12) d = v3(0, 0, -1);
    l->n = v3norm(d);
}

bool ls_inspect_set(SceneDesc *d, int sel_light, int sel_prim,
                    FieldId id, double v) {
    if (sel_light >= 0 && sel_light < d->nlights) {
        Light *l = &d->lights[sel_light];
        switch (id) {
            case FLD_L_KIND: {
                LightKind k = (LightKind)(int)(v + 0.5);
                if (k == l->kind) return false;
                l->kind = k;
                /* Give a newly-area light a sane extent rather than zero, which
                 * would make its radiance infinite. */
                if (l->radius <= 0.0) l->radius = 0.05;
                if (v3len2(l->ex) <= 0.0) l->ex = v3(0.05, 0, 0);
                if (v3len2(l->ey) <= 0.0) l->ey = v3(0, -0.05, 0);
                if (k == LS_LIGHT_RECT) l->n = v3norm(v3cross(l->ex, l->ey));
                if (k == LS_LIGHT_SPOT && l->beam_k <= 0.0)
                    l->beam_k = ls_light_k_from_beam_angle(40.0);
                break;
            }
            case FLD_L_LM:
                l->flux_in_lumens = true;
                l->flux_authored = v;
                l->phi_e = ls_watts_from_lumens(v, &l->s_hat);
                break;
            case FLD_L_EV:
                /* Lux in, W/m^2 stored. Same conversion as lumens -> watts:
                 * units.h owns it, so no photometric value is ever kept. */
                l->e_perp = ls_watts_from_lumens(v, &l->s_hat);
                break;
            case FLD_L_W:
                if (l->kind == LS_LIGHT_DIRECTIONAL) { l->e_perp = v; }
                else {
                    l->flux_in_lumens = false;
                    l->flux_authored = v;
                    l->phi_e = v;
                }
                break;
            case FLD_L_BEAM: l->beam_k = ls_light_k_from_beam_angle(v); break;
            case FLD_L_CCT:
                l->spd_a = v;
                /* SIMPLE never meets the distinction: a Planckian radiator below
                 * about 4000 K, a CIE daylight illuminant above, which is the
                 * usual convention for "warm white" versus "daylight". */
                if (l->spd_kind != LS_SPD_BLACKBODY && l->spd_kind != LS_SPD_DAYLIGHT)
                    l->spd_kind = LS_SPD_DAYLIGHT;
                if (l->spd_kind == LS_SPD_DAYLIGHT && v < 4000.0) l->spd_kind = LS_SPD_BLACKBODY;
                if (l->spd_kind == LS_SPD_BLACKBODY && v >= 4000.0) l->spd_kind = LS_SPD_DAYLIGHT;
                respec(l);
                break;
            case FLD_L_SPD:
                l->spd_kind = (LsSpdKind)(int)(v + 0.5);
                /* spd_a means a different quantity per kind -- a peak
                 * wavelength in nm for an LED, a colour temperature in K for
                 * the others -- so changing the kind has to re-home it, and the
                 * width separately. Gating the wavelength on "the width is
                 * unset", as this once did, only re-homed it the FIRST time a
                 * light became an LED: cycle DAYLIGHT -> LED -> ... -> LED and
                 * the second pass left 4000 in a nanometre field, which is a
                 * Gaussian with no power in the band and no way to normalise. */
                if (l->spd_kind == LS_SPD_LED) {
                    if (!(l->spd_a >= 380.0 && l->spd_a <= 780.0)) l->spd_a = 550.0;
                    if (!(l->spd_b > 0.0))                         l->spd_b = 40.0;
                } else if (l->spd_kind == LS_SPD_BLACKBODY ||
                           l->spd_kind == LS_SPD_DAYLIGHT) {
                    if (!(l->spd_a >= 1200.0 && l->spd_a <= 20000.0)) l->spd_a = 4000.0;
                }
                respec(l);
                break;
            case FLD_L_LED_C: l->spd_a = v; respec(l); break;
            case FLD_L_LED_W: l->spd_b = v; respec(l); break;
            case FLD_L_X: l->p.x = v; break;
            case FLD_L_Y: l->p.y = v; break;
            case FLD_L_Z: l->p.z = v; break;
            case FLD_L_AIMX: set_aim(l, 0, v); break;
            case FLD_L_AIMY: set_aim(l, 1, v); break;
            case FLD_L_AIMZ: set_aim(l, 2, v); break;
            case FLD_L_SIZEU: l->ex = v3scale(v3norm(l->ex), v * 0.5); break;
            case FLD_L_SIZEV: l->ey = v3scale(v3norm(l->ey), v * 0.5); break;
            case FLD_L_RADIUS: l->radius = v; break;
            default: return false;
        }
        if (l->kind == LS_LIGHT_RECT) l->n = v3norm(v3cross(l->ex, l->ey));
        ls_scene_update_light(d, sel_light);
        return true;
    }

    if (sel_prim >= 0 && sel_prim < d->nprims) {
        Prim *p = &d->prims[sel_prim];
        Material *m = &d->mats[p->mat_id];
        switch (id) {
            case FLD_P_X: p->c.x = v; break;
            case FLD_P_Y: p->c.y = v; break;
            case FLD_P_Z: p->c.z = v; break;
            case FLD_P_RADIUS: p->r = v; break;
            case FLD_P_SIZEU: p->ex = v3scale(v3norm(p->ex), v * 0.5); break;
            case FLD_P_SIZEV: p->ey = v3scale(v3norm(p->ey), v * 0.5); break;
            case FLD_M_KIND: {
                int k = (int)(v + 0.5);
                if (k == 1) {
                    m->bsdf.kind = LS_BSDF_CONDUCTOR;
                    if (m->bsdf.alpha <= 0.0) m->bsdf.alpha = 0.1;
                    if (!m->metal[0]) snprintf(m->metal, sizeof m->metal, "al");
                    ls_metal_aluminium(&m->bsdf.eta, &m->bsdf.kappa);
                } else {
                    m->bsdf.kind = LS_BSDF_LAMBERT;
                    if (ls_spectrum_mean(&m->bsdf.rho) <= 0.0)
                        m->bsdf.rho = ls_spectrum_const(0.5);
                }
                break;
            }
            case FLD_M_ALBEDO: m->bsdf.rho = ls_spectrum_const(v); break;
            case FLD_M_ROUGH:  m->bsdf.alpha = v; break;
            case FLD_M_METAL: {
                int i = (int)(v + 0.5);
                const char *nm = (i == 1) ? "cu" : (i == 2) ? "au" : "al";
                snprintf(m->metal, sizeof m->metal, "%s", nm);
                if (i == 1)      ls_metal_copper(&m->bsdf.eta, &m->bsdf.kappa);
                else if (i == 2) ls_metal_gold(&m->bsdf.eta, &m->bsdf.kappa);
                else             ls_metal_aluminium(&m->bsdf.eta, &m->bsdf.kappa);
                break;
            }
            default: return false;
        }
        /* A part that is an area light's own geometry must not drift from the
         * light; push the change back through the light instead. */
        if (p->light_id >= 0) {
            d->lights[p->light_id].p = p->c;
            ls_scene_update_light(d, p->light_id);
        }
        ls_scene_rebuild(d);
        return true;
    }
    return false;
}
