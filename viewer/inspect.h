/* inspect.h — the property inspector's model.
 *
 * Free of SDL so the field lists and the edit semantics can be exercised
 * headlessly; viewer/main.c only draws what this produces and feeds edits back.
 *
 * TWO INDEPENDENT AXES. Conflating them was a real bug: a viewer switched to
 * radiometric still authored its lights in lumens, because the only way to
 * reach watts was to raise the tier.
 *
 * The TIER decides how much is revealed. It never changes what is stored or
 * simulated:
 *
 *   SIMPLE      what is printed on an off-the-shelf luminaire: flux, beam
 *               angle, colour temperature, where it is.
 *   ADVANCED    the source itself: emitter kind and size, aim, the spectral
 *               model and its parameters, surface reflectivity.
 *   SCIENTIFIC  everything derived: both flux systems side by side, efficacy,
 *               radiant and luminous intensity, radiance, the beam exponent
 *               and its solid angle.
 *
 * The UNIT SYSTEM decides which of the two equivalent numbers leads, at every
 * tier: lumens and lux when photometric, watts and W/m^2 when radiometric.
 * Editing either one is editing the same stored radiant flux -- units.h keeps
 * the photometric value unreachable outside its own layer, so what a light
 * holds is always watts, whichever row you typed into.
 */
#ifndef LIGHTSIM_VIEWER_INSPECT_H
#define LIGHTSIM_VIEWER_INSPECT_H

#include "lightsim/sceneedit.h"
#include "lightsim/units.h"

typedef enum { LS_TIER_SIMPLE, LS_TIER_ADVANCED, LS_TIER_SCIENTIFIC } LsTier;

typedef enum {
    FLD_NONE = 0,
    /* ---- light ---- */
    FLD_L_KIND, FLD_L_LM, FLD_L_W, FLD_L_EV,
    FLD_L_BEAM, FLD_L_FIELD,
    FLD_L_CCT, FLD_L_SPD, FLD_L_LED_C, FLD_L_LED_W,
    FLD_L_X, FLD_L_Y, FLD_L_Z,
    FLD_L_AIMX, FLD_L_AIMY, FLD_L_AIMZ,
    FLD_L_SIZEU, FLD_L_SIZEV, FLD_L_RADIUS,
    FLD_L_EFFICACY, FLD_L_I0_CD, FLD_L_I0_W, FLD_L_RADIANCE,
    FLD_L_K, FLD_L_OMEGA,
    /* ---- part ---- */
    FLD_P_KIND, FLD_P_X, FLD_P_Y, FLD_P_Z,
    FLD_P_RADIUS, FLD_P_SIZEU, FLD_P_SIZEV,
    FLD_P_TRIS, FLD_P_AREA,
    FLD_M_KIND, FLD_M_ALBEDO, FLD_M_ROUGH, FLD_M_METAL,
    FLD_COUNT
} FieldId;

typedef struct {
    FieldId     id;
    const char *label;
    const char *unit;        /* "" when the label carries it */
    double      value;
    double      lo, hi;      /* clamped on set */
    bool        readonly;    /* derived; shown, not edited */
    bool        is_enum;
    const char *const *names;
    int         nnames;
    bool        heading;     /* a section rule, not a value */
} Field;

#define LS_INSPECT_MAX 32

/* Build the visible field list for the current selection, tier and unit system.
 * `sel_light` and `sel_prim` are indices, or -1. Returns the count. */
int ls_inspect_fields(const SceneDesc *d, int sel_light, int sel_prim,
                      LsUnitSystem units,
                      LsTier tier, Field *out, int max);

/* Apply an edit. Returns true if anything actually changed, so the caller knows
 * whether to re-solve. The caller is responsible for having taken an undo
 * snapshot first, and for calling ls_scene_update_light afterwards -- which this
 * does internally for light fields, since a light edit is meaningless without
 * its geometry following. */
bool ls_inspect_set(SceneDesc *d, int sel_light, int sel_prim,
                    FieldId id, double v);

/* Format a field's value for display. `buf` must hold at least 32 bytes. */
void ls_inspect_format(const Field *f, char *buf, size_t n);

#endif /* LIGHTSIM_VIEWER_INSPECT_H */
