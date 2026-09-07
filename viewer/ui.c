#include "ui.h"
#include <stddef.h>

static void add(Toolbar *t, UiAction a, const char *label, const char *hint,
                bool sep, int *y) {
    UiButton *b = &t->buttons[t->count++];
    b->action = a;
    b->label = label;
    b->hint = hint;
    b->separator_above = sep;
    b->enabled = true;
    b->active = false;
    if (sep) *y += UI_GROUP_GAP;
    b->rect.x = UI_MARGIN;
    b->rect.y = *y;
    b->rect.w = UI_BUTTON_W;
    b->rect.h = UI_BUTTON_H;
    *y += UI_BUTTON_H + UI_BUTTON_GAP;
}

void ui_init(Toolbar *t) {
    t->count = 0;
    t->hover = -1;
    t->pressed = -1;
    int y = UI_TOP_MARGIN;
    add(t, UI_MODE_GRID,   "FIELD MAP", "1", false, &y);
    add(t, UI_MODE_RENDER, "RENDER",    "2", false, &y);
    add(t, UI_UNITS,       "LUX",       "U", true,  &y);
    add(t, UI_TRANSPORT,   "FULL",      "T", false, &y);
    add(t, UI_QUALITY,     "DRAFT",     "Q", false, &y);
    add(t, UI_TIER,        "SIMPLE",    "3", false, &y);
    add(t, UI_ADD_LIGHT,   "ADD LIGHT", "A", true,  &y);
    add(t, UI_ADD_PART,    "ADD PART",  "P", false, &y);
    add(t, UI_DUPLICATE,   "DUPLICATE", "D", false, &y);
    add(t, UI_DELETE,      "DELETE",    "DEL", false, &y);
    add(t, UI_UNDO,        "UNDO",      "^Z", true,  &y);
    add(t, UI_REDO,        "REDO",      "^Y", false, &y);
    add(t, UI_SOLVE,       "SOLVE",     "R", true,  &y);
    add(t, UI_SAVE,        "SAVE PPM",  "S", false, &y);
    add(t, UI_SAVE_SCENE,  "SAVE SCENE","W", false, &y);
    add(t, UI_BLENDER,     "BLENDER",   "B", false, &y);
}

int ui_hit_test(const Toolbar *t, int x, int y) {
    for (int i = 0; i < t->count; ++i) {
        const UiRect *r = &t->buttons[i].rect;
        if (x >= r->x && x < r->x + r->w && y >= r->y && y < r->y + r->h) return i;
    }
    return -1;
}

bool ui_contains(const Toolbar *t, int x, int y) {
    (void)t;
    return x >= 0 && x < UI_TOOLBAR_W && y >= 0;
}

const char *ui_label(const Toolbar *t, int index, UiState s) {
    switch (t->buttons[index].action) {
        case UI_UNITS:     return s.photometric ? "LUX" : "WATT/M2";
        case UI_TRANSPORT: return s.direct_only ? "DIRECT" : "FULL";
        case UI_QUALITY:   return s.high_quality ? "FINE" : "DRAFT";
        case UI_SOLVE:     return s.solving ? "SOLVING" : "SOLVE";
        case UI_TIER:      return s.tier == 0 ? "SIMPLE"
                                : s.tier == 1 ? "ADVANCED" : "SCIENTIFIC";
        default:           return t->buttons[index].label;
    }
}

void ui_apply_state(Toolbar *t, UiState s) {
    for (int i = 0; i < t->count; ++i) {
        UiButton *b = &t->buttons[i];
        switch (b->action) {
            case UI_MODE_GRID:
                b->active = s.grid_mode;
                b->enabled = s.has_grid;
                break;
            case UI_MODE_RENDER:
                b->active = !s.grid_mode;
                b->enabled = s.has_camera;
                break;
            case UI_UNITS:
                b->active = s.photometric;
                b->enabled = true;
                break;
            case UI_TRANSPORT:
                b->active = !s.direct_only;
                b->enabled = true;
                break;
            case UI_QUALITY:
                b->active = s.high_quality;
                b->enabled = true;
                break;
            case UI_SOLVE:
                b->active = s.solving;
                /* A solve only means something for the field map; the render
                 * accumulates continuously on its own. */
                b->enabled = s.has_grid && !s.solving && s.grid_mode;
                break;
            case UI_TIER:
                b->enabled = true;
                b->active = (s.tier > 0);
                break;
            case UI_ADD_LIGHT:
                /* Placement needs somewhere to click, which means the 3D view. */
                b->enabled = s.has_camera && !s.grid_mode;
                b->active = (s.tool == UI_TOOL_LIGHT);
                break;
            case UI_ADD_PART:
                b->enabled = s.has_camera && !s.grid_mode;
                b->active = (s.tool == UI_TOOL_PART);
                break;
            case UI_DUPLICATE:
            case UI_DELETE:
                b->enabled = s.has_selection;
                b->active = false;
                break;
            case UI_UNDO:
                b->enabled = s.can_undo;
                b->active = false;
                break;
            case UI_REDO:
                b->enabled = s.can_redo;
                b->active = false;
                break;
            case UI_SAVE:
            case UI_SAVE_SCENE:
                b->enabled = true;
                b->active = false;
                break;
            case UI_BLENDER:
                b->enabled = s.has_grid;
                b->active = false;
                break;
            default:
                break;
        }
    }
}
