#include "ui.h"
#include <stddef.h>

/* The toolbar's contents, in order. A `true` group flag starts a new group,
 * drawn with a separator rule above it and a little extra breathing room. */
typedef struct {
    UiAction    action;
    const char *label;
    const char *hint;
    const char *tip;        /* hover text: what it does, and what it wants */
    bool        group_start;
} ButtonSpec;

static const ButtonSpec BUTTON_SPECS[] = {
    { UI_VIEW_3D,    "3D VIEW",   "1",
      "Look at the scene in perspective. Drag to orbit, wheel to zoom.", false },
    { UI_VIEW_TOP,   "TOP VIEW",  "2",
      "Look straight down, orthographic and to scale, framed on the measurement "
      "grid so geometry and result line up.", false },
    { UI_VIEW_HEAT,  "HEAT MAP",  "3",
      "Colour every surface by the illuminance arriving at it instead of the "
      "light leaving it. Not a third camera: it works through whichever view "
      "is showing.", false },
    { UI_UNITS,      "LUX",       "U",
      "Switch between photometric units, weighted by the eye's response, and "
      "raw radiometric watts per square metre.", true  },
    { UI_TRANSPORT,  "FULL",      "T",
      "Count only light straight from the sources, or every bounce as well. "
      "The difference is the contribution of the room itself.", false },
    { UI_QUALITY,    "DRAFT",     "Q",
      "Trade noise for speed. DRAFT refines while you work; FINE costs more "
      "per pass and is what a reported number should come from.", false },
    { UI_TIER,       "SIMPLE",    "V",
      "How much of the physics the inspector exposes. SIMPLE shows what you "
      "place and aim, SCIENTIFIC every spectral and sampling parameter.", false },
    { UI_ADD_LIGHT,  "ADD LIGHT", "A",
      "Arm the light tool, then click a surface to drop a luminaire on it. "
      "Escape disarms.", true  },
    { UI_ADD_PART,   "ADD PART",  "P",
      "Arm the part tool, then click a surface to stand a block on it. "
      "Escape disarms.", false },
    { UI_IMPORT,     "IMPORT",    "I",
      "Bring in geometry from an OBJ or STL file. Dropping the file on the "
      "window does the same thing.", false },
    { UI_DUPLICATE,  "DUPLICATE", "D",
      "Copy whatever is selected, offset a little so you can see both. "
      "Select a light or a part first.", false },
    { UI_DELETE,     "DELETE",    "DEL",
      "Remove whatever is selected. Select a light or a part first.", false },
    { UI_UNDO,       "UNDO",      "^Z",
      "Step back through your edits.", true  },
    { UI_REDO,       "REDO",      "^Y",
      "Step forward again through undone edits.", false },
    { UI_SOLVE,      "SOLVE",     "R",
      "Re-measure the field on the grid at the current quality. Needs a scene "
      "with a measurement grid; the render refines on its own.", true  },
    { UI_SAVE,       "SAVE PPM",  "S",
      "Write the picture and the field out: the render as a PPM, the grid as "
      "CSV and JSON, alongside a false-colour map.", false },
    { UI_SAVE_SCENE, "SAVE SCENE","W",
      "Save the scene itself -- lights, parts, camera and grid -- so you can "
      "come back to it. Unlike SAVE PPM, which writes results.", false },
    { UI_BLENDER,    "BLENDER",   "B",
      "Write a Blender script that rebuilds this scene with the measured field "
      "draped over the grid. Needs a measurement grid.", false },
    { UI_HELP,       "HELP",      "H",
      "List every key and what it does. Everything here has a button as well, "
      "and every button shows its key on the right.", true  },
};
#define BUTTON_SPEC_COUNT ((int)(sizeof BUTTON_SPECS / sizeof BUTTON_SPECS[0]))

void ui_init(Toolbar *t, int strip_height) {
    t->count = 0;
    t->hover = -1;
    t->pressed = -1;

    int n = BUTTON_SPEC_COUNT < UI_ACTION_COUNT ? BUTTON_SPEC_COUNT : UI_ACTION_COUNT;
    int groups = 0;
    for (int i = 1; i < n; ++i) if (BUTTON_SPECS[i].group_start) groups++;

    /* Squeeze to fit rather than running off the bottom of the window, where a
     * button can neither be read nor pressed. Spacing goes first, since it
     * costs the least, and only then the buttons themselves -- down to a floor
     * where the label would stop being legible. */
    int button_h = UI_BUTTON_H, gap = UI_BUTTON_GAP, group_gap = UI_GROUP_GAP;
    int needed = UI_TOP_MARGIN + n * (button_h + gap) + groups * group_gap;
    if (strip_height > 0 && needed > strip_height) {
        gap = 2;
        group_gap = 4;
        needed = UI_TOP_MARGIN + n * (button_h + gap) + groups * group_gap;
        if (needed > strip_height) {
            int room = strip_height - UI_TOP_MARGIN - groups * group_gap - n * gap;
            button_h = room / n;
            if (button_h < UI_BUTTON_MIN_H) button_h = UI_BUTTON_MIN_H;
        }
    }

    int y = UI_TOP_MARGIN;
    for (int i = 0; i < n; ++i) {
        const ButtonSpec *spec = &BUTTON_SPECS[i];
        if (spec->group_start && i > 0) y += group_gap;

        UiButton *b = &t->buttons[t->count++];
        b->action = spec->action;
        b->label = spec->label;
        b->hint = spec->hint;
        b->tip = spec->tip;
        b->rect = (UiRect){ UI_MARGIN, y, UI_BUTTON_W, button_h };
        b->enabled = true;
        b->active = false;
        b->separator_above = spec->group_start && i > 0;

        y += button_h + gap;
    }
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

static const UiButton *button_for_action(const Toolbar *t, UiAction action) {
    for (int i = 0; i < t->count; ++i)
        if (t->buttons[i].action == action) return &t->buttons[i];
    return NULL;
}

bool ui_action_enabled(const Toolbar *t, UiAction action) {
    const UiButton *b = button_for_action(t, action);
    return b && b->enabled;
}

const char *ui_action_tip(const Toolbar *t, UiAction action) {
    const UiButton *b = button_for_action(t, action);
    return b ? b->tip : "";
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
            case UI_VIEW_3D:
                b->active = (s.view == 0);
                b->enabled = true;
                break;
            case UI_VIEW_TOP:
                b->active = (s.view == 1);
                b->enabled = true;
                break;
            case UI_VIEW_HEAT:
                /* Not a third camera: a different quantity to show through
                 * whichever camera is active, on every surface in the scene. */
                b->active = s.shade_heat;
                b->enabled = true;
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
                b->enabled = s.has_grid && !s.solving;
                break;
            case UI_TIER:
                b->enabled = true;
                b->active = (s.tier > 0);
                break;
            case UI_ADD_LIGHT:
                /* A placement tool, not an operation on a selection: always
                 * offered, and lit while it is armed. */
                b->enabled = true;
                b->active = (s.tool == UI_TOOL_LIGHT);
                break;
            case UI_ADD_PART:
                b->enabled = true;
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
            case UI_IMPORT:
                /* Always available: an import needs no selection and no grid,
                 * and dropping a file on the window does the same thing. */
                b->enabled = true;
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
            case UI_HELP:
                b->enabled = true;
                b->active = s.help_open;
                break;
            case UI_NONE:
            case UI_ACTION_COUNT:
                break;
        }
    }
}
