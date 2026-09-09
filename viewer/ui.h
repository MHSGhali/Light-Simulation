/* ui.h — toolbar model: layout, hit-testing, and per-button state.
 * Deliberately free of SDL so the rules can be exercised headlessly; drawing
 * lives in draw.c. */
#ifndef LIGHTSIM_VIEWER_UI_H
#define LIGHTSIM_VIEWER_UI_H

#include <stdbool.h>

/* Wide enough that the longest label ("SCIENTIFIC", "SAVE SCENE") still clears
 * the hotkey hint drawn against the button's right edge. test_ui holds this to
 * a real gap, not merely to non-overlap. */
#define UI_TOOLBAR_W    168
#define UI_BUTTON_W     148
#define UI_BUTTON_H      30
#define UI_BUTTON_MIN_H  18     /* how far buttons squeeze to fit a short window */
#define UI_BUTTON_GAP     4
#define UI_GROUP_GAP     10
#define UI_MARGIN        10
#define UI_TOP_MARGIN    12

typedef enum {
    UI_NONE = 0,
    UI_VIEW_3D,
    UI_VIEW_TOP,
    UI_VIEW_HEAT,
    UI_UNITS,
    UI_TRANSPORT,
    UI_QUALITY,
    UI_TIER,
    UI_ADD_LIGHT,
    UI_ADD_PART,
    UI_DUPLICATE,
    UI_DELETE,
    UI_UNDO,
    UI_REDO,
    UI_SOLVE,
    UI_SAVE,
    UI_SAVE_SCENE,
    UI_BLENDER,
    UI_IMPORT,
    UI_HELP,
    UI_ACTION_COUNT
} UiAction;

/* Which placement tool is armed, if any. Arming is a toggle and the armed
 * button lights in the same colour as a selection, so "this tool is armed" and
 * "this object is selected" read as one visual language. */
typedef enum { UI_TOOL_NONE = 0, UI_TOOL_LIGHT, UI_TOOL_PART } UiTool;

typedef struct { int x, y, w, h; } UiRect;

typedef struct {
    UiAction    action;
    const char *label;
    const char *hint;            /* hotkey reminder ('^' draws a caret) */
    const char *tip;             /* what it does and what it needs, on hover */
    UiRect      rect;
    bool        enabled;
    bool        active;          /* toggled on / currently selected */
    bool        separator_above;
} UiButton;

typedef struct {
    UiButton buttons[UI_ACTION_COUNT];
    int count;
    int hover;                   /* index under the cursor, or -1 */
    int pressed;                 /* index the mouse went down on, or -1 */
} Toolbar;

/* Everything ui_apply_state needs, so ui.c never reaches into the app. */
typedef struct {
    int    view;                 /* 0 = 3D perspective, 1 = top plan */
    bool   shade_heat;           /* false-colour illuminance on every surface */
    bool   photometric;
    bool   direct_only;
    bool   high_quality;
    bool   solving;              /* a solve is in flight */
    bool   has_grid;             /* the scene declares a measurement grid */
    bool   has_camera;
    int    tier;                 /* 0 simple, 1 advanced, 2 scientific */
    UiTool tool;                 /* armed placement tool */
    bool   has_selection;
    bool   can_undo, can_redo;
    bool   help_open;            /* the key list is showing */
} UiState;

/* Builds the button layout to fit a strip `strip_height` tall (the window's
 * height). Call at startup and again on resize: with a fixed layout the
 * buttons below the fold are simply unreachable. A height of 0 means "do not
 * squeeze", for headless callers that only want the rules. */
void ui_init(Toolbar *t, int strip_height);

int  ui_hit_test(const Toolbar *t, int x, int y);
bool ui_contains(const Toolbar *t, int x, int y);
void ui_apply_state(Toolbar *t, UiState s);

/* Label a button should currently show (toggles change wording). */
const char *ui_label(const Toolbar *t, int index, UiState s);

/* Whether `action`'s button is enabled right now, and its hover text. The
 * hotkeys go through these so a key can never disagree with the button beside
 * it -- and so a refused key can explain itself in the same words the
 * button's tooltip uses. Unknown actions read as disabled. */
bool        ui_action_enabled(const Toolbar *t, UiAction action);
const char *ui_action_tip(const Toolbar *t, UiAction action);

#endif /* LIGHTSIM_VIEWER_UI_H */
