/* ui.h — toolbar model: layout, hit-testing, and per-button state.
 * Deliberately free of SDL so the rules can be exercised headlessly; drawing
 * lives in draw.c. */
#ifndef LIGHTSIM_VIEWER_UI_H
#define LIGHTSIM_VIEWER_UI_H

#include <stdbool.h>

#define UI_TOOLBAR_W  152
#define UI_BUTTON_W   132
#define UI_BUTTON_H    28
#define UI_BUTTON_GAP   4
#define UI_GROUP_GAP   12
#define UI_MARGIN      10
#define UI_TOP_MARGIN  12

typedef enum {
    UI_NONE = 0,
    UI_MODE_GRID,
    UI_MODE_RENDER,
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
    const char *hint;            /* hotkey reminder */
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
    bool   grid_mode;            /* showing the field map rather than the render */
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
} UiState;

void ui_init(Toolbar *t);
int  ui_hit_test(const Toolbar *t, int x, int y);
bool ui_contains(const Toolbar *t, int x, int y);
void ui_apply_state(Toolbar *t, UiState s);

/* Label a button should currently show (toggles change wording). */
const char *ui_label(const Toolbar *t, int index, UiState s);

#endif /* LIGHTSIM_VIEWER_UI_H */
