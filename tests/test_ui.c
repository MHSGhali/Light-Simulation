/* The toolbar model is kept free of SDL precisely so its rules can be checked
 * headlessly. These run in the normal suite, with no window and no SDL2. */
#include "test.h"
#include "tests.h"
#include "../viewer/ui.h"
#include <string.h>
#include "../viewer/font.h"

void test_ui(void) {
    SECTION("toolbar layout and rules");
    {
        Toolbar t;
        ui_init(&t);
        CHECK(t.count == UI_ACTION_COUNT - 1);   /* every action but UI_NONE */

        /* Buttons stack without overlapping, and stay inside the strip. */
        for (int i = 0; i < t.count; ++i) {
            CHECK(t.buttons[i].rect.x >= 0);
            CHECK(t.buttons[i].rect.x + t.buttons[i].rect.w <= UI_TOOLBAR_W);
            if (i > 0) {
                const UiRect *a = &t.buttons[i-1].rect, *b = &t.buttons[i].rect;
                CHECK(b->y >= a->y + a->h);
            }
        }

        /* Hit-testing finds each button at its own centre, and nothing in the
         * gap above the first one. */
        for (int i = 0; i < t.count; ++i) {
            const UiRect *r = &t.buttons[i].rect;
            CHECK(ui_hit_test(&t, r->x + r->w / 2, r->y + r->h / 2) == i);
        }
        CHECK(ui_hit_test(&t, UI_MARGIN, 0) == -1);
        CHECK(ui_hit_test(&t, UI_TOOLBAR_W + 50, 40) == -1);
        CHECK(ui_contains(&t, 10, 400));
        CHECK(!ui_contains(&t, UI_TOOLBAR_W + 1, 400));
    }

    SECTION("toolbar state");
    {
        Toolbar t;
        ui_init(&t);
        UiState s;
        memset(&s, 0, sizeof s);
        s.has_grid = true; s.has_camera = true; s.grid_mode = true;
        s.photometric = true;
        ui_apply_state(&t, s);

        for (int i = 0; i < t.count; ++i) {
            const UiButton *b = &t.buttons[i];
            if (b->action == UI_MODE_GRID)   CHECK(b->active);
            if (b->action == UI_MODE_RENDER) CHECK(!b->active);
            /* Solving the field only means something while it is on screen. */
            if (b->action == UI_SOLVE)       CHECK(b->enabled);
        }

        /* In render mode SOLVE has nothing to do -- the render accumulates on
         * its own -- so it must go disabled rather than silently no-op. */
        s.grid_mode = false;
        ui_apply_state(&t, s);
        for (int i = 0; i < t.count; ++i)
            if (t.buttons[i].action == UI_SOLVE) CHECK(!t.buttons[i].enabled);

        /* A scene with no grid cannot show the field map or export one. */
        s.has_grid = false; s.grid_mode = false;
        ui_apply_state(&t, s);
        for (int i = 0; i < t.count; ++i) {
            const UiButton *b = &t.buttons[i];
            if (b->action == UI_MODE_GRID) CHECK(!b->enabled);
            if (b->action == UI_BLENDER)   CHECK(!b->enabled);
        }

        /* Toggle labels track the state they report. */
        s.has_grid = true; s.grid_mode = true;
        ui_apply_state(&t, s);
        for (int i = 0; i < t.count; ++i) {
            if (t.buttons[i].action == UI_UNITS) {
                s.photometric = true;
                CHECK(strcmp(ui_label(&t, i, s), "LUX") == 0);
                s.photometric = false;
                CHECK(strcmp(ui_label(&t, i, s), "WATT/M2") == 0);
            }
            if (t.buttons[i].action == UI_TRANSPORT) {
                s.direct_only = false;
                CHECK(strcmp(ui_label(&t, i, s), "FULL") == 0);
                s.direct_only = true;
                CHECK(strcmp(ui_label(&t, i, s), "DIRECT") == 0);
            }
        }
    }

    SECTION("font coverage");
    {
        /* Every character the viewer actually draws must have a glyph, or a
         * label silently renders as blanks. */
        const char *used = "FIELD MAP RENDER LUX WATT/M2 FULL DIRECT DRAFT FINE "
                           "SOLVE SOLVING SAVE PPM BLENDER PROBE ILLUMINANCE "
                           "IRRADIANCE MIN MAX MEAN STDDEV U0 UD CONTRAST "
                           "CROSS SECTION AT Y M PASSES SPP DRAG TO ORBIT "
                           "WHEEL ZOOM 0123456789.-";
        for (const char *p = used; *p; ++p) {
            const unsigned char *g = font_glyph(*p);
            bool blank = (g[0] | g[1] | g[2] | g[3] | g[4]) == 0;
            if (*p != ' ') CHECK(!blank);
        }
        CHECK(font_text_width("ABC", 2) == (3 * 6 - 1) * 2);
    }
}
