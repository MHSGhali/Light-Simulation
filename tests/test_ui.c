/* The toolbar model is kept free of SDL precisely so its rules can be checked
 * headlessly. These run in the normal suite, with no window and no SDL2. */
#include "test.h"
#include "tests.h"
#include "../viewer/ui.h"
#include "../viewer/status.h"
#include <string.h>
#include "../viewer/font.h"

void test_ui(void) {
    SECTION("toolbar layout and rules");
    {
        Toolbar t;
        ui_init(&t, 0);
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

        /* Every button explains itself on hover, and every one of those words
         * is drawable -- a tip with no glyph renders as a blank panel. */
        for (int i = 0; i < t.count; ++i) {
            CHECK(t.buttons[i].tip != NULL && t.buttons[i].tip[0] != '\0');
            CHECK(t.buttons[i].hint != NULL && t.buttons[i].hint[0] != '\0');
            for (const char *c = t.buttons[i].tip; *c; ++c) {
                const unsigned char *g = font_glyph(*c);
                if (*c != ' ') CHECK((g[0]|g[1]|g[2]|g[3]|g[4]) != 0);
            }
        }

        /* The widest label still clears the hotkey hint beside it, with space
         * to spare. Merely not overlapping is not enough: two runs of text that
         * touch read as one unpronounceable word. */
        UiState widest;
        memset(&widest, 0, sizeof widest);
        widest.tier = 2;                          /* "SCIENTIFIC", the longest */
        for (int i = 0; i < t.count; ++i) {
            int label_end = 9 + font_text_width(ui_label(&t, i, widest), 2);
            int hint_start = UI_BUTTON_W - 9 - font_text_width(t.buttons[i].hint, 1);
            CHECK(hint_start - label_end >= 6);
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

    SECTION("toolbar fits the window it is drawn in");
    {
        /* A layout taller than the strip puts buttons below the fold, where
         * they can be neither read nor pressed. */
        Toolbar full, squeezed;
        ui_init(&full, 0);
        int needed = full.buttons[full.count-1].rect.y +
                     full.buttons[full.count-1].rect.h;

        ui_init(&squeezed, needed / 2);
        CHECK(squeezed.count == full.count);
        int bottom = squeezed.buttons[squeezed.count-1].rect.y +
                     squeezed.buttons[squeezed.count-1].rect.h;
        CHECK(bottom < needed);
        for (int i = 0; i < squeezed.count; ++i)
            CHECK(squeezed.buttons[i].rect.h >= UI_BUTTON_MIN_H);

        /* Never squeezed past legibility, even for a window far too short --
         * better to run off the bottom than to draw a row of slivers. */
        Toolbar tiny;
        ui_init(&tiny, 40);
        for (int i = 0; i < tiny.count; ++i)
            CHECK(tiny.buttons[i].rect.h == UI_BUTTON_MIN_H);

        /* A strip with room to spare is left alone. */
        Toolbar roomy;
        ui_init(&roomy, needed * 2);
        for (int i = 0; i < roomy.count; ++i)
            CHECK(roomy.buttons[i].rect.h == UI_BUTTON_H);
    }

    SECTION("toolbar state");
    {
        Toolbar t;
        ui_init(&t, 0);
        UiState s;
        memset(&s, 0, sizeof s);
        s.has_grid = true; s.has_camera = true; s.view = 0;
        ui_apply_state(&t, s);

        /* Exactly one CAMERA button is active at a time. Heat is not a third
         * camera -- it is a different quantity shown through whichever camera
         * is active, on every surface -- so it toggles independently. */
        for (int v = 0; v < 2; ++v) {
            s.view = v;
            ui_apply_state(&t, s);
            int active = 0;
            for (int i = 0; i < t.count; ++i) {
                const UiButton *b = &t.buttons[i];
                if (b->action == UI_VIEW_3D || b->action == UI_VIEW_TOP) {
                    if (b->active) active++;
                    int want = (b->action == UI_VIEW_3D) ? 0 : 1;
                    CHECK(b->active == (want == v));
                }
            }
            CHECK(active == 1);
        }

        /* The heat toggle tracks its own flag and never disables the cameras. */
        for (int hflag = 0; hflag < 2; ++hflag) {
            s.shade_heat = (hflag != 0);
            ui_apply_state(&t, s);
            for (int i = 0; i < t.count; ++i) {
                const UiButton *b = &t.buttons[i];
                if (b->action == UI_VIEW_HEAT) {
                    CHECK(b->active == s.shade_heat);
                    CHECK(b->enabled);
                }
                /* Placement stays available whichever quantity is on screen,
                 * because both views still have geometry to click on. */
                if (b->action == UI_ADD_LIGHT || b->action == UI_ADD_PART)
                    CHECK(b->enabled);
            }
        }
        s.shade_heat = false;

        /* A scene with no measurement grid cannot export one. Heat shading is
         * unaffected: it measures the surfaces themselves, not the grid. */
        s.has_grid = false;
        ui_apply_state(&t, s);
        for (int i = 0; i < t.count; ++i) {
            const UiButton *b = &t.buttons[i];
            if (b->action == UI_BLENDER)   CHECK(!b->enabled);
            if (b->action == UI_VIEW_HEAT) CHECK(b->enabled);
        }

        /* Edit actions need something selected; history needs history. */
        s.has_grid = true; s.has_selection = false;
        s.can_undo = false; s.can_redo = false;
        ui_apply_state(&t, s);
        for (int i = 0; i < t.count; ++i) {
            const UiButton *b = &t.buttons[i];
            if (b->action == UI_DELETE || b->action == UI_DUPLICATE) CHECK(!b->enabled);
            if (b->action == UI_UNDO || b->action == UI_REDO)        CHECK(!b->enabled);
        }
        s.has_selection = true; s.can_undo = true;
        ui_apply_state(&t, s);
        for (int i = 0; i < t.count; ++i) {
            const UiButton *b = &t.buttons[i];
            if (b->action == UI_DELETE || b->action == UI_DUPLICATE) CHECK(b->enabled);
            if (b->action == UI_UNDO) CHECK(b->enabled);
            if (b->action == UI_REDO) CHECK(!b->enabled);
        }

        /* Toggle labels track the state they report. */
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
            if (t.buttons[i].action == UI_TIER) {
                s.tier = 0; CHECK(strcmp(ui_label(&t, i, s), "SIMPLE") == 0);
                s.tier = 2; CHECK(strcmp(ui_label(&t, i, s), "SCIENTIFIC") == 0);
            }
        }

        /* HELP lights while the key list is showing, and is never refused --
         * the one button that has to work when nothing else does. */
        for (int hflag = 0; hflag < 2; ++hflag) {
            s.help_open = (hflag != 0);
            ui_apply_state(&t, s);
            for (int i = 0; i < t.count; ++i)
                if (t.buttons[i].action == UI_HELP) {
                    CHECK(t.buttons[i].active == s.help_open);
                    CHECK(t.buttons[i].enabled);
                }
        }
        s.help_open = false;

        /* Looking a rule up by action agrees with the button carrying it, so a
         * hotkey and the button beside it can never disagree. */
        s.has_grid = false;
        ui_apply_state(&t, s);
        CHECK(!ui_action_enabled(&t, UI_BLENDER));
        CHECK(ui_action_enabled(&t, UI_HELP));
        CHECK(ui_action_tip(&t, UI_BLENDER)[0] != '\0');
        /* An action with no button reads as refused rather than crashing. */
        CHECK(!ui_action_enabled(&t, UI_NONE));
        CHECK(ui_action_tip(&t, UI_NONE)[0] == '\0');
    }

    SECTION("message log");
    {
        StatusLog log;
        status_init(&log);
        const StatusMessage *seen[STATUS_HISTORY];
        CHECK(status_visible(&log, 0, seen, STATUS_HISTORY) == 0);
        CHECK(!log.sticky_set);

        /* Newest first, and only as many as are kept. */
        status_push(&log, STATUS_INFO, 1000, "FIRST");
        status_push(&log, STATUS_WARN, 1000, "SECOND");
        CHECK(status_visible(&log, 1000, seen, STATUS_HISTORY) == 2);
        CHECK(strcmp(seen[0]->text, "SECOND") == 0);
        CHECK(seen[0]->level == STATUS_WARN);
        CHECK(strcmp(seen[1]->text, "FIRST") == 0);

        /* The same thing said twice running is one thing that is still true:
         * it restamps rather than scrolling the line away. */
        status_push(&log, STATUS_WARN, 4000, "SECOND");
        CHECK(status_visible(&log, 4000, seen, STATUS_HISTORY) == 2);
        CHECK(seen[0]->stamp_ms == 4000);

        /* Oldest falls off the end rather than growing without bound. */
        status_push(&log, STATUS_INFO, 4000, "THIRD");
        status_push(&log, STATUS_INFO, 4000, "FOURTH");
        int n = status_visible(&log, 4000, seen, STATUS_HISTORY);
        CHECK(n == STATUS_HISTORY);
        CHECK(strcmp(seen[0]->text, "FOURTH") == 0);
        for (int i = 0; i < n; ++i) CHECK(strcmp(seen[i]->text, "FIRST") != 0);

        /* Solid, then fading, then gone -- and gone means invisible. */
        status_init(&log);
        status_push(&log, STATUS_INFO, 1000, "HELLO");
        CHECK(status_alpha(&log.recent[0], 1000) == 255);
        CHECK(status_alpha(&log.recent[0], 1000 + STATUS_FADE_MS
                                               - STATUS_FADE_TAIL_MS) == 255);
        int mid = status_alpha(&log.recent[0],
                              1000 + STATUS_FADE_MS - STATUS_FADE_TAIL_MS / 2);
        CHECK(mid > 0 && mid < 255);
        CHECK(status_alpha(&log.recent[0], 1000 + STATUS_FADE_MS) == 0);
        CHECK(status_visible(&log, 1000 + STATUS_FADE_MS, seen, STATUS_HISTORY) == 0);

        /* Unsigned arithmetic, so a message does not become immortal when
         * SDL_GetTicks wraps past 2^32. */
        status_init(&log);
        status_push(&log, STATUS_INFO, 0xFFFFFF00u, "NEAR THE WRAP");
        CHECK(status_alpha(&log.recent[0], 0xFFFFFF00u + 100u) == 255);
        CHECK(status_alpha(&log.recent[0],
                           0xFFFFFF00u + STATUS_FADE_MS + 1u) == 0);

        /* The banner is a condition, so it does not fade and it can be
         * withdrawn when the condition stops holding. */
        status_set_sticky(&log, STATUS_ERROR, "JAMMED");
        CHECK(log.sticky_set);
        CHECK(log.sticky_level == STATUS_ERROR);
        CHECK(strcmp(log.sticky, "JAMMED") == 0);
        status_set_sticky(&log, STATUS_INFO, "");   /* empty clears it */
        CHECK(!log.sticky_set);
        status_set_sticky(&log, STATUS_WARN, "HOLDING");
        status_clear_sticky(&log);
        CHECK(!log.sticky_set);

        /* Nothing to say is not a message. */
        status_init(&log);
        status_push(&log, STATUS_INFO, 10, "");
        status_push(&log, STATUS_INFO, 10, NULL);
        CHECK(log.count == 0);
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
