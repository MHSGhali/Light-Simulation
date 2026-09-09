#include "draw.h"
#include "font.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

/* Grounds and rules. */
const Col COL_BG        = { 0x14, 0x14, 0x14 };
const Col COL_PANEL     = { 0x1E, 0x1E, 0x22 };
const Col COL_RULE      = { 0x37, 0x39, 0x40 };
/* Button faces: resting, under the cursor, held down, and refused. A pressed
 * button is lighter than a hovered one, so a click that lands is felt. */
const Col COL_FACE      = { 0x2D, 0x2F, 0x34 };
const Col COL_FACE_HOT  = { 0x3E, 0x41, 0x48 };
const Col COL_FACE_DOWN = { 0x4E, 0x52, 0x5C };
const Col COL_FACE_OFF  = { 0x22, 0x23, 0x27 };
const Col COL_EDGE      = { 0x46, 0x49, 0x52 };
const Col COL_EDGE_OFF  = { 0x30, 0x31, 0x36 };
/* Text. A disabled button dims its ink rather than its alpha: a half-
 * transparent label over a dark panel reads as blurred, not as unavailable. */
const Col COL_INK       = { 0xCD, 0xD2, 0xDC };
const Col COL_HEAD      = { 0x96, 0x9C, 0xAC };
const Col COL_MUTED     = { 0x82, 0x88, 0x94 };
const Col COL_DIM       = { 0x5C, 0x5F, 0x66 };
const Col COL_ACCENT    = { 0xFF, 0xE1, 0x46 };
const Col COL_ACCENT_DIM= { 0xC8, 0xB4, 0x46 };
const Col COL_WARN      = { 0xF5, 0xC3, 0x5A };
const Col COL_ERR       = { 0xF0, 0x6E, 0x5F };
const Col COL_S1        = { 0x50, 0xC8, 0xDC };
const Col COL_S2        = { 0xF0, 0xAA, 0x50 };
const Col COL_AXIS_X    = { 0xE3, 0x49, 0x48 };
const Col COL_AXIS_Y    = { 0x1B, 0xAF, 0x7A };
const Col COL_AXIS_Z    = { 0x39, 0x87, 0xE5 };

Col draw_level_col(StatusLevel level) {
    switch (level) {
        case STATUS_WARN:  return COL_WARN;
        case STATUS_ERROR: return COL_ERR;
        case STATUS_INFO:
        default:           return COL_INK;
    }
}

void draw_rect_fill(SDL_Renderer *ren, int x, int y, int w, int h, Col c, Uint8 a) {
    SDL_SetRenderDrawColor(ren, c.r, c.g, c.b, a);
    SDL_Rect r = { x, y, w, h };
    SDL_RenderFillRect(ren, &r);
}

void draw_rect_line(SDL_Renderer *ren, int x, int y, int w, int h, Col c, Uint8 a) {
    SDL_SetRenderDrawColor(ren, c.r, c.g, c.b, a);
    SDL_Rect r = { x, y, w, h };
    SDL_RenderDrawRect(ren, &r);
}

void draw_line(SDL_Renderer *ren, int x0, int y0, int x1, int y1, Col c, Uint8 a) {
    SDL_SetRenderDrawColor(ren, c.r, c.g, c.b, a);
    SDL_RenderDrawLine(ren, x0, y0, x1, y1);
}

void draw_marker(SDL_Renderer *ren, int x, int y, int r, Col c, Uint8 a) {
    SDL_SetRenderDrawColor(ren, c.r, c.g, c.b, a);
    SDL_RenderDrawLine(ren, x - r, y, x + r, y);
    SDL_RenderDrawLine(ren, x, y - r, x, y + r);
    SDL_Rect box = { x - r, y - r, 2 * r, 2 * r };
    SDL_RenderDrawRect(ren, &box);
}

void draw_text(SDL_Renderer *ren, int x, int y, int scale, const char *s, Col c, Uint8 a) {
    SDL_SetRenderDrawColor(ren, c.r, c.g, c.b, a);
    int pen = x;
    for (const char *p = s; *p; ++p) {
        for (int col = 0; col < FONT_W; ++col)
            for (int row = 0; row < FONT_H; ++row)
                if (font_pixel(*p, col, row)) {
                    SDL_Rect px = { pen + col * scale, y + row * scale, scale, scale };
                    SDL_RenderFillRect(ren, &px);
                }
        pen += (FONT_W + 1) * scale;
    }
}

void draw_text_right(SDL_Renderer *ren, int right_x, int y, int scale, const char *s, Col c, Uint8 a) {
    draw_text(ren, right_x - font_text_width(s, scale), y, scale, s, c, a);
}

void draw_text_mid(SDL_Renderer *ren, int mid_x, int y, int scale, const char *s, Col c, Uint8 a) {
    draw_text(ren, mid_x - font_text_width(s, scale) / 2, y, scale, s, c, a);
}

void draw_viridis(double t, Uint8 *r, Uint8 *g, Uint8 *b) {
    static const double V[11][3] = {
        {0.267004,0.004874,0.329415},{0.282623,0.140926,0.457517},
        {0.253935,0.265254,0.529983},{0.206756,0.371758,0.553117},
        {0.163625,0.471133,0.558148},{0.127568,0.566949,0.550556},
        {0.134692,0.658636,0.517649},{0.266941,0.748751,0.440573},
        {0.477504,0.821444,0.318195},{0.741388,0.873449,0.149561},
        {0.993248,0.906157,0.143936}
    };
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    double s = t * 10.0;
    int i = (int)s;
    if (i > 9) i = 9;
    double u = s - i;
    *r = (Uint8)(255.0 * (V[i][0] + u * (V[i+1][0] - V[i][0])) + 0.5);
    *g = (Uint8)(255.0 * (V[i][1] + u * (V[i+1][1] - V[i][1])) + 0.5);
    *b = (Uint8)(255.0 * (V[i][2] + u * (V[i+1][2] - V[i][2])) + 0.5);
}

/* -------------------------------------------------------------- toolbar --- */
#define LABEL_SCALE 2
#define HINT_SCALE  1
#define TEXT_PAD    9

void draw_toolbar(SDL_Renderer *ren, const Toolbar *t, UiState st, int strip_h) {
    draw_rect_fill(ren, 0, 0, UI_TOOLBAR_W, strip_h, COL_PANEL, 255);
    draw_line(ren, UI_TOOLBAR_W - 1, 0, UI_TOOLBAR_W - 1, strip_h, COL_RULE, 255);

    for (int i = 0; i < t->count; ++i) {
        const UiButton *b = &t->buttons[i];
        const UiRect *r = &b->rect;

        if (b->separator_above) {
            int sy = r->y - UI_GROUP_GAP / 2 - UI_BUTTON_GAP / 2;
            draw_line(ren, r->x + 8, sy, r->x + r->w - 8, sy, COL_RULE, 255);
        }

        Col face = !b->enabled      ? COL_FACE_OFF
                 : (t->pressed == i) ? COL_FACE_DOWN
                 : (t->hover == i)   ? COL_FACE_HOT
                                     : COL_FACE;
        draw_rect_fill(ren, r->x, r->y, r->w, r->h, face, 255);

        /* A toggled-on button borrows the canvas's selection yellow, as an
         * outline rather than a fill: the label stays legible and the button
         * still reads as a button. */
        Col edge = b->active && b->enabled ? COL_ACCENT
                 : b->enabled              ? COL_EDGE
                                           : COL_EDGE_OFF;
        draw_rect_line(ren, r->x, r->y, r->w, r->h, edge, 255);

        Col ink  = !b->enabled ? COL_DIM : b->active ? COL_ACCENT     : COL_INK;
        Col hint = !b->enabled ? COL_EDGE : b->active ? COL_ACCENT_DIM : COL_MUTED;

        int lh = FONT_H * LABEL_SCALE, hh = FONT_H * HINT_SCALE;
        draw_text(ren, r->x + TEXT_PAD, r->y + (r->h - lh) / 2, LABEL_SCALE,
                  ui_label(t, i, st), ink, 255);
        if (b->hint && b->hint[0])
            draw_text_right(ren, r->x + r->w - TEXT_PAD, r->y + (r->h - hh) / 2,
                            HINT_SCALE, b->hint, hint, 255);
    }
}

/* ---------------------------------------------------- wrapped text boxes --- */

int draw_wrap(int scale, int max_w, const char *text,
              char lines[][DRAW_WRAP_LINE_CHARS], int max_lines, int *widest_out) {
    int count = 0, widest = 0;
    if (!text || !*text || max_lines <= 0) {
        if (widest_out) *widest_out = 0;
        return 0;
    }
    const char *p = text;
    while (*p && count < max_lines) {
        while (*p == ' ') ++p;
        if (!*p) break;
        int len = 0, last_space = -1;
        char buf[DRAW_WRAP_LINE_CHARS];
        while (p[len] && len < (int)sizeof buf - 1) {
            buf[len] = p[len];
            if (p[len] == ' ') last_space = len;
            buf[len + 1] = '\0';
            if (font_text_width(buf, scale) > max_w) {
                if (last_space > 0) { len = last_space; buf[len] = '\0'; }
                break;
            }
            ++len;
        }
        buf[len] = '\0';
        if (len == 0) break;
        snprintf(lines[count], DRAW_WRAP_LINE_CHARS, "%s", buf);
        int w = font_text_width(lines[count], scale);
        if (w > widest) widest = w;
        ++count;
        p += len;
    }
    if (widest_out) *widest_out = widest;
    return count;
}

#define TIP_SCALE   1
#define TIP_LINE_H  13
#define TIP_PAD     8
#define TIP_MAX_W   268
#define TIP_LINES   8

void draw_tooltip(SDL_Renderer *ren, UiRect anchor, const char *text, int win_w, int win_h) {
    if (!text || !*text) return;

    char lines[TIP_LINES][DRAW_WRAP_LINE_CHARS];
    int widest = 0;
    int n = draw_wrap(TIP_SCALE, TIP_MAX_W, text, lines, TIP_LINES, &widest);
    if (n == 0) return;

    UiRect panel = {
        anchor.x + anchor.w + 10,
        anchor.y,
        widest + 2 * TIP_PAD,
        (n - 1) * TIP_LINE_H + FONT_H * TIP_SCALE + 2 * TIP_PAD,
    };
    /* Keep it on screen: a tooltip that runs off the edge tells you nothing. */
    if (panel.x + panel.w > win_w - 6) panel.x = win_w - 6 - panel.w;
    if (panel.x < 6) panel.x = 6;
    if (panel.y + panel.h > win_h - 6) panel.y = win_h - 6 - panel.h;
    if (panel.y < 6) panel.y = 6;

    draw_rect_fill(ren, panel.x, panel.y, panel.w, panel.h, COL_PANEL, 245);
    draw_rect_line(ren, panel.x, panel.y, panel.w, panel.h, COL_EDGE, 255);
    for (int i = 0; i < n; ++i)
        draw_text(ren, panel.x + TIP_PAD, panel.y + TIP_PAD + i * TIP_LINE_H,
                  TIP_SCALE, lines[i], COL_INK, 255);
}

/* ------------------------------------------------------------- messages --- */
#define STATUS_SCALE  1
#define STATUS_LINE_H 13
#define STATUS_PAD    10

void draw_status(SDL_Renderer *ren, const StatusLog *log, UiRect canvas, unsigned now_ms) {
    /* A condition that is still true goes across the top, where it cannot be
     * mistaken for something that has already happened. */
    if (log->sticky_set) {
        Col c = draw_level_col(log->sticky_level);
        /* The banner is the loudest thing the viewer says, so it is set at the
         * toolbar's own size -- unless the message is too long for the canvas
         * at that size, where being readable beats being loud. */
        int scale = 2;
        if (font_text_width(log->sticky, scale) > canvas.w - 2 * STATUS_PAD) scale = 1;
        draw_rect_fill(ren, canvas.x, canvas.y, canvas.w, DRAW_BANNER_H, COL_PANEL, 235);
        draw_line(ren, canvas.x, canvas.y + DRAW_BANNER_H - 1,
                  canvas.x + canvas.w, canvas.y + DRAW_BANNER_H - 1, c, 255);
        draw_text_mid(ren, canvas.x + canvas.w / 2,
                      canvas.y + (DRAW_BANNER_H - FONT_H * scale) / 2, scale,
                      log->sticky, c, 255);
    }

    /* Everything else stacks up from the bottom-left corner, newest lowest --
     * nearest the mouse, and out of the way of the picture itself. */
    const StatusMessage *shown[STATUS_HISTORY];
    int n = status_visible(log, now_ms, shown, STATUS_HISTORY);
    int max_w = canvas.w - 2 * STATUS_PAD;
    int y = canvas.y + canvas.h - STATUS_PAD - FONT_H * STATUS_SCALE;

    for (int i = 0; i < n; ++i) {
        Col c = draw_level_col(shown[i]->level);
        int alpha = status_alpha(shown[i], now_ms);
        /* Older lines dim further, so the newest reads first. */
        if (i > 0) alpha = (alpha * 55) / 100;
        if (alpha <= 0) continue;

        char lines[DRAW_WRAP_MAX_LINES][DRAW_WRAP_LINE_CHARS];
        int lc = draw_wrap(STATUS_SCALE, max_w, shown[i]->text, lines,
                           DRAW_WRAP_MAX_LINES, NULL);
        y -= (lc - 1) * STATUS_LINE_H;
        if (y < canvas.y + DRAW_BANNER_H) break;

        for (int k = 0; k < lc; ++k)
            draw_text(ren, canvas.x + STATUS_PAD, y + k * STATUS_LINE_H,
                      STATUS_SCALE, lines[k], c, (Uint8)alpha);
        y -= STATUS_LINE_H;
    }
}

/* ----------------------------------------------------------------- help --- */

/* The keys and gestures that have no button, so the help can list them too.
 * Everything else in the overlay is read straight off the toolbar, which is
 * what keeps the two from drifting apart. */
typedef struct { const char *keys, *what; } KeyNote;

static const KeyNote EXTRA_KEYS[] = {
    { "CLICK",       "SELECT THE LIGHT OR PART UNDER THE CURSOR" },
    { "DRAG IT",     "SLIDE THE SELECTION ALONG THE SURFACE BEHIND IT" },
    { "DRAG GIZMO",  "MOVE ALONG ONE AXIS, OR TURN ABOUT IT" },
    { "DRAG EMPTY",  "ORBIT THE VIEW" },
    { "WHEEL",       "ZOOM" },
    { "TAB",         "STEP THROUGH THE LIGHTS, THEN THE PARTS" },
    { "DRAG A ROW",  "SCRUB THAT INSPECTOR VALUE" },
    { "TYPE",        "SET THE FOCUSED ROW, ENTER COMMITS" },
    { "DROP A FILE", "IMPORT OBJ OR STL GEOMETRY" },
    { "F",           "DRAPE THE MEASURED FIELD OVER THE GEOMETRY" },
    { "ESC",         "DISARM, DESELECT, THEN QUIT -- OR CLOSE THIS" },
};
#define EXTRA_KEY_COUNT ((int)(sizeof EXTRA_KEYS / sizeof EXTRA_KEYS[0]))

#define HELP_SCALE   1
#define HELP_ROW_H   15
#define HELP_KEY_W   72      /* the hotkey column, wide enough for "DROP A FILE" */
#define HELP_WHAT_W  280     /* and what it does beside it */
#define HELP_COL_W   (HELP_KEY_W + HELP_WHAT_W)
#define HELP_PAD     20
#define HELP_TOP     44      /* room for the title above the first row */

void draw_help(SDL_Renderer *ren, const Toolbar *t, UiState st, UiRect canvas) {
    int total = t->count + 1 + EXTRA_KEY_COUNT;   /* +1 for the divider */
    /* Two balanced columns, rather than one long one beside empty space. */
    int per_col = (total + 1) / 2;

    UiRect panel = { 0, 0, 2 * HELP_COL_W + 3 * HELP_PAD,
                           HELP_TOP + per_col * HELP_ROW_H + HELP_PAD };
    if (panel.w > canvas.w - 40) panel.w = canvas.w - 40;
    if (panel.h > canvas.h - 32) {
        panel.h = canvas.h - 32;
        per_col = (panel.h - HELP_TOP - HELP_PAD) / HELP_ROW_H;
        if (per_col < 1) per_col = 1;
    }
    if (panel.w < 240 || panel.h < 120) return;
    panel.x = canvas.x + (canvas.w - panel.w) / 2;
    panel.y = canvas.y + (canvas.h - panel.h) / 2;

    draw_rect_fill(ren, panel.x, panel.y, panel.w, panel.h, COL_PANEL, 246);
    draw_rect_line(ren, panel.x, panel.y, panel.w, panel.h, COL_EDGE, 255);
    draw_text_mid(ren, panel.x + panel.w / 2, panel.y + 14, 2,
                  "EVERY COMMAND. H OR ESC CLOSES.", COL_INK, 255);

    int col_w = (panel.w - 3 * HELP_PAD) / 2;
    for (int i = 0; i < total; ++i) {
        int col = i / per_col, row = i % per_col;
        if (col > 1) break;                   /* two columns is all there is room for */
        int x = panel.x + HELP_PAD + col * (col_w + HELP_PAD);
        int y = panel.y + HELP_TOP + row * HELP_ROW_H;

        const char *keys, *what;
        if (i < t->count) {
            /* Read straight off the toolbar, so a command cannot be added
             * without appearing here -- and so a toggle's key is listed
             * against the label it currently shows. */
            keys = t->buttons[i].hint;
            what = ui_label(t, i, st);
        } else if (i == t->count) {
            draw_text(ren, x, y, HELP_SCALE, "MOUSE AND THE REST", COL_HEAD, 255);
            continue;
        } else {
            const KeyNote *k = &EXTRA_KEYS[i - t->count - 1];
            keys = k->keys;
            what = k->what;
        }
        draw_text(ren, x, y, HELP_SCALE, keys, COL_ACCENT, 255);
        draw_text(ren, x + HELP_KEY_W, y, HELP_SCALE, what, COL_INK, 255);
    }
}

/* ------------------------------------------------------------- readouts --- */

void draw_colorbar(SDL_Renderer *ren, int x, int y, int w, int h,
                   double lo, double hi, const char *unit) {
    for (int i = 0; i < w; ++i) {
        Uint8 r, g, b;
        draw_viridis((double)i / (w > 1 ? w - 1 : 1), &r, &g, &b);
        SDL_SetRenderDrawColor(ren, r, g, b, 255);
        SDL_Rect px = { x + i, y, 1, h };
        SDL_RenderFillRect(ren, &px);
    }
    draw_rect_line(ren, x, y, w, h, COL_RULE, 255);
    char buf[64];
    for (int k = 0; k <= 4; ++k) {
        double v = lo + (hi - lo) * k / 4.0;
        if (v >= 100.0) snprintf(buf, sizeof buf, "%.0f", v);
        else if (v >= 1.0) snprintf(buf, sizeof buf, "%.1f", v);
        else snprintf(buf, sizeof buf, "%.2f", v);
        int tx = x + (w - 1) * k / 4;
        if (k == 0)      draw_text(ren, tx, y + h + 5, 1, buf, COL_MUTED, 255);
        else if (k == 4) draw_text_right(ren, tx, y + h + 5, 1, buf, COL_MUTED, 255);
        else             draw_text_mid(ren, tx, y + h + 5, 1, buf, COL_MUTED, 255);
    }
    draw_text_right(ren, x + w, y - 10, 1, unit, COL_MUTED, 255);
}

void draw_plot(SDL_Renderer *ren, int x, int y, int w, int h,
               const double *full, const double *direct, int n,
               double ymax, const char *unit, const char *caption) {
    draw_rect_fill(ren, x, y, w, h, COL_PANEL, 255);
    draw_rect_line(ren, x, y, w, h, COL_RULE, 255);

    const int ml = 54, mr = 12, mt = 20, mb = 20;
    int px0 = x + ml, py0 = y + mt;
    int pw = w - ml - mr, ph = h - mt - mb;
    if (pw <= 4 || ph <= 4 || n < 2) return;
    if (ymax <= 0.0) ymax = 1.0;

    char buf[64];
    for (int k = 0; k <= 3; ++k) {
        double v = ymax * k / 3.0;
        int gy = py0 + ph - ph * k / 3;
        draw_line(ren, px0, gy, px0 + pw, gy, COL_RULE, 255);
        if (v >= 100.0) snprintf(buf, sizeof buf, "%.0f", v);
        else snprintf(buf, sizeof buf, "%.2f", v);
        draw_text_right(ren, px0 - 6, gy - 3, 1, buf, COL_MUTED, 255);
    }

    for (int pass = 0; pass < 2; ++pass) {
        const double *s = pass ? full : direct;
        Col c = pass ? COL_S1 : COL_S2;
        if (!s) continue;
        SDL_SetRenderDrawColor(ren, c.r, c.g, c.b, 255);
        for (int i = 1; i < n; ++i) {
            int x0 = px0 + (i - 1) * pw / (n - 1);
            int x1 = px0 + i * pw / (n - 1);
            int y0 = py0 + ph - (int)(ph * (s[i-1] / ymax));
            int y1 = py0 + ph - (int)(ph * (s[i]   / ymax));
            SDL_RenderDrawLine(ren, x0, y0, x1, y1);
            SDL_RenderDrawLine(ren, x0, y0 + 1, x1, y1 + 1);   /* 2px stroke */
        }
    }

    draw_text(ren, x + 10, y + 6, 1, caption, COL_HEAD, 255);
    draw_text_right(ren, x + w - 10, y + 6, 1, unit, COL_MUTED, 255);
    /* Legend, so the two series are never identified by colour alone. */
    draw_rect_fill(ren, x + w - 150, y + h - 14, 12, 2, COL_S1, 255);
    draw_text(ren, x + w - 134, y + h - 18, 1, "FULL", COL_MUTED, 255);
    draw_rect_fill(ren, x + w - 92, y + h - 14, 12, 2, COL_S2, 255);
    draw_text(ren, x + w - 76, y + h - 18, 1, "DIRECT", COL_MUTED, 255);
}
