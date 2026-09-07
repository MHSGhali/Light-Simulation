#include "draw.h"
#include "font.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

const Col COL_BG     = { 0x12, 0x18, 0x1A };
const Col COL_PANEL  = { 0x1A, 0x22, 0x24 };
const Col COL_RULE   = { 0x2A, 0x35, 0x38 };
const Col COL_INK    = { 0xE7, 0xED, 0xEC };
const Col COL_MUTED  = { 0x8B, 0x9B, 0x9A };
const Col COL_ACCENT = { 0x54, 0xC3, 0xB4 };
const Col COL_S1     = { 0x39, 0x87, 0xE5 };
const Col COL_S2     = { 0xD9, 0x59, 0x26 };
const Col COL_WARN   = { 0xED, 0xA1, 0x00 };

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

void draw_toolbar(SDL_Renderer *ren, const Toolbar *t, UiState st, int strip_h) {
    draw_rect_fill(ren, 0, 0, UI_TOOLBAR_W, strip_h, COL_PANEL, 255);
    SDL_SetRenderDrawColor(ren, COL_RULE.r, COL_RULE.g, COL_RULE.b, 255);
    SDL_RenderDrawLine(ren, UI_TOOLBAR_W - 1, 0, UI_TOOLBAR_W - 1, strip_h);

    for (int i = 0; i < t->count; ++i) {
        const UiButton *b = &t->buttons[i];
        if (b->separator_above) {
            SDL_SetRenderDrawColor(ren, COL_RULE.r, COL_RULE.g, COL_RULE.b, 255);
            int sy = b->rect.y - UI_GROUP_GAP / 2 - UI_BUTTON_GAP / 2;
            SDL_RenderDrawLine(ren, UI_MARGIN, sy, UI_MARGIN + UI_BUTTON_W, sy);
        }
        bool hot = (t->hover == i && b->enabled);
        Col fill = b->active ? COL_ACCENT : (hot ? COL_RULE : COL_PANEL);
        Uint8 alpha = b->enabled ? 255 : 90;
        draw_rect_fill(ren, b->rect.x, b->rect.y, b->rect.w, b->rect.h, fill, alpha);
        draw_rect_line(ren, b->rect.x, b->rect.y, b->rect.w, b->rect.h,
                       b->active ? COL_ACCENT : COL_RULE, alpha);
        Col ink = b->active ? COL_BG : COL_INK;
        draw_text(ren, b->rect.x + 9, b->rect.y + 10, 2, ui_label(t, i, st), ink,
                  b->enabled ? 255 : 110);
        if (b->hint)
            draw_text_right(ren, b->rect.x + b->rect.w - 8, b->rect.y + 11, 1,
                            b->hint, b->active ? COL_BG : COL_MUTED,
                            b->enabled ? 200 : 90);
    }
}

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
        else             draw_text(ren, tx - font_text_width(buf, 1) / 2, y + h + 5, 1,
                                   buf, COL_MUTED, 255);
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
        int gy = py0 + ph - (int)(ph * k / 3.0);
        SDL_SetRenderDrawColor(ren, COL_RULE.r, COL_RULE.g, COL_RULE.b, 255);
        SDL_RenderDrawLine(ren, px0, gy, px0 + pw, gy);
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

    draw_text(ren, x + 10, y + 6, 1, caption, COL_MUTED, 255);
    draw_text_right(ren, x + w - 10, y + 6, 1, unit, COL_MUTED, 255);
    /* Legend, so the two series are never identified by colour alone. */
    draw_rect_fill(ren, x + w - 150, y + h - 14, 12, 2, COL_S1, 255);
    draw_text(ren, x + w - 134, y + h - 18, 1, "FULL", COL_MUTED, 255);
    draw_rect_fill(ren, x + w - 92, y + h - 14, 12, 2, COL_S2, 255);
    draw_text(ren, x + w - 76, y + h - 18, 1, "DIRECT", COL_MUTED, 255);
}
