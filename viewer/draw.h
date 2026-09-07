/* draw.h — SDL drawing for the viewer: text, panels, toolbar, colour bar, plot. */
#ifndef LIGHTSIM_VIEWER_DRAW_H
#define LIGHTSIM_VIEWER_DRAW_H

#include <SDL2/SDL.h>
#include "ui.h"

/* Instrument palette. Dark by design: a false-colour ramp reads better against
 * a dark ground, and the viewer is looked at for long stretches. */
typedef struct { Uint8 r, g, b; } Col;
extern const Col COL_BG, COL_PANEL, COL_RULE, COL_INK, COL_MUTED,
                 COL_ACCENT, COL_S1, COL_S2, COL_WARN;

void draw_rect_fill(SDL_Renderer *ren, int x, int y, int w, int h, Col c, Uint8 a);
void draw_rect_line(SDL_Renderer *ren, int x, int y, int w, int h, Col c, Uint8 a);
void draw_line(SDL_Renderer *ren, int x0, int y0, int x1, int y1, Col c, Uint8 a);
/* A small cross-in-a-box marker, for objects that project to a single point. */
void draw_marker(SDL_Renderer *ren, int x, int y, int r, Col c, Uint8 a);
void draw_text(SDL_Renderer *ren, int x, int y, int scale, const char *s, Col c, Uint8 a);
void draw_text_right(SDL_Renderer *ren, int right_x, int y, int scale, const char *s, Col c, Uint8 a);

void draw_toolbar(SDL_Renderer *ren, const Toolbar *t, UiState st, int strip_h);

/* Horizontal viridis ramp with `lo`..`hi` tick labels beneath it. */
void draw_colorbar(SDL_Renderer *ren, int x, int y, int w, int h,
                   double lo, double hi, const char *unit);

/* Two-series cross-section. `n` samples each; NULL series is skipped. */
void draw_plot(SDL_Renderer *ren, int x, int y, int w, int h,
               const double *full, const double *direct, int n,
               double ymax, const char *unit, const char *caption);

/* Viridis lookup, matching the ramp the PPM writer and the Blender export use. */
void draw_viridis(double t, Uint8 *r, Uint8 *g, Uint8 *b);

#endif /* LIGHTSIM_VIEWER_DRAW_H */
