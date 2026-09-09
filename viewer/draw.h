/* draw.h — SDL drawing for the viewer: text, panels, toolbar, colour bar, plot. */
#ifndef LIGHTSIM_VIEWER_DRAW_H
#define LIGHTSIM_VIEWER_DRAW_H

#include <SDL2/SDL.h>
#include "ui.h"
#include "status.h"

/* Instrument palette. Neutral graphite, with ONE accent -- the selection
 * yellow -- so "this is the thing you are acting on" is never said two
 * different ways. Dark by design: a false-colour ramp reads better against a
 * dark ground, and the viewer is looked at for long stretches. */
typedef struct { Uint8 r, g, b; } Col;

/* Grounds and rules. */
extern const Col COL_BG, COL_PANEL, COL_RULE;
/* Button faces, by state. */
extern const Col COL_FACE, COL_FACE_HOT, COL_FACE_DOWN, COL_FACE_OFF;
extern const Col COL_EDGE, COL_EDGE_OFF;
/* Text, brightest to dimmest. */
extern const Col COL_INK, COL_HEAD, COL_MUTED, COL_DIM;
/* The accent, and its dimmed companion for an active button's hotkey hint. */
extern const Col COL_ACCENT, COL_ACCENT_DIM;
/* Message levels, and the two plot series. */
extern const Col COL_WARN, COL_ERR, COL_S1, COL_S2;
/* Gizmo axis colours, taken from the validated categorical palette so the three
 * stay distinguishable to a colour-blind reader as well. */
extern const Col COL_AXIS_X, COL_AXIS_Y, COL_AXIS_Z;

/* The colour a message of `level` is drawn in. */
Col draw_level_col(StatusLevel level);

void draw_rect_fill(SDL_Renderer *ren, int x, int y, int w, int h, Col c, Uint8 a);
void draw_rect_line(SDL_Renderer *ren, int x, int y, int w, int h, Col c, Uint8 a);
void draw_line(SDL_Renderer *ren, int x0, int y0, int x1, int y1, Col c, Uint8 a);
/* A small cross-in-a-box marker, for objects that project to a single point. */
void draw_marker(SDL_Renderer *ren, int x, int y, int r, Col c, Uint8 a);
void draw_text(SDL_Renderer *ren, int x, int y, int scale, const char *s, Col c, Uint8 a);
void draw_text_right(SDL_Renderer *ren, int right_x, int y, int scale, const char *s, Col c, Uint8 a);
void draw_text_mid(SDL_Renderer *ren, int mid_x, int y, int scale, const char *s, Col c, Uint8 a);

void draw_toolbar(SDL_Renderer *ren, const Toolbar *t, UiState st, int strip_h);

/* ---- wrapped text, shared by the tooltip, the status log and the help ---- */
#define DRAW_WRAP_MAX_LINES 10
#define DRAW_WRAP_LINE_CHARS 96

/* Breaks `text` into lines no wider than `max_w` at `scale`, on word
 * boundaries. Returns how many lines were written, and the widest through
 * `widest_out` when that is not NULL. */
int draw_wrap(int scale, int max_w, const char *text,
              char lines[][DRAW_WRAP_LINE_CHARS], int max_lines, int *widest_out);

/* Hover text for one toolbar button, parked beside `anchor` and kept on
 * screen. Draws nothing for an empty tip. */
void draw_tooltip(SDL_Renderer *ren, UiRect anchor, const char *text, int win_w, int win_h);

/* A condition that still holds goes across the top of `canvas`; everything
 * that has merely happened stacks up from its bottom-left corner, newest
 * lowest and older lines dimmer. */
#define DRAW_BANNER_H 26
void draw_status(SDL_Renderer *ren, const StatusLog *log, UiRect canvas, unsigned now_ms);

/* Every command in one place, generated from the toolbar itself plus the keys
 * and gestures that have no button. */
void draw_help(SDL_Renderer *ren, const Toolbar *t, UiState st, UiRect canvas);

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
