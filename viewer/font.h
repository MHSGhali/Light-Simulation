/* font.h — a 5x7 bitmap font covering ASCII 32..95 (space through underscore,
 * which includes the digits and uppercase). Lowercase is upcased on the way in.
 *
 * Deliberately not SDL_ttf: the viewer needs six words and a lot of numbers,
 * and a table costs less than a font dependency. Kept free of SDL so the glyph
 * data can be checked headlessly.
 */
#ifndef LIGHTSIM_FONT_H
#define LIGHTSIM_FONT_H

#include <stdbool.h>

#define FONT_W 5
#define FONT_H 7

/* Column bitmap for `c`: 5 bytes, bit 0 is the top row. Returns the blank
 * glyph for anything outside the covered range. */
const unsigned char *font_glyph(char c);

/* True if the pixel at (col,row) of `c`'s glyph is set. */
bool font_pixel(char c, int col, int row);

/* Advance width of `text` at `scale`, including the 1px inter-glyph gap. */
int font_text_width(const char *text, int scale);

#endif /* LIGHTSIM_FONT_H */
