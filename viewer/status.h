/* status.h — the viewer's message channel.
 *
 * Everything the viewer has to say -- a refused command, a file written, an
 * import that failed -- goes through here so it can be shown ON THE CANVAS.
 * It used to go to stdout plus one transient line, which is invisible to
 * anyone who launched the app from the Finder, or whose window covers the
 * terminal it was started from.
 *
 * Two kinds of message, because there are two kinds of thing to say:
 *   - a moment that has passed ("wrote out/render.ppm"), which fades out;
 *   - a condition that is still true ("solving the field"), which stays up
 *     until whatever caused it goes away.
 *
 * Free of any SDL dependency so the rules can be exercised headlessly;
 * drawing lives in draw.c (draw_status).
 */
#ifndef LIGHTSIM_VIEWER_STATUS_H
#define LIGHTSIM_VIEWER_STATUS_H

#include <stdbool.h>

#define STATUS_MAX_TEXT       192
#define STATUS_HISTORY        3       /* how many recent lines stay on screen */
#define STATUS_FADE_MS        6000u   /* how long one stays before fading */
#define STATUS_FADE_TAIL_MS   900u    /* the fading part at the end of that */

typedef enum { STATUS_INFO, STATUS_WARN, STATUS_ERROR } StatusLevel;

typedef struct {
    char        text[STATUS_MAX_TEXT];
    StatusLevel level;
    unsigned    stamp_ms;
} StatusMessage;

typedef struct {
    StatusMessage recent[STATUS_HISTORY];   /* oldest first; [count-1] newest */
    int           count;

    char        sticky[STATUS_MAX_TEXT];
    StatusLevel sticky_level;
    bool        sticky_set;
} StatusLog;

void status_init(StatusLog *log);

/* Adds a message stamped `now_ms`. Repeating the newest message only restamps
 * it, so holding a key that refuses does not scroll the same line three
 * times. */
void status_push(StatusLog *log, StatusLevel level, unsigned now_ms, const char *text);

/* The condition banner: set while it holds, cleared when it stops holding. */
void status_set_sticky(StatusLog *log, StatusLevel level, const char *text);
void status_clear_sticky(StatusLog *log);

/* Messages not yet faded out, NEWEST FIRST. Returns how many were written. */
int status_visible(const StatusLog *log, unsigned now_ms,
                   const StatusMessage **out, int max_out);

/* 0..255 alpha for a message at `now_ms`: solid, then fading over the last
 * STATUS_FADE_TAIL_MS, then zero. */
int status_alpha(const StatusMessage *msg, unsigned now_ms);

#endif /* LIGHTSIM_VIEWER_STATUS_H */
