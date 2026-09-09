#include "status.h"

#include <stdio.h>
#include <string.h>

void status_init(StatusLog *log) {
    log->count = 0;
    log->sticky[0] = '\0';
    log->sticky_level = STATUS_INFO;
    log->sticky_set = false;
    for (int i = 0; i < STATUS_HISTORY; ++i) {
        log->recent[i].text[0] = '\0';
        log->recent[i].level = STATUS_INFO;
        log->recent[i].stamp_ms = 0;
    }
}

void status_push(StatusLog *log, StatusLevel level, unsigned now_ms, const char *text) {
    if (!text || !*text) return;

    /* The same thing said twice running is one thing that is still true, not
     * two things that happened -- restamp it rather than repeating it. */
    if (log->count > 0) {
        StatusMessage *newest = &log->recent[log->count - 1];
        if (newest->level == level && strncmp(newest->text, text, STATUS_MAX_TEXT) == 0) {
            newest->stamp_ms = now_ms;
            return;
        }
    }

    if (log->count >= STATUS_HISTORY) {
        for (int i = 1; i < STATUS_HISTORY; ++i) log->recent[i - 1] = log->recent[i];
        log->count = STATUS_HISTORY - 1;
    }
    StatusMessage *m = &log->recent[log->count++];
    snprintf(m->text, sizeof m->text, "%s", text);
    m->level = level;
    m->stamp_ms = now_ms;
}

void status_set_sticky(StatusLog *log, StatusLevel level, const char *text) {
    if (!text || !*text) { status_clear_sticky(log); return; }
    snprintf(log->sticky, sizeof log->sticky, "%s", text);
    log->sticky_level = level;
    log->sticky_set = true;
}

void status_clear_sticky(StatusLog *log) {
    log->sticky[0] = '\0';
    log->sticky_set = false;
}

int status_alpha(const StatusMessage *msg, unsigned now_ms) {
    /* Unsigned subtraction, so this stays right across the SDL_GetTicks wrap. */
    unsigned age = now_ms - msg->stamp_ms;
    if (age >= STATUS_FADE_MS) return 0;
    unsigned solid = STATUS_FADE_MS - STATUS_FADE_TAIL_MS;
    if (age <= solid) return 255;
    unsigned into_tail = age - solid;
    return (int)(255u - (into_tail * 255u) / STATUS_FADE_TAIL_MS);
}

int status_visible(const StatusLog *log, unsigned now_ms,
                   const StatusMessage **out, int max_out) {
    int n = 0;
    for (int i = log->count - 1; i >= 0 && n < max_out; --i) {
        if (status_alpha(&log->recent[i], now_ms) <= 0) continue;
        out[n++] = &log->recent[i];
    }
    return n;
}
