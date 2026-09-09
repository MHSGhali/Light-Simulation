#include "capture.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* One line of the script. `frames` is how many presented frames the step
 * occupies; 0 means it takes effect between frames and costs no time, which
 * keeps `record` from adding a stray leading frame to a recording. */
typedef enum {
    STEP_WAIT, STEP_KEY, STEP_MOVE, STEP_CLICK, STEP_DRAG, STEP_DROP,
    STEP_SHOT, STEP_RECORD, STEP_STOP
} StepKind;

typedef struct {
    StepKind    kind;
    int         frames;
    int         x0, y0, x1, y1;
    SDL_Keycode key;
    Uint16      mod;
    char        name[256];      /* shot stem, record prefix, or dropped path */
} Step;

struct Capture {
    Step  *steps;
    int    nsteps;
    int    step;            /* the step being executed */
    int    frame;           /* frames already spent on it */

    char   out_dir[512];
    char   rec[128];        /* recording prefix; empty when not recording */
    int    rec_every;       /* frames between written frames, >= 1 */
    int    rec_countdown;
    int    rec_n;
    char   shot[128];       /* file for THIS frame; empty when there is none */

    bool   mod_held;        /* a modifier is set for a key event still in flight */
    Uint32 frame_start_ms;
    int    written;
};

/* Milliseconds as the script writes them, rounded to whole frames. Never
 * zero: a step asking for 5 ms still wants to happen. */
static int ms_to_frames(int ms) {
    int f = (ms + CAPTURE_FRAME_MS / 2) / CAPTURE_FRAME_MS;
    return f < 1 ? 1 : f;
}

/* ------------------------------------------------------------- parsing --- */

/* Splits on whitespace, in place, stopping at a '#'. Returns the token count. */
static int tokenize(char *line, char **tok, int max_tok) {
    int n = 0;
    char *p = line;
    while (*p && n < max_tok) {
        while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') ++p;
        if (!*p || *p == '#') break;
        tok[n++] = p;
        while (*p && *p != ' ' && *p != '\t' && *p != '\r' && *p != '\n') ++p;
        if (*p) *p++ = '\0';
    }
    return n;
}

/* Modifier prefixes on a `key` line. cmd and ctrl are one thing here: the
 * viewer accepts either for undo, so a script should not have to know which
 * platform it is describing. */
static bool parse_mod(const char *word, Uint16 *mod) {
    if (!strcmp(word, "cmd") || !strcmp(word, "ctrl")) { *mod |= KMOD_LGUI | KMOD_LCTRL; return true; }
    if (!strcmp(word, "shift")) { *mod |= KMOD_LSHIFT; return true; }
    return false;
}

static bool parse_step(char **tok, int n, Step *s, char *err, size_t errn) {
    memset(s, 0, sizeof *s);
    const char *cmd = tok[0];

    if (!strcmp(cmd, "wait") && n == 2) {
        s->kind = STEP_WAIT;
        s->frames = ms_to_frames(atoi(tok[1]));
        return true;
    }
    if (!strcmp(cmd, "key") && n >= 2) {
        int i = 1;
        while (i < n - 1 && parse_mod(tok[i], &s->mod)) ++i;
        if (i != n - 1) {
            snprintf(err, errn, "key: expected modifiers then one key name");
            return false;
        }
        /* SDL's own key names, so the script says "Tab" and "Escape" rather
         * than carrying a second table that can drift from SDL's. */
        s->key = SDL_GetKeyFromName(tok[i]);
        if (s->key == SDLK_UNKNOWN) {
            snprintf(err, errn, "key: SDL does not know a key named '%s'", tok[i]);
            return false;
        }
        s->kind = STEP_KEY;
        s->frames = 1;
        return true;
    }
    if (!strcmp(cmd, "move") && n == 3) {
        s->kind = STEP_MOVE;
        s->x0 = s->x1 = atoi(tok[1]);
        s->y0 = s->y1 = atoi(tok[2]);
        s->frames = 1;
        return true;
    }
    if (!strcmp(cmd, "click") && n == 3) {
        s->kind = STEP_CLICK;
        s->x0 = s->x1 = atoi(tok[1]);
        s->y0 = s->y1 = atoi(tok[2]);
        s->frames = 1;
        return true;
    }
    if (!strcmp(cmd, "drag") && n == 6) {
        s->kind = STEP_DRAG;
        s->x0 = atoi(tok[1]); s->y0 = atoi(tok[2]);
        s->x1 = atoi(tok[3]); s->y1 = atoi(tok[4]);
        s->frames = ms_to_frames(atoi(tok[5]));
        return true;
    }
    if (!strcmp(cmd, "drop") && n == 2) {
        s->kind = STEP_DROP;
        snprintf(s->name, sizeof s->name, "%s", tok[1]);
        s->frames = 1;
        return true;
    }
    if (!strcmp(cmd, "shot") && n == 2) {
        s->kind = STEP_SHOT;
        snprintf(s->name, sizeof s->name, "%s", tok[1]);
        s->frames = 1;
        return true;
    }
    if (!strcmp(cmd, "record") && (n == 2 || n == 3)) {
        s->kind = STEP_RECORD;
        snprintf(s->name, sizeof s->name, "%s", tok[1]);
        /* A gif wants tens of frames, not the hundreds a ten-second sequence
         * would produce at the capture rate -- and each one is a full-window
         * PPM on disk. The stride is in milliseconds so it also states the
         * playback rate the frames were sampled at. */
        s->x0 = (n == 3) ? ms_to_frames(atoi(tok[2])) : 1;
        s->frames = 0;
        return true;
    }
    if (!strcmp(cmd, "stop") && n == 1) {
        s->kind = STEP_STOP;
        s->frames = 0;
        return true;
    }
    snprintf(err, errn, "unknown command '%s' with %d argument(s)", cmd, n - 1);
    return false;
}

Capture *capture_open(const char *script_path, const char *out_dir,
                      char *err, size_t errn) {
    FILE *f = fopen(script_path, "r");
    if (!f) {
        snprintf(err, errn, "cannot open %s", script_path);
        return NULL;
    }
    Capture *c = calloc(1, sizeof *c);
    if (!c) { fclose(f); snprintf(err, errn, "out of memory"); return NULL; }
    snprintf(c->out_dir, sizeof c->out_dir, "%s", out_dir);

    int cap = 64;
    c->steps = calloc((size_t)cap, sizeof *c->steps);
    if (!c->steps) { free(c); fclose(f); snprintf(err, errn, "out of memory"); return NULL; }

    char line[512];
    int lineno = 0;
    while (fgets(line, sizeof line, f)) {
        char *tok[8];
        ++lineno;
        int n = tokenize(line, tok, 8);
        if (n == 0) continue;

        if (c->nsteps == cap) {
            int grown = cap * 2;
            Step *bigger = realloc(c->steps, (size_t)grown * sizeof *c->steps);
            if (!bigger) { capture_free(c); fclose(f); snprintf(err, errn, "out of memory"); return NULL; }
            c->steps = bigger;
            cap = grown;
        }
        char why[160];
        if (!parse_step(tok, n, &c->steps[c->nsteps], why, sizeof why)) {
            snprintf(err, errn, "%s:%d: %s", script_path, lineno, why);
            capture_free(c);
            fclose(f);
            return NULL;
        }
        ++c->nsteps;
    }
    fclose(f);
    if (c->nsteps == 0) {
        snprintf(err, errn, "%s is empty", script_path);
        capture_free(c);
        return NULL;
    }
    return c;
}

void capture_free(Capture *c) {
    if (!c) return;
    free(c->steps);
    free(c);
}

int capture_written(const Capture *c) { return c ? c->written : 0; }

/* -------------------------------------------------------------- input --- */

static void push_motion(int x, int y) {
    SDL_Event e;
    memset(&e, 0, sizeof e);
    e.type = SDL_MOUSEMOTION;
    e.motion.x = x;
    e.motion.y = y;
    SDL_PushEvent(&e);
}

static void push_button(int x, int y, bool down) {
    SDL_Event e;
    memset(&e, 0, sizeof e);
    e.type = down ? SDL_MOUSEBUTTONDOWN : SDL_MOUSEBUTTONUP;
    e.button.button = SDL_BUTTON_LEFT;
    e.button.state = down ? SDL_PRESSED : SDL_RELEASED;
    e.button.clicks = 1;
    e.button.x = x;
    e.button.y = y;
    SDL_PushEvent(&e);
}

/* SDL owns the string on a real drop and the handler frees it with SDL_free,
 * so a scripted drop has to allocate the same way. */
static void push_drop(const char *path) {
    SDL_Event e;
    memset(&e, 0, sizeof e);
    e.type = SDL_DROPFILE;
    e.drop.file = SDL_strdup(path);
    if (!e.drop.file) return;
    if (SDL_PushEvent(&e) <= 0) SDL_free(e.drop.file);
}

static void push_key(SDL_Keycode key, Uint16 mod) {
    SDL_Event e;
    memset(&e, 0, sizeof e);
    e.type = SDL_KEYDOWN;
    e.key.state = SDL_PRESSED;
    e.key.keysym.sym = key;
    e.key.keysym.scancode = SDL_GetScancodeFromKey(key);
    e.key.keysym.mod = mod;
    SDL_PushEvent(&e);
}

/* The undo path asks SDL for the live modifier state rather than reading the
 * event, which a pushed event cannot answer. Setting it globally does, and it
 * is cleared on the next frame -- by then the event has been polled. */
static void hold_mods(Capture *c, Uint16 mod) {
    if (!mod) return;
    SDL_SetModState((SDL_Keymod)mod);
    c->mod_held = true;
}

static int lerp_i(int a, int b, int i, int n) {
    if (n <= 1) return b;
    return a + (b - a) * i / (n - 1);
}

bool capture_begin_frame(Capture *c) {
    c->frame_start_ms = SDL_GetTicks();
    c->shot[0] = '\0';
    if (c->mod_held) { SDL_SetModState(KMOD_NONE); c->mod_held = false; }

    /* Steps that cost no time run here and hand over to the next one, so a
     * `record`/`stop` pair brackets exactly the frames between them. */
    for (;;) {
        if (c->step >= c->nsteps) return false;
        Step *s = &c->steps[c->step];
        if (s->frames != 0) break;
        if (s->kind == STEP_RECORD) {
            snprintf(c->rec, sizeof c->rec, "%s", s->name);
            c->rec_every = s->x0 > 0 ? s->x0 : 1;
            c->rec_countdown = 0;
            c->rec_n = 0;
        } else if (s->kind == STEP_STOP) {
            c->rec[0] = '\0';
        }
        ++c->step;
    }

    Step *s = &c->steps[c->step];
    int i = c->frame, last = s->frames - 1;

    switch (s->kind) {
    case STEP_KEY:
        push_key(s->key, s->mod);
        hold_mods(c, s->mod);
        break;
    case STEP_MOVE:
        push_motion(s->x0, s->y0);
        break;
    case STEP_CLICK:
        push_motion(s->x0, s->y0);
        push_button(s->x0, s->y0, true);
        push_button(s->x0, s->y0, false);
        break;
    case STEP_DRAG:
        /* The press lands on the start point before any movement, because the
         * viewer decides between "orbit" and "move the selection" from what
         * was under the cursor when the button went down. */
        if (i == 0) {
            push_motion(s->x0, s->y0);
            push_button(s->x0, s->y0, true);
        }
        push_motion(lerp_i(s->x0, s->x1, i, s->frames),
                    lerp_i(s->y0, s->y1, i, s->frames));
        if (i == last) push_button(s->x1, s->y1, false);
        break;
    case STEP_DROP:
        push_drop(s->name);
        break;
    case STEP_SHOT:
        snprintf(c->shot, sizeof c->shot, "%s", s->name);
        break;
    case STEP_WAIT:
    case STEP_RECORD:
    case STEP_STOP:
        break;
    }

    if (++c->frame >= s->frames) { c->frame = 0; ++c->step; }
    return true;
}

/* ------------------------------------------------------------- output --- */

static void write_ppm(Capture *c, SDL_Renderer *ren, const char *stem) {
    int w = 0, h = 0;
    if (SDL_GetRendererOutputSize(ren, &w, &h) != 0 || w <= 0 || h <= 0) {
        fprintf(stderr, "capture: renderer size: %s\n", SDL_GetError());
        return;
    }
    size_t pitch = (size_t)w * 3u;
    unsigned char *px = malloc(pitch * (size_t)h);
    if (!px) return;
    if (SDL_RenderReadPixels(ren, NULL, SDL_PIXELFORMAT_RGB24, px, (int)pitch) != 0) {
        fprintf(stderr, "capture: read pixels: %s\n", SDL_GetError());
        free(px);
        return;
    }
    char path[640];
    snprintf(path, sizeof path, "%s/%s.ppm", c->out_dir, stem);
    FILE *f = fopen(path, "wb");
    if (!f) {
        fprintf(stderr, "capture: cannot write %s\n", path);
        free(px);
        return;
    }
    fprintf(f, "P6\n%d %d\n255\n", w, h);
    fwrite(px, 1, pitch * (size_t)h, f);
    fclose(f);
    free(px);
    ++c->written;
}

void capture_end_frame(Capture *c, SDL_Renderer *ren) {
    if (c->shot[0]) write_ppm(c, ren, c->shot);
    if (c->rec[0]) {
        if (c->rec_countdown <= 0) {
            char stem[160];
            snprintf(stem, sizeof stem, "%s_%04d", c->rec, c->rec_n++);
            write_ppm(c, ren, stem);
            c->rec_countdown = c->rec_every;
        }
        --c->rec_countdown;
    }
    /* Pace to a fixed frame, so "wait 2000" is two seconds of settling here
     * and two seconds on any other machine. Writing a frame can overrun it,
     * in which case the next frame simply starts late rather than compounding. */
    Uint32 spent = SDL_GetTicks() - c->frame_start_ms;
    if (spent < CAPTURE_FRAME_MS) SDL_Delay(CAPTURE_FRAME_MS - spent);
}
