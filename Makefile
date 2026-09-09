# lightsim — spectral radiometric lighting simulator
# GNU Make 3.81 compatible (macOS ships 3.81): no .ONESHELL, no $(file ...).

CC      := cc
CSTD    := -std=c11 -D_DEFAULT_SOURCE
WARN    := -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wdouble-promotion \
           -Wstrict-prototypes -Wmissing-prototypes

# -ffp-contract=off keeps FMA contraction from shifting results at the 1e-9
# level, which would make the deterministic tests flaky across -O levels.
# -fno-fast-math is policy, not preference: fast-math would license
# reassociation of the spectral sums and disable the NaN/Inf semantics the
# estimators rely on to fail loudly.
FPFLAGS := -ffp-contract=off -fno-fast-math

# Each configuration builds into its own directory. Sharing one caused
# sanitizer-instrumented objects to be linked into a plain build.
VARIANT ?= release
ifeq ($(VARIANT),release)
  OPT  := -O2 $(FPFLAGS)
  SAN  :=
else ifeq ($(VARIANT),debug)
  OPT  := -O0 -g3 $(FPFLAGS)
  SAN  :=
else ifeq ($(VARIANT),asan)
  OPT  := -O1 -g -fno-omit-frame-pointer $(FPFLAGS)
  SAN  := -fsanitize=address,undefined -fno-sanitize-recover=all
else ifeq ($(VARIANT),ubsan)
  OPT  := -O1 -g -fno-omit-frame-pointer $(FPFLAGS)
  SAN  := -fsanitize=undefined -fno-sanitize-recover=all
endif

# -MMD -MP emits a .d file of header prerequisites next to every object.
# Without this a change to a struct in a header recompiles only the .c files
# that changed, and the rest of the build keeps the OLD struct layout -- which
# does not fail to link, it just reads garbage at runtime.
CFLAGS  := $(CSTD) $(WARN) $(OPT) $(SAN) -Iinclude -MMD -MP
LDLIBS  := -lm -lpthread

BUILD   := build/$(VARIANT)
SRC     := $(wildcard src/*.c)
OBJ     := $(patsubst src/%.c,$(BUILD)/%.o,$(SRC))
TESTSRC := $(wildcard tests/*.c)
TESTOBJ := $(patsubst tests/%.c,$(BUILD)/tests_%.o,$(TESTSRC))
APPSRC  := $(wildcard apps/*.c)
APPOBJ  := $(patsubst apps/%.c,$(BUILD)/apps_%.o,$(APPSRC))

.PHONY: all test debug clean test-asan test-ubsan check-purity view

# One build command. The viewer is included whenever SDL2 is available; without
# it the library, CLI and tests still build exactly as before.
SDL_PROBE := $(shell pkg-config --exists sdl2 2>/dev/null && echo yes)
ifeq ($(SDL_PROBE),yes)
  DEFAULT_TARGETS := lightsim lightsim-view
else
  DEFAULT_TARGETS := lightsim
endif

all: $(DEFAULT_TARGETS)
	@if [ "$(SDL_PROBE)" != "yes" ]; then \
	    echo "note: SDL2 not found, so the interactive viewer was skipped."; \
	    echo "      brew install sdl2   then re-run make"; \
	fi

lightsim: $(OBJ) $(APPOBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(APPOBJ) $(LDLIBS)

$(BUILD)/%.o: src/%.c | $(BUILD)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILD)/tests_%.o: tests/%.c | $(BUILD)
	$(CC) $(CFLAGS) -Itests -c -o $@ $<

$(BUILD)/apps_%.o: apps/%.c | $(BUILD)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILD):
	mkdir -p $(BUILD)

# ---- interactive viewer (optional: needs SDL2) ----------------------------
# Kept out of `all` and out of $(SRC) so the library, the CLI and the tests all
# build with no SDL2 installed.
SDL_CFLAGS := $(shell pkg-config --cflags sdl2 2>/dev/null)
SDL_LIBS   := $(shell pkg-config --libs sdl2 2>/dev/null)
VIEWSRC    := $(wildcard viewer/*.c)
VIEWOBJ    := $(patsubst viewer/%.c,$(BUILD)/viewer_%.o,$(VIEWSRC))

$(BUILD)/viewer_%.o: viewer/%.c | $(BUILD)
	$(CC) $(CFLAGS) -Iviewer $(SDL_CFLAGS) -c -o $@ $<

view: lightsim-view

lightsim-view: $(OBJ) $(VIEWOBJ)
ifeq ($(strip $(SDL_LIBS)),)
	@echo "lightsim-view needs SDL2. Install it with:  brew install sdl2"
	@exit 1
else
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(VIEWOBJ) $(SDL_LIBS) $(LDLIBS)
endif

# ui.c, font.c and status.c carry no SDL dependency, so the headless suite can
# exercise the toolbar rules, the glyphs and the message log without a window.
HEADLESS_VIEW := $(BUILD)/viewer_ui.o $(BUILD)/viewer_font.o $(BUILD)/viewer_inspect.o \
                 $(BUILD)/viewer_status.o

$(BUILD)/run_tests: $(OBJ) $(TESTOBJ) $(HEADLESS_VIEW)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(TESTOBJ) $(HEADLESS_VIEW) $(LDLIBS)

test: check-purity $(BUILD)/run_tests
	$(BUILD)/run_tests

debug:
	$(MAKE) VARIANT=debug all

test-asan:
	$(MAKE) VARIANT=asan test

test-ubsan:
	$(MAKE) VARIANT=ubsan test

# The units invariant from the design, enforced by the build rather than by
# discipline: photometry (Km, V(lambda), lux, candela) may appear ONLY in the
# colorimetry/units layer. If it reaches the transport core, a lumen can land in
# an accumulator that is later scaled by a radiometric BSDF.
PURE_SRC := $(filter-out src/units.c src/color.c src/cie_data.c src/spectrum.c,$(SRC))
# Word-anchored so "flux" does not match "lux".
PURE_PAT := \b683\b|ybar|ls_photometric|\blux\b|\bcandela\b|\blumens?\b
check-purity:
	@if [ -n "$(PURE_SRC)" ] && grep -nE '$(PURE_PAT)' $(PURE_SRC) 2>/dev/null; then \
	    echo "FAIL: photometry leaked out of the units layer (above)"; exit 1; \
	fi
	@echo "check-purity: transport core is free of photometric constants"

-include $(OBJ:.o=.d) $(TESTOBJ:.o=.d) $(APPOBJ:.o=.d) $(VIEWOBJ:.o=.d)

clean:
	rm -rf build lightsim lightsim-view run_tests
