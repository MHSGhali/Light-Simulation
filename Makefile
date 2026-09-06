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

CFLAGS  := $(CSTD) $(WARN) $(OPT) $(SAN) -Iinclude
LDLIBS  := -lm -lpthread

BUILD   := build/$(VARIANT)
SRC     := $(wildcard src/*.c)
OBJ     := $(patsubst src/%.c,$(BUILD)/%.o,$(SRC))
TESTSRC := $(wildcard tests/*.c)
TESTOBJ := $(patsubst tests/%.c,$(BUILD)/tests_%.o,$(TESTSRC))
APPSRC  := $(wildcard apps/*.c)
APPOBJ  := $(patsubst apps/%.c,$(BUILD)/apps_%.o,$(APPSRC))

.PHONY: all test debug clean test-asan test-ubsan check-purity

all: lightsim

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

$(BUILD)/run_tests: $(OBJ) $(TESTOBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(TESTOBJ) $(LDLIBS)

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

clean:
	rm -rf build lightsim run_tests
