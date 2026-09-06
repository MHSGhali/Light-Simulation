# lightsim — spectral radiometric lighting simulator
# GNU Make 3.81 compatible (macOS ships 3.81): no .ONESHELL, no $(file ...).

CC      := cc
CSTD    := -std=c11 -D_DEFAULT_SOURCE
WARN    := -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wdouble-promotion \
           -Wstrict-prototypes -Wmissing-prototypes
# -ffp-contract=off keeps FMA contraction from shifting results at the 1e-9
# level, which would make the deterministic tests flaky across -O levels.
# -fno-fast-math is a policy, not a preference: fast-math licenses reassociation
# of the spectral sums and disables the NaN/Inf semantics the estimators rely on
# to fail loudly rather than quietly.
OPT     := -O3 -ffp-contract=off -fno-fast-math
CFLAGS  := $(CSTD) $(WARN) $(OPT) -Iinclude
LDLIBS  := -lm -lpthread

BUILD   := build
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

run_tests: $(OBJ) $(TESTOBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(TESTOBJ) $(LDLIBS)

test: check-purity run_tests
	./run_tests

# The units invariant from the design, enforced by the build rather than by
# discipline: photometry (Km, V(lambda), lux, candela) may appear ONLY in the
# colorimetry/units layer. If it shows up in the transport core, a lumen can
# reach an accumulator that is later scaled by a radiometric BSDF.
PURE_SRC := $(filter-out src/units.c src/color.c src/cie_data.c src/spectrum.c,$(SRC))
check-purity:
	@if [ -n "$(PURE_SRC)" ] && \
	    grep -n -e '\b683\b' -e 'ybar' -e 'photometric' -e 'lux' -e 'candela' \
	         $(PURE_SRC) 2>/dev/null; then \
	    echo "FAIL: photometry leaked out of the units layer (above)"; exit 1; \
	fi
	@echo "check-purity: transport core is free of photometric constants"



debug: CFLAGS := $(CSTD) $(WARN) -O0 -g3 -Iinclude
debug: clean lightsim

# Sanitizers rebuild from scratch at -O1 so the traps are meaningful.
test-asan:
	$(MAKE) clean
	$(MAKE) OPT="-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer" test

test-ubsan:
	$(MAKE) clean
	$(MAKE) OPT="-O1 -g -fsanitize=undefined -fno-sanitize-recover=all" test

clean:
	rm -rf $(BUILD) lightsim run_tests
