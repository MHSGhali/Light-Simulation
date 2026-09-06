# Light-Simulation

A physically-based **spectral** ray tracer in C that computes radiometric
quantities from first principles and reports any of them in either radiometric
or photometric units.

Radiance is carried as 95 bins across 360–830 nm. Photometric quantities come
from an exact `∫ Φ(λ)·V(λ) dλ` against the CIE 1931 observer, so switching units
re-projects one stored spectrum rather than applying a correction factor.

## Build

```sh
make            # library + CLI      (no dependencies)
make view       # interactive viewer (needs SDL2: brew install sdl2)
make test       # validation suite
```

The CLI, the library and the tests build with no SDL2 installed; the viewer is
the only thing that needs it.

## Run

**Characterise a source spectrum**

```sh
./lightsim source blackbody 2856 --units photometric
./lightsim source daylight 6504 --units radiometric
./lightsim source led 620 25
./lightsim source mono 555                # 683.0000 lm, exactly
```

**Measure a scene**

```sh
./lightsim grid scenes/workcell.scene --units photometric \
    --spp 256 --indirect 256 --depth 5 \
    --out out/wc --blender out/wc.py
```

Writes a CSV of every measurement point, a false-colour PPM, a JSON summary,
and optionally a Blender script. `--indirect 0 --depth 0` gives direct light
only, for comparison against a closed-form model.

**Render a view**

```sh
./lightsim render scenes/workcell.scene --spp 256 --depth 6 --out out/render.ppm
```

Writes a tone-mapped PPM plus a PFM holding raw radiance in physical units.

**Interactive viewer**

```sh
./lightsim-view scenes/workcell.scene            # field map
./lightsim-view scenes/workcell.scene --render   # progressive render
```

| key | action | | key | action |
|---|---|---|---|---|
| `1` | field map | | `R` | re-solve |
| `2` | progressive render | | `S` | save PPM |
| `U` | lux ↔ W/m² | | `B` | export Blender script |
| `T` | full ↔ direct-only | | `Esc` | quit |
| `Q` | draft ↔ fine | | | |

Hover the field map to probe a point; the cross-section follows the row under
the cursor. In render mode, drag to orbit and scroll to zoom — the film resets
and re-accumulates.

## Blender export

`--blender FILE.py`, or `B` in the viewer, writes a self-contained script.
Open it in Blender's Scripting tab and press Run. It builds:

- the scene geometry with matching albedo,
- real Blender lights carrying the **same radiant flux in watts** the simulation
  used — Cycles reads area-light power as radiant flux, so a Cycles render of
  the exported scene is directly comparable with the simulated field,
- the measurement grid as a mesh, one face per point, false-coloured through the
  same viridis ramp, with the raw value kept in a `Value` face attribute so the
  numbers survive the trip.

## Scene files

```
material floor lambert 0.20
material wall  lambert 0.80
material steel metal al 0.15          # spectral complex IOR

quad   floor  0 0 0    0 0 1   0.3 0 0   0 0.3 0
sphere steel  0.06 0.02 0.09   0.05

light rect -0.12 -0.12 0.45  0.05 0 0  0 -0.05 0  lm 200 daylight 5000
light spot 0 0 0.4  0 0 -1  30 20  W 5 blackbody 3000

grid   -0.3 -0.3 0.001   0.6 0 0   0 0.6 0   64 64
camera 0 -0.72 0.30      0 0 0.06   42 480 360
```

Lamp flux may be given in `lm` or `W`; lumens are converted to watts once, at
load, so no photometric value is ever stored internally.

Lights: `point`, `spot`, `rect`, `disk`, `sphere`, `sun`.
Spectra: `flat`, `blackbody <K>`, `daylight <CCT>`, `led <centre_nm> <fwhm_nm>`.

## Validation

`make test` runs 5 237 assertions, each pinned to a closed form or a published
CIE constant, clean under AddressSanitizer and UBSan at both 5 nm and 10 nm bins.

| check | result |
|---|---|
| 1 W at 555 nm | 683.000 lm exactly |
| ∫V(λ)dλ | 106.857 nm |
| blackbody 2856 K, LER vs total radiant power | 16.45 lm/W |
| blackbody LER peak | 95.43 lm/W at 6628 K |
| illuminant A chromaticity | (0.44753, 0.40743) |
| uniform sphere ≡ point source of equal flux | 0.015 % apart |
| furnace, ρ = 0.9 → `L = Le/(1−ρ)` | 10.008 vs 10 |
| NEE vs BSDF-sampling vs MIS, same scene | agree to 0.5 % |

The furnace test is the one that matters most: it fails loudly on
energy-conservation bugs, missing cosine factors, wrong PDFs, and bad Russian
roulette — the errors that are otherwise invisible in a picture that looks fine.

## Licence

MIT. See `LICENSE`.
