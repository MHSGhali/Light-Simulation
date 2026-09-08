#!/usr/bin/env python3
"""Regenerate the README's images from the simulator itself.

    make && python3 tools/make_docs.py

Every picture in the README is a real run: the renders come from `lightsim
render`, the field maps are coloured from the CSV `lightsim grid` writes, and
the ramp is the viridis table viewer/draw.c uses, so the README and the tool
agree. Nothing here is drawn by hand or touched up.

The one screenshot, docs/viewer.png, is not generated -- it is a capture of the
running viewer.

Needs Pillow.
"""
import csv, json, math, os, subprocess
from PIL import Image, ImageDraw, ImageFont

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "docs")
TMP = os.path.join(REPO, "out", "docgen")
SCENE = "scenes/workcell.scene"
os.makedirs(OUT, exist_ok=True)
os.makedirs(TMP, exist_ok=True)

# viewer/draw.c draw_viridis(), verbatim, so the README matches the app.
V = [(0.267004, 0.004874, 0.329415), (0.282623, 0.140926, 0.457517),
     (0.253935, 0.265254, 0.529983), (0.206756, 0.371758, 0.553117),
     (0.163625, 0.471133, 0.558148), (0.127568, 0.566949, 0.550556),
     (0.134692, 0.658636, 0.517649), (0.266941, 0.748751, 0.440573),
     (0.477504, 0.821444, 0.318195), (0.741388, 0.873449, 0.149561),
     (0.993248, 0.906157, 0.143936)]


def viridis(t):
    t = 0.0 if t < 0 else (1.0 if t > 1 else t)
    s = t * 10.0
    i = min(int(s), 9)
    u = s - i
    return tuple(int(255 * (V[i][k] + u * (V[i + 1][k] - V[i][k])) + 0.5)
                 for k in range(3))


def run(args):
    subprocess.run(args, cwd=REPO, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def load_grid(prefix):
    """(values[j][i], nu, nv, unit) from the CSV `lightsim grid` writes."""
    with open(prefix + ".csv") as f:
        rows = [r for r in f if not r.startswith("#")]
    vals, unit = {}, None
    for r in csv.DictReader(rows):
        key = [k for k in r if k.startswith("value_")][0]
        unit = key.split("value_", 1)[1]
        vals[(int(r["j"]), int(r["i"]))] = float(r[key])
    nu = max(k[1] for k in vals) + 1
    nv = max(k[0] for k in vals) + 1
    return [[vals[(j, i)] for i in range(nu)] for j in range(nv)], nu, nv, unit


HEAD = 30      # caption strip, so a label never sits on top of the data


def field_image(vals, nu, nv, lo, hi, px, caption, unit, sub=None):
    """A false-coloured field map with a colour bar and a caption strip.

    `lo`/`hi` are passed in rather than taken per-image: a sequence of frames
    renormalised individually would animate the colour map instead of the
    physics, which is the one thing a measurement picture must not do."""
    cell = max(1, px // max(nu, nv))
    w, h = nu * cell, nv * cell
    img = Image.new("RGB", (w, HEAD + h + 34), (18, 24, 26))
    d = ImageDraw.Draw(img)
    f = ImageFont.load_default()

    d.text((4, 4), caption, font=f, fill=(240, 245, 244))
    if sub:
        d.text((4, 16), sub, font=f, fill=(84, 195, 180))

    span = (hi - lo) if hi > lo else 1.0
    for j in range(nv):
        for i in range(nu):
            c = viridis((vals[nv - 1 - j][i] - lo) / span)      # +y up
            d.rectangle([i * cell, HEAD + j * cell,
                         (i + 1) * cell - 1, HEAD + (j + 1) * cell - 1], fill=c)

    by = HEAD + h + 8
    for x in range(w):
        d.line([(x, by), (x, by + 9)], fill=viridis(x / max(1, w - 1)))
    d.rectangle([0, by, w - 1, by + 9], outline=(60, 74, 78))
    d.text((2, by + 13), f"{lo:.0f}", font=f, fill=(160, 175, 174))
    d.text((w - 52, by + 13), f"{hi:.0f} {unit}", font=f, fill=(160, 175, 174))
    return img


def scene_with(replace):
    """The scene with whole lines rewritten by `replace(line) -> line or None`."""
    out = []
    for line in open(os.path.join(REPO, SCENE)):
        line = line.rstrip("\n")
        out.append(replace(line) or line)
    return "\n".join(out) + "\n"


def hero():
    print("hero render")
    run(["./lightsim", "render", SCENE, "--spp", "900", "--depth", "6",
         "--out", f"{TMP}/hero.ppm"])
    Image.open(f"{TMP}/hero.ppm").convert("RGB").save(f"{OUT}/render.png")


def orbit(n=48):
    print(f"orbit, {n} frames")
    tgt, rad, zoff = (0.0, 0.0, 0.06), 0.72, 0.24
    frames = []
    for k in range(n):
        az = -math.pi / 2 + 2 * math.pi * k / n
        ex, ey = tgt[0] + rad * math.cos(az), tgt[1] + rad * math.sin(az)

        def rep(line, ex=ex, ey=ey):
            if line.startswith("camera"):
                return (f"camera {ex:.5f} {ey:.5f} {tgt[2] + zoff:.5f}  "
                        f"{tgt[0]} {tgt[1]} {tgt[2]}  42 384 288")
            return None

        open(f"{TMP}/orb.scene", "w").write(scene_with(rep))
        run(["./lightsim", "render", f"{TMP}/orb.scene", "--spp", "110",
             "--depth", "5", "--out", f"{TMP}/orb.ppm"])
        frames.append(Image.open(f"{TMP}/orb.ppm").convert("RGB"))
    frames[0].save(f"{OUT}/orbit.gif", save_all=True, append_images=frames[1:],
                   duration=70, loop=0, optimize=True)


def mounting_height():
    # The panels must stay inside the enclosure: the ceiling quad is at z = 0.5,
    # and a panel above it lights the outside of the ceiling and nothing else.
    print("mounting-height sweep")
    heights = [0.15 + 0.02 * k for k in range(17)]
    runs = []
    for z in heights + heights[::-1][1:-1]:          # ping-pong
        def rep(line, z=z):
            if line.startswith("light rect"):
                t = line.split()
                t[4] = f"{z:.3f}"
                return " ".join(t)
            return None

        open(f"{TMP}/h.scene", "w").write(scene_with(rep))
        # U0 is min/mean -- one worst point -- so it is the noisiest number
        # here and needs the samples. A grid is 30 ms; there is no reason to
        # skimp, and at low counts the sweep shows noise as if it were physics.
        run(["./lightsim", "grid", f"{TMP}/h.scene", "--spp", "1600",
             "--indirect", "640", "--depth", "5", "--out", f"{TMP}/h"])
        vals, nu, nv, unit = load_grid(f"{TMP}/h")
        st = json.load(open(f"{TMP}/h.json"))["stats"]
        runs.append((z, vals, nu, nv, unit, st))
        print(f"   {z*100:.0f} cm  mean {st['mean']:7.1f}  U0 {st['u0']:.3f}")

    lo = min(min(min(r) for r in v[1]) for v in runs)
    hi = max(max(max(r) for r in v[1]) for v in runs)
    imgs = [field_image(v, nu, nv, lo, hi, 300, f"panels at {z*100:.0f} cm",
                        unit, f"mean {st['mean']:.0f} lx   U0 {st['u0']:.2f}")
            for (z, v, nu, nv, unit, st) in runs]
    imgs[0].save(f"{OUT}/mounting-height.gif", save_all=True,
                 append_images=imgs[1:], duration=150, loop=0, optimize=True)


def direct_vs_full():
    print("direct vs full transport")
    run(["./lightsim", "grid", SCENE, "--spp", "512", "--indirect", "0",
         "--depth", "0", "--out", f"{TMP}/d"])
    run(["./lightsim", "grid", SCENE, "--spp", "512", "--indirect", "256",
         "--depth", "5", "--out", f"{TMP}/f"])
    dv, nu, nv, unit = load_grid(f"{TMP}/d")
    fv, _, _, _ = load_grid(f"{TMP}/f")
    ds = json.load(open(f"{TMP}/d.json"))["stats"]
    fs = json.load(open(f"{TMP}/f.json"))["stats"]
    lo = min(min(min(r) for r in dv), min(min(r) for r in fv))
    hi = max(max(max(r) for r in dv), max(max(r) for r in fv))
    a = field_image(dv, nu, nv, lo, hi, 300, "direct only", unit,
                    f"mean {ds['mean']:.0f}  U0 {ds['u0']:.2f}")
    b = field_image(fv, nu, nv, lo, hi, 300, "full transport", unit,
                    f"mean {fs['mean']:.0f}  U0 {fs['u0']:.2f}")
    side = Image.new("RGB", (a.width + b.width + 12, a.height), (18, 24, 26))
    side.paste(a, (0, 0))
    side.paste(b, (a.width + 12, 0))
    side.save(f"{OUT}/direct-vs-full.png")
    print(f"   interreflection adds {100*(fs['mean']/ds['mean']-1):.0f}% "
          f"to the mean; {ds['occluded']} points get no direct light at all")


if __name__ == "__main__":
    hero()
    orbit()
    mounting_height()
    direct_vs_full()
    print("wrote docs/")
