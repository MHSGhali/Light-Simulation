/* import.h — bringing geometry in from elsewhere.
 *
 * WHAT COMES BACK, AND WHAT DOES NOT
 *   Blender writes OBJ natively, so that is the road in. An OBJ carries
 *   vertices, faces and a material name; it does not carry spectra, lights,
 *   cameras or units. So this module reads geometry and a base colour, and
 *   everything physical stays authored in the .scene file beside it.
 *
 *   Deliberately refused rather than guessed:
 *     Ke  emissive materials. A mesh emitter is found by BSDF sampling but is
 *         NOT in the light list, so ls_estimate_irradiance would miss it
 *         entirely and ls_estimate_irradiance_full would discard its first-hit
 *         emission as "already counted by direct sampling" -- which for a mesh
 *         it never was. The grid would come out systematically low with a
 *         plausible-looking number and no test firing. Use `light rect`.
 *     Ks/Ns  speculars. LS_BSDF_CONDUCTOR takes measured spectral eta/kappa;
 *         there is no honest map from a Phong exponent to it. Write
 *         `material steel metal al 0.15` by hand instead.
 *     vn  vertex normals. geom.h requires the raw geometric normal, and a
 *         shading normal that disagrees with it breaks energy conservation at
 *         grazing angles. This simulator's output is a number, not a picture.
 *
 *   Kd becomes a spectrum through ls_spectrum_from_rgb_reflectance, and the
 *   authored triple is kept on the Material so the scene writer can emit it
 *   back as `material <name> rgb r g b`.
 *
 * STL
 *   Also read, in both its ASCII and binary forms. It is a simpler format and
 *   a poorer one: no materials at all (not one -- zero), no names, no groups,
 *   so an STL arrives as a single mesh with a default grey to be replaced by
 *   hand.
 *
 *   AND NO UNITS. This is the trap. STL is unitless by convention and CAD
 *   tools almost universally mean millimetres, while Blender writes metres.
 *   Nothing in the file distinguishes them, so nothing can detect it: a part
 *   exported from a CAD tool arrives 1000x too large and swallows the room.
 *
 *   Rather than guess, `scale` is explicit and defaults to 1 (metres, matching
 *   every other length in this engine), and the importer prints the bounding
 *   box of what it read. A units mistake then reads as "40 x 25 x 12 m" on the
 *   way in rather than as a mysteriously black measurement afterwards.
 */
#ifndef LIGHTSIM_IMPORT_H
#define LIGHTSIM_IMPORT_H

#include "scenefile.h"

/* Directory containing `path`, with the trailing separator removed. "" when
 * `path` has no directory part. */
void ls_dir_of(const char *path, char *out, size_t n);

/* `rel` resolved against `base_dir`, unless it is already absolute. */
void ls_resolve_path(const char *base_dir, const char *rel, char *out, size_t n);

/* One `usemtl` group of an OBJ (or the whole file when `group` is NULL) as a
 * built, BVH'd mesh with one reference. NULL on failure, with why in `err`.
 * Vertices are taken as authored -- Blender's exporter bakes the object's world
 * matrix, so the caller places the result at the identity. */
Mesh *ls_obj_load(const char *path, const char *group, ls_real scale,
                  char err[256]);

/* An STL, ASCII or binary -- the form is detected from the content, since a
 * binary STL may also begin with the bytes "solid". One mesh, no material. */
Mesh *ls_stl_load(const char *path, ls_real scale, char err[256]);

/* Add every usemtl group of `path` to `d` as its own Mesh and Prim, with
 * materials from the sidecar `mtllib` where there is one. One prim per group,
 * because a Prim carries exactly one material.
 *
 * Returns the number of prims added, or -1 with why in `d->err`. */
int ls_scene_import_obj(SceneDesc *d, const char *path, ls_real scale);
int ls_scene_import_obj_at(SceneDesc *d, const char *path, ls_real scale,
                           const vec3 *at);

/* Import any supported geometry file, picking the reader by extension.
 * `scale` multiplies every vertex; 1.0 means the file is already in metres.
 *
 * `at` places the result: the footprint is centred on (at.x, at.y) and the
 * BOTTOM of the geometry is set to at.z, so `at = (0,0,0)` stands the part on
 * the floor at the origin. NULL leaves the coordinates exactly as authored.
 *
 * The distinction matters more than it looks. A CAD tool lays parts out on a
 * build plate, so an STL's coordinates are typically a few hundred millimetres
 * from the origin in x and y -- import one as-authored into a room centred on
 * the origin and it lands outside the walls, traced perfectly and invisible. */
int ls_scene_import_at(SceneDesc *d, const char *path, ls_real scale,
                       const vec3 *at);

/* As authored. */
int ls_scene_import(SceneDesc *d, const char *path, ls_real scale);

/* The scale a file of `file_extent` metres most likely wants, to join a scene
 * spanning `scene_extent` metres: 1, or 0.001 when the geometry is so far out
 * of proportion that it can only be a millimetre file read as metres.
 *
 * A guess, and only defensible where it is announced and undoable -- the
 * viewer, which has no command line to take a scale on. The CLI asks instead. */
ls_real ls_import_units_hint(ls_real file_extent, ls_real scene_extent);

#endif /* LIGHTSIM_IMPORT_H */
