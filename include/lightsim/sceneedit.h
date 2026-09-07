/* sceneedit.h — mutating a loaded scene.
 *
 * The parser builds a SceneDesc once; these operations let an editor change it
 * afterwards. Two invariants are maintained by every one of them, and are the
 * reason this is a module rather than a handful of inline array pokes:
 *
 *   1. `Light.index` always equals the light's position in the array, because
 *      that index is the contribution-matrix column and is what `Prim.light_id`
 *      refers to.
 *
 *   2. Every rect, disk and sphere light has exactly one paired Prim carrying an
 *      emissive Material, so the source is visible to the camera and findable by
 *      BSDF sampling. If that pairing drifts, the render and the field map
 *      disagree -- the render shows a light the NEE path no longer samples, or
 *      the reverse -- and neither picture looks obviously wrong.
 *      ls_scene_sync_light_geom() is the ONLY place the pairing is established.
 *
 * Every mutator calls ls_scene_rebuild() before returning, so `d->scene` never
 * points at a stale array after a realloc.
 */
#ifndef LIGHTSIM_SCENEEDIT_H
#define LIGHTSIM_SCENEEDIT_H

#include "scenefile.h"

/* Refresh d->scene from the backing arrays. Cheap; call after any mutation. */
void ls_scene_rebuild(SceneDesc *d);

/* Deep copy, for undo snapshots. `dst` is overwritten and takes ownership of
 * its own allocations. Returns false only on allocation failure. */
bool ls_scene_clone(const SceneDesc *src, SceneDesc *dst);

/* Create or refresh the emissive geometry paired to light `index`, and bring it
 * into line with the light's current kind, placement, extent and radiance.
 * Delta lights (point, spot, sun) have no geometry; any stale prim is removed. */
void ls_scene_sync_light_geom(SceneDesc *d, int index);

/* Append a light. It is finalized (which asserts its flux normalisation) and
 * its geometry synced. Returns its index, or -1. */
int  ls_scene_add_light(SceneDesc *d, Light l);

/* Remove a light and its paired geometry, compacting the arrays and remapping
 * every surviving Light.index and Prim.light_id. */
void ls_scene_remove_light(SceneDesc *d, int index);

/* Re-finalize a light after its fields were edited in place, and resync its
 * geometry. Call this rather than touching d->lights[i] and hoping. */
void ls_scene_update_light(SceneDesc *d, int index);

/* Duplicate a light, offset by `offset`. Returns the new index, or -1. */
int  ls_scene_duplicate_light(SceneDesc *d, int index, vec3 offset);

int  ls_scene_add_material(SceneDesc *d, Material m, const char *name);
int  ls_scene_add_prim(SceneDesc *d, Prim p);
void ls_scene_remove_prim(SceneDesc *d, int index);

/* Index of the prim paired to light `index`, or -1. */
int  ls_scene_light_prim(const SceneDesc *d, int index);

#endif /* LIGHTSIM_SCENEEDIT_H */
