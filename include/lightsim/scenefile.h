/* scenefile.h — a small line-oriented scene description.
 *
 *   # comment
 *   camera   <eye xyz> <target xyz> <fov_y_deg> <width> <height>
 *   grid     <origin xyz> <edge_u xyz> <edge_v xyz> <nu> <nv>
 *   material <name> lambert <albedo>
 *   material <name> metal   <al|cu|au> <alpha>
 *   material <name> emit    <radiance>
 *   plane    <mat> <centre xyz> <normal xyz>
 *   quad     <mat> <centre xyz> <normal xyz> <half_u xyz> <half_v xyz>
 *   sphere   <mat> <centre xyz> <radius>
 *   light point  <xyz> <W|lm> <value> <spd>
 *   light spot   <xyz> <dir xyz> <total_deg> <falloff_deg> <W|lm> <value> <spd>
 *   light rect   <centre xyz> <half_u xyz> <half_v xyz> <W|lm> <value> <spd>
 *   light disk   <centre xyz> <normal xyz> <radius> <W|lm> <value> <spd>
 *   light sphere <centre xyz> <radius> <W|lm> <value> <spd>
 *   light sun    <dir xyz> <irradiance W/m2> <spd>
 *
 *   spd := flat | blackbody <K> | daylight <K> | led <centre_nm> <fwhm_nm>
 *
 * A `rect`, `disk` or `sphere` light also inserts matching emissive geometry,
 * so the source is visible to the camera and to BSDF sampling, and is bound to
 * the light so MIS can weight the two strategies against each other.
 */
#ifndef LIGHTSIM_SCENEFILE_H
#define LIGHTSIM_SCENEFILE_H

#include "scene.h"
#include "camera.h"

typedef struct {
    Scene    scene;
    Camera   camera;
    bool     has_camera;
    /* The camera as authored. Camera itself stores an orthonormal basis, from
     * which the eye/target pair that produced it cannot be recovered uniquely,
     * so the writer keeps the original numbers. */
    vec3     cam_eye, cam_target;
    ls_real  cam_fov_deg;
    vec3     grid_o, grid_u, grid_v;
    int      grid_nu, grid_nv;
    bool     has_grid;

    Prim     *prims;   int nprims,  cap_prims;
    Material *mats;    int nmats,   cap_mats;
    char    (*names)[32];
    Light    *lights;  int nlights, cap_lights;
    char      err[256];
} SceneDesc;

bool ls_scene_load(SceneDesc *d, const char *path);
void ls_scene_desc_free(SceneDesc *d);

/* Re-emit the scene in the same format ls_scene_load reads, so an edited
 * session round-trips, stays diffable, and can be re-run from the CLI.
 *
 * Emissive geometry paired to an area light is deliberately NOT written: it is
 * regenerated on load by the same code that created it in the first place, so
 * writing it would duplicate every area light on the next load. */
bool ls_scene_save(const SceneDesc *d, const char *path);

#endif /* LIGHTSIM_SCENEFILE_H */
