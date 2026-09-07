/* camera.h — pinhole camera and ray generation. */
#ifndef LIGHTSIM_CAMERA_H
#define LIGHTSIM_CAMERA_H

#include "geom.h"

typedef struct {
    vec3    eye;
    vec3    fwd, right, up;   /* orthonormal */
    ls_real fov_y;            /* radians, vertical */
    int     width, height;
} Camera;

Camera ls_camera_look_at(vec3 eye, vec3 target, vec3 up_hint,
                         ls_real fov_y_deg, int width, int height);

/* Ray through pixel (x,y) offset by (jx,jy) in [0,1) for antialiasing. */
Ray ls_camera_ray(const Camera *c, int x, int y, ls_real jx, ls_real jy);

/* Ray through a continuous pixel coordinate, for picking. Same mapping as
 * ls_camera_ray with the jitter folded into the fractional part. */
Ray ls_camera_pick_ray(const Camera *c, ls_real sx, ls_real sy);

/* World point -> continuous pixel coordinate: the inverse of ls_camera_pick_ray.
 * Returns false when the point is at or behind the eye plane, where no pixel
 * corresponds to it. Used to place overlay gizmos and to pick delta lights,
 * which have no geometry to intersect. */
bool ls_camera_project(const Camera *c, vec3 p, ls_real *sx, ls_real *sy);

#endif /* LIGHTSIM_CAMERA_H */
