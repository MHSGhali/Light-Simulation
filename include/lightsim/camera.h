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

#endif /* LIGHTSIM_CAMERA_H */
