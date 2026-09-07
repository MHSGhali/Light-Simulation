#include "lightsim/camera.h"
#include <string.h>

Camera ls_camera_look_at(vec3 eye, vec3 target, vec3 up_hint,
                         ls_real fov_y_deg, int width, int height) {
    /* Zeroed first: every field must be set, and `ortho` in particular is a
     * _Bool whose garbage value would make a perspective camera behave as an
     * orthographic one at random. Adding a field to Camera must not be able to
     * leave it uninitialised here. */
    Camera c;
    memset(&c, 0, sizeof c);
    c.eye = eye;
    c.fwd = v3norm(v3sub(target, eye));
    c.right = v3norm(v3cross(c.fwd, up_hint));
    /* Rebuild up from the orthogonalised right so the basis stays orthonormal
     * even when up_hint is not perpendicular to the view direction. */
    c.up = v3cross(c.right, c.fwd);
    c.fov_y = fov_y_deg * LS_PI / 180.0;
    c.width = width;
    c.height = height;
    c.ortho = false;
    c.ortho_height = 0.0;
    return c;
}

Camera ls_camera_ortho(vec3 eye, vec3 target, vec3 up_hint,
                       ls_real height_world, int width, int height) {
    Camera c = ls_camera_look_at(eye, target, up_hint, 45.0, width, height);
    c.ortho = true;
    c.ortho_height = height_world;
    return c;
}

ls_real ls_camera_world_per_pixel(const Camera *c, vec3 p) {
    if (c->ortho) return c->ortho_height / (ls_real)c->height;
    ls_real z = v3dot(v3sub(p, c->eye), c->fwd);
    if (z < 1e-6) z = 1e-6;
    return z * 2.0 * tan(0.5 * c->fov_y) / (ls_real)c->height;
}

Ray ls_camera_ray(const Camera *c, int x, int y, ls_real jx, ls_real jy) {
    return ls_camera_pick_ray(c, (ls_real)x + jx, (ls_real)y + jy);
}

Ray ls_camera_pick_ray(const Camera *c, ls_real sx, ls_real sy) {
    ls_real aspect = (ls_real)c->width / (ls_real)c->height;
    ls_real half_h = c->ortho ? 0.5 * c->ortho_height : tan(0.5 * c->fov_y);
    ls_real half_w = half_h * aspect;
    /* Pixel centre in NDC, y flipped so row 0 is the top of the image. */
    ls_real nx = (2.0 * (sx / (ls_real)c->width)  - 1.0) * half_w;
    ls_real ny = (1.0 - 2.0 * (sy / (ls_real)c->height)) * half_h;
    Ray r;
    if (c->ortho) {
        /* Parallel rays: the film position moves the ORIGIN, not the direction. */
        r.o = v3add(c->eye, v3add(v3scale(c->right, nx), v3scale(c->up, ny)));
        r.d = c->fwd;
    } else {
        r.o = c->eye;
        r.d = v3norm(v3add(c->fwd, v3add(v3scale(c->right, nx), v3scale(c->up, ny))));
    }
    r.tmin = 0.0;
    r.tmax = HUGE_VAL;
    return r;
}

bool ls_camera_project(const Camera *c, vec3 p, ls_real *sx, ls_real *sy) {
    vec3 v = v3sub(p, c->eye);
    ls_real z = v3dot(v, c->fwd);
    if (z <= 1e-9) return false;               /* at or behind the eye plane */
    ls_real aspect = (ls_real)c->width / (ls_real)c->height;
    ls_real half_h = c->ortho ? 0.5 * c->ortho_height : tan(0.5 * c->fov_y);
    ls_real half_w = half_h * aspect;
    /* Orthographic: no divide by depth, which is exactly what makes the view
     * measurable. */
    ls_real nx = v3dot(v, c->right) / (c->ortho ? 1.0 : z);
    ls_real ny = v3dot(v, c->up)    / (c->ortho ? 1.0 : z);
    *sx = (nx / half_w + 1.0) * 0.5 * (ls_real)c->width;
    *sy = (1.0 - ny / half_h) * 0.5 * (ls_real)c->height;
    return true;
}
