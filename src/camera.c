#include "lightsim/camera.h"

Camera ls_camera_look_at(vec3 eye, vec3 target, vec3 up_hint,
                         ls_real fov_y_deg, int width, int height) {
    Camera c;
    c.eye = eye;
    c.fwd = v3norm(v3sub(target, eye));
    c.right = v3norm(v3cross(c.fwd, up_hint));
    /* Rebuild up from the orthogonalised right so the basis stays orthonormal
     * even when up_hint is not perpendicular to the view direction. */
    c.up = v3cross(c.right, c.fwd);
    c.fov_y = fov_y_deg * LS_PI / 180.0;
    c.width = width;
    c.height = height;
    return c;
}

Ray ls_camera_ray(const Camera *c, int x, int y, ls_real jx, ls_real jy) {
    ls_real aspect = (ls_real)c->width / (ls_real)c->height;
    ls_real half_h = tan(0.5 * c->fov_y);
    ls_real half_w = half_h * aspect;
    /* Pixel centre in NDC, y flipped so row 0 is the top of the image. */
    ls_real sx = (2.0 * (((ls_real)x + jx) / (ls_real)c->width) - 1.0) * half_w;
    ls_real sy = (1.0 - 2.0 * (((ls_real)y + jy) / (ls_real)c->height)) * half_h;
    Ray r;
    r.o = c->eye;
    r.d = v3norm(v3add(c->fwd, v3add(v3scale(c->right, sx), v3scale(c->up, sy))));
    r.tmin = 0.0;
    r.tmax = HUGE_VAL;
    return r;
}
