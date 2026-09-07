/* export.h — writes a self-contained Blender Python script that rebuilds the
 * simulated scene AND carries the simulation result into it.
 *
 * Two things go across:
 *   1. The physical setup — geometry with matching albedo, and real Blender
 *      lights whose power is the same radiant flux in watts the simulator used,
 *      tinted by the source's own spectrum. Cycles treats area-light power as
 *      radiant flux in W, so phi_e maps across with no fudge factor.
 *   2. The measured field — the illuminance/irradiance grid as a mesh with one
 *      face per measurement point, false-coloured through the same viridis ramp
 *      the PPM writer uses, on an emission shader so it reads in the viewport.
 *
 * Convention: 1 simulator world unit = 1 metre, which is Blender's internal
 * unit, so coordinates cross unchanged.
 */
#ifndef LIGHTSIM_EXPORT_H
#define LIGHTSIM_EXPORT_H

#include "scenefile.h"

bool ls_export_blender(const SceneDesc *d,
                       const ls_real *grid, int nu, int nv,
                       ls_real lo, ls_real hi, const char *unit_name,
                       const char *quantity_name, const char *path);

#endif /* LIGHTSIM_EXPORT_H */
