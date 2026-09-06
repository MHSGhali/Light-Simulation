/* analysis.h — descriptive statistics over a measurement set.
 *
 * These are the numbers a lighting study actually reports: the spread of the
 * distribution and how uniform it is. Kept unit-agnostic -- feed it lux or
 * W/m^2 and the ratios come out the same. */
#ifndef LIGHTSIM_ANALYSIS_H
#define LIGHTSIM_ANALYSIS_H

#include "core.h"

typedef struct {
    ls_real min, max, mean, stddev;   /* stddev uses the N-1 denominator */
    ls_real michelson;                /* (max-min)/(max+min)             */
    ls_real u0;                       /* min/mean, the CIE uniformity U0 */
    ls_real ud;                       /* min/max, the diversity ratio    */
    int     count;
} LsStats;

LsStats ls_stats(const ls_real *v, int n);

#endif /* LIGHTSIM_ANALYSIS_H */
