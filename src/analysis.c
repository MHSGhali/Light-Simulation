#include "lightsim/analysis.h"
#include <math.h>

LsStats ls_stats(const ls_real *v, int n) {
    LsStats s;
    s.min = s.max = s.mean = s.stddev = s.michelson = s.u0 = s.ud = 0.0;
    s.count = n;
    if (n <= 0) return s;

    s.min = s.max = v[0];
    ls_real sum = 0.0;
    for (int i = 0; i < n; ++i) {
        if (v[i] < s.min) s.min = v[i];
        if (v[i] > s.max) s.max = v[i];
        sum += v[i];
    }
    s.mean = sum / (ls_real)n;

    if (n > 1) {
        ls_real acc = 0.0;
        for (int i = 0; i < n; ++i) acc += ls_sqr(v[i] - s.mean);
        s.stddev = sqrt(acc / (ls_real)(n - 1));
    }
    ls_real denom = s.max + s.min;
    s.michelson = (denom > 0.0) ? (s.max - s.min) / denom : 0.0;
    s.u0 = (s.mean > 0.0) ? s.min / s.mean : 0.0;
    s.ud = (s.max  > 0.0) ? s.min / s.max  : 0.0;
    return s;
}
