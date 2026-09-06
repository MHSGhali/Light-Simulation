/* thread.h — pthreads work distribution.
 *
 * Work items are handed out by an atomic counter, and callers seed their RNG
 * from the ITEM index rather than a thread id. That makes results bit-identical
 * regardless of thread count or scheduling order, which is what lets a reported
 * number be defensible and lets the tests compare runs.
 */
#ifndef LIGHTSIM_THREAD_H
#define LIGHTSIM_THREAD_H

#include "core.h"

typedef void (*LsWorkFn)(int item, void *user);

/* Run fn(0..n-1) across `nthreads` workers. nthreads <= 0 means "one per core".
 * Returns the number of workers actually used. */
int ls_parallel_for(int n, int nthreads, LsWorkFn fn, void *user);

/* Cores available, for defaulting and for reporting. */
int ls_hardware_threads(void);

#endif /* LIGHTSIM_THREAD_H */
