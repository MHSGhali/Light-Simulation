#include "lightsim/thread.h"
#include <pthread.h>
#include <stdatomic.h>
#include <stdlib.h>
#include <unistd.h>

int ls_hardware_threads(void) {
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return n > 0 ? (int)n : 1;
}

typedef struct {
    atomic_int next;
    int        n;
    LsWorkFn   fn;
    void      *user;
} Pool;

static void *worker(void *arg) {
    Pool *p = arg;
    for (;;) {
        int i = atomic_fetch_add(&p->next, 1);
        if (i >= p->n) break;
        p->fn(i, p->user);
    }
    return NULL;
}

int ls_parallel_for(int n, int nthreads, LsWorkFn fn, void *user) {
    if (n <= 0) return 0;
    if (nthreads <= 0) nthreads = ls_hardware_threads();
    if (nthreads > n) nthreads = n;

    Pool p;
    atomic_init(&p.next, 0);
    p.n = n; p.fn = fn; p.user = user;

    if (nthreads <= 1) { worker(&p); return 1; }

    pthread_t *th = malloc((size_t)nthreads * sizeof *th);
    if (!th) { worker(&p); return 1; }

    int started = 0;
    for (int i = 0; i < nthreads; ++i)
        if (pthread_create(&th[i], NULL, worker, &p) == 0) started++;
    if (started == 0) { worker(&p); free(th); return 1; }
    /* If some threads failed to start, this one helps drain the queue. */
    if (started < nthreads) worker(&p);
    for (int i = 0; i < started; ++i) pthread_join(th[i], NULL);
    free(th);
    return started;
}
