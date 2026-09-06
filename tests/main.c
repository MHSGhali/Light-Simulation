#include "test.h"

int ls_test_failures = 0;
int ls_test_count = 0;

#include "tests.h"

int main(void) {
    printf("lightsim validation suite\n=========================\n");
    test_spectral();
    test_transport();
    test_bsdf();
    test_furnace();
    printf("\n-------------------------\n%d checks, %d failure%s\n",
           ls_test_count, ls_test_failures, ls_test_failures == 1 ? "" : "s");
    return ls_test_failures == 0 ? 0 : 1;
}
