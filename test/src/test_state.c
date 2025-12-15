#ifndef Q_TEST_H
#include "test.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

int main(void) {
    cplx_t **cb = test_cb_pure(2);
    cplx_t **xb = test_xb_pure(2);
    cplx_t **yb = test_yb_pure(2);

    // Print the Pauli-Z, -X, -Y basis states
    const unsigned dim = 1 << 2;

    printf("Computational basis states:\n");
    printf("---------------------------\n");
    for (unsigned i = 0; i < dim; i++) {
        vprint(cb[i]);
    }

    printf("Hadamard basis states:\n");
    printf("----------------------\n");
    for (unsigned i = 0; i < dim; i++) {
        vprint(xb[i]);
    }

    printf("Circular basis states:\n");
    printf("---------------------\n");
    for (unsigned i = 0; i < dim; i++) {
        vprint(yb[i]);
    }
}
