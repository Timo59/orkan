#ifndef Q_TEST_H
#include "test.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

int main(void) {
    const unsigned nqubits = 2;
    cplx_t **vecs = test_gen_states_pure(nqubits);

    // Print the Pauli-Z, -X, -Y basis states
    const unsigned dim = 1 << nqubits;

    printf("Computational basis states:\n");
    printf("---------------------------\n");
    for (unsigned i = 0; i < dim; i++) {
        vprint(vecs[i], dim);
    }
    printf("\n");

    printf("Hadamard basis states:\n");
    printf("----------------------\n");
    for (unsigned i = 0; i < dim; i++) {
        vprint(vecs[dim + i], dim);
    }
    printf("\n");

    printf("Circular basis states:\n");
    printf("---------------------\n");
    for (unsigned i = 0; i < dim; i++) {
        vprint(vecs[2 * dim + i], dim);
    }

    test_rm_states_pure(nqubits, vecs);
}
