//
// Created by Timo Ziegler on 31.08.24.
//

#ifndef TIMING_H
#include "timing.h"
#endif

#define MAXQUBITS   16

int main(void) {
    state_t state;
    FILE* fptr;

    char* file = "timePaulis/timePaulis.dat";
    fptr = fopen(file, "w");
    fprintf(fptr, "Pauli-X\n%-10s%-16s%-16s%-16s%-16s%-16s\n", "qubits", "old", "blas", "omp", "omp_blas", "gcd");

    for (qubit_t qubits = 10; qubits <= MAXQUBITS; ++qubits) {
        fprintf(fptr, "%-10d", qubits);
        timeQubitGate(&state, applyX, qubits, fptr);
        timeQubitGate(&state, applyX_blas, qubits, fptr);
        timeQubitGate(&state, applyX_omp, qubits, fptr);
        timeQubitGate(&state, applyX_blas_omp, qubits, fptr);
        timeQubitGate(&state, applyX_gcd, qubits, fptr);

        fprintf(fptr, "\n");
    }

    fprintf(fptr, "\nPauli-Y\n%-10s%-16s%-16s%-16s%-16s%-16s\n", "qubits", "old", "blas", "omp", "omp_blas", "gcd");

    for (qubit_t qubits = 10; qubits <= MAXQUBITS; ++qubits) {
        fprintf(fptr, "%-10d", qubits);
        timeQubitGate(&state, applyY, qubits, fptr);
        timeQubitGate(&state, applyY_blas, qubits, fptr);
        timeQubitGate(&state, applyY_omp, qubits, fptr);
        timeQubitGate(&state, applyY_blas_omp, qubits, fptr);
        timeQubitGate(&state, applyY_gcd, qubits, fptr);

        fprintf(fptr, "\n");
    }

    fprintf(fptr, "\nPauli-Z\n%-10s%-16s%-16s%-16s%-16s%-16s\n", "qubits", "old", "blas", "omp", "omp_blas", "gcd");

    for (qubit_t qubits = 10; qubits <= MAXQUBITS; ++qubits) {
        fprintf(fptr, "%-10d", qubits);
        timeQubitGate(&state, applyZ, qubits, fptr);
        timeQubitGate(&state, applyZblas, qubits, fptr);
        timeQubitGate(&state, applyZomp, qubits, fptr);
        timeQubitGate(&state, applyZomp_blas, qubits, fptr);
        timeQubitGate(&state, applyZgcd, qubits, fptr);

        fprintf(fptr, "\n");
    }

    fclose(fptr);

    return 0;
}