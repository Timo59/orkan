//
// Created by Timo Ziegler on 26.08.24.
//

#include "timing.h"

#define MAXQUBITS               24

double angles[MAXQUBITS] = {7.94578, -1.44730, 7.05692, 8.45498, -2.21193, 4.61068, 2.33483, -2.06368, 5.52334,
                            -8.51154, -5.35689, -0.27744, 2.92276, -8.99122, -0.75470, 7.95795, 9.44236, -7.86635,
                            9.26042, 0.42754, 3.28808, 8.28486, -0.85614, -6.73783};

int main() {
    state_t state;
    FILE *fptr;

    fptr = fopen("timeEvoString.dat", "w");
    fprintf(fptr, "OLD GATES\n");
    fprintf(fptr, "%-10s%-16s%-16s%-16s\n", "qubits", "old", "blas", "omp");

    for (qubit_t qubits = 10; qubits <= MAXQUBITS; ++qubits) {
        fprintf(fptr, "%-10d", qubits);
        pauli_t* comps = compsRand(qubits);
        timeEvoPauli(&state, evolvePauliStr, comps, angles, qubits, fptr);
        timeEvoPauli(&state, evolvePauliStr_omp, comps, angles, qubits, fptr);
        timeEvoPauli(&state, evolvePauliStr_blas, comps, angles, qubits, fptr);
        fprintf(fptr, "\n");
        free(comps);
    }
    fprintf(fptr, "\n");
    fclose(fptr);

#define BLAS_GATES
    fptr = fopen("timeEvoString.dat", "a");
    fprintf(fptr, "BLAS GATES\n");
    fprintf(fptr, "%-10s%-16s%-16s%-16s\n", "qubits", "old", "blas", "omp");

    for (qubit_t qubits = 10; qubits <= MAXQUBITS; ++qubits) {
        fprintf(fptr, "%-10d", qubits);
        pauli_t* comps = compsRand(qubits);
        timeEvoPauli(&state, evolvePauliStr, comps, angles, qubits, fptr);
        timeEvoPauli(&state, evolvePauliStr_omp, comps, angles, qubits, fptr);
        timeEvoPauli(&state, evolvePauliStr_blas, comps, angles, qubits, fptr);
        fprintf(fptr, "\n");
        free(comps);
    }
    fprintf(fptr, "\n");
    fclose(fptr);
#undef BLAS_GATES
#define BLAS_OMP_GATES
    fptr = fopen("timeEvoString.dat", "a");
    fprintf(fptr, "BLAS OMP GATES\n");
    fprintf(fptr, "%-10s%-16s%-16s%-16s\n", "qubits", "old", "blas", "omp");

    for (qubit_t qubits = 10; qubits <= MAXQUBITS; ++qubits) {
        fprintf(fptr, "%-10d", qubits);
        pauli_t* comps = compsRand(qubits);
        timeEvoPauli(&state, evolvePauliStr, comps, angles, qubits, fptr);
        timeEvoPauli(&state, evolvePauliStr_omp, comps, angles, qubits, fptr);
        timeEvoPauli(&state, evolvePauliStr_blas, comps, angles, qubits, fptr);
        fprintf(fptr, "\n");
        free(comps);
    }
    fprintf(fptr, "\n");
    fclose(fptr);
#undef BLAS_OMP_GATES
#define OMP_GATES
    fptr = fopen("timeEvoString.dat", "a");
    fprintf(fptr, "OMP GATES\n");
    fprintf(fptr, "%-10s%-16s%-16s%-16s\n", "qubits", "old", "blas", "omp");

    for (qubit_t qubits = 10; qubits <= MAXQUBITS; ++qubits) {
        fprintf(fptr, "%-10d", qubits);
        pauli_t* comps = compsRand(qubits);
        timeEvoPauli(&state, evolvePauliStr, comps, angles, qubits, fptr);
        timeEvoPauli(&state, evolvePauliStr_omp, comps, angles, qubits, fptr);
        timeEvoPauli(&state, evolvePauliStr_blas, comps, angles, qubits, fptr);
        fprintf(fptr, "\n");
        free(comps);
    }
    fprintf(fptr, "\n");
    fclose(fptr);
#undef OMP_GATES
    return 0;
}