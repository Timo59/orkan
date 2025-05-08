//
// Created by Timo Ziegler on 26.08.24.
//


#ifndef TIMING_H
#include "timing.h"
#endif

#define MAXQUBITS               24

double coeffsB[MAXQUBITS] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                            1.};
double coeffsR[MAXQUBITS] = {-3.40697, -9.56289, 3.20127, 6.16643, 1.66549, 6.99230, 0.46997, 4.51850, 2.80129,
                             -1.25172, -2.65258, 9.66951, 1.52371, 9.40358, 2.13824, 5.84664, 0.66015, 2.17288,
                             -7.35707, -4.24229, -3.75625, 7.95174, -8.86732, -7.60499};

int main() {
    state_t state;
    FILE* fptr;

    char* Bfile = "timeObsB.dat";
    char* Randfile = "timeObsRand.dat";

    fptr = fopen(Bfile, "w");
    fprintf(fptr, "OLD GATES \n");
    fprintf(fptr, "%-10s%-16s%-16s%-16s%-16s\n", "qubits", "old", "blas", "omp", "blas_omp");

    for (qubit_t qubits = 10; qubits <= MAXQUBITS; ++qubits) {
        fprintf(fptr, "%-10d", qubits);
        pauliObs_t obsB;
        pauli_t* pauliB = compsB(qubits);
        init_obs(qubits, &obsB, pauliB, coeffsB);
        timeObsPauli(&state, applyObsPauli, &obsB, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_blas, &obsB, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_omp, &obsB, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_blas_omp, &obsB, qubits, fptr);
        fprintf(fptr, "\n");
        free(pauliB);
    }
    fprintf(fptr, "\n");
    fclose(fptr);

    fptr = fopen(Randfile, "w");
    fprintf(fptr, "OLD GATES \n");
    fprintf(fptr, "%-10s%-16s%-16s%-16s%-16s\n", "qubits", "old", "blas", "omp", "omp_blas");

    for (qubit_t qubits = 10; qubits <= MAXQUBITS; ++qubits) {
        fprintf(fptr, "%-10d", qubits);
        pauliObs_t obsRand;
        pauli_t* pauliRand = compsRand(qubits);
        init_obs(qubits, &obsRand, pauliRand, coeffsR);
        timeObsPauli(&state, applyObsPauli, &obsRand, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_blas, &obsRand, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_omp, &obsRand, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_blas_omp, &obsRand, qubits, fptr);
        fprintf(fptr, "\n");
        free(pauliRand);
    }
    fprintf(fptr, "\n");
    fclose(fptr);

#define BLAS_GATES
    fptr = fopen(Bfile, "a");
    fprintf(fptr, "BLAS GATES \n");
    fprintf(fptr, "%-10s%-16s%-16s%-16s%-16s\n", "qubits", "old", "blas", "omp", "blas_omp");

    for (qubit_t qubits = 10; qubits <= MAXQUBITS; ++qubits) {
        fprintf(fptr, "%-10d", qubits);
        pauliObs_t obsB;
        pauli_t* pauliB = compsB(qubits);
        init_obs(qubits, &obsB, pauliB, coeffsB);
        timeObsPauli(&state, applyObsPauli, &obsB, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_blas, &obsB, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_omp, &obsB, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_blas_omp, &obsB, qubits, fptr);
        fprintf(fptr, "\n");
        free(pauliB);
    }
    fprintf(fptr, "\n");
    fclose(fptr);

    fptr = fopen(Randfile, "a");
    fprintf(fptr, "BLAS GATES \n");
    fprintf(fptr, "%-10s%-16s%-16s%-16s%-16s\n", "qubits", "old", "blas", "omp", "omp_blas");

    for (qubit_t qubits = 10; qubits <= MAXQUBITS; ++qubits) {
        fprintf(fptr, "%-10d", qubits);
        pauliObs_t obsRand;
        pauli_t* pauliRand = compsRand(qubits);
        init_obs(qubits, &obsRand, pauliRand, coeffsR);
        timeObsPauli(&state, applyObsPauli, &obsRand, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_blas, &obsRand, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_omp, &obsRand, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_blas_omp, &obsRand, qubits, fptr);
        fprintf(fptr, "\n");
        free(pauliRand);
    }
    fprintf(fptr, "\n");
    fclose(fptr);

#undef BLAS_GATES
#define BLAS_OMP_GATES
    fptr = fopen(Bfile, "a");
    fprintf(fptr, "BLAS OMP GATES \n");
    fprintf(fptr, "%-10s%-16s%-16s%-16s%-16s\n", "qubits", "old", "blas", "omp", "blas_omp");

    for (qubit_t qubits = 10; qubits <= MAXQUBITS; ++qubits) {
        fprintf(fptr, "%-10d", qubits);
        pauliObs_t obsB;
        pauli_t* pauliB = compsB(qubits);
        init_obs(qubits, &obsB, pauliB, coeffsB);
        timeObsPauli(&state, applyObsPauli, &obsB, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_blas, &obsB, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_omp, &obsB, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_blas_omp, &obsB, qubits, fptr);
        fprintf(fptr, "\n");
        free(pauliB);
    }
    fprintf(fptr, "\n");
    fclose(fptr);

    fptr = fopen(Randfile, "a");
    fprintf(fptr, "BLAS OMP GATES \n");
    fprintf(fptr, "%-10s%-16s%-16s%-16s%-16s\n", "qubits", "old", "blas", "omp", "omp_blas");

    for (qubit_t qubits = 10; qubits <= MAXQUBITS; ++qubits) {
        fprintf(fptr, "%-10d", qubits);
        pauliObs_t obsRand;
        pauli_t* pauliRand = compsRand(qubits);
        init_obs(qubits, &obsRand, pauliRand, coeffsR);
        timeObsPauli(&state, applyObsPauli, &obsRand, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_blas, &obsRand, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_omp, &obsRand, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_blas_omp, &obsRand, qubits, fptr);
        fprintf(fptr, "\n");
        free(pauliRand);
    }
    fprintf(fptr, "\n");
    fclose(fptr);
#undef BLAS_OMP_GATES
#define OMP_GATES
    fptr = fopen(Bfile, "a");
    fprintf(fptr, "OMP GATES \n");
    fprintf(fptr, "%-10s%-16s%-16s%-16s%-16s\n", "qubits", "old", "blas", "omp", "blas_omp");

    for (qubit_t qubits = 10; qubits <= MAXQUBITS; ++qubits) {
        fprintf(fptr, "%-10d", qubits);
        pauliObs_t obsB;
        pauli_t* pauliB = compsB(qubits);
        init_obs(qubits, &obsB, pauliB, coeffsB);
        timeObsPauli(&state, applyObsPauli, &obsB, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_blas, &obsB, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_omp, &obsB, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_blas_omp, &obsB, qubits, fptr);
        fprintf(fptr, "\n");
        free(pauliB);
    }
    fprintf(fptr, "\n");
    fclose(fptr);

    fptr = fopen(Randfile, "a");
    fprintf(fptr, "OMP GATES \n");
    fprintf(fptr, "%-10s%-16s%-16s%-16s%-16s\n", "qubits", "old", "blas", "omp", "omp_blas");

    for (qubit_t qubits = 10; qubits <= MAXQUBITS; ++qubits) {
        fprintf(fptr, "%-10d", qubits);
        pauliObs_t obsRand;
        pauli_t* pauliRand = compsRand(qubits);
        init_obs(qubits, &obsRand, pauliRand, coeffsR);
        timeObsPauli(&state, applyObsPauli, &obsRand, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_blas, &obsRand, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_omp, &obsRand, qubits, fptr);
        timeObsPauli(&state, applyObsPauli_blas_omp, &obsRand, qubits, fptr);
        fprintf(fptr, "\n");
        free(pauliRand);
    }
    fprintf(fptr, "\n");
    fclose(fptr);
#undef OMP_GATES

    return 0;
}