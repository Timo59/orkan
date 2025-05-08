//
// Created by Timo Ziegler on 17.09.24.
//

#include "timing.h"

#define CIRCDEPTH   4
double par[CIRCDEPTH] = {-5.24058, -2.96042, -9.73666, -0.06649};

int main(int argc, char* argv[]) {

    if (argc != 2) {
        printf("Usage: %s 'number of qubits'\n", argv[0]);
        return 0;
    }

    state_t state;
    unsigned int qubits = atoi(argv[1]);
    srand(mach_absolute_time());
    printf("Measure time of fundamental functionalities on %d qubits...\n", qubits);
    FILE* fptr;
    char* file = "timing.dat";

    /*
     * Timing of gate operations
     */
    printf("Gate operations...\n");
    fptr = fopen(file, "w");
    fprintf(fptr, "TIME GATESET\n");
    fprintf(fptr, "%-10s%-16s%-16s%-16s\n", "", "Pauli-X", "Pauli-Y", "Pauli-Z");

    fprintf(fptr, "%-10s", "Old");
    timeQubitGate(&state, applyX, qubits, fptr);
    timeQubitGate(&state, applyY, qubits, fptr);
    timeQubitGate(&state, applyZ, qubits, fptr);
    fprintf(fptr, "\n");

    fprintf(fptr, "%-10s", "BLAS");
    timeQubitGate(&state, applyX_blas, qubits, fptr);
    timeQubitGate(&state, applyY_blas, qubits, fptr);
    timeQubitGate(&state, applyZ_blas, qubits, fptr);
    fprintf(fptr, "\n");

    fprintf(fptr, "%-10s", "BLAS OMP");
    timeQubitGate(&state, applyX_blas_omp, qubits, fptr);
    timeQubitGate(&state, applyY_blas_omp, qubits, fptr);
    timeQubitGate(&state, applyZ_blas_omp, qubits, fptr);
    fprintf(fptr, "\n");

    fprintf(fptr, "%-10s", "OMP");
    timeQubitGate(&state, applyX_omp, qubits, fptr);
    timeQubitGate(&state, applyY_omp, qubits, fptr);
    timeQubitGate(&state, applyZ_omp, qubits, fptr);
    fprintf(fptr, "\n");

    fprintf(fptr, "%-10s", "GCD");
    timeQubitGate(&state, applyX_gcd, qubits, fptr);
    timeQubitGate(&state, applyY_gcd, qubits, fptr);
    timeQubitGate(&state, applyZ_gcd, qubits, fptr);
    fprintf(fptr, "\n");

    fclose(fptr);
    printf("DONE\n");

    /*
     * Timing of observable application
     */
    printf("Observable applications...\n");
    pauliObs_t observable;
    pauli_t* Bcomps = compsB(qubits);
    double coeffsB[qubits];
    for (unsigned int i = 0; i < qubits; ++i) {
        coeffsB[i] = 1.;
    }
    init_obs(qubits, &observable, Bcomps, coeffsB);
    fptr = fopen(file, "a");
    fprintf(fptr, "\nTIME B-OBSERVABLE\n");
    fprintf(fptr, "%-16s%-16s%-16s%-16s%-16s\n", "", "Old", "BLAS", "BLAS OMP", "OMP");
    fprintf(fptr, "%-16s", "Old gates");
    timeObsPauli(&state, applyObsPauli_old, &observable, qubits, fptr, applyPauliStr);
    timeObsPauli(&state, applyObsPauli_blas, &observable, qubits, fptr, applyPauliStr);
    timeObsPauli(&state, applyObsPauli_omp, &observable, qubits, fptr, applyPauliStr);
    timeObsPauli(&state, applyObsPauli_blas_omp, &observable, qubits, fptr, applyPauliStr);
    fprintf(fptr, "\n");

    fprintf(fptr, "%-16s", "BLAS gates");
    timeObsPauli(&state, applyObsPauli_old, &observable, qubits, fptr, applyPauliStr_blas);
    timeObsPauli(&state, applyObsPauli_blas, &observable, qubits, fptr, applyPauliStr_blas);
    timeObsPauli(&state, applyObsPauli_omp, &observable, qubits, fptr, applyPauliStr_blas);
    timeObsPauli(&state, applyObsPauli_blas_omp, &observable, qubits, fptr, applyPauliStr_blas);
    fprintf(fptr, "\n");

    fprintf(fptr, "%-16s", "BLAS OMP gates");
    timeObsPauli(&state, applyObsPauli_old, &observable, qubits, fptr, applyPauliStr_blas_omp);
    timeObsPauli(&state, applyObsPauli_blas, &observable, qubits, fptr, applyPauliStr_blas_omp);
    timeObsPauli(&state, applyObsPauli_omp, &observable, qubits, fptr, applyPauliStr_blas_omp);
    timeObsPauli(&state, applyObsPauli_blas_omp, &observable, qubits, fptr, applyPauliStr_blas_omp);
    fprintf(fptr, "\n");

    fprintf(fptr, "%-16s", "OMP gates");
    timeObsPauli(&state, applyObsPauli_old, &observable, qubits, fptr, applyPauliStr_omp);
    timeObsPauli(&state, applyObsPauli_blas, &observable, qubits, fptr, applyPauliStr_omp);
    timeObsPauli(&state, applyObsPauli_omp, &observable, qubits, fptr, applyPauliStr_omp);
    timeObsPauli(&state, applyObsPauli_blas_omp, &observable, qubits, fptr, applyPauliStr_omp);
    fprintf(fptr, "\n");

    fclose(fptr);

    pauli_t* compsR = compsRand(qubits);
    double coeffsR[qubits];
    for (unsigned int i = 0; i < qubits; ++i) {
        coeffsR[i] = (double) rand() / (double) (RAND_MAX / 10);
    }
    init_obs(qubits, &observable, compsR, coeffsR);
    fptr = fopen(file, "a");
    fprintf(fptr, "\nTIME RANDOM PAULI OBSERVABLE\n");
    fprintf(fptr, "%-16s%-16s%-16s%-16s%-16s\n", "", "Old", "BLAS", "BLAS OMP", "OMP");
    fprintf(fptr, "%-16s", "Old gates");
    timeObsPauli(&state, applyObsPauli_old, &observable, qubits, fptr, applyPauliStr);
    timeObsPauli(&state, applyObsPauli_blas, &observable, qubits, fptr, applyPauliStr);
    timeObsPauli(&state, applyObsPauli_omp, &observable, qubits, fptr, applyPauliStr);
    timeObsPauli(&state, applyObsPauli_blas_omp, &observable, qubits, fptr, applyPauliStr);
    fprintf(fptr, "\n");

    fprintf(fptr, "%-16s", "BLAS gates");
    timeObsPauli(&state, applyObsPauli_old, &observable, qubits, fptr, applyPauliStr_blas);
    timeObsPauli(&state, applyObsPauli_blas, &observable, qubits, fptr, applyPauliStr_blas);
    timeObsPauli(&state, applyObsPauli_omp, &observable, qubits, fptr, applyPauliStr_blas);
    timeObsPauli(&state, applyObsPauli_blas_omp, &observable, qubits, fptr, applyPauliStr_blas);
    fprintf(fptr, "\n");

    fprintf(fptr, "%-16s", "BLAS OMP gates");
    timeObsPauli(&state, applyObsPauli_old, &observable, qubits, fptr, applyPauliStr_blas_omp);
    timeObsPauli(&state, applyObsPauli_blas, &observable, qubits, fptr, applyPauliStr_blas_omp);
    timeObsPauli(&state, applyObsPauli_omp, &observable, qubits, fptr, applyPauliStr_blas_omp);
    timeObsPauli(&state, applyObsPauli_blas_omp, &observable, qubits, fptr, applyPauliStr_blas_omp);
    fprintf(fptr, "\n");

    fprintf(fptr, "%-16s", "OMP gates");
    timeObsPauli(&state, applyObsPauli_old, &observable, qubits, fptr, applyPauliStr_omp);
    timeObsPauli(&state, applyObsPauli_blas, &observable, qubits, fptr, applyPauliStr_omp);
    timeObsPauli(&state, applyObsPauli_omp, &observable, qubits, fptr, applyPauliStr_omp);
    timeObsPauli(&state, applyObsPauli_blas_omp, &observable, qubits, fptr, applyPauliStr_omp);
    fprintf(fptr, "\n");

    fclose(fptr);
    printf("DONE\n");

    /*
     * Timing of evolution with Pauli strings
     */
    printf("Pauli string evolution...\n");
    init_obs(qubits, &observable, Bcomps, coeffsB);
    fptr = fopen(file, "a");
    fprintf(fptr, "\nTIME EVOLUTION B OBSERVABLE\n");
    fprintf(fptr, "%-16s%-16s%-16s%-16s\n", "", "Old", "BLAS", "OMP");
    fprintf(fptr, "%-16s", "Old gates");
    timeEvoPauli(&state, evolvePauliStr_old, Bcomps, coeffsR, qubits, fptr, applyPauliStr);
    timeEvoPauli(&state, evolvePauliStr_omp, Bcomps, coeffsR, qubits, fptr, applyPauliStr);
    timeEvoPauli(&state, evolvePauliStr_blas, Bcomps, coeffsR, qubits, fptr, applyPauliStr);
    fprintf(fptr, "\n");

    fprintf(fptr, "%-16s", "BLAS gates");
    timeEvoPauli(&state, evolvePauliStr_old, Bcomps, coeffsR, qubits, fptr, applyPauliStr_blas);
    timeEvoPauli(&state, evolvePauliStr_omp, Bcomps, coeffsR, qubits, fptr, applyPauliStr_blas);
    timeEvoPauli(&state, evolvePauliStr_blas, Bcomps, coeffsR, qubits, fptr, applyPauliStr_blas);
    fprintf(fptr, "\n");

    fprintf(fptr, "%-16s", "BLAS OMP gates");
    timeEvoPauli(&state, evolvePauliStr_old, Bcomps, coeffsR, qubits, fptr, applyPauliStr_blas_omp);
    timeEvoPauli(&state, evolvePauliStr_omp, Bcomps, coeffsR, qubits, fptr, applyPauliStr_blas_omp);
    timeEvoPauli(&state, evolvePauliStr_blas, Bcomps, coeffsR, qubits, fptr, applyPauliStr_blas_omp);
    fprintf(fptr, "\n");

    fprintf(fptr, "%-16s", "OMP gates");
    timeEvoPauli(&state, evolvePauliStr_old, Bcomps, coeffsR, qubits, fptr, applyPauliStr_omp);
    timeEvoPauli(&state, evolvePauliStr_omp, Bcomps, coeffsR, qubits, fptr, applyPauliStr_omp);
    timeEvoPauli(&state, evolvePauliStr_blas, Bcomps, coeffsR, qubits, fptr, applyPauliStr_omp);
    fprintf(fptr, "\n");

    fclose(fptr);

    init_obs(qubits, &observable, compsR, coeffsR);
    fptr = fopen(file, "a");
    fprintf(fptr, "\nTIME EVOLUTION RANDOM PAULI OBSERVABLE\n");
    fprintf(fptr, "%-16s%-16s%-16s%-16s\n", "", "Old", "BLAS", "OMP");
    fprintf(fptr, "%-16s", "Old gates");
    timeEvoPauli(&state, evolvePauliStr_old, compsR, coeffsR, qubits, fptr, applyPauliStr);
    timeEvoPauli(&state, evolvePauliStr_omp, compsR, coeffsR, qubits, fptr, applyPauliStr);
    timeEvoPauli(&state, evolvePauliStr_blas, compsR, coeffsR, qubits, fptr, applyPauliStr);
    fprintf(fptr, "\n");

    fprintf(fptr, "%-16s", "BLAS gates");
    timeEvoPauli(&state, evolvePauliStr_old, compsR, coeffsR, qubits, fptr, applyPauliStr_blas);
    timeEvoPauli(&state, evolvePauliStr_omp, compsR, coeffsR, qubits, fptr, applyPauliStr_blas);
    timeEvoPauli(&state, evolvePauliStr_blas, compsR, coeffsR, qubits, fptr, applyPauliStr_blas);
    fprintf(fptr, "\n");

    fprintf(fptr, "%-16s", "BLAS OMP gates");
    timeEvoPauli(&state, evolvePauliStr_old, compsR, coeffsR, qubits, fptr, applyPauliStr_blas_omp);
    timeEvoPauli(&state, evolvePauliStr_omp, compsR, coeffsR, qubits, fptr, applyPauliStr_blas_omp);
    timeEvoPauli(&state, evolvePauliStr_blas, compsR, coeffsR, qubits, fptr, applyPauliStr_blas_omp);
    fprintf(fptr, "\n");

    fprintf(fptr, "%-16s", "OMP gates");
    timeEvoPauli(&state, evolvePauliStr_old, compsR, coeffsR, qubits, fptr, applyPauliStr_omp);
    timeEvoPauli(&state, evolvePauliStr_omp, compsR, coeffsR, qubits, fptr, applyPauliStr_omp);
    timeEvoPauli(&state, evolvePauliStr_blas, compsR, coeffsR, qubits, fptr, applyPauliStr_omp);
    fprintf(fptr, "\n");

    fclose(fptr);
    printf("DONE\n");

    /*
     * Timing of observable's expectation value
     */
    printf("Expectation value observable...\n");
    fptr = fopen(file, "a");
    fprintf(fptr, "\nTIME MEAN RANDOM PAULI OBSERVABLE\n");
    fprintf(fptr, "%-16s%-16s%-16s%-16s%-16s\n", "", "Old", "BLAS", "BLAS OMP", "OMP");
    fprintf(fptr, "%-16s", "Old gates");
    timeMeanObs(&state, expValObsPauli_old, &observable, qubits, fptr, applyPauliStr);
    timeMeanObs(&state, expValObsPauli_blas, &observable, qubits, fptr, applyPauliStr);
    timeMeanObs(&state, expValObsPauli_blas_omp, &observable, qubits, fptr, applyPauliStr);
    timeMeanObs(&state, expValObsPauli_omp, &observable, qubits, fptr, applyPauliStr);
    fprintf(fptr, "\n");
#define BLAS_GATES
    fprintf(fptr, "%-16s", "BLAS gates");
    timeMeanObs(&state, expValObsPauli_old, &observable, qubits, fptr, applyPauliStr_blas);
    timeMeanObs(&state, expValObsPauli_blas, &observable, qubits, fptr, applyPauliStr_blas);
    timeMeanObs(&state, expValObsPauli_blas_omp, &observable, qubits, fptr, applyPauliStr_blas);
    timeMeanObs(&state, expValObsPauli_omp, &observable, qubits, fptr, applyPauliStr_blas);
    fprintf(fptr, "\n");
#undef BLAS_GATES
#define BLAS_OMP_GATES
    fprintf(fptr, "%-16s", "BLAS OMP gates");
    timeMeanObs(&state, expValObsPauli_old, &observable, qubits, fptr, applyPauliStr_blas_omp);
    timeMeanObs(&state, expValObsPauli_blas, &observable, qubits, fptr, applyPauliStr_blas_omp);
    timeMeanObs(&state, expValObsPauli_blas_omp, &observable, qubits, fptr, applyPauliStr_blas_omp);
    timeMeanObs(&state, expValObsPauli_omp, &observable, qubits, fptr, applyPauliStr_blas_omp);
    fprintf(fptr, "\n");
#undef BLAS_OMP_GATES
#define OMP_GATES
    fprintf(fptr, "%-16s", "OMP gates");
    timeMeanObs(&state, expValObsPauli_old, &observable, qubits, fptr, applyPauliStr_omp);
    timeMeanObs(&state, expValObsPauli_blas, &observable, qubits, fptr, applyPauliStr_omp);
    timeMeanObs(&state, expValObsPauli_blas_omp, &observable, qubits, fptr, applyPauliStr_omp);
    timeMeanObs(&state, expValObsPauli_omp, &observable, qubits, fptr, applyPauliStr_omp);
    fprintf(fptr, "\n");
#undef OMP_GATES

    fclose(fptr);
    printf("DONE\n");

/*    *//*
     * Timing of observable's gradient wrt to a PQC
     *//*
    printf("Gradient PQC...\n");
    obs_t* evoOps[CIRCDEPTH];                                       // Declare array of evolutional operators
    pauliObs_t* evoPauli[CIRCDEPTH];                                // Declare array of evolutional Pauli obs
    for (depth_t i = 0; i < CIRCDEPTH; ++i) {
        evoPauli[i] = malloc(sizeof(*(evoPauli[i])));           // Initialize one evolutional Pauli obs
        evoPauli[i]->comps = compsRand(qubits);
        evoPauli[i]->coeffs = coeffsR;
        evoPauli[i]->qubits = qubits;
        evoPauli[i]->length = (complength_t) qubits;

        evoOps[i] = malloc(sizeof(*(evoOps[i])));               // Initialize one evolutional operator
        evoOps[i]->type = PAULI;
        evoOps[i]->pauliObs = evoPauli[i];
    }
    obs_t obs;
    obs.type = PAULI;
    obs.pauliObs = &observable;
    fptr = fopen(file, "a");
    fprintf(fptr, "\nTIME GRADIENT RANDOM PAULI OBSERVABLE\n");
    fprintf(fptr, "\nOLD CIRC\n");
    fprintf(fptr, "%-16s%-16s%-16s%-16s%-16s\n", "", "Old", "BLAS", "BLAS OMP", "OMP");
    fprintf(fptr, "%-16s", "Old gates");
    timeGradientPQC(&state, gradPQC, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    fprintf(fptr, "\n");
#define BLAS_GATES
    fprintf(fptr, "%-16s", "BLAS gates");
    timeGradientPQC(&state, gradPQC, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    fprintf(fptr, "\n");
#undef BLAS_GATES
#define BLAS_OMP_GATES
    fprintf(fptr, "%-16s", "BLAS OMP gates");
    timeGradientPQC(&state, gradPQC, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    fprintf(fptr, "\n");
#undef BLAS_OMP_GATES
#define OMP_GATES
    fprintf(fptr, "%-16s", "OMP gates");
    timeGradientPQC(&state, gradPQC, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    fprintf(fptr, "\n");
#undef OMP_GATES

#define BLAS_CIRC
    fprintf(fptr, "\nBLAS CIRC\n");
    fprintf(fptr, "%-16s%-16s%-16s%-16s%-16s\n", "", "Old", "BLAS", "BLAS OMP", "OMP");
    fprintf(fptr, "%-16s", "Old gates");
    timeGradientPQC(&state, gradPQC, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    fprintf(fptr, "\n");
#define BLAS_GATES
    fprintf(fptr, "%-16s", "BLAS gates");
    timeGradientPQC(&state, gradPQC, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    fprintf(fptr, "\n");
#undef BLAS_GATES
#define BLAS_OMP_GATES
    fprintf(fptr, "%-16s", "BLAS OMP gates");
    timeGradientPQC(&state, gradPQC, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    fprintf(fptr, "\n");
#undef BLAS_OMP_GATES
#define OMP_GATES
    fprintf(fptr, "%-16s", "OMP gates");
    timeGradientPQC(&state, gradPQC, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    fprintf(fptr, "\n");
#undef OMP_GATES
#undef BLAS_CIRC

#define BLAS_OMP_CIRC
    fprintf(fptr, "\nBLAS OMP CIRC\n");
    fprintf(fptr, "%-16s%-16s%-16s%-16s%-16s\n", "", "Old", "BLAS", "BLAS OMP", "OMP");
    fprintf(fptr, "%-16s", "Old gates");
    timeGradientPQC(&state, gradPQC, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    fprintf(fptr, "\n");
#define BLAS_GATES
    fprintf(fptr, "%-16s", "BLAS gates");
    timeGradientPQC(&state, gradPQC, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    fprintf(fptr, "\n");
#undef BLAS_GATES
#define BLAS_OMP_GATES
    fprintf(fptr, "%-16s", "BLAS OMP gates");
    timeGradientPQC(&state, gradPQC, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    fprintf(fptr, "\n");
#undef BLAS_OMP_GATES
#define OMP_GATES
    fprintf(fptr, "%-16s", "OMP gates");
    timeGradientPQC(&state, gradPQC, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    fprintf(fptr, "\n");
#undef OMP_GATES
#undef BLAS_OMP_CIRC

#define OMP_CIRC
    fprintf(fptr, "\nOMP CIRC\n");
    fprintf(fptr, "%-16s%-16s%-16s%-16s%-16s\n", "", "Old", "BLAS", "BLAS OMP", "OMP");
    fprintf(fptr, "%-16s", "Old gates");
    timeGradientPQC(&state, gradPQC, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    fprintf(fptr, "\n");
#define BLAS_GATES
    fprintf(fptr, "%-16s", "BLAS gates");
    timeGradientPQC(&state, gradPQC, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    fprintf(fptr, "\n");
#undef BLAS_GATES
#define BLAS_OMP_GATES
    fprintf(fptr, "%-16s", "BLAS OMP gates");
    timeGradientPQC(&state, gradPQC, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    fprintf(fptr, "\n");
#undef BLAS_OMP_GATES
#define OMP_GATES
    fprintf(fptr, "%-16s", "OMP gates");
    timeGradientPQC(&state, gradPQC, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    timeGradientPQC(&state, gradPQC_blas_omp, par, &obs, evoOps, CIRCDEPTH, qubits, fptr);
    fprintf(fptr, "\n");
#undef OMP_GATES
#undef OMP_CIRC

    fclose(fptr);
    printf("DONE\n");
    for (unsigned int i = 0; i < CIRCDEPTH; ++i) {
        free(evoOps[i]->pauliObs->comps);
        free(evoOps[i]->pauliObs);
        free(evoOps[i]);
    }*/
    free(Bcomps);
    free(compsR);
    return 1;
}