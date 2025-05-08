//
// Created by Timo Ziegler on 08.10.24.
//

#include "timing.h"

#define MAXQUBITS   24
#define MAXDEPTH    8

double par[MAXDEPTH] = {9.74095, -0.93873, 3.10498, 8.29608, -6.84950, -9.80926, -3.81814, 6.05164};
double coeffs[MAXQUBITS] = {4.23369, -3.31494, 5.11566, 7.09852, 8.81981, 7.05529, 3.21002, -6.99860, -8.08379,
                              5.50778, -5.52321, 9.49860, 1.21603, 3.95303, -9.37972, -9.87905, -0.04753, 9.63584,
                              -8.52856, -5.91887, 3.37913, -2.52934, 4.84267, 5.21956};

int main() {
    state_t state;
    FILE* fptr;
    char* file = "walltimeGradPQC.dat";
    fptr = fopen(file, "w");

    for (depth_t circdepth = 1; circdepth <= MAXDEPTH; circdepth *= 2) {
        fprintf(fptr, "Circdepth: %d\n%-10s%-16s%-16s%-16s%-16s%-16s\n",
                circdepth, "qubits", "old", "blas", "omp", "omp_blas", "gcd");
        for (qubit_t qubits = 10; qubits <= MAXQUBITS; ++qubits) {
            fprintf(fptr, "%-10d", qubits);

            pauliObs_t obsPauli;                                            // Declare Pauli observable
            pauli_t* compsObs = compsRand(qubits);                          // Random comps for Pauli observable
            init_obs(qubits, &obsPauli, compsObs, coeffs);      // Initialize Pauli observable

            obs_t obs;                                                      // Declare observable
            obs.type = PAULI;                                               // Initialize observable
            obs.pauliObs = &obsPauli;

            obs_t* evoOps[circdepth];                                       // Declare array of evolutional operators
            pauliObs_t* evoPauli[circdepth];                                // Declare array of evolutional Pauli obs
            for (depth_t i = 0; i < circdepth; ++i) {
                evoPauli[i] = malloc(sizeof(*(evoPauli[i])));           // Initialize one evolutional Pauli obs
                evoPauli[i]->comps = compsRand(qubits);
                evoPauli[i]->coeffs = coeffs;
                evoPauli[i]->qubits = qubits;
                evoPauli[i]->length = (complength_t) qubits;

                evoOps[i] = malloc(sizeof(*(evoOps[i])));               // Initialize one evolutional operator
                evoOps[i]->type = PAULI;
                evoOps[i]->pauliObs = evoPauli[i];
            }

            timeGradientPQC(&state, gradPQC, par, &obs, evoOps, circdepth, qubits, fptr);
            timeGradientPQC(&state, gradPQCblas, par, &obs, evoOps, circdepth, qubits, fptr);
            timeGradientPQC(&state, gradPQComp, par, &obs, evoOps, circdepth, qubits, fptr);
            timeGradientPQC(&state, gradPQCblas_omp, par, &obs, evoOps, circdepth, qubits, fptr);

            for (depth_t i = 0; i < circdepth; ++i) {
                free(evoPauli[i]->comps);
                free(evoPauli[i]);
                free(evoOps[i]);
            }
            free(compsObs);
            fprintf(fptr, "\n");
        }
        fprintf(fptr, "\n");
    }

    fclose(fptr);
    return 0;
}
