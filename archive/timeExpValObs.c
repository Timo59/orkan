//
// Created by Timo Ziegler on 08.10.24.
//

#include "timing.h"
#include "profile_exact.h"

#define QUBITS               24

double coeffs[QUBITS] = {-9.82646, -6.20157, 6.46052, -5.89780, 9.10301, -4.63317, -1.72492, -4.27666, 6.15689,
                            7.29567, 7.06829, -1.78119, -1.44909, -7.92318, -2.14942, -4.07057, -3.78105, 9.93711,
                            5.58031, -7.02385, 0.02740, 1.68612, 5.73233, 7.03682};

int main() {
    state_t state;
    FILE* fptr;

    char* file = "timeExpValObs.dat";

    fptr = fopen(file, "w");
    fprintf(fptr, "%d qubits", QUBITS);
    fprintf(fptr, "OLD EXACT");
    fprintf(fptr, "%-16s%-16s%-16s%-16s%-16s\n", "", "old circ", "blas circ", "blas-omp circ", "omp circ");

    pauliObs_t obs;
    pauli_t* comps = compsRand(QUBITS);
    init_obs(QUBITS, &obs, comps, coeffs);

    fprintf(fptr, "%-16s", "old gates");
#define OLD_CIRC
    timeMeanObs(&state, expValObsPauli, &obs, QUBITS, fptr);
#undef OLD_CIRC
#define BLAS_CIRC
    timeMeanObs(&state, expValObsPauli, &obs, QUBITS, fptr);
#undef BLAS_CIRC
#define BLAS_OMP_CIRC
    timeMeanObs(&state, expValObsPauli, &obs, QUBITS, fptr);
#undef BLAS_OMP_CIRC
#define OMP_CIRC
    timeMeanObs(&state, expValObsPauli, &obs, QUBITS, fptr);
#undef OMP_CIRC
    fprintf(fptr, "\n")

#define

    timeMeanObs(&state, expValObsPauli_blas, &obs, qubits, fptr);
    timeMeanObs(&state, expValObsPauli_omp, &obs, qubits, fptr);
    timeMeanObs(&state, expValObsPauli_blas_omp, &obs, qubits, fptr);
    fprintf(fptr, "\n");
    fclose(fptr);
#undef OLD_GATES

    return 0;
}