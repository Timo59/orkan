//
// Created by Timo Ziegler on 13.09.24.
//

#ifndef TIMING_H
#define TIMING_H

#ifndef QLIB_H
#include "qLib.h"
#endif

#ifndef _MACH_MACH_TIME_H_
#include <mach/mach_time.h>
#endif

#ifndef _STDIO_H_
#include <stdio.h>
#endif

/*
 * =====================================================================================================================
 *                                              type definitions
 * =====================================================================================================================
 */
typedef void (*qubitGate)(state_t* state, qubit_t target);
typedef void (*obsPauli)(state_t* state, const pauliObs_t* obs, ptrStr applyStr);
typedef void (*evoPauli)(state_t* state, const pauli_t paulistr[], double angle, ptrStr applyStr);
typedef double (*meanObs)(const state_t* state, const pauliObs_t* obs, ptrStr applyStr);

/*
 * =====================================================================================================================
 *                                              helper functions
 * =====================================================================================================================
 */

uint64_t get_time() {
    static mach_timebase_info_data_t info = {0, 0};
    if (info.denom == 0) {
        mach_timebase_info(&info);
    }
    uint64_t now = mach_absolute_time();
    now *= info.numer;
    now /= info.denom;
    return now;
}

pauli_t* compsB(qubit_t qubits) {                                               // Creates the Pauli components for the
    pauli_t* out = calloc(qubits * qubits, sizeof(pauli_t ));       // B observable on qubits qubits
    for (qubit_t i = 0; i < qubits; ++i) {
        out[i * qubits + i] = X;
    }
    return out;
}

pauli_t* compsRand(qubit_t qubits) {                                            // Creates random Pauli components for
    srand(mach_absolute_time());                                                // for qubits qubits
    pauli_t* out = malloc(qubits * qubits * sizeof(pauli_t));
    for (int i = 0; i < qubits * qubits; ++i) {
        out[i] = rand() % 4;
    }
    return out;
}

void init_obs(qubit_t qubits, pauliObs_t* obs, pauli_t* comps, double coeffs[]) {
    obs->comps = comps;
    obs->coeffs = coeffs;
    complength_t length = qubits;
    obs->length = length;
    obs->qubits = qubits;
}

/*
 * =====================================================================================================================
 *                                              time measurements
 * =====================================================================================================================
 */

void timeQubitGate(state_t* state, qubitGate applyGate, qubit_t qubits, FILE* fptr) {
    stateInitPlus(state, qubits);                                                   // Initialize state to plus state
    uint64_t start = get_time();                                                    // Start timing
    for (qubit_t target = 0; target < qubits; ++target) {                           // Iterate all qubits
        applyGate(state, target);                                                   // Apply gate to the current qubit
    }
    uint64_t end = get_time();                                                      // End timing
    uint64_t time_elapsed = end - start;                                            // Difference between times
    time_elapsed /= (long long) qubits;                                             // Divide difference by number of
                                                                                    // executions
    fprintf(fptr, "%-16.9f", (double) time_elapsed / NSEC_PER_SEC);                 // Add time in seconds to file
    stateFreeVector(state);
}

void timeObsPauli(state_t* state, obsPauli applyPauliObs, pauliObs_t* obs, qubit_t qubits, FILE* fptr, ptrStr applyStr)
{
    stateInitPlus(state, qubits);                                                   // Initialize state to plus state
    uint64_t start = get_time();                                                    // Start timing
    applyPauliObs(state, obs, applyStr);                                            // Apply observable
    uint64_t end = get_time();                                                      // End timing
    uint64_t time_elapsed = end - start;                                            // Difference between times
    fprintf(fptr, "%-16.9f", (double) time_elapsed / NSEC_PER_SEC);                 // Add time in seconds to file
    stateFreeVector(state);
}

void timeEvoPauli(state_t* state,
                  evoPauli evolution,
                  pauli_t* paulistrs,
                  double angles[],
                  qubit_t qubits,
                  FILE* fptr,
                  ptrStr applyStr)
                  {
    stateInitPlus(state, qubits);                                                   // Initialize state to plus state
    uint64_t start = get_time();                                                    // Start timing
    for (qubit_t i = 0; i < qubits; ++i) {                                          // Iterate all Pauli strings
        evolution(state, paulistrs + i * qubits, angles[i], applyStr);// Apply string evolution for current
    }                                                                               // Pauli string
    uint64_t end = get_time();                                                      // End timing
    uint64_t time_elapsed = end - start;                                            // Difference between times
    time_elapsed /= (long long) qubits;                                             // Divide difference by number of
                                                                                    // executions
    fprintf(fptr, "%-16.9f", (double) time_elapsed / NSEC_PER_SEC);                 // Add time in seconds to file
    stateFreeVector(state);
}

void timeMeanObs(state_t* state, meanObs mean, const pauliObs_t* obs, qubit_t qubits, FILE* fptr, ptrStr applyStr) {
    stateInitPlus(state, qubits);                                                   // Initialize state to plus state
    uint64_t start = get_time();                                                    // Start timing
    mean(state, obs, applyStr);                                                     // Apply observable
    uint64_t end = get_time();                                                      // End timing
    uint64_t time_elapsed = end - start;                                            // Difference between times
    fprintf(fptr, "%-16.9f", (double) time_elapsed / NSEC_PER_SEC);                 // Add time in seconds to file
    stateFreeVector(state);
}

void timeGradientPQC(state_t* state,
                     gradientPQC gradient,
                     const double par[],
                     const obs_t* obs,
                     const obs_t* evoOps[],
                     depth_t circdepth,
                     qubit_t qubits,
                     FILE* fptr)
{
    stateInitPlus(state, qubits);                                                   // Initialize state to plus state
    uint64_t start = get_time();                                                    // Start timing
    gradient(state, par, obs, evoOps, circdepth);                                   // Apply observable
    uint64_t end = get_time();                                                      // End timing
    uint64_t time_elapsed = end - start;                                            // Difference between times
    fprintf(fptr, "%-16.9f", (double) time_elapsed / NSEC_PER_SEC);                 // Add time in seconds to file
    stateFreeVector(state);
}

#endif /* TIMING_H */