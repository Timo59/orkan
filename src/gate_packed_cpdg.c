/**
 * @file gate_packed_cpdg.c
 * @brief CPdg (controlled-phase-dagger) gate for mixed states in packed lower-triangular storage
 */

#include "gate.h"

void cpdg_packed(state_t *state, const qubit_t control, const qubit_t target, const double theta) {
    GATE_VALIDATE(0, "cpdg: packed not yet implemented");
}
