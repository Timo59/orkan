/**
 * @file gate_packed_cp.c
 * @brief CP (controlled-phase) gate for mixed states in packed lower-triangular storage
 */

#include "gate.h"

void cp_packed(state_t *state, const qubit_t control, const qubit_t target, const double theta) {
    GATE_VALIDATE(0, "cp: packed not yet implemented");
}
