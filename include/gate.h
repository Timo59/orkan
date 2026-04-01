// gate.h   - Functions representing the action of local quantum gates to a general quantum state
//
// It includes <stdio.h> and <stdlib.h> because GATE_VALIDATE needs fprintf() and exit().

#ifndef GATE_H
#define GATE_H

/*
 * =====================================================================================================================
 * Includes
 * =====================================================================================================================
 */
#include "state.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/*
 * =====================================================================================================================
 * C++ check
 * =====================================================================================================================
 */
#ifdef __cplusplus
extern "C" {
#endif


/*
 * =====================================================================================================================
 * Error handling
 * =====================================================================================================================
 *
 * GATE_VALIDATE calls exit() on failure. This is a deliberate design choice for a
 * scientific computing library where correctness is paramount: passing an out-of-range
 * qubit index or null pointer indicates a programming error (not a recoverable runtime
 * condition), and continuing would silently produce wrong simulation results. The
 * fail-fast behavior ensures such bugs surface immediately during development and testing.
 */

#define GATE_VALIDATE(cond, msg) do {                           \
    if (!(cond)) {                                              \
        fprintf(stderr, "qlib: gate: %s\n", (msg));            \
        exit(EXIT_FAILURE);                                     \
    }                                                           \
} while(0)

/*
 * =====================================================================================================================
 * Pauli gates
 * =====================================================================================================================
 */

void x(state_t* state, qubit_t target);
void y(state_t* state, qubit_t target);
void z(state_t* state, qubit_t target);

/*
 * =====================================================================================================================
 * Clifford gates
 * =====================================================================================================================
 */

void h(state_t* state, qubit_t target);
void s(state_t* state, qubit_t target);
void sdg(state_t* state, qubit_t target);
void t(state_t* state, qubit_t target);
void tdg(state_t* state, qubit_t target);
void hy(state_t* state, qubit_t target);

/*
 * =====================================================================================================================
 * Rotation gates
 * =====================================================================================================================
 */

void rx(state_t* state, qubit_t target, double theta);
void ry(state_t* state, qubit_t target, double theta);
void rz(state_t* state, qubit_t target, double theta);
void p(state_t* state, qubit_t target, double theta);

/*
 * =====================================================================================================================
 * Two-qubit gates
 * =====================================================================================================================
 */

void cx(state_t* state, qubit_t control, qubit_t target);
void cy(state_t* state, qubit_t control, qubit_t target);
void cz(state_t* state, qubit_t control, qubit_t target);
void swap_gate(state_t* state, qubit_t q1, qubit_t q2);

/*
 * =====================================================================================================================
 * Three-qubit gates
 * =====================================================================================================================
 */

void ccx(state_t* state, qubit_t ctrl1, qubit_t ctrl2, qubit_t target);

#ifdef __cplusplus
}
#endif

#endif // GATE_H
