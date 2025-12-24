// mhipster.c - Functions representing fundamental quantum gates applied to mixed states, i.e., density matrices on
// Hilbert spaces

#ifndef GATE_H
#include "gate.h"
#endif

#include <math.h>

#include <stdio.h>

#if defined(__APPLE__)
    #include <vecLib/cblas_new.h>
    #include <vecLib/lapack.h>
#elif defined(__linux__)
    #include <cblas.h>
#endif

/*
 * =====================================================================================================================
 * Pauli gates
 * =====================================================================================================================
 */

void x_mixed(state_t *state, qubit_t target) {

}