// channel.c - Superoperator construction and channel dispatch for mixed states

#include "channel.h"

#include <stdio.h>
#include <stdlib.h>

#define CHANNEL_VALIDATE(cond, msg) do {                        \
    if (!(cond)) {                                              \
        fprintf(stderr, "orkan: channel: %s\n", (msg));        \
        exit(EXIT_FAILURE);                                     \
    }                                                           \
} while(0)

/*
 * =====================================================================================================================
 * Compute superoperator from Kraus representation
 * =====================================================================================================================
 */

superop_t kraus_to_superop(const kraus_t *kraus) {
  CHANNEL_VALIDATE(kraus && kraus->data, "kraus_to_superop: null kraus or data");

  const idx_t N = (idx_t)1 << kraus->n_qubits;  // Hilbert space dimension
  const idx_t M = N * N;                        // Superoperator dimension

  superop_t out = {0};
  out.n_qubits = kraus->n_qubits;
  out.data = calloc(M * M, sizeof(cplx_t));
  CHANNEL_VALIDATE(out.data, "kraus_to_superop: allocation failed");

  for (uint64_t k = 0; k < kraus->n_terms; ++k) {
    /* Pointer to first element in k-th Kraus operator */
    cplx_t * restrict op = kraus->data + k * M;

    for (idx_t block_row = 0; block_row < N; ++block_row) {
      for (idx_t row = 0; row < N; ++row) {
        for (idx_t block_col = 0; block_col < N; ++block_col) {
          idx_t offset = (block_row * N + row) * M + block_col * N;

          for (idx_t col = 0; col < N; ++col) {
            out.data[col + offset] += op[block_col + block_row * N] * conj(op[col + row * N]);
          }
        }
      }
    }
  }

  return out;
}

/*
 * =====================================================================================================================
 * Backend declarations
 * =====================================================================================================================
 */

extern void channel_packed_1q(state_t *state, const cplx_t *sop, const qubit_t target);
extern void channel_tiled_1q(state_t *state, const cplx_t *sop, const qubit_t target);

/*
 * =====================================================================================================================
 * Dispatcher
 * =====================================================================================================================
 */

void channel_1q(state_t *state, const superop_t *sop, const qubit_t target) {
  CHANNEL_VALIDATE(state && state->data,  "channel_1q: null state or data pointer");
  CHANNEL_VALIDATE(target < state->qubits, "channel_1q: target qubit out of range");
  CHANNEL_VALIDATE(sop && sop->data,       "channel_1q: null superoperator or data pointer");
  CHANNEL_VALIDATE(sop->n_qubits == 1,     "channel_1q: superoperator is not 1-local");

  switch (state->type) {
    case PURE:
      CHANNEL_VALIDATE(0, "channel_1q: not available for pure states");
    case MIXED_PACKED:
      channel_packed_1q(state, sop->data, target);
      break;
    case MIXED_TILED:
      channel_tiled_1q(state, sop->data, target);
      break;
    default:
      CHANNEL_VALIDATE(0, "channel_1q: unknown state type");
  }
}
