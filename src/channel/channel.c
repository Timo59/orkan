// channel.c  - Channel dispatchers with input validation and two way type dispatch

#include "channel.h"

/*
 * =====================================================================================================================
 * Compute superoperator from Kraus representation
 * =====================================================================================================================
 */

superop_t kraus_to_superop(const kraus_t *kraus) {
  if (!kraus || !kraus->data) {
    fprintf(stderr, "kraus_to_superop: null kraus or data\n");
    exit(EXIT_FAILURE);
  }

  const idx_t N = (idx_t)1 << kraus->n_qubits;  // Hilbert space dimension
  const idx_t M = N * N;                        // Superoperator dimension

  superop_t out = {0};
  out.n_qubits = kraus->n_qubits;
  out.data = calloc(M * M, sizeof(cplx_t));
  if (!out.data) {
    fprintf(stderr, "kraus_to_superop: data allocation failed\n");
    exit(EXIT_FAILURE);
  }

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
 * Packed mixed state channel implementations
 * =====================================================================================================================
 */

extern void channel_packed_1q(state_t *state, const cplx_t *sop, const qubit_t target);

/*
 * =====================================================================================================================
 * Tiled mixed state channel implementations 
 * =====================================================================================================================
 */

extern void channel_tiled_1q(state_t *state, const cplx_t *sop, const qubit_t target);

/*
 * =====================================================================================================================
 * Dispatchers
 * =====================================================================================================================
 */

void channel_1q(state_t *state, const superop_t *sop, const qubit_t target) {
  if (!state || !state->data) {
    fprintf(stderr, "channel_1q: null state or data pointer\n");
    exit(EXIT_FAILURE);
  }
  if (target >= state->qubits) {
    fprintf(stderr, "channel_1q: target qubit out of range\n");
    exit(EXIT_FAILURE);
  }
  if (!sop || !sop->data) {
    fprintf(stderr, "channel_1q: null superoperator or data pointer\n");
    exit(EXIT_FAILURE);
  }
  if (sop->n_qubits != 1) {
    fprintf(stderr, "channel_1q: number of targeted qubits is not 1\n");
    exit(EXIT_FAILURE);
  }

  switch (state->type) {
    case PURE:
      fprintf(stderr, "channel_1q: Kraus operation is not available for pure states\n");
      exit(EXIT_FAILURE);
    case MIXED_PACKED:
      channel_packed_1q(state, sop->data, target);
      break;
    case MIXED_TILED:
      channel_tiled_1q(state, sop->data, target);
      break;
  }
}

