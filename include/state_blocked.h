#ifndef STATE_BLOCKED_H
#define STATE_BLOCKED_H

/**
 * @file state_blocked.h
 * @brief Blocked memory layout for density matrices
 *
 * Implements a tile-based blocking of density matrices optimized for:
 * - Cache locality (64x64 tiles fit in L2 cache)
 * - SIMD vectorization (fixed stride within tiles)
 * - OpenMP parallelization (independent tile processing)
 *
 * Recommended for n >= 10 qubits where memory overhead is < 6%.
 */

#include "q_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Elements per block dimension (64x64 x 16 bytes = 64KB per tile) */
#define BLOCK_DIM 64

/**
 * @brief Blocked state structure (parallel to state_t for mixed states)
 */
typedef struct state_blocked {
    cplx_t *data;       /* Blocked storage */
    qubit_t qubits;     /* Number of qubits */
    dim_t n_blocks;     /* Number of tile rows/cols */
    dim_t padded_dim;   /* Dimension padded to BLOCK_DIM multiple */
} state_blocked_t;

/*
 * =====================================================================================================================
 * State management
 * =====================================================================================================================
 */

/**
 * @brief Calculate storage length for blocked format
 *
 * @param[in] qubits Number of qubits
 * @return Number of complex elements needed for blocked storage
 */
dim_t state_blocked_len(qubit_t qubits);

/**
 * @brief Initialize a blocked state (allocates memory, sets to zero)
 *
 * @param[out] state State to initialize
 * @param[in] qubits Number of qubits
 */
void state_blocked_init(state_blocked_t *state, qubit_t qubits);

/**
 * @brief Free blocked state memory
 *
 * @param[in,out] state State to free
 */
void state_blocked_free(state_blocked_t *state);

/**
 * @brief Initialize blocked state to uniform superposition |+⟩⟨+|
 *
 * @param[out] state State to initialize
 * @param[in] qubits Number of qubits
 */
void state_blocked_plus(state_blocked_t *state, qubit_t qubits);

/*
 * =====================================================================================================================
 * Index helpers
 * =====================================================================================================================
 */

/**
 * @brief Get tile offset in storage (column-major lower triangular)
 *
 * @param[in] tile_row Row index of tile
 * @param[in] tile_col Column index of tile (must be <= tile_row)
 * @param[in] n_blocks Number of tile rows/cols
 * @return Offset to start of tile in storage array
 */
dim_t tile_offset(dim_t tile_row, dim_t tile_col, dim_t n_blocks);

/**
 * @brief Get pointer to tile data
 *
 * @param[in] state Blocked state
 * @param[in] tile_row Row index of tile
 * @param[in] tile_col Column index of tile (must be <= tile_row)
 * @return Pointer to first element of tile
 */
cplx_t* get_tile(const state_blocked_t *state, dim_t tile_row, dim_t tile_col);

/*
 * =====================================================================================================================
 * Element access (for debugging and measurement)
 * =====================================================================================================================
 */

/**
 * @brief Get element from blocked state (handles Hermitian symmetry)
 *
 * @param[in] state Blocked state
 * @param[in] row Row index in density matrix
 * @param[in] col Column index in density matrix
 * @return Element value (conjugated if row < col)
 */
cplx_t state_blocked_get(const state_blocked_t *state, dim_t row, dim_t col);

/**
 * @brief Set element in blocked state (only lower triangle)
 *
 * @param[in,out] state Blocked state
 * @param[in] row Row index (must be >= col)
 * @param[in] col Column index
 * @param[in] val Value to set
 */
void state_blocked_set(state_blocked_t *state, dim_t row, dim_t col, cplx_t val);

/**
 * @brief Get probability of measuring a basis state
 *
 * @param[in] state Blocked state
 * @param[in] basis_state Basis state index
 * @return Probability (diagonal element)
 */
double state_blocked_prob(const state_blocked_t *state, dim_t basis_state);

/**
 * @brief Compute trace of density matrix
 *
 * @param[in] state Blocked state
 * @return Trace (should be 1.0 for valid density matrix)
 */
double state_blocked_trace(const state_blocked_t *state);

#ifdef __cplusplus
}
#endif

#endif /* STATE_BLOCKED_H */
