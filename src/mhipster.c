/**
 * @file mhipster.c
 * @brief Quantum gates for mixed states (density matrices in packed lower-triangular storage)
 *
 * Mixed states are stored in LAPACK packed lower-triangular column-major format.
 * For an N×N density matrix, we store N(N+1)/2 complex elements.
 *
 * Single-qubit gates transform ρ → UρU†. For target qubit t, this mixes elements
 * based on bit t of their row/column indices. The traversal visits all elements
 * systematically, processing four logical quadrants:
 *   - (0,0): row bit t = 0, col bit t = 0
 *   - (1,1): row bit t = 1, col bit t = 1
 *   - (1,0): row bit t = 1, col bit t = 0  (in lower triangle)
 *   - (0,1): row bit t = 0, col bit t = 1  (in upper triangle, accessed via conjugate)
 */

#include "gate.h"
#include <complex.h>

#if defined(__APPLE__)
    #include <vecLib/cblas_new.h>
#elif defined(__linux__)
    #include <cblas.h>
#endif


/*
 * =====================================================================================================================
 * Single-qubit gate traversal for packed lower-triangular density matrices
 * =====================================================================================================================
 *
 * This macro generates the loop structure for applying any single-qubit gate to a mixed state.
 * The caller provides operation macros that are expanded at three key points:
 *
 *   DIAG_OP(data, n, stride):
 *       Process n elements from (0,0) quadrant starting at data[0]
 *       paired with n elements from (1,1) quadrant starting at data[stride].
 *       Called once per column in the first row-block, and once per row-block thereafter.
 *
 *   CONJ_OP(d10, count, start_k, dim, col_block):
 *       Process the triangular region where (1,0) elements pair with (0,1) elements.
 *       The (0,1) elements are in the upper triangle, so we access them as conjugates
 *       of their transpose. d10 points to first (1,0) element; count = number of pairs.
 *
 *   OFFDIAG_OP(d10, d01, n):
 *       Process n elements from (1,0) quadrant starting at d10
 *       paired with n elements from (0,1) quadrant starting at d01.
 *       Unlike CONJ_OP, both regions are in lower triangle here.
 *
 * Variables available to operation macros:
 *   dim     - Hilbert space dimension (2^n_qubits)
 *   incr    - Distance between indices differing only in target bit (2^target)
 *   subdim  - Invariant subspace dimension (2^(target+1))
 */
#define TRAVERSE_MIXED_1Q(state, target, DIAG_OP, CONJ_OP, OFFDIAG_OP)                                                 \
do {                                                                                                                   \
    const dim_t dim = (dim_t)1 << (state)->qubits;                                                                     \
    const dim_t incr = (dim_t)1 << (target);                                                                           \
    const dim_t subdim = (dim_t)1 << ((target) + 1);                                                                   \
                                                                                                                       \
    /* Iterate column blocks: columns where bit t = 0 */                                                               \
    for (dim_t col_block = 0; col_block < dim; col_block += subdim) {                                                  \
        /* Offset to start of col_block in packed array */                                                             \
        dim_t offset = col_block * (2 * dim - col_block + 1) / 2;                                                      \
        /* Distance from (i,col) to (i+incr, col+incr) in packed format */                                             \
        dim_t stride = incr * (2 * (dim - col_block) - incr + 1) / 2;                                                  \
                                                                                                                       \
        /* Process each column where bit t = 0 */                                                                      \
        for (dim_t col = col_block; col < col_block + incr; ++col) {                                                   \
            cplx_t *data = (state)->data + offset;                                                                     \
                                                                                                                       \
            /* --- First row block: rows [col, col_block + incr) --- */                                                \
            const dim_t n_first = incr - (col - col_block);                                                            \
                                                                                                                       \
            /* (0,0) ↔ (1,1) pairs in first block */                                                                   \
            DIAG_OP(data, n_first, stride);                                                                            \
                                                                                                                       \
            /* (1,0) ↔ (0,1) pairs where (0,1) is in upper triangle (access via conjugate) */                          \
            cplx_t *d10 = data + n_first;                                                                              \
            CONJ_OP(d10, n_first, col - col_block, dim, col_block);                                             \
                                                                                                                       \
            /* --- Remaining row blocks --- */                                                                         \
            cplx_t *row_data = d10 + incr;                                                                             \
            for (dim_t row_block = col_block + subdim; row_block < dim; row_block += subdim) {                         \
                /* (0,0) ↔ (1,1) pairs */                                                                              \
                DIAG_OP(row_data, incr, stride);                                                                       \
                                                                                                                       \
                /* (1,0) ↔ (0,1) pairs (both in lower triangle) */                                                     \
                OFFDIAG_OP(row_data + incr, row_data + incr + stride - subdim, incr);                                  \
                                                                                                                       \
                row_data += subdim;                                                                                    \
            }                                                                                                          \
                                                                                                                       \
            offset += dim - col;                                                                                       \
            stride -= incr;                                                                                            \
        }                                                                                                              \
    }                                                                                                                  \
} while (0)


/*
 * =====================================================================================================================
 * Gate-specific inline operations
 * =====================================================================================================================
 */

/**
 * @brief Swap n complex elements: data[0:n] ↔ data[stride:stride+n]
 */
static inline void op_swap_stride(cplx_t *data, dim_t n, dim_t stride) {
    cblas_zswap(n, data, 1, data + stride, 1);
}

/**
 * @brief Swap n complex elements: a[0:n] ↔ b[0:n]
 */
static inline void op_swap(cplx_t *a, cplx_t *b, dim_t n) {
    cblas_zswap(n, a, 1, b, 1);
}

/**
 * @brief Swap and negate n complex elements: a ↔ b, then negate both
 */
static inline void op_swap_neg(cplx_t *a, cplx_t *b, dim_t n) {
    const cplx_t neg = -1.0;
    cblas_zswap(n, a, 1, b, 1);
    cblas_zscal(n, &neg, a, 1);
    cblas_zscal(n, &neg, b, 1);
}

/**
 * @brief Negate n complex elements at two locations
 */
static inline void op_neg_both(cplx_t *a, cplx_t *b, dim_t n) {
    const cplx_t neg = -1.0;
    cblas_zscal(n, &neg, a, 1);
    cblas_zscal(n, &neg, b, 1);
}

/**
 * @brief Swap and conjugate pairs in the triangular (1,0)↔(0,1) region.
 *
 * Elements (i,j) with row bit t=1 and col bit t=0 must be swapped with (j',i')
 * where j' has bit t=1 and i' has bit t=0. When j' > i', this element is in the
 * upper triangle; we access it as conj(ρ[i',j']) from the lower triangle.
 *
 * The step between d10[k] and its pair element increases as we iterate because
 * the paired element is progressively further down in the packed array.
 *
 * @param d10        Pointer to first (1,0) element in this region
 * @param count      Number of pairs to process
 * @param start_k    Starting index offset (loop accesses d10[start_k], d10[start_k+1], ...)
 * @param dim        Hilbert space dimension
 * @param col_block  Current column block start
 */
static inline void op_swap_conj(cplx_t *d10, dim_t count, dim_t start_k, dim_t dim, dim_t col_block) {
    dim_t step = 0;
    for (dim_t k = start_k; k < start_k + count; ++k) {
        cplx_t tmp = d10[k];
        d10[k] = conj(d10[k + step]);
        d10[k + step] = conj(tmp);
        step += dim - (k + col_block + 2);
    }
}

/**
 * @brief Swap, conjugate, and negate pairs in the triangular region (for Y gate).
 */
static inline void op_swap_conj_neg(cplx_t *d10, dim_t count, dim_t start_k, dim_t dim, dim_t col_block) {
    dim_t step = 0;
    for (dim_t k = start_k; k < start_k + count; ++k) {
        cplx_t tmp = d10[k];
        d10[k] = -conj(d10[k + step]);
        d10[k + step] = -conj(tmp);
        step += dim - (k + col_block + 2);
    }
}

/**
 * @brief Negate n elements starting at d10 (for Z gate CONJ region).
 *
 * Unlike X/Y which process n_first elements with variable-stride pairs,
 * Z negates incr contiguous elements (the full (1,0) region in this column).
 */
static inline void op_neg_n(cplx_t *d10, dim_t n) {
    const cplx_t neg = -1.0;
    cblas_zscal(n, &neg, d10, 1);
}

/** @brief No-op for diagonal blocks (used by Z gate) */
static inline void op_noop_diag(cplx_t *data, dim_t n, dim_t stride) {
    (void)data; (void)n; (void)stride;
}


/*
 * =====================================================================================================================
 * Pauli-X gate: ρ → XρX
 * =====================================================================================================================
 *
 * X = |0⟩⟨1| + |1⟩⟨0| swaps amplitudes where target bit = 0 with target bit = 1.
 * For density matrix: swaps (0,0)↔(1,1) blocks and swaps (0,1)↔(1,0) blocks.
 */

#define X_DIAG_OP(data, n, stride)                      op_swap_stride(data, n, stride)
#define X_CONJ_OP(d10, count, start_k, dim, col_block)  op_swap_conj(d10, count, start_k, dim, col_block)
#define X_OFFDIAG_OP(d10, d01, n)                       op_swap(d10, d01, n)

void x_mixed(state_t *state, const qubit_t target) {
    TRAVERSE_MIXED_1Q(state, target, X_DIAG_OP, X_CONJ_OP, X_OFFDIAG_OP);
}


/*
 * =====================================================================================================================
 * Pauli-Y gate: ρ → YρY†
 * =====================================================================================================================
 *
 * Y = -i|0⟩⟨1| + i|1⟩⟨0|. The transformation YρY† swaps (0,0)↔(1,1) like X,
 * but the off-diagonal blocks (0,1)↔(1,0) also get negated because (-i)(i) = 1
 * on diagonal but (-i)(i)* = -1 on off-diagonal crossings.
 */

#define Y_DIAG_OP(data, n, stride)                      op_swap_stride(data, n, stride)
#define Y_CONJ_OP(d10, count, start_k, dim, col_block)  op_swap_conj_neg(d10, count, start_k, dim, col_block)
#define Y_OFFDIAG_OP(d10, d01, n)                       op_swap_neg(d10, d01, n)

void y_mixed(state_t *state, const qubit_t target) {
    TRAVERSE_MIXED_1Q(state, target, Y_DIAG_OP, Y_CONJ_OP, Y_OFFDIAG_OP);
}


/*
 * =====================================================================================================================
 * Pauli-Z gate: ρ → ZρZ
 * =====================================================================================================================
 *
 * Z = |0⟩⟨0| - |1⟩⟨1| is diagonal, so ZρZ only affects off-diagonal blocks:
 *   - (0,0) and (1,1) blocks: unchanged (phase factors cancel)
 *   - (0,1) and (1,0) blocks: negated (mixed phase factors: 1 × -1 = -1)
 */

#define Z_DIAG_OP(data, n, stride)                      op_noop_diag(data, n, stride)
#define Z_CONJ_OP(d10, count, start_k, dim, col_block)  op_neg_n(d10, incr)  /* uses incr from outer scope */
#define Z_OFFDIAG_OP(d10, d01, n)                       op_neg_both(d10, d01, n)

void z_mixed(state_t *state, const qubit_t target) {
    TRAVERSE_MIXED_1Q(state, target, Z_DIAG_OP, Z_CONJ_OP, Z_OFFDIAG_OP);
}
