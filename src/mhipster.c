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
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(__APPLE__)
    #include <vecLib/cblas_new.h>
#elif defined(__linux__)
    #include <cblas.h>
#endif

/* Minimum dimension to enable OpenMP parallelization (avoid thread overhead for small systems)
 * Mixed states: O(dim²) work per gate due to density matrix size, so parallelization is
 * beneficial at much smaller dim than for pure states. At dim=64 (n=6), the packed array
 * has ~2000 elements — sufficient work to offset thread overhead.
 */
#define OMP_THRESHOLD 64


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
 *
 * OpenMP parallelization:
 *   The outer col_block loop iterates over independent column blocks. Each block accesses
 *   disjoint memory regions, so threads won't conflict. Parallelization is enabled when
 *   dim >= OMP_THRESHOLD and there are multiple col_blocks (subdim < dim).
 */
#define TRAVERSE_MIXED_1Q(state, target, DIAG_OP, CONJ_OP, OFFDIAG_OP)                                                 \
do {                                                                                                                   \
    const dim_t dim = (dim_t)1 << (state)->qubits;                                                                     \
    const dim_t incr = (dim_t)1 << (target);                                                                           \
    const dim_t subdim = (dim_t)1 << ((target) + 1);                                                                   \
    const int parallel_outer = (subdim < dim) && (dim >= OMP_THRESHOLD);                                               \
                                                                                                                       \
    /* Iterate column blocks: columns where bit t = 0 */                                                               \
    _Pragma("omp parallel for if(parallel_outer)")                                                                     \
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
            CONJ_OP(d10, n_first, col - col_block, dim, col_block);                                                    \
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


/*
 * =====================================================================================================================
 * Hadamard gate: ρ → HρH (H = H†)
 * =====================================================================================================================
 *
 * Hadamard is a "mixing" gate: all 4 quadrants of each 2×2 block mix together.
 * This requires reading all 4 elements before computing the transformation.
 *
 * For each 2×2 block indexed by (r, c) with bit t = 0 in both:
 *   ρ'_00 = (ρ_00 + ρ_01 + ρ_10 + ρ_11) / 2
 *   ρ'_11 = (ρ_00 - ρ_01 - ρ_10 + ρ_11) / 2
 *   ρ'_01 = (ρ_00 - ρ_01 + ρ_10 - ρ_11) / 2
 *   ρ'_10 = (ρ_00 + ρ_01 - ρ_10 - ρ_11) / 2
 */

/** @brief Compute packed array index for element (r, c) where r >= c */
static inline dim_t pack_idx(dim_t dim, dim_t r, dim_t c) {
    return c * (2 * dim - c + 1) / 2 + (r - c);
}

void h_mixed(state_t *state, const qubit_t target) {
    const dim_t dim = (dim_t)1 << state->qubits;
    const dim_t incr = (dim_t)1 << target;
    cplx_t *data = state->data;

    /*
     * Iterate over all 2×2 blocks indexed by (r, c) with r >= c, bit t = 0 in both.
     * Each 2×2 block transformation is independent, so we can parallelize.
     * Use dynamic scheduling since the inner loop length varies with c.
     */
    #pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
    for (dim_t c = 0; c < dim; ++c) {
        if (c & incr) continue;
        for (dim_t r = c; r < dim; ++r) {
            if (r & incr) continue;

            dim_t r0 = r, c0 = c, r1 = r | incr, c1 = c | incr;

            /* Packed indices for the 3 always-lower-triangle elements */
            dim_t idx00 = pack_idx(dim, r0, c0);
            dim_t idx10 = pack_idx(dim, r1, c0);
            dim_t idx11 = pack_idx(dim, r1, c1);

            /* Read 4 elements (element (r0, c1) may be in upper triangle) */
            cplx_t rho00 = data[idx00];
            cplx_t rho10 = data[idx10];
            cplx_t rho11 = data[idx11];
            cplx_t rho01;
            dim_t idx01;
            int r0_ge_c1 = (r0 >= c1);
            if (r0_ge_c1) {
                idx01 = pack_idx(dim, r0, c1);
                rho01 = data[idx01];
            } else {
                idx01 = pack_idx(dim, c1, r0);
                rho01 = conj(data[idx01]);
            }

            /* Hadamard transformation: HρH */
            cplx_t new00 = 0.5 * (rho00 + rho01 + rho10 + rho11);
            cplx_t new11 = 0.5 * (rho00 - rho01 - rho10 + rho11);
            cplx_t new01 = 0.5 * (rho00 - rho01 + rho10 - rho11);
            cplx_t new10 = 0.5 * (rho00 + rho01 - rho10 - rho11);

            /* Write to lower triangle positions */
            data[idx00] = new00;
            data[idx10] = new10;
            data[idx11] = new11;
            if (r0_ge_c1) {
                data[idx01] = new01;
            } else {
                data[idx01] = conj(new01);
            }
        }
    }
}


/*
 * =====================================================================================================================
 * S gate (Phase gate): ρ → SρS†
 * =====================================================================================================================
 *
 * S = diag(1, i), S† = diag(1, -i)
 * Diagonal blocks (0,0), (1,1) unchanged: phases cancel
 * Off-diagonal: (1,0) block *= i, (0,1) block *= -i
 */

/** @brief Multiply n elements by i */
static inline void op_scale_i(cplx_t *a, dim_t n) {
    const cplx_t alpha = I;
    cblas_zscal(n, &alpha, a, 1);
}

/** @brief Multiply n elements by -i */
static inline void op_scale_neg_i(cplx_t *a, dim_t n) {
    const cplx_t alpha = -I;
    cblas_zscal(n, &alpha, a, 1);
}

/** @brief Scale d10 by i and d01 by -i */
static inline void op_scale_i_neg_i(cplx_t *d10, cplx_t *d01, dim_t n) {
    const cplx_t alpha_i = I;
    const cplx_t alpha_neg_i = -I;
    cblas_zscal(n, &alpha_i, d10, 1);
    cblas_zscal(n, &alpha_neg_i, d01, 1);
}

#define S_DIAG_OP(data, n, stride)                      op_noop_diag(data, n, stride)
#define S_CONJ_OP(d10, count, start_k, dim, col_block)  op_scale_i(d10, incr)
#define S_OFFDIAG_OP(d10, d01, n)                       op_scale_i_neg_i(d10, d01, n)

void s_mixed(state_t *state, const qubit_t target) {
    TRAVERSE_MIXED_1Q(state, target, S_DIAG_OP, S_CONJ_OP, S_OFFDIAG_OP);
}


/*
 * =====================================================================================================================
 * S-dagger gate: ρ → S†ρS
 * =====================================================================================================================
 *
 * S† = diag(1, -i), S = diag(1, i)
 * Diagonal blocks unchanged
 * Off-diagonal: (1,0) block *= -i, (0,1) block *= i
 */

#define SDG_DIAG_OP(data, n, stride)                      op_noop_diag(data, n, stride)
#define SDG_CONJ_OP(d10, count, start_k, dim, col_block)  op_scale_neg_i(d10, incr)
#define SDG_OFFDIAG_OP(d10, d01, n)                       op_scale_i_neg_i(d01, d10, n)  /* swapped: d01 by i, d10 by -i */

void sdg_mixed(state_t *state, const qubit_t target) {
    TRAVERSE_MIXED_1Q(state, target, SDG_DIAG_OP, SDG_CONJ_OP, SDG_OFFDIAG_OP);
}


/*
 * =====================================================================================================================
 * T gate (π/8 gate): ρ → TρT†
 * =====================================================================================================================
 *
 * T = diag(1, e^(iπ/4)), T† = diag(1, e^(-iπ/4))
 * Diagonal blocks unchanged: phases cancel
 * Off-diagonal: (1,0) block *= e^(iπ/4), (0,1) block *= e^(-iπ/4)
 */

/** @brief Scale d10 by e^(iπ/4) and d01 by e^(-iπ/4) */
static inline void op_scale_t_phases(cplx_t *d10, cplx_t *d01, dim_t n) {
    const cplx_t alpha_t = M_SQRT1_2 + M_SQRT1_2 * I;      /* e^(iπ/4) */
    const cplx_t alpha_tdg = M_SQRT1_2 - M_SQRT1_2 * I;    /* e^(-iπ/4) */
    cblas_zscal(n, &alpha_t, d10, 1);
    cblas_zscal(n, &alpha_tdg, d01, 1);
}

static inline void op_scale_t(cplx_t *a, dim_t n) {
    const cplx_t alpha = M_SQRT1_2 + M_SQRT1_2 * I;
    cblas_zscal(n, &alpha, a, 1);
}

#define T_DIAG_OP(data, n, stride)                      op_noop_diag(data, n, stride)
#define T_CONJ_OP(d10, count, start_k, dim, col_block)  op_scale_t(d10, incr)
#define T_OFFDIAG_OP(d10, d01, n)                       op_scale_t_phases(d10, d01, n)

void t_mixed(state_t *state, const qubit_t target) {
    TRAVERSE_MIXED_1Q(state, target, T_DIAG_OP, T_CONJ_OP, T_OFFDIAG_OP);
}


/*
 * =====================================================================================================================
 * T-dagger gate: ρ → T†ρT
 * =====================================================================================================================
 *
 * T† = diag(1, e^(-iπ/4)), T = diag(1, e^(iπ/4))
 * Diagonal blocks unchanged
 * Off-diagonal: (1,0) block *= e^(-iπ/4), (0,1) block *= e^(iπ/4)
 */

static inline void op_scale_tdg(cplx_t *a, dim_t n) {
    const cplx_t alpha = M_SQRT1_2 - M_SQRT1_2 * I;
    cblas_zscal(n, &alpha, a, 1);
}

#define TDG_DIAG_OP(data, n, stride)                      op_noop_diag(data, n, stride)
#define TDG_CONJ_OP(d10, count, start_k, dim, col_block)  op_scale_tdg(d10, incr)
#define TDG_OFFDIAG_OP(d10, d01, n)                       op_scale_t_phases(d01, d10, n)  /* swapped phases */

void tdg_mixed(state_t *state, const qubit_t target) {
    TRAVERSE_MIXED_1Q(state, target, TDG_DIAG_OP, TDG_CONJ_OP, TDG_OFFDIAG_OP);
}
