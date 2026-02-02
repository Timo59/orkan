/**
 * @file mhipster.c
 * @brief Quantum gates for mixed states (density matrices in packed lower-triangular storage)
 *
 * Mixed states are stored in LAPACK packed lower-triangular column-major format.
 * For an N×N density matrix, we store N(N+1)/2 complex elements.
 *
 * Single-qubit gates transform ρ → UρU†. This implementation uses a unified
 * 2×2 block processing approach where all gates iterate over blocks indexed
 * by (r0, c0) with bit t=0 in both indices. Each block contains:
 *   - (0,0) at (r0, c0)
 *   - (1,1) at (r1, c1) where r1=r0|incr, c1=c0|incr
 *   - (1,0) at (r1, c0)
 *   - (0,1) at (r0, c1) - may be in upper triangle
 *
 * Key optimizations:
 *   - Direct enumeration using insertBit0() - no wasted iterations
 *   - Zero-multiplication for phase gates (S, T, etc.)
 *   - All operations inline - no BLAS overhead
 *   - Auto-vectorization-friendly code structure
 */

#include "gate.h"
#include <complex.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Minimum dimension to enable OpenMP parallelization (avoid thread overhead for small systems)
 * Mixed states: O(dim²) work per gate due to density matrix size, so parallelization is
 * beneficial at much smaller dim than for pure states. At dim=64 (n=6), the packed array
 * has ~2000 elements — sufficient work to offset thread overhead.
 */
#define OMP_THRESHOLD 64


/*
 * =====================================================================================================================
 * Index computation helpers
 * =====================================================================================================================
 */

/**
 * @brief Insert a 0 bit at position `pos` in value `val`
 *
 * This allows direct enumeration of indices with bit t=0, avoiding wasted
 * iterations and conditional skipping. For k from 0 to dim/2-1,
 * insertBit0(k, target) produces all indices with bit t=0.
 *
 * Example: pos=2, val=0b101 -> 0b1001 (insert 0 at bit 2)
 */
static inline dim_t insertBit0(dim_t val, qubit_t pos) {
    dim_t mask = ((dim_t)1 << pos) - 1;
    return (val & mask) | ((val & ~mask) << 1);
}

/**
 * @brief Compute packed array index for element (r, c) where r >= c
 */
static inline dim_t pack_idx(dim_t dim, dim_t r, dim_t c) {
    return c * (2 * dim - c + 1) / 2 + (r - c);
}


/*
 * =====================================================================================================================
 * Unified 2×2 block traversal macro
 * =====================================================================================================================
 *
 * Iterates over all 2×2 blocks in the density matrix that mix under the target qubit.
 * Each block is indexed by (r0, c0) where bit t=0 in both.
 *
 * For each block:
 *   - idx00: packed index of (r0, c0) - always in lower triangle
 *   - idx11: packed index of (r1, c1) - always in lower triangle
 *   - idx10: packed index of (r1, c0) - always in lower triangle
 *   - idx01: packed index of (r0, c1) or (c1, r0) if in upper triangle
 *   - lower_01: 1 if (r0, c1) is in lower triangle, 0 if we read conjugate
 *
 * The BLOCK_OP macro receives these indices and a pointer to the data array.
 */
#define TRAVERSE_PACKED_BLOCKS(state, target, BLOCK_OP)                         \
do {                                                                            \
    const dim_t dim = (dim_t)1 << (state)->qubits;                              \
    const dim_t incr = (dim_t)1 << (target);                                    \
    const dim_t half_dim = dim >> 1;                                            \
    cplx_t * restrict data = (state)->data;                                     \
                                                                                \
    _Pragma("omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)")      \
    for (dim_t bc = 0; bc < half_dim; ++bc) {                                   \
        dim_t c0 = insertBit0(bc, target);                                      \
        dim_t c1 = c0 | incr;                                                   \
                                                                                \
        for (dim_t br = bc; br < half_dim; ++br) {                              \
            dim_t r0 = insertBit0(br, target);                                  \
            dim_t r1 = r0 | incr;                                               \
                                                                                \
            /* Compute packed indices for the 4 elements */                     \
            dim_t idx00 = pack_idx(dim, r0, c0);                                \
            dim_t idx11 = pack_idx(dim, r1, c1);                                \
            dim_t idx10 = pack_idx(dim, r1, c0);                                \
                                                                                \
            /* (r0, c1): in lower triangle if r0 >= c1 */                       \
            int lower_01 = (r0 >= c1);                                          \
            dim_t idx01 = lower_01 ? pack_idx(dim, r0, c1)                       \
                                   : pack_idx(dim, c1, r0);                     \
                                                                                \
            /* On diagonal blocks (bc==br), idx01==idx10 points to same elem */ \
            int diag_block = (bc == br);                                        \
                                                                                \
            BLOCK_OP(data, idx00, idx01, idx10, idx11, lower_01, diag_block);   \
        }                                                                       \
    }                                                                           \
} while (0)


/*
 * =====================================================================================================================
 * Pauli gates
 * =====================================================================================================================
 */

/*
 * =====================================================================================================================
 * Pauli-X gate: ρ → XρX
 * =====================================================================================================================
 *
 * X = |0⟩⟨1| + |1⟩⟨0| swaps amplitudes where target bit = 0 with target bit = 1.
 * For density matrix:
 *   - Swap (0,0) ↔ (1,1)
 *   - Swap (0,1) ↔ (1,0)
 *
 * Zero-multiplication gate: only swaps, no arithmetic.
 */

#define X_BLOCK_OP(data, idx00, idx01, idx10, idx11, lower_01, diag_block) do { \
    /* Swap (0,0) <-> (1,1) */                                                  \
    cplx_t tmp = data[idx00];                                                   \
    data[idx00] = data[idx11];                                                  \
    data[idx11] = tmp;                                                          \
                                                                                \
    /* Swap (0,1) <-> (1,0), handling upper triangle */                         \
    if (diag_block) {                                                           \
        /* On diag blocks: stored = rho10, rho01 = conj(rho10) */               \
        /* After swap: new_rho10 = rho01 = conj(stored) */                      \
        /* We store new_rho10 directly (it's still at lower triangle pos) */    \
        data[idx10] = conj(data[idx10]);                                        \
    } else {                                                                    \
        cplx_t rho01 = lower_01 ? data[idx01] : conj(data[idx01]);              \
        cplx_t rho10 = data[idx10];                                             \
        data[idx10] = rho01;                                                    \
        data[idx01] = lower_01 ? rho10 : conj(rho10);                           \
    }                                                                           \
} while (0)

void x_mixed(state_t *state, const qubit_t target) {
    TRAVERSE_PACKED_BLOCKS(state, target, X_BLOCK_OP);
}


/*
 * =====================================================================================================================
 * Pauli-Y gate: ρ → YρY†
 * =====================================================================================================================
 *
 * Y = -i|0⟩⟨1| + i|1⟩⟨0|. The transformation YρY† swaps (0,0)↔(1,1) like X,
 * but the off-diagonal blocks get negated: (0,1) ↔ -(1,0).
 *
 * Zero-multiplication gate: swaps and negations only.
 */

#define Y_BLOCK_OP(data, idx00, idx01, idx10, idx11, lower_01, diag_block) do { \
    /* Swap (0,0) <-> (1,1) */                                                  \
    cplx_t tmp = data[idx00];                                                   \
    data[idx00] = data[idx11];                                                  \
    data[idx11] = tmp;                                                          \
                                                                                \
    /* Swap (0,1) <-> (1,0) with negation */                                    \
    /* On diagonal blocks, rho01 = conj(rho10), so -conj(rho10) is stored */    \
    if (diag_block) {                                                           \
        data[idx10] = -conj(data[idx10]);                                       \
    } else {                                                                    \
        cplx_t rho01 = lower_01 ? data[idx01] : conj(data[idx01]);              \
        cplx_t rho10 = data[idx10];                                             \
        data[idx10] = -rho01;                                                   \
        data[idx01] = lower_01 ? -rho10 : conj(-rho10);                         \
    }                                                                           \
} while (0)

void y_mixed(state_t *state, const qubit_t target) {
    TRAVERSE_PACKED_BLOCKS(state, target, Y_BLOCK_OP);
}


/*
 * =====================================================================================================================
 * Pauli-Z gate: ρ → ZρZ
 * =====================================================================================================================
 *
 * Z = |0⟩⟨0| - |1⟩⟨1| is diagonal, so ZρZ only affects off-diagonal blocks:
 *   - (0,0) and (1,1) blocks: unchanged (phase factors cancel: 1*1=1, (-1)*(-1)=1)
 *   - (0,1) and (1,0) blocks: negated (mixed phase factors: 1 × -1 = -1)
 *
 * Zero-multiplication gate: only negations.
 */

#define Z_BLOCK_OP(data, idx00, idx01, idx10, idx11, lower_01, diag_block) do { \
    (void)idx00; (void)idx11; (void)lower_01;                                   \
    /* Diagonal unchanged, negate off-diagonal */                               \
    data[idx10] = -data[idx10];                                                 \
    /* On diagonal blocks, idx01 == idx10, so skip to avoid double-negation */  \
    if (!diag_block) {                                                          \
        data[idx01] = -data[idx01];                                             \
    }                                                                           \
} while (0)

void z_mixed(state_t *state, const qubit_t target) {
    TRAVERSE_PACKED_BLOCKS(state, target, Z_BLOCK_OP);
}


/*
 * =====================================================================================================================
 * Clifford gates
 * =====================================================================================================================
 */

/*
 * =====================================================================================================================
 * Hadamard gate: ρ → HρH (H = H†)
 * =====================================================================================================================
 *
 * Hadamard is a "mixing" gate: all 4 quadrants of each 2×2 block mix together.
 * The transformation is:
 *   ρ'_00 = (ρ_00 + ρ_01 + ρ_10 + ρ_11) / 2
 *   ρ'_11 = (ρ_00 - ρ_01 - ρ_10 + ρ_11) / 2
 *   ρ'_01 = (ρ_00 - ρ_01 + ρ_10 - ρ_11) / 2
 *   ρ'_10 = (ρ_00 + ρ_01 - ρ_10 - ρ_11) / 2
 */

#define H_BLOCK_OP(data, idx00, idx01, idx10, idx11, lower_01, diag_block) do { \
    cplx_t rho00 = data[idx00];                                                 \
    cplx_t rho11 = data[idx11];                                                 \
    cplx_t rho10 = data[idx10];                                                 \
    /* On diagonal blocks, rho01 = conj(rho10) by hermiticity */                \
    cplx_t rho01 = diag_block ? conj(rho10)                                     \
                  : (lower_01 ? data[idx01] : conj(data[idx01]));               \
                                                                                \
    cplx_t new00 = 0.5 * (rho00 + rho01 + rho10 + rho11);                       \
    cplx_t new11 = 0.5 * (rho00 - rho01 - rho10 + rho11);                       \
    cplx_t new01 = 0.5 * (rho00 - rho01 + rho10 - rho11);                       \
    cplx_t new10 = 0.5 * (rho00 + rho01 - rho10 - rho11);                       \
                                                                                \
    data[idx00] = new00;                                                        \
    data[idx11] = new11;                                                        \
    data[idx10] = new10;                                                        \
    /* On diagonal blocks, new01 = conj(new10) so store conj(new01) = new10 */  \
    if (!diag_block) {                                                          \
        data[idx01] = lower_01 ? new01 : conj(new01);                           \
    }                                                                           \
} while (0)

void h_mixed(state_t *state, const qubit_t target) {
    TRAVERSE_PACKED_BLOCKS(state, target, H_BLOCK_OP);
}



/*
 * =====================================================================================================================
 * S gate (Phase gate): ρ → SρS†
 * =====================================================================================================================
 *
 * S = diag(1, i), S† = diag(1, -i)
 * For ρ' = SρS†:
 *   - (0,0): unchanged (1 * ρ_00 * 1 = ρ_00)
 *   - (1,1): unchanged (i * ρ_11 * (-i) = ρ_11)
 *   - (1,0): *= i (i * ρ_10 * 1 = i*ρ_10)
 *   - (0,1): *= -i (1 * ρ_01 * (-i) = -i*ρ_01)
 *
 * Zero-multiplication: multiply by i is (a+bi) -> (-b+ai)
 *                      multiply by -i is (a+bi) -> (b-ai)
 */

#define S_BLOCK_OP(data, idx00, idx01, idx10, idx11, lower_01, diag_block) do { \
    (void)idx00; (void)idx11;                                                   \
    /* (1,0) *= i: (a+bi) -> (-b+ai) */                                         \
    double re10 = creal(data[idx10]), im10 = cimag(data[idx10]);                \
    data[idx10] = CMPLX(-im10, re10);                                           \
                                                                                \
    /* On diagonal blocks, idx01 == idx10, so skip (already processed) */       \
    if (diag_block) break;                                                      \
                                                                                \
    /* (0,1) *= -i: if in lower, (a+bi) -> (b-ai) */                            \
    /* if in upper (stored as conj), we read conj and write conj(-i * val) */   \
    double re01 = creal(data[idx01]), im01 = cimag(data[idx01]);                \
    if (lower_01) {                                                             \
        data[idx01] = CMPLX(im01, -re01);                                       \
    } else {                                                                    \
        /* stored = conj(rho01), rho01 = conj(stored) = (re01, -im01) */        \
        /* new_rho01 = -i * rho01 = -i * (re01 - i*im01) = -im01 - i*re01 */    \
        /* store conj(new_rho01) = (-im01, re01) */                             \
        data[idx01] = CMPLX(-im01, re01);                                       \
    }                                                                           \
} while (0)

void s_mixed(state_t *state, const qubit_t target) {
    TRAVERSE_PACKED_BLOCKS(state, target, S_BLOCK_OP);
}


/*
 * =====================================================================================================================
 * S-dagger gate: ρ → S†ρS
 * =====================================================================================================================
 *
 * S† = diag(1, -i), S = diag(1, i)
 * For ρ' = S†ρS:
 *   - (0,0): unchanged
 *   - (1,1): unchanged
 *   - (1,0): *= -i
 *   - (0,1): *= i
 */

#define SDG_BLOCK_OP(data, idx00, idx01, idx10, idx11, lower_01, diag_block) do {\
    (void)idx00; (void)idx11;                                                   \
    /* (1,0) *= -i: (a+bi) -> (b-ai) */                                         \
    double re10 = creal(data[idx10]), im10 = cimag(data[idx10]);                \
    data[idx10] = CMPLX(im10, -re10);                                           \
                                                                                \
    /* On diagonal blocks, idx01 == idx10, so skip (already processed) */       \
    if (diag_block) break;                                                      \
                                                                                \
    /* (0,1) *= i */                                                            \
    double re01 = creal(data[idx01]), im01 = cimag(data[idx01]);                \
    if (lower_01) {                                                             \
        /* (a+bi)*i = (-b+ai) */                                                \
        data[idx01] = CMPLX(-im01, re01);                                       \
    } else {                                                                    \
        /* stored = conj(rho01), rho01 = (re01, -im01) */                       \
        /* new_rho01 = i * rho01 = i * (re01 - i*im01) = im01 + i*re01 */       \
        /* store conj(new_rho01) = (im01, -re01) */                             \
        data[idx01] = CMPLX(im01, -re01);                                       \
    }                                                                           \
} while (0)

void sdg_mixed(state_t *state, const qubit_t target) {
    TRAVERSE_PACKED_BLOCKS(state, target, SDG_BLOCK_OP);
}


/*
 * =====================================================================================================================
 * T gate (π/8 gate): ρ → TρT†
 * =====================================================================================================================
 *
 * T = diag(1, e^(iπ/4)), T† = diag(1, e^(-iπ/4))
 * For ρ' = TρT†:
 *   - (0,0): unchanged
 *   - (1,1): unchanged (e^(iπ/4) * ρ_11 * e^(-iπ/4) = ρ_11)
 *   - (1,0): *= e^(iπ/4)
 *   - (0,1): *= e^(-iπ/4)
 *
 * Minimal multiplication:
 *   e^(iπ/4) = (1+i)/√2, so (a+bi)*e^(iπ/4) = ((a-b) + (a+b)i)/√2
 *   e^(-iπ/4) = (1-i)/√2, so (a+bi)*e^(-iπ/4) = ((a+b) + (b-a)i)/√2
 */

#define T_BLOCK_OP(data, idx00, idx01, idx10, idx11, lower_01, diag_block) do { \
    (void)idx00; (void)idx11;                                                   \
    /* (1,0) *= e^(iπ/4): (a+bi) -> ((a-b)*√½, (a+b)*√½) */                     \
    double re10 = creal(data[idx10]), im10 = cimag(data[idx10]);                \
    data[idx10] = CMPLX((re10 - im10) * M_SQRT1_2,                              \
                        (re10 + im10) * M_SQRT1_2);                             \
                                                                                \
    /* On diagonal blocks, idx01 == idx10, so skip (already processed) */       \
    if (diag_block) break;                                                      \
                                                                                \
    /* (0,1) *= e^(-iπ/4): (a+bi) -> ((a+b)*√½, (b-a)*√½) */                    \
    double re01 = creal(data[idx01]), im01 = cimag(data[idx01]);                \
    if (lower_01) {                                                             \
        data[idx01] = CMPLX((re01 + im01) * M_SQRT1_2,                          \
                            (im01 - re01) * M_SQRT1_2);                         \
    } else {                                                                    \
        /* stored = conj(rho01), rho01 = (re01, -im01) */                       \
        /* new = e^(-iπ/4) * (re01 - i*im01) = ((re01-im01) + (-im01-re01)i)/√2*/\
        /* store conj = ((re01-im01)/√2, (im01+re01)/√2) */                     \
        data[idx01] = CMPLX((re01 - im01) * M_SQRT1_2,                          \
                            (im01 + re01) * M_SQRT1_2);                         \
    }                                                                           \
} while (0)

void t_mixed(state_t *state, const qubit_t target) {
    TRAVERSE_PACKED_BLOCKS(state, target, T_BLOCK_OP);
}


/*
 * =====================================================================================================================
 * T-dagger gate: ρ → T†ρT
 * =====================================================================================================================
 *
 * T† = diag(1, e^(-iπ/4)), T = diag(1, e^(iπ/4))
 * For ρ' = T†ρT:
 *   - (0,0): unchanged
 *   - (1,1): unchanged
 *   - (1,0): *= e^(-iπ/4)
 *   - (0,1): *= e^(iπ/4)
 */

#define TDG_BLOCK_OP(data, idx00, idx01, idx10, idx11, lower_01, diag_block) do {\
    (void)idx00; (void)idx11;                                                   \
    /* (1,0) *= e^(-iπ/4): (a+bi) -> ((a+b)*√½, (b-a)*√½) */                    \
    double re10 = creal(data[idx10]), im10 = cimag(data[idx10]);                \
    data[idx10] = CMPLX((re10 + im10) * M_SQRT1_2,                              \
                        (im10 - re10) * M_SQRT1_2);                             \
                                                                                \
    /* On diagonal blocks, idx01 == idx10, so skip (already processed) */       \
    if (diag_block) break;                                                      \
                                                                                \
    /* (0,1) *= e^(iπ/4): (a+bi) -> ((a-b)*√½, (a+b)*√½) */                     \
    double re01 = creal(data[idx01]), im01 = cimag(data[idx01]);                \
    if (lower_01) {                                                             \
        data[idx01] = CMPLX((re01 - im01) * M_SQRT1_2,                          \
                            (re01 + im01) * M_SQRT1_2);                         \
    } else {                                                                    \
        /* stored = conj(rho01), rho01 = (re01, -im01) */                       \
        /* new = e^(iπ/4) * (re01 - i*im01) = ((re01+im01) + (re01-im01)i)/√2 */ \
        /* store conj = ((re01+im01)/√2, (im01-re01)/√2) */                     \
        data[idx01] = CMPLX((re01 + im01) * M_SQRT1_2,                          \
                            (im01 - re01) * M_SQRT1_2);                         \
    }                                                                           \
} while (0)

void tdg_mixed(state_t *state, const qubit_t target) {
    TRAVERSE_PACKED_BLOCKS(state, target, TDG_BLOCK_OP);
}
