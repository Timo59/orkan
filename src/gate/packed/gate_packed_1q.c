/**
 * @file gate_packed_1q.c
 * @brief Single-qubit gates for mixed states (density matrices in packed lower-triangular storage)
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
 *   - Auto-vectorization-friendly code structure
 */

#include "gate.h"
#include "index.h"
#include <complex.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Minimum dimension to enable OpenMP parallelization (avoid thread overhead for small systems)
 * Mixed states: O(dim²) work per gate due to density matrix size, so parallelization is
 * beneficial at much smaller dim than for pure states. At dim=512 (n=9), the packed array
 * has ~131k elements — sufficient work to offset thread overhead. Below this, thread-launch
 * overhead (~30 µs) dominates the actual computation.
 */


/*
 * =====================================================================================================================
 * Index computation helpers
 * =====================================================================================================================
 */

/* insertBit0() and pack_idx() are defined in gate.h as static inlines shared across all gate backends. */


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
    const idx_t dim = (idx_t)1 << (state)->qubits;                    \
    const idx_t incr = (idx_t)1 << (target);                          \
    const idx_t half_dim = dim >> 1;                                       \
    cplx_t * restrict data = (state)->data;                                     \
                                                                                \
    _Pragma("omp parallel for schedule(static) if(dim >= OMP_THRESHOLD_MIXED)")       \
    for (idx_t bc = 0; bc < half_dim; ++bc) {                              \
        idx_t c0 = insertBit0(bc, target);                                 \
        idx_t c1 = c0 | incr;                                              \
                                                                                \
        /* Precompute column base for c1: pack_idx(dim, c1, c1) */             \
        idx_t col1_base = c1 * (2 * dim - c1 + 1) / 2;                    \
        /* Precompute column base for c0: pack_idx(dim, c0, c0) */             \
        idx_t col0_base = c0 * (2 * dim - c0 + 1) / 2;                    \
                                                                                \
        for (idx_t br = bc; br < half_dim; ++br) {                         \
            idx_t r0 = insertBit0(br, target);                             \
            idx_t r1 = r0 | incr;                                          \
                                                                                \
            /* Compute idx00 = pack_idx(dim, r0, c0) from col0_base */         \
            idx_t idx00 = col0_base + (r0 - c0);                           \
            /* idx10 = pack_idx(dim, r1, c0): same column, row offset by incr */\
            idx_t idx10 = idx00 + incr;                                    \
            /* idx11 = pack_idx(dim, r1, c1) from col1_base */                 \
            idx_t idx11 = col1_base + (r1 - c1);                           \
                                                                                \
            /* (r0, c1): in lower triangle if r0 >= c1 */                       \
            int lower_01 = (r0 >= c1);                                          \
            idx_t idx01 = lower_01 ? (col1_base + (r0 - c1))               \
                                        : pack_idx(dim, c1, r0);               \
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

void x_packed(state_t *state, const qubit_t target) {
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

void y_packed(state_t *state, const qubit_t target) {
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

void z_packed(state_t *state, const qubit_t target) {
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

void h_packed(state_t *state, const qubit_t target) {
    TRAVERSE_PACKED_BLOCKS(state, target, H_BLOCK_OP);
}


/*
 * =====================================================================================================================
 * Hy (Hadamard-Y) gate: ρ → Hy ρ Hy†
 * =====================================================================================================================
 *
 * Hy = [[1,-1],[1,1]]/√2, Hy† = [[1,1],[-1,1]]/√2  (NOT self-adjoint)
 *
 * For ρ' = Hy ρ Hy†:
 *   ρ'_00 = (ρ00 - ρ01 - ρ10 + ρ11) / 2
 *   ρ'_01 = (ρ00 + ρ01 - ρ10 - ρ11) / 2
 *   ρ'_10 = (ρ00 - ρ01 + ρ10 - ρ11) / 2
 *   ρ'_11 = (ρ00 + ρ01 + ρ10 + ρ11) / 2
 */

#define HY_BLOCK_OP(data, idx00, idx01, idx10, idx11, lower_01, diag_block) do { \
    cplx_t rho00 = data[idx00];                                                 \
    cplx_t rho11 = data[idx11];                                                 \
    cplx_t rho10 = data[idx10];                                                 \
    cplx_t rho01 = diag_block ? conj(rho10)                                     \
                  : (lower_01 ? data[idx01] : conj(data[idx01]));               \
                                                                                \
    cplx_t new00 = 0.5 * (rho00 - rho01 - rho10 + rho11);                       \
    cplx_t new11 = 0.5 * (rho00 + rho01 + rho10 + rho11);                       \
    cplx_t new01 = 0.5 * (rho00 + rho01 - rho10 - rho11);                       \
    cplx_t new10 = 0.5 * (rho00 - rho01 + rho10 - rho11);                       \
                                                                                \
    data[idx00] = new00;                                                        \
    data[idx11] = new11;                                                        \
    data[idx10] = new10;                                                        \
    if (!diag_block) {                                                          \
        data[idx01] = lower_01 ? new01 : conj(new01);                           \
    }                                                                           \
} while (0)

void hy_packed(state_t *state, const qubit_t target) {
    TRAVERSE_PACKED_BLOCKS(state, target, HY_BLOCK_OP);
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
    if (!diag_block) {                                                          \
        /* (0,1) *= -i: if in lower, (a+bi) -> (b-ai) */                        \
        /* if in upper (stored as conj), we read conj and write conj(-i*val) */ \
        double re01 = creal(data[idx01]), im01 = cimag(data[idx01]);            \
        if (lower_01) {                                                         \
            data[idx01] = CMPLX(im01, -re01);                                   \
        } else {                                                                \
            /* stored = conj(rho01), rho01 = conj(stored) = (re01, -im01) */    \
            /* new_rho01 = -i * rho01 = -i*(re01 - i*im01) = -im01 - i*re01 */ \
            /* store conj(new_rho01) = (-im01, re01) */                         \
            data[idx01] = CMPLX(-im01, re01);                                   \
        }                                                                       \
    }                                                                           \
} while (0)

void s_packed(state_t *state, const qubit_t target) {
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
    if (!diag_block) {                                                          \
        /* (0,1) *= i */                                                        \
        double re01 = creal(data[idx01]), im01 = cimag(data[idx01]);            \
        if (lower_01) {                                                         \
            /* (a+bi)*i = (-b+ai) */                                            \
            data[idx01] = CMPLX(-im01, re01);                                   \
        } else {                                                                \
            /* stored = conj(rho01), rho01 = (re01, -im01) */                   \
            /* new_rho01 = i*(re01 - i*im01) = im01 + i*re01 */                \
            /* store conj(new_rho01) = (im01, -re01) */                         \
            data[idx01] = CMPLX(im01, -re01);                                   \
        }                                                                       \
    }                                                                           \
} while (0)

void sdg_packed(state_t *state, const qubit_t target) {
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
    /* (1,0) *= e^(iπ/4): (a+bi) -> ((a-b)*sqrt(1/2), (a+b)*sqrt(1/2)) */      \
    double re10 = creal(data[idx10]), im10 = cimag(data[idx10]);                \
    data[idx10] = CMPLX((re10 - im10) * M_SQRT1_2,                              \
                        (re10 + im10) * M_SQRT1_2);                             \
                                                                                \
    /* On diagonal blocks, idx01 == idx10, so skip (already processed) */       \
    if (!diag_block) {                                                          \
        /* (0,1) *= e^(-iπ/4): (a+bi) -> ((a+b)*sqrt(1/2), (b-a)*sqrt(1/2)) */\
        double re01 = creal(data[idx01]), im01 = cimag(data[idx01]);            \
        if (lower_01) {                                                         \
            data[idx01] = CMPLX((re01 + im01) * M_SQRT1_2,                      \
                                (im01 - re01) * M_SQRT1_2);                     \
        } else {                                                                \
            /* stored = conj(rho01), rho01 = (re01, -im01) */                   \
            /* new = e^(-iπ/4)*(re01 - i*im01) = ((re01-im01)+(-im01-re01)i)/sqrt(2) */\
            /* store conj = ((re01-im01)/sqrt(2), (im01+re01)/sqrt(2)) */       \
            data[idx01] = CMPLX((re01 - im01) * M_SQRT1_2,                      \
                                (im01 + re01) * M_SQRT1_2);                     \
        }                                                                       \
    }                                                                           \
} while (0)

void t_packed(state_t *state, const qubit_t target) {
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
    /* (1,0) *= e^(-iπ/4): (a+bi) -> ((a+b)*sqrt(1/2), (b-a)*sqrt(1/2)) */     \
    double re10 = creal(data[idx10]), im10 = cimag(data[idx10]);                \
    data[idx10] = CMPLX((re10 + im10) * M_SQRT1_2,                              \
                        (im10 - re10) * M_SQRT1_2);                             \
                                                                                \
    /* On diagonal blocks, idx01 == idx10, so skip (already processed) */       \
    if (!diag_block) {                                                          \
        /* (0,1) *= e^(iπ/4): (a+bi) -> ((a-b)*sqrt(1/2), (a+b)*sqrt(1/2)) */ \
        double re01 = creal(data[idx01]), im01 = cimag(data[idx01]);            \
        if (lower_01) {                                                         \
            data[idx01] = CMPLX((re01 - im01) * M_SQRT1_2,                      \
                                (re01 + im01) * M_SQRT1_2);                     \
        } else {                                                                \
            /* stored = conj(rho01), rho01 = (re01, -im01) */                   \
            /* new = e^(iπ/4)*(re01 - i*im01) = ((re01+im01)+(re01-im01)i)/sqrt(2) */\
            /* store conj = ((re01+im01)/sqrt(2), (im01-re01)/sqrt(2)) */       \
            data[idx01] = CMPLX((re01 + im01) * M_SQRT1_2,                      \
                                (im01 - re01) * M_SQRT1_2);                     \
        }                                                                       \
    }                                                                           \
} while (0)

void tdg_packed(state_t *state, const qubit_t target) {
    TRAVERSE_PACKED_BLOCKS(state, target, TDG_BLOCK_OP);
}


/*
 * P(θ) gate: ρ → P(θ)ρP(θ)†. P(θ) = diag(1, e^(iθ)).
 *   - ρ'_00 = ρ00 (unchanged)
 *   - ρ'_11 = ρ11 (unchanged)
 *   - ρ'_01 = e^(-iθ)·ρ01
 *   - ρ'_10 = e^(iθ)·ρ10
 *
 * This is a diagonal gate: only applies phase to off-diagonal elements.
 * e^(-iθ)·(a+bi) = (a·cos(θ) + b·sin(θ)) + i·(b·cos(θ) - a·sin(θ))
 * e^(iθ)·(a+bi) = (a·cos(θ) - b·sin(θ)) + i·(b·cos(θ) + a·sin(θ))
 * The density matrix conjugation is identical to Rz(θ):
 *   ρ'_00 = ρ00, ρ'_11 = ρ11, ρ'_10 = e^(iθ)ρ10, ρ'_01 = e^(-iθ)ρ01
 */
void p_packed(state_t *state, const qubit_t target, const double theta) {
    const double c = cos(theta);
    const double s = sin(theta);

    const idx_t dim = (idx_t)1 << state->qubits;
    const idx_t incr = (idx_t)1 << target;
    const idx_t half_dim = dim >> 1;
    cplx_t * restrict data = state->data;

    #pragma omp parallel for schedule(static) if(dim >= OMP_THRESHOLD_MIXED)
    for (idx_t bc = 0; bc < half_dim; ++bc) {
        idx_t c0 = insertBit0(bc, target);
        idx_t c1 = c0 | incr;

        for (idx_t br = bc; br < half_dim; ++br) {
            idx_t r0 = insertBit0(br, target);
            idx_t r1 = r0 | incr;

            idx_t idx10 = pack_idx(dim, r1, c0);

            int lower_01 = (r0 >= c1);
            idx_t idx01 = lower_01 ? pack_idx(dim, r0, c1) : pack_idx(dim, c1, r0);
            int diag_block = (bc == br);

            /* (1,0) *= e^(iθ): (a+bi) -> (a·c - b·s) + i·(b·c + a·s) */
            double re10 = creal(data[idx10]), im10 = cimag(data[idx10]);
            data[idx10] = CMPLX(re10 * c - im10 * s, im10 * c + re10 * s);

            if (diag_block) continue;

            /* (0,1) *= e^(-iθ): (a+bi) -> (a·c + b·s) + i·(b·c - a·s) */
            double re01 = creal(data[idx01]), im01 = cimag(data[idx01]);
            if (lower_01) {
                data[idx01] = CMPLX(re01 * c + im01 * s, im01 * c - re01 * s);
            } else {
                /* stored = conj(rho01), rho01 = (re01, -im01) */
                /* new = e^(-iθ) * (re01 - i·im01) = (re01·c - im01·s) + i·(-im01·c - re01·s) */
                /* store conj = (re01·c - im01·s) + i·(im01·c + re01·s) */
                data[idx01] = CMPLX(re01 * c - im01 * s, im01 * c + re01 * s);
            }
        }
    }
}


/*
 * =====================================================================================================================
 * Rotation gates
 * =====================================================================================================================
 *
 * For rotation gates Rx, Ry, Rz with angle θ, the transformation is ρ' = U ρ U†.
 * Unlike fixed-angle gates, these require the angle parameter passed at runtime.
 */

/*
 * =====================================================================================================================
 * Rx gate: ρ → Rx(θ) ρ Rx(θ)†
 * =====================================================================================================================
 *
 * Rx(θ) = [[c, -is], [-is, c]] where c = cos(θ/2), s = sin(θ/2)
 *
 * For ρ' = Rx ρ Rx†:
 *   - ρ'_00 = c²·ρ00 + s²·ρ11 + ics·(ρ01 - ρ10)
 *   - ρ'_11 = s²·ρ00 + c²·ρ11 - ics·(ρ01 - ρ10)
 *   - ρ'_01 = ics·(ρ00 - ρ11) + c²·ρ01 + s²·ρ10
 *   - ρ'_10 = -ics·(ρ00 - ρ11) + c²·ρ10 + s²·ρ01
 *
 * Note: On diagonal blocks (bc==br), ρ00 and ρ11 are real diagonal elements and
 * ρ01 - ρ10 is purely imaginary, so ics·(ρ01 - ρ10) is real. On cross-blocks,
 * all elements can be complex.
 */
void rx_packed(state_t *state, const qubit_t target, const double theta) {
    const double c = cos(theta / 2.0);
    const double s = sin(theta / 2.0);
    const double c2 = c * c;
    const double s2 = s * s;
    const double cs = c * s;

    const idx_t dim = (idx_t)1 << state->qubits;
    const idx_t incr = (idx_t)1 << target;
    const idx_t half_dim = dim >> 1;
    cplx_t * restrict data = state->data;

    #pragma omp parallel for schedule(static) if(dim >= OMP_THRESHOLD_MIXED)
    for (idx_t bc = 0; bc < half_dim; ++bc) {
        idx_t c0 = insertBit0(bc, target);
        idx_t c1 = c0 | incr;

        for (idx_t br = bc; br < half_dim; ++br) {
            idx_t r0 = insertBit0(br, target);
            idx_t r1 = r0 | incr;

            idx_t idx00 = pack_idx(dim, r0, c0);
            idx_t idx11 = pack_idx(dim, r1, c1);
            idx_t idx10 = pack_idx(dim, r1, c0);

            int lower_01 = (r0 >= c1);
            idx_t idx01 = lower_01 ? pack_idx(dim, r0, c1) : pack_idx(dim, c1, r0);
            int diag_block = (bc == br);

            cplx_t rho00 = data[idx00];
            cplx_t rho11 = data[idx11];
            cplx_t rho10 = data[idx10];
            cplx_t rho01 = diag_block ? conj(rho10)
                          : (lower_01 ? data[idx01] : conj(data[idx01]));

            /* i*cs*(a+bi) = -cs*b + i*cs*a */
            cplx_t diff_01_10 = rho01 - rho10;
            cplx_t ics_diff_01_10 = CMPLX(-cs * cimag(diff_01_10), cs * creal(diff_01_10));

            cplx_t diff_diag = rho00 - rho11;
            cplx_t ics_diff_diag = CMPLX(-cs * cimag(diff_diag), cs * creal(diff_diag));

            cplx_t new00 = c2 * rho00 + s2 * rho11 + ics_diff_01_10;
            cplx_t new11 = s2 * rho00 + c2 * rho11 - ics_diff_01_10;
            cplx_t new01 = ics_diff_diag + c2 * rho01 + s2 * rho10;
            cplx_t new10 = -ics_diff_diag + c2 * rho10 + s2 * rho01;

            data[idx00] = new00;
            data[idx11] = new11;
            data[idx10] = new10;
            if (!diag_block) {
                data[idx01] = lower_01 ? new01 : conj(new01);
            }
        }
    }
}


/*
 * =====================================================================================================================
 * Ry gate: ρ → Ry(θ) ρ Ry(θ)†
 * =====================================================================================================================
 *
 * Ry(θ) = [[c, -s], [s, c]] where c = cos(θ/2), s = sin(θ/2)
 *
 * For ρ' = Ry ρ Ry†:
 *   - ρ'_00 = c²·ρ00 + s²·ρ11 - cs·(ρ01 + ρ10)
 *   - ρ'_11 = s²·ρ00 + c²·ρ11 + cs·(ρ01 + ρ10)
 *   - ρ'_01 = cs·(ρ00 - ρ11) + c²·ρ01 - s²·ρ10
 *   - ρ'_10 = cs·(ρ00 - ρ11) + c²·ρ10 - s²·ρ01
 *
 * Note: On diagonal blocks, ρ01 + ρ10 = 2·Re(ρ01) for Hermitian matrices.
 * On cross-blocks, all elements can be complex.
 */
void ry_packed(state_t *state, const qubit_t target, const double theta) {
    const double c = cos(theta / 2.0);
    const double s = sin(theta / 2.0);
    const double c2 = c * c;
    const double s2 = s * s;
    const double cs = c * s;

    const idx_t dim = (idx_t)1 << state->qubits;
    const idx_t incr = (idx_t)1 << target;
    const idx_t half_dim = dim >> 1;
    cplx_t * restrict data = state->data;

    #pragma omp parallel for schedule(static) if(dim >= OMP_THRESHOLD_MIXED)
    for (idx_t bc = 0; bc < half_dim; ++bc) {
        idx_t c0 = insertBit0(bc, target);
        idx_t c1 = c0 | incr;

        for (idx_t br = bc; br < half_dim; ++br) {
            idx_t r0 = insertBit0(br, target);
            idx_t r1 = r0 | incr;

            idx_t idx00 = pack_idx(dim, r0, c0);
            idx_t idx11 = pack_idx(dim, r1, c1);
            idx_t idx10 = pack_idx(dim, r1, c0);

            int lower_01 = (r0 >= c1);
            idx_t idx01 = lower_01 ? pack_idx(dim, r0, c1) : pack_idx(dim, c1, r0);
            int diag_block = (bc == br);

            cplx_t rho00 = data[idx00];
            cplx_t rho11 = data[idx11];
            cplx_t rho10 = data[idx10];
            cplx_t rho01 = diag_block ? conj(rho10)
                          : (lower_01 ? data[idx01] : conj(data[idx01]));

            cplx_t sum_01_10 = rho01 + rho10;
            cplx_t diff_diag = rho00 - rho11;

            cplx_t new00 = c2 * rho00 + s2 * rho11 - cs * sum_01_10;
            cplx_t new11 = s2 * rho00 + c2 * rho11 + cs * sum_01_10;
            cplx_t new01 = cs * diff_diag + c2 * rho01 - s2 * rho10;
            cplx_t new10 = cs * diff_diag + c2 * rho10 - s2 * rho01;

            data[idx00] = new00;
            data[idx11] = new11;
            data[idx10] = new10;
            if (!diag_block) {
                data[idx01] = lower_01 ? new01 : conj(new01);
            }
        }
    }
}


/*
 * =====================================================================================================================
 * Rz gate: ρ → Rz(θ) ρ Rz(θ)†
 * =====================================================================================================================
 *
 * Rz(θ) = diag(e^(-iθ/2), e^(iθ/2))
 *
 * For ρ' = Rz ρ Rz†:
 *   - ρ'_00 = ρ00 (unchanged)
 *   - ρ'_11 = ρ11 (unchanged)
 *   - ρ'_01 = e^(-iθ)·ρ01
 *   - ρ'_10 = e^(iθ)·ρ10
 *
 * This is a diagonal gate: only applies phase to off-diagonal elements.
 * e^(-iθ)·(a+bi) = (a·cos(θ) + b·sin(θ)) + i·(b·cos(θ) - a·sin(θ))
 * e^(iθ)·(a+bi) = (a·cos(θ) - b·sin(θ)) + i·(b·cos(θ) + a·sin(θ))
 * Its effect is identical to the P(θ) gate.
 */
void rz_packed(state_t *state, const qubit_t target, const double theta) {
    p_packed(state, target, theta);
}
