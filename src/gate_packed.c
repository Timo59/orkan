/**
 * @file gate_packed.c
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

/* insertBit0() is defined in gate.h as a static inline shared across all gate backends. */

/**
 * @brief Compute packed array index for element (r, c) where r >= c
 */
static inline gate_idx_t pack_idx(gate_idx_t dim, gate_idx_t r, gate_idx_t c) {
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
    const gate_idx_t dim = (gate_idx_t)1 << (state)->qubits;                    \
    const gate_idx_t incr = (gate_idx_t)1 << (target);                          \
    const gate_idx_t half_dim = dim >> 1;                                       \
    cplx_t * restrict data = (state)->data;                                     \
                                                                                \
    _Pragma("omp parallel for schedule(static) if(dim >= OMP_THRESHOLD)")       \
    for (gate_idx_t bc = 0; bc < half_dim; ++bc) {                              \
        gate_idx_t c0 = insertBit0(bc, target);                                 \
        gate_idx_t c1 = c0 | incr;                                              \
                                                                                \
        /* Precompute column base for c1: pack_idx(dim, c1, c1) */             \
        gate_idx_t col1_base = c1 * (2 * dim - c1 + 1) / 2;                    \
        /* Precompute column base for c0: pack_idx(dim, c0, c0) */             \
        gate_idx_t col0_base = c0 * (2 * dim - c0 + 1) / 2;                    \
                                                                                \
        for (gate_idx_t br = bc; br < half_dim; ++br) {                         \
            gate_idx_t r0 = insertBit0(br, target);                             \
            gate_idx_t r1 = r0 | incr;                                          \
                                                                                \
            /* Compute idx00 = pack_idx(dim, r0, c0) from col0_base */         \
            gate_idx_t idx00 = col0_base + (r0 - c0);                           \
            /* idx10 = pack_idx(dim, r1, c0): same column, row offset by incr */\
            gate_idx_t idx10 = idx00 + incr;                                    \
            /* idx11 = pack_idx(dim, r1, c1) from col1_base */                 \
            gate_idx_t idx11 = col1_base + (r1 - c1);                           \
                                                                                \
            /* (r0, c1): in lower triangle if r0 >= c1 */                       \
            int lower_01 = (r0 >= c1);                                          \
            gate_idx_t idx01 = lower_01 ? (col1_base + (r0 - c1))               \
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

    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    const gate_idx_t incr = (gate_idx_t)1 << target;
    const gate_idx_t half_dim = dim >> 1;
    cplx_t * restrict data = state->data;

    #pragma omp parallel for schedule(static) if(dim >= OMP_THRESHOLD)
    for (gate_idx_t bc = 0; bc < half_dim; ++bc) {
        gate_idx_t c0 = insertBit0(bc, target);
        gate_idx_t c1 = c0 | incr;

        for (gate_idx_t br = bc; br < half_dim; ++br) {
            gate_idx_t r0 = insertBit0(br, target);
            gate_idx_t r1 = r0 | incr;

            gate_idx_t idx10 = pack_idx(dim, r1, c0);

            int lower_01 = (r0 >= c1);
            gate_idx_t idx01 = lower_01 ? pack_idx(dim, r0, c1) : pack_idx(dim, c1, r0);
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

    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    const gate_idx_t incr = (gate_idx_t)1 << target;
    const gate_idx_t half_dim = dim >> 1;
    cplx_t * restrict data = state->data;

    #pragma omp parallel for schedule(static) if(dim >= OMP_THRESHOLD)
    for (gate_idx_t bc = 0; bc < half_dim; ++bc) {
        gate_idx_t c0 = insertBit0(bc, target);
        gate_idx_t c1 = c0 | incr;

        for (gate_idx_t br = bc; br < half_dim; ++br) {
            gate_idx_t r0 = insertBit0(br, target);
            gate_idx_t r1 = r0 | incr;

            gate_idx_t idx00 = pack_idx(dim, r0, c0);
            gate_idx_t idx11 = pack_idx(dim, r1, c1);
            gate_idx_t idx10 = pack_idx(dim, r1, c0);

            int lower_01 = (r0 >= c1);
            gate_idx_t idx01 = lower_01 ? pack_idx(dim, r0, c1) : pack_idx(dim, c1, r0);
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

    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    const gate_idx_t incr = (gate_idx_t)1 << target;
    const gate_idx_t half_dim = dim >> 1;
    cplx_t * restrict data = state->data;

    #pragma omp parallel for schedule(static) if(dim >= OMP_THRESHOLD)
    for (gate_idx_t bc = 0; bc < half_dim; ++bc) {
        gate_idx_t c0 = insertBit0(bc, target);
        gate_idx_t c1 = c0 | incr;

        for (gate_idx_t br = bc; br < half_dim; ++br) {
            gate_idx_t r0 = insertBit0(br, target);
            gate_idx_t r1 = r0 | incr;

            gate_idx_t idx00 = pack_idx(dim, r0, c0);
            gate_idx_t idx11 = pack_idx(dim, r1, c1);
            gate_idx_t idx10 = pack_idx(dim, r1, c0);

            int lower_01 = (r0 >= c1);
            gate_idx_t idx01 = lower_01 ? pack_idx(dim, r0, c1) : pack_idx(dim, c1, r0);
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


/*
 * =====================================================================================================================
 * Two-qubit gates
 * =====================================================================================================================
 */

/*
 * =====================================================================================================================
 * CNOT gate: ρ → CX · ρ · CX†  (CX = CX† for CNOT)
 * =====================================================================================================================
 *
 * CNOT is a permutation gate that swaps basis states |ctrl=1,tgt=0> <-> |ctrl=1,tgt=1>.
 * For density matrices ρ' = CX·ρ·CX†, the transformation permutes rows and columns:
 *   new_rho[i,j] = rho[perm(i), perm(j)]
 * where perm(x) = x XOR (1<<target) if (x & (1<<control)), else x.
 *
 * This is a "zero-flop" gate in the sense that CNOT involves no matrix multiplication,
 * only permutation of elements. However, for packed Hermitian storage, some elements
 * require conjugation when they cross the diagonal.
 */
/*
 * In-place CNOT for packed storage.
 *
 * CNOT is an involution: perm(perm(x)) = x, where perm(x) = x ^ target_bit
 * when control_bit is set. Since the density matrix transformation is
 * new_rho[r,c] = old_rho[perm(r), perm(c)], and the permutation is its own
 * inverse, every stored element is either a fixed point or part of a 2-cycle.
 *
 * For each stored (r,c) with r >= c, we compute the destination (r',c') =
 * (perm(r), perm(c)) canonicalized to lower-triangle form. If this differs
 * from (r,c), we swap the two packed elements. To avoid double-swapping,
 * we only swap when the source packed index < destination packed index.
 *
 * When the source and destination are on opposite sides of the diagonal
 * (one has r>=c, the other has r<c), the values must be conjugated.
 */
void cx_packed(state_t *state, const qubit_t control, const qubit_t target) {
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    cplx_t * restrict data = state->data;

    const gate_idx_t incr_ctrl = (gate_idx_t)1 << control;
    const gate_idx_t incr_tgt = (gate_idx_t)1 << target;

    #pragma omp parallel for schedule(static) if(dim >= OMP_THRESHOLD)
    for (gate_idx_t c = 0; c < dim; ++c) {
        for (gate_idx_t r = c; r < dim; ++r) {
            /* Compute permuted row and column */
            gate_idx_t pr = (r & incr_ctrl) ? (r ^ incr_tgt) : r;
            gate_idx_t pc = (c & incr_ctrl) ? (c ^ incr_tgt) : c;

            /* Fixed point: nothing to do */
            if (pr == r && pc == c) continue;

            /* Canonicalize destination to lower triangle */
            int dst_conj = (pr < pc);  /* Need conjugation if dest was upper-tri */
            gate_idx_t dr = dst_conj ? pc : pr;
            gate_idx_t dc = dst_conj ? pr : pc;

            /* Compute packed indices */
            gate_idx_t src_idx = pack_idx(dim, r, c);
            gate_idx_t dst_idx = pack_idx(dim, dr, dc);

            /* Self-conjugation: element stays in place but crosses the diagonal.
             * This happens when (r,c) maps to (pr,pc) with pr < pc, and after
             * canonicalization (dr,dc) == (r,c). The stored value must be conjugated. */
            if (src_idx == dst_idx) {
                if (dst_conj) data[src_idx] = conj(data[src_idx]);
                continue;
            }

            /* Only swap once per 2-cycle: process when src < dst */
            if (src_idx > dst_idx) continue;

            /* Swap, applying conjugation if the permutation crosses the diagonal.
             *
             * new_rho[r,c] = old_rho[pr,pc]. In packed storage:
             *   - src stores rho[r,c] directly (r >= c, lower triangle)
             *   - dst stores rho[dr,dc] directly (dr >= dc, lower triangle)
             *
             * If dst_conj: the logical source for (r,c) is rho[pr,pc] where pr<pc,
             * which is conj(rho[pc,pr]) = conj(data[dst_idx]). And the logical
             * source for (dr,dc)=(pc,pr) is rho[r,c] = data[src_idx], but we need
             * to store it in lower-tri form, which is conj(rho[r,c]) since the
             * destination represents the transposed position.
             */
            cplx_t src_val = data[src_idx];
            cplx_t dst_val = data[dst_idx];
            data[src_idx] = dst_conj ? conj(dst_val) : dst_val;
            data[dst_idx] = dst_conj ? conj(src_val) : src_val;
        }
    }
}





void cy_packed(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "cy: packed not yet implemented");
}

void cz_packed(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "cz: packed not yet implemented");
}

void cs_packed(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "cs: packed not yet implemented");
}

void csdg_packed(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "csdg: packed not yet implemented");
}

void ch_packed(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "ch: packed not yet implemented");
}

void chy_packed(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "chy: packed not yet implemented");
}

void ct_packed(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "ct: packed not yet implemented");
}

void ctdg_packed(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "ctdg: packed not yet implemented");
}

void cp_packed(state_t *state, const qubit_t control, const qubit_t target, const double theta) {
    GATE_VALIDATE(0, "cp: packed not yet implemented");
}

void cpdg_packed(state_t *state, const qubit_t control, const qubit_t target, const double theta) {
    GATE_VALIDATE(0, "cpdg: packed not yet implemented");
}

/*
 * In-place SWAP for packed storage.
 *
 * SWAP permutes basis states by exchanging bits `lo` and `hi`:
 *   perm(x) = x with bit lo and bit hi swapped
 * For density matrices: new_rho[r,c] = old_rho[perm(r), perm(c)].
 *
 * Like CNOT, SWAP is a zero-flop permutation gate — no arithmetic, only
 * element moves. The permutation maps each 4×4 block (indexed by the two-bit
 * pattern at positions lo,hi) into swap pairs:
 *
 *   Pair 1: (r01, c00) <-> (r10, c00)   always lower triangle
 *   Pair 2: (r01, c01) <-> (r10, c10)   always lower triangle
 *   Pair 3: (r11, c01) <-> (r11, c10)   always lower triangle
 *   Pair 4: (r10, c01) <-> (r01, c10)   may cross diagonal
 *   Pair 5: (r00, c01) <-> (r00, c10)   may cross diagonal
 *   Pair 6: (r10, c11) <-> (r01, c11)   may cross diagonal
 *
 * We iterate only the lower-triangular half of the block grid (br >= bc).
 * When a swap partner falls in the upper triangle (r < c), its packed
 * representative is at the transposed position with a conjugate. These
 * elements belong to block (bc, br) which is never iterated, so both
 * directions of the swap must be written in the current iteration.
 */
void swap_packed(state_t *state, const qubit_t q1, const qubit_t q2) {
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    const qubit_t hi = q1 > q2 ? q1 : q2;
    const qubit_t lo = q1 < q2 ? q1 : q2;
    const gate_idx_t incr_hi = (gate_idx_t)1 << hi;
    const gate_idx_t incr_lo = (gate_idx_t)1 << lo;
    cplx_t * restrict data = state->data;
    const gate_idx_t quarter_dim = dim >> 2;

    #pragma omp parallel for schedule(guided) if(dim >= OMP_THRESHOLD)
    for (gate_idx_t bc = 0; bc < quarter_dim; ++bc) {
        const gate_idx_t c00 = insertBits2_0(bc, lo, hi);
        const gate_idx_t c01 = c00 | incr_lo, c10 = c00 | incr_hi, c11 = c10 | incr_lo;

        /* Precompute offsets of columns */
        gate_idx_t offset_c00 = c00 * (2 * dim - c00 + 1) / 2;
        gate_idx_t offset_c01 = c01 * (2 * dim - c01 + 1) / 2;
        gate_idx_t offset_c10 = c10 * (2 * dim - c10 + 1) / 2;
        gate_idx_t offset_c11 = c11 * (2 * dim - c11 + 1) / 2;


        for (gate_idx_t br = bc; br < quarter_dim; ++br) {
            const gate_idx_t r00 = insertBits2_0(br, lo, hi);
            const gate_idx_t r01 = r00 | incr_lo, r10 = r00 | incr_hi, r11 = r10 | incr_lo;

            /* Pair 1: (r01, c00) <-> (r10, c00) — both always in lower triangle */
            gate_idx_t idx_a = (r01 - c00) + offset_c00;
            gate_idx_t idx_b = (r10 - c00) + offset_c00;
            cplx_t tmp = data[idx_a];
            data[idx_a] = data[idx_b];
            data[idx_b] = tmp;

            /* Pair 2: (r01, c01) <-> (r10, c10) — both always in lower triangle */
            idx_a = (r01 - c01) + offset_c01;
            idx_b = (r10 - c10) + offset_c10;
            tmp = data[idx_a];
            data[idx_a] = data[idx_b];
            data[idx_b] = tmp;

            /* Pair 3: (r11, c01) <-> (r11, c10) — both always in lower triangle */
            idx_a = (r11 - c01) + offset_c01;
            idx_b = (r11 - c10) + offset_c10;
            tmp = data[idx_a];
            data[idx_a] = data[idx_b];
            data[idx_b] = tmp;

            /* Pair 4: (r10, c01) <-> (r01, c10)
             * (r10, c01) is always lower-tri. (r01, c10) may cross the diagonal:
             *   r01 > c10: both lower-tri, plain swap.
             *   r01 <= c10: (r01, c10) is upper-tri; its packed representative is
             *     conj(data) at (c10, r01). That position belongs to block (bc, br),
             *     never iterated, so we must write both directions here. */
            idx_a = (r10 - c01) + offset_c01;
            tmp = data[idx_a];
            if (r01 > c10) {
                idx_b = (r01 - c10) + offset_c10;
                data[idx_a] = data[idx_b];
                data[idx_b] = tmp;
            }
            else {
                idx_b = pack_idx(dim, c10, r01);
                data[idx_a] = conj(data[idx_b]);
                data[idx_b] = conj(tmp);
            }

            /* Pairs 5+6: (r00, c01) <-> (r00, c10) and (r10, c11) <-> (r01, c11)
             *
             * Cases based on how r00 compares to c01 and c10:
             *
             * (a) r00 > c10 > c01: all elements in lower triangle, plain swaps.
             * (b) r00 > c01 but r00 <= c10: (r00, c01) is lower-tri but (r00, c10)
             *     is upper-tri. Swap via conjugated packed representatives at
             *     (c10, r00) and (c11, r01). Same logic as Pair 4.
             * (c) r00 <= c01 (only for br > bc, starting at 4+ qubits): both
             *     elements of each pair are upper-tri. Their packed representatives
             *     (c01, r00) <-> (c10, r00) and (c11, r10) <-> (c11, r01) are in
             *     the lower triangle but in unreachable block (bc, br). Swap directly.
             */
            if (r00 > c01) {
                idx_a = (r00 - c01) + offset_c01;
                const gate_idx_t idx_aa = (r10 - c11) + offset_c11;

                if (r00 > c10) {
                    /* (a) Both lower-tri: plain swap */
                    idx_b = (r00 - c10) + offset_c10;
                    tmp = data[idx_a];
                    data[idx_a] = data[idx_b];
                    data[idx_b] = tmp;

                    idx_b = (r01 - c11) + offset_c11;
                    tmp = data[idx_aa];
                    data[idx_aa] = data[idx_b];
                    data[idx_b] = tmp;
                }
                else {
                    /* (b) Cross-diagonal: conjugated swap via (c10,r00) and (c11,r01) */
                    idx_b = pack_idx(dim, c10, r00);
                    tmp = data[idx_a];
                    data[idx_a] = conj(data[idx_b]);
                    data[idx_b] = conj(tmp);

                    idx_b = pack_idx(dim, c11, r01);
                    tmp = data[idx_aa];
                    data[idx_aa] = conj(data[idx_b]);
                    data[idx_b] = conj(tmp);
                }
            }
            else if (bc != br) {
                /* (c) Both upper-tri: swap packed representatives directly */
                idx_a = pack_idx(dim, c01, r00);
                idx_b = pack_idx(dim, c10, r00);
                tmp = data[idx_a];
                data[idx_a] = data[idx_b];
                data[idx_b] = tmp;

                idx_a = pack_idx(dim, c11, r10);
                idx_b = pack_idx(dim, c11, r01);
                tmp = data[idx_a];
                data[idx_a] = data[idx_b];
                data[idx_b] = tmp;
            }
        }
    }
}

void ccx_packed(state_t *state, const qubit_t ctrl1, const qubit_t ctrl2, const qubit_t target) {
    GATE_VALIDATE(0, "ccx: packed not yet implemented");
}
