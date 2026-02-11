/**
 * @file gate_tiled.c
 * @brief Quantum gates for mixed states (density matrices in tiled lower-triangular storage)
 *
 * Tiled storage partitions the density matrix into TILE_DIM x TILE_DIM tiles.
 * Only lower-triangular tiles are stored (at tile level). Within each tile,
 * elements are stored row-major.
 *
 * Single-qubit gates: ρ → UρU†. Two traversal modes:
 *
 *   Within-tile (target < LOG_TILE_DIM):
 *     The target bit affects only local coordinates within a tile. All 4 butterfly
 *     elements lie within the same tile. No Hermitian conjugation branches needed.
 *
 *   Cross-tile (target >= LOG_TILE_DIM):
 *     The target bit affects tile coordinates. The 4 elements span 4 tiles:
 *     T(tr0,tc0), T(tr0,tc1), T(tr1,tc0), T(tr1,tc1). Only lower-triangular
 *     tiles are stored, so upper-triangle tiles require conjugation.
 */

#include "gate.h"
#include <complex.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Minimum dimension to enable OpenMP parallelization */
#define OMP_THRESHOLD 64

/*
 * =====================================================================================================================
 * Tiled index helpers
 * =====================================================================================================================
 */

/* insertBit0() is defined in gate.h as a static inline shared across all gate backends. */

static inline gate_idx_t tile_off(gate_idx_t tr, gate_idx_t tc) {
    return (tr * (tr + 1) / 2 + tc) * TILE_SIZE;
}

static inline gate_idx_t elem_off(gate_idx_t lr, gate_idx_t lc) {
    return lr * TILE_DIM + lc;
}

/*
 * =====================================================================================================================
 * Within-tile traversal macro (target < LOG_TILE_DIM)
 * =====================================================================================================================
 *
 * Iterates over every stored tile, and within each tile processes all 2×2 blocks
 * that mix under the target qubit. Since both row and column indices only differ
 * in a local (within-tile) bit, all 4 elements of each block are in the same tile.
 *
 * For each block:
 *   - a00 = &tile[elem_off(lr0, lc0)]  -- local row/col both have target bit 0
 *   - a01 = &tile[elem_off(lr0, lc1)]  -- local col has target bit 1
 *   - a10 = &tile[elem_off(lr1, lc0)]  -- local row has target bit 1
 *   - a11 = &tile[elem_off(lr1, lc1)]  -- both have target bit 1
 *
 * Tiles store full TILE_DIM × TILE_DIM elements (including diagonal tiles), so
 * all 4 elements are directly accessible with no conjugation logic needed.
 */
#define TRAVERSE_TILED_WITHIN(state, target, BLOCK_OP)                          \
do {                                                                            \
    const gate_idx_t dim = (gate_idx_t)1 << (state)->qubits;                    \
    const gate_idx_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;            \
    const gate_idx_t stride = (gate_idx_t)1 << (target);                        \
    const gate_idx_t tile_dim = (dim < TILE_DIM) ? dim : TILE_DIM;              \
    const gate_idx_t half_td = tile_dim >> 1;                                   \
    cplx_t *data = (state)->data;                                               \
                                                                                \
    _Pragma("omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)")      \
    for (gate_idx_t tr = 0; tr < n_tiles; ++tr) {                               \
        for (gate_idx_t tc = 0; tc <= tr; ++tc) {                               \
            cplx_t *tile = data + tile_off(tr, tc);                             \
                                                                                \
            for (gate_idx_t br = 0; br < half_td; ++br) {                       \
                gate_idx_t lr0 = insertBit0(br, (target));                      \
                gate_idx_t lr1 = lr0 | stride;                                  \
                for (gate_idx_t bc = 0; bc < half_td; ++bc) {                   \
                    gate_idx_t lc0 = insertBit0(bc, (target));                  \
                    gate_idx_t lc1 = lc0 | stride;                              \
                    cplx_t *a00 = tile + elem_off(lr0, lc0);                   \
                    cplx_t *a01 = tile + elem_off(lr0, lc1);                   \
                    cplx_t *a10 = tile + elem_off(lr1, lc0);                   \
                    cplx_t *a11 = tile + elem_off(lr1, lc1);                   \
                    BLOCK_OP(a00, a01, a10, a11);                               \
                }                                                               \
            }                                                                   \
        }                                                                       \
    }                                                                           \
} while(0)

/*
 * =====================================================================================================================
 * Cross-tile traversal macro (target >= LOG_TILE_DIM)
 * =====================================================================================================================
 *
 * The target bit affects tile coordinates. For tile-row index tr, bits above
 * LOG_TILE_DIM correspond to the tile-row, so the target bit (relative to tile
 * indexing) is at position (target - LOG_TILE_DIM).
 *
 * We enumerate tile-row pairs (tr0, tr1) where tr1 = tr0 | (1 << tile_bit),
 * and tile-col pairs (tc0, tc1) similarly. The 4 blocks are:
 *   T(tr0,tc0), T(tr0,tc1), T(tr1,tc0), T(tr1,tc1)
 *
 * Since only lower-triangular tiles are stored:
 *   - T(tr0,tc0): stored if tr0 >= tc0 (always, by construction)
 *   - T(tr1,tc1): stored if tr1 >= tc1 (always, since tr1>tr0, tc1>tc0)
 *   - T(tr1,tc0): stored if tr1 >= tc0 (always, since tr1>tr0>=tc0)
 *   - T(tr0,tc1): stored only if tr0 >= tc1. If not, use conj of T(tc1,tr0).
 *
 * Within each tile-block, the local coordinates are identical for all 4 tiles.
 * We iterate over all (lr, lc) within the tile; for diagonal tile-blocks
 * (tc0 == tr0 at tile level), only the lower triangle needs processing.
 */
/*
 * Cross-tile traversal for diagonal tiles:
 *
 * When is_diag (i_tc == i_tr, so tr0 == tc0 and tr1 == tc1), tiles t00 and t11
 * are diagonal tiles that store only the lower triangle (lr >= lc). To avoid
 * reading uninitialized upper-triangle data, we iterate lc from 0 to lr
 * (inclusive), ensuring lr >= lc. The BLOCK_OP then processes the lower-triangle
 * stored values directly.
 *
 * For t10 = T(tr1, tc0): off-diagonal tile, all elements stored -- lr,lc fine.
 * For t01: when is_diag, tr0 == tc0 < tc1, so lower_01 = 0, meaning
 *          t01 = T(tc1, tr0) which is off-diagonal -- elem at (lc, lr).
 *          But with the flipped iteration (lr >= lc), t01 elem is at (lc, lr).
 */
#define TRAVERSE_TILED_CROSS(state, target, BLOCK_OP)                           \
do {                                                                            \
    const gate_idx_t dim = (gate_idx_t)1 << (state)->qubits;                    \
    const gate_idx_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;            \
    const qubit_t tile_bit = (target) - LOG_TILE_DIM;                           \
    const gate_idx_t tile_stride = (gate_idx_t)1 << tile_bit;                   \
    const gate_idx_t half_nt = n_tiles >> 1;                                    \
    cplx_t *data = (state)->data;                                               \
                                                                                \
    /* Iterate tile-row pairs (tr0, tr1) with tile_bit=0 */                     \
    _Pragma("omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)")      \
    for (gate_idx_t i_tr = 0; i_tr < half_nt; ++i_tr) {                        \
        gate_idx_t tr0 = insertBit0(i_tr, tile_bit);                            \
        gate_idx_t tr1 = tr0 | tile_stride;                                     \
                                                                                \
        /* Iterate tile-column pairs (tc0, tc1) with tile_bit=0 */             \
        /* Up to i_tr so that tr0 >= tc0 */                                    \
        for (gate_idx_t i_tc = 0; i_tc <= i_tr; ++i_tc) {                      \
            gate_idx_t tc0 = insertBit0(i_tc, tile_bit);                        \
            gate_idx_t tc1 = tc0 | tile_stride;                                 \
                                                                                \
            /* Tile pointers: T00, T11, T10 always in lower triangle */        \
            cplx_t *t00 = data + tile_off(tr0, tc0);                           \
            cplx_t *t11 = data + tile_off(tr1, tc1);                           \
            cplx_t *t10 = data + tile_off(tr1, tc0);                           \
                                                                                \
            /* T01: T(tr0, tc1) - only stored if tr0 >= tc1 */                 \
            int lower_01 = (tr0 >= tc1);                                       \
            cplx_t *t01 = lower_01 ? data + tile_off(tr0, tc1)                 \
                                   : data + tile_off(tc1, tr0);                \
                                                                                \
            int is_diag = (i_tc == i_tr);                                      \
                                                                                \
            /* Iterate all (lr, lc) within tiles.                              \
             * For is_diag: iterate lr >= lc (lower triangle) because t00     \
             * and t11 are diagonal tiles that only store lower-triangle       \
             * elements. For !is_diag: iterate all elements. */                \
            for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {                     \
                gate_idx_t lc_end = is_diag ? lr + 1 : TILE_DIM;              \
                for (gate_idx_t lc = 0; lc < lc_end; ++lc) {                  \
                    cplx_t *a00 = t00 + elem_off(lr, lc);                      \
                    cplx_t *a11 = t11 + elem_off(lr, lc);                      \
                    cplx_t *a10 = t10 + elem_off(lr, lc);                      \
                    /* For t01: if stored as T(tr0,tc1), elem is (lr,lc)      \
                     * If stored as T(tc1,tr0), elem is at (lc,lr) */          \
                    cplx_t *a01 = lower_01 ? t01 + elem_off(lr, lc)            \
                                           : t01 + elem_off(lc, lr);           \
                    int diag_block = is_diag && (lr == lc);                     \
                                                                                \
                    BLOCK_OP(a00, a01, a10, a11, lower_01, diag_block);        \
                    /* Mirror upper triangle of diagonal tiles t00, t11 */      \
                    if (is_diag && lr != lc) {                                  \
                        *(t00 + elem_off(lc, lr)) = conj(*a00);                \
                        *(t11 + elem_off(lc, lr)) = conj(*a11);                \
                    }                                                           \
                }                                                               \
            }                                                                   \
        }                                                                       \
    }                                                                           \
} while(0)

/*
 * =====================================================================================================================
 * Block operation macros: Within-tile (direct access, no conjugation)
 * =====================================================================================================================
 *
 * These macros receive pointers to the 4 elements of a 2×2 block.
 * Since tiles store full TILE_DIM × TILE_DIM data, all 4 elements are
 * directly accessible with no conjugation logic needed.
 */

/*
 * Off-diagonal tile ops: full tile stored, all 4 elements directly accessible.
 * No conjugation needed. Only 4 args: a00, a01, a10, a11.
 */
#define X_OFFDIAG_OP(a00, a01, a10, a11) do { \
    cplx_t tmp = *(a00);                       \
    *(a00) = *(a11);                            \
    *(a11) = tmp;                               \
    tmp = *(a01);                               \
    *(a01) = *(a10);                            \
    *(a10) = tmp;                               \
} while(0)

#define Y_OFFDIAG_OP(a00, a01, a10, a11) do { \
    cplx_t tmp = *(a00);                       \
    *(a00) = *(a11);                            \
    *(a11) = tmp;                               \
    tmp = *(a01);                               \
    *(a01) = -*(a10);                           \
    *(a10) = -tmp;                              \
} while(0)

#define Z_OFFDIAG_OP(a00, a01, a10, a11) do { \
    (void)(a00); (void)(a11);                  \
    *(a01) = -*(a01);                          \
    *(a10) = -*(a10);                          \
} while(0)

#define H_OFFDIAG_OP(a00, a01, a10, a11) do { \
    cplx_t v00 = *(a00), v01 = *(a01);        \
    cplx_t v10 = *(a10), v11 = *(a11);        \
    *(a00) = 0.5 * (v00 + v01 + v10 + v11);   \
    *(a01) = 0.5 * (v00 - v01 + v10 - v11);   \
    *(a10) = 0.5 * (v00 + v01 - v10 - v11);   \
    *(a11) = 0.5 * (v00 - v01 - v10 + v11);   \
} while(0)

#define S_OFFDIAG_OP(a00, a01, a10, a11) do { \
    (void)(a00); (void)(a11);                        \
    *(a01) = CMPLX(cimag(*(a01)), -creal(*(a01)));  \
    *(a10) = CMPLX(-cimag(*(a10)), creal(*(a10)));  \
} while(0)

#define SDG_OFFDIAG_OP(a00, a01, a10, a11) do {\
    (void)(a00); (void)(a11);                        \
    *(a01) = CMPLX(-cimag(*(a01)), creal(*(a01)));  \
    *(a10) = CMPLX(cimag(*(a10)), -creal(*(a10)));  \
} while(0)

#define T_OFFDIAG_OP(a00, a01, a10, a11) do { \
    (void)(a00); (void)(a11);                         \
    double r01 = creal(*(a01)), i01 = cimag(*(a01)); \
    double r10 = creal(*(a10)), i10 = cimag(*(a10)); \
    *(a01) = CMPLX((r01 + i01) * M_SQRT1_2,          \
                   (i01 - r01) * M_SQRT1_2);          \
    *(a10) = CMPLX((r10 - i10) * M_SQRT1_2,          \
                   (r10 + i10) * M_SQRT1_2);          \
} while(0)

#define TDG_OFFDIAG_OP(a00, a01, a10, a11) do {\
    (void)(a00); (void)(a11);                         \
    double r01 = creal(*(a01)), i01 = cimag(*(a01)); \
    double r10 = creal(*(a10)), i10 = cimag(*(a10)); \
    *(a01) = CMPLX((r01 - i01) * M_SQRT1_2,          \
                   (r01 + i01) * M_SQRT1_2);          \
    *(a10) = CMPLX((r10 + i10) * M_SQRT1_2,          \
                   (i10 - r10) * M_SQRT1_2);          \
} while(0)

#define HY_OFFDIAG_OP(a00, a01, a10, a11) do { \
    cplx_t v00 = *(a00), v01 = *(a01);          \
    cplx_t v10 = *(a10), v11 = *(a11);          \
    *(a00) = 0.5 * (v00 - v01 - v10 + v11);     \
    *(a01) = 0.5 * (v00 + v01 - v10 - v11);     \
    *(a10) = 0.5 * (v00 - v01 + v10 - v11);     \
    *(a11) = 0.5 * (v00 + v01 + v10 + v11);     \
} while(0)

/*
 * =====================================================================================================================
 * Block operation macros: Cross-tile (with conjugation handling)
 * =====================================================================================================================
 *
 * Cross-tile macros receive pointers and flags:
 *   lower_01: if 1, a01 points to T(tr0,tc1) directly; if 0, it points to
 *             T(tc1,tr0) and the stored value is the conjugate of ρ(r0,c1).
 *   diag_block: if 1, the (0,1) and (1,0) elements are a conjugate pair.
 *
 * The pattern mirrors the packed BLOCK_OP macros but reads/writes through
 * tile pointers instead of packed indices.
 */

#define X_CROSS_OP(a00, a01, a10, a11, lower_01, diag_block) do { \
    cplx_t tmp = *a00; *a00 = *a11; *a11 = tmp;                    \
    if (diag_block) {                                               \
        /* a10 is stored directly, a01 = conj(a10) */              \
        *a10 = conj(*a10);                                         \
        *a01 = lower_01 ? conj(*a10) : *a10;                      \
    } else {                                                        \
        cplx_t rho01 = lower_01 ? *a01 : conj(*a01);               \
        cplx_t rho10 = *a10;                                        \
        *a10 = rho01;                                               \
        *a01 = lower_01 ? rho10 : conj(rho10);                     \
    }                                                               \
} while(0)

#define Y_CROSS_OP(a00, a01, a10, a11, lower_01, diag_block) do { \
    cplx_t tmp = *a00; *a00 = *a11; *a11 = tmp;                    \
    if (diag_block) {                                               \
        *a10 = -conj(*a10);                                        \
        *a01 = lower_01 ? conj(*a10) : *a10;                      \
    } else {                                                        \
        cplx_t rho01 = lower_01 ? *a01 : conj(*a01);               \
        cplx_t rho10 = *a10;                                        \
        *a10 = -rho01;                                              \
        *a01 = lower_01 ? -rho10 : conj(-rho10);                   \
    }                                                               \
} while(0)

#define Z_CROSS_OP(a00, a01, a10, a11, lower_01, diag_block) do { \
    (void)(a00); (void)(a11);                                       \
    *a10 = -*a10;                                                   \
    if (!diag_block) *a01 = -*a01;                                  \
} while(0)

#define H_CROSS_OP(a00, a01, a10, a11, lower_01, diag_block) do {  \
    cplx_t rho01 = lower_01 ? *a01 : conj(*a01);                    \
    cplx_t v00 = *a00, v01 = rho01, v10 = *a10, v11 = *a11;        \
    *a00 = 0.5 * (v00 + v01 + v10 + v11);                           \
    *a11 = 0.5 * (v00 - v01 - v10 + v11);                           \
    cplx_t new10 = 0.5 * (v00 + v01 - v10 - v11);                   \
    *a10 = new10;                                                    \
    if (diag_block) {                                                \
        /* On diagonal blocks, a01 = conj(a10) by Hermiticity */    \
        *a01 = lower_01 ? conj(new10) : new10;                      \
    } else {                                                         \
        cplx_t new01 = 0.5 * (v00 - v01 + v10 - v11);              \
        *a01 = lower_01 ? new01 : conj(new01);                      \
    }                                                                \
} while(0)

#define S_CROSS_OP(a00, a01, a10, a11, lower_01, diag_block) do {  \
    (void)(a00); (void)(a11);                                        \
    /* rho10 *= i, rho01 *= -i */                                   \
    *a10 = CMPLX(-cimag(*a10), creal(*a10));                        \
    if (!diag_block) {                                               \
        if (lower_01) {                                              \
            *a01 = CMPLX(cimag(*a01), -creal(*a01));                \
        } else {                                                     \
            /* stored = conj(rho01); new = conj(rho01 * (-i))     */\
            /* = stored * i */                                      \
            *a01 = CMPLX(-cimag(*a01), creal(*a01));                \
        }                                                            \
    }                                                                \
} while(0)

#define SDG_CROSS_OP(a00, a01, a10, a11, lower_01, diag_block) do {\
    (void)(a00); (void)(a11);                                        \
    /* rho10 *= -i, rho01 *= i */                                   \
    *a10 = CMPLX(cimag(*a10), -creal(*a10));                        \
    if (!diag_block) {                                               \
        if (lower_01) {                                              \
            *a01 = CMPLX(-cimag(*a01), creal(*a01));                \
        } else {                                                     \
            /* stored = conj(rho01); new = conj(rho01 * i)        */\
            /* = stored * (-i) */                                   \
            *a01 = CMPLX(cimag(*a01), -creal(*a01));                \
        }                                                            \
    }                                                                \
} while(0)

#define T_CROSS_OP(a00, a01, a10, a11, lower_01, diag_block) do {  \
    (void)(a00); (void)(a11);                                        \
    double r10 = creal(*a10), i10 = cimag(*a10);                    \
    /* rho10 *= e^(iπ/4) = (1+i)/√2 */                             \
    *a10 = CMPLX((r10 - i10) * M_SQRT1_2,                           \
                 (r10 + i10) * M_SQRT1_2);                           \
    if (!diag_block) {                                               \
        if (lower_01) {                                              \
            double r01 = creal(*a01), i01 = cimag(*a01);            \
            /* rho01 *= e^(-iπ/4) = (1-i)/√2 */                    \
            *a01 = CMPLX((r01 + i01) * M_SQRT1_2,                   \
                         (i01 - r01) * M_SQRT1_2);                   \
        } else {                                                     \
            double r = creal(*a01), im = cimag(*a01);               \
            /* stored = conj(rho01); new = stored * (1+i)/√2 */    \
            *a01 = CMPLX((r - im) * M_SQRT1_2,                      \
                         (r + im) * M_SQRT1_2);                      \
        }                                                            \
    }                                                                \
} while(0)

#define TDG_CROSS_OP(a00, a01, a10, a11, lower_01, diag_block) do {\
    (void)(a00); (void)(a11);                                        \
    double r10 = creal(*a10), i10 = cimag(*a10);                    \
    /* rho10 *= e^(-iπ/4) = (1-i)/√2 */                            \
    *a10 = CMPLX((r10 + i10) * M_SQRT1_2,                           \
                 (i10 - r10) * M_SQRT1_2);                           \
    if (!diag_block) {                                               \
        if (lower_01) {                                              \
            double r01 = creal(*a01), i01 = cimag(*a01);            \
            /* rho01 *= e^(iπ/4) = (1+i)/√2 */                     \
            *a01 = CMPLX((r01 - i01) * M_SQRT1_2,                   \
                         (r01 + i01) * M_SQRT1_2);                   \
        } else {                                                     \
            double r = creal(*a01), im = cimag(*a01);               \
            /* stored = conj(rho01); new = stored * (1-i)/√2 */    \
            *a01 = CMPLX((r + im) * M_SQRT1_2,                      \
                         (im - r) * M_SQRT1_2);                      \
        }                                                            \
    }                                                                \
} while(0)

#define HY_CROSS_OP(a00, a01, a10, a11, lower_01, diag_block) do { \
    cplx_t rho01 = lower_01 ? *a01 : conj(*a01);                    \
    cplx_t v00 = *a00, v01 = rho01, v10 = *a10, v11 = *a11;        \
    *a00 = 0.5 * (v00 - v01 - v10 + v11);                           \
    *a11 = 0.5 * (v00 + v01 + v10 + v11);                           \
    cplx_t new10 = 0.5 * (v00 - v01 + v10 - v11);                   \
    *a10 = new10;                                                    \
    if (diag_block) {                                                \
        /* On diagonal blocks, a01 = conj(a10) by Hermiticity */    \
        *a01 = lower_01 ? conj(new10) : new10;                      \
    } else {                                                         \
        cplx_t new01 = 0.5 * (v00 + v01 - v10 - v11);              \
        *a01 = lower_01 ? new01 : conj(new01);                      \
    }                                                                \
} while(0)

/*
 * =====================================================================================================================
 * Dispatch: within-tile vs cross-tile
 * =====================================================================================================================
 */

#define DISPATCH_TILED_1Q(state, target, WITHIN_OP, CROSS_OP)      \
do {                                                                \
    if ((target) < LOG_TILE_DIM) {                                  \
        TRAVERSE_TILED_WITHIN(state, target, WITHIN_OP);            \
    } else {                                                        \
        TRAVERSE_TILED_CROSS(state, target, CROSS_OP);              \
    }                                                               \
} while(0)

/*
 * =====================================================================================================================
 * Pauli gates
 * =====================================================================================================================
 */

void x_tiled(state_t *state, const qubit_t target) {
    DISPATCH_TILED_1Q(state, target, X_OFFDIAG_OP, X_CROSS_OP);
}

void y_tiled(state_t *state, const qubit_t target) {
    DISPATCH_TILED_1Q(state, target, Y_OFFDIAG_OP, Y_CROSS_OP);
}

void z_tiled(state_t *state, const qubit_t target) {
    DISPATCH_TILED_1Q(state, target, Z_OFFDIAG_OP, Z_CROSS_OP);
}

/*
 * =====================================================================================================================
 * Clifford gates
 * =====================================================================================================================
 */

void h_tiled(state_t *state, const qubit_t target) {
    DISPATCH_TILED_1Q(state, target, H_OFFDIAG_OP, H_CROSS_OP);
}

void s_tiled(state_t *state, const qubit_t target) {
    DISPATCH_TILED_1Q(state, target, S_OFFDIAG_OP, S_CROSS_OP);
}

void sdg_tiled(state_t *state, const qubit_t target) {
    DISPATCH_TILED_1Q(state, target, SDG_OFFDIAG_OP, SDG_CROSS_OP);
}

void t_tiled(state_t *state, const qubit_t target) {
    DISPATCH_TILED_1Q(state, target, T_OFFDIAG_OP, T_CROSS_OP);
}

void tdg_tiled(state_t *state, const qubit_t target) {
    DISPATCH_TILED_1Q(state, target, TDG_OFFDIAG_OP, TDG_CROSS_OP);
}

void hy_tiled(state_t *state, const qubit_t target) {
    DISPATCH_TILED_1Q(state, target, HY_OFFDIAG_OP, HY_CROSS_OP);
}

/*
 * =====================================================================================================================
 * Rotation gates - Traversal macros to eliminate boilerplate
 * =====================================================================================================================
 *
 * Mixing rotation gates (Rx, Ry) share identical traversal structure: all 4 quadrants
 * of each 2x2 block mix together. Only the computation formula differs. These macros
 * factor out the common traversal, reading, and write-back logic.
 *
 * COMPUTE_OP(v00, v01, v10, v11, n00, n01, n10, n11): given input values v00..v11,
 * compute and assign output values n00..n11. Captures trig constants from enclosing scope.
 */

/*
 * Within-tile traversal for mixing rotation gates.
 * target < LOG_TILE_DIM: all 4 butterfly elements within the same tile.
 * Tiles store full TILE_DIM × TILE_DIM data, so all elements are directly accessible.
 */
#define TRAVERSE_TILED_WITHIN_ROT(state, target, COMPUTE_OP)                    \
do {                                                                            \
    const gate_idx_t dim = (gate_idx_t)1 << (state)->qubits;                    \
    const gate_idx_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;            \
    const gate_idx_t stride = (gate_idx_t)1 << (target);                        \
    const gate_idx_t tile_dim = (dim < TILE_DIM) ? dim : TILE_DIM;              \
    const gate_idx_t half_td = tile_dim >> 1;                                   \
    cplx_t *data = (state)->data;                                               \
                                                                                \
    _Pragma("omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)")      \
    for (gate_idx_t tr = 0; tr < n_tiles; ++tr) {                               \
        for (gate_idx_t tc = 0; tc <= tr; ++tc) {                               \
            cplx_t *tile = data + tile_off(tr, tc);                             \
                                                                                \
            for (gate_idx_t br = 0; br < half_td; ++br) {                       \
                gate_idx_t lr0 = insertBit0(br, (target));                      \
                gate_idx_t lr1 = lr0 | stride;                                  \
                for (gate_idx_t bc = 0; bc < half_td; ++bc) {                   \
                    gate_idx_t lc0 = insertBit0(bc, (target));                  \
                    gate_idx_t lc1 = lc0 | stride;                              \
                                                                                \
                    cplx_t v00 = tile[elem_off(lr0, lc0)];                      \
                    cplx_t v01 = tile[elem_off(lr0, lc1)];                      \
                    cplx_t v10 = tile[elem_off(lr1, lc0)];                      \
                    cplx_t v11 = tile[elem_off(lr1, lc1)];                      \
                                                                                \
                    cplx_t n00, n01, n10, n11;                                  \
                    COMPUTE_OP(v00, v01, v10, v11, n00, n01, n10, n11);         \
                                                                                \
                    tile[elem_off(lr0, lc0)] = n00;                             \
                    tile[elem_off(lr0, lc1)] = n01;                             \
                    tile[elem_off(lr1, lc0)] = n10;                             \
                    tile[elem_off(lr1, lc1)] = n11;                             \
                }                                                               \
            }                                                                   \
        }                                                                       \
    }                                                                           \
} while(0)

/*
 * Cross-tile traversal for mixing rotation gates.
 * target >= LOG_TILE_DIM: the 4 elements span 4 tiles.
 *
 * NOTE: When is_diag (i_tc == i_tr), tiles t00 and t11 are diagonal tiles
 * that only store the lower triangle (lr >= lc). We iterate lc from 0 to lr
 * (inclusive) to ensure we only read valid lower-triangle positions.
 */
#define TRAVERSE_TILED_CROSS_ROT(state, target, COMPUTE_OP)                     \
do {                                                                            \
    const gate_idx_t dim = (gate_idx_t)1 << (state)->qubits;                    \
    const gate_idx_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;            \
    const qubit_t tile_bit = (target) - LOG_TILE_DIM;                           \
    const gate_idx_t tile_stride = (gate_idx_t)1 << tile_bit;                   \
    const gate_idx_t half_nt = n_tiles >> 1;                                    \
    cplx_t *data = (state)->data;                                               \
                                                                                \
    _Pragma("omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)")      \
    for (gate_idx_t i_tr = 0; i_tr < half_nt; ++i_tr) {                        \
        gate_idx_t tr0 = insertBit0(i_tr, tile_bit);                            \
        gate_idx_t tr1 = tr0 | tile_stride;                                     \
        for (gate_idx_t i_tc = 0; i_tc <= i_tr; ++i_tc) {                      \
            gate_idx_t tc0 = insertBit0(i_tc, tile_bit);                        \
            gate_idx_t tc1 = tc0 | tile_stride;                                 \
                                                                                \
            cplx_t *t00 = data + tile_off(tr0, tc0);                           \
            cplx_t *t11 = data + tile_off(tr1, tc1);                           \
            cplx_t *t10 = data + tile_off(tr1, tc0);                           \
            int lower_01 = (tr0 >= tc1);                                       \
            cplx_t *t01 = lower_01 ? data + tile_off(tr0, tc1)                 \
                                   : data + tile_off(tc1, tr0);                \
            int is_diag = (i_tc == i_tr);                                      \
                                                                                \
            for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {                     \
                gate_idx_t lc_end = is_diag ? lr + 1 : TILE_DIM;              \
                for (gate_idx_t lc = 0; lc < lc_end; ++lc) {                  \
                    int diag_block = is_diag && (lr == lc);                     \
                                                                                \
                    cplx_t v00 = t00[elem_off(lr, lc)];                        \
                    cplx_t v11 = t11[elem_off(lr, lc)];                        \
                    cplx_t v10 = t10[elem_off(lr, lc)];                        \
                    cplx_t v01 = lower_01 ? t01[elem_off(lr, lc)]              \
                                          : conj(t01[elem_off(lc, lr)]);       \
                                                                                \
                    cplx_t n00, n01, n10, n11;                                  \
                    COMPUTE_OP(v00, v01, v10, v11, n00, n01, n10, n11);         \
                                                                                \
                    t00[elem_off(lr, lc)] = n00;                               \
                    t11[elem_off(lr, lc)] = n11;                               \
                    t10[elem_off(lr, lc)] = n10;                               \
                    if (diag_block) {                                            \
                        if (lower_01) t01[elem_off(lr, lc)] = conj(n10);       \
                        else          t01[elem_off(lc, lr)] = n10;             \
                    } else {                                                     \
                        if (lower_01) t01[elem_off(lr, lc)] = n01;             \
                        else          t01[elem_off(lc, lr)] = conj(n01);       \
                    }                                                            \
                    /* Mirror upper triangle of diagonal tiles t00, t11 */      \
                    if (is_diag && lr != lc) {                                  \
                        t00[elem_off(lc, lr)] = conj(n00);                     \
                        t11[elem_off(lc, lr)] = conj(n11);                     \
                    }                                                            \
                }                                                               \
            }                                                                   \
        }                                                                       \
    }                                                                           \
} while(0)

/*
 * =====================================================================================================================
 * Rx gate (tiled): ρ → Rx(θ) ρ Rx(θ)†
 * =====================================================================================================================
 *
 * Rx(θ) = [[c, -is], [-is, c]] where c=cos(θ/2), s=sin(θ/2)
 *
 *   ρ'_00 = c²ρ00 + s²ρ11 + ics(ρ01 - ρ10)
 *   ρ'_11 = s²ρ00 + c²ρ11 - ics(ρ01 - ρ10)
 *   ρ'_01 = ics(ρ00 - ρ11) + c²ρ01 + s²ρ10
 *   ρ'_10 = -ics(ρ00 - ρ11) + c²ρ10 + s²ρ01
 */
#define RX_COMPUTE(v00, v01, v10, v11, n00, n01, n10, n11) do {                \
    cplx_t _d01 = (v01) - (v10);                                               \
    cplx_t _ics = CMPLX(-cs * cimag(_d01), cs * creal(_d01));                  \
    cplx_t _dd  = (v00) - (v11);                                               \
    cplx_t _icd = CMPLX(-cs * cimag(_dd), cs * creal(_dd));                    \
    (n00) = c2 * (v00) + s2 * (v11) + _ics;                                    \
    (n11) = s2 * (v00) + c2 * (v11) - _ics;                                    \
    (n10) = -_icd + c2 * (v10) + s2 * (v01);                                   \
    (n01) =  _icd + c2 * (v01) + s2 * (v10);                                   \
} while(0)

void rx_tiled(state_t *state, const qubit_t target, const double theta) {
    const double co = cos(theta / 2.0);
    const double si = sin(theta / 2.0);
    const double c2 = co * co, s2 = si * si, cs = co * si;

    if (target < LOG_TILE_DIM) {
        TRAVERSE_TILED_WITHIN_ROT(state, target, RX_COMPUTE);
    } else {
        TRAVERSE_TILED_CROSS_ROT(state, target, RX_COMPUTE);
    }
}
#undef RX_COMPUTE

/*
 * =====================================================================================================================
 * Ry gate (tiled): ρ → Ry(θ) ρ Ry(θ)†
 * =====================================================================================================================
 *
 * Ry(θ) = [[c, -s], [s, c]] where c=cos(θ/2), s=sin(θ/2)
 *
 *   ρ'_00 = c²ρ00 + s²ρ11 - cs(ρ01 + ρ10)
 *   ρ'_11 = s²ρ00 + c²ρ11 + cs(ρ01 + ρ10)
 *   ρ'_01 = cs(ρ00 - ρ11) + c²ρ01 - s²ρ10
 *   ρ'_10 = cs(ρ00 - ρ11) + c²ρ10 - s²ρ01
 */
#define RY_COMPUTE(v00, v01, v10, v11, n00, n01, n10, n11) do {                \
    cplx_t _sum = (v01) + (v10);                                                \
    cplx_t _dd  = (v00) - (v11);                                               \
    (n00) = c2 * (v00) + s2 * (v11) - cs * _sum;                               \
    (n11) = s2 * (v00) + c2 * (v11) + cs * _sum;                               \
    (n10) = cs * _dd + c2 * (v10) - s2 * (v01);                                \
    (n01) = cs * _dd + c2 * (v01) - s2 * (v10);                                \
} while(0)

void ry_tiled(state_t *state, const qubit_t target, const double theta) {
    const double co = cos(theta / 2.0);
    const double si = sin(theta / 2.0);
    const double c2 = co * co, s2 = si * si, cs = co * si;

    if (target < LOG_TILE_DIM) {
        TRAVERSE_TILED_WITHIN_ROT(state, target, RY_COMPUTE);
    } else {
        TRAVERSE_TILED_CROSS_ROT(state, target, RY_COMPUTE);
    }
}
#undef RY_COMPUTE

/*
 * =====================================================================================================================
 * Rz gate (tiled): ρ → Rz(θ) ρ Rz(θ)†
 * =====================================================================================================================
 *
 * Rz(θ) = diag(e^(-iθ/2), e^(iθ/2)). Diagonal gate:
 *   ρ'_00 = ρ00, ρ'_11 = ρ11
 *   ρ'_10 = e^(iθ) ρ10, ρ'_01 = e^(-iθ) ρ01
 *
 * Rz is a diagonal gate that only applies phase to off-diagonal elements.
 * The within-tile and cross-tile paths differ from the mixing gates (Rx, Ry)
 * because diagonal elements are unchanged and only rho10/rho01 are modified.
 */
void rz_tiled(state_t *state, const qubit_t target, const double theta) {
    const double c = cos(theta), sn = sin(theta);

    if (target < LOG_TILE_DIM) {
        const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
        const gate_idx_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
        const gate_idx_t stride = (gate_idx_t)1 << target;
        const gate_idx_t tile_dim = (dim < TILE_DIM) ? dim : TILE_DIM;
        const gate_idx_t half_td = tile_dim >> 1;
        cplx_t *data = state->data;

        #pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t tr = 0; tr < n_tiles; ++tr) {
            for (gate_idx_t tc = 0; tc <= tr; ++tc) {
                cplx_t *tile = data + tile_off(tr, tc);

                for (gate_idx_t br = 0; br < half_td; ++br) {
                    gate_idx_t lr0 = insertBit0(br, target);
                    gate_idx_t lr1 = lr0 | stride;
                    for (gate_idx_t bc = 0; bc < half_td; ++bc) {
                        gate_idx_t lc0 = insertBit0(bc, target);
                        gate_idx_t lc1 = lc0 | stride;

                        /* rho10 *= e^(iθ) */
                        cplx_t *p10 = tile + elem_off(lr1, lc0);
                        double r10 = creal(*p10), i10 = cimag(*p10);
                        *p10 = CMPLX(r10 * c - i10 * sn, i10 * c + r10 * sn);

                        /* rho01 *= e^(-iθ) */
                        cplx_t *p01 = tile + elem_off(lr0, lc1);
                        double r01 = creal(*p01), i01 = cimag(*p01);
                        *p01 = CMPLX(r01 * c + i01 * sn, i01 * c - r01 * sn);
                    }
                }
            }
        }
    } else {
        /* Cross-tile traversal: iterate lr >= lc for diagonal tile-blocks */
        const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
        const gate_idx_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
        const qubit_t tile_bit = target - LOG_TILE_DIM;
        const gate_idx_t tile_stride = (gate_idx_t)1 << tile_bit;
        const gate_idx_t half_nt = n_tiles >> 1;
        cplx_t *data = state->data;

        #pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t i_tr = 0; i_tr < half_nt; ++i_tr) {
            gate_idx_t tr0 = insertBit0(i_tr, tile_bit);
            gate_idx_t tr1 = tr0 | tile_stride;
            for (gate_idx_t i_tc = 0; i_tc <= i_tr; ++i_tc) {
                gate_idx_t tc0 = insertBit0(i_tc, tile_bit);
                gate_idx_t tc1 = tc0 | tile_stride;

                cplx_t *t10 = data + tile_off(tr1, tc0);

                int lower_01 = (tr0 >= tc1);
                cplx_t *t01 = lower_01 ? data + tile_off(tr0, tc1)
                                       : data + tile_off(tc1, tr0);
                int is_diag = (i_tc == i_tr);

                for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                    gate_idx_t lc_end = is_diag ? lr + 1 : TILE_DIM;
                    for (gate_idx_t lc = 0; lc < lc_end; ++lc) {
                        int diag_block = is_diag && (lr == lc);

                        /* rho10 *= e^(iθ) */
                        cplx_t *p10 = t10 + elem_off(lr, lc);
                        double r10 = creal(*p10), i10 = cimag(*p10);
                        *p10 = CMPLX(r10 * c - i10 * sn, i10 * c + r10 * sn);

                        if (diag_block) continue;

                        /* rho01 *= e^(-iθ) */
                        cplx_t *p01;
                        double r01, i01;
                        if (lower_01) {
                            p01 = t01 + elem_off(lr, lc);
                            r01 = creal(*p01); i01 = cimag(*p01);
                            *p01 = CMPLX(r01 * c + i01 * sn, i01 * c - r01 * sn);
                        } else {
                            p01 = t01 + elem_off(lc, lr);
                            r01 = creal(*p01); i01 = cimag(*p01);
                            /* stored = conj(rho01); new = e^(iθ)*stored */
                            *p01 = CMPLX(r01 * c - i01 * sn, i01 * c + r01 * sn);
                        }
                    }
                }
            }
        }
    }
}

/*
 * =====================================================================================================================
 * P gate (tiled): ρ → P(θ) ρ P(θ)†
 * =====================================================================================================================
 *
 * P(θ) = diag(1, e^(iθ)). Diagonal gate:
 *   ρ'_00 = ρ00, ρ'_11 = ρ11
 *   ρ'_10 = e^(iθ) ρ10, ρ'_01 = e^(-iθ) ρ01
 *
 * Note: Same structure as Rz but with phase e^(iθ) instead of e^(iθ) (Rz uses full θ).
 * Actually P(θ) has the same conjugation structure as Rz(θ) since both are diagonal.
 */
void p_tiled(state_t *state, const qubit_t target, const double theta) {
    /* P(θ) conjugation is identical to Rz(θ) -- same e^(iθ) on (1,0), e^(-iθ) on (0,1) */
    rz_tiled(state, target, theta);
}

/*
 * =====================================================================================================================
 * Two-qubit gate helpers
 * =====================================================================================================================
 */

/*
 * Helpers for reading/writing from diagonal tiles through the canonical
 * (lower-triangle) position. In diagonal tiles (tr == tc), elements at
 * (lr, lc) and (lc, lr) are conjugates. During swap operations, one
 * element of a conjugate pair may be modified before the other is read,
 * causing stale data. These helpers always access the lower-triangle
 * position, avoiding stale upper-triangle reads.
 */
static inline cplx_t dtile_read(cplx_t *tile, gate_idx_t lr, gate_idx_t lc) {
    return (lr >= lc) ? tile[elem_off(lr, lc)] : conj(tile[elem_off(lc, lr)]);
}

static inline void dtile_write(cplx_t *tile, gate_idx_t lr, gate_idx_t lc, cplx_t val) {
    if (lr >= lc) {
        tile[elem_off(lr, lc)] = val;
    } else {
        tile[elem_off(lc, lr)] = conj(val);
    }
}

/*
 * =====================================================================================================================
 * Two-qubit gates
 * =====================================================================================================================
 */

/*
 * CX (CNOT) tiled: ρ' = CX·ρ·CX, zero-flop permutation gate (CX is self-adjoint).
 *
 * perm(x) = x ^ (1<<target) when bit control is set, else x.
 * ρ'[r,c] = ρ[perm(r), perm(c)].  Involution: each element is fixed or in a 2-cycle.
 *
 * For each block (br, bc), row indices r00,r01,r10,r11 and col indices c00,c01,c10,c11
 * (where bit pattern ab means control=a, target=b), the 6 swap pairs are:
 *   Pair A: (r10,c00) <-> (r11,c00)   -- always lower-tri
 *   Pair E: (r10,c10) <-> (r11,c11)   -- always lower-tri
 *   Pair B: (r11,c01) <-> (r10,c01)   -- may cross diagonal
 *   Pair F: (r11,c10) <-> (r10,c11)   -- may cross diagonal
 *   Pair C: (r00,c10) <-> (r00,c11)   -- may cross diagonal
 *   Pair D: (r01,c10) <-> (r01,c11)   -- may cross diagonal
 *
 * Three cases based on where lo=min(ctrl,tgt) and hi=max(ctrl,tgt) fall
 * relative to LOG_TILE_DIM, following the same structure as swap_tiled().
 */
void cx_tiled(state_t *state, const qubit_t control, const qubit_t target) {
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    cplx_t * restrict data = state->data;

    const qubit_t lo = (control < target) ? control : target;
    const qubit_t hi = (control > target) ? control : target;
    const gate_idx_t incr_ctrl = (gate_idx_t)1 << control;
    const gate_idx_t incr_tgt  = (gate_idx_t)1 << target;

    if (hi < LOG_TILE_DIM) {
        /*
         * Case 1: Both qubits within a tile.
         *
         * Both bits lo and hi are in local coordinates (0..LOG_TILE_DIM-1).
         * All elements of each swap pair are in the same tile.
         * Tiles store full TILE_DIM x TILE_DIM elements, so no conjugation needed.
         */
        const gate_idx_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
        const gate_idx_t tile_dim = (dim < TILE_DIM) ? dim : TILE_DIM;
        const gate_idx_t quarter_td = tile_dim >> 2;

        #pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t tr = 0; tr < n_tiles; ++tr) {
            for (gate_idx_t tc = 0; tc <= tr; ++tc) {
                cplx_t *tile = data + tile_off(tr, tc);

                for (gate_idx_t br = 0; br < quarter_td; ++br) {
                    gate_idx_t lr00 = insertBits2_0(br, lo, hi);
                    gate_idx_t lr01 = lr00 | incr_tgt;
                    gate_idx_t lr10 = lr00 | incr_ctrl;
                    gate_idx_t lr11 = lr10 | incr_tgt;

                    for (gate_idx_t bc = 0; bc < quarter_td; ++bc) {
                        gate_idx_t lc00 = insertBits2_0(bc, lo, hi);
                        gate_idx_t lc01 = lc00 | incr_tgt;
                        gate_idx_t lc10 = lc00 | incr_ctrl;
                        gate_idx_t lc11 = lc10 | incr_tgt;

                        cplx_t tmp;

                        /* Pair A: (r10,c00) <-> (r11,c00) */
                        tmp = tile[elem_off(lr10, lc00)];
                        tile[elem_off(lr10, lc00)] = tile[elem_off(lr11, lc00)];
                        tile[elem_off(lr11, lc00)] = tmp;

                        /* Pair E: (r10,c10) <-> (r11,c11) */
                        tmp = tile[elem_off(lr10, lc10)];
                        tile[elem_off(lr10, lc10)] = tile[elem_off(lr11, lc11)];
                        tile[elem_off(lr11, lc11)] = tmp;

                        /* Pair B: (r11,c01) <-> (r10,c01) */
                        tmp = tile[elem_off(lr11, lc01)];
                        tile[elem_off(lr11, lc01)] = tile[elem_off(lr10, lc01)];
                        tile[elem_off(lr10, lc01)] = tmp;

                        /* Pair F: (r11,c10) <-> (r10,c11) */
                        tmp = tile[elem_off(lr11, lc10)];
                        tile[elem_off(lr11, lc10)] = tile[elem_off(lr10, lc11)];
                        tile[elem_off(lr10, lc11)] = tmp;

                        /* Pair C: (r00,c10) <-> (r00,c11) */
                        tmp = tile[elem_off(lr00, lc10)];
                        tile[elem_off(lr00, lc10)] = tile[elem_off(lr00, lc11)];
                        tile[elem_off(lr00, lc11)] = tmp;

                        /* Pair D: (r01,c10) <-> (r01,c11) */
                        tmp = tile[elem_off(lr01, lc10)];
                        tile[elem_off(lr01, lc10)] = tile[elem_off(lr01, lc11)];
                        tile[elem_off(lr01, lc11)] = tmp;
                    }
                }
            }
        }
    }
    else if (lo < LOG_TILE_DIM) {
        /*
         * Case 2: lo within tile, hi across tiles.
         *
         * Bit lo is in local coordinates, bit hi is in tile coordinates
         * (tile_bit = hi - LOG_TILE_DIM). The CX permutation mixes local
         * and tile coordinates depending on which qubit is control vs target.
         *
         * Row/col quartets split across tiles:
         *   r00 = (tr0, lr0), r01 = (tr0, lr1) or (tr1, lr0) depending on
         *         whether lo==target or lo==control.
         *
         * If control > target (hi==control, lo==target):
         *   incr_ctrl = 1<<control affects tile coord, incr_tgt = 1<<target affects local coord
         *   r00 = (tr0, lr0), r01 = (tr0, lr1), r10 = (tr1, lr0), r11 = (tr1, lr1)
         *   c00 = (tc0, lc0), c01 = (tc0, lc1), c10 = (tc1, lc0), c11 = (tc1, lc1)
         *
         * If control < target (hi==target, lo==control):
         *   incr_tgt = 1<<target affects tile coord, incr_ctrl = 1<<control affects local coord
         *   r00 = (tr0, lr0), r01 = (tr1, lr0), r10 = (tr0, lr1), r11 = (tr1, lr1)
         *   c00 = (tc0, lc0), c01 = (tc1, lc0), c10 = (tc0, lc1), c11 = (tc1, lc1)
         *
         * We use _tile and _loc suffixes to disambiguate which increment
         * goes to tile coordinates and which goes to local coordinates.
         */
        const qubit_t tile_bit = hi - LOG_TILE_DIM;
        const gate_idx_t tile_stride = (gate_idx_t)1 << tile_bit;
        const gate_idx_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
        const gate_idx_t half_nt = n_tiles >> 1;
        const gate_idx_t tile_dim = (dim < TILE_DIM) ? dim : TILE_DIM;
        const gate_idx_t half_td = tile_dim >> 1;
        const gate_idx_t incr_loc = (gate_idx_t)1 << lo;

        /*
         * Determine whether hi is control or target.
         * hi_is_ctrl: hi==control, lo==target
         *   tile0 = ctrl=0, tile1 = ctrl=1
         *   loc0  = tgt=0,  loc1  = tgt=1
         * !hi_is_ctrl: hi==target, lo==control
         *   tile0 = tgt=0, tile1 = tgt=1
         *   loc0  = ctrl=0, loc1 = ctrl=1
         */
        const int hi_is_ctrl = (control > target);

        #pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t i_tr = 0; i_tr < half_nt; ++i_tr) {
            gate_idx_t tr0 = insertBit0(i_tr, tile_bit);
            gate_idx_t tr1 = tr0 | tile_stride;

            for (gate_idx_t i_tc = 0; i_tc <= i_tr; ++i_tc) {
                gate_idx_t tc0 = insertBit0(i_tc, tile_bit);
                gate_idx_t tc1 = tc0 | tile_stride;

                /* Tile pointers: T00, T11, T10 always in lower triangle */
                cplx_t *t00 = data + tile_off(tr0, tc0);
                cplx_t *t11 = data + tile_off(tr1, tc1);
                cplx_t *t10 = data + tile_off(tr1, tc0);

                /* T01 = T(tr0, tc1): only stored if tr0 >= tc1 */
                int lower_01 = (tr0 >= tc1);
                cplx_t *t01 = lower_01 ? data + tile_off(tr0, tc1)
                                       : data + tile_off(tc1, tr0);

                int is_diag = (i_tc == i_tr);

                /*
                 * Map the CX 6 swap pairs to (tile, local) coordinates.
                 *
                 * Notation: T[tile_row, tile_col][local_row, local_col]
                 *
                 * For hi_is_ctrl (tile bit = ctrl bit, local bit = tgt bit):
                 *   r_ab = (tr_a, lr_b), c_ab = (tc_a, lc_b) where a=ctrl, b=tgt
                 *   Pair A: T10[lr0,lc0] <-> T11[lr0,lc0]  -- t10 <-> T(tr1,tc0+tgt)
                 *     But wait: r10=(tr1,lr0), c00=(tc0,lc0) => T(tr1,tc0)[lr0,lc0]
                 *               r11=(tr1,lr1), c00=(tc0,lc0) => T(tr1,tc0)[lr1,lc0]
                 *   So Pair A: T(tr1,tc0)[lr0,lc0] <-> T(tr1,tc0)[lr1,lc0]  (within t10)
                 *   Pair E: T(tr1,tc1)[lr0,lc0] <-> T(tr1,tc1)[lr1,lc1]     (within t11... wait)
                 *     r10=(tr1,lr0), c10=(tc1,lc0) => T(tr1,tc1)[lr0,lc0]
                 *     r11=(tr1,lr1), c11=(tc1,lc1) => T(tr1,tc1)[lr1,lc1]
                 *   Pair B: T(tr1,tc0)[lr1,lc1] <-> T(tr1,tc0)[lr0,lc1]     (within t10)
                 *   Pair F: T(tr1,tc1)[lr1,lc0] <-> T(tr1,tc1)[lr0,lc1]     (within t11... wait)
                 *     r11=(tr1,lr1), c10=(tc1,lc0) => T(tr1,tc1)[lr1,lc0]
                 *     r10=(tr1,lr0), c11=(tc1,lc1) => T(tr1,tc1)[lr0,lc1]
                 *   Pair C: T(tr0,tc1)[lr0,lc0] <-> T(tr0,tc1)[lr0,lc1]     (t01)
                 *     r00=(tr0,lr0), c10=(tc1,lc0) => T(tr0,tc1)[lr0,lc0]
                 *     r00=(tr0,lr0), c11=(tc1,lc1) => T(tr0,tc1)[lr0,lc1]
                 *   Pair D: T(tr0,tc1)[lr1,lc0] <-> T(tr0,tc1)[lr1,lc1]     (t01)
                 *     r01=(tr0,lr1), c10=(tc1,lc0) => T(tr0,tc1)[lr1,lc0]
                 *     r01=(tr0,lr1), c11=(tc1,lc1) => T(tr0,tc1)[lr1,lc1]
                 *
                 * So for hi_is_ctrl: Pairs A,B within t10; Pairs E,F within t11;
                 *   Pairs C,D within t01 (which may be upper-tri).
                 *   t00 is UNTOUCHED.
                 *
                 * For !hi_is_ctrl (tile bit = tgt bit, local bit = ctrl bit):
                 *   r_ab = (tr_b, lr_a), c_ab = (tc_b, lc_a)
                 *   r00=(tr0,lr0), r01=(tr1,lr0), r10=(tr0,lr1), r11=(tr1,lr1)
                 *   c00=(tc0,lc0), c01=(tc1,lc0), c10=(tc0,lc1), c11=(tc1,lc1)
                 *
                 *   Pair A: (r10,c00)<->(r11,c00) = T(tr0,tc0)[lr1,lc0]<->T(tr1,tc0)[lr1,lc0]
                 *     i.e. t00[lr1,lc0] <-> t10[lr1,lc0]
                 *   Pair E: (r10,c10)<->(r11,c11) = T(tr0,tc0)[lr1,lc1]<->T(tr1,tc1)[lr1,lc1]
                 *     i.e. t00[lr1,lc1] <-> t11[lr1,lc1]
                 *   Pair B: (r11,c01)<->(r10,c01) = T(tr1,tc1)[lr1,lc0]<->T(tr0,tc1)[lr1,lc0]
                 *     i.e. t11[lr1,lc0] <-> t01[lr1,lc0]
                 *   Pair F: (r11,c10)<->(r10,c11) = T(tr1,tc0)[lr1,lc1]<->T(tr0,tc1)[lr1,lc1]
                 *     i.e. t10[lr1,lc1] <-> t01[lr1,lc1]
                 *   Pair C: (r00,c10)<->(r00,c11) = T(tr0,tc0)[lr0,lc1]<->T(tr0,tc1)[lr0,lc1]
                 *     i.e. t00[lr0,lc1] <-> t01[lr0,lc1]
                 *   Pair D: (r01,c10)<->(r01,c11) = T(tr1,tc0)[lr0,lc1]<->T(tr1,tc1)[lr0,lc1]
                 *     i.e. t10[lr0,lc1] <-> t11[lr0,lc1]
                 *
                 * So for !hi_is_ctrl: Pairs A,E touch t00; Pairs B,F touch t01;
                 *   all 6 pairs span 2 different tiles each.
                 */

                if (hi_is_ctrl) {
                    /*
                     * hi == control, lo == target.
                     * Pairs A,B are within t10 (always lower-tri off-diag tile).
                     * Pairs E,F are within t11.
                     * Pairs C,D are within t01 (may be upper-tri).
                     *
                     * When is_diag: t00 and t11 are diagonal tiles.
                     * CRITICAL: t01 and t10 point to the SAME physical tile
                     * (T(tc1,tr0) = T(tr1,tr0) = t10 when tr0==tc0).
                     * Pair C's swap t01[lc0,lr0]<->t01[lc1,lr0] is the
                     * transpose of Pair A's swap t10[lr0,lc0]<->t10[lr1,lc0]
                     * from the (b_lc,b_lr) iteration. Iterating all (b_lr,b_lc)
                     * for Pairs A,B already covers the work of Pairs C,D.
                     * So: skip Pairs C,D entirely when is_diag.
                     *
                     * Pairs E,F operate on t11 (diagonal tile): use
                     * dtile_read/dtile_write. Because t11 is diagonal,
                     * dtile helpers canonicalize through the lower triangle.
                     * Iterations (b_lr,b_lc) and (b_lc,b_lr) then access the
                     * SAME physical positions, causing each swap to be done
                     * twice (i.e. undone). Guard with b_lc <= b_lr so each
                     * canonical pair is processed exactly once.
                     */
                    for (gate_idx_t b_lr = 0; b_lr < half_td; ++b_lr) {
                        gate_idx_t lr0 = insertBit0(b_lr, lo);
                        gate_idx_t lr1 = lr0 | incr_loc;

                        for (gate_idx_t b_lc = 0; b_lc < half_td; ++b_lc) {
                            gate_idx_t lc0 = insertBit0(b_lc, lo);
                            gate_idx_t lc1 = lc0 | incr_loc;

                            cplx_t tmp;

                            /* Pair A: t10[lr0,lc0] <-> t10[lr1,lc0] */
                            tmp = t10[elem_off(lr0, lc0)];
                            t10[elem_off(lr0, lc0)] = t10[elem_off(lr1, lc0)];
                            t10[elem_off(lr1, lc0)] = tmp;

                            /* Pair B: t10[lr1,lc1] <-> t10[lr0,lc1] */
                            tmp = t10[elem_off(lr1, lc1)];
                            t10[elem_off(lr1, lc1)] = t10[elem_off(lr0, lc1)];
                            t10[elem_off(lr0, lc1)] = tmp;

                            if (is_diag && b_lc <= b_lr) {
                                /* t11 is a diagonal tile; use dtile helpers.
                                 * Restricted to b_lc <= b_lr: see comment above. */

                                /* Pair E: t11[lr0,lc0] <-> t11[lr1,lc1] */
                                cplx_t vE_a = dtile_read(t11, lr0, lc0);
                                cplx_t vE_b = dtile_read(t11, lr1, lc1);
                                dtile_write(t11, lr0, lc0, vE_b);
                                dtile_write(t11, lr1, lc1, vE_a);

                                /* Pair F: t11[lr1,lc0] <-> t11[lr0,lc1] */
                                cplx_t vF_a = dtile_read(t11, lr1, lc0);
                                cplx_t vF_b = dtile_read(t11, lr0, lc1);
                                dtile_write(t11, lr1, lc0, vF_b);
                                dtile_write(t11, lr0, lc1, vF_a);

                                /* Pairs C,D: SKIPPED when is_diag.
                                 * t01 aliases t10, so Pair C's transposed swap
                                 * on t01 is already covered by Pair A operating
                                 * on t10 across all (b_lr,b_lc) iterations. */
                            } else if (!is_diag) {
                                /* Pair E: t11[lr0,lc0] <-> t11[lr1,lc1] */
                                tmp = t11[elem_off(lr0, lc0)];
                                t11[elem_off(lr0, lc0)] = t11[elem_off(lr1, lc1)];
                                t11[elem_off(lr1, lc1)] = tmp;

                                /* Pair F: t11[lr1,lc0] <-> t11[lr0,lc1] */
                                tmp = t11[elem_off(lr1, lc0)];
                                t11[elem_off(lr1, lc0)] = t11[elem_off(lr0, lc1)];
                                t11[elem_off(lr0, lc1)] = tmp;

                                /* Pairs C,D: within t01 (may be upper-tri) */
                                if (lower_01) {
                                    /* Pair C: t01[lr0,lc0] <-> t01[lr0,lc1] */
                                    tmp = t01[elem_off(lr0, lc0)];
                                    t01[elem_off(lr0, lc0)] = t01[elem_off(lr0, lc1)];
                                    t01[elem_off(lr0, lc1)] = tmp;

                                    /* Pair D: t01[lr1,lc0] <-> t01[lr1,lc1] */
                                    tmp = t01[elem_off(lr1, lc0)];
                                    t01[elem_off(lr1, lc0)] = t01[elem_off(lr1, lc1)];
                                    t01[elem_off(lr1, lc1)] = tmp;
                                } else {
                                    /* t01 = T(tc1, tr0), access via transposed+conj.
                                     * T(tr0,tc1)[lr,lc] = conj(t01[lc,lr])
                                     * Swap conj(t01[lc0,lr0]) <-> conj(t01[lc1,lr0])
                                     * => swap t01[lc0,lr0] <-> t01[lc1,lr0] */
                                    /* Pair C: t01[lc0,lr0] <-> t01[lc1,lr0] */
                                    tmp = t01[elem_off(lc0, lr0)];
                                    t01[elem_off(lc0, lr0)] = t01[elem_off(lc1, lr0)];
                                    t01[elem_off(lc1, lr0)] = tmp;

                                    /* Pair D: t01[lc0,lr1] <-> t01[lc1,lr1] */
                                    tmp = t01[elem_off(lc0, lr1)];
                                    t01[elem_off(lc0, lr1)] = t01[elem_off(lc1, lr1)];
                                    t01[elem_off(lc1, lr1)] = tmp;
                                }
                            }
                        }
                    }

                    /* Mirror upper triangle of diagonal tiles. When is_diag,
                     * t11 is a diagonal tile whose lower triangle was modified
                     * by pairs E,F. t00 is untouched by CX when hi==control.
                     * dtile_write maintained both triangles, but do a full
                     * mirror pass for consistency. */
                    if (is_diag) {
                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            for (gate_idx_t lc = lr + 1; lc < TILE_DIM; ++lc) {
                                t11[elem_off(lr, lc)] = conj(t11[elem_off(lc, lr)]);
                            }
                        }
                    }
                } else {
                    /*
                     * hi == target, lo == control.
                     * Pair A: t00[lr1,lc0] <-> t10[lr1,lc0]
                     * Pair E: t00[lr1,lc1] <-> t11[lr1,lc1]
                     * Pair B: t11[lr1,lc0] <-> t01[lr1,lc0]
                     * Pair F: t10[lr1,lc1] <-> t01[lr1,lc1]
                     * Pair C: t00[lr0,lc1] <-> t01[lr0,lc1]
                     * Pair D: t10[lr0,lc1] <-> t11[lr0,lc1]
                     *
                     * Here lr0/lr1 refer to ctrl=0/1 in local coord.
                     * Pairs A,E modify t00 which may be diagonal.
                     * Pairs B,C,F involve t01 which may be upper-tri.
                     */
                    if (is_diag) {
                        /* t00 and t11 are diagonal tiles. t01==t10 by storage
                         * (when is_diag, lower_01==0, t01=T(tc1,tr0)=T(tr1,tr0)=t10).
                         *
                         * Blocks (b_lr, b_lc) and (b_lc, b_lr) produce the same
                         * physical swaps due to Hermiticity, so process only b_lc<=b_lr. */
                        for (gate_idx_t b_lr = 0; b_lr < half_td; ++b_lr) {
                            gate_idx_t lr0 = insertBit0(b_lr, lo);
                            gate_idx_t lr1 = lr0 | incr_loc;

                            for (gate_idx_t b_lc = 0; b_lc <= b_lr; ++b_lc) {
                                gate_idx_t lc0 = insertBit0(b_lc, lo);
                                gate_idx_t lc1 = lc0 | incr_loc;

                                /* Read all values needed before writing.
                                 * t01 = T(tr1,tr0) = t10, accessed as conj(t10[lc,lr])
                                 * for logical element T(tr0,tc1)[lr,lc]. */
                                cplx_t pA_a = dtile_read(t00, lr1, lc0);
                                cplx_t pA_b = t10[elem_off(lr1, lc0)];

                                cplx_t pE_a = dtile_read(t00, lr1, lc1);
                                cplx_t pE_b = dtile_read(t11, lr1, lc1);

                                /* t01 logical [lr1,lc0] = conj(t10[lc0,lr1]) */
                                cplx_t pB_a = dtile_read(t11, lr1, lc0);
                                cplx_t pB_b = conj(t10[elem_off(lc0, lr1)]);

                                /* t01 logical [lr1,lc1] = conj(t10[lc1,lr1]) */
                                cplx_t pF_a = t10[elem_off(lr1, lc1)];
                                cplx_t pF_b = conj(t10[elem_off(lc1, lr1)]);

                                /* t01 logical [lr0,lc1] = conj(t10[lc1,lr0]) */
                                cplx_t pC_a = dtile_read(t00, lr0, lc1);
                                cplx_t pC_b = conj(t10[elem_off(lc1, lr0)]);

                                cplx_t pD_a = t10[elem_off(lr0, lc1)];
                                cplx_t pD_b = dtile_read(t11, lr0, lc1);

                                /* Write back */
                                /* Pair A */
                                dtile_write(t00, lr1, lc0, pA_b);
                                t10[elem_off(lr1, lc0)] = pA_a;

                                /* Pair E */
                                dtile_write(t00, lr1, lc1, pE_b);
                                dtile_write(t11, lr1, lc1, pE_a);

                                /* Pair B: t11[lr1,lc0] <-> t01[lr1,lc0] */
                                dtile_write(t11, lr1, lc0, pB_b);
                                t10[elem_off(lc0, lr1)] = conj(pB_a);

                                /* Pair F: t10[lr1,lc1] <-> t01[lr1,lc1] */
                                t10[elem_off(lr1, lc1)] = pF_b;
                                t10[elem_off(lc1, lr1)] = conj(pF_a);

                                /* Pair C: t00[lr0,lc1] <-> t01[lr0,lc1] */
                                dtile_write(t00, lr0, lc1, pC_b);
                                t10[elem_off(lc1, lr0)] = conj(pC_a);

                                /* Pair D */
                                t10[elem_off(lr0, lc1)] = pD_b;
                                dtile_write(t11, lr0, lc1, pD_a);
                            }
                        }

                        /* Mirror upper triangle of diagonal tiles */
                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            for (gate_idx_t lc = lr + 1; lc < TILE_DIM; ++lc) {
                                t00[elem_off(lr, lc)] = conj(t00[elem_off(lc, lr)]);
                                t11[elem_off(lr, lc)] = conj(t11[elem_off(lc, lr)]);
                            }
                        }
                    } else {
                        for (gate_idx_t b_lr = 0; b_lr < half_td; ++b_lr) {
                            gate_idx_t lr0 = insertBit0(b_lr, lo);
                            gate_idx_t lr1 = lr0 | incr_loc;

                            for (gate_idx_t b_lc = 0; b_lc < half_td; ++b_lc) {
                                gate_idx_t lc0 = insertBit0(b_lc, lo);
                                gate_idx_t lc1 = lc0 | incr_loc;

                                cplx_t tmp;

                                /* Pair A: t00[lr1,lc0] <-> t10[lr1,lc0] */
                                tmp = t00[elem_off(lr1, lc0)];
                                t00[elem_off(lr1, lc0)] = t10[elem_off(lr1, lc0)];
                                t10[elem_off(lr1, lc0)] = tmp;

                                /* Pair E: t00[lr1,lc1] <-> t11[lr1,lc1] */
                                tmp = t00[elem_off(lr1, lc1)];
                                t00[elem_off(lr1, lc1)] = t11[elem_off(lr1, lc1)];
                                t11[elem_off(lr1, lc1)] = tmp;

                                /* Pair D: t10[lr0,lc1] <-> t11[lr0,lc1] */
                                tmp = t10[elem_off(lr0, lc1)];
                                t10[elem_off(lr0, lc1)] = t11[elem_off(lr0, lc1)];
                                t11[elem_off(lr0, lc1)] = tmp;

                                /* Pairs B, C, F involve t01 */
                                if (lower_01) {
                                    /* Pair B: t11[lr1,lc0] <-> t01[lr1,lc0] */
                                    tmp = t11[elem_off(lr1, lc0)];
                                    t11[elem_off(lr1, lc0)] = t01[elem_off(lr1, lc0)];
                                    t01[elem_off(lr1, lc0)] = tmp;

                                    /* Pair F: t10[lr1,lc1] <-> t01[lr1,lc1] */
                                    tmp = t10[elem_off(lr1, lc1)];
                                    t10[elem_off(lr1, lc1)] = t01[elem_off(lr1, lc1)];
                                    t01[elem_off(lr1, lc1)] = tmp;

                                    /* Pair C: t00[lr0,lc1] <-> t01[lr0,lc1] */
                                    tmp = t00[elem_off(lr0, lc1)];
                                    t00[elem_off(lr0, lc1)] = t01[elem_off(lr0, lc1)];
                                    t01[elem_off(lr0, lc1)] = tmp;
                                } else {
                                    /* t01 = T(tc1, tr0): element [lr,lc] of
                                     * logical T(tr0,tc1) stored as conj(t01[lc,lr]) */

                                    /* Pair B: t11[lr1,lc0] <-> conj(t01[lc0,lr1]) */
                                    {
                                        cplx_t val_a = t11[elem_off(lr1, lc0)];
                                        cplx_t val_b = conj(t01[elem_off(lc0, lr1)]);
                                        t11[elem_off(lr1, lc0)] = val_b;
                                        t01[elem_off(lc0, lr1)] = conj(val_a);
                                    }

                                    /* Pair F: t10[lr1,lc1] <-> conj(t01[lc1,lr1]) */
                                    {
                                        cplx_t val_a = t10[elem_off(lr1, lc1)];
                                        cplx_t val_b = conj(t01[elem_off(lc1, lr1)]);
                                        t10[elem_off(lr1, lc1)] = val_b;
                                        t01[elem_off(lc1, lr1)] = conj(val_a);
                                    }

                                    /* Pair C: t00[lr0,lc1] <-> conj(t01[lc1,lr0]) */
                                    {
                                        cplx_t val_a = t00[elem_off(lr0, lc1)];
                                        cplx_t val_b = conj(t01[elem_off(lc1, lr0)]);
                                        t00[elem_off(lr0, lc1)] = val_b;
                                        t01[elem_off(lc1, lr0)] = conj(val_a);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else {
        /*
         * Case 3: Both qubits across tiles.
         *
         * Both bits lo and hi are in tile coordinates. The CX permutation
         * only affects tile indices, not local coordinates within tiles.
         *
         * tile_ctrl = control - LOG_TILE_DIM, tile_tgt = target - LOG_TILE_DIM.
         * tile_lo = lo - LOG_TILE_DIM, tile_hi = hi - LOG_TILE_DIM.
         *
         * We iterate tile-row quartets and tile-col quartets.
         * Row/col tile quartets: tr00, tr01, tr10, tr11 (bit pattern = ctrl,tgt).
         *
         * The 6 swap pairs at tile level:
         *   Pair A: T(tr10,tc00) <-> T(tr11,tc00)
         *   Pair E: T(tr10,tc10) <-> T(tr11,tc11)
         *   Pair B: T(tr11,tc01) <-> T(tr10,tc01)
         *   Pair F: T(tr11,tc10) <-> T(tr10,tc11)
         *   Pair C: T(tr00,tc10) <-> T(tr00,tc11)
         *   Pair D: T(tr01,tc10) <-> T(tr01,tc11)
         */
        const qubit_t tile_lo = lo - LOG_TILE_DIM;
        const qubit_t tile_hi = hi - LOG_TILE_DIM;
        const qubit_t tile_ctrl = control - LOG_TILE_DIM;
        const qubit_t tile_tgt  = target - LOG_TILE_DIM;
        const gate_idx_t incr_tile_ctrl = (gate_idx_t)1 << tile_ctrl;
        const gate_idx_t incr_tile_tgt  = (gate_idx_t)1 << tile_tgt;
        const gate_idx_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
        const gate_idx_t quarter_nt = n_tiles >> 2;

        #pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t b_tr = 0; b_tr < quarter_nt; ++b_tr) {
            gate_idx_t tr00 = insertBits2_0(b_tr, tile_lo, tile_hi);
            gate_idx_t tr01 = tr00 | incr_tile_tgt;
            gate_idx_t tr10 = tr00 | incr_tile_ctrl;
            gate_idx_t tr11 = tr10 | incr_tile_tgt;

            for (gate_idx_t b_tc = 0; b_tc <= b_tr; ++b_tc) {
                gate_idx_t tc00 = insertBits2_0(b_tc, tile_lo, tile_hi);
                gate_idx_t tc01 = tc00 | incr_tile_tgt;
                gate_idx_t tc10 = tc00 | incr_tile_ctrl;
                gate_idx_t tc11 = tc10 | incr_tile_tgt;

                int is_diag = (b_tc == b_tr);

                /* ---- Pair A: T(tr10,tc00) <-> T(tr11,tc00) ----
                 * Both always lower-tri (tr10 > tr00 >= tc00, tr11 > tr00 >= tc00). */
                {
                    cplx_t *tA = data + tile_off(tr10, tc00);
                    cplx_t *tB = data + tile_off(tr11, tc00);
                    for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                        cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                    }
                }

                /* ---- Pair E: T(tr10,tc10) <-> T(tr11,tc11) ----
                 * Both always lower-tri. */
                {
                    cplx_t *tA = data + tile_off(tr10, tc10);
                    cplx_t *tB = data + tile_off(tr11, tc11);
                    for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                        cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                    }
                }

                /* ---- Pair B: T(tr11,tc01) <-> T(tr10,tc01) ----
                 * T(tr11,tc01): always lower-tri (tr11 has both ctrl+tgt bits set).
                 * T(tr10,tc01): lower-tri iff tr10 >= tc01.
                 *   When tile_ctrl < tile_tgt: tr10 = tr00|ctrl_bit (smaller),
                 *   tc01 = tc00|tgt_bit (larger), so tr10 < tc01 is possible.
                 *
                 * When is_diag and !lower_b: T(tr10,tc01) stored as T(tc01,tr10).
                 * Unlike Pair F, these are two distinct tiles (not the same physical
                 * tile), so this is NOT hermitian transpose in place. */
                {
                    cplx_t *tA = data + tile_off(tr11, tc01);
                    int lower_b = (tr10 >= tc01);
                    if (lower_b) {
                        cplx_t *tB = data + tile_off(tr10, tc01);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                        }
                    } else {
                        /* T(tr10,tc01) upper-tri, stored as T(tc01, tr10) */
                        cplx_t *tB = data + tile_off(tc01, tr10);
                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                cplx_t val_a = tA[elem_off(lr, lc)];
                                cplx_t val_b = conj(tB[elem_off(lc, lr)]);
                                tA[elem_off(lr, lc)] = val_b;
                                tB[elem_off(lc, lr)] = conj(val_a);
                            }
                        }
                    }
                }

                /* ---- Pair F: T(tr11,tc10) <-> T(tr10,tc11) ----
                 * T(tr11,tc10): always lower-tri.
                 * T(tr10,tc11): lower-tri iff tr10 >= tc11.
                 *
                 * When is_diag (tr00==tc00): tr10=tc10, tc11=tr11, so
                 * T(tr10,tc11)=T(tc10,tr11). Since tc10<tr11, T(tr10,tr11)
                 * is upper-tri, stored as T(tr11,tr10). The swap pair
                 * becomes T(tr11,tc10) <-> T(tc10,tr11) where both
                 * physically point to T(tr11,tr10). This is Hermitian
                 * transpose in place. */
                {
                    cplx_t *tA = data + tile_off(tr11, tc10);
                    int lower_f = (tr10 >= tc11);
                    if (lower_f) {
                        cplx_t *tB = data + tile_off(tr10, tc11);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                        }
                    } else if (is_diag) {
                        /* tA = T(tr11, tc10) = T(tr11, tr10). Hermitian transpose in place.
                         * Logical swap: tA[lr,lc] <-> conj(tA[lc,lr]) */
                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            tA[elem_off(lr, lr)] = conj(tA[elem_off(lr, lr)]);
                            for (gate_idx_t lc = lr + 1; lc < TILE_DIM; ++lc) {
                                cplx_t val_a = tA[elem_off(lr, lc)];
                                cplx_t val_b = tA[elem_off(lc, lr)];
                                tA[elem_off(lr, lc)] = conj(val_b);
                                tA[elem_off(lc, lr)] = conj(val_a);
                            }
                        }
                    } else {
                        /* T(tr10,tc11) upper-tri, stored as T(tc11, tr10) */
                        cplx_t *tB = data + tile_off(tc11, tr10);
                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                cplx_t val_a = tA[elem_off(lr, lc)];
                                cplx_t val_b = conj(tB[elem_off(lc, lr)]);
                                tA[elem_off(lr, lc)] = val_b;
                                tB[elem_off(lc, lr)] = conj(val_a);
                            }
                        }
                    }
                }

                /* ---- Pairs C+D ----
                 * Pair C: T(tr00,tc10) <-> T(tr00,tc11)
                 * Pair D: T(tr01,tc10) <-> T(tr01,tc11)
                 *
                 * Lower-tri conditions:
                 *   T(tr00,tc10): iff tr00 >= tc10    [lo_10]
                 *   T(tr00,tc11): iff tr00 >= tc11    [lo_11]
                 *   T(tr01,tc10): iff tr01 >= tc10    [lo2_10]
                 *   T(tr01,tc11): iff tr01 >= tc11    [lo2_11]
                 *
                 * When is_diag: tr00==tc00, so tc10>tr00 and tc11>tr00,
                 * meaning lo_10==0 and lo_11==0.
                 *
                 * Pairs C+D, when is_diag, map to the same physical
                 * storage as Pairs A+B (the stored lower-tri reps coincide),
                 * so executing them would undo A+B. Skip.
                 */
                if (!is_diag) {
                    int lo_10 = (tr00 >= tc10);
                    int lo_11 = (tr00 >= tc11);

                    /* Pair C: T(tr00,tc10) <-> T(tr00,tc11) */
                    if (lo_10 && lo_11) {
                        cplx_t *tCa = data + tile_off(tr00, tc10);
                        cplx_t *tCb = data + tile_off(tr00, tc11);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tCa[i]; tCa[i] = tCb[i]; tCb[i] = tmp;
                        }
                    }
                    else if (lo_10 && !lo_11) {
                        /* T(tr00,tc10) lower, T(tr00,tc11) upper -> T(tc11,tr00) */
                        cplx_t *tCa = data + tile_off(tr00, tc10);
                        cplx_t *tCb = data + tile_off(tc11, tr00);
                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                cplx_t val_a = tCa[elem_off(lr, lc)];
                                cplx_t val_b = conj(tCb[elem_off(lc, lr)]);
                                tCa[elem_off(lr, lc)] = val_b;
                                tCb[elem_off(lc, lr)] = conj(val_a);
                            }
                        }
                    }
                    else {
                        /* Both upper-tri: swap stored reps directly.
                         * T(tr00,tc10) stored as T(tc10,tr00).
                         * T(tr00,tc11) stored as T(tc11,tr00). */
                        cplx_t *tCa = data + tile_off(tc10, tr00);
                        cplx_t *tCb = data + tile_off(tc11, tr00);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tCa[i]; tCa[i] = tCb[i]; tCb[i] = tmp;
                        }
                    }

                    /* Pair D: T(tr01,tc10) <-> T(tr01,tc11) */
                    int lo2_10 = (tr01 >= tc10);
                    int lo2_11 = (tr01 >= tc11);

                    if (lo2_10 && lo2_11) {
                        cplx_t *tDa = data + tile_off(tr01, tc10);
                        cplx_t *tDb = data + tile_off(tr01, tc11);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tDa[i]; tDa[i] = tDb[i]; tDb[i] = tmp;
                        }
                    }
                    else if (lo2_10 && !lo2_11) {
                        cplx_t *tDa = data + tile_off(tr01, tc10);
                        cplx_t *tDb = data + tile_off(tc11, tr01);
                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                cplx_t val_a = tDa[elem_off(lr, lc)];
                                cplx_t val_b = conj(tDb[elem_off(lc, lr)]);
                                tDa[elem_off(lr, lc)] = val_b;
                                tDb[elem_off(lc, lr)] = conj(val_a);
                            }
                        }
                    }
                    else {
                        cplx_t *tDa = data + tile_off(tc10, tr01);
                        cplx_t *tDb = data + tile_off(tc11, tr01);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tDa[i]; tDa[i] = tDb[i]; tDb[i] = tmp;
                        }
                    }
                }
            }
        }
    }
}

void cy_tiled(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "cy: tiled not yet implemented");
}

void cz_tiled(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "cz: tiled not yet implemented");
}

void cs_tiled(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "cs: tiled not yet implemented");
}

void csdg_tiled(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "csdg: tiled not yet implemented");
}

void ch_tiled(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "ch: tiled not yet implemented");
}

void chy_tiled(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "chy: tiled not yet implemented");
}

void ct_tiled(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "ct: tiled not yet implemented");
}

void ctdg_tiled(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "ctdg: tiled not yet implemented");
}

void cp_tiled(state_t *state, const qubit_t control, const qubit_t target, const double theta) {
    GATE_VALIDATE(0, "cp: tiled not yet implemented");
}

void cpdg_tiled(state_t *state, const qubit_t control, const qubit_t target, const double theta) {
    GATE_VALIDATE(0, "cpdg: tiled not yet implemented");
}

/*
 * =====================================================================================================================
 * SWAP gate (tiled): rho -> SWAP * rho * SWAP
 * =====================================================================================================================
 *
 * SWAP on qubits (lo, hi) permutes basis states by swapping the bit values at
 * positions lo and hi. For density matrix element rho[r,c]:
 *   rho'[perm(r), perm(c)] = rho[r,c]
 * where perm(x) exchanges bits lo and hi in x. SWAP is Hermitian and unitary.
 *
 * Since SWAP is a self-inverse permutation, every element is either a fixed point
 * (bits lo and hi have the same value in the index) or part of a 2-cycle. We only
 * need to swap each 2-cycle pair once.
 *
 * The SWAP permutation acts on 4x4 blocks indexed by the 2-bit pattern at positions
 * (lo, hi). Of the 16 element positions, 4 are fixed (00->00, 11->11) and 12 form
 * 6 swap pairs (see swap_packed for the derivation):
 *   Pair 1: (r01, c00) <-> (r10, c00)
 *   Pair 2: (r01, c01) <-> (r10, c10)
 *   Pair 3: (r11, c01) <-> (r11, c10)
 *   Pair 4: (r10, c01) <-> (r01, c10)  -- may cross the diagonal
 *   Pair 5: (r00, c01) <-> (r00, c10)  -- may cross the diagonal
 *   Pair 6: (r10, c11) <-> (r01, c11)  -- may cross the diagonal
 *
 * Three cases based on where lo and hi fall relative to LOG_TILE_DIM:
 *   Case 1: hi < LOG_TILE_DIM   -- both qubits within a tile (local coordinates)
 *   Case 2: lo < LOG_TILE_DIM, hi >= LOG_TILE_DIM -- lo within tile, hi across tiles
 *   Case 3: lo >= LOG_TILE_DIM  -- both qubits across tiles
 */

void swap_tiled(state_t *state, const qubit_t q1, const qubit_t q2) {
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    cplx_t * restrict data = state->data;

    const qubit_t lo = (q1 < q2) ? q1 : q2;
    const qubit_t hi = (q1 > q2) ? q1 : q2;

    if (hi < LOG_TILE_DIM) {
        /*
         * Case 1: Both qubits within a tile.
         *
         * Both bits lo and hi are in local coordinates (0..LOG_TILE_DIM-1).
         * The permutation only affects local row/col within each tile.
         * All 4 elements of each swap pair are in the same tile.
         * Tiles store full TILE_DIM x TILE_DIM elements, so no conjugation needed.
         *
         * For diagonal tiles: the swaps modify both the lower and upper triangle
         * positions. Since SWAP preserves Hermiticity and both conjugate positions
         * are swapped consistently (both are within the same tile), the result is
         * automatically Hermitian-consistent.
         */
        const gate_idx_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
        const gate_idx_t tile_dim = (dim < TILE_DIM) ? dim : TILE_DIM;
        const gate_idx_t quarter_td = tile_dim >> 2;
        const gate_idx_t incr_lo = (gate_idx_t)1 << lo;
        const gate_idx_t incr_hi = (gate_idx_t)1 << hi;

        #pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t tr = 0; tr < n_tiles; ++tr) {
            for (gate_idx_t tc = 0; tc <= tr; ++tc) {
                cplx_t *tile = data + tile_off(tr, tc);

                for (gate_idx_t br = 0; br < quarter_td; ++br) {
                    gate_idx_t lr00 = insertBits2_0(br, lo, hi);
                    gate_idx_t lr01 = lr00 | incr_lo;
                    gate_idx_t lr10 = lr00 | incr_hi;
                    gate_idx_t lr11 = lr10 | incr_lo;

                    for (gate_idx_t bc = 0; bc < quarter_td; ++bc) {
                        gate_idx_t lc00 = insertBits2_0(bc, lo, hi);
                        gate_idx_t lc01 = lc00 | incr_lo;
                        gate_idx_t lc10 = lc00 | incr_hi;
                        gate_idx_t lc11 = lc10 | incr_lo;

                        cplx_t tmp;

                        /* Pair 1: (r01,c00) <-> (r10,c00) */
                        tmp = tile[elem_off(lr01, lc00)];
                        tile[elem_off(lr01, lc00)] = tile[elem_off(lr10, lc00)];
                        tile[elem_off(lr10, lc00)] = tmp;

                        /* Pair 2: (r01,c01) <-> (r10,c10) */
                        tmp = tile[elem_off(lr01, lc01)];
                        tile[elem_off(lr01, lc01)] = tile[elem_off(lr10, lc10)];
                        tile[elem_off(lr10, lc10)] = tmp;

                        /* Pair 3: (r11,c01) <-> (r11,c10) */
                        tmp = tile[elem_off(lr11, lc01)];
                        tile[elem_off(lr11, lc01)] = tile[elem_off(lr11, lc10)];
                        tile[elem_off(lr11, lc10)] = tmp;

                        /* Pair 4: (r10,c01) <-> (r01,c10) */
                        tmp = tile[elem_off(lr10, lc01)];
                        tile[elem_off(lr10, lc01)] = tile[elem_off(lr01, lc10)];
                        tile[elem_off(lr01, lc10)] = tmp;

                        /* Pair 5: (r00,c01) <-> (r00,c10) */
                        tmp = tile[elem_off(lr00, lc01)];
                        tile[elem_off(lr00, lc01)] = tile[elem_off(lr00, lc10)];
                        tile[elem_off(lr00, lc10)] = tmp;

                        /* Pair 6: (r01,c11) <-> (r10,c11) */
                        tmp = tile[elem_off(lr01, lc11)];
                        tile[elem_off(lr01, lc11)] = tile[elem_off(lr10, lc11)];
                        tile[elem_off(lr10, lc11)] = tmp;
                    }
                }
            }
        }
    }
    else if (lo < LOG_TILE_DIM) {
        /*
         * Case 2: lo within tile, hi across tiles.
         *
         * Bit lo is in local coordinates, bit hi is in tile coordinates
         * (tile_bit = hi - LOG_TILE_DIM). The permutation swaps a local bit
         * with a tile-coordinate bit, so elements move between tiles.
         *
         * For row index r with tile-row tr and local-row lr:
         *   perm(r) swaps bit lo in lr with bit tile_bit in tr
         *
         * Non-trivial action: elements in (tr0, lr|lo=1) swap with
         * elements in (tr1, lr&~lo), and similarly for columns.
         *
         * We iterate tile-row pairs (tr0, tr1) and tile-col pairs (tc0, tc1)
         * based on tile_bit, and within each 4-tile group, iterate local
         * coordinate half-pairs based on lo.
         *
         * The 6 swap pairs, with row/col split into (tile, local):
         *   Pair 1: T(tr0,tc0)[lr1,lc0] <-> T(tr1,tc0)[lr0,lc0]
         *   Pair 2: T(tr0,tc0)[lr1,lc1] <-> T(tr1,tc1)[lr0,lc0]
         *   Pair 3: T(tr1,tc0)[lr1,lc1] <-> T(tr1,tc1)[lr1,lc0]
         *   Pair 4: T(tr1,tc0)[lr0,lc1] <-> T(tr0,tc1)[lr1,lc0]
         *   Pair 5: T(tr0,tc0)[lr0,lc1] <-> T(tr0,tc1)[lr0,lc0]
         *   Pair 6: T(tr1,tc1)[lr0,lc1] <-> T(tr0,tc1)[lr1,lc1]
         *
         * Pairs 1-3: both tiles always in lower triangle.
         * Pairs 4-6: T(tr0,tc1) may be upper-tri, accessed via T(tc1,tr0)+conj.
         *
         * Diagonal tile-block note (is_diag, i.e., tr0==tc0, tr1==tc1):
         *   T(tr0,tr0) and T(tr1,tr1) are diagonal tiles where elem(lr,lc)
         *   and elem(lc,lr) are conjugates. Swaps modify some positions in
         *   these tiles. To avoid reading stale upper-triangle values, we
         *   use dtile_read/dtile_write which canonicalize through the
         *   lower-triangle position. After all swaps, we mirror the upper
         *   triangle of diagonal tiles to restore consistency.
         */
        const qubit_t tile_bit = hi - LOG_TILE_DIM;
        const gate_idx_t tile_stride = (gate_idx_t)1 << tile_bit;
        const gate_idx_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
        const gate_idx_t half_nt = n_tiles >> 1;
        const gate_idx_t tile_dim = (dim < TILE_DIM) ? dim : TILE_DIM;
        const gate_idx_t half_td = tile_dim >> 1;
        const gate_idx_t incr_lo = (gate_idx_t)1 << lo;

        #pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t i_tr = 0; i_tr < half_nt; ++i_tr) {
            gate_idx_t tr0 = insertBit0(i_tr, tile_bit);
            gate_idx_t tr1 = tr0 | tile_stride;

            for (gate_idx_t i_tc = 0; i_tc <= i_tr; ++i_tc) {
                gate_idx_t tc0 = insertBit0(i_tc, tile_bit);
                gate_idx_t tc1 = tc0 | tile_stride;

                /* Tile pointers: T00, T11, T10 always in lower triangle */
                cplx_t *t00 = data + tile_off(tr0, tc0);
                cplx_t *t11 = data + tile_off(tr1, tc1);
                cplx_t *t10 = data + tile_off(tr1, tc0);

                /* T01 = T(tr0, tc1): only stored if tr0 >= tc1 */
                int lower_01 = (tr0 >= tc1);
                cplx_t *t01 = lower_01 ? data + tile_off(tr0, tc1)
                                       : data + tile_off(tc1, tr0);

                int is_diag = (i_tc == i_tr);

                for (gate_idx_t b_lr = 0; b_lr < half_td; ++b_lr) {
                    gate_idx_t lr0 = insertBit0(b_lr, lo);
                    gate_idx_t lr1 = lr0 | incr_lo;

                    for (gate_idx_t b_lc = 0; b_lc < half_td; ++b_lc) {
                        gate_idx_t lc0 = insertBit0(b_lc, lo);
                        gate_idx_t lc1 = lc0 | incr_lo;

                        if (is_diag) {
                            /* t01 == t10 (same memory), and t00/t11 are diagonal
                             * tiles where dtile canonicalization creates further
                             * aliasing. Blocks (b_lr, b_lc) and (b_lc, b_lr)
                             * write identical values to the same physical
                             * locations, so process only b_lc <= b_lr. */
                            if (b_lc > b_lr) break;
                            cplx_t p1a = dtile_read(t00, lr1, lc0);
                            cplx_t p2a = dtile_read(t00, lr1, lc1);
                            cplx_t p5a = dtile_read(t00, lr0, lc1);
                            cplx_t p2b = dtile_read(t11, lr0, lc0);
                            cplx_t p3b = dtile_read(t11, lr1, lc0);
                            cplx_t p6a = dtile_read(t11, lr0, lc1);
                            cplx_t p1b = t10[elem_off(lr0, lc0)];
                            cplx_t p3a = t10[elem_off(lr1, lc1)];
                            cplx_t p4a = t10[elem_off(lr0, lc1)];
                            cplx_t p4b = conj(t10[elem_off(lc0, lr1)]);
                            cplx_t p5b = conj(t10[elem_off(lc0, lr0)]);
                            cplx_t p6b = conj(t10[elem_off(lc1, lr1)]);

                            /* Pair 1 */ dtile_write(t00, lr1, lc0, p1b);
                                         t10[elem_off(lr0, lc0)] = p1a;
                            /* Pair 2 */ dtile_write(t00, lr1, lc1, p2b);
                                         dtile_write(t11, lr0, lc0, p2a);
                            /* Pair 3 */ t10[elem_off(lr1, lc1)] = p3b;
                                         dtile_write(t11, lr1, lc0, p3a);
                            /* Pair 4 */ t10[elem_off(lr0, lc1)] = p4b;
                                         t10[elem_off(lc0, lr1)] = conj(p4a);
                            /* Pair 5 */ dtile_write(t00, lr0, lc1, p5b);
                                         t10[elem_off(lc0, lr0)] = conj(p5a);
                            /* Pair 6 */ dtile_write(t11, lr0, lc1, p6b);
                                         t10[elem_off(lc1, lr1)] = conj(p6a);
                        } else {
                            cplx_t tmp;

                            /* Pair 1: T(tr0,tc0)[lr1,lc0] <-> T(tr1,tc0)[lr0,lc0] */
                            tmp = t00[elem_off(lr1, lc0)];
                            t00[elem_off(lr1, lc0)] = t10[elem_off(lr0, lc0)];
                            t10[elem_off(lr0, lc0)] = tmp;

                            /* Pair 2: T(tr0,tc0)[lr1,lc1] <-> T(tr1,tc1)[lr0,lc0] */
                            tmp = t00[elem_off(lr1, lc1)];
                            t00[elem_off(lr1, lc1)] = t11[elem_off(lr0, lc0)];
                            t11[elem_off(lr0, lc0)] = tmp;

                            /* Pair 3: T(tr1,tc0)[lr1,lc1] <-> T(tr1,tc1)[lr1,lc0] */
                            tmp = t10[elem_off(lr1, lc1)];
                            t10[elem_off(lr1, lc1)] = t11[elem_off(lr1, lc0)];
                            t11[elem_off(lr1, lc0)] = tmp;

                            /* Pair 4: T(tr1,tc0)[lr0,lc1] <-> T(tr0,tc1)[lr1,lc0] */
                            {
                                cplx_t val_a = t10[elem_off(lr0, lc1)];
                                cplx_t val_b;
                                if (lower_01) {
                                    val_b = t01[elem_off(lr1, lc0)];
                                    t10[elem_off(lr0, lc1)] = val_b;
                                    t01[elem_off(lr1, lc0)] = val_a;
                                } else {
                                    val_b = conj(t01[elem_off(lc0, lr1)]);
                                    t10[elem_off(lr0, lc1)] = val_b;
                                    t01[elem_off(lc0, lr1)] = conj(val_a);
                                }
                            }

                            /* Pair 5: T(tr0,tc0)[lr0,lc1] <-> T(tr0,tc1)[lr0,lc0] */
                            {
                                cplx_t val_a = t00[elem_off(lr0, lc1)];
                                cplx_t val_b;
                                if (lower_01) {
                                    val_b = t01[elem_off(lr0, lc0)];
                                    t01[elem_off(lr0, lc0)] = val_a;
                                } else {
                                    val_b = conj(t01[elem_off(lc0, lr0)]);
                                    t01[elem_off(lc0, lr0)] = conj(val_a);
                                }
                                t00[elem_off(lr0, lc1)] = val_b;
                            }

                            /* Pair 6: T(tr1,tc1)[lr0,lc1] <-> T(tr0,tc1)[lr1,lc1] */
                            {
                                cplx_t val_a = t11[elem_off(lr0, lc1)];
                                cplx_t val_b;
                                if (lower_01) {
                                    val_b = t01[elem_off(lr1, lc1)];
                                    t01[elem_off(lr1, lc1)] = val_a;
                                } else {
                                    val_b = conj(t01[elem_off(lc1, lr1)]);
                                    t01[elem_off(lc1, lr1)] = conj(val_a);
                                }
                                t11[elem_off(lr0, lc1)] = val_b;
                            }
                        }
                    }
                }

                /* Mirror upper triangle of diagonal tiles from lower triangle.
                 * When is_diag, T(tr0,tr0) and T(tr1,tr1) are diagonal tiles
                 * whose lower-triangle elements are now correct. Propagate to
                 * upper triangle for storage consistency. */
                if (is_diag) {
                    for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                        for (gate_idx_t lc = lr + 1; lc < TILE_DIM; ++lc) {
                            t00[elem_off(lr, lc)] = conj(t00[elem_off(lc, lr)]);
                            t11[elem_off(lr, lc)] = conj(t11[elem_off(lc, lr)]);
                        }
                    }
                }
            }
        }
    }
    else {
        /*
         * Case 3: Both qubits across tiles.
         *
         * Both bits lo and hi are in tile coordinates. The permutation only
         * affects tile indices, not local coordinates within tiles. For each
         * local position (lr, lc), the swap operates on tile-level quartets.
         *
         * tile_lo = lo - LOG_TILE_DIM, tile_hi = hi - LOG_TILE_DIM.
         *
         * We iterate tile-row quartets and tile-col quartets (lower triangle
         * at the quartet base level). Since local coordinates are unaffected,
         * the swaps operate on entire tiles (TILE_SIZE elements at a time).
         *
         * The 6 swap pairs (at the tile level):
         *   Pair 1: T(tr01,tc00) <-> T(tr10,tc00)
         *   Pair 2: T(tr01,tc01) <-> T(tr10,tc10)
         *   Pair 3: T(tr11,tc01) <-> T(tr11,tc10)
         *   Pair 4: T(tr10,tc01) <-> T(tr01,tc10)
         *   Pair 5: T(tr00,tc01) <-> T(tr00,tc10)
         *   Pair 6: T(tr10,tc11) <-> T(tr01,tc11)
         *
         * Tile lower-triangle status:
         *   Pairs 1-3: both tiles always in lower triangle.
         *   Pair 4: T(tr01,tc10) may be upper-tri if tr01 < tc10.
         *   Pairs 5,6: may involve upper-tri tiles depending on whether
         *     tr00 >= tc01 and tr00 >= tc10.
         */
        const qubit_t tile_lo = lo - LOG_TILE_DIM;
        const qubit_t tile_hi = hi - LOG_TILE_DIM;
        const gate_idx_t incr_lo = (gate_idx_t)1 << tile_lo;
        const gate_idx_t incr_hi = (gate_idx_t)1 << tile_hi;
        const gate_idx_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
        const gate_idx_t quarter_nt = n_tiles >> 2;

        #pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t b_tr = 0; b_tr < quarter_nt; ++b_tr) {
            gate_idx_t tr00 = insertBits2_0(b_tr, tile_lo, tile_hi);
            gate_idx_t tr01 = tr00 | incr_lo;
            gate_idx_t tr10 = tr00 | incr_hi;
            gate_idx_t tr11 = tr10 | incr_lo;

            for (gate_idx_t b_tc = 0; b_tc <= b_tr; ++b_tc) {
                gate_idx_t tc00 = insertBits2_0(b_tc, tile_lo, tile_hi);
                gate_idx_t tc01 = tc00 | incr_lo;
                gate_idx_t tc10 = tc00 | incr_hi;
                gate_idx_t tc11 = tc10 | incr_lo;

                int is_diag = (b_tc == b_tr);

                /* ---- Pair 1: T(tr01,tc00) <-> T(tr10,tc00) ----
                 * Both always lower-tri (tr01 >= tr00 >= tc00, tr10 > tr00 >= tc00). */
                {
                    cplx_t *tA = data + tile_off(tr01, tc00);
                    cplx_t *tB = data + tile_off(tr10, tc00);
                    for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                        cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                    }
                }

                /* ---- Pair 2: T(tr01,tc01) <-> T(tr10,tc10) ----
                 * Both always lower-tri (differences equal tr00-tc00 >= 0). */
                {
                    cplx_t *tA = data + tile_off(tr01, tc01);
                    cplx_t *tB = data + tile_off(tr10, tc10);
                    for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                        cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                    }
                }

                /* ---- Pair 3: T(tr11,tc01) <-> T(tr11,tc10) ----
                 * Both always lower-tri. */
                {
                    cplx_t *tA = data + tile_off(tr11, tc01);
                    cplx_t *tB = data + tile_off(tr11, tc10);
                    for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                        cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                    }
                }

                /* ---- Pair 4: T(tr10,tc01) <-> T(tr01,tc10) ----
                 * T(tr10,tc01): always lower-tri.
                 * T(tr01,tc10): lower-tri iff tr01 >= tc10.
                 *
                 * When is_diag (tr00==tc00): tr01=tc01, tc10=tr10, so this becomes
                 * T(tr10,tr01) <-> T(tr01,tr10). Since tr10 > tr01, T(tr01,tr10) is
                 * upper-tri, stored as T(tr10,tr01) = tA. So tA == tB and the "swap"
                 * is the in-place Hermitian transpose of a single tile. */
                {
                    cplx_t *tA = data + tile_off(tr10, tc01);
                    int lower_b = (tr01 >= tc10);
                    if (lower_b) {
                        cplx_t *tB = data + tile_off(tr01, tc10);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                        }
                    } else if (is_diag) {
                        /* tA points to T(tr10, tr01). The logical swap maps
                         * rho[tr10+lr, tr01+lc] <-> rho[tr01+lr, tr10+lc].
                         * In storage: tA[lr,lc] <-> conj(tA[lc,lr]).
                         * Hermitian-transpose the tile in place. */
                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            tA[elem_off(lr, lr)] = conj(tA[elem_off(lr, lr)]);
                            for (gate_idx_t lc = lr + 1; lc < TILE_DIM; ++lc) {
                                cplx_t val_a = tA[elem_off(lr, lc)];
                                cplx_t val_b = tA[elem_off(lc, lr)];
                                tA[elem_off(lr, lc)] = conj(val_b);
                                tA[elem_off(lc, lr)] = conj(val_a);
                            }
                        }
                    } else {
                        /* tA = T(tr10,tc01), tB = T(tc10,tr01). Distinct tiles.
                         * Logical: tA[lr,lc] <-> conj(tB[lc,lr]). */
                        cplx_t *tB = data + tile_off(tc10, tr01);
                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                cplx_t val_a = tA[elem_off(lr, lc)];
                                cplx_t val_b = conj(tB[elem_off(lc, lr)]);
                                tA[elem_off(lr, lc)] = val_b;
                                tB[elem_off(lc, lr)] = conj(val_a);
                            }
                        }
                    }
                }

                /* ---- Pairs 5+6 ----
                 * Pair 5: T(tr00,tc01) <-> T(tr00,tc10)
                 * Pair 6: T(tr10,tc11) <-> T(tr01,tc11)
                 *
                 * Lower-tri conditions (using incr_hi > incr_lo):
                 *   T(tr00,tc01): iff tr00 >= tc01     [lo_01]
                 *   T(tr00,tc10): iff tr00 >= tc10     [lo_10]
                 *   T(tr10,tc11): iff tr10 >= tc11, equiv tr00 >= tc01  [lo_01]
                 *   T(tr01,tc11): iff tr01 >= tc11, equiv tr00 >= tc10  [lo_10]
                 *
                 * Three sub-cases for !is_diag:
                 * (a) lo_01 && lo_10: all 4 tiles lower-tri -> plain swap.
                 * (b) lo_01 && !lo_10: mixed -> conjugate access for 2 tiles.
                 * (c) !lo_01: all 4 upper-tri -> swap stored reps directly.
                 * For is_diag: tr00==tc00 implies tc01>tr00, so !lo_01 always.
                 * Handled separately after the if/else chain.
                 */
                {
                    int lo_01 = (tr00 >= tc01);
                    int lo_10 = (tr00 >= tc10);

                    if (lo_01 && lo_10) {
                        /* (a) All tiles in lower triangle: plain tile swaps */
                        cplx_t *t5a = data + tile_off(tr00, tc01);
                        cplx_t *t5b = data + tile_off(tr00, tc10);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = t5a[i]; t5a[i] = t5b[i]; t5b[i] = tmp;
                        }
                        cplx_t *t6a = data + tile_off(tr10, tc11);
                        cplx_t *t6b = data + tile_off(tr01, tc11);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = t6a[i]; t6a[i] = t6b[i]; t6b[i] = tmp;
                        }
                    }
                    else if (lo_01) {
                        /* (b) T(tr00,tc01) lower, T(tr00,tc10) upper -> T(tc10,tr00)
                         *     T(tr10,tc11) lower, T(tr01,tc11) upper -> T(tc11,tr01) */
                        cplx_t *t5a = data + tile_off(tr00, tc01);
                        cplx_t *t5b = data + tile_off(tc10, tr00);
                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                cplx_t val_a = t5a[elem_off(lr, lc)];
                                cplx_t val_b = conj(t5b[elem_off(lc, lr)]);
                                t5a[elem_off(lr, lc)] = val_b;
                                t5b[elem_off(lc, lr)] = conj(val_a);
                            }
                        }
                        cplx_t *t6a = data + tile_off(tr10, tc11);
                        cplx_t *t6b = data + tile_off(tc11, tr01);
                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                cplx_t val_a = t6a[elem_off(lr, lc)];
                                cplx_t val_b = conj(t6b[elem_off(lc, lr)]);
                                t6a[elem_off(lr, lc)] = val_b;
                                t6b[elem_off(lc, lr)] = conj(val_a);
                            }
                        }
                    }
                    else if (!is_diag) {
                        /* (c) All 4 tiles upper-tri. Swap their stored lower-tri
                         * representatives directly (conjugation on both sides
                         * cancels for a plain swap between transposed tiles). */
                        cplx_t *t5a = data + tile_off(tc01, tr00);
                        cplx_t *t5b = data + tile_off(tc10, tr00);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = t5a[i]; t5a[i] = t5b[i]; t5b[i] = tmp;
                        }
                        cplx_t *t6a = data + tile_off(tc11, tr10);
                        cplx_t *t6b = data + tile_off(tc11, tr01);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = t6a[i]; t6a[i] = t6b[i]; t6b[i] = tmp;
                        }
                    }
                    /* is_diag && !lo_01: When is_diag, tr==tc so Pairs 5+6
                     * operate on the same physical tiles as Pairs 1+3
                     * (the stored lower-tri representatives coincide).
                     * Executing them would undo the Pair 1+3 swaps. Skip. */
                }
            }
        }
    }
}

void ccx_tiled(state_t *state, const qubit_t ctrl1, const qubit_t ctrl2, const qubit_t target) {
    GATE_VALIDATE(0, "ccx: tiled not yet implemented");
}
