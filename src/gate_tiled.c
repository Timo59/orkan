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

static inline dim_t tile_off(dim_t tr, dim_t tc) {
    return (tr * (tr + 1) / 2 + tc) * TILE_SIZE;
}

static inline dim_t elem_off(dim_t lr, dim_t lc) {
    return lr * TILE_DIM + lc;
}

/**
 * @brief Insert a 0 bit at position `pos` in value `val`
 */
static inline dim_t insertBit0(dim_t val, qubit_t pos) {
    dim_t mask = ((dim_t)1 << pos) - 1;
    return (val & mask) | ((val & ~mask) << 1);
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
 * Diagonal tiles (tr == tc): Only process blocks where lr0 >= lc0 (lower triangle
 * within the tile). Also skip the block at lr0==lc0 for the (0,1) element since
 * it is the conjugate of the (1,0) element.
 *
 * Off-diagonal tiles (tr > tc): Process all blocks (full tile stored).
 */
#define TRAVERSE_TILED_WITHIN(state, target, OFFDIAG_OP, DIAG_OP)               \
do {                                                                            \
    const dim_t dim = POW2((state)->qubits, dim_t);                             \
    const dim_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;                \
    const qubit_t tgt = (target);                                               \
    const dim_t stride = (dim_t)1 << tgt;                                       \
    const dim_t tile_dim = (dim < TILE_DIM) ? dim : TILE_DIM;                    \
    const dim_t half_td = tile_dim >> 1;                                        \
    cplx_t *data = (state)->data;                                               \
                                                                                \
    /* Iterate lower-triangular tiles */                                        \
    _Pragma("omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)")      \
    for (dim_t tr = 0; tr < n_tiles; ++tr) {                                    \
        for (dim_t tc = 0; tc <= tr; ++tc) {                                    \
            cplx_t *tile = data + tile_off(tr, tc);                             \
            int is_diag = (tr == tc);                                           \
                                                                                \
            /* Iterate blocks within the tile */                                \
            for (dim_t br = 0; br < half_td; ++br) {                            \
                dim_t lr0 = insertBit0(br, tgt);                                \
                dim_t lr1 = lr0 | stride;                                       \
                                                                                \
                for (dim_t bc = 0; bc < half_td; ++bc) {                        \
                    dim_t lc0 = insertBit0(bc, tgt);                            \
                    dim_t lc1 = lc0 | stride;                                   \
                                                                                \
                    if (is_diag) {                                               \
                        /* Only process blocks where lr0 >= lc0 */              \
                        if (lr0 < lc0) continue;                                \
                        /* On diagonal tile, a01 at (lr0,lc1) may be upper */   \
                        /* triangle. If lr0 < lc1, read conj from (lc1,lr0) */ \
                        int lower_01 = (lr0 >= lc1);                            \
                        cplx_t *a01_ptr = lower_01                              \
                            ? tile + elem_off(lr0, lc1)                         \
                            : tile + elem_off(lc1, lr0);                        \
                        int diag_block = (lr0 == lc0);                          \
                        cplx_t *a00 = tile + elem_off(lr0, lc0);               \
                        cplx_t *a10 = tile + elem_off(lr1, lc0);               \
                        cplx_t *a11 = tile + elem_off(lr1, lc1);               \
                        DIAG_OP(a00, a01_ptr, a10, a11,                         \
                                lower_01, diag_block);                          \
                    } else {                                                     \
                        /* Off-diagonal tile: full tile stored, no conj */       \
                        cplx_t *a00 = tile + elem_off(lr0, lc0);               \
                        cplx_t *a01 = tile + elem_off(lr0, lc1);               \
                        cplx_t *a10 = tile + elem_off(lr1, lc0);               \
                        cplx_t *a11 = tile + elem_off(lr1, lc1);               \
                        OFFDIAG_OP(a00, a01, a10, a11);                         \
                    }                                                            \
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
#define TRAVERSE_TILED_CROSS(state, target, BLOCK_OP)                           \
do {                                                                            \
    const dim_t dim = POW2((state)->qubits, dim_t);                             \
    const dim_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;                \
    const qubit_t tile_bit = (target) - LOG_TILE_DIM;                           \
    const dim_t tile_stride = (dim_t)1 << tile_bit;                             \
    const dim_t half_nt = n_tiles >> 1;                                         \
    cplx_t *data = (state)->data;                                               \
                                                                                \
    /* Iterate tile-column pairs (tc0, tc1) with tile_bit=0 */                  \
    _Pragma("omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)")      \
    for (dim_t i_tc = 0; i_tc < half_nt; ++i_tc) {                             \
        dim_t tc0 = insertBit0(i_tc, tile_bit);                                \
        dim_t tc1 = tc0 | tile_stride;                                         \
                                                                                \
        /* Iterate tile-row pairs (tr0, tr1) with tile_bit=0 */                \
        /* Start from i_tc so that tr0 >= tc0 */                               \
        for (dim_t i_tr = i_tc; i_tr < half_nt; ++i_tr) {                     \
            dim_t tr0 = insertBit0(i_tr, tile_bit);                            \
            dim_t tr1 = tr0 | tile_stride;                                     \
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
            /* Iterate all (lr, lc) within tiles */                            \
            for (dim_t lr = 0; lr < TILE_DIM; ++lr) {                          \
                dim_t lc_start = is_diag ? lr : 0;                             \
                for (dim_t lc = lc_start; lc < TILE_DIM; ++lc) {              \
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
 * These macros receive pointers to the 4 elements and a diag_block flag.
 * On diagonal blocks, a01 and a10 point to the same conjugate pair.
 * Since tiles store full (non-packed) data, a01 and a10 are independent
 * EXCEPT on diagonal blocks of diagonal tiles where a01==conj(a10) by
 * Hermitian symmetry. But since we store both explicitly in the tile,
 * we operate on both and the result maintains Hermitian structure.
 *
 * For diagonal blocks: after the operation, we enforce Hermiticity by
 * setting a01 = conj(a10).
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
    *(a01) = CMPLX(cimag(*(a01)), -creal(*(a01)));  \
    *(a10) = CMPLX(-cimag(*(a10)), creal(*(a10)));  \
} while(0)

#define SDG_OFFDIAG_OP(a00, a01, a10, a11) do {\
    *(a01) = CMPLX(-cimag(*(a01)), creal(*(a01)));  \
    *(a10) = CMPLX(cimag(*(a10)), -creal(*(a10)));  \
} while(0)

#define T_OFFDIAG_OP(a00, a01, a10, a11) do { \
    double r01 = creal(*(a01)), i01 = cimag(*(a01)); \
    double r10 = creal(*(a10)), i10 = cimag(*(a10)); \
    *(a01) = CMPLX((r01 + i01) * M_SQRT1_2,          \
                   (i01 - r01) * M_SQRT1_2);          \
    *(a10) = CMPLX((r10 - i10) * M_SQRT1_2,          \
                   (r10 + i10) * M_SQRT1_2);          \
} while(0)

#define TDG_OFFDIAG_OP(a00, a01, a10, a11) do {\
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
    cplx_t new01 = 0.5 * (v00 - v01 + v10 - v11);                   \
    *a10 = new10;                                                    \
    *a01 = lower_01 ? new01 : conj(new01);                           \
    if (diag_block) {                                                \
        *a01 = lower_01 ? conj(*a10) : *a10;                        \
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
    cplx_t new01 = 0.5 * (v00 + v01 - v10 - v11);                   \
    *a10 = new10;                                                    \
    *a01 = lower_01 ? new01 : conj(new01);                           \
    if (diag_block) {                                                \
        *a01 = lower_01 ? conj(*a10) : *a10;                        \
    }                                                                \
} while(0)

/*
 * =====================================================================================================================
 * Dispatch: within-tile vs cross-tile
 * =====================================================================================================================
 */

#define DISPATCH_TILED_1Q(state, target, OFFDIAG_OP, CROSS_OP)      \
do {                                                                \
    if ((target) < LOG_TILE_DIM) {                                  \
        TRAVERSE_TILED_WITHIN(state, target, OFFDIAG_OP, CROSS_OP);\
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
 * Rotation gates
 * =====================================================================================================================
 */

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
void rx_tiled(state_t *state, const qubit_t target, const double theta) {
    const double co = cos(theta / 2.0);
    const double si = sin(theta / 2.0);
    const double c2 = co * co, s2 = si * si, cs = co * si;

    if (target < LOG_TILE_DIM) {
        /* Within-tile traversal */
        const dim_t dim = POW2(state->qubits, dim_t);
        const dim_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
        const dim_t stride = (dim_t)1 << target;
        const dim_t tile_dim = (dim < TILE_DIM) ? dim : TILE_DIM;
        const dim_t half_td = tile_dim >> 1;
        cplx_t *data = state->data;

        #pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
        for (dim_t tr = 0; tr < n_tiles; ++tr) {
            for (dim_t tc = 0; tc <= tr; ++tc) {
                cplx_t *tile = data + tile_off(tr, tc);
                int is_diag = (tr == tc);

                for (dim_t br = 0; br < half_td; ++br) {
                    dim_t lr0 = insertBit0(br, target);
                    dim_t lr1 = lr0 | stride;
                    for (dim_t bc = 0; bc < half_td; ++bc) {
                        dim_t lc0 = insertBit0(bc, target);
                        dim_t lc1 = lc0 | stride;

                        if (is_diag && lr0 < lc0) continue;
                        int diag_block = is_diag && (lr0 == lc0);

                        cplx_t v00, v01, v10, v11;
                        v00 = tile[elem_off(lr0, lc0)];
                        v11 = tile[elem_off(lr1, lc1)];
                        v10 = tile[elem_off(lr1, lc0)];

                        if (is_diag) {
                            int lower_01 = (lr0 >= lc1);
                            v01 = lower_01 ? tile[elem_off(lr0, lc1)]
                                           : conj(tile[elem_off(lc1, lr0)]);
                        } else {
                            v01 = tile[elem_off(lr0, lc1)];
                        }

                        cplx_t diff_01_10 = v01 - v10;
                        cplx_t ics_diff = CMPLX(-cs * cimag(diff_01_10), cs * creal(diff_01_10));
                        cplx_t diff_diag = v00 - v11;
                        cplx_t ics_ddiag = CMPLX(-cs * cimag(diff_diag), cs * creal(diff_diag));

                        cplx_t n00 = c2 * v00 + s2 * v11 + ics_diff;
                        cplx_t n11 = s2 * v00 + c2 * v11 - ics_diff;
                        cplx_t n10 = -ics_ddiag + c2 * v10 + s2 * v01;
                        cplx_t n01 = ics_ddiag + c2 * v01 + s2 * v10;

                        tile[elem_off(lr0, lc0)] = n00;
                        tile[elem_off(lr1, lc1)] = n11;
                        tile[elem_off(lr1, lc0)] = n10;

                        if (is_diag) {
                            int lower_01 = (lr0 >= lc1);
                            if (diag_block) {
                                if (lower_01)
                                    tile[elem_off(lr0, lc1)] = conj(n10);
                                else
                                    tile[elem_off(lc1, lr0)] = n10;
                            } else {
                                if (lower_01)
                                    tile[elem_off(lr0, lc1)] = n01;
                                else
                                    tile[elem_off(lc1, lr0)] = conj(n01);
                            }
                        } else {
                            tile[elem_off(lr0, lc1)] = n01;
                        }
                    }
                }
            }
        }
    } else {
        /* Cross-tile traversal */
        const dim_t dim = POW2(state->qubits, dim_t);
        const dim_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
        const qubit_t tile_bit = target - LOG_TILE_DIM;
        const dim_t tile_stride = (dim_t)1 << tile_bit;
        const dim_t half_nt = n_tiles >> 1;
        cplx_t *data = state->data;

        #pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
        for (dim_t i_tc = 0; i_tc < half_nt; ++i_tc) {
            dim_t tc0 = insertBit0(i_tc, tile_bit);
            dim_t tc1 = tc0 | tile_stride;
            for (dim_t i_tr = i_tc; i_tr < half_nt; ++i_tr) {
                dim_t tr0 = insertBit0(i_tr, tile_bit);
                dim_t tr1 = tr0 | tile_stride;

                cplx_t *t00 = data + tile_off(tr0, tc0);
                cplx_t *t11 = data + tile_off(tr1, tc1);
                cplx_t *t10 = data + tile_off(tr1, tc0);

                int lower_01 = (tr0 >= tc1);
                cplx_t *t01 = lower_01 ? data + tile_off(tr0, tc1)
                                       : data + tile_off(tc1, tr0);
                int is_diag = (i_tc == i_tr);

                for (dim_t lr = 0; lr < TILE_DIM; ++lr) {
                    dim_t lc_start = is_diag ? lr : 0;
                    for (dim_t lc = lc_start; lc < TILE_DIM; ++lc) {
                        int diag_block = is_diag && (lr == lc);

                        cplx_t v00 = t00[elem_off(lr, lc)];
                        cplx_t v11 = t11[elem_off(lr, lc)];
                        cplx_t v10 = t10[elem_off(lr, lc)];
                        cplx_t rho01 = lower_01 ? t01[elem_off(lr, lc)]
                                                : conj(t01[elem_off(lc, lr)]);

                        cplx_t diff_01_10 = rho01 - v10;
                        cplx_t ics_diff = CMPLX(-cs * cimag(diff_01_10), cs * creal(diff_01_10));
                        cplx_t diff_diag = v00 - v11;
                        cplx_t ics_ddiag = CMPLX(-cs * cimag(diff_diag), cs * creal(diff_diag));

                        cplx_t n00 = c2 * v00 + s2 * v11 + ics_diff;
                        cplx_t n11 = s2 * v00 + c2 * v11 - ics_diff;
                        cplx_t n10 = -ics_ddiag + c2 * v10 + s2 * rho01;
                        cplx_t n01 = ics_ddiag + c2 * rho01 + s2 * v10;

                        t00[elem_off(lr, lc)] = n00;
                        t11[elem_off(lr, lc)] = n11;
                        t10[elem_off(lr, lc)] = n10;
                        if (diag_block) {
                            if (lower_01)
                                t01[elem_off(lr, lc)] = conj(n10);
                            else
                                t01[elem_off(lc, lr)] = n10;
                        } else {
                            if (lower_01)
                                t01[elem_off(lr, lc)] = n01;
                            else
                                t01[elem_off(lc, lr)] = conj(n01);
                        }
                    }
                }
            }
        }
    }
}

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
void ry_tiled(state_t *state, const qubit_t target, const double theta) {
    const double co = cos(theta / 2.0);
    const double si = sin(theta / 2.0);
    const double c2 = co * co, s2 = si * si, cs_val = co * si;

    if (target < LOG_TILE_DIM) {
        const dim_t dim = POW2(state->qubits, dim_t);
        const dim_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
        const dim_t stride = (dim_t)1 << target;
        const dim_t tile_dim = (dim < TILE_DIM) ? dim : TILE_DIM;
        const dim_t half_td = tile_dim >> 1;
        cplx_t *data = state->data;

        #pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
        for (dim_t tr = 0; tr < n_tiles; ++tr) {
            for (dim_t tc = 0; tc <= tr; ++tc) {
                cplx_t *tile = data + tile_off(tr, tc);
                int is_diag = (tr == tc);

                for (dim_t br = 0; br < half_td; ++br) {
                    dim_t lr0 = insertBit0(br, target);
                    dim_t lr1 = lr0 | stride;
                    for (dim_t bc = 0; bc < half_td; ++bc) {
                        dim_t lc0 = insertBit0(bc, target);
                        dim_t lc1 = lc0 | stride;

                        if (is_diag && lr0 < lc0) continue;
                        int diag_block = is_diag && (lr0 == lc0);

                        cplx_t v00, v01, v10, v11;
                        v00 = tile[elem_off(lr0, lc0)];
                        v11 = tile[elem_off(lr1, lc1)];
                        v10 = tile[elem_off(lr1, lc0)];

                        if (is_diag) {
                            int lower_01 = (lr0 >= lc1);
                            v01 = lower_01 ? tile[elem_off(lr0, lc1)]
                                           : conj(tile[elem_off(lc1, lr0)]);
                        } else {
                            v01 = tile[elem_off(lr0, lc1)];
                        }

                        cplx_t sum_01_10 = v01 + v10;
                        cplx_t diff_diag = v00 - v11;

                        cplx_t n00 = c2 * v00 + s2 * v11 - cs_val * sum_01_10;
                        cplx_t n11 = s2 * v00 + c2 * v11 + cs_val * sum_01_10;
                        cplx_t n10 = cs_val * diff_diag + c2 * v10 - s2 * v01;
                        cplx_t n01 = cs_val * diff_diag + c2 * v01 - s2 * v10;

                        tile[elem_off(lr0, lc0)] = n00;
                        tile[elem_off(lr1, lc1)] = n11;
                        tile[elem_off(lr1, lc0)] = n10;

                        if (is_diag) {
                            int lower_01 = (lr0 >= lc1);
                            if (diag_block) {
                                if (lower_01)
                                    tile[elem_off(lr0, lc1)] = conj(n10);
                                else
                                    tile[elem_off(lc1, lr0)] = n10;
                            } else {
                                if (lower_01)
                                    tile[elem_off(lr0, lc1)] = n01;
                                else
                                    tile[elem_off(lc1, lr0)] = conj(n01);
                            }
                        } else {
                            tile[elem_off(lr0, lc1)] = n01;
                        }
                    }
                }
            }
        }
    } else {
        const dim_t dim = POW2(state->qubits, dim_t);
        const dim_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
        const qubit_t tile_bit = target - LOG_TILE_DIM;
        const dim_t tile_stride = (dim_t)1 << tile_bit;
        const dim_t half_nt = n_tiles >> 1;
        cplx_t *data = state->data;

        #pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
        for (dim_t i_tc = 0; i_tc < half_nt; ++i_tc) {
            dim_t tc0 = insertBit0(i_tc, tile_bit);
            dim_t tc1 = tc0 | tile_stride;
            for (dim_t i_tr = i_tc; i_tr < half_nt; ++i_tr) {
                dim_t tr0 = insertBit0(i_tr, tile_bit);
                dim_t tr1 = tr0 | tile_stride;

                cplx_t *t00 = data + tile_off(tr0, tc0);
                cplx_t *t11 = data + tile_off(tr1, tc1);
                cplx_t *t10 = data + tile_off(tr1, tc0);

                int lower_01 = (tr0 >= tc1);
                cplx_t *t01 = lower_01 ? data + tile_off(tr0, tc1)
                                       : data + tile_off(tc1, tr0);
                int is_diag = (i_tc == i_tr);

                for (dim_t lr = 0; lr < TILE_DIM; ++lr) {
                    dim_t lc_start = is_diag ? lr : 0;
                    for (dim_t lc = lc_start; lc < TILE_DIM; ++lc) {
                        int diag_block = is_diag && (lr == lc);

                        cplx_t v00 = t00[elem_off(lr, lc)];
                        cplx_t v11 = t11[elem_off(lr, lc)];
                        cplx_t v10 = t10[elem_off(lr, lc)];
                        cplx_t rho01 = lower_01 ? t01[elem_off(lr, lc)]
                                                : conj(t01[elem_off(lc, lr)]);

                        cplx_t sum_01_10 = rho01 + v10;
                        cplx_t diff_diag = v00 - v11;

                        cplx_t n00 = c2 * v00 + s2 * v11 - cs_val * sum_01_10;
                        cplx_t n11 = s2 * v00 + c2 * v11 + cs_val * sum_01_10;
                        cplx_t n10 = cs_val * diff_diag + c2 * v10 - s2 * rho01;
                        cplx_t n01 = cs_val * diff_diag + c2 * rho01 - s2 * v10;

                        t00[elem_off(lr, lc)] = n00;
                        t11[elem_off(lr, lc)] = n11;
                        t10[elem_off(lr, lc)] = n10;
                        if (diag_block) {
                            if (lower_01)
                                t01[elem_off(lr, lc)] = conj(n10);
                            else
                                t01[elem_off(lc, lr)] = n10;
                        } else {
                            if (lower_01)
                                t01[elem_off(lr, lc)] = n01;
                            else
                                t01[elem_off(lc, lr)] = conj(n01);
                        }
                    }
                }
            }
        }
    }
}

/*
 * =====================================================================================================================
 * Rz gate (tiled): ρ → Rz(θ) ρ Rz(θ)†
 * =====================================================================================================================
 *
 * Rz(θ) = diag(e^(-iθ/2), e^(iθ/2)). Diagonal gate:
 *   ρ'_00 = ρ00, ρ'_11 = ρ11
 *   ρ'_10 = e^(iθ) ρ10, ρ'_01 = e^(-iθ) ρ01
 */
void rz_tiled(state_t *state, const qubit_t target, const double theta) {
    const double c = cos(theta), sn = sin(theta);

    if (target < LOG_TILE_DIM) {
        const dim_t dim = POW2(state->qubits, dim_t);
        const dim_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
        const dim_t stride = (dim_t)1 << target;
        const dim_t tile_dim = (dim < TILE_DIM) ? dim : TILE_DIM;
        const dim_t half_td = tile_dim >> 1;
        cplx_t *data = state->data;

        #pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
        for (dim_t tr = 0; tr < n_tiles; ++tr) {
            for (dim_t tc = 0; tc <= tr; ++tc) {
                cplx_t *tile = data + tile_off(tr, tc);
                int is_diag = (tr == tc);

                for (dim_t br = 0; br < half_td; ++br) {
                    dim_t lr0 = insertBit0(br, target);
                    dim_t lr1 = lr0 | stride;
                    for (dim_t bc = 0; bc < half_td; ++bc) {
                        dim_t lc0 = insertBit0(bc, target);
                        dim_t lc1 = lc0 | stride;

                        if (is_diag && lr0 < lc0) continue;
                        int diag_block = is_diag && (lr0 == lc0);

                        /* rho10 *= e^(iθ) */
                        cplx_t *p10 = tile + elem_off(lr1, lc0);
                        double r10 = creal(*p10), i10 = cimag(*p10);
                        *p10 = CMPLX(r10 * c - i10 * sn, i10 * c + r10 * sn);

                        if (diag_block) continue;

                        /* rho01 *= e^(-iθ) */
                        if (is_diag) {
                            int lower_01 = (lr0 >= lc1);
                            cplx_t *p01 = lower_01 ? tile + elem_off(lr0, lc1)
                                                   : tile + elem_off(lc1, lr0);
                            double r01 = creal(*p01), i01 = cimag(*p01);
                            if (lower_01) {
                                *p01 = CMPLX(r01 * c + i01 * sn, i01 * c - r01 * sn);
                            } else {
                                /* stored = conj(rho01); new_stored = conj(e^(-iθ) * rho01)
                                 * = e^(iθ) * conj(rho01) = e^(iθ) * stored */
                                *p01 = CMPLX(r01 * c - i01 * sn, i01 * c + r01 * sn);
                            }
                        } else {
                            cplx_t *p01 = tile + elem_off(lr0, lc1);
                            double r01 = creal(*p01), i01 = cimag(*p01);
                            *p01 = CMPLX(r01 * c + i01 * sn, i01 * c - r01 * sn);
                        }
                    }
                }
            }
        }
    } else {
        const dim_t dim = POW2(state->qubits, dim_t);
        const dim_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
        const qubit_t tile_bit = target - LOG_TILE_DIM;
        const dim_t tile_stride = (dim_t)1 << tile_bit;
        const dim_t half_nt = n_tiles >> 1;
        cplx_t *data = state->data;

        #pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
        for (dim_t i_tc = 0; i_tc < half_nt; ++i_tc) {
            dim_t tc0 = insertBit0(i_tc, tile_bit);
            dim_t tc1 = tc0 | tile_stride;
            for (dim_t i_tr = i_tc; i_tr < half_nt; ++i_tr) {
                dim_t tr0 = insertBit0(i_tr, tile_bit);
                dim_t tr1 = tr0 | tile_stride;

                cplx_t *t10 = data + tile_off(tr1, tc0);

                int lower_01 = (tr0 >= tc1);
                cplx_t *t01 = lower_01 ? data + tile_off(tr0, tc1)
                                       : data + tile_off(tc1, tr0);
                int is_diag = (i_tc == i_tr);

                for (dim_t lr = 0; lr < TILE_DIM; ++lr) {
                    dim_t lc_start = is_diag ? lr : 0;
                    for (dim_t lc = lc_start; lc < TILE_DIM; ++lc) {
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
                            /* stored = conj(rho01); new = conj(e^(-iθ)*rho01) = e^(iθ)*stored */
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
 * Two-qubit gates (stubs)
 * =====================================================================================================================
 */

void cx_tiled(state_t *state, const qubit_t control, const qubit_t target) {
    GATE_VALIDATE(0, "cx: tiled not yet implemented");
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

void swap_tiled(state_t *state, const qubit_t q1, const qubit_t q2) {
    GATE_VALIDATE(0, "swap: tiled not yet implemented");
}

void ccx_tiled(state_t *state, const qubit_t ctrl1, const qubit_t ctrl2, const qubit_t target) {
    GATE_VALIDATE(0, "ccx: tiled not yet implemented");
}
