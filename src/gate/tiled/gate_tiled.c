/**
 * @file gate_tiled.c
 * @brief Single-qubit and rotation gates for tiled density matrices
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
 *
 * Two-qubit gates are in separate files: gate_tiled_cx.c, gate_tiled_swap.c.
 */

#include "gate_tiled.h"

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
//    DISPATCH_TILED_1Q(state, target, X_OFFDIAG_OP, X_CROSS_OP);
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    const gate_idx_t n_tiles = (dim + TILE_DIM - (gate_idx_t)1) >> LOG_TILE_DIM;
    const gate_idx_t dim_tile = dim < TILE_DIM ? dim : TILE_DIM;
    const gate_idx_t incr = (gate_idx_t)1 << target;
    cplx_t * restrict data = state->data;

    if (target < LOG_TILE_DIM) {
        const gate_idx_t n_base = dim_tile >> 1;

        #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t tr = 0; tr < n_tiles; ++tr) {
            for (gate_idx_t tc = 0; tc <= tr; ++tc) {
                gate_idx_t offset_tile = (tr * (tr + 1) / 2 + tc) * TILE_SIZE;

                for (gate_idx_t br = 0; br < n_base; ++br) {
                    const gate_idx_t r0 = insertBit0(br, target);
                    const gate_idx_t r1 = r0 | incr;

                    const gate_idx_t offset_r0 = r0 * TILE_DIM + offset_tile;
                    const gate_idx_t offset_r1 = r1 * TILE_DIM + offset_tile;

                    for (gate_idx_t bc = 0; bc < n_base; ++bc) {
                        const gate_idx_t c0 = insertBit0(bc, target);
                        const gate_idx_t c1 = c0 | incr;

                        cplx_t tmp = data[c0 + offset_r0];
                        data[c0 + offset_r0] = data[c1 + offset_r1];
                        data[c1 + offset_r1] = tmp;

                        tmp = data[c0 + offset_r1];
                        data[c0 + offset_r1] = data[c1 + offset_r0];
                        data[c1 + offset_r0] = tmp;
                    }
                }
            }
        }  
    }

    else {
        const gate_idx_t n_base = n_tiles >> 1;
        const qubit_t target_tile = target - LOG_TILE_DIM;
        const gate_idx_t incr_tile = (gate_idx_t)1 << target_tile;

        #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t btr = 0; btr < n_base; ++btr) {
            const gate_idx_t tr0 = insertBit0(btr, target_tile);
            const gate_idx_t tr1 = tr0 | incr_tile;

            for (gate_idx_t btc = 0; btc <= btr; ++btc) {
                const gate_idx_t tc0 = insertBit0(btc, target_tile);
                const gate_idx_t tc1 = tc0 | incr_tile;

                const gate_idx_t offset_t00 = (tr0 * (tr0 + 1) / 2 + tc0) * TILE_SIZE;
                const gate_idx_t offset_t10 = (tr1 * (tr1 + 1) / 2 + tc0) * TILE_SIZE;
                const gate_idx_t offset_t11 = (tr1 * (tr1 + 1) / 2 + tc1) * TILE_SIZE;

                if (tr0 > tc1) {
                    const gate_idx_t offset_t01 = (tr0 * (tr0 + 1) / 2 + tc1) * TILE_SIZE;

                    for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                        const gate_idx_t idx00 = lr * TILE_DIM + offset_t00; 
                        const gate_idx_t idx01 = lr * TILE_DIM + offset_t01;
                        const gate_idx_t idx10 = lr * TILE_DIM + offset_t10;
                        const gate_idx_t idx11 = lr * TILE_DIM + offset_t11;

                        for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                            cplx_t tmp = data[idx00 + lc];
                            data[idx00 + lc] = data[idx11 + lc];
                            data[idx11 + lc] = tmp;

                            tmp = data[idx01 + lc];
                            data[idx01 + lc] = data[idx10 + lc];
                            data[idx10 + lc] = tmp;
                        }
                    }
                }
                else if (btr != btc) {
                    const gate_idx_t offset_t01 = (tc1 * (tc1 + 1) / 2 + tr0) * TILE_SIZE;

                    for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                        const gate_idx_t idx00 = lr * TILE_DIM + offset_t00; 
                        const gate_idx_t idx10 = lr * TILE_DIM + offset_t10;
                        const gate_idx_t idx11 = lr * TILE_DIM + offset_t11;

                        for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                            cplx_t tmp = data[idx00 + lc];
                            data[idx00 + lc] = data[idx11 + lc];
                            data[idx11 + lc] = tmp;

                            tmp = data[lr + lc * TILE_DIM + offset_t01];
                            data[lr + lc * TILE_DIM + offset_t01] = conj(data[idx10 + lc]);
                            data[idx10 + lc] = conj(tmp);
                        }
                    }
                }
                else {

                    for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                        const gate_idx_t idx00 = lr * TILE_DIM + offset_t00; 
                        const gate_idx_t idx10 = lr * TILE_DIM + offset_t10;
                        const gate_idx_t idx11 = lr * TILE_DIM + offset_t11;

                        for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                            cplx_t tmp = data[idx00 + lc];
                            data[idx00 + lc] = data[idx11 + lc];
                            data[idx11 + lc] = tmp;
                        }

                        for (gate_idx_t lc = 0; lc <= lr; ++lc) {
                            cplx_t tmp = data[idx10 + lc];
                            data[idx10 + lc] = conj(data[lr + lc * TILE_DIM + offset_t10]);
                            data[lr + lc * TILE_DIM + offset_t10] = conj(tmp);
                        }
                    }
                }
            }
        }
    }
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

/* Two-qubit gates: see gate_tiled_cx.c, gate_tiled_swap.c */
