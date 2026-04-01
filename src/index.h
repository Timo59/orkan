/**
 * @file index.h
 * @brief Shared inline helpers for index computation across all backends
 *
 * Provides bit-insertion and storage-offset functions used by the gate, channel,
 * and state modules. All functions operate on idx_t (uint64_t) from q_types.h.
 */

#ifndef INDEX_H
#define INDEX_H

#include "q_types.h"

/*
 * =====================================================================================================================
 * Bit-insertion helpers
 * =====================================================================================================================
 */

/**
 * @brief Insert a 0 bit at position @p pos in value @p val
 *
 * For k from 0 to dim/2-1, insertBit0(k, target) produces all indices with bit
 * @p pos equal to 0. This enables direct enumeration without wasted iterations.
 *
 * Example: pos=2, val=0b101 -> 0b1001 (insert 0 at bit 2)
 */
static inline idx_t insertBit0(idx_t val, qubit_t pos) {
    idx_t mask = ((idx_t)1 << pos) - 1;
    return (val & mask) | ((val & ~mask) << 1);
}

/**
 * @brief Insert 0 bits at two positions in value @p val
 *
 * For two-qubit gates, we iterate over indices where both qubit bits are 0.
 * This helper inserts 0 bits at positions @p lo and @p hi (lo < hi).
 */
static inline idx_t insertBits2_0(idx_t val, qubit_t lo, qubit_t hi) {
    idx_t mask_lo = ((idx_t)1 << lo) - 1;
    val = (val & mask_lo) | ((val & ~mask_lo) << 1);
    idx_t mask_hi = ((idx_t)1 << hi) - 1;
    return (val & mask_hi) | ((val & ~mask_hi) << 1);
}

/*
 * =====================================================================================================================
 * Packed storage helpers
 * =====================================================================================================================
 */

/**
 * @brief Packed lower-triangle column-major index for element (r, c) where r >= c
 *
 * Formula: c * (2*dim - c + 1) / 2 + (r - c)
 */
static inline idx_t pack_idx(idx_t dim, idx_t r, idx_t c) {
    return c * (2 * dim - c + 1) / 2 + (r - c);
}

/**
 * @brief Column offset in packed lower-triangle column-major format
 *
 * Returns the index of the first element of column @p col.
 */
static inline idx_t col_off(idx_t dim, idx_t col) {
    return col * (2 * dim - col + 1) / 2;
}

/*
 * =====================================================================================================================
 * Tiled storage helpers
 * =====================================================================================================================
 */

/**
 * @brief Base offset of tile (tr, tc) in the data array (tr >= tc)
 */
static inline idx_t tile_off(idx_t tr, idx_t tc) {
    return (tr * (tr + 1) / 2 + tc) * TILE_SIZE;
}

/**
 * @brief Element offset within a tile (row-major)
 */
static inline idx_t elem_off(idx_t lr, idx_t lc) {
    return lr * TILE_DIM + lc;
}

#endif /* INDEX_H */
