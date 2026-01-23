# Implementation Plan: mhipster_block.c - Blocked Memory Layout for Mixed State Gates

## Overview

Implement a blocked memory layout alternative to the packed lower-triangular storage in `mhipster.c`, optimized for cache locality, SIMD vectorization, and better parallel decomposition.

## Design Decision: Tile-Based Blocking

**Chosen approach**: 2D tile blocking of the density matrix, storing tiles in row-major order within each tile, tiles ordered to preserve lower-triangular structure.

**Rationale**:
- Enables BLAS-3 operations (zgemm) for Hadamard-type gates
- Fixed stride within tiles enables SIMD vectorization
- Independent tiles enable embarrassingly parallel processing
- Cache-friendly: 64x64 complex tile = 64KB fits in L2

**Block size**: 64 complex elements per dimension (configurable via `BLOCK_DIM`)
- 64x64 x 16 bytes = 64KB per full tile
- Fits in L2 cache (128KB on M4), leaves room for working data
- Power-of-2 for efficient index calculations

## File Structure

```
src/
├── mhipster.c          (UNCHANGED)
├── state.c             (UNCHANGED)
├── state_blocked.c     (NEW - blocked state management)
├── mhipster_block.c    (NEW - blocked gate implementations)

include/
├── state.h             (UNCHANGED)
├── state_blocked.h     (NEW - blocked state type + functions)
├── mhipster_block.h    (NEW - blocked gate declarations)

test/
├── src/
│   ├── test_mhipster.c       (UNCHANGED)
│   ├── test_mhipster_block.c (NEW - blocked gate tests)
│   └── test_helpers.c        (NEW - shared helpers extracted from test_mhipster.c)
├── include/
│   └── test_helpers.h        (NEW - multiQubitGate declaration)

benchmark/
├── src/
│   └── bench_mixed.c   (MODIFY - only file modified, add blocked variant)
```

**Key design decisions**:
1. **Completely parallel infrastructure** - No modifications to existing files except bench_mixed.c
2. **Independent format** - States start blocked and stay blocked (conversion helpers provided for debugging)
3. **Extract test helpers** - Move `multiQubitGate()` to separate `test_helpers.c` to share between test files (avoids duplicate main() symbols)

## Phase 1: Foundation (Test-Driven)

### 1.1 Blocked State Type (Parallel to Existing)

**File**: `include/state_blocked.h` (NEW - does not modify state.h)

```c
#ifndef STATE_BLOCKED_H
#define STATE_BLOCKED_H

#include "q_types.h"

#define BLOCK_DIM 64  // Elements per block dimension (64x64 = 64KB per tile)

// Blocked state structure (parallel to state_t)
typedef struct state_blocked {
    cplx_t *data;       // Blocked storage
    qubit_t qubits;     // Number of qubits
    dim_t n_blocks;     // Number of tile rows/cols
    dim_t padded_dim;   // Dimension padded to BLOCK_DIM multiple
} state_blocked_t;

// State management
void state_blocked_init(state_blocked_t *state, qubit_t qubits);
void state_blocked_free(state_blocked_t *state);
void state_blocked_plus(state_blocked_t *state, qubit_t qubits);
dim_t state_blocked_len(qubit_t qubits);

// Index helpers
dim_t tile_offset(dim_t tile_row, dim_t tile_col, dim_t n_blocks);
cplx_t* get_tile(state_blocked_t *state, dim_t tile_row, dim_t tile_col);

// Element access (for debugging and measurement)
cplx_t state_blocked_get(const state_blocked_t *state, dim_t row, dim_t col);
void state_blocked_set(state_blocked_t *state, dim_t row, dim_t col, cplx_t val);

// Measurement support
double state_blocked_prob(const state_blocked_t *state, dim_t basis_state);
double state_blocked_trace(const state_blocked_t *state);

#endif
```

### 1.2 Blocked Gate Declarations

**File**: `include/mhipster_block.h` (NEW)

```c
#ifndef MHIPSTER_BLOCK_H
#define MHIPSTER_BLOCK_H

#include "state_blocked.h"

// Blocked gate declarations (separate API from packed)
void x_blocked(state_blocked_t *state, qubit_t target);
void y_blocked(state_blocked_t *state, qubit_t target);
void z_blocked(state_blocked_t *state, qubit_t target);
void h_blocked(state_blocked_t *state, qubit_t target);
void s_blocked(state_blocked_t *state, qubit_t target);
void sdg_blocked(state_blocked_t *state, qubit_t target);
void t_blocked(state_blocked_t *state, qubit_t target);
void tdg_blocked(state_blocked_t *state, qubit_t target);

#endif
```

### 1.3 Blocked Memory Layout

**Hermitian Symmetry**: We exploit Hermitian symmetry at the **tile level**:
- Only store tiles from the lower triangle: tile(i,j) where i >= j
- For n_blocks tiles per dimension, store n_blocks(n_blocks+1)/2 tiles (not n_blocks^2)
- Upper triangle elements accessed as conjugate of lower triangle counterpart

**Layout**: Lower-triangular tiles stored contiguously, each tile in row-major order.

```
For 4x4 blocks (N=4 tiles):
Tile(0,0) |
Tile(1,0) | Tile(1,1) |
Tile(2,0) | Tile(2,1) | Tile(2,2) |
Tile(3,0) | Tile(3,1) | Tile(3,2) | Tile(3,3) |

Storage: [Tile(0,0)][Tile(1,0)][Tile(1,1)][Tile(2,0)]...

Within each tile: FULL BLOCK_DIM x BLOCK_DIM, row-major for SIMD-friendly operations
```

**Memory comparison**:
- Packed format: N(N+1)/2 elements (optimal for Hermitian)
- Blocked format: n_tiles x BLOCK_DIM^2 where n_tiles = n_blocks(n_blocks+1)/2
- Trade-off: Blocked uses more memory (due to full tiles + padding) but enables SIMD

**Index functions**:
```c
// Tile index in storage (lower triangular)
static inline dim_t tile_offset(dim_t tile_row, dim_t tile_col, dim_t n_blocks) {
    // col tiles before this column + row offset within column
    return tile_col * (2 * n_blocks - tile_col + 1) / 2 + (tile_row - tile_col);
}

// Element index within tile (row-major)
static inline dim_t elem_in_tile(dim_t local_row, dim_t local_col) {
    return local_row * BLOCK_DIM + local_col;
}
```

## Phase 2: Blocked State Infrastructure (Test-Driven)

### 2.1 Test File: `test/src/test_mhipster_block.c` (State Tests Section)

```c
// Test: blocked state allocation succeeds
void test_state_init_blocked(void);

// Test: blocked state length calculation correct
void test_state_len_blocked(void);

// Test: plus state initialization in blocked format
void test_state_plus_blocked(void);

// Test: memory layout is correct (elements accessible via tile indexing)
void test_blocked_memory_layout(void);
```

### 2.2 Implement State Functions: `src/state_blocked.c`

1. `state_blocked_len()` - Calculate storage size for blocked format
2. `state_blocked_init()` - Allocate blocked storage (with padding to BLOCK_DIM)
3. `state_blocked_plus()` - Initialize uniform superposition in blocked format
4. Helper: `pack_idx_blocked()` - Index calculation for blocked layout

**Storage formula**:
```c
// Number of tiles per dimension (rounded up)
n_blocks = (dim + BLOCK_DIM - 1) / BLOCK_DIM;

// Total tiles in lower triangle
n_tiles = n_blocks * (n_blocks + 1) / 2;

// Total elements (each tile is BLOCK_DIM x BLOCK_DIM)
storage_len = n_tiles * BLOCK_DIM * BLOCK_DIM;
```

## Phase 3: Gate Implementation (Incremental)

### 3.1 Z Gate First (Simplest - Diagonal)

**Why Z first**: Only affects off-diagonal blocks, no mixing, tests basic traversal.

**Test**: `test/src/test_mhipster_block.c`
```c
void test_z_blocked(void) {
    // For each qubit count and target:
    // 1. Initialize state directly in blocked format (state_blocked_plus)
    // 2. Apply z_blocked()
    // 3. Extract result values via tile indexing
    // 4. Compare with reference: multiQubitGate() on full matrix
}
```

**Test state generation**: Initialize directly in blocked format, validate against full-matrix reference implementation (not packed format).

**Implementation pattern**:
```c
void z_blocked(state_blocked_t *state, qubit_t target) {
    const dim_t dim = (dim_t)1 << state->qubits;
    const dim_t n_blocks = (dim + BLOCK_DIM - 1) / BLOCK_DIM;
    const dim_t n_tiles = n_blocks * (n_blocks + 1) / 2;
    const dim_t incr = (dim_t)1 << target;

    // Linearized loop over lower-triangular tiles
    // NOTE: Cannot use collapse(2) because inner loop bound depends on outer variable
    #pragma omp parallel for schedule(dynamic)
    for (dim_t t = 0; t < n_tiles; ++t) {
        // Derive tile coordinates from linear index (column-major lower triangular)
        dim_t tb_col = (dim_t)(sqrt(0.25 + 2.0 * t) - 0.5);
        dim_t col_start = tb_col * (tb_col + 1) / 2;
        if (col_start + tb_col < t) ++tb_col;  // Adjust for rounding
        col_start = tb_col * (tb_col + 1) / 2;
        dim_t tb_row = t - col_start + tb_col;

        // Process tile (tb_row, tb_col)
        cplx_t *tile = get_tile(state, tb_row, tb_col);
        z_tile(tile, tb_row, tb_col, incr, dim);
    }
}
```

**Alternative** (simpler but less parallel granularity):
```c
    // Parallelize outer loop only - safe with triangular iteration
    #pragma omp parallel for schedule(dynamic)
    for (dim_t tb_row = 0; tb_row < n_blocks; ++tb_row) {
        for (dim_t tb_col = 0; tb_col <= tb_row; ++tb_col) {
            cplx_t *tile = get_tile(state, tb_row, tb_col);
            z_tile(tile, tb_row, tb_col, incr, dim);
        }
    }
```

### 3.2 X Gate (Permutation - Tests SIMD Swap)

**Implementation**: Block-aware swap with SIMD optimization.

```c
// SIMD swap within tile using NEON/AVX2
static void swap_tile_elements_simd(cplx_t *a, cplx_t *b, dim_t n);
```

### 3.3 Y Gate (Permutation + Negation)

Similar to X with additional negation.

### 3.4 H Gate (Mixing - Main Benefit of Blocking)

**Key optimization**: Use BLAS-3 (zgemm) for tile-level Hadamard transform.

```c
void h_blocked(state_blocked_t *state, qubit_t target) {
    // For each 2x2 group of tiles affected by target qubit:
    // Load 4 tiles, apply blocked Hadamard transform, store back

    // Within tiles: exploit SIMD for parallel element processing
}
```

### 3.5 S, Sdg, T, Tdg Gates (Diagonal Phase)

Similar structure to Z, different phase factors.

## Phase 4: SIMD Optimization

### 4.1 Platform Detection

```c
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    #include <arm_neon.h>
    #define USE_NEON 1
#elif defined(__AVX2__)
    #include <immintrin.h>
    #define USE_AVX2 1
#endif
```

### 4.2 SIMD Primitives

**File**: `src/mhipster_block.c`

```c
// Swap n complex elements (NEON version)
// Uses vld1q_f64 for direct load (no deinterleaving overhead)
static inline void simd_swap_neon(cplx_t *a, cplx_t *b, dim_t n) {
    dim_t i = 0;
    // Process 2 complex numbers at a time (4 doubles = 32 bytes)
    for (; i + 2 <= n; i += 2) {
        float64x2_t va0 = vld1q_f64((double*)(a + i));      // a[i]
        float64x2_t va1 = vld1q_f64((double*)(a + i + 1));  // a[i+1]
        float64x2_t vb0 = vld1q_f64((double*)(b + i));      // b[i]
        float64x2_t vb1 = vld1q_f64((double*)(b + i + 1));  // b[i+1]
        vst1q_f64((double*)(a + i), vb0);
        vst1q_f64((double*)(a + i + 1), vb1);
        vst1q_f64((double*)(b + i), va0);
        vst1q_f64((double*)(b + i + 1), va1);
    }
    // Scalar remainder
    for (; i < n; ++i) {
        cplx_t tmp = a[i];
        a[i] = b[i];
        b[i] = tmp;
    }
}

// Conjugate n complex elements (NEON version)
static inline void simd_conj_neon(cplx_t *a, dim_t n) {
    const float64x2_t neg_mask = {1.0, -1.0};  // Preserves real, negates imag
    for (dim_t i = 0; i < n; ++i) {
        float64x2_t v = vld1q_f64((double*)(a + i));
        v = vmulq_f64(v, neg_mask);
        vst1q_f64((double*)(a + i), v);
    }
}

// Scale n complex elements by real scalar (NEON version)
static inline void simd_scale_real_neon(cplx_t *a, double alpha, dim_t n) {
    float64x2_t valpha = vdupq_n_f64(alpha);
    for (dim_t i = 0; i < n; ++i) {
        float64x2_t v = vld1q_f64((double*)(a + i));
        v = vmulq_f64(v, valpha);
        vst1q_f64((double*)(a + i), v);
    }
}
```

### 4.3 Fallback to BLAS/Scalar

When SIMD unavailable or for small n, fall back to:
- `cblas_zswap`, `cblas_zscal` for BLAS-1
- Scalar loops for remainder elements

## Phase 5: OpenMP Parallelization

### 5.1 Tile-Level Parallelism

**Important**: Cannot use `collapse(2)` with non-rectangular (triangular) iteration spaces.

```c
// Option A: Linearized tile index (maximum parallelism)
const dim_t n_tiles = n_blocks * (n_blocks + 1) / 2;
#pragma omp parallel for schedule(dynamic)
for (dim_t t = 0; t < n_tiles; ++t) {
    dim_t tb_row, tb_col;
    tile_coords_from_index(t, &tb_row, &tb_col);  // Helper function
    process_tile(tb_row, tb_col);
}

// Option B: Parallelize outer loop only (simpler, sufficient for most cases)
#pragma omp parallel for schedule(dynamic)
for (dim_t tb_row = 0; tb_row < n_blocks; ++tb_row) {
    for (dim_t tb_col = 0; tb_col <= tb_row; ++tb_col) {
        process_tile(tb_row, tb_col);
    }
}
```

### 5.2 Nested Parallelism for Large Tiles

For very large systems (n > 14), enable nested parallelism within tiles:

```c
#pragma omp parallel for if(BLOCK_DIM >= 128)
for (dim_t row = 0; row < BLOCK_DIM; ++row) {
    process_tile_row(tile, row);
}
```

## Phase 6: Benchmark Integration

### 6.1 Modify `benchmark/src/bench_mixed.c`

Add blocked variant:

```c
#include "state_blocked.h"
#include "mhipster_block.h"

bench_result_t bench_qlib_blocked(qubit_t qubits, const char *gate_name,
                                   void (*gate_fn)(state_blocked_t*, qubit_t),
                                   int iterations, int warmup) {
    state_blocked_t state = {0};
    state_blocked_plus(&state, qubits);

    // Warm-up
    for (int i = 0; i < warmup; ++i) {
        gate_fn(&state, 0);
    }

    // Timed iterations
    double start = get_time_ns();
    for (int i = 0; i < iterations; ++i) {
        gate_fn(&state, i % qubits);
    }
    double elapsed = get_time_ns() - start;

    state_blocked_free(&state);  // Always clean up

    return (bench_result_t){
        .avg_ns = elapsed / iterations,
        .iterations = iterations
    };
}
```

### 6.2 Gate Definition Extension

```c
typedef struct gate_def {
    const char *name;
    qs_error_t (*fn_packed)(state_t*, qubit_t);
    void (*fn_blocked)(state_blocked_t*, qubit_t);  // NEW - uses blocked state type
    const cplx_t *mat;
} gate_def_t;

static const gate_def_t gates[] = {
    {"X", x, x_blocked, XMAT},
    {"H", h, h_blocked, HMAT},
    {"Z", z, z_blocked, ZMAT},
    {NULL, NULL, NULL, NULL}
};
```

### 6.3 New Benchmark Output

```
qubits: 12, dim: 4096

Gate: H
  qlib_packed     10.5 ms (95,238 ops/sec)
  qlib_blocked     7.2 ms (138,889 ops/sec)   <<< NEW
  blas_dense      850.0 ms (1,176 ops/sec)
  Speedup: 1.46x blocked vs packed, 118x vs dense
```

## Phase 7: CMake Integration

### 7.1 Modify `src/CMakeLists.txt` (build config only)

```cmake
set(Q_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/state.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/state_blocked.c"    # NEW
    "${CMAKE_CURRENT_SOURCE_DIR}/gate.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/qhipster.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/mhipster.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/mhipster_block.c"   # NEW
)
```

### 7.2 Modify `test/CMakeLists.txt` (build config only)

```cmake
# Shared test helper library (extracted from test_mhipster.c)
add_library(test_helpers STATIC
    src/test_helpers.c
    src/gatemat.c
)
target_include_directories(test_helpers PUBLIC include)
target_link_libraries(test_helpers PRIVATE q)

# Existing test (update to use shared helpers)
add_executable(test_mhipster src/test_mhipster.c)
target_link_libraries(test_mhipster PRIVATE q unity test_helpers)

# Add blocked gate tests
add_executable(test_mhipster_block src/test_mhipster_block.c)
target_link_libraries(test_mhipster_block PRIVATE q unity test_helpers)
add_test(NAME test_mhipster_block COMMAND test_mhipster_block)
```

**Note**: CMakeLists.txt modifications are build configuration, not source code changes. The constraint "do not edit existing files except bench_mixed.c" applies to source/header files, not build scripts.

**Helper extraction**: Move `multiQubitGate()` and related reference implementations from `test_mhipster.c` to `test_helpers.c`. This avoids duplicate `main()` symbols when linking multiple test executables.

## Implementation Order (Tight TDD - Test Immediately Followed by Implementation)

| Phase | Task | Files |
|-------|------|-------|
| 1 | Create headers | `include/state_blocked.h`, `include/mhipster_block.h` |
| 2 | State tests + impl | `test/src/test_mhipster_block.c` (state section), `src/state_blocked.c` |
| 3 | Z gate: test -> impl | `test/src/test_mhipster_block.c` (test_z), `src/mhipster_block.c` (z_blocked) |
| 4 | X gate: test -> impl | Same files (test_x -> x_blocked with SIMD) |
| 5 | Y gate: test -> impl | Same files (test_y -> y_blocked) |
| 6 | S gate: test -> impl | Same files (test_s -> s_blocked) |
| 7 | Sdg gate: test -> impl | Same files (test_sdg -> sdg_blocked) |
| 8 | T gate: test -> impl | Same files (test_t -> t_blocked) |
| 9 | Tdg gate: test -> impl | Same files (test_tdg -> tdg_blocked) |
| 10 | H gate: test -> impl | Same files (test_h -> h_blocked, key optimization target) |
| 11 | SIMD tuning | `src/mhipster_block.c` (NEON/AVX2 paths) |
| 12 | OpenMP tuning | `src/mhipster_block.c` (parallel regions) |
| 13 | Benchmark | `benchmark/src/bench_mixed.c` (add blocked variant) |

**TDD Rule**: Never implement a gate without its test passing first. Each gate has:
1. Test written and failing (references `multiQubitGate()` from test_mhipster.c)
2. Implementation written
3. Test passing
4. Move to next gate

## Success Criteria

1. **Correctness**: All blocked gates produce identical results to packed gates (tolerance: 1e-12)
2. **Performance**: H gate blocked >= 1.3x faster than packed for n >= 10 qubits
3. **Benchmarks**: `./benchmark/bench_mixed` runs with `--blocked` option showing comparison
4. **Tests pass**: `ctest` includes blocked gate tests
5. **Documentation**: Clear guidance that blocked format is recommended for n >= 10 qubits (memory overhead < 6%)

## Risk Mitigation

1. **Memory overhead**: Blocked layout has variable overhead depending on qubit count:

   | Qubits | Dim   | Packed Elements | Blocked Elements | Overhead |
   |--------|-------|-----------------|------------------|----------|
   | 5      | 32    | 528             | 4,096            | **7.8x** |
   | 8      | 256   | 32,896          | 40,960           | 24%      |
   | 10     | 1024  | 524,800         | 557,056          | 6%       |
   | 12     | 4096  | 8,386,560       | 8,519,680        | 1.6%     |
   | 14     | 16384 | 134,225,920     | 135,274,496      | 0.8%     |

   Overhead is severe for small systems (n < 8) due to tile padding, but negligible (<2%) for n ≥ 12 where blocked format provides performance benefits. Mitigation: Document that blocked format is recommended only for n ≥ 10 qubits; use packed format for smaller systems.

2. **SIMD portability**: ARM NEON vs x86 AVX2 differences. Mitigation: Abstract SIMD behind macros with scalar/BLAS fallback when intrinsics unavailable.

3. **Test complexity**: Validating blocked against packed requires careful test state generation. Mitigation: Use reference implementation (multiQubitGate in test harness) rather than comparing blocked vs packed directly.

## Files to Create/Modify

**New source/header files** (completely parallel infrastructure):
- `include/state_blocked.h` - Blocked state type and functions
- `include/mhipster_block.h` - Blocked gate declarations
- `src/state_blocked.c` - Blocked state implementation
- `src/mhipster_block.c` - Blocked gate implementations
- `test/src/test_mhipster_block.c` - Blocked format tests
- `test/src/test_helpers.c` - Shared test helpers (extracted from test_mhipster.c)
- `test/include/test_helpers.h` - Header for multiQubitGate() and related functions

**Modified files** (source code):
- `benchmark/src/bench_mixed.c` - Add `bench_qlib_blocked()` function (ONLY source file modified)
- `test/src/test_mhipster.c` - Remove helper functions (moved to test_helpers.c), add include

**Modified files** (build configuration only):
- `src/CMakeLists.txt` - Add state_blocked.c and mhipster_block.c to Q_SOURCES
- `test/CMakeLists.txt` - Add test_helpers library and test_mhipster_block target
- `benchmark/CMakeLists.txt` - Ensure blocked sources are linked
