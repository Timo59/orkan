# Implementation Plan: Optimizing `h_tiled()` -- Optimizations 7 and 8

**Date**: 2026-02-14
**Target file**: `/Users/timo/Code/orkan/src/gate/tiled/gate_tiled_h.c`
**Profiling harness**: `/Users/timo/Code/orkan/profile/src/profile_h_tiled.c`

**Build configurations**:
- Debug: `cmake-build-debug` with `LOG_TILE_DIM=3` (8x8 tiles, `TILE_DIM=8`, `TILE_SIZE=64`)
- Release: `cmake-build-release` with `LOG_TILE_DIM=5` (32x32 tiles, `TILE_DIM=32`, `TILE_SIZE=1024`)

**Critical note**: Correctness tests run in Debug mode (8x8 tiles). Profiling runs in Release mode (32x32 tiles). The code must work correctly for ALL values of `LOG_TILE_DIM`.

---

## Background

Two performance optimizations were identified during code review of the tiled gate implementations:

- **Optimization 7**: Restructure the within-tile `bc` loop to expose contiguous memory runs for SIMD vectorization
- **Optimization 8**: Use a stack-allocated scratch buffer in cross-tile sub-case (b) to eliminate stride-512 column access

Both optimizations apply to all single-qubit tiled gates. This plan uses `h_tiled()` as the exemplar.

---

## Step 1: Create the Profiling Harness

### 1.1 Create `profile/src/profile_h_tiled.c`

This file follows the exact pattern of `profile/src/profile_x_tiled.c`, substituting `h_tiled`/`h_packed` for `x_tiled`/`x_packed`.

```c
/**
 * @file profile_h_tiled.c
 * @brief Micro-profiling harness for h_tiled()
 *
 * Profiles the tiled Hadamard gate across multiple qubit counts and target
 * qubits, separating within-tile (target < LOG_TILE_DIM) from cross-tile
 * (target >= LOG_TILE_DIM) code paths.
 *
 * The Hadamard gate is a non-permutation, non-diagonal gate that mixes all
 * 4 elements in each butterfly block. Unlike X (swap-only), H performs
 * arithmetic (4 additions + 1 multiply per output element), making it
 * compute-bound for small problem sizes and memory-bound for large ones.
 *
 * Usage:
 *   ./profile_h_tiled [options]
 *
 * Options:
 *   --min-qubits N    Minimum qubit count (default: LOG_TILE_DIM)
 *   --max-qubits N    Maximum qubit count (default: 12)
 *   --warmup N        Warmup iterations (default: 200)
 *   --iters N         Measurement iterations (default: 2000)
 *   --csv FILE        Write CSV output to FILE
 *   --target N        Profile only target qubit N (default: all valid targets)
 *   --compare         Also profile h_packed() for comparison
 *   -v, --verbose     Print extra info (state sizes, tile config)
 */

#include "profile.h"
#include <string.h>
#include <getopt.h>

/* Forward declarations for backend functions (not in public gate.h) */
extern void h_tiled(state_t *state, const qubit_t target);
extern void h_packed(state_t *state, const qubit_t target);

/*
 * =====================================================================================================================
 * Command-line options
 * =====================================================================================================================
 */

typedef struct {
    qubit_t min_qubits;
    qubit_t max_qubits;
    int     warmup;
    int     iters;
    char    csv_path[256];
    int     specific_target;  /* -1 = all targets */
    int     compare;
    int     verbose;
} options_t;

static options_t parse_args(int argc, char *argv[]) {
    options_t o = {
        .min_qubits = (LOG_TILE_DIM < 2) ? 2 : LOG_TILE_DIM,
        .max_qubits = PROF_DEFAULT_MAX_QUBITS,
        .warmup     = PROF_DEFAULT_WARMUP,
        .iters      = PROF_DEFAULT_ITERS,
        .csv_path   = "",
        .specific_target = -1,
        .compare    = 0,
        .verbose    = 0,
    };

    static struct option long_opts[] = {
        {"min-qubits", required_argument, NULL, 'm'},
        {"max-qubits", required_argument, NULL, 'M'},
        {"warmup",     required_argument, NULL, 'w'},
        {"iters",      required_argument, NULL, 'i'},
        {"csv",        required_argument, NULL, 'c'},
        {"target",     required_argument, NULL, 't'},
        {"compare",    no_argument,       NULL, 'C'},
        {"verbose",    no_argument,       NULL, 'v'},
        {"help",       no_argument,       NULL, 'h'},
        {NULL, 0, NULL, 0}
    };

    int ch;
    while ((ch = getopt_long(argc, argv, "m:M:w:i:c:t:Cvh", long_opts, NULL)) != -1) {
        switch (ch) {
            case 'm': o.min_qubits = (qubit_t)atoi(optarg); break;
            case 'M': o.max_qubits = (qubit_t)atoi(optarg); break;
            case 'w': o.warmup = atoi(optarg); break;
            case 'i': o.iters  = atoi(optarg); break;
            case 'c': strncpy(o.csv_path, optarg, sizeof(o.csv_path) - 1); break;
            case 't': o.specific_target = atoi(optarg); break;
            case 'C': o.compare = 1; break;
            case 'v': o.verbose = 1; break;
            case 'h':
                printf("Usage: %s [--min-qubits N] [--max-qubits N] [--warmup N] "
                       "[--iters N] [--csv FILE] [--target N] [--compare] [-v]\n", argv[0]);
                exit(0);
            default: break;
        }
    }

    if (o.min_qubits < 2) o.min_qubits = 2;
    return o;
}

/*
 * =====================================================================================================================
 * Main
 * =====================================================================================================================
 */

int main(int argc, char *argv[]) {
    options_t opts = parse_args(argc, argv);

    printf("=== h_tiled() profiling ===\n");
    printf("Tile config: TILE_DIM=%d, LOG_TILE_DIM=%d, TILE_SIZE=%d\n",
           TILE_DIM, LOG_TILE_DIM, TILE_SIZE);
    printf("Warmup: %d, Iterations: %d\n", opts.warmup, opts.iters);
    printf("Qubit range: %d to %d\n\n", opts.min_qubits, opts.max_qubits);

    FILE *csv = NULL;
    if (opts.csv_path[0]) {
        csv = fopen(opts.csv_path, "w");
        if (!csv) {
            fprintf(stderr, "Cannot open %s for writing\n", opts.csv_path);
            return 1;
        }
        prof_print_csv_header(csv);
    }

    /* Print table header */
    prof_print_header();

    for (qubit_t n = opts.min_qubits; n <= opts.max_qubits; n++) {
        qubit_t max_target = n - 1;

        if (opts.verbose) {
            uint64_t dim = (uint64_t)1 << n;
            uint64_t nt = (dim + TILE_DIM - 1) / TILE_DIM;
            printf("\n--- %d qubits: dim=%llu, n_tiles=%llu, "
                   "data=%.1f KiB ---\n",
                   n, (unsigned long long)dim, (unsigned long long)nt,
                   (double)prof_tiled_bytes(n) / 2048.0);
        }

        for (qubit_t t = 0; t <= max_target; t++) {
            if (opts.specific_target >= 0 && t != (qubit_t)opts.specific_target)
                continue;

            /* Profile h_tiled */
            prof_result_t r = prof_run_1q(
                h_tiled, "h_tiled", MIXED_TILED,
                n, t, opts.warmup, opts.iters);
            prof_print_result(&r);
            if (csv) prof_print_csv_row(csv, &r);

            /* Optionally compare with h_packed */
            if (opts.compare) {
                prof_result_t rp = prof_run_1q(
                    h_packed, "h_packed", MIXED_PACKED,
                    n, t, opts.warmup, opts.iters);
                prof_print_result(&rp);
                if (csv) prof_print_csv_row(csv, &rp);
            }
        }
    }

    if (csv) {
        fclose(csv);
        printf("\nCSV written to: %s\n", opts.csv_path);
    }

    /*
     * =====================================================================================================================
     * Scaling analysis: fix target, vary qubits
     * =====================================================================================================================
     */
    printf("\n\n=== Scaling analysis (fixed target, varying qubits) ===\n");
    prof_print_header();

    /* Within-tile target (target=0, always within-tile) */
    for (qubit_t n = opts.min_qubits; n <= opts.max_qubits; n++) {
        prof_result_t r = prof_run_1q(
            h_tiled, "h_tiled", MIXED_TILED,
            n, 0, opts.warmup, opts.iters);
        prof_print_result(&r);
        if (csv) prof_print_csv_row(csv, &r);
    }

    /* Cross-tile target (target=LOG_TILE_DIM, when n > LOG_TILE_DIM) */
    if (opts.max_qubits > LOG_TILE_DIM) {
        printf("\n");
        for (qubit_t n = LOG_TILE_DIM + 1; n <= opts.max_qubits; n++) {
            prof_result_t r = prof_run_1q(
                h_tiled, "h_tiled", MIXED_TILED,
                n, LOG_TILE_DIM, opts.warmup, opts.iters);
            prof_print_result(&r);
            if (csv) prof_print_csv_row(csv, &r);
        }
    }

    printf("\n=== Profiling complete ===\n");
    return 0;
}
```

### 1.2 Add Build Target to `profile/CMakeLists.txt`

Append after the existing `profile_x_tiled` target (after line 22):

```cmake
# --- profile_h_tiled ---
add_executable(profile_h_tiled "${PROF_SRC_DIR}/profile_h_tiled.c")

set_target_properties(profile_h_tiled PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/profile"
)

target_include_directories(profile_h_tiled PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_link_libraries(profile_h_tiled PRIVATE q)
```

### 1.3 Build and Verify

```bash
cd /Users/timo/Code/orkan/cmake-build-release
cmake .. && cmake --build .
./profile/profile_h_tiled --min-qubits 5 --max-qubits 5 --iters 10
```

**Success criteria**: Binary runs, prints a table header followed by 5 rows (targets 0-4), each showing `h_tiled`, qubits `5`, path `within-tile`.

---

## Step 2: Baseline Measurement

```bash
cd /Users/timo/Code/orkan/cmake-build-release
./profile/profile_h_tiled \
    --min-qubits 5 --max-qubits 12 --iters 5000 \
    --csv /Users/timo/Code/orkan/profile/results_h_tiled_baseline.csv \
    --compare -v
```

**Key data points to record**:
- Within-tile targets 1-4 (affected by Optimization 7)
- Cross-tile targets at LOG_TILE_DIM for n >= 7 (affected by Optimization 8)
- `h_packed` numbers as a fixed reference baseline

---

## Step 3: Implement Optimization 7 -- Vectorized Within-Tile Inner Loop

### 3.1 Index Math

`insertBit0(bc, target)` produces `c0` values in contiguous runs of length `2^target`:

- **target=0**: runs of 1 (c0 = 0, 2, 4, ...)
- **target=1**: runs of 2 (c0 = 0, 1, 4, 5, 8, 9, ...)
- **target=2**: runs of 4 (c0 = 0, 1, 2, 3, 8, 9, 10, 11, ...)
- **target=k**: runs of `2^k`, base = `g << (k+1)`, c1 = c0 + `2^k`

The key identity:
```
insertBit0(bc, target) = (g << (target + 1)) + off
```
where `g = bc >> target` and `off = bc & ((1 << target) - 1)`.

**Concrete example** (TILE_DIM=32, target=2):
- `incr = 4`, `n_base = 16`, `group_len = 4`, `n_groups = 4`
- Group 0: c0 = {0,1,2,3}, c1 = {4,5,6,7}
- Group 1: c0 = {8,9,10,11}, c1 = {12,13,14,15}
- Group 2: c0 = {16,17,18,19}, c1 = {20,21,22,23}
- Group 3: c0 = {24,25,26,27}, c1 = {28,29,30,31}

All accesses within a group are contiguous, enabling vectorization.

### 3.2 Code Change

**File**: `src/gate/tiled/gate_tiled_h.c`

**Replace lines 48-68** (the `br`/`bc` nested loop inside `if (target < LOG_TILE_DIM)`) with:

```c
                /* Number of contiguous-column groups and elements per group.
                 * insertBit0(bc, target) produces c0 values in contiguous runs
                 * of length group_len = 2^target, starting at base = g << (target+1).
                 * Within a run: c0 = base + off, c1 = base + incr + off.
                 * Decomposing bc into (group, offset) exposes contiguous inner
                 * loops that the compiler can vectorize. */
                const gate_idx_t group_len = incr;              /* 1 << target */
                const gate_idx_t n_groups  = n_base >> target;  /* dim_tile >> (target+1) */

                for (gate_idx_t br = 0; br < n_base; ++br) {
                    const gate_idx_t r0 = insertBit0(br, target);
                    const gate_idx_t r1 = r0 | incr;

                    const gate_idx_t offset_r0 = r0 * TILE_DIM + offset_tile;
                    const gate_idx_t offset_r1 = r1 * TILE_DIM + offset_tile;

                    cplx_t * restrict p_r0 = data + offset_r0;
                    cplx_t * restrict p_r1 = data + offset_r1;

                    for (gate_idx_t g = 0; g < n_groups; ++g) {
                        const gate_idx_t c0_base = g << (target + 1);
                        const gate_idx_t c1_base = c0_base + incr;

                        #pragma clang loop vectorize(enable)
                        for (gate_idx_t off = 0; off < group_len; ++off) {
                            cplx_t v00 = p_r0[c0_base + off];
                            cplx_t v01 = p_r0[c1_base + off];
                            cplx_t v10 = p_r1[c0_base + off];
                            cplx_t v11 = p_r1[c1_base + off];

                            p_r0[c0_base + off] = 0.5 * (v00 + v01 + v10 + v11);
                            p_r0[c1_base + off] = 0.5 * (v00 - v01 + v10 - v11);
                            p_r1[c0_base + off] = 0.5 * (v00 + v01 - v10 - v11);
                            p_r1[c1_base + off] = 0.5 * (v00 - v01 - v10 + v11);
                        }
                    }
                }
```

### 3.3 Why `restrict` Helps

The original code uses `data[c0 + offset_r0]` where the compiler cannot prove the four butterfly addresses don't alias. By introducing `p_r0 = data + offset_r0` and `p_r1 = data + offset_r1` with `restrict`, we tell the compiler these row regions don't overlap. This is correct because `r0 != r1` (they differ in bit `target`).

### 3.4 Edge Cases

- **target=0**: `group_len=1`, `n_groups=n_base`. Inner loop has trip count 1 -- degenerates to original scalar behavior. No regression.
- **target=1**: `group_len=2`. Minimal vectorization benefit (32 bytes per iteration).
- **target=4** (TILE_DIM=32 only): `group_len=16`, `n_groups=1`. Best case for vectorization.
- **Debug mode** (TILE_DIM=8): target range 0-2 for within-tile. target=2 gives `group_len=4`, `n_groups=1`. Works correctly.
- **dim < TILE_DIM**: `dim_tile = dim`, `n_base = dim/2`. Values scale proportionally.

---

## Step 4: Implement Optimization 8 -- Scratch Buffer for Sub-case (b)

### 4.1 Problem

In sub-case (b), tile T(tr0,tc1) is accessed via its transpose T(tc1,tr0):

```c
data[lr + lc * TILE_DIM + offset_t01]   // stride TILE_DIM in lc
```

For TILE_DIM=32, each `lc` increment jumps `32 * 16 = 512 bytes`. Every access touches a different cache line, causing severe L1 pressure.

### 4.2 Strategy

1. Load-transpose the t01 tile into a stack scratch buffer (applying conjugation)
2. Process the butterfly with all four tile pointers having contiguous row-major access
3. Write-back-transpose the modified scratch into storage (applying conjugation)

**Scratch size**: `TILE_DIM * TILE_DIM * sizeof(cplx_t)` = 16 KiB (Release) or 1 KiB (Debug). Fits easily on the stack (OpenMP default stack is 2-8 MiB).

### 4.3 Code Change

**File**: `src/gate/tiled/gate_tiled_h.c`

**Replace the entire `else if (btr != btc)` block** (lines 124-155) with:

```c
                else if (btr != btc) {
                    /*
                     * Sub-case b: off-diagonal (btr != btc), tr0 <= tc1.
                     * T(tr0,tc1) not in lower triangle; access via T(tc1,tr0).
                     * stored[lc,lr] = conj(rho01[lr,lc]).
                     *
                     * Optimization: transpose-load the t01 tile into a stack
                     * scratch buffer so the inner loop accesses all 4 tiles
                     * with contiguous (row-major) stride. After processing,
                     * transpose-write the scratch back to storage.
                     *
                     * Scratch size: TILE_DIM * TILE_DIM * sizeof(cplx_t)
                     *   = 32*32*16 = 16 KiB (Release) or 8*8*16 = 1 KiB (Debug).
                     */
                    const gate_idx_t offset_t01 = tile_off(tc1, tr0);
                    cplx_t scratch01[TILE_SIZE];

                    /* Load-transpose: scratch01[lr][lc] = rho01[lr,lc]
                     *                                   = conj(stored[lc,lr])
                     *                                   = conj(data[lc*TILE_DIM + lr + offset_t01]) */
                    for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                        for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                            scratch01[lr * TILE_DIM + lc] =
                                conj(data[lc * TILE_DIM + lr + offset_t01]);
                        }
                    }

                    /* Process: all 4 tiles now have contiguous row-major access */
                    for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                        cplx_t * restrict row00 = data + lr * TILE_DIM + offset_t00;
                        cplx_t * restrict row01 = scratch01 + lr * TILE_DIM;
                        cplx_t * restrict row10 = data + lr * TILE_DIM + offset_t10;
                        cplx_t * restrict row11 = data + lr * TILE_DIM + offset_t11;

                        #pragma clang loop vectorize(enable)
                        for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                            cplx_t v00 = row00[lc];
                            cplx_t v01 = row01[lc];
                            cplx_t v10 = row10[lc];
                            cplx_t v11 = row11[lc];

                            row00[lc] = 0.5 * (v00 + v01 + v10 + v11);
                            row01[lc] = 0.5 * (v00 - v01 + v10 - v11);
                            row10[lc] = 0.5 * (v00 + v01 - v10 - v11);
                            row11[lc] = 0.5 * (v00 - v01 - v10 + v11);
                        }
                    }

                    /* Write-back-transpose: stored[lc,lr] = conj(scratch01[lr][lc])
                     * i.e. data[lc*TILE_DIM + lr + offset_t01] = conj(scratch01[lr*TILE_DIM + lc]) */
                    for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                        for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                            data[lc * TILE_DIM + lr + offset_t01] =
                                conj(scratch01[lr * TILE_DIM + lc]);
                        }
                    }
                }
```

### 4.4 Correctness Argument

- **Load-transpose**: `stored[lc,lr] = conj(rho01[lr,lc])` by the Hermitian storage convention, so `scratch01[lr*TILE_DIM + lc] = conj(stored[lc,lr]) = rho01[lr,lc]`. Correct.
- **Processing**: Identical butterfly arithmetic to the original. `row01[lc]` writes go to scratch (not storage). The processing loop is structurally identical to sub-case (a).
- **Write-back-transpose**: `stored[lc,lr] = conj(new_rho01[lr,lc]) = conj(scratch01[lr*TILE_DIM + lc])`. Matches the original code's `data[lr + lc * TILE_DIM + offset_t01] = conj(new01)`.
- **Read-before-write safety**: All reads of t01 happen in the load step (before any writes). All writes happen in the write-back step (after all processing). No read-after-write hazard.
- **OpenMP safety**: `scratch01` is a local variable inside the OMP parallel loop body. Each thread gets its own stack copy.

### 4.5 Performance Characteristics

**Net effect**: Two bad-stride passes (load + write-back) of 16 KiB each, versus the original 32 bad-stride passes (one per `lr`) each touching 32 elements at stride 512. The main processing loop runs at full SIMD throughput with contiguous access.

---

## Step 5: Verify Correctness

Build and test in Debug mode (8x8 tiles):

```bash
cd /Users/timo/Code/orkan/cmake-build-debug
cmake --build . && ctest --output-on-failure
```

**Success criteria**: All tests pass, 0 failures. The gate tests cover:
- Every qubit count from 2 to ~7
- Every valid target qubit for each qubit count
- Random density matrix inputs verified against a Kronecker-product reference
- Both within-tile (target 0 to LOG_TILE_DIM-1) and cross-tile (LOG_TILE_DIM+) targets

### If Tests Fail

**Optimization 7 debug**: Print `c0_base + off` values for a small case (n=3, target=1) and compare against `insertBit0(bc, target)` for all `bc` in `[0, n_base)`. Every value must appear exactly once.

**Optimization 8 debug**: After load-transpose, assert `scratch01[lr * TILE_DIM + lc] == conj(data[lc * TILE_DIM + lr + offset_t01])`. Common mistake: swapping `lr`/`lc` in the transpose indexing, or forgetting `conj()`.

---

## Step 6: Post-Optimization Measurement

```bash
cd /Users/timo/Code/orkan/cmake-build-release
cmake --build .
./profile/profile_h_tiled \
    --min-qubits 5 --max-qubits 12 --iters 5000 \
    --csv /Users/timo/Code/orkan/profile/results_h_tiled_optimized.csv \
    --compare -v
```

### What to Look For

**Optimization 7** (within-tile, targets 1-4):
- target=0: minimal change (group_len=1, no vectorization gain)
- target=2 (group_len=4): measurable speedup (64 bytes per vector iteration)
- target=4 (group_len=16): best speedup
- Most visible at n=9-11 where data fits in cache (compute-bound regime)
- Expected: 10-30% improvement on within-tile paths for targets >= 2

**Optimization 8** (cross-tile sub-case b):
- Cross-tile targets at n=10-11 where data fits in L2/SLC
- Expected: 15-40% improvement on sub-case (b) tile-blocks
- Overall cross-tile improvement smaller (sub-cases a and c unaffected)

### Optional: Isolated A/B Comparison

1. Apply only Optimization 7, profile as `results_h_tiled_opt7.csv`
2. Apply Optimization 8 on top, profile as `results_h_tiled_opt78.csv`

---

## Step 7: Write Profiling Report

Create `profile/REPORT_h_tiled.md` following the structure of `profile/REPORT_x_tiled.md`:

1. Executive Summary (key findings, overall speedup table)
2. Methodology (harness, iterations, tools)
3. Baseline Results (scaling table, within-tile vs cross-tile, h_tiled vs h_packed)
4. Optimization 7 Results (before/after table for within-tile targets)
5. Optimization 8 Results (before/after table for cross-tile paths)
6. Combined Results (final comparison, h_tiled vs h_packed after optimization)
7. Remaining Optimization Opportunities (sub-case c, blocked transpose, NEON conj)

---

## Files Summary

| Action | File | Step |
|--------|------|------|
| CREATE | `profile/src/profile_h_tiled.c` | 1.1 |
| MODIFY | `profile/CMakeLists.txt` (append after line 22) | 1.2 |
| MODIFY | `src/gate/tiled/gate_tiled_h.c` (replace lines 48-68) | 3 |
| MODIFY | `src/gate/tiled/gate_tiled_h.c` (replace lines 124-155) | 4 |
| CREATE | `profile/REPORT_h_tiled.md` | 7 |

## Commands Summary

```bash
# Step 1: Build profiling harness (after creating files)
cd /Users/timo/Code/orkan/cmake-build-release && cmake .. && cmake --build .

# Step 2: Baseline
./profile/profile_h_tiled --min-qubits 5 --max-qubits 12 --iters 5000 \
    --csv /Users/timo/Code/orkan/profile/results_h_tiled_baseline.csv --compare -v

# Steps 3+4: Apply optimizations to gate_tiled_h.c

# Step 5: Correctness (Debug build)
cd /Users/timo/Code/orkan/cmake-build-debug && cmake --build . && ctest --output-on-failure

# Step 6: Post-optimization (Release build)
cd /Users/timo/Code/orkan/cmake-build-release && cmake --build .
./profile/profile_h_tiled --min-qubits 5 --max-qubits 12 --iters 5000 \
    --csv /Users/timo/Code/orkan/profile/results_h_tiled_optimized.csv --compare -v
```

## Invariants to Preserve

1. **Hadamard butterfly**: `rho'_00 = 0.5*(+v00+v01+v10+v11)`, `rho'_01 = 0.5*(+v00-v01+v10-v11)`, `rho'_10 = 0.5*(+v00+v01-v10-v11)`, `rho'_11 = 0.5*(+v00-v01-v10+v11)`
2. **Hermitian storage**: Reading T(r,c) with r < c requires `conj(data[lc*TILE_DIM + lr + tile_off(c,r)])`. Writing requires storing `conj(value)`.
3. **Do not change sub-case (a) or (c)**: Only modify the within-tile `bc` loop and sub-case (b).
4. **Do not change OpenMP structure**: The `#pragma omp parallel for` directives must remain unchanged.

## Applicability to Other Gates

Both optimizations apply to all single-qubit tiled gates. After validating on H:

| Group | Gates | Opt 7 | Opt 8 |
|-------|-------|-------|-------|
| Mixing (4 outputs) | H, HY, Rx, Ry | Same loop structure, different butterfly | Same scratch pattern, different butterfly |
| Diagonal (2 outputs) | Rz, S, Sdg, T, Tdg, Z, P | Same loop structure, only a01/a10 modified | Simpler scratch (only t01 read+write) |
| Swap-based | X, Y | Same loop structure, swap/phase-swap | Same scratch pattern |
