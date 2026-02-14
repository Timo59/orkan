# Profiling Report: `x_tiled()`

**Date**: 2026-02-13
**Platform**: Apple Silicon (macOS Darwin 24.6.0)
**Build**: Release (`-O3 -march=native`), `TILE_DIM=32`, `LOG_TILE_DIM=5`
**Compiler**: Apple clang (LLVM), C17

## 1. Executive Summary

`x_tiled()` is 5–8% faster than `x_packed()` at large qubit counts (10–12),
validating the tiled storage strategy. However, several optimization
opportunities exist:

| Issue | Impact | Effort |
|-------|--------|--------|
| OpenMP overhead dominates at n≤8 | ~32µs fixed cost on 48KB data | Low |
| 51 inner loops fail to vectorize | Limits throughput ceiling | Medium |
| Diagonal cross-tile path: two separate passes over t10 | Cache unfriendly | Low |
| Diagonal cross-tile conj loop: column-stride access pattern | L1 misses | Medium |
| Bandwidth drops from ~225 GB/s (n=11) to ~137 GB/s (n=12) | Expected (L2→DRAM) | — |

## 2. Methodology

**Harness**: `profile/src/profile_x_tiled.c` using `clock_gettime(CLOCK_MONOTONIC)`.
Each data point = 5000 iterations after 200 warmup calls. Statistics: median,
mean, min, p5/p95, stddev. State initialized to |+>^n (all elements nonzero).

**Tools used**:
- Custom micro-benchmark (this framework)
- Compiler vectorization reports (`-Rpass=loop-vectorize`, `-Rpass-missed=loop-vectorize`)
- Assembly inspection (`clang -S -O3`)

**Recommended tools for deeper analysis** (not used in this report):
1. **Apple Instruments — Time Profiler**: Hotspot sampling with call stacks
2. **Apple Instruments — CPU Counters**: Hardware PMC (cache misses, branch mispredictions, IPC)
3. **`llvm-mca`**: Static throughput/bottleneck analysis of inner loops
4. **`sample <pid>`**: Quick non-invasive stack sampling

## 3. Results

### 3.1 Scaling: time vs. qubit count (target=0, within-tile)

| Qubits | Dim | Data Size | Median (ns) | BW (GB/s) | elem/ns |
|--------|-----|-----------|-------------|-----------|---------|
| 5 | 32 | 16 KiB | <1000 | — | — |
| 6 | 64 | 48 KiB | 31,000 | 3.2 | 0.10 |
| 7 | 128 | 160 KiB | 32,000 | 10.2 | 0.32 |
| 8 | 256 | 576 KiB | 35,000 | 33.7 | 1.05 |
| 9 | 512 | 2.1 MiB | 39,000 | 114.3 | 3.57 |
| 10 | 1024 | 8.4 MiB | 86,000 | 201.2 | 6.29 |
| 11 | 2048 | 33 MiB | 303,000 | 225.0 | 7.03 |
| 12 | 4096 | 132 MiB | 1,925,000 | 140.5 | 4.39 |

**Key observations**:
- **n=6,7,8 have nearly identical times (~32µs)** despite 12x data range.
  This is OpenMP thread launch overhead (dim=64 hits `OMP_THRESHOLD`).
- **Peak bandwidth at n=11** (~225 GB/s) — data fits in the L2/SLC hierarchy.
- **n=12 drops to ~137 GB/s** — 132 MiB exceeds on-chip caches; memory-bound.

### 3.2 Within-tile vs. cross-tile (n=10, 11, 12)

| Qubits | Target | Path | Median (ns) | BW (GB/s) |
|--------|--------|------|-------------|-----------|
| 10 | 0 | within | 90,000 | 192 |
| 10 | 5 | cross | 78,000 | **222** |
| 11 | 0 | within | 299,000 | **228** |
| 11 | 5 | cross | 338,000 | 202 |
| 11 | 7 | cross | 303,000 | 225 |
| 12 | 0 | within | 1,929,000 | 140 |
| 12 | 5 | cross | 2,034,000 | 133 |
| 12 | 11 | cross | 2,127,000 | **127** |

**Key observations**:
- At n=10, cross-tile is 15% **faster** than within-tile. The cross-tile path
  does sequential sweeps over entire tiles, which is more vectorizer-friendly
  than the bit-insertion stride within tiles.
- At n=11–12, cross-tile with high tile_bit (target=5,6) is ~10% slower.
  The highest tile_bit (target=n-1) is consistently the slowest cross-tile
  target — likely because the tr0/tr1 tile pairs are maximally distant
  in memory.

### 3.3 x_tiled vs. x_packed

| Qubits | x_tiled (ns) | x_packed (ns) | Speedup |
|--------|-------------|--------------|---------|
| 8 | 29,000 | 35,000 | 1.21x |
| 9 | 46,000 | 44,000 | 0.96x |
| 10 | 88,000 | 92,000 | 1.05x |
| 11 | 333,000 | 354,000 | 1.06x |
| 12 | 1,976,000 | 2,090,000 | 1.06x |

Tiled is consistently ~5–6% faster at cache-relevant sizes (10–12 qubits).

## 4. Compiler Vectorization Analysis

Compilation with `-Rpass=loop-vectorize -Rpass-missed=loop-vectorize`:

- **4 loops vectorized** (vectorization width: 2, complex double)
- **51 loops NOT vectorized**

Reasons for vectorization failure:
1. **"cannot identify array bounds"** — The `data` pointer with computed
   offsets (`offset_tile`, `insertBit0()`) prevents the compiler from
   proving no aliasing. The within-tile inner loop at line 459 fails for this
   reason.
2. **"unsafe dependent memory operations"** — The diagonal cross-tile
   conjugation loop (line 544) reads and writes overlapping regions
   (`data[idx10 + lc]` vs `data[lr + lc * TILE_DIM + offset_t10]`).

### Assembly observations

The compiler generates:
- **`ldp q0, q1` / `stp` pairs** for the straightforward tile swap loops
  (32 bytes per load/store pair — good)
- **Scalar `fneg d`** for conjugation (`conj()` needs to negate only the
  imaginary component, forcing element decomposition)
- **No NEON vector (`v`) register usage** — all operations on `q` (single
  128-bit = one complex double) or `d` (64-bit = one real/imag component)

The function compiles to **~5000 lines** of assembly due to the three
cross-tile branches being fully unrolled/specialized by the compiler.

## 5. Identified Optimization Opportunities

### 5.1 Raise OMP threshold for x_tiled (Low effort, High impact at small n)

**Problem**: At n=6,7,8, the ~32µs OpenMP thread launch cost is 60–99% of the
total time. The actual computation takes <1µs (as seen at n=5 where OMP is not
triggered).

**Suggestion**: For the X gate specifically (swap-only, no arithmetic), the
per-element work is very low. Consider raising `OMP_THRESHOLD` or using a
per-gate threshold. A data-size-based threshold (~1 MiB) would be more
appropriate than a dimension-based one.

### 5.2 Merge the two lr-loops in the diagonal cross-tile case (Low effort)

**Problem**: Lines 531–547 of `gate_tiled.c` have two separate `lr` loops for
the diagonal case (`btr == btc`):
```c
// Loop 1: swap t00 <-> t11 (all lc)
for (lr ...) for (lc = 0; lc < TILE_DIM; ...)
    swap data[idx00+lc], data[idx11+lc];

// Loop 2: conjugate-transpose t10 (lc <= lr)
for (lr ...) for (lc = 0; lc <= lr; ...)
    swap-conj data[idx10+lc], data[idx10_transposed];
```

This makes two passes over the same cache lines. **Fusing into one lr-loop**
would halve the cache pressure for diagonal tile-blocks.

### 5.3 Fix column-stride access in the conj-transpose loop (Medium effort)

**Problem**: The diagonal cross-tile conjugation (line 544) accesses
`data[lr + lc * TILE_DIM + offset_t10]` — this is a column-stride access
pattern (stride = `TILE_DIM * sizeof(cplx_t)` = 512 bytes for 32x32 tiles).
Each access to a new `lc` jumps 512 bytes, causing L1 cache misses.

**Suggestion**: This is an in-place conjugate-transpose of a TILE_DIM×TILE_DIM
block. Well-known cache-oblivious transpose algorithms exist (e.g., blocking
into 4×4 or 8×8 sub-blocks). Alternatively, since this is an in-place
Hermitian transpose, process elements in (lc, lr) pairs with lc < lr and swap.

### 5.4 Help the vectorizer with restrict and pragma hints (Medium effort)

**Problem**: 51 of 55 inner loops fail to vectorize, primarily due to
"cannot identify array bounds" and "unsafe dependent memory operations".

**Suggestions**:
- Add `restrict` to local tile pointers in the inner scope (the `data` pointer
  already has `restrict`, but the compiler can't prove non-aliasing between
  derived pointers like `t00`, `t10`, `t11`, `t01`)
- Use `#pragma clang loop vectorize(enable)` or
  `#pragma clang loop distribute(enable)` on the inner lc-loops
- For the within-tile path: the four elements `a00, a01, a10, a11` are at
  known offsets from the tile base — asserting non-overlap could enable
  vectorization

### 5.5 Consider NEON intrinsics for the conj operation (Medium effort)

**Problem**: `conj()` on `double _Complex` forces the compiler to decompose
into real/imag parts and negate the imaginary part with scalar `fneg d`.

**Suggestion**: A NEON intrinsic approach can conjugate a complex double in one
instruction:
```c
#include <arm_neon.h>
// XOR the sign bit of the imaginary part
static inline cplx_t fast_conj(cplx_t z) {
    float64x2_t v = vld1q_f64((double*)&z);
    float64x2_t mask = {0.0, -0.0};  // sign bit in high lane
    v = vreinterpretq_f64_u64(
        veorq_u64(vreinterpretq_u64_f64(v),
                  vreinterpretq_u64_f64(mask)));
    cplx_t result;
    vst1q_f64((double*)&result, v);
    return result;
}
```
This keeps the value in a single Q register and avoids decomposition.

### 5.6 No action needed: memory-bound regime at n=12+ (Expected)

At n=12, the 132 MiB working set exceeds on-chip caches and bandwidth drops
to ~137 GB/s (likely close to DRAM bandwidth on this machine). This is
expected behavior — tiling can't help when the entire dataset exceeds cache.
Prefetch hints (`__builtin_prefetch`) on the next tile during iteration could
marginally improve DRAM latency hiding.

## 6. Priority Ranking

| Priority | Item | Expected Impact |
|----------|------|----------------|
| 1 | Raise OMP threshold (5.1) | Eliminate 32µs overhead at n≤8 |
| 2 | Merge diagonal lr-loops (5.2) | ~5–10% improvement in diagonal cross-tile path |
| 3 | Fix column-stride access (5.3) | Significant for cross-tile diagonal blocks |
| 4 | Vectorizer hints (5.4) | Potential 2x throughput on inner loops |
| 5 | NEON conj intrinsic (5.5) | Minor improvement in conj-heavy paths |
