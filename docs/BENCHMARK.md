# Mixed State Gate Performance Analysis

This document summarizes performance analysis of mixed-state gate operations, including memory layout considerations, profiling results, and optimization recommendations.

The two mixed-state implementations are:
- **Packed** (`gate_packed.c`): LAPACK packed lower-triangular column-major storage
- **Tiled** (`gate_tiled.c`): Cache-friendly tiled lower-triangular storage

## Executive Summary

The packed implementation with BLAS-1 operations achieves **51-56 GB/s average bandwidth** on Apple M4, which is approximately 50% of peak memory bandwidth. Three distinct bottleneck regimes exist depending on target qubit:

| Regime | Target Qubits | Bottleneck | Achieved Bandwidth |
|--------|---------------|------------|-------------------|
| Low | 0-3 | BLAS call overhead (millions of tiny calls) | 17-72 GB/s |
| Mid | 4-5 | Balanced (sweet spot) | **90-93 GB/s** |
| High | 6+ | CONJ_OP variable-stride access | 58-79 GB/s |

The tiled implementation (`MIXED_TILED`) addresses the CONJ_OP bottleneck by reorganizing data into fixed-size tiles, enabling better cache locality and SIMD vectorization for high target qubits.

---

## 1. Memory Layouts

### Packed Lower-Triangular Storage (`MIXED_PACKED`)

Density matrices are stored in LAPACK packed lower-triangular column-major format:

```
Matrix (N=4):                  Linear storage:
┌──────────────────────┐
│ ρ00                  │       [0]: ρ00  ┐
│ ρ10  ρ11             │       [1]: ρ10  │ Column 0
│ ρ20  ρ21  ρ22        │       [2]: ρ20  │
│ ρ30  ρ31  ρ32  ρ33   │       [3]: ρ30  ┘
└──────────────────────┘       [4]: ρ11  ┐
                               [5]: ρ21  │ Column 1
                               [6]: ρ31  ┘
                               ...
```

**Index formula**: `(i, j) → j*N - j*(j-1)/2 + (i-j)` for i >= j, where N = 2^n

### Tiled Lower-Triangular Storage (`MIXED_TILED`)

The lower triangle is partitioned into `TILE_DIM x TILE_DIM` tiles (default TILE_DIM = 32). Tiles are stored in row-major order within each tile, and tile pairs are indexed as `(tile_row, tile_col)` with `tile_row >= tile_col`.

```
Tile grid (N=64, TILE_DIM=32):

     col 0    col 1
   ┌────────┬────────┐
   │ T(0,0) │        │    Storage: T(0,0), T(1,0), T(1,1)
   ├────────┼────────┤    Each tile: TILE_DIM² contiguous elements
   │ T(1,0) │ T(1,1) │
   └────────┴────────┘
```

**Storage**: `n_tiles * (n_tiles + 1) / 2 * TILE_SIZE` elements, where `n_tiles = ceil(2^n / TILE_DIM)`.

**WARNING**: Intra-tile layout is row-major, NOT column-major like LAPACK packed. Do NOT pass tile pointers to LAPACK routines expecting column-major storage.

### Gate Traversal Structure (Packed)

Single-qubit gates in `gate_packed.c` use the `TRAVERSE_MIXED_1Q` macro with three operation types:

1. **DIAG_OP**: Process (0,0)↔(1,1) quadrant pairs (contiguous, uses `cblas_zswap`)
2. **CONJ_OP**: Process (1,0)↔(0,1) pairs where (0,1) is in upper triangle (variable stride, scalar loop)
3. **OFFDIAG_OP**: Process (1,0)↔(0,1) pairs both in lower triangle (contiguous, uses `cblas_zswap`)

---

## 2. Profiling Results (Packed)

### Test Environment
- **CPU**: Apple M4
- **Theoretical bandwidth**: ~120 GB/s
- **Test sizes**: 10-14 qubits (8 MB - 2 GB packed arrays)

### Overall Time Distribution (Pauli-X gate, all targets)

| n_qubits | DIAG_OP | CONJ_OP | OFFDIAG_OP |
|----------|---------|---------|------------|
| 10 | 51.1% | 2.2% | 46.6% |
| 12 | 51.4% | 2.6% | 46.0% |
| 14 | 49.8% | 3.5% | 46.7% |

**CONJ_OP is only ~3% of total time** when averaged across all target qubits.

### Per-Target Breakdown

CONJ_OP becomes significant only for the highest target qubits:

| Target (n=14) | CONJ % of that target's time |
|---------------|------------------------------|
| 0-5 | <1% |
| 10 | 6.5% |
| 11 | 15.2% |
| 12 | 35.6% |
| 13 (highest) | **68.5%** |

### BLAS Call Analysis (n=12 qubits)

| Target | incr | Total BLAS calls | Avg length | CONJ pairs |
|--------|------|------------------|------------|------------|
| 0 | 1 | 4,194,304 | 1.0 | 2,048 |
| 1 | 2 | 2,097,152 | 2.0 | 3,072 |
| 4 | 16 | 262,144 | 15.9 | 17,408 |
| 7 | 128 | 32,768 | 124.0 | 132,096 |
| 11 | 2048 | 2,048 | 1024.5 | 2,098,176 |

### Bandwidth by Target Qubit (n=12)

```
Target:   0    1    2    3    4    5    6    7    8    9   10   11
BW GB/s: 17   32   49   72   92   93   59   58   79   76   73   66
         └─── BLAS overhead ───┘   └ peak ┘   └── CONJ_OP overhead ──┘
```

---

## 3. Tiled Layout Analysis

### Design

The tiled layout reorganizes the density matrix into fixed-size tiles where elements processed together are stored contiguously:

```
Tile(i,j) contains elements with row indices in [i*TILE_DIM, (i+1)*TILE_DIM)
                              and col indices in [j*TILE_DIM, (j+1)*TILE_DIM)
```

### Benefits

| Aspect | Benefit |
|--------|---------|
| Cache locality | Elements in same tile fit in L1/L2 |
| CONJ_OP elimination | Fixed stride within tiles avoids the variable-stride problem |
| Parallelization | Independent tiles are embarrassingly parallel |
| SIMD friendliness | Contiguous intra-tile data enables vectorization |

### Trade-offs

| Factor | Impact |
|--------|--------|
| **BLAS-1 vectors** | Tiling fragments long contiguous vectors into tile-sized chunks |
| **Pauli gates** | X, Y, Z only need swaps; tiling adds indexing overhead for small systems |
| **Diagonal gates** | RZZ and phase gates see no benefit from tiling |
| **Memory overhead** | Padding for non-power-of-TILE_DIM dimensions |

### When to Prefer Each Layout

| Gate Type | Packed | Tiled |
|-----------|--------|-------|
| Z, S, T, Rz (diagonal) | Preferred | No benefit |
| X, Y (permutation) | Good for mid qubits | Better for high qubits |
| H, Rx, Ry (mixing) | CONJ_OP bottleneck at high qubits | Avoids CONJ_OP |
| QAOA circuits (RZZ-dominated) | Preferred | No benefit |
| General circuits (mixed gates) | Good | Better cache behavior at scale |

---

## 4. CONJ_OP Deep Dive (Packed)

### The Problem

```c
// Variable stride that grows each iteration - terrible for SIMD/prefetcher
for (dim_t k = start_k; k < start_k + count; ++k) {
    cplx_t tmp = d10[k];
    d10[k] = conj(d10[k + step]);
    d10[k + step] = conj(tmp);
    step += dim - (k + col_block + 2);  // Growing step!
}
```

### Micro-benchmark Results (n=14, target=13)

| Implementation | Time | Speedup |
|----------------|------|---------|
| Scalar (packed) | 0.170 ms | 1x |
| Precomputed indices | 0.014 ms | 12x |
| Unrolled + prefetch | 0.012 ms | 15x |
| Tiled (SIMD) | 0.003 ms | **50x** |

### Why CONJ_OP is Slow in Packed

- Achieves only **1.5-11 GB/s** vs 100+ GB/s theoretical
- Variable stride defeats prefetcher
- No SIMD vectorization possible
- Each iteration has data dependency on previous `step`

### Impact Assessment

```
CONJ_OP contribution (packed layout):
- Overall (all targets): ~3% of total time
- Highest target only: ~70% of that target's time
- Fixing CONJ_OP: ~1.03x overall, ~2x for highest qubit only
```

The tiled layout eliminates this problem entirely by ensuring contiguous access within tiles.

---

## 5. Optimization Recommendations (Packed)

### Priority 1: Low-Qubit BLAS Batching (High Impact)

**Problem**: Targets 0-3 make millions of single-element BLAS calls.

**Solution**: Skip BLAS for tiny operations, use direct SIMD:

```c
#define BLAS_THRESHOLD 16

static inline void swap_elements(cplx_t *a, cplx_t *b, dim_t n) {
    if (n >= BLAS_THRESHOLD) {
        cblas_zswap((int)n, a, 1, b, 1);
    } else {
        // Direct swap - avoid BLAS call overhead
        for (dim_t i = 0; i < n; i++) {
            cplx_t tmp = a[i];
            a[i] = b[i];
            b[i] = tmp;
        }
    }
}
```

**Expected improvement**: 2x bandwidth for targets 0-3.

### Priority 2: CONJ_OP Optimization (Medium Impact)

**Option A**: Precompute step values to enable loop unrolling.

**Option B**: For highest 1-2 qubits only, use the tiled layout. The three-way dispatch in `gate.c` makes it possible to select storage format per-circuit based on gate composition.

### Priority 3: OpenMP Tuning (Low Impact)

Current parallelization is at `col_block` level. For highest qubits, there's only 1 col_block (no parallelism).

### Priority 4: Memory Alignment (Minor)

Ensure 64-byte alignment for cache line optimization. (Already done in current implementation.)

---

## 6. Recommendations by Use Case

### QAOA / VQE (Diagonal-Dominated)

**Recommendation**: Use `MIXED_PACKED`.

- RZZ gates are diagonal (no CONJ_OP)
- Current bandwidth is good (50-90 GB/s)
- Tiled layout provides no benefit for diagonal gates

### General Circuits (Mixed Gates)

**Recommendation**: Benchmark both layouts; use `MIXED_TILED` for n >= 8 qubits.

- Tiled avoids CONJ_OP overhead for high target qubits
- Packed is faster for small systems due to lower indexing overhead
- The `bench_mixed` tool compares both layouts head-to-head

### Hadamard-Heavy Circuits

**Recommendation**: Use `MIXED_TILED`.

- Hadamard/rotation gates trigger CONJ_OP at high target qubits
- Tiled layout eliminates the variable-stride bottleneck
- Benefit increases with qubit count

---

## 7. Benchmark Commands

```bash
# Build (from project root)
mkdir -p cmake-build-debug && cd cmake-build-debug
cmake ..
cmake --build . --target bench_mixed

# Run comparative benchmark
./benchmark/bench_mixed                              # Default: 2-12 qubits, 1000 iterations
./benchmark/bench_mixed --min-qubits 4 --max-qubits 14 --iterations 100
./benchmark/bench_mixed --csv                        # CSV output for analysis
./benchmark/bench_mixed --pgfplots                   # pgfplots-compatible .dat output

# Legacy micro-benchmarks (standalone, not CMake-managed)
clang -O3 -march=native test/bench_bandwidth.c -o test/bench_bandwidth -framework Accelerate
clang -O3 -march=native test/bench_gate_bw.c -o test/bench_gate_bw -framework Accelerate
clang -O3 -march=native test/bench_conj_op.c -o test/bench_conj_op
```

---

## 8. Summary

| Metric | Packed | Tiled | Notes |
|--------|--------|-------|-------|
| Avg bandwidth | 51 GB/s | TBD | Packed profiled on M4 |
| Low-qubit bandwidth | 17-72 GB/s | Similar | Both limited by call overhead |
| Mid-qubit bandwidth | 90-93 GB/s | Similar | Both near-optimal |
| High-qubit bandwidth | 58-79 GB/s | Expected higher | Tiled eliminates CONJ_OP |
| CONJ_OP (isolated) | 1.5 GB/s | N/A | Tiled avoids variable-stride |
| Memory overhead | Minimal | +padding | Tiled pads to TILE_DIM boundary |

**Bottom line**: Both layouts are available via three-way dispatch (`PURE`, `MIXED_PACKED`, `MIXED_TILED`). Use packed for diagonal-dominated workloads (QAOA) and tiled for general circuits with many Hadamard/rotation gates. Run `bench_mixed --pgfplots` to generate comparative data for your target hardware.
