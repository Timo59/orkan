# Mixed State Gate Performance Analysis

This document summarizes performance analysis of mixed-state gate operations in `mhipster.c`, including memory layout considerations, profiling results, and optimization recommendations.

## Executive Summary

The current packed lower-triangular storage with BLAS-1 operations achieves **51-56 GB/s average bandwidth** on Apple M4, which is approximately 50% of peak memory bandwidth. Three distinct bottleneck regimes exist depending on target qubit:

| Regime | Target Qubits | Bottleneck | Achieved Bandwidth |
|--------|---------------|------------|-------------------|
| Low | 0-3 | BLAS call overhead (millions of tiny calls) | 17-72 GB/s |
| Mid | 4-5 | Balanced (sweet spot) | **90-93 GB/s** |
| High | 6+ | CONJ_OP variable-stride access | 58-79 GB/s |

**Key finding**: Blocked memory layout is **not recommended** for this codebase. The current layout is near-optimal for QAOA workloads where diagonal gates (RZZ) dominate.

---

## 1. Current Memory Layout

### Packed Lower-Triangular Storage

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

**Index formula**: `(i, j) → j*N - j*(j-1)/2 + (i-j)` for i ≥ j, where N = 2ⁿ

### Gate Traversal Structure

Single-qubit gates use the `TRAVERSE_MIXED_1Q` macro with three operation types:

1. **DIAG_OP**: Process (0,0)↔(1,1) quadrant pairs (contiguous, uses `cblas_zswap`)
2. **CONJ_OP**: Process (1,0)↔(0,1) pairs where (0,1) is in upper triangle (variable stride, scalar loop)
3. **OFFDIAG_OP**: Process (1,0)↔(0,1) pairs both in lower triangle (contiguous, uses `cblas_zswap`)

---

## 2. Profiling Results

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

## 3. Blocked Memory Layout Analysis

### Concept

Reorganize density matrix into cache-friendly blocks where elements processed together are stored together:

```
Block(i,j) contains elements with row indices in [i*B, (i+1)*B)
                          and col indices in [j*B, (j+1)*B)
```

### Potential Benefits

| Aspect | Benefit |
|--------|---------|
| Cache locality | Elements in same block fit in L1/L2 |
| BLAS-3 for Hadamard | Enable matrix-matrix operations |
| Parallelization | Independent blocks → embarrassingly parallel |
| SIMD for CONJ_OP | Fixed stride within blocks |

### Why Blocking is NOT Recommended

| Factor | Impact |
|--------|--------|
| **BLAS-1 efficiency** | Blocking fragments long vectors into short ones, *reducing* BLAS efficiency |
| **Pauli gates** | X, Y, Z only need BLAS-1; blocking adds overhead without enabling BLAS-3 |
| **QAOA workload** | Dominated by RZZ (diagonal) which doesn't benefit from blocking |
| **Implementation complexity** | Two code paths, conversion overhead, index complexity |
| **Memory overhead** | Block offset tables, potential padding |

### Blocking Trade-off Summary

| Gate Type | Current Efficiency | Blocking Benefit |
|-----------|-------------------|------------------|
| Z, S, T (diagonal) | High | None |
| X, Y (permutation) | Medium | Small |
| H, Rx, Ry, Rz (mixing) | Lower | Moderate (BLAS-3) |

**Verdict**: Only worthwhile if workload is dominated by Hadamard/rotation gates AND you can stay in blocked layout for entire circuits.

---

## 4. CONJ_OP Deep Dive

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
| Scalar (current) | 0.170 ms | 1× |
| Precomputed indices | 0.014 ms | 12× |
| Unrolled + prefetch | 0.012 ms | 15× |
| Blocked (SIMD) | 0.003 ms | **50×** |

### Why CONJ_OP is Slow

- Achieves only **1.5-11 GB/s** vs 100+ GB/s theoretical
- Variable stride defeats prefetcher
- No SIMD vectorization possible
- Each iteration has data dependency on previous `step`

### Impact Assessment

```
CONJ_OP contribution:
- Overall (all targets): ~3% of total time
- Highest target only: ~70% of that target's time
- Fixing CONJ_OP: ~1.03× overall, ~2× for highest qubit only
```

---

## 5. Optimization Recommendations

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

**Expected improvement**: 2× bandwidth for targets 0-3.

### Priority 2: CONJ_OP Optimization (Medium Impact)

**Option A**: Precompute step values to enable loop unrolling:

```c
static void op_swap_conj_unrolled(cplx_t *d10, dim_t count, dim_t start_k,
                                   dim_t dim, dim_t col_block) {
    dim_t step = 0;
    dim_t k = start_k;

    // Unroll by 4 with prefetch
    for (; k + 4 <= start_k + count; k += 4) {
        __builtin_prefetch(&d10[k + 8], 1, 3);

        dim_t s0 = step; step += dim - (k + col_block + 2);
        dim_t s1 = step; step += dim - (k + col_block + 3);
        dim_t s2 = step; step += dim - (k + col_block + 4);
        dim_t s3 = step; step += dim - (k + col_block + 5);

        // Process 4 pairs
        cplx_t t0 = d10[k-start_k+0], t1 = d10[k-start_k+1];
        cplx_t t2 = d10[k-start_k+2], t3 = d10[k-start_k+3];

        d10[k-start_k+0] = conj(d10[k-start_k+0+s0]);
        d10[k-start_k+1] = conj(d10[k-start_k+1+s1]);
        d10[k-start_k+2] = conj(d10[k-start_k+2+s2]);
        d10[k-start_k+3] = conj(d10[k-start_k+3+s3]);

        d10[k-start_k+0+s0] = conj(t0);
        d10[k-start_k+1+s1] = conj(t1);
        d10[k-start_k+2+s2] = conj(t2);
        d10[k-start_k+3+s3] = conj(t3);
    }
    // ... remainder loop
}
```

**Expected improvement**: 2-3× for CONJ_OP regions.

**Option B**: For highest 1-2 qubits only, consider blocked approach for CONJ_OP.

### Priority 3: OpenMP Tuning (Low Impact)

Current parallelization is at `col_block` level. For highest qubits, there's only 1 col_block (no parallelism).

**Solution**: Add inner parallelization for high qubits:

```c
if (target >= n_qubits - 2) {
    // Parallelize at finer granularity for highest qubits
    #pragma omp parallel for
    for (dim_t col = col_block; col < col_block + incr; ++col) {
        // ... process column
    }
}
```

### Priority 4: Memory Alignment (Minor)

Ensure 64-byte alignment for cache line optimization:

```c
state->data = aligned_alloc(64, size * sizeof(cplx_t));
```

(Already done in current implementation.)

---

## 6. Recommendations by Use Case

### QAOA / VQE (Diagonal-Dominated)

**Recommendation**: Keep current implementation.

- RZZ gates are diagonal (no CONJ_OP)
- Current bandwidth is good (50-90 GB/s)
- Blocked layout provides no benefit

### General Circuits (Mixed Gates)

**Recommendation**: Implement Priority 1 + 2 optimizations.

- Low-qubit batching: +30% overall
- CONJ_OP unrolling: +5% overall
- Total expected improvement: ~35%

### Hadamard-Heavy Circuits

**Recommendation**: Consider blocked layout only if:
1. >50% of gates are Hadamard/rotations
2. Can stay in blocked layout for entire circuit
3. n ≥ 8 qubits (blocking overhead amortized)

---

## 7. Benchmark Commands

```bash
# Build benchmarks
clang -O3 -march=native test/bench_bandwidth.c -o test/bench_bandwidth -framework Accelerate
clang -O3 -march=native test/bench_gate_bw.c -o test/bench_gate_bw -framework Accelerate
clang -O3 -march=native test/bench_conj_op.c -o test/bench_conj_op

# Run benchmarks
./test/bench_bandwidth 22 30      # Memory bandwidth baseline (64MB arrays)
./test/bench_gate_bw 12 200       # Gate bandwidth by target qubit
./test/bench_conj_op 14 200       # CONJ_OP micro-benchmark
./test/analyze_gate_calls 12      # BLAS call pattern analysis
```

---

## 8. Summary

| Metric | Current | With Optimizations | Notes |
|--------|---------|-------------------|-------|
| Avg bandwidth | 51 GB/s | ~70 GB/s | +35% |
| Low-qubit bandwidth | 17-72 GB/s | 50-90 GB/s | BLAS batching |
| High-qubit bandwidth | 58-79 GB/s | 70-85 GB/s | CONJ_OP unrolling |
| Peak bandwidth | 93 GB/s | 93 GB/s | Already optimal |
| CONJ_OP (isolated) | 1.5 GB/s | 50+ GB/s | With blocking/SIMD |

**Bottom line**: The current implementation is well-designed. Targeted optimizations (BLAS batching for low qubits, CONJ_OP unrolling for high qubits) can improve overall performance by ~35% without architectural changes.
