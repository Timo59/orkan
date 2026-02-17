# Mixed State Benchmarks

Benchmarks comparing qlib's packed sparse density matrix implementation against BLAS dense, QuEST, and Qulacs.

## How to Run

```bash
# From project root
cd cmake-build-debug
cmake --build . --target bench_mixed
./benchmark/bench_mixed
```

### Options

```
--min-qubits N   Minimum qubit count (default: 2)
--max-qubits N   Maximum qubit count (default: 12)
--step N         Qubit count increment (default: 1)
--iterations N   Number of iterations (default: 1000)
--warmup N       Warm-up iterations (default: 100)
--csv            Output CSV format
--pgfplots       Output pgfplots-compatible .dat format
--verbose        Show per-method memory usage inline
--help           Show help
```

### Examples

```bash
./benchmark/bench_mixed --max-qubits 8           # Quick test (every qubit: 2,3,4,...,8)
./benchmark/bench_mixed --step 2                 # Even qubits only (2,4,6,...,12)
./benchmark/bench_mixed --csv > results.csv      # Export CSV data
./benchmark/bench_mixed --pgfplots > results.dat # Export for pgfplots/LaTeX
./benchmark/bench_mixed --iterations 5000        # Higher precision
```

## Results (10 qubits, 1024-dim)

| Method | X gate | H gate | Z gate |
|--------|--------|--------|--------|
| qlib packed | 130 ms | 125 ms | 79 ms |
| qlib blocked | 104 ms | 179 ms | 77 ms |
| QuEST | 524 ms | 520 ms | 392 ms |
| Qulacs | 549 ms | 548 ms | 552 ms |

**Speedup (packed vs others):** qlib is 4-5x faster than QuEST and Qulacs.

**Blocked vs Packed:** X gate 1.26x faster, Z gate 1.03x faster at 10 qubits.

**Memory (10 qubits):** packed=8.0 MB, blocked=8.5 MB, QuEST=16.9 MB, Qulacs=16.0 MB

### Blocked Format Performance by Qubit Count

| Qubits | X blocked vs packed | Z blocked vs packed |
|--------|---------------------|---------------------|
| 6 | **13.9x faster** | **22.9x faster** |
| 8 | 0.86x (near parity) | 0.91x (near parity) |
| 10 | **1.32x faster** | **1.08x faster** |

The blocked format excels at small sizes (avoiding BLAS overhead) and at larger sizes (cache-friendly tile structure with SIMD).

## Benchmark Methodology

### Timing Strategy

Each benchmark measures the time to apply a gate to every qubit in the system:

1. **Warm-up phase**: Run `warmup` iterations (default: 100) to stabilize caches and JIT
2. **Timed phase**: Run `iterations` iterations (default: 1000), applying the gate to each qubit
3. **Metrics**: Total time and ops/sec (one op = one gate application)

The benchmark applies gates to all target qubits per iteration to measure realistic workloads where gates are applied across the register.

### Memory Strategy

Memory is measured as the RSS (resident set size) delta capturing the full allocation from density matrix initialization through gate simulation setup:

1. Record RSS before any allocation
2. Initialize the density matrix (random values for benchmarking)
3. Build any required gate matrices (for dense methods, this includes full N×N unitary for each target qubit)
4. Record RSS after setup
5. Report the delta as actual memory used

This captures real memory usage including:
- The density matrix storage
- Gate matrices and temporary buffers
- Framework-internal allocations

### What We Compare

- **Time**: Wall-clock time for gate applications (lower is better)
- **Memory**: Actual RSS used for state representation and operations (lower is better)
- **Speedup**: Ratio of other methods' time to qlib time
- **Memory savings**: Percentage reduction vs dense representation

## Why Qulacs Peaks at 6 Qubits

Qulacs performs well relative to qlib at 6 qubits but falls behind at all other tested sizes (2, 4, 8, 10, 12). The cause is cache fitting combined with data volume differences.

**Root cause:** Qulacs' `dm_single_qubit_dense_matrix_gate` iterates over the full N×N density matrix (row-major), while qlib's packed format stores only the lower triangle (~N(N+1)/2 elements). At 6 qubits the full matrix fits in L1 cache, neutralizing this 2× data advantage.

| Qubits | dim | Qulacs matrix size | Why Qulacs loses |
|--------|-----|-------------------|------------------|
| 2, 4 | 4–16 | 256 B – 4 KB | Fork/pipe overhead (~1-5 ms) dwarfs sub-µs gate work |
| **6** | **64** | **64 KB** | **Fits L1 cache; simple row-major loop runs at full bandwidth** |
| 8 | 256 | 1 MB | Exceeds L1; `insert_zero_to_basis_index` strides cause cache thrashing |
| 10 | 1024 | 16 MB | L3/main memory; OpenMP enables at dim ≥ 1024 adding thread overhead |
| 12 | 4096 | 256 MB | Main memory only; 2× data volume over packed compounds with cache misses |

**Key code references:**
- Qulacs OpenMP threshold: `set_qulacs_num_threads(dim, 10)` enables threads when `dim >= 1024` (qulacs `utility.cpp:134`)
- Qulacs dm gate loop: `(dim/2)²` iterations × 4 elements = full dim² traversal (`update_ops_dm.cpp:59-97`)
- Qulacs benchmark forks per call: `bench_qulacs.cpp:146` — constant process creation cost penalizes small sizes

## Comparison Methods

| Method | Description | Complexity |
|--------|-------------|------------|
| qlib packed | Packed lower-triangular with stride traversal (OpenMP) | O(N) per gate |
| qlib blocked | Tile-based blocking with NEON SIMD (OpenMP) | O(N) per gate |
| BLAS dense | Full matrix UρU† via zgemm | O(N³) |
| naive loop | Element-by-element UρU† computation | O(N³) |
| QuEST | Multithreaded quantum simulator (OpenMP) | O(N²) |
| Qulacs | High-performance simulator with SIMD/OpenMP (QunaSys) | O(N²) |

Note: BLAS dense and naive loop are only run for ≤8 qubits due to O(N³) performance.

Memory is measured as actual runtime RSS (resident set size) delta, not theoretical sizes. This captures real allocations including temporary buffers.

### Memory Measurement Methodology

Memory is measured as RSS (resident set size) delta capturing the full allocation cost.

**Process isolation**: QuEST and Qulacs benchmarks run in forked child processes to ensure accurate memory measurement. This avoids memory reuse between benchmarks and captures the true cost of each framework's initialization.

**Full initialization**: Each framework's memory measurement includes all setup:
- **QuEST**: Measures `initQuESTEnv()` + `createDensityQureg()` + state initialization (forked process, RSS delta)
- **Qulacs**: Reports theoretical matrix size (`dim² × 16 bytes`) since RSS delta after fork() is unreliable on macOS
- **BLAS/naive**: Measures density matrix + gate matrix allocations (RSS delta)
- **qlib**: Measures packed storage allocation (RSS delta)

**Page-touching**: To avoid lazy allocation where the OS defers page mapping:
- **QuEST**: Uses `initDebugState()` to touch all state vector pages
- **Qulacs**: Uses `memset()` to initialize matrix (memory reported as theoretical size)
- **BLAS/naive/qlib**: Fill with random values, naturally touching all pages

**Theoretical sizes** (for reference):
- Dense: `2^n × 2^n × 16 bytes` (complex double)
- Packed (qlib): `2^n × (2^n + 1) / 2 × 16 bytes` (lower triangle)

## Scientific Rigor Considerations

When citing these benchmarks, document the following:

### Environment
- CPU model, core count, clock speed
- RAM size and speed
- OS and kernel version
- Compiler and version (e.g., clang 15.0, gcc 13.2)
- BLAS implementation (OpenBLAS, Accelerate, MKL)
- Library versions: QuEST, Qulacs

### Build Configuration
- Optimization level (-O0, -O2, -O3)
- Architecture flags (-march=native)
- Debug vs Release build
- Link-time optimization (LTO) enabled?

### Statistical Validity
- Number of independent runs (current: single run)
- Variance and standard deviation not reported
- No confidence intervals
- Results may vary between runs due to system load

### Timing
- Clock source: `CLOCK_MONOTONIC` via `clock_gettime()`
- Timer resolution: nanosecond precision
- Wall-clock time (includes potential OS scheduling noise)

### Fairness Considerations
- QuEST uses OpenMP multithreading; thread count affects results
- Qulacs uses OpenMP for parallelization (enabled by default)
- BLAS may use multiple threads depending on implementation
- qlib uses OpenMP parallelization (enabled by default via `ENABLE_OPENMP` CMake option)
- All frameworks now have comparable multithreading capabilities

### Limitations
- Only single-qubit gates benchmarked
- No multi-qubit gate operations (CNOT, CZ)
- No measurement operations
- Random state initialization (not physically meaningful states)
- Does not measure gate matrix construction time separately
- Memory measurement granularity limited by OS page size (typically 4KB/16KB)

## External Frameworks

QuEST and Qulacs are included as git submodules in `extern/`. They are auto-detected during cmake configuration. Check the cmake output for status:

```
-- Benchmark configuration:
--   QuEST:     TRUE
--   Qulacs:    TRUE
--   OpenMP:    TRUE
```

### Setup (if submodules not initialized)

```bash
git submodule update --init extern/QuEST extern/qulacs
brew install libomp boost  # macOS dependencies

# Build QuEST
cd extern/QuEST && mkdir build && cd build
OpenMP_ROOT=$(brew --prefix)/opt/libomp cmake ..
cmake --build .

# Build Qulacs (csim_static library)
# IMPORTANT: Use -DCMAKE_BUILD_TYPE=Release for optimized builds (-O3)
cd extern/qulacs && mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DUSE_PYTHON=No -DUSE_TEST=No -DUSE_OMP=Yes \
  -DOpenMP_C_FLAGS="-Xclang -fopenmp -I$(brew --prefix libomp)/include" \
  -DOpenMP_CXX_FLAGS="-Xclang -fopenmp -I$(brew --prefix libomp)/include" \
  -DOpenMP_C_LIB_NAMES="omp" -DOpenMP_CXX_LIB_NAMES="omp" \
  -DOpenMP_omp_LIBRARY="$(brew --prefix libomp)/lib/libomp.dylib"
cmake --build . --target csim_static
```

### Custom paths

If frameworks are installed elsewhere, pass hints to cmake:

```bash
cmake -DQuEST_DIR=/path/to/quest -DQULACS_DIR=/path/to/qulacs ..
```

## Files

```
benchmark/
├── CMakeLists.txt
├── include/bench.h      # Timing utilities, result types
└── src/
    ├── bench_mixed.c    # Main benchmark + BLAS/naive implementations
    ├── bench_quest.c    # QuEST wrapper
    └── bench_qulacs.cpp # Qulacs wrapper
```
