# Mixed State Benchmarks

Benchmarks comparing qlib's packed sparse density matrix implementation against BLAS dense, QuEST, and Quantum++.

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
--iterations N   Number of iterations (default: 1000)
--warmup N       Warm-up iterations (default: 100)
--csv            Output CSV format
--verbose        Show per-method memory usage inline
--help           Show help
```

### Examples

```bash
./benchmark/bench_mixed --max-qubits 8           # Quick test
./benchmark/bench_mixed --csv > results.csv      # Export data
./benchmark/bench_mixed --iterations 5000        # Higher precision
```

## Results (8 qubits, 256-dim)

| Method | X gate | H gate | Z gate |
|--------|--------|--------|--------|
| qlib packed | 215 ms | 196 ms | 108 ms |
| BLAS dense | 5848 ms | 5821 ms | 5852 ms |
| QuEST | 725 ms | 727 ms | 641 ms |
| Quantum++ | 1760 ms | 1694 ms | 1789 ms |

**Speedup:** qlib is 27-54x faster than BLAS, 3-6x faster than QuEST, 8-17x faster than Quantum++.

**Memory (8 qubits, runtime RSS):** qlib=0.5 MB, dense=9.3 MB, QuEST=2.3 MB, Qpp=1.0 MB

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

## Comparison Methods

| Method | Description | Complexity |
|--------|-------------|------------|
| qlib packed | Packed lower-triangular with stride traversal | O(N) per gate |
| BLAS dense | Full matrix UρU† via zgemm | O(N³) |
| naive loop | Element-by-element UρU† computation | O(N³) |
| QuEST | Multithreaded quantum simulator | O(N²) |
| Quantum++ | Eigen3-based C++17 library | O(N²) |

Note: naive loop is only run for ≤8 qubits due to performance.

Memory is measured as actual runtime RSS (resident set size) delta, not theoretical sizes. This captures real allocations including temporary buffers.

### Memory Measurement Methodology

Memory is measured as RSS (resident set size) delta capturing the full allocation cost.

**Process isolation**: QuEST and Quantum++ benchmarks run in forked child processes to ensure accurate memory measurement. This avoids memory reuse between benchmarks and captures the true cost of each framework's initialization.

**Full initialization**: Each framework's memory measurement includes all setup:
- **QuEST**: Measures `initQuESTEnv()` + `createDensityQureg()` + state initialization (forked process)
- **Quantum++**: Measures Eigen matrix allocation (forked process)
- **BLAS/naive**: Measures density matrix + gate matrix allocations
- **qlib**: Measures packed storage allocation

**Page-touching**: To avoid lazy allocation where the OS defers page mapping:
- **QuEST**: Uses `initDebugState()` to touch all state vector pages
- **Quantum++**: Initializes matrices with `Ones()` before setting values
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
- Library versions: QuEST, Quantum++, Eigen

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
- Quantum++ single-threaded by default
- BLAS may use multiple threads depending on implementation
- qlib is single-threaded
- For fair comparison, consider `OMP_NUM_THREADS=1`

### Limitations
- Only single-qubit gates benchmarked
- No multi-qubit gate operations (CNOT, CZ)
- No measurement operations
- Random state initialization (not physically meaningful states)
- Does not measure gate matrix construction time separately
- Memory measurement granularity limited by OS page size (typically 4KB/16KB)

## External Frameworks

QuEST and Quantum++ are included as git submodules in `extern/`. They are auto-detected during cmake configuration. Check the cmake output for status:

```
-- Benchmark configuration:
--   QuEST:     TRUE
--   Quantum++: TRUE
```

### Setup (if submodules not initialized)

```bash
git submodule update --init extern/QuEST extern/qpp
brew install eigen libomp  # macOS dependencies

# Build QuEST
cd extern/QuEST && mkdir build && cd build
OpenMP_ROOT=$(brew --prefix)/opt/libomp cmake ..
cmake --build .
```

### Custom paths

If frameworks are installed elsewhere, pass hints to cmake:

```bash
cmake -DQuEST_DIR=/path/to/quest -DQPP_DIR=/path/to/qpp ..
```

## Files

```
benchmark/
├── CMakeLists.txt
├── include/bench.h      # Timing utilities, result types
└── src/
    ├── bench_mixed.c    # Main benchmark + BLAS/naive implementations
    ├── bench_quest.c    # QuEST wrapper
    └── bench_qpp.cpp    # Quantum++ wrapper
```
