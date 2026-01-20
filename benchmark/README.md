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

## Comparison Methods

| Method | Description | Complexity |
|--------|-------------|------------|
| qlib packed | Packed lower-triangular with stride traversal | O(N) per gate |
| BLAS dense | Full matrix UρU† via zgemm | O(N³) |
| QuEST | Multithreaded quantum simulator | O(N²) |
| Quantum++ | Eigen3-based C++17 library | O(N²) |

## External Frameworks

QuEST and Quantum++ are included as git submodules in `extern/`. They are auto-detected during cmake configuration.

### Setup (if submodules not initialized)

```bash
git submodule update --init extern/QuEST extern/qpp
brew install eigen libomp  # macOS dependencies

# Build QuEST
cd extern/QuEST && mkdir build && cd build
OpenMP_ROOT=$(brew --prefix)/opt/libomp cmake ..
cmake --build .
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
