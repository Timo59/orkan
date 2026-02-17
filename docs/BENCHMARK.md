# Benchmark Suite

Reference for the `bench_mixed` benchmark, which compares qlib's mixed-state gate implementations against each other and against external frameworks.

## File Layout

```
benchmark/
├── CMakeLists.txt              # Build config; auto-detects QuEST + Qulacs
├── include/
│   └── bench.h                 # Types, timing helpers, memory measurement, function declarations
└── src/
    ├── bench_mixed.c           # Main benchmark (gate definitions, CLI, output, benchmark loops)
    ├── bench_quest.c           # QuEST adapter (compiled only when WITH_QUEST is set)
    └── bench_qulacs.cpp        # Qulacs adapter (compiled only when WITH_QULACS is set)
```

## Build

Benchmarks are **Release-only** (guarded in `CMakeLists.txt`). Debug builds skip them.

```bash
# From project root
cd cmake-build-release          # or create: mkdir cmake-build-release && cd cmake-build-release
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .                 # builds bench_mixed (and all other targets)
```

## Running

```bash
# From the Release build directory
./benchmark/bench_mixed                                           # Default: 2-12 qubits, 1000 iter, 100 warmup
./benchmark/bench_mixed --max-qubits 10 --iterations 100          # Quick test run
./benchmark/bench_mixed --csv                                     # CSV output
./benchmark/bench_mixed --pgfplots                                # pgfplots-compatible .dat tables
./benchmark/bench_mixed --verbose                                 # Extra detail
```

All CLI options (see `--help`):

| Flag | Default | Purpose |
|------|---------|---------|
| `--min-qubits N` | 2 | Minimum qubit count |
| `--max-qubits N` | 12 | Maximum qubit count |
| `--iterations N` | 1000 | Timed iterations per gate/method |
| `--warmup N` | 100 | Warm-up iterations (not timed) |
| `--csv` | off | One CSV row per measurement |
| `--pgfplots` | off | One `.dat` table per gate (for LaTeX plots) |
| `--verbose` | off | Show additional detail in console output |

Qubit counts step by 2 (2, 4, 6, ...). Dense BLAS and naive methods are skipped above 8 qubits (O(N^3) too slow).

## Methods Compared

Six methods are indexed in `bench_mixed.c` (enum `M_QLIB` ... `M_QULACS`):

| Index | Method | What it does | Source |
|-------|--------|-------------|--------|
| 0 | `qlib_packed` | Production packed lower-triangular (BLAS-1) | `src/gate/packed/` |
| 1 | `qlib_tiled` | Production tiled lower-triangular | `src/gate/tiled/` |
| 2 | `blas_dense` | Reference: full N x N matrix, `zgemm` U rho U† | `bench_mixed.c` |
| 3 | `naive_loop` | Reference: element-by-element loop | `bench_mixed.c` |
| 4 | `quest` | QuEST density-matrix backend | `bench_quest.c` |
| 5 | `qulacs` | Qulacs csim density-matrix backend | `bench_qulacs.cpp` |

QuEST and Qulacs are auto-detected from `extern/` by CMake and conditionally compiled (`WITH_QUEST`, `WITH_QULACS` defines).

## Gates Benchmarked

Defined in `bench_mixed.c` around line 678:

**Single-qubit** (`gates[]`): X, H, Z — always benchmarked with all six methods.

**Two-qubit** (`gates_2q[]`): CX, SWAP — each entry carries a `has_tiled` flag that controls whether the tiled method is included.

| Gate | `has_tiled` | Notes |
|------|-------------|-------|
| CX | 1 | Tiled implementation in `gate_tiled_cx.c` |
| SWAP | 1 | Tiled implementation in `gate_tiled_swap.c` |

Two-qubit gates do not run the naive method. The `has_tiled` flag is the single switch for enabling/disabling tiled benchmarking per 2Q gate. To add a new 2Q gate, append to `gates_2q[]` and set the flag.

## Key Data Structures

Defined in `bench.h`:

- **`bench_result_t`**: Single measurement (qubits, dim, gate name, method name, time_ms, ops/sec, memory bytes, iterations).
- **`bench_options_t`**: Parsed CLI options.
- **`pgfplots_data_t`** (in `bench_mixed.c`): Accumulates all results for `--pgfplots` output. Indexed as `time_ms[gate][qubit_config][method]`.

## External Dependencies

Both are optional. If absent, their columns are simply omitted.

| Dependency | Location | Build |
|------------|----------|-------|
| QuEST v4 | `extern/QuEST/` | Pre-built; CMake finds `quest/include/quest.h` + `libQuEST` |
| Qulacs | `extern/qulacs/` | Built in-tree as `csim_static` from `src/csim/*.cpp`; needs Eigen |

## Memory Layout Background

The benchmark exists to compare two production storage formats for density matrices:

**Packed** (`MIXED_PACKED`): LAPACK packed lower-triangular column-major. Index formula: `(i,j) -> j*N - j*(j-1)/2 + (i-j)` for `i >= j`. Favors BLAS-1 operations but suffers variable-stride access (CONJ_OP) for high target qubits.

**Tiled** (`MIXED_TILED`): Lower triangle partitioned into `TILE_DIM x TILE_DIM` tiles (default 32). Intra-tile layout is **row-major** (not column-major). Tile index: `(tr*(tr+1)/2 + tc) * TILE_SIZE + lr*TILE_DIM + lc`. Eliminates CONJ_OP bottleneck via contiguous cache-friendly access.

See `docs/STATE_MODULE.md` for the full memory layout specification.
