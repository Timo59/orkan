# Benchmark Module

Technical specification for `bench_gate`, which measures single-gate throughput of qlib's mixed-state implementations and compares them against external frameworks. Circuit-level benchmarking is out of scope; a future `bench_circuit` executable is reserved for that purpose.

## File Layout

```
benchmark/
├── CMakeLists.txt              # Release-only build; auto-detects QuEST + Qulacs
├── README.md                   # Developer convenience: setup instructions, sample results, methodology
├── include/
│   └── bench.h                 # Types, timing helpers, statistical helpers, memory utilities, shared constants, function declarations
└── src/
    ├── bench_mixed.c           # qlib adapters (packed/tiled, 1Q and 2Q) + shared build_all_pairs() utility
    ├── bench_baselines.c       # BLAS/naive baselines: gate matrices, bench_blas_dense/naive_loop/blas_dense_2q
    ├── bench_main.c            # Orchestration: gate tables, CLI parsing, output formatting, main()
    ├── bench_quest.c           # QuEST adapter (compiled only when WITH_QUEST is defined)
    └── bench_qulacs.cpp        # Qulacs adapter (compiled only when WITH_QULACS is defined)
```

`README.md` is a developer convenience file. This spec is the authoritative reference.

`bench.h` declares all public types, constants, inline utilities, and function prototypes — including `build_all_pairs()` and `BENCH_MAX_RUNS`. Output functions (`bench_print_result`, `bench_print_csv`, `bench_print_csv_header`, `bench_parse_options`) are implemented in `bench_main.c`. Baseline functions (`bench_blas_dense`, `bench_naive_loop`, `bench_blas_dense_2q`) are implemented in `bench_baselines.c`.

## Dependencies

**qlib module (`libq`):** `bench_gate` links against `libq` and exercises the public gate API directly — `src/gate/packed/` (packed lower-triangular) and `src/gate/tiled/` (tiled lower-triangular). BLAS and OpenMP linkage are inherited transitively via `libq`'s public CMake dependencies.

**External frameworks** (optional, auto-detected at configure time):

| Dependency | Location | Detection |
|------------|----------|-----------|
| QuEST v4 | `extern/QuEST/` | `find_path(quest/include/quest.h)` + `find_library(QuEST)` |
| Qulacs (`csim`) | `extern/qulacs/src/csim/` | Presence of `update_ops_dm.hpp` |
| Eigen | System or Homebrew | `find_path(Eigen/Core)` |
| OpenMP | System or Homebrew libomp | `find_package(OpenMP)` (macOS: auto-configured from `brew --prefix libomp`) |

If a framework is absent, its columns are omitted from all output; the binary still builds.

## Build Integration

### Requirements

Benchmarks require a **Release build**. `benchmark/CMakeLists.txt` calls `return()` immediately (with a status message) if `CMAKE_BUILD_TYPE != Release`.

### Commands

```bash
# From project root
cmake -B cmake-build-release -DCMAKE_BUILD_TYPE=Release -DBUILD_BENCHMARKS=ON
cmake --build cmake-build-release --target bench_gate
# Binary lands at: cmake-build-release/benchmark/bench_gate
```

### CMake targets produced

| Target | Type | Condition |
|--------|------|-----------|
| `bench_gate` | Executable | Always (Release only) |
| `csim_static` | Static library | `QULACS_FOUND` — compiled with `-w $<$<CONFIG:Release>:-march=native>` to match qlib's vectorisation flags |

### Compile definitions

`WITH_QUEST` is defined when QuEST is found. `WITH_QULACS` is defined when Qulacs is found. Both are set automatically by CMake — there are no manual `-DWITH_QUEST=ON` or `-DWITH_QULACS=ON` flags.

### Custom install paths

Pass `-DQuEST_DIR=/path` to cmake for a non-default QuEST location. There is no `QULACS_DIR` hint; Qulacs must reside at `extern/qulacs/`.

### Configuration summary

CMake prints at configure time:

```
-- Benchmark configuration:
--   QuEST:     TRUE/FALSE
--   Qulacs:    TRUE/FALSE
--   OpenMP:    TRUE/FALSE
```

The OpenMP line reflects `OpenMP_CXX_FOUND`.

## CLI Reference

```bash
./benchmark/bench_gate                                    # All defaults
./benchmark/bench_gate --min-qubits 4 --max-qubits 8
./benchmark/bench_gate --step 2
./benchmark/bench_gate --iterations 5000
./benchmark/bench_gate --runs 20
./benchmark/bench_gate --csv > results.csv
./benchmark/bench_gate --pgfplots > results.dat
./benchmark/bench_gate --verbose
./benchmark/bench_gate --help
```

All flags:

| Flag | Default | Valid range | Purpose |
|------|---------|-------------|---------|
| `--min-qubits N` | 2 | 1–30 | Minimum qubit count |
| `--max-qubits N` | 12 | min–30 | Maximum qubit count |
| `--step N` | 1 | ≥1 | Qubit count increment |
| `--iterations N` | 1000 | ≥1 | Gate calls per timed run |
| `--warmup N` | 100 | ≥0 | Untimed warm-up iterations before each run series |
| `--runs N` | 10 | 1–200 | Independent timing runs per (gate, qubits, method) |
| `--csv` | off | — | One CSV row per measurement; output to stdout |
| `--pgfplots` | off | — | Five normalized timing tables per gate + one memory table; output to stdout |
| `--verbose` | off | — | Append min, median, and memory bytes to each console result line |
| `--help` | — | — | Print usage and exit |

All values are validated at startup; invalid values print to stderr and exit with code 1. Defaults are defined as macros in `bench.h` (`BENCH_DEFAULT_*`). The maximum runs count (`BENCH_MAX_RUNS = 200`) is also defined in `bench.h` and used by both the validation check and the stats computation.

**`blas_dense` and `naive_loop`** are skipped for qubit counts above 8 (hardcoded: `run_dense = (qubits <= 8)`). This cutoff is not derived from `--max-qubits` and is not configurable.

**`naive_loop`** runs at 1/10th the iteration and warmup counts relative to `--iterations`/`--warmup` (hardcoded in `bench_main.c`). The reported speedup is scaled by ×10 to make it comparable to other methods.

## Statistical Model

### Structure of a single measurement

For each `(gate, qubits, method)` triple:

1. `warmup` untimed gate applications are performed, cycling through all qubit targets (`i % qubits`) or all qubit pairs (`i % num_pairs`), to bring code paths and data into cache.
2. `runs` independent timed runs are performed. Each run applies the gate `iterations × sweep_size` times and is measured as a single wall-clock interval.
3. `bench_compute_stats()` computes the five statistics below from the `runs` wall-clock measurements.

### `sweep_size`

`sweep_size` is the number of gate applications per timed iteration:

- **1Q gates**: `sweep_size = qubits` (one application per qubit, target 0..n-1)
- **2Q gates**: `sweep_size = qubits * (qubits - 1) / 2` (one application per ordered pair `(q1, q2)` with `q1 < q2`)

The 2Q all-pairs sweep is deliberate: it exercises within-tile, cross-tile, and tile-boundary code paths in a single run.

### Statistics collected

All statistics are computed by `bench_compute_stats()` from the array of `runs` wall-clock measurements (each in ms):

| Statistic | `bench_result_t` field | Definition |
|-----------|------------------------|------------|
| Mean | `time_ms` | Arithmetic mean of `runs` wall-clock measurements (ms) |
| Std dev | `time_ms_std` | Sample standard deviation, Bessel-corrected: `sqrt(Σ(xᵢ−x̄)² / (n−1))` (ms); 0 if `runs == 1` |
| Minimum | `time_ms_min` | Minimum of `runs` measurements — best-case hardware throughput (ms) |
| Median | `time_ms_median` | Median of `runs` measurements — robust to scheduler spikes (ms) |
| CV | `time_ms_cv` | Coefficient of variation: `std_dev / mean × 100` (%) |

`ops_per_sec` is derived from the mean: `(iterations × sweep_size) / (time_ms / 1000)`.

`time_per_gate_ms` (CSV only) is the mean normalized per gate application: `time_ms / (iterations × sweep_size)`.

### Reproducibility

State data is filled with random values via `rand()` at state initialisation. No fixed seed is set; the PRNG state is whatever it happens to be at the time of the call. Since state values do not affect gate execution timing (no data-dependent branches in the arithmetic), this has no impact on measurement fairness or repeatability.

## Output Formats

### Console (default)

For each qubit count, for each gate: one result line per method showing method name, `time_ms ± time_ms_std`, CV, and ops/sec. Format:

```
  %-14s %8.3f ± %6.3f ms (CV %4.1f%%)  %12.0f ops/sec
```

With `--verbose`, each result line is extended with:

```
  [min %8.3f ms  med %8.3f ms  mem %zu B]
```

After each gate block: a speedup summary line (packed vs other methods) and a memory comparison line (values auto-scaled to KB or MB, with packed vs dense savings percentage).

### CSV (`--csv`)

One row per `(gate, qubits, method)` triple, written to stdout. Header:

```
qubits,dim,gate,method,runs,iterations,sweep_size,time_ms_mean,time_ms_std,time_ms_min,time_ms_median,time_ms_cv,time_per_gate_ms,ops_per_sec,memory_bytes
```

Row format: `%u,%lld,%s,%s,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.3f,%.9f,%.0f,%zu`.

`sweep_size` and `time_per_gate_ms` are computed at emit time via the `SET_SWEEP` macro in `bench_main.c` and are not stored in `bench_result_t` between benchmarks.

### pgfplots (`--pgfplots`)

Output to stdout. For each of the 5 gates (X, H, Z, CX, SWAP), five tables are emitted in this order:

| Table | Column header comment | Unit | Definition |
|-------|-----------------------|------|------------|
| Mean per-gate time | `time_per_gate_us` | µs | `time_ms × 1000 / (iterations × sweep_size)` |
| 95% CI half-width | `ci95_per_gate` | µs | `t_crit(runs−1) × (std / sqrt(runs)) × 1000 / (iterations × sweep_size)` |
| Min per-gate time | `min_per_gate` | µs | `time_ms_min × 1000 / (iterations × sweep_size)` |
| Median per-gate time | `median_per_gate` | µs | `time_ms_median × 1000 / (iterations × sweep_size)` |
| CV | `cv` | % | `time_ms_cv` (dimensionless; invariant under normalization) |

All timing values are in **µs** (microseconds), not ms. CV is reproduced directly from the stored percentage value without renormalization.

Each table has a leading `qubits` column followed by one column per method. Five methods appear: `qlib_packed`, `qlib_tiled`, `blas`, `quest`, `qulacs`. `naive_loop` is neither executed nor stored when `--pgfplots` is active. Missing values (framework not built, or qubit count exceeds the dense cutoff) are emitted as `nan`.

After all per-gate tables, one memory table is emitted with values in MB:

```
# Memory usage (MB)
qubits  qlib_packed     qlib_tiled      blas            quest           qulacs
```

Memory values are taken from gate index 0 (X gate) for all methods. Missing values are `nan`.

#### Confidence interval formula

```
CI₉₅ = t_crit(df) × (σ / √n) × 1000 / (iterations × sweep_size)
```

where `df = runs − 1`, `σ = time_ms_std`, and `n = runs`. The t-critical value is looked up from a hard-coded table covering df = 1..49 exactly (49-entry array in `bench_main.c`). For df ≥ 50, the normal approximation 1.96 is used. If `runs < 2`, the CI is emitted as `nan`.

## Key Data Structures

All types are declared in `bench.h`. See that file for the authoritative field-level documentation.

**`bench_result_t`** — result for one `(gate, qubits, method)` triple. Carries timing statistics (`time_ms`, `time_ms_std`, `time_ms_min`, `time_ms_median`, `time_ms_cv`, `ops_per_sec`), memory footprint (`memory_bytes`), run metadata (`iterations`, `runs`), and identification (`qubits`, `dim`, `gate_name`, `method`). `sweep_size` and `time_per_gate_ms` are populated at emit time only (CSV path, via `SET_SWEEP`); they are zero-initialized and not set by the benchmark functions themselves.

**`bench_options_t`** — parsed CLI options.

**`bench_run_stats_t`** — internal helper returned by `bench_compute_stats()`. Not stored persistently; values are immediately copied into `bench_result_t`.

**`pgfplots_data_t`** (defined in `bench_main.c`) — heap-allocated only when `--pgfplots` is passed. Stores six 3-D arrays indexed as `[gate_index][qubit_config_index][method_index]` for all timing statistics and memory. `MAX_QUBIT_CONFIGS = 32`. Normalization to per-gate µs values is computed at emit time in `print_pgfplots_output()`.

**`bench_summary_t`** (defined in `bench_main.c`) — collects all method results for one `(gate, qubit-count)` combination and is passed to `print_speedup_and_memory()`.

## Implementations Benchmarked

Six method slots are defined (enum in `bench_main.c`):

| Index | Enum | Method string | Description | Source |
|-------|------|---------------|-------------|--------|
| 0 | `M_QLIB` | `qlib_packed` | Production packed lower-triangular | `src/gate/packed/` |
| 1 | `M_QLIB_TILED` | `qlib_tiled` | Production tiled lower-triangular | `src/gate/tiled/` |
| 2 | `M_BLAS` | `blas_dense` | Full N×N density matrix, `zgemm` for UρU† | `bench_baselines.c` |
| 3 | `M_NAIVE` | `naive_loop` | Element-by-element UρU† loop | `bench_baselines.c` |
| 4 | `M_QUEST` | `quest` | QuEST density-matrix backend | `bench_quest.c` |
| 5 | `M_QULACS` | `qulacs` | Qulacs csim density-matrix backend | `bench_qulacs.cpp` |

QuEST and Qulacs are conditionally compiled. When not present, their slots contain zero-initialized results and are omitted from output.

## Gates Benchmarked

**Single-qubit** (all methods): X, H, Z.

**Two-qubit** (no `naive_loop`):

| Gate | Tiled | Notes |
|------|-------|-------|
| CX | Yes | |
| SWAP | Yes | |

Each 2Q gate iterates over all `n*(n-1)/2` ordered pairs `(q1, q2)` with `q1 < q2`. This all-pairs sweep is intentional: it exercises within-tile, cross-tile, and tile-boundary code paths in a single run, providing a representative average over the gate's full qubit-pair distribution.

## Timing and Memory Methodology

### Timing

All methods — qlib, `blas_dense`, `naive_loop`, and both external framework adapters — use `clock_gettime(CLOCK_MONOTONIC)` via `bench_time_ns()` defined in `bench.h`.

### Memory reporting

All methods report only the **theoretical** state storage footprint, not process or framework overhead. Scratch buffers allocated outside the hot loop are not included.

| Method | Measurement |
|--------|-------------|
| `qlib_packed` | `dim*(dim+1)/2 × sizeof(cplx_t)` via `bench_packed_size()` |
| `qlib_tiled` | tile-pair count × `TILE_SIZE × sizeof(cplx_t)` via `bench_tiled_size()` |
| `blas_dense` | `dim² × sizeof(cplx_t)` via `bench_dense_size()` |
| `naive_loop` | same as `blas_dense` |
| QuEST | `dim² × sizeof(double) × 2` (full density matrix, complex double) |
| Qulacs | `dim² × sizeof(CTYPE)` |

Note: `blas_dense` and `naive_loop` pre-allocate all per-target gate matrices (one full N×N matrix per qubit or qubit-pair) outside the timed loop. At 12 qubits these matrices dominate memory; the reported `memory_bytes` covers only the density matrix, not these auxiliary buffers.

### External framework adapters

**QuEST** uses an init-once pattern. `bench_quest_init()` calls `initQuESTEnv()` once at program start. Each `bench_quest()` call creates a fresh `Qureg` with `createDensityQureg`, calls `initDebugState` to touch all pages (ensuring page faults occur before timing begins), runs warmup + timed rounds, then destroys the `Qureg`. `bench_quest_cleanup()` calls `finalizeQuESTEnv()` once at program exit. This avoids forking while respecting QuEST's init/finalize constraint.

**Qulacs** runs directly in-process. `bench_qulacs()` allocates the density matrix with `malloc` and initialises it with non-trivial values to fault all pages before the timed region, runs warmup + timed rounds, then frees with `free`. No global state; no init/finalize.

## Status

**In Progress.** Single-gate throughput benchmarking is functional for 1Q gates (X, H, Z) and 2Q gates (CX, SWAP). External framework adapters (QuEST, Qulacs) are complete.

Known gaps:
- Gate matrix construction time is not separated from gate application time.
- Measurement (non-unitary) operations are not benchmarked.
- No thread-count sweep; OpenMP uses all available threads by default (control via `OMP_NUM_THREADS`).
- Circuit-level benchmarking is reserved for a future `bench_circuit` executable.
