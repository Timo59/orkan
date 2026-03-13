# Benchmark Module

Technical specification for `bench_gate`, which measures single-gate throughput of qlib's mixed-state implementations and compares them against external frameworks in two modes: **sweep mode** (aggregate throughput over all qubit targets per qubit count) and **per-qubit mode** (per-target breakdown enabling position-by-position cost analysis). Per-qubit mode disaggregates what sweep mode averages — it measures the same gate on each individual target qubit (or pair) separately rather than aggregating across all targets. Circuit-level benchmarking is out of scope.

## File Layout

```
benchmark/
├── CMakeLists.txt              # Release-only build; auto-detects QuEST + Qulacs + Aer DM
├── README.md                   # Developer convenience: setup instructions, sample results, methodology
├── include/
│   └── bench.h                 # Types, timing helpers, statistical helpers, memory utilities, shared constants, function declarations
└── src/
    ├── bench_mixed.c           # qlib adapters (packed/tiled, 1Q and 2Q) + shared build_all_pairs() and build_all_ordered_pairs() utilities
    ├── bench_baselines.c       # BLAS/naive baselines: gate matrices, bench_blas_dense/naive_loop/blas_dense_2q
    ├── bench_main.c            # Orchestration: gate tables, CLI parsing, output formatting, main()
    ├── bench_quest.c           # QuEST adapter (compiled only when WITH_QUEST is defined)
    ├── bench_qulacs.cpp        # Qulacs adapter (compiled only when WITH_QULACS is defined)
    └── bench_aer_dm.cpp        # Qiskit-Aer DM adapter (compiled only when WITH_AER_DM is defined)
```

`README.md` is a developer convenience file. This spec is the authoritative reference.

`bench.h` declares all public types, constants, inline utilities, and function prototypes — including `build_all_pairs()`, `build_all_ordered_pairs()`, and `BENCH_MAX_RUNS`. Output functions (`bench_print_result`, `bench_print_csv`, `bench_print_csv_header`, `bench_parse_options`, `bench_print_perq_csv_header`, `bench_print_perq_csv`, `bench_print_perq_console`) are implemented in `bench_main.c`. Baseline functions (`bench_blas_dense`, `bench_naive_loop`, `bench_blas_dense_2q`) are implemented in `bench_baselines.c`.

## Dependencies

**qlib module (`libq`):** `bench_gate` links against `libq` and exercises the public gate API directly — `src/gate/packed/` (packed lower-triangular) and `src/gate/tiled/` (tiled lower-triangular). BLAS and OpenMP linkage are inherited transitively via `libq`'s public CMake dependencies.

**External frameworks** (optional, auto-detected at configure time):

| Dependency | Location | Detection |
|------------|----------|-----------|
| QuEST v4 | `extern/QuEST/` | `find_path(quest/include/quest.h)` + `find_library(QuEST)` |
| Qulacs (`csim`) | `extern/qulacs/src/csim/` | Presence of `update_ops_dm.hpp` |
| Qiskit-Aer DM kernel | `extern/aer-dm/include/` | Presence of `simulators/density_matrix/densitymatrix.hpp` |
| nlohmann/json | System or Homebrew | `find_path(nlohmann/json.hpp)` (macOS: auto-configured from `brew --prefix nlohmann-json`); required by `extern/aer-dm/include/framework/json.hpp` |
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
| `bench_gate` | Executable (`EXCLUDE_FROM_ALL`) | Always (Release only) |
| `csim_static` | Static library | `QULACS_FOUND` — compiled with `-w $<$<CONFIG:Release>:-march=native>` to match qlib's vectorisation flags |
| `aer_dm_static` | Interface library | `AER_DM_FOUND` — carries SYSTEM include paths for vendored Aer headers and nlohmann/json |
| `bench_aer_only` | Executable (`EXCLUDE_FROM_ALL`) | `AER_DM_FOUND` — Aer-only data-collection variant; see below |

`bench_gate` and `bench_aer_only` are both `EXCLUDE_FROM_ALL`; they must be built with an explicit `--target` argument.

**`bench_aer_only`** uses the same source list as `bench_gate` with `BENCH_AER_ONLY_MODE=1` defined. This compile definition gates out all non-Aer method call sites in `bench_main.c` while preserving the full CLI, output formatting, and statistical infrastructure.

### Compile definitions

`WITH_QUEST` is defined when QuEST is found. `WITH_QULACS` is defined when Qulacs is found. Both are set automatically by CMake — there are no manual `-DWITH_QUEST=ON` or `-DWITH_QULACS=ON` flags.

`WITH_AER_DM` is defined on `bench_gate` and `bench_aer_only` when headers are found at `extern/aer-dm/include/`. Set automatically; no manual flag required. When headers are present, the `qiskit_aer_dm` column appears in all output modes.

### Custom install paths

Pass `-DQuEST_DIR=/path` to cmake for a non-default QuEST location. There is no `QULACS_DIR` hint; Qulacs must reside at `extern/qulacs/`.

### Configuration summary

CMake prints at configure time:

```
-- Benchmark configuration:
--   QuEST:     TRUE/FALSE
--   Qulacs:    TRUE/FALSE
--   Aer DM:    TRUE/FALSE
--   OpenMP:    TRUE/FALSE
```

The OpenMP line reflects `OpenMP_CXX_FOUND`.

## Key Data Structures

All types are declared in `bench.h`. Field-level documentation is authoritative there; this section summarises what each type represents.

**`bench_result_t`** — result for one `(gate, qubits, method)` triple. Carries timing statistics (`time_ms`, `time_ms_std`, `time_ms_min`, `time_ms_median`, `time_ms_cv`, `ops_per_sec`), memory footprint (`memory_bytes`), run metadata (`iterations`, `runs`, `sweep_size`), and identification (`qubits`, `dim`, `gate_name`, `method`). The fields `sweep_size` and `time_per_gate_ms` are zero-initialized by adapter benchmark functions and are not populated during timing. They are set by the `SET_SWEEP(res, sw)` macro (defined in `bench_main.c`), which is invoked by the orchestration loop in `bench_main.c` for each result struct immediately before passing that struct to the output functions (`bench_print_csv()`, pgfplots emitters). `SET_SWEEP` writes `sweep_size` directly into the struct and derives `time_per_gate_ms` as `time_ms / (iterations × sweep_size)`.

**`bench_result_perq_t`** — result for one `(gate, qubits, target, method)` tuple. Extends the sweep result with per-gate microsecond timings (`time_per_gate_us`, `time_per_gate_us_median`, `time_per_gate_us_min`) and explicit target fields (`target`, `target2`, `is_2q`). `target2` is `QUBIT_NONE` (defined in `bench.h` as `(qubit_t)-1`, the maximum value of `qubit_t`) for 1Q gates. See `bench.h` for the full field list.

**`bench_options_t`** — parsed CLI options. Includes `runs_explicit` and `iterations_explicit` flags that track whether `--runs`/`--iterations` were supplied on the command line (used to apply per-qubit mode default overrides).

**`bench_run_stats_t`** — public typedef in `bench.h`. Returned by `bench_compute_stats()`. Contains five fields: `mean`, `std_dev`, `min`, `median`, and `cv` — all derived from an array of K wall-clock run times in ms. Not stored persistently; values are immediately copied into the result struct by the caller.

**`bench_harness_t`** — encapsulates an opaque gate callback, its context, and loop counts for the per-qubit timed-run harness (`bench_run_timed()`).

**`bench_run_timed()`** — public `static inline` function in `bench.h`. Standard timing harness used by all `_at` adapter implementations. Takes a `const bench_harness_t *h`, an output array `double *out` of at least `h->runs` elements, and `int qubits` (controls warmup batch size). Performs time-based warmup (until `BENCH_PQ_WARMUP_MS` has elapsed), then records `h->runs` independent timing samples in milliseconds.

**`bench_fill_perq_stats()`** — public `static inline` function in `bench.h`. Signature: `void bench_fill_perq_stats(bench_result_perq_t *r, const bench_run_stats_t *s, int iterations)`. Populates all timing fields of `*r` (`time_ms_mean`, `time_ms_std`, `time_ms_min`, `time_ms_median`, `time_ms_cv`, `ops_per_sec`, `time_per_gate_us`, `time_per_gate_us_median`, `time_per_gate_us_min`, `iterations`) from a completed `bench_run_stats_t` and the gate application count. Because `bench_run_stats_t` carries no run-count field, callers must set `r->runs` separately after calling this function. This is the canonical pattern that every `_at` adapter implementation must follow: call `bench_run_timed()` to collect samples, `bench_compute_stats()` to summarize them, `bench_fill_perq_stats()` to populate the result, then set `r->runs` explicitly.

## CLI Reference

### Sweep Mode (default)

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

### Per-Qubit Mode

```bash
./benchmark/bench_gate --per-qubit
./benchmark/bench_gate --per-qubit --gate X
./benchmark/bench_gate --per-qubit --gate CX --max-qubits 10
./benchmark/bench_gate --per-qubit --csv > perq.csv
./benchmark/bench_gate --per-qubit --runs 30 --iterations 500
```

### All Flags

| Flag | Default | Valid range | Purpose |
|------|---------|-------------|---------|
| `--min-qubits N` | 2 | 1–30 | Minimum qubit count |
| `--max-qubits N` | 12 | min–30 | Maximum qubit count |
| `--step N` | 1 | ≥1 | Qubit count increment |
| `--iterations N` | 1000 (sweep) / 100 (per-qubit) | ≥1 | Gate calls per timed run; overridden by per-qubit calibration |
| `--warmup N` | 100 | ≥0 | Untimed warm-up iterations before each run series (sweep mode only; per-qubit uses time-based warmup) |
| `--runs N` | 10 (sweep) / 15 (per-qubit) | 1–200 | Independent timing runs per measurement tuple |
| `--per-qubit` | off | — | Enable per-qubit (per-target) benchmark mode |
| `--gate NAME` | (all gates) | Gate name string | Restrict to a single named gate; only meaningful with `--per-qubit` |
| `--csv` | off | — | Emit CSV to stdout (format differs between modes; see Output Formats) |
| `--pgfplots` | off | — | Five normalized timing tables per gate + one memory table; stdout (sweep mode only) |
| `--verbose` | off | — | Append min, median, and memory bytes to each console result line (sweep mode only) |
| `--help` | — | — | Print usage and exit |

**Default overrides in per-qubit mode:** When `--per-qubit` is active and `--runs`/`--iterations` were not explicitly supplied, defaults are overridden to `BENCH_PQ_DEFAULT_RUNS = 15` and `BENCH_PQ_DEFAULT_ITERATIONS = 100`. Explicitly supplied values are honoured. These overrides are tracked via `runs_explicit` and `iterations_explicit` in `bench_options_t`.

**`--gate` without `--per-qubit`:** The `--gate` flag has no effect in sweep mode. If specified without `--per-qubit`, a warning is printed to stderr and the flag is silently ignored.

**Gate name matching:** `--gate NAME` is a case-insensitive match against the gate name strings in the gate table: `"X"`, `"H"`, `"Z"`, `"CX"`, `"SWAP"` (uses `strcasecmp`).

**Mutual exclusion:** `--per-qubit` and `--pgfplots` are mutually exclusive. Passing both prints an error to stderr and exits with code 1.

All values are validated at startup; invalid values print to stderr and exit with code 1. Sweep-mode defaults are defined as `BENCH_DEFAULT_*` macros in `bench.h`. `BENCH_MAX_RUNS = 200` is defined in `bench.h` and used by both the validation check and stats computation.

**`blas_dense` and `naive_loop`** are skipped for qubit counts above 8 (hardcoded: `run_dense = (qubits <= 8)`). This cutoff is not derived from `--max-qubits` and is not configurable.

**`naive_loop`** runs at 1/10th the iteration and warmup counts relative to `--iterations`/`--warmup` (hardcoded in `bench_main.c`; 1Q gates only — `naive_loop` is not run for 2Q gates). The reported speedup is scaled by ×10 to make it comparable to other methods (1Q gates only).

## Implementations Benchmarked

Seven method slots are defined (enum in `bench_main.c`):

| Index | Enum | Method string | Description | Source |
|-------|------|---------------|-------------|--------|
| 0 | `M_QLIB` | `qlib_packed` | Production packed lower-triangular | `src/gate/packed/` |
| 1 | `M_QLIB_TILED` | `qlib_tiled` | Production tiled lower-triangular | `src/gate/tiled/` |
| 2 | `M_BLAS` | `blas_dense` | Full N×N density matrix, `zgemm` for UρU† | `bench_baselines.c` |
| 3 | `M_NAIVE` | `naive_loop` | Element-by-element UρU† loop | `bench_baselines.c` |
| 4 | `M_QUEST` | `quest` | QuEST density-matrix backend | `bench_quest.c` |
| 5 | `M_QULACS` | `qulacs` | Qulacs csim density-matrix backend | `bench_qulacs.cpp` |
| 6 | `M_AER_DM` | `qiskit_aer_dm` | Qiskit-Aer internal C++ density matrix kernel | `bench_aer_dm.cpp` |

QuEST, Qulacs, and Aer DM are conditionally compiled. When not present, their slots contain zero-initialized results and are omitted from output.

**Adding a new method:** Implement a sweep function returning `bench_result_t` and per-qubit `_at` variants returning `bench_result_perq_t`, following the signatures in `bench.h`. Register the method in the enum and gate table in `bench_main.c`. `bench_aer_only` (with `BENCH_AER_ONLY_MODE=1`) gates out all non-Aer call sites in `bench_main.c` via compile-time checks; new methods must respect that pattern if they are not Aer-only.

## Gates Benchmarked

**Single-qubit** (all methods): X, H, Z.

**Two-qubit** (no `naive_loop`):

| Gate | Tiled | Notes |
|------|-------|-------|
| CX | Yes | |
| SWAP | Yes | |

## Statistical Model

### Sweep Mode: Structure of a Single Measurement

For each `(gate, qubits, method)` triple:

1. `warmup` untimed gate applications are performed, cycling through all qubit targets (`i % qubits`) or all qubit pairs (`i % num_pairs`), to bring code paths and data into cache.
2. `runs` independent timed runs are performed. Each run applies the gate `iterations × sweep_size` times and is measured as a single wall-clock interval.
3. `bench_compute_stats()` computes statistics from the `runs` wall-clock measurements.

### Per-Qubit Mode: Structure of a Single Measurement

For each `(gate, qubits, target, method)` tuple:

1. **Calibration (two-phase adaptive probe):** Applied identically to all methods (`qlib_packed`, `qlib_tiled`, `blas_dense`, `quest`, `qulacs`, `qiskit_aer_dm`), using target qubit 0 (or pair `(q1s[0], q2s[0])` for 2Q gates). The resulting iteration count applies to all target positions for that method.
   - **Phase 1 — single-shot probe:** The gate is invoked once (1 iteration, 1 run). If the returned `time_per_gate_us / 1000.0 >= MIN_MEASUREMENT_MS = 1.0 ms`, the gate is slow enough on its own; `probe_iters` is set to 1. Otherwise `probe_iters = BENCH_PQ_PROBE_ITERS = 100`.
   - **Phase 2 — multi-run probe:** `BENCH_PQ_PROBE_RUNS = 5` runs of `probe_iters` iterations each are timed; the minimum elapsed time across those runs is taken. `calibrate_iterations()` then scales the iteration count up until a single run would exceed `MIN_MEASUREMENT_MS`. If the probe elapsed time rounds to 0 ns, `BENCH_PQ_FALLBACK_ITERATIONS = 10000` is used. The final calibrated count is capped at `BENCH_PQ_MAX_CALIBRATED_ITERATIONS = 1,000,000`.
2. **Warmup:** The gate callback is invoked in batches until at least `BENCH_PQ_WARMUP_MS = 100 ms` of wall-clock time has elapsed (time-based, not fixed-count). Batch size is hardcoded in `bench.h` and scales with qubit count to bound latency at large state sizes:

   | Qubit count | Warmup batch size |
   |-------------|-------------------|
   | ≤ 4 | 128 calls per batch |
   | ≤ 7 | 16 calls per batch |
   | ≤ 10 | 4 calls per batch |
   | > 10 | 1 call per batch |

3. `runs` independent timed runs are recorded. Each run applies the gate `iterations` times to the same target.
4. `bench_compute_stats()` computes statistics from the `runs` measurements.

### `sweep_size` (Sweep Mode)

`sweep_size` is the number of gate applications per timed iteration:

- **1Q gates**: `sweep_size = qubits` (one application per qubit, target 0..n-1)
- **2Q gates**: `sweep_size = qubits * (qubits - 1) / 2` (one application per unordered pair `(q1, q2)` with `q1 < q2`)

The 2Q all-pairs sweep exercises within-tile, cross-tile, and tile-boundary code paths in a single run, providing a representative average over the gate's full qubit-pair distribution.

In per-qubit mode, each `_at` call targets a single qubit or qubit pair; `sweep_size` is always 1.

### 2Q Pair-Building in Per-Qubit Mode

In per-qubit mode the pair list for a 2Q gate is built differently depending on whether the gate is symmetric, controlled by the `is_symmetric` field in `gates_2q_perq[]` (defined in `bench_main.c`):

| Gate | `is_symmetric` | Function | Pair count | Rationale |
|------|---------------|----------|------------|-----------|
| CX | 0 | `build_all_ordered_pairs()` | `N*(N-1)` | Control and target are distinct roles; `(q0,q1)` and `(q1,q0)` are different gates |
| SWAP | 1 | `build_all_pairs()` | `N*(N-1)/2` | Order is irrelevant; `(q0,q1)` and `(q1,q0)` are identical up to labelling |

In **sweep mode**, `SET_SWEEP` always uses `sw2q = N*(N-1)/2` for all 2Q gates regardless of `is_symmetric` — sweep mode does not use `build_all_ordered_pairs`.

**Adding a new 2Q gate:** Set `is_symmetric = 0` for gates where control/target order matters (CX-like), or `is_symmetric = 1` for gates where it does not (SWAP-like). An incorrect setting will silently produce the wrong number of per-qubit result rows and a misleading per-pair coverage claim.

### Statistics Collected

`bench_compute_stats()` operates on an array of `runs` wall-clock measurements (ms):

| Statistic | Field in `bench_result_t` (sweep) | Field in `bench_result_perq_t` (per-qubit) | Definition |
|-----------|-----------------------------------|---------------------------------------------|------------|
| Mean | `time_ms` | `time_ms_mean` | Arithmetic mean (ms) |
| Std dev | `time_ms_std` | `time_ms_std` | Sample standard deviation, Bessel-corrected: `sqrt(Σ(xᵢ−x̄)² / (n−1))`; 0 if `runs == 1` |
| Minimum | `time_ms_min` | `time_ms_min` | Minimum — best-case hardware throughput (ms) |
| Median | `time_ms_median` | `time_ms_median` | Median — robust to scheduler spikes (ms) |
| CV | `time_ms_cv` | `time_ms_cv` | `std_dev / mean × 100` (%) |

`ops_per_sec` is derived from the mean. In sweep mode: `(iterations × sweep_size) / (time_ms / 1000)`. In per-qubit mode: `iterations / (time_ms_mean / 1000)`.

Per-qubit mode additionally reports:
- `time_per_gate_us`: `time_ms_mean / iterations × 1000` (µs)
- `time_per_gate_us_median`: `time_ms_median / iterations × 1000` (µs)
- `time_per_gate_us_min`: `time_ms_min / iterations × 1000` (µs)

### Per-Qubit Confidence Interval Note

At the default `BENCH_PQ_DEFAULT_RUNS = 15` runs, the 95% CI on the Bessel-corrected standard deviation spans approximately [0.74σ, 1.83σ] (chi-squared distribution, 14 degrees of freedom). For publication-quality intervals, use at least 30–50 runs. Results with CV > 50% are flagged as `(noisy)` in console output; no IQR-based filtering is applied.

### Reproducibility

State data is filled with random values via `bench_fill_random()` (xoshiro256** seeded from OS entropy) at state initialisation. No fixed seed is set. Since gate arithmetic has no data-dependent branches, state values do not affect timing.

## Timing and Memory Methodology

### Timing

All methods use `clock_gettime(CLOCK_MONOTONIC)` via `bench_time_ns()` defined in `bench.h`. In per-qubit mode, `bench_timing_barrier()` is called three times per timed run within `bench_run_timed()`: immediately before the start-clock read, immediately after the start-clock read, and immediately after the timed loop (before reading elapsed time). This pattern ensures the clock reads themselves are separated from the measured region by serializing instructions. The barrier is architecture-aware (serializing instruction on ARM64 and x86-64; compiler-only barrier on other platforms) and is defined as a static inline in `bench.h`.

**Linux note:** For reliable results, set the CPU governor to `performance` mode (`cpupower frequency-set -g performance`). Variable-frequency governors can introduce systematic bias when warmup duration does not fully stabilize the clock.

### Memory Reporting

All methods report only the **theoretical** state storage footprint, not process or framework overhead. Scratch buffers allocated outside the hot loop are not included.

| Method | Measurement |
|--------|-------------|
| `qlib_packed` | `dim*(dim+1)/2 × sizeof(cplx_t)` via `bench_packed_size()` |
| `qlib_tiled` | tile-pair count × `TILE_SIZE × sizeof(cplx_t)` via `bench_tiled_size()` |
| `blas_dense` | `dim² × sizeof(cplx_t)` via `bench_dense_size()` |
| `naive_loop` | same as `blas_dense` |
| QuEST | `dim² × sizeof(double) × 2` (full density matrix, complex double) |
| Qulacs | `dim² × sizeof(CTYPE)` |
| `qiskit_aer_dm` | `dim² × sizeof(std::complex<double>)` |

Note: `blas_dense` and `naive_loop` pre-allocate all per-target gate matrices outside the timed loop. At 12 qubits these matrices dominate memory; the reported `memory_bytes` covers only the density matrix, not these auxiliary buffers.

### Huge Pages

State allocations at or above `BENCH_HUGEPAGE_THRESHOLD = 1 GB` use `bench_alloc_huge()` / `bench_free_huge()`, which request huge-page-backed memory via `mmap` on Linux (with `MADV_HUGEPAGE` for transparent huge pages) and macOS. Below the threshold, `bench_alloc_huge()` calls `malloc`. `bench_free_huge()` must be called with the same `n` as `bench_alloc_huge()`.

### External Framework Adapters

**QuEST** uses an init-once pattern: `bench_quest_init()` calls `initQuESTEnv()` once at program start; `bench_quest_cleanup()` calls `finalizeQuESTEnv()` once at exit. This is required because QuEST's global environment (`initQuESTEnv`/`finalizeQuESTEnv`) must be called exactly once per process. Each per-call or per-`_at` invocation creates a fresh `Qureg`, touches all pages before timing begins, runs warmup + timed rounds, then destroys the `Qureg`.

**Qulacs** and **Qiskit-Aer DM** run directly in-process with no global state or init/finalize. Each call allocates (or constructs) the state object, fills its internal buffer to fault all pages before timing begins, runs warmup + timed rounds, then frees (or destroys) the object. The Aer DM adapter provides explicit `bench_aer_dm_alloc()` / `bench_aer_dm_reinit()` / `bench_aer_dm_free()` helpers because `AER::QV::DensityMatrix<double>` is a C++ object that cannot be manipulated directly from C callers. Exceptions from C++ are caught at the `extern "C"` boundary; a zero `time_ms` result indicates failure (e.g., `std::bad_alloc` at large qubit counts).

## Adapter Function Signatures

Per-qubit measurement uses `_at` function variants. Each takes a pre-allocated state object, applies the gate to the specified target, and returns a `bench_result_perq_t`. The caller owns state allocation and lifetime; no allocation occurs inside the `_at` function.

**qlib backends (`bench_mixed.c`):**
```c
bench_result_perq_t bench_qlib_packed_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t),
    qubit_t target, int iterations, int runs, state_t *state);

bench_result_perq_t bench_qlib_tiled_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t),
    qubit_t target, int iterations, int runs, state_t *state);

bench_result_perq_t bench_qlib_packed_2q_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t, qubit_t),
    qubit_t q1, qubit_t q2, int iterations, int runs, state_t *state);

bench_result_perq_t bench_qlib_tiled_2q_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t, qubit_t),
    qubit_t q1, qubit_t q2, int iterations, int runs, state_t *state);
```

**BLAS baseline (`bench_baselines.c`):**
```c
bench_result_perq_t bench_blas_dense_at(qubit_t qubits, const char *gate_name,
    const cplx_t gate_mat[4],
    qubit_t target, int iterations, int runs, cplx_t *rho);

bench_result_perq_t bench_blas_dense_2q_at(qubit_t qubits, const char *gate_name,
    const cplx_t gate_mat[16],
    qubit_t q1, qubit_t q2, int iterations, int runs, cplx_t *rho);
```

**QuEST (`bench_quest.c` — exposed only when `WITH_QUEST` is defined and `quest/include/quest.h` has been included before `bench.h` in the same C translation unit; C++ inclusion is not supported due to QuEST v4 header incompatibilities that cause overload-resolution errors when included transitively):**
```c
bench_result_perq_t bench_quest_at(qubit_t qubits, const char *gate_name,
    qubit_t target, int iterations, int runs, Qureg *qureg);

bench_result_perq_t bench_quest_2q_at(qubit_t qubits, const char *gate_name,
    qubit_t q1, qubit_t q2, int iterations, int runs, Qureg *qureg);
```

**Qulacs (`bench_qulacs.cpp` — `WITH_QULACS` only):**
```c
bench_result_perq_t bench_qulacs_at(qubit_t qubits, const char *gate_name,
    qubit_t target, int iterations, int runs, void *rho);

bench_result_perq_t bench_qulacs_2q_at(qubit_t qubits, const char *gate_name,
    qubit_t q1, qubit_t q2, int iterations, int runs, void *rho);
```

**Qiskit-Aer DM (`bench_aer_dm.cpp` — `WITH_AER_DM` only):**
```c
bench_result_perq_t bench_aer_dm_at(qubit_t qubits, const char *gate_name,
    qubit_t target, int iterations, int runs, void *dm);

bench_result_perq_t bench_aer_dm_2q_at(qubit_t qubits, const char *gate_name,
    qubit_t q1, qubit_t q2, int iterations, int runs, void *dm);

void *bench_aer_dm_alloc(qubit_t qubits);
void  bench_aer_dm_reinit(void *dm, qubit_t qubits);
void  bench_aer_dm_free(void *dm);
```

## Constants

All constants are defined in `bench.h`.

| Constant | Value | Applies to | Meaning |
|----------|-------|------------|---------|
| `BENCH_DEFAULT_MIN_QUBITS` | 2 | Sweep | Default `--min-qubits` |
| `BENCH_DEFAULT_MAX_QUBITS` | 12 | Sweep | Default `--max-qubits` |
| `BENCH_DEFAULT_STEP` | 1 | Sweep | Default `--step` |
| `BENCH_DEFAULT_ITERATIONS` | 1000 | Sweep | Default `--iterations` in sweep mode |
| `BENCH_DEFAULT_WARMUP` | 100 | Sweep | Default `--warmup` in sweep mode |
| `BENCH_DEFAULT_RUNS` | 10 | Sweep | Default `--runs` in sweep mode |
| `BENCH_MAX_RUNS` | 200 | Both | Hard upper bound on `--runs`; also bounds the stats sort buffer |
| `BENCH_PQ_DEFAULT_RUNS` | 15 | Per-qubit | Default `--runs` override in per-qubit mode |
| `BENCH_PQ_DEFAULT_ITERATIONS` | 100 | Per-qubit | Starting iteration count before calibration |
| `BENCH_PQ_PROBE_ITERS` | 100 | Per-qubit | Iterations per probe run during calibration |
| `BENCH_PQ_PROBE_RUNS` | 5 | Per-qubit | Probe runs during calibration; minimum elapsed time is used |
| `BENCH_PQ_WARMUP_MS` | 100.0 | Per-qubit | Minimum warmup duration (ms) before timed runs begin |
| `MIN_MEASUREMENT_MS` | 1.0 | Per-qubit | Target minimum duration for a single timed run; calibration scales iterations until met |
| `BENCH_PQ_FALLBACK_ITERATIONS` | 10000 | Per-qubit | Fallback iteration count when probe elapsed time rounds to 0 ns |
| `BENCH_PQ_MAX_CALIBRATED_ITERATIONS` | 1,000,000 | Per-qubit | Upper bound on calibrated iteration count |
| `QUBIT_NONE` | `(qubit_t)-1` | Per-qubit | Sentinel for absent second qubit in 1Q gate results |
| `MAX_PAIRS` | 128 | Both | Maximum unordered pairs (q1 < q2) for `build_all_pairs()`; sufficient for ≤16 qubits |
| `MAX_TARGETS` | 256 | Per-qubit | `2 × MAX_PAIRS`; upper bound on per-qubit result array size (accommodates ordered pairs from `build_all_ordered_pairs()`) |
| `MAX_QUBITS_FOR_ORDERED_PAIRS` | 16 | Per-qubit | Maximum qubit count accepted by `build_all_ordered_pairs()`; returns 0 if exceeded |
| `BENCH_HUGEPAGE_THRESHOLD` | `1 << 30` (1 GB) | Both | State allocation size above which `bench_alloc_huge()` requests huge-page-backed memory |

## Output Formats

### Console (Sweep Mode, default)

For each qubit count, for each gate: one result line per method showing method name, `time_ms ± time_ms_std`, CV, and ops/sec. Columns: method name (14 chars), mean time (ms), std dev (ms), CV (%), ops/sec.

With `--verbose`, each result line is extended with: min time (ms), median time (ms), memory footprint (bytes).

After each gate block: a speedup summary line (packed vs other methods) and a memory comparison line (values auto-scaled to KB or MB, with packed vs dense savings percentage).

### Console (Per-Qubit Mode)

`bench_print_perq_console()` in `bench_main.c` is called once per method per `(gate, qubits)` block, receiving all targets for that method. Output differs by gate arity:

**1Q gates:** A header line `--- <gate> | <method> | <N> qubits ---` followed by a column header and one row per target qubit. Columns: `target` (qubit index), `mean(us/g)` (`time_per_gate_us`), `median(us/g)` (`time_per_gate_us_median`), `min(us/g)` (`time_per_gate_us_min`), `CV(%)`, `ops/sec`. Results with CV > 50% are suffixed `(noisy)`. Failed measurements (`time_ms_mean <= 0`) are silently skipped.

**2Q gates:** A header line followed by a single aggregated summary line: `min=X.XXX us (qA,qB)  mean=X.XXX us  max=X.XXX us (qC,qD)`, identifying the pair coordinates of the minimum and maximum. A note `(full per-pair data in CSV)` is appended. No per-pair rows are printed to console.

### CSV (Sweep Mode: `--csv`)

One row per `(gate, qubits, method)` triple, written to stdout.

Header columns:
`qubits`, `dim`, `gate`, `method`, `runs`, `iterations`, `sweep_size`, `time_ms_mean`, `time_ms_std`, `time_ms_min`, `time_ms_median`, `time_ms_cv`, `time_per_gate_ms`, `ops_per_sec`, `memory_bytes`

(`sweep_size` and `time_per_gate_ms` are set by the `SET_SWEEP` macro before this function is called; see Key Data Structures.)

### CSV (Per-Qubit Mode: `--per-qubit --csv`)

One row per `(gate, qubits, target, method)` tuple, written to stdout. Uses a distinct header from sweep CSV.

Header columns:
`qubits`, `dim`, `gate`, `method`, `runs`, `iterations`, `target`, `target2`, `is_2q`, `time_ms_mean`, `time_ms_std`, `time_ms_min`, `time_ms_median`, `time_ms_cv`, `time_per_gate_us`, `time_per_gate_us_median`, `time_per_gate_us_min`, `ops_per_sec`, `memory_bytes`

Notes:
- `time_per_gate_us` and related fields are in **microseconds**. Microseconds are used because per-gate times at small qubit counts are sub-millisecond and would lose significant digits in ms representation.
- `target2` is emitted as `-1` for 1Q gates (when `target2 == QUBIT_NONE`).
- Rows where `time_ms_mean <= 0.0` are silently skipped (indicates a failed measurement, e.g., `std::bad_alloc` in Aer DM at large qubit counts).

**Bandwidth not reported.** No bytes-per-second column is included. The bandwidth formula differs by a different factor for every backend (packed lower-triangular, tiled, dense full-matrix, and external framework internals each imply different effective access fractions), making a bandwidth column misleading for cross-backend comparison.

### pgfplots (`--pgfplots`, Sweep Mode Only)

Output to stdout. For each of the 5 gates (X, H, Z, CX, SWAP), five tables are emitted in this order:

| Table | Column name | Unit | Definition |
|-------|-------------|------|------------|
| Mean per-gate time | `time_per_gate_us` | µs | `time_ms × 1000 / (iterations × sweep_size)` |
| 95% CI half-width | `ci95_per_gate` | µs | `t_crit(runs−1) × (std / sqrt(runs)) × 1000 / (iterations × sweep_size)` |
| Min per-gate time | `min_per_gate` | µs | `time_ms_min × 1000 / (iterations × sweep_size)` |
| Median per-gate time | `median_per_gate` | µs | `time_ms_median × 1000 / (iterations × sweep_size)` |
| CV | `cv` | % | `time_ms_cv` (dimensionless; invariant under normalization) |

All timing values are in **µs**. CV is reproduced directly from the stored percentage value.

Each table has a leading `qubits` column followed by one column per method. Six methods appear: `qlib_packed`, `qlib_tiled`, `blas`, `quest`, `qulacs`, `qiskit_aer_dm`. `naive_loop` is not called when `--pgfplots` is active: the `bench_naive_loop()` call site is guarded by `if (!opts.pgfplots_output)` (inside the `if (run_dense)` block), so the function is never invoked and no result is produced. Separately, `M_NAIVE` is excluded from the pgfplots method index list regardless. Missing values (framework not built, or qubit count exceeds the dense cutoff) are emitted as `nan`.

After all per-gate tables, one memory table is emitted (values in MB). Memory values are taken from gate index 0 (X gate) for all methods.

#### Confidence Interval Formula

```
CI₉₅ = t_crit(df) × (σ / √n) × 1000 / (iterations × sweep_size)
```

where `df = runs − 1`, `σ = time_ms_std`, and `n = runs`. The t-critical value is looked up from a hard-coded table covering df = 1..49 (49-entry array in `bench_main.c`). For df ≥ 50, the normal approximation 1.96 is used. If `runs < 2`, CI is emitted as `nan`.

## Limitations and Known Gaps

| ID | Limitation |
|----|------------|
| L1 | **State initialization differs between backends.** qlib states are filled with `bench_fill_random()` (xoshiro256** from OS entropy); QuEST uses its own debug-state initializer; Qulacs and Aer DM use direct buffer fills. Initial state values do not affect timing (no data-dependent branches in gate arithmetic), but initialization cost is excluded from all timed regions. |
| L2 | **Per-qubit mode measures steady-state single-target throughput, not circuit-representative throughput.** Each benchmark applies the same gate to the same target in a tight loop — this saturates cache and hides realistic qubit-to-qubit state transition overhead. Results are not directly comparable to sweep-mode `time_per_gate_ms`, which averages over all targets in a single timed interval. |
| L3 | **NUMA and thread affinity not controlled.** On multi-socket systems, memory allocation and thread scheduling may land on different NUMA nodes across runs, introducing variance not captured by the CV metric. |
| L4 | **No IQR outlier filtering.** Only a CV > 50% flag is emitted. Scheduler spikes that fall within the distribution but inflate the mean are not removed. The minimum and median columns are more robust in the presence of such spikes. |
| L5 | **Variable iteration counts across methods (per-qubit mode).** Calibration targets `MIN_MEASUREMENT_MS = 1 ms` per run but runs independently per method, so different methods may execute different iteration counts. Use `time_per_gate_us` columns for cross-method comparison; do not compare raw `time_ms_mean` values across methods. |
| L6 | **Calibration-on-target-0 bias (per-qubit mode).** Iteration count calibration uses target qubit 0. For backends where per-qubit cost varies with target index (e.g., tile-boundary effects in the tiled backend), the calibrated count may result in timed runs that are too short for high-index targets. |
| L7 | **Gate matrix construction time is not separated from gate application time.** |
| L8 | **Measurement (non-unitary) operations are not benchmarked.** |
| L9 | **No thread-count sweep.** OpenMP uses all available threads by default; control via `OMP_NUM_THREADS`. |
| L10 | **Circuit-level benchmarking is not implemented.** |
