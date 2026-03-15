# Benchmark Module

Technical specification for the benchmark module, which contains two independent executables:

- **`bench_gate`** — measures single-gate throughput of qlib's mixed-state implementations in sweep and per-qubit modes, compared against external frameworks.
- **`bench_circuit`** — measures end-to-end circuit execution throughput across all backends on three application-motivated circuits (QAOA p=1, VQE-HEA, Quantum Volume).

Both targets share `bench.h` infrastructure (types, timing helpers, statistics) and the same external framework dependencies.

## File Layout

```
benchmark/
├── CMakeLists.txt              # Release-only build; auto-detects QuEST + Qulacs + Aer DM
├── README.md                   # Developer convenience: setup instructions, sample results, methodology
├── CIRCUIT_PLAN.md             # Implementation plan for bench_circuit (design decisions, rationale)
├── include/
│   ├── bench.h                 # Types, timing helpers, statistical helpers, memory utilities, shared constants, function declarations
│   └── bench_circuit.h         # Circuit abstraction (circuit_gate_t, circuit_op_t, circuit_t), runner typedef, context structs, circuit constructors, runner declarations
└── src/
    ├── bench_mixed.c           # qlib adapters (packed/tiled, 1Q and 2Q) + build_all_pairs() + build_all_ordered_pairs()
    ├── bench_baselines.c       # BLAS/naive baseline: gate matrices, bench_blas_dense/naive_loop/blas_dense_2q
    ├── bench_main.c            # bench_gate orchestration: gate tables, CLI parsing, output formatting, main()
    ├── bench_quest.c           # bench_gate QuEST adapter: sweep + per-qubit (WITH_QUEST only); owns all Qureg lifecycle code so bench_main.c never needs quest.h
    ├── bench_qulacs.cpp        # bench_gate Qulacs adapter (WITH_QULACS only)
    ├── bench_aer_dm.cpp        # bench_gate Qiskit-Aer DM adapter (WITH_AER_DM only)
    ├── bench_circuit.c         # Circuit lifecycle, append helpers, QAOA/VQE-HEA/QV constructors
    ├── bench_circuit_main.c    # bench_circuit orchestration: CLI parsing, calibration, output, main()
    ├── bench_circuit_qlib.c    # qlib packed/tiled circuit runners + context alloc/free
    ├── bench_circuit_quest.c   # QuEST circuit runner + context alloc/free (WITH_QUEST only)
    ├── bench_circuit_qulacs.cpp# Qulacs circuit runner + context alloc/free (WITH_QULACS only)
    ├── bench_circuit_aer_dm.cpp# Aer-DM circuit runner + context alloc/free (WITH_AER_DM only)
    └── bench_util.c            # Shared utilities: timing barrier, RNG, statistics, timed-run harness, huge-page allocator
```

`README.md` and `CIRCUIT_PLAN.md` are developer convenience files. This spec is the authoritative reference.

## Dependencies

**qlib module (`libq`):** Both targets link against `libq` and exercise the public gate and state API. BLAS and OpenMP linkage are inherited transitively via `libq`'s public CMake dependencies.

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
cmake --build cmake-build-release --target bench_circuit
# Binaries land at: cmake-build-release/benchmark/bench_gate
#                   cmake-build-release/benchmark/bench_circuit
```

### CMake targets produced

| Target | Type | Condition |
|--------|------|-----------|
| `bench_gate` | Executable (`EXCLUDE_FROM_ALL`) | Always (Release only) |
| `bench_circuit` | Executable (`EXCLUDE_FROM_ALL`) | Always (Release only) |
| `csim_static` | Static library | `QULACS_FOUND` — compiled with `-w $<$<CONFIG:Release>:-march=native>` |
| `aer_dm_static` | Interface library | `AER_DM_FOUND` — carries SYSTEM include paths for vendored Aer headers and nlohmann/json |
| `bench_aer_only` | Executable (`EXCLUDE_FROM_ALL`) | `AER_DM_FOUND` — Aer-only data-collection variant; see below |

All three executables are `EXCLUDE_FROM_ALL`; they must be built with an explicit `--target` argument.

**`bench_aer_only`** uses the same source list as `bench_gate` with `BENCH_AER_ONLY_MODE=1` defined. This compile definition gates out all non-Aer method call sites in `bench_main.c` while preserving the full CLI, output formatting, and statistical infrastructure. It is a temporary data-collection tool.

### Compile definitions

`WITH_QUEST` is defined on both `bench_gate` and `bench_circuit` when QuEST is found. `WITH_QULACS` is defined when Qulacs is found. `WITH_AER_DM` is defined when headers are found at `extern/aer-dm/include/`. All are set automatically by CMake — there are no manual `-DWITH_*=ON` flags.

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
--   bench_circuit: enabled (circuit-level benchmarks)
```

The OpenMP line reflects `OpenMP_CXX_FOUND`.

---

## bench_gate — Single-Gate Throughput Benchmark

`bench_gate` measures single-gate throughput in two modes: **sweep mode** (aggregate throughput over all qubit targets per qubit count) and **per-qubit mode** (per-target breakdown enabling position-by-position cost analysis). Circuit-level benchmarking is out of scope for `bench_gate`; see `bench_circuit` below.

### Key Data Structures (bench_gate)

All types are declared in `bench.h`. Field-level documentation is authoritative there; this section summarises what each type represents.

**`bench_result_t`** — result for one `(gate, qubits, method)` triple. Carries timing statistics (`time_ms`, `time_ms_std`, `time_ms_min`, `time_ms_median`, `time_ms_cv`, `ops_per_sec`), memory footprint (`memory_bytes`), run metadata (`iterations`, `runs`, `sweep_size`), and identification (`qubits`, `dim`, `gate_name`, `method`). The fields `sweep_size` and `time_per_gate_ms` are zero-initialized by adapter benchmark functions and are not populated during timing. They are set by the `SET_SWEEP(res, sw)` macro (defined in `bench_main.c`), which is invoked by the orchestration loop in `bench_main.c` for each result struct immediately before passing that struct to the output functions (`bench_print_csv()`, pgfplots emitters). `SET_SWEEP` writes `sweep_size` directly into the struct and derives `time_per_gate_ms` as `time_ms / (iterations × sweep_size)`.

In `bench_circuit`, `sweep_size` is repurposed to hold `n_ops` (gate count of the circuit), and `time_per_gate_ms` is computed as `time_ms / (iterations × n_ops)`. No `SET_SWEEP` macro is used; circuit runners set these fields directly.

**`bench_result_perq_t`** — result for one `(gate, qubits, target, method)` tuple. Extends the sweep result with per-gate microsecond timings (`time_per_gate_us`, `time_per_gate_us_median`, `time_per_gate_us_min`) and explicit target fields (`target`, `target2`, `is_2q`). `target2` is `QUBIT_NONE` (defined in `bench.h` as `(qubit_t)-1`, the maximum value of `qubit_t`) for 1Q gates. See `bench.h` for the full field list.

**`bench_options_t`** — parsed CLI options. Includes `runs_explicit` and `iterations_explicit` flags that track whether `--runs`/`--iterations` were supplied on the command line (used to apply per-qubit mode default overrides).

**`bench_run_stats_t`** — public typedef in `bench.h`. Returned by `bench_compute_stats()`. Contains five fields: `mean`, `std_dev`, `min`, `median`, and `cv` — all derived from an array of K wall-clock run times in ms. Not stored persistently; values are immediately copied into the result struct by the caller.

**`bench_harness_t`** — encapsulates an opaque gate callback, its context, and loop counts for the per-qubit timed-run harness (`bench_run_timed()`).

**`bench_run_timed()`** — public `static inline` function in `bench.h`. Standard timing harness used by all `_at` adapter implementations. Takes a `const bench_harness_t *h`, an output array `double *out` of at least `h->runs` elements, and `int qubits` (controls warmup batch size). Performs time-based warmup (until `BENCH_PQ_WARMUP_MS` has elapsed), then records `h->runs` independent timing samples in milliseconds.

**`bench_fill_perq_stats()`** — public function in `bench_util.c`. Signature: `void bench_fill_perq_stats(bench_result_perq_t *r, const bench_run_stats_t *s, int iterations, int runs)`. Populates all timing fields of `*r` from a completed `bench_run_stats_t`, the gate application count, and the run count. This is the canonical pattern that every `_at` adapter implementation must follow: call `bench_run_timed()` to collect samples, `bench_compute_stats()` to summarize them, then `bench_fill_perq_stats()` to populate the result.

### CLI Reference (bench_gate)

#### Sweep Mode (default)

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

#### Per-Qubit Mode

```bash
./benchmark/bench_gate --per-qubit
./benchmark/bench_gate --per-qubit --gate X
./benchmark/bench_gate --per-qubit --gate CX --max-qubits 10
./benchmark/bench_gate --per-qubit --csv > perq.csv
./benchmark/bench_gate --per-qubit --runs 30 --iterations 500
```

#### All Flags

| Flag | Default | Valid range | Purpose |
|------|---------|-------------|---------|
| `--min-qubits N` | 2 | 1–30 | Minimum qubit count |
| `--max-qubits N` | 12 | min–30 | Maximum qubit count |
| `--step N` | 1 | ≥1 | Qubit count increment |
| `--iterations N` | 1000 (sweep) / 100 (per-qubit) | ≥1 | Gate calls per timed run; overridden by per-qubit calibration |
| `--warmup N` | 100 | ≥0 | Untimed warm-up iterations before each run series (sweep mode only; per-qubit uses time-based warmup) |
| `--runs N` | 10 (sweep) / 15 (per-qubit) | 1–200 | Independent timing runs per measurement tuple |
| `--per-qubit` | off | — | Single-target throughput mode: benchmark each qubit target independently |
| `--gate NAME` | all | gate name | Restrict `--per-qubit` run to one named gate (e.g. `X`, `CX`); no effect without `--per-qubit` |
| `--csv` | off | — | Emit CSV to stdout (format differs between modes; see Output Formats) |
| `--pgfplots` | off | — | Five normalized timing tables per gate + one memory table; stdout. Combined with `--per-qubit`: emits per-target .dat tables aggregated over targets. |
| `--verbose` | off | — | Append min, median, and memory bytes to each console result line (sweep mode only) |
| `--help` | — | — | Print usage and exit |

**Default overrides in per-qubit mode:** When `--per-qubit` is active and `--runs`/`--iterations` were not explicitly supplied, defaults are overridden to `BENCH_PQ_DEFAULT_RUNS = 15` and `BENCH_PQ_DEFAULT_ITERATIONS = 100`. Explicitly supplied values are honoured. These overrides are tracked via `runs_explicit` and `iterations_explicit` in `bench_options_t`.

**`--gate` without `--per-qubit`:** The `--gate` flag has no effect in sweep mode. If specified without `--per-qubit`, a warning is printed to stderr and the flag is silently ignored.

**Gate name matching:** `--gate NAME` is a case-insensitive match against the gate name strings in the gate table: `"X"`, `"H"`, `"Z"`, `"CX"`, `"SWAP"` (uses `strcasecmp`).

**`--per-qubit --pgfplots` combined mode:** When both flags are set, console and CSV output are suppressed; instead, `pgfplots_perq_data_t` accumulates per-target results and `print_pgfplots_perq_output()` emits all tables at the end. If `--csv` is also set, it is ignored with a warning to stderr. See the Per-Qubit Benchmark section below for full details on the emitted tables.

All values are validated at startup; invalid values print to stderr and exit with code 1. Sweep-mode defaults are defined as `BENCH_DEFAULT_*` macros in `bench.h`. `BENCH_MAX_RUNS = 200` is defined in `bench.h` and used by both the validation check and stats computation.

**`blas_dense` and `naive_loop`** are skipped for qubit counts above 8 (hardcoded: `run_dense = (qubits <= 8)`). This cutoff is not derived from `--max-qubits` and is not configurable.

**`naive_loop`** runs at 1/10th the iteration and warmup counts relative to `--iterations`/`--warmup` (hardcoded in `bench_main.c`; 1Q gates only — `naive_loop` is not run for 2Q gates). The reported speedup is scaled by ×10 to make it comparable to other methods (1Q gates only).

### Implementations Benchmarked (bench_gate)

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

### Gates Benchmarked (bench_gate)

**Single-qubit** (all methods): X, H, Z.

**Two-qubit** (no `naive_loop`):

| Gate | Tiled | Notes |
|------|-------|-------|
| CX | Yes | |
| SWAP | Yes | |

### Statistical Model (bench_gate)

#### Sweep Mode: Structure of a Single Measurement

For each `(gate, qubits, method)` triple:

1. `warmup` untimed gate applications are performed, cycling through all qubit targets (`i % qubits`) or all qubit pairs (`i % num_pairs`), to bring code paths and data into cache.
2. `runs` independent timed runs are performed. Each run applies the gate `iterations × sweep_size` times and is measured as a single wall-clock interval.
3. `bench_compute_stats()` computes statistics from the `runs` wall-clock measurements.

#### Per-Qubit Mode: Structure of a Single Measurement

For each `(gate, qubits, target, method)` tuple:

1. **Calibration (two-phase adaptive probe):** Applied identically to all methods, using target qubit 0 (or pair `(q1s[0], q2s[0])` for 2Q gates). The resulting iteration count applies to all target positions for that method.
   - **Phase 1 — single-shot probe:** The gate is invoked once (1 iteration, 1 run). If the returned `time_per_gate_us / 1000.0 >= MIN_MEASUREMENT_MS = 1.0 ms`, the gate is slow enough on its own; `probe_iters` is set to 1. Otherwise `probe_iters = BENCH_PQ_PROBE_ITERS = 100`.
   - **Phase 2 — multi-run probe:** `BENCH_PQ_PROBE_RUNS = 5` runs of `probe_iters` iterations each are timed; the minimum elapsed time across those runs is taken. `calibrate_iterations()` then scales the iteration count up until a single run would exceed `MIN_MEASUREMENT_MS`. If the probe elapsed time rounds to 0 ns, `BENCH_PQ_FALLBACK_ITERATIONS = 10000` is used. The final calibrated count is capped at `BENCH_PQ_MAX_CALIBRATED_ITERATIONS = 1,000,000`.
2. **Warmup:** The gate callback is invoked in batches until at least `BENCH_PQ_WARMUP_MS = 100 ms` of wall-clock time has elapsed. Batch size scales with qubit count:

   | Qubit count | Warmup batch size |
   |-------------|-------------------|
   | ≤ 4 | 128 calls per batch |
   | ≤ 7 | 16 calls per batch |
   | ≤ 10 | 4 calls per batch |
   | > 10 | 1 call per batch |

3. `runs` independent timed runs are recorded. Each run applies the gate `iterations` times to the same target.
4. `bench_compute_stats()` computes statistics from the `runs` measurements.

#### `sweep_size` (Sweep Mode)

`sweep_size` is the number of gate applications per timed iteration:

- **1Q gates**: `sweep_size = qubits` (one application per qubit, target 0..n-1)
- **2Q gates**: `sweep_size = qubits * (qubits - 1) / 2` (one application per unordered pair `(q1, q2)` with `q1 < q2`)

In per-qubit mode, each `_at` call targets a single qubit or qubit pair; `sweep_size` is always 1.

#### 2Q Pair-Building in Per-Qubit Mode

| Gate | `is_symmetric` | Function | Pair count | Rationale |
|------|---------------|----------|------------|-----------|
| CX | 0 | `build_all_ordered_pairs()` | `N*(N-1)` | Control and target are distinct roles |
| SWAP | 1 | `build_all_pairs()` | `N*(N-1)/2` | Order is irrelevant |

In **sweep mode**, `SET_SWEEP` always uses `sw2q = N*(N-1)/2` for all 2Q gates regardless of `is_symmetric`.

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

### Adapter Function Signatures (bench_gate)

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

**QuEST (`bench_quest.c` — `WITH_QUEST` only):**

The per-qubit public API exposed in `bench.h` uses opaque state management: `bench_quest_perq_1q` and `bench_quest_perq_2q` create and destroy a `Qureg` internally, so callers in `bench_main.c` never need `Qureg` or `quest.h` in scope. This avoids QuEST v4 header C++ incompatibilities when `quest.h` is included transitively.

The lower-level `bench_quest_at` / `bench_quest_2q_at` functions (which take an explicit `Qureg*`) remain internal to `bench_quest.c` and are not declared in `bench.h`.

```c
/* Probe-calibrate-measure loop over all targets — Qureg managed internally */
void bench_quest_perq_1q(qubit_t qubits, const char *gate_name,
                          int iters_explicit, int explicit_iters, int runs,
                          bench_result_perq_t *out);

void bench_quest_perq_2q(qubit_t qubits, const char *gate_name,
                          const qubit_t *q1s, const qubit_t *q2s, int n_pairs,
                          int iters_explicit, int explicit_iters, int runs,
                          bench_result_perq_t *out);
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

### Output Formats (bench_gate)

#### Console (Sweep Mode, default)

For each qubit count, for each gate: one result line per method showing method name, `time_ms ± time_ms_std`, CV, and ops/sec. Columns: method name (14 chars), mean time (ms), std dev (ms), CV (%), ops/sec.

With `--verbose`, each result line is extended with: min time (ms), median time (ms), memory footprint (bytes).

After each gate block: a speedup summary line (packed vs other methods) and a memory comparison line (values auto-scaled to KB or MB, with packed vs dense savings percentage).

#### Console (Per-Qubit Mode)

`bench_print_perq_console()` in `bench_main.c` is called once per method per `(gate, qubits)` block.

**1Q gates:** Header `--- <gate> | <method> | <N> qubits ---` followed by one row per target qubit. Columns: `target`, `mean(us/g)`, `median(us/g)`, `min(us/g)`, `CV(%)`, `ops/sec`. Results with CV > 50% are suffixed `(noisy)`.

**2Q gates:** A single aggregated summary line: `min=X.XXX us (qA,qB)  mean=X.XXX us  max=X.XXX us (qC,qD)`. No per-pair rows are printed to console.

#### CSV (Sweep Mode: `--csv`)

One row per `(gate, qubits, method)` triple.

Header: `qubits,dim,gate,method,runs,iterations,sweep_size,time_ms_mean,time_ms_std,time_ms_min,time_ms_median,time_ms_cv,time_per_gate_ms,ops_per_sec,memory_bytes`

#### CSV (Per-Qubit Mode: `--per-qubit --csv`)

One row per `(gate, qubits, target, method)` tuple.

Header: `qubits,dim,gate,method,runs,iterations,target,target2,is_2q,time_ms_mean,time_ms_std,time_ms_min,time_ms_median,time_ms_cv,time_per_gate_us,time_per_gate_us_median,time_per_gate_us_min,ops_per_sec,memory_bytes`

Notes:
- `time_per_gate_us` and related fields are in **microseconds**.
- `target2` is emitted as `-1` for 1Q gates (when `target2 == QUBIT_NONE`).
- Rows where `time_ms_mean <= 0.0` are silently skipped.

**Bandwidth not reported.** No bytes-per-second column is included; the bandwidth formula differs by a different factor for every backend.

#### pgfplots (`--pgfplots`, Sweep Mode Only)

For each of the 5 gates, five tables are emitted: mean per-gate time, 95% CI half-width, min per-gate time, median per-gate time, CV (all timing in µs). After all per-gate tables, one memory table in MB. `naive_loop` is never invoked and is excluded from the pgfplots method index. Missing values are emitted as `nan`.

The 95% CI formula: `t_crit(df) × (σ / √n) × 1000 / (iterations × sweep_size)`, where t-critical is from a hard-coded 49-entry table; for df ≥ 50, 1.96 is used.

---

## bench_circuit — Circuit-Level Benchmark

`bench_circuit` measures end-to-end density matrix circuit execution throughput. All backends receive the identical gate sequence (same `circuit_op_t` array), operate in density matrix mode, and start each timed run from the maximally mixed state ρ = I/dim.

### Circuit Abstraction

Defined in `bench_circuit.h` and implemented in `bench_circuit.c`.

#### `circuit_gate_t`

Enum of gate types used by the three circuit constructors. Only gates actually used are defined:

| Value | Name | Description |
|-------|------|-------------|
| 0 | `GATE_H` | Hadamard |
| 1 | `GATE_X` | Pauli X |
| 2 | `GATE_Z` | Pauli Z |
| 3 | `GATE_RX` | X-rotation (parametric) |
| 4 | `GATE_RY` | Y-rotation (parametric) |
| 5 | `GATE_RZ` | Z-rotation (parametric) |
| 6 | `GATE_CX` | Controlled-NOT |

All dispatch switches must handle all 7 cases; unimplemented cases print to stderr.

#### `circuit_op_t`

A single gate operation:

```c
typedef struct {
    circuit_gate_t gate_type;
    qubit_t        q0;      /* target qubit (or control for CX) */
    qubit_t        q1;      /* second qubit for 2Q gates; 0 otherwise */
    double         angle;   /* rotation angle; 0.0 for non-parametric gates */
} circuit_op_t;
```

#### `circuit_t`

A flat, dynamically-allocated gate sequence:

```c
typedef struct {
    circuit_op_t *ops;       /* heap-allocated array */
    int           n_ops;     /* number of gates in the circuit */
    int           cap;       /* allocated capacity (internal) */
    int           n_qubits;  /* circuit width */
    const char   *name;      /* static string, not copied */
} circuit_t;
```

The ops array grows with doubling realloc. `name` points to a static string (e.g., `"qaoa_p1"`); callers must not free it.

#### `runner_fn` typedef

```c
typedef void (*runner_fn)(const circuit_t *c, void *ctx, int iterations);
```

Used by `calibrate_iterations()` in `bench_circuit_main.c`. A runner receives a circuit, a backend-specific context (cast to the appropriate `*_circuit_ctx_t`), and executes the circuit `iterations` times. The runner does not reinitialize state between iterations. Context is allocated and owned by the caller via the `*_circuit_ctx_alloc()` / `*_circuit_ctx_free()` pairs.

#### Context Structs

| Struct | Backend | Key fields |
|--------|---------|------------|
| `qlib_circuit_ctx_t` | qlib packed/tiled | `state_t *s`, `state_type_t type` |
| `quest_circuit_ctx_t` | QuEST (`WITH_QUEST`) | `void *qureg` (opaque `Qureg*`) |
| `qulacs_circuit_ctx_t` | Qulacs (`WITH_QULACS`) | `void *rho` (opaque `CTYPE*`), `size_t dim` |
| `aer_dm_circuit_ctx_t` | Aer-DM (`WITH_AER_DM`) | `void *dm` (opaque `DensityMatrix<double>*`), `int n_qubits` |

All context types are defined in `bench_circuit.h`. The opaque pointer fields avoid pulling backend-specific C++ headers into C translation units.

### Circuit Constructors

Implemented in `bench_circuit.c`. All constructors return a `circuit_t` with `n_ops == 0` on allocation failure.

#### `circuit_qaoa_p1(int n, double gamma, double beta)` — QAOA p=1 on cycle graph

Structure: `n` H gates + `n` RZZ(gamma) decompositions + `n` RX(2·beta) gates.

RZZ(gamma) is decomposed as: `CX(i, i+1 mod n)` · `Rz(i+1 mod n, 2·gamma)` · `CX(i, i+1 mod n)` — three gates per edge.

**Total gate count: 5n.** Enforced by `assert(c.n_ops == 5 * n)`.

Default parameters: `gamma = π/4`, `beta = π/4` (overridable via `--gamma` and `--beta` CLI flags).

#### `circuit_vqe_hea(int n, int depth)` — Hardware-Efficient VQE Ansatz

`depth` layers of: `n` Ry(π/4) gates + CX ladder (0→1, 1→2, ..., n-2→n-1). If `depth < 0`, defaults to `floor(n/2)`.

**Total gate count: `depth × (2n − 1)`.** Enforced by assert.

The VQE angle (π/4) is hardcoded for reproducibility. It is not CLI-overridable.

#### `circuit_qv(int n, int seed)` — Quantum Volume

`n` must be even. Constructs `n` layers; each layer: Fisher-Yates random permutation of n qubits (LCG-seeded), then `n/2` qubit pairs each receiving:

- 3 Ry(random angle) on `qa`
- 3 Ry(random angle) on `qb`
- CX(qa, qb) · CX(qb, qa) · CX(qa, qb)
- 3 Rz(random angle) on `qa`
- 3 Rz(random angle) on `qb`

15 gates per pair. **Total gate count: `15 × (n/2) × n = 15n²/2`.** Enforced by assert.

This exercises the same gate mix as a KAK decomposition (9 single-qubit rotations + 3 CNOT gates per pair) but angles are LCG-generated, not KAK-computed interaction parameters.

The RNG is a 64-bit LCG with constants `(6364136223846793005, 1442695040888963407)`. Angles are drawn in `[0, π]`. The same seed produces the same circuit across runs and machines.

### Backends Benchmarked (bench_circuit)

Five competitors, all in density matrix mode:

| Method string | Source file | Notes |
|---------------|-------------|-------|
| `qlib_packed` | `bench_circuit_qlib.c` | Always present |
| `qlib_tiled` | `bench_circuit_qlib.c` | Always present |
| `quest` | `bench_circuit_quest.c` | `WITH_QUEST` only |
| `qulacs` | `bench_circuit_qulacs.cpp` | `WITH_QULACS` only; uses csim low-level `dm_*_gate` functions |
| `qiskit_aer_dm` | `bench_circuit_aer_dm.cpp` | `WITH_AER_DM` only; rotation gates construct per-op matrices |

Each backend implements the same `bench_circuit_<name>()` function signature:
```c
bench_result_t bench_circuit_<name>(const circuit_t *c, int iterations, int warmup, int runs);
```

The `warmup` parameter is accepted but ignored by all backends; warmup is time-based (see Measurement Methodology).

### Measurement Methodology (bench_circuit)

All backends follow the same protocol:

1. **State initialization:** Allocate the density matrix once before the warmup loop. Initialize to the maximally mixed state ρ = I/dim: zero the buffer, then set diagonal entries to `1/dim`. This is done via `memset` + a loop over diagonal indices; it is not |+><+|/dim.

2. **Time-based warmup:** Run full circuit passes continuously until `BENCH_PQ_WARMUP_MS = 100 ms` of wall-clock time has elapsed. State is not reinitialized during warmup.

3. **Per-run reinitialization:** Before each of the `runs` timed runs, the state is reset to ρ = I/dim. State is NOT reinitialized between the `iterations` circuit applications within a single run. This measures **sustained throughput** (steady-state gate-application rate), not single-shot latency.

4. **Timing barriers:** Three `bench_timing_barrier()` calls per timed run (immediately before the start-clock read, immediately after the start-clock read, and immediately after the timed loop), matching the `bench_gate` per-qubit protocol.

5. **Statistics:** `bench_compute_stats()` on the `runs` wall-clock measurements, producing mean, std dev, min, median, CV.

6. **Iteration calibration:** By default (`--iterations 0`), each backend independently calibrates its iteration count before the main benchmark loop. The calibration runs `BENCH_PQ_PROBE_ITERS = 100` iterations for `BENCH_PQ_PROBE_RUNS = 5` probe rounds, takes the minimum elapsed time, and scales to target at least `MIN_MEASUREMENT_MS = 1.0 ms` per timed run. Result is clamped to `[1, BENCH_PQ_MAX_CALIBRATED_ITERATIONS]`. If the probe elapsed time is < 1e-6 ms, `BENCH_PQ_FALLBACK_ITERATIONS = 10000` is used.

   For `qlib_packed` and `qlib_tiled`, calibration uses the `runner_fn` / context pattern. For QuEST, Qulacs, and Aer-DM, calibration calls the respective `bench_circuit_*()` function directly for `BENCH_PQ_PROBE_RUNS` probe runs with `BENCH_PQ_PROBE_ITERS` iterations each.

   Calibration is per-backend and per-`(circuit, qubit count)`. Different backends may calibrate to different iteration counts; compare using `time_per_gate_ms`, not raw `time_ms`.

7. **`sweep_size` field:** In `bench_result_t`, `sweep_size` is set to `c->n_ops` (the circuit's gate count) by circuit runners. `time_per_gate_ms` = `time_ms / (iterations × n_ops)`. `ops_per_sec` counts individual gate applications per second.

### Memory Reporting (bench_circuit)

All backends report the theoretical state storage footprint only, not framework overhead.

| Method | Reported `memory_bytes` |
|--------|------------------------|
| `qlib_packed` | `bench_packed_size(qubits)` = `dim*(dim+1)/2 × sizeof(cplx_t)` |
| `qlib_tiled` | `bench_tiled_size(qubits)` = tile-pair count × `TILE_SIZE × sizeof(cplx_t)` |
| `quest` | `dim² × sizeof(double) × 2` |
| `qulacs` | `dim² × sizeof(CTYPE)` |
| `qiskit_aer_dm` | `dim² × sizeof(std::complex<double>)` |

### Correctness Verification (bench_circuit, debug builds only)

Active when `NDEBUG` is not defined (i.e., Debug builds) and `n_qubits <= 4`. Called automatically after each `run_all_methods()` invocation.

Each included backend applies the circuit once from ρ = I/dim and extracts the resulting full density matrix as a `dim × dim` array of `cplx_t`. All pairs of extracted matrices are compared using the **Frobenius norm**, tolerance `1e-6`. Results are printed to stderr (to avoid polluting stdout CSV/pgfplots output).

**Layout note for QuEST:** `cpuAmps` is stored column-major (`cpuAmps[col * dim + row]`); the extraction transposes into row-major before comparison.

**Qulacs limitation:** The `csim` gate functions are declared in C++ headers and cannot be called from the C translation unit where verification runs. Full-matrix cross-check for Qulacs is skipped; a TODO note is present in `bench_circuit_main.c` indicating a `qulacs_apply_circuit_once()` helper must be added to `bench_circuit_qulacs.cpp` to enable it.

**Aer-DM limitation:** `AER::QV::DensityMatrix<double>` is opaque from C. Full-matrix cross-check for Aer-DM is skipped; an `aer_dm_extract_matrix()` helper would need to be added to `bench_circuit_aer_dm.cpp` to enable it.

Only `qlib_packed`, `qlib_tiled`, and (when `WITH_QUEST`) `quest` currently participate in pairwise Frobenius comparison.

### CLI Reference (bench_circuit)

```bash
./benchmark/bench_circuit                          # All circuits, defaults
./benchmark/bench_circuit --qaoa                   # QAOA p=1 only
./benchmark/bench_circuit --vqe                    # VQE-HEA only
./benchmark/bench_circuit --qv                     # Quantum Volume only
./benchmark/bench_circuit --qubits-min 4 --qubits-max 12 --qubits-step 2
./benchmark/bench_circuit --runs 20
./benchmark/bench_circuit --iterations 50          # Fixed; skips calibration
./benchmark/bench_circuit --gamma 0.5 --beta 0.3   # QAOA parameters
./benchmark/bench_circuit --seed 123               # QV RNG seed
./benchmark/bench_circuit --csv > results.csv
./benchmark/bench_circuit --pgfplots > results.dat
./benchmark/bench_circuit --verbose
./benchmark/bench_circuit --help
```

All flags:

| Flag | Default | Purpose |
|------|---------|---------|
| `--qubits-min N` | 4 | Minimum qubit count |
| `--qubits-max N` | 12 | Maximum qubit count |
| `--qubits-step N` | 2 | Qubit count increment |
| `--runs N` | 10 | Independent timing runs (1..`BENCH_MAX_RUNS`) |
| `--seed N` | 42 | RNG seed for QV circuit (must be ≥ 0) |
| `--qaoa` | — | Run only QAOA p=1 |
| `--vqe` | — | Run only VQE-HEA |
| `--qv` | — | Run only Quantum Volume |
| `--all` | — | Run all circuits (same as default) |
| `--gamma FLOAT` | π/4 | QAOA cost Hamiltonian parameter |
| `--beta FLOAT` | π/4 | QAOA mixer parameter |
| `--iterations N` | 0 | Gate applications per timed run; 0 = calibrate |
| `--csv` | off | CSV output to stdout |
| `--pgfplots` | off | pgfplots-compatible output to stdout |
| `--verbose` | off | Extended per-result info (min, median, iterations, memory) |
| `--help` | — | Print usage and exit |

**Default circuit selection:** If none of `--qaoa`, `--vqe`, `--qv`, `--all` is specified, all three circuits run (equivalent to `--all`).

**QV and odd qubit counts:** `circuit_qv` requires even n. When `--qv` is active and the current qubit count is odd, QV is skipped with a warning to stderr; other circuits at that qubit count still run.

**`--qubits-min` validation:** Values below 2 are clamped to 2. (`n=1` would cause `circuit_rzz(0, 0)` which sets control == target for CX, invoking undefined behavior.)

### CSV Output Format (bench_circuit)

One row per `(circuit, method, n_qubits)` triple.

Header:
```
circuit,method,n_qubits,n_ops,iterations,runs,mean_ms,ci95_ms,min_ms,median_ms,cv
```

Fields:
- `circuit`: circuit name string (`qaoa_p1`, `vqe_hea`, `qv`)
- `method`: backend method string
- `n_qubits`: qubit count
- `n_ops`: total gate count of the circuit (`sweep_size` in `bench_result_t`)
- `iterations`: calibrated (or fixed) circuit applications per timed run
- `runs`: number of timed runs
- `mean_ms`: mean run time (ms)
- `ci95_ms`: 95% CI half-width = `1.96 × std / sqrt(runs)` (0 if `runs == 1`)
- `min_ms`: minimum run time (ms)
- `median_ms`: median run time (ms)
- `cv`: coefficient of variation (%)

Note: `ci95_ms` uses the normal approximation (1.96) rather than t-critical lookup. No `ops_per_sec` or `time_per_gate_ms` column is emitted; derive them as needed from `mean_ms`, `n_ops`, `iterations`.

### pgfplots Output Format (bench_circuit)

One table per circuit (preceded by `# CIRCUIT: <name>`). Each row is one `(n_qubits, method)` pair.

Header:
```
qubits method time_us_per_gate_mean time_us_per_gate_min time_us_per_gate_median time_us_per_gate_cv_pct memory_bytes
```

`time_us_per_gate_mean` = `time_per_gate_ms × 1000` (µs). `time_us_per_gate_min` and `_median` are derived analogously from `time_ms_min` and `time_ms_median` divided by `iterations × n_ops`.

### Console Output (bench_circuit)

For each `(circuit, n_qubits)` block:
- Header line: `Circuit: <name>  qubits: N  n_ops: M`
- One result line per backend: `  <method>  mean ± std ms (CV%)  ops/sec`
- With `--verbose`: additionally `[min X ms  med X ms  iters N  mem B bytes]`
- Speedup line: `Speedup: X.XXx tiled vs packed` (console only, not CSV/pgfplots)
- Separator: `---...---`

In Debug builds, correctness verification output appears on stderr after the speedup line.

The header block (printed once at startup in non-CSV/non-pgfplots mode) includes: date/time, build type (Debug/Release), OpenMP thread count, qubit range, and active backends.

---

## Shared Infrastructure

### Timing and Memory Methodology

All timing uses `clock_gettime(CLOCK_MONOTONIC)` via `bench_time_ns()`. The `bench_timing_barrier()` call is architecture-aware: serializing instruction on ARM64 (`isb`) and x86-64 (`lfence`); compiler-only barrier (`asm volatile("" ::: "memory")`) on other platforms.

**Linux note:** For reliable results, set the CPU governor to `performance` mode (`cpupower frequency-set -g performance`).

### Memory Reporting (bench_gate)

| Method | Measurement |
|--------|-------------|
| `qlib_packed` | `dim*(dim+1)/2 × sizeof(cplx_t)` via `bench_packed_size()` |
| `qlib_tiled` | tile-pair count × `TILE_SIZE × sizeof(cplx_t)` via `bench_tiled_size()` |
| `blas_dense` | `dim² × sizeof(cplx_t)` via `bench_dense_size()` |
| `naive_loop` | same as `blas_dense` |
| QuEST | `dim² × sizeof(double) × 2` |
| Qulacs | `dim² × sizeof(CTYPE)` |
| `qiskit_aer_dm` | `dim² × sizeof(std::complex<double>)` |

### Huge Pages

State allocations at or above `BENCH_HUGEPAGE_THRESHOLD = 1 GB` use `bench_alloc_huge()` / `bench_free_huge()`, which request huge-page-backed memory via `mmap` on Linux (with `MADV_HUGEPAGE`) and macOS. Below the threshold, `bench_alloc_huge()` calls `malloc`. `bench_free_huge()` must be called with the same `n` as `bench_alloc_huge()`.

### External Framework Adapter Notes

**QuEST** uses an init-once pattern: `bench_quest_init()` calls `initQuESTEnv()` once at program start; `bench_quest_cleanup()` calls `finalizeQuESTEnv()` once at exit. This applies to both `bench_gate` and `bench_circuit`. QuEST v4 API differences from v3: rotation gates renamed to `applyRotateX/Y/Z`; `createDensityQureg()` takes only `numQubits` (no `QuESTEnv` argument); no `initMixedState()` for ρ = I/dim (implemented manually via `initBlankState()` + diagonal writes); `cpuAmps` layout is column-major.

**Qulacs** has no global state or init/finalize. `bench_circuit_qulacs.cpp` uses the low-level `csim` `dm_*_gate` functions directly, bypassing the higher-level `cppsim` gate object layer.

**Aer-DM** has no global state. The `bench_circuit_aer_dm.cpp` adapter precomputes fixed-gate matrices (`H`, `X`, `Z`, `CX`) once per benchmark call. Rotation gates (`Rx`, `Ry`, `Rz`) construct their 2×2 matrices per operation because angles vary per gate instance. C++ exceptions are caught at the `extern "C"` boundary; a zero `time_ms` result indicates failure (e.g., `std::bad_alloc` at large qubit counts).

## Per-Qubit Benchmark

`--per-qubit` mode measures steady-state single-target throughput: it applies the same gate to the same fixed target qubit many times on a single state that is re-initialized before each target's warmup phase. This is a deliberate isolation measurement that exposes whether gate cost depends on target position in the memory layout of each representation. It is **not** circuit-representative throughput; do not compare `time_per_gate_us` (per-qubit mode, µs) directly to `time_per_gate_ms` (sweep mode, ms) without qualification.

### CLI Reference

| Flag | Default | Notes |
|------|---------|-------|
| `--per-qubit` | off | Enables single-target throughput mode |
| `--gate NAME` | all gates | Restricts run to one named gate; no effect without `--per-qubit` |
| `--runs N` | `BENCH_PQ_DEFAULT_RUNS` | Override auto-default for per-qubit mode |
| `--iterations N` | `BENCH_PQ_DEFAULT_ITERATIONS` | Override auto-default for per-qubit mode |

**`--pgfplots` with `--per-qubit`:** Supported. Emits pgfplots `.dat` tables to stdout. Results are grouped by `(gate, qubit_count)`: within each block the x-axis is target position (1Q: integer qubit index; 2Q: `q1_q2` pair label), not qubit count. Five table types are emitted per block: `time_per_gate_us` (mean), `ci95_per_gate` (95% CI half-width, `t_crit(runs-1) × std / sqrt(runs) / iters × 1000`), `min_per_gate`, `median_per_gate`, and `cv`. Label and method columns are all 14 characters wide. If `--csv` is also set, it is ignored with a warning to stderr.

**Implementation:** `print_pgfplots_perq_output()` calls `emit_perq_table()` five times (once per field). Gate symmetry (`is_symmetric`) is read directly from the file-scope `gates_2q_perq[]` table. Pair reconstruction is handled by `reconstruct_pair()`, whose return value is asserted at each call site.

### Known Limitations

| ID | Description |
|----|-------------|
| L1 | **State initialization differs between backends.** qlib and BLAS use `bench_fill_random()`; Qulacs and Aer-DM use `std::mt19937_64`. QuEST's `initDebugState` produces a deterministic state unrelated to the others. At low qubit counts where the matrix fits in L1/L2 cache, different data values can affect SIMD fast-paths. Accept this for cross-framework comparison; document it when citing. |
| L2 | **Measurement semantics vs. circuit workloads.** Per-qubit `time_per_gate_us` is not directly comparable to sweep `time_per_gate_ms`. The sweep benchmark is a better proxy for circuit-like workloads. |
| L3 | **Bandwidth not reported.** The `theoretical_bw_gbs` field is removed because the formula `2 * memory_bytes / time_per_gate` is wrong by different factors for each backend. If bandwidth analysis is needed, measure it with STREAM on the same hardware and compare independently. |
| L4 | **NUMA and thread affinity.** On multi-socket machines, target-qubit cost differences may reflect NUMA topology, not representation cost. Use `numactl --interleave=all` on Linux for publication results. Not applicable on Apple Silicon. |
| L5 | **No IQR outlier filtering.** With 15 runs, IQR trimming is feasible but not implemented. A `(noisy)` flag at CV > 50% is the only outlier indicator. |
| L6 | **Variable iteration counts across methods.** Per-qubit mode calibrates once per (gate, qubits, method). Different methods will have different `iterations` counts. Use `time_per_gate_us` (not total run time) for all cross-method comparisons. |
| L7 | **Calibration-on-target-0 bias.** The probe runs on target 0 only. For packed storage representations, gates on higher-index targets touch more of the state buffer, so the calibrated iteration count may be too high for high targets. Runs for high-index targets will be shorter in wall time than runs for target 0. |
| L8 | **No IQR outlier filtering in pgfplots mode.** CI95 is computed from std dev and run count stored in `perq_cell_t`; outlier trimming is not applied. Use `--per-qubit --csv` for full per-target statistical detail. |

## Status

### Reproducibility

`bench_gate` state data is filled with random values via `bench_fill_random()` (xoshiro256** seeded from OS entropy). No fixed seed is set; gate arithmetic has no data-dependent branches, so state values do not affect timing.

`bench_circuit` state is initialized to ρ = I/dim (maximally mixed), not random. QV circuit angles are determined entirely by `--seed` (LCG); the same seed produces the same circuit deterministically.

### Per-Qubit Confidence Interval Note

At the default `BENCH_PQ_DEFAULT_RUNS = 15` runs, the 95% CI on the Bessel-corrected standard deviation spans approximately [0.74σ, 1.83σ] (chi-squared distribution, 14 degrees of freedom). For publication-quality intervals, use at least 30–50 runs. Results with CV > 50% are flagged as `(noisy)` in console output; no IQR-based filtering is applied.

### Constants

All constants are defined in `bench.h`.

| Constant | Value | Applies to | Meaning |
|----------|-------|------------|---------|
| `BENCH_DEFAULT_MIN_QUBITS` | 2 | bench_gate sweep | Default `--min-qubits` |
| `BENCH_DEFAULT_MAX_QUBITS` | 12 | bench_gate sweep | Default `--max-qubits` |
| `BENCH_DEFAULT_STEP` | 1 | bench_gate sweep | Default `--step` |
| `BENCH_DEFAULT_ITERATIONS` | 1000 | bench_gate sweep | Default `--iterations` in sweep mode |
| `BENCH_DEFAULT_WARMUP` | 100 | bench_gate sweep | Default `--warmup` in sweep mode |
| `BENCH_DEFAULT_RUNS` | 10 | bench_gate sweep | Default `--runs` in sweep mode |
| `BENCH_MAX_RUNS` | 200 | Both | Hard upper bound on `--runs`; also bounds the stats sort buffer |
| `BENCH_PQ_DEFAULT_RUNS` | 15 | bench_gate per-qubit | Default `--runs` override in per-qubit mode |
| `BENCH_PQ_DEFAULT_ITERATIONS` | 100 | bench_gate per-qubit | Starting iteration count before calibration |
| `BENCH_PQ_PROBE_ITERS` | 100 | Both | Iterations per probe run during calibration |
| `BENCH_PQ_PROBE_RUNS` | 5 | Both | Probe runs during calibration; minimum elapsed time is used |
| `BENCH_PQ_WARMUP_MS` | 100.0 | Both | Minimum warmup duration (ms) before timed runs begin |
| `MIN_MEASUREMENT_MS` | 1.0 | Both | Target minimum duration for a single timed run |
| `BENCH_PQ_FALLBACK_ITERATIONS` | 10000 | Both | Fallback iteration count when probe elapsed time rounds to 0 ns |
| `BENCH_PQ_MAX_CALIBRATED_ITERATIONS` | 1,000,000 | Both | Upper bound on calibrated iteration count |
| `QUBIT_NONE` | `(qubit_t)-1` | bench_gate per-qubit | Sentinel for absent second qubit in 1Q gate results |
| `MAX_PAIRS` | 256 | bench_gate | Maximum unordered/ordered pairs for `build_all_pairs()` / `build_all_ordered_pairs()`; covers 16×15 = 240 ordered pairs |
| `MAX_TARGETS` | 512 | bench_gate per-qubit | `2 × MAX_PAIRS`; upper bound on per-qubit result array size |
| `MAX_QUBITS_FOR_ORDERED_PAIRS` | 16 | bench_gate per-qubit | Maximum qubit count accepted by `build_all_ordered_pairs()` |
| `BENCH_HUGEPAGE_THRESHOLD` | `1 << 30` (1 GB) | Both | Size above which `bench_alloc_huge()` requests huge-page-backed memory |

---

## Limitations and Known Gaps

| ID | Limitation |
|----|------------|
| L1 | **State initialization differs between backends (bench_gate).** qlib states are filled with `bench_fill_random()` (xoshiro256** from OS entropy); QuEST uses its own debug-state initializer; Qulacs and Aer DM use direct buffer fills. Initial state values do not affect timing (no data-dependent branches in gate arithmetic), but initialization cost is excluded from all timed regions. |
| L2 | **Per-qubit mode measures steady-state single-target throughput, not circuit-representative throughput.** Each benchmark applies the same gate to the same target in a tight loop. Results are not directly comparable to sweep-mode `time_per_gate_ms`, which averages over all targets. |
| L3 | **NUMA and thread affinity not controlled.** On multi-socket systems, memory allocation and thread scheduling may land on different NUMA nodes across runs, introducing variance not captured by CV. |
| L4 | **No IQR outlier filtering.** Only a CV > 50% flag is emitted. The minimum and median columns are more robust in the presence of scheduler spikes. |
| L5 | **Variable iteration counts across methods (per-qubit mode and bench_circuit).** Calibration runs independently per method. Use `time_per_gate_us` (per-qubit) or `time_per_gate_ms` (circuit) for cross-method comparison; do not compare raw `time_ms_mean` values across methods. |
| L6 | **Calibration-on-target-0 bias (bench_gate per-qubit mode).** Iteration count calibration uses target qubit 0. For backends where per-qubit cost varies with target index (e.g., tile-boundary effects in the tiled backend), the calibrated count may result in timed runs that are too short for high-index targets. |
| L7 | **Gate matrix construction time is not separated from gate application time.** In the Aer-DM circuit runner, rotation gate matrices are reconstructed per operation. This overhead is inherent to the Aer-DM API and is measured as part of backend cost. |
| L8 | **Measurement (non-unitary) operations are not benchmarked.** |
| L9 | **No thread-count sweep.** OpenMP uses all available threads by default; control via `OMP_NUM_THREADS`. |
| L10 | **Qulacs and Aer-DM excluded from bench_circuit correctness verification.** Qulacs csim headers are C++ only; Aer-DM DensityMatrix is opaque from C. Only qlib_packed, qlib_tiled, and QuEST (when compiled in) participate in pairwise Frobenius norm checks. Adding helpers (`qulacs_apply_circuit_once()` in `bench_circuit_qulacs.cpp`, `aer_dm_extract_matrix()` in `bench_circuit_aer_dm.cpp`) would enable full cross-check. |
| L11 | **bench_circuit CSV omits ops_per_sec and time_per_gate_ms columns.** Derive from `mean_ms / (iterations × n_ops)` if needed. |
| L12 | **bench_circuit ci95_ms uses normal approximation (1.96) rather than t-critical.** For small run counts (< 30), the true t-critical value is larger; CI widths are understated. |
