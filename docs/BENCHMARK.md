# Benchmark Module

Technical spec for `bench_mixed`, which compares qlib's mixed-state gate implementations against each other and against external frameworks.

## File Layout

```
benchmark/
├── CMakeLists.txt              # Release-only build config; auto-detects QuEST + Qulacs
├── README.md                   # Local notes: setup instructions, sample results, methodology
├── include/
│   └── bench.h                 # Types (bench_result_t, bench_options_t), timing helpers, memory measurement, function declarations
└── src/
    ├── bench_mixed.c           # Main benchmark: gate definitions, CLI parsing, output formatting, benchmark loops, blas_dense and naive_loop implementations
    ├── bench_quest.c           # QuEST adapter (compiled only when WITH_QUEST is defined)
    └── bench_qulacs.cpp        # Qulacs adapter (compiled only when WITH_QULACS is defined)
```

`README.md` is a developer convenience file (setup commands, sample results, methodology notes). This spec is the authoritative reference.

## Dependencies (Inputs)

**qlib module (`libq`)**: `bench_mixed` links against `libq` and exercises the public gate API directly — specifically `src/gate/packed/` (packed lower-triangular) and `src/gate/tiled/` (tiled lower-triangular). It inherits BLAS and OpenMP linkage transitively via `libq`'s public CMake dependencies.

**External frameworks** (optional, auto-detected):

| Dependency | Location | Detection | Build |
|------------|----------|-----------|-------|
| QuEST v4 | `extern/QuEST/` | `find_path(quest/include/quest.h)` + `find_library(QuEST)` | Pre-built; see `benchmark/README.md` for setup |
| Qulacs (`csim`) | `extern/qulacs/src/csim/` | Checks for `update_ops_dm.hpp` | Built in-tree as `csim_static`; requires Eigen |
| Eigen | System or Homebrew | `find_path(Eigen/Core)` | Header-only; needed by `csim_static` |
| OpenMP | System or Homebrew libomp | `find_package(OpenMP)` | macOS: auto-configured from `brew --prefix libomp` |

If a framework is absent its columns are omitted from output; the binary still builds.

## Build Integration

Benchmarks require a **Release build**. `CMakeLists.txt` returns immediately (with a status message) if `CMAKE_BUILD_TYPE != Release`.

```bash
# From project root
cmake -B cmake-build-release -DCMAKE_BUILD_TYPE=Release -DBUILD_BENCHMARKS=ON
cmake --build cmake-build-release --target bench_mixed
# Binary lands at: cmake-build-release/benchmark/bench_mixed
```

The `bench` convenience target (requires `BUILD_TESTS=ON`) runs the full test suite first, then the benchmark:

```bash
cmake --build cmake-build-release --target bench
```

CMake prints a configuration summary at the end of the benchmark configure step:

```
-- Benchmark configuration:
--   QuEST:     TRUE/FALSE
--   Qulacs:    TRUE/FALSE
--   OpenMP:    TRUE/FALSE
```

Note: the OpenMP line reflects `OpenMP_CXX_FOUND`, not the generic `OpenMP_FOUND`.

**CMake targets produced:**

| Target | Type | Condition |
|--------|------|-----------|
| `bench_mixed` | Executable | Always (Release only) |
| `csim_static` | Static library | `QULACS_FOUND` |
| `bench` | Custom (ctest --output-on-failure + bench_mixed) | `BUILD_TESTS=ON` |

**Compile definitions:** `_GNU_SOURCE` always; `WITH_QUEST` if QuEST found; `WITH_QULACS` if Qulacs found.

**Custom install paths:** Pass `-DQuEST_DIR=/path` to cmake for a non-default QuEST location. There is no equivalent `QULACS_DIR` hint — the Qulacs path is hardcoded to `${CMAKE_SOURCE_DIR}/extern/qulacs`.

**`bench` target dependencies:** The `bench` custom target depends on `test_state` and `test_gate` specifically (not all test targets). If new test targets are added, this list must be updated manually in `benchmark/CMakeLists.txt`.

## Running

```bash
# From the Release build directory
./benchmark/bench_mixed                                    # Defaults: qubits 2–12, step 1, 1000 iter, 100 warmup
./benchmark/bench_mixed --min-qubits 4 --max-qubits 8     # Qubit range
./benchmark/bench_mixed --step 2                           # Even qubit counts only (2, 4, 6, …)
./benchmark/bench_mixed --iterations 5000                  # Higher precision
./benchmark/bench_mixed --csv > results.csv                # CSV export
./benchmark/bench_mixed --pgfplots > results.dat           # pgfplots export (stdout)
./benchmark/bench_mixed --verbose                          # Per-method memory inline
./benchmark/bench_mixed --help                             # Full usage
```

All CLI flags:

| Flag | Default | Valid range | Purpose |
|------|---------|-------------|---------|
| `--min-qubits N` | 2 | 1–30 | Minimum qubit count |
| `--max-qubits N` | 12 | min–30 | Maximum qubit count |
| `--step N` | 1 | ≥1 | Qubit count increment |
| `--iterations N` | 1000 | ≥1 | Timed iterations per gate/method |
| `--warmup N` | 100 | ≥0 | Warm-up iterations (untimed) |
| `--csv` | off | — | One CSV row per measurement |
| `--pgfplots` | off | — | One timing table per gate + one memory table (all to stdout) |
| `--verbose` | off | — | Show per-method memory inline |

All arguments are validated at startup; invalid values print an error and exit with code 1. Defaults are defined as macros in `bench.h` (`BENCH_DEFAULT_*`).

`blas_dense` and `naive_loop` are skipped for qubit counts above 8 (hardcoded `run_dense = (qubits <= 8)`; not derived from `--max-qubits`, not configurable). `qlib_packed` and `qlib_tiled` run for all qubit counts up to `--max-qubits` with no additional size guard.

`naive_loop` runs at 1/10th the iteration and warmup counts relative to `--iterations`/`--warmup` (hardcoded in `bench_mixed.c`) to keep runtimes tractable for its O(N³) cost.

## Output Format

**Default (console):** For each qubit count, for each gate: one line per method showing method name, time_ms, and ops/sec. After each gate block: a speedup summary line and a memory comparison line (values auto-scaled to KB or MB). Memory is always shown in the default mode; `--verbose` adds per-method memory inline on the method line itself.

**`--csv`:** One row per `(gate, qubits, method)` triple. Header and columns:

```
qubits,dim,gate,method,time_ms,ops_per_sec,memory_bytes,iterations
```

Values formatted as: `%u,%ld,%s,%s,%.6f,%.0f,%zu,%d`. Output to stdout; redirect to file as needed.

**`--pgfplots`:** Output to stdout. Contains one timing table per gate (columns: `qubits`, then one column per method), followed by one memory table (columns: `qubits`, then one column per method, values in MB). `naive_loop` is excluded from pgfplots output (5 methods shown: `qlib`, `qlib_tiled`, `blas`, `quest`, `qulacs`). Missing values (not-built or skipped) are emitted as `nan`.

## Methods Compared

Six methods are indexed in `bench_mixed.c` (enum `M_QLIB` … `M_QULACS`):

| Index | Enum | Method string | Description | Source |
|-------|------|---------------|-------------|--------|
| 0 | `M_QLIB` | `qlib_packed` | Production packed lower-triangular (BLAS-1) | `src/gate/packed/` |
| 1 | `M_QLIB_TILED` | `qlib_tiled` | Production tiled lower-triangular | `src/gate/tiled/` |
| 2 | `M_BLAS` | `blas_dense` | Full N×N matrix, `zgemm` UρU† | `bench_mixed.c` |
| 3 | `M_NAIVE` | `naive_loop` | Element-by-element UρU† loop | `bench_mixed.c` |
| 4 | `M_QUEST` | `quest` | QuEST density-matrix backend | `bench_quest.c` |
| 5 | `M_QULACS` | `qulacs` | Qulacs csim density-matrix backend | `bench_qulacs.cpp` |

QuEST and Qulacs are conditionally compiled; their slots are absent from output when not built.

## Gates Benchmarked

**Single-qubit** (`gates[]` in `bench_mixed.c`): X, H, Z — benchmarked with all available methods.

**Two-qubit** (`gates_2q[]`): Each entry carries a `has_tiled` flag controlling whether the tiled method runs.

| Gate | `has_tiled` | Notes |
|------|-------------|-------|
| CX | 1 | Tiled implementation: `gate_tiled_cx.c` |
| SWAP | 1 | Tiled implementation: `gate_tiled_swap.c` |

Two-qubit gates do not run `naive_loop`. Each 2Q gate iterates over all `n*(n-1)/2` ordered qubit pairs `(q1, q2)` with `q1 < q2` (not just adjacent pairs), so every code path — within-tile, mixed, and tile-level — is exercised. To add a new 2Q gate: append to `gates_2q[]` and set `has_tiled`.

## Key Data Structures

All defined in `bench.h` unless noted:

- **`bench_result_t`**: Single measurement record — qubits, dim, gate name, method name, time_ms, ops/sec, memory bytes, iterations.
- **`bench_options_t`**: Parsed CLI options.
- **`bench_summary_t`** (defined in `bench_mixed.c`): Collects all method results for one `(gate, qubit-count)` combination. Used by `print_speedup_and_memory()` to emit the speedup and memory lines — shared between 1Q and 2Q loops.
- **`pgfplots_data_t`** (defined in `bench_mixed.c`, not `bench.h`): Accumulates results for `--pgfplots` output. Heap-allocated only when `--pgfplots` is passed (zero BSS impact otherwise). Indexed as `time_ms[gate][qubit_config][method]`.

## Timing and Memory Methodology

**Timing:** Wall-clock time. qlib, blas_dense, and naive_loop use `clock_gettime(CLOCK_MONOTONIC)` via `bench_time_ns()` in `bench.h`. The Qulacs adapter uses `std::chrono::high_resolution_clock` instead.

Each measurement:
1. `warmup` untimed iterations, cycling through all qubit targets (`i % qubits`) or all qubit pairs (`i % num_pairs`) so that every code path is cache-warm before timing starts.
2. `iterations` timed iterations: single-qubit gates apply to every qubit 0..n-1; two-qubit gates apply to every ordered pair `(q1, q2)` with `q1 < q2` (`n*(n-1)/2` pairs per iteration).
3. `srand(42)` is called before each `bench_*` invocation for a given `(gate, qubits)` combination, so all methods initialise their state from the same random number sequence.

ops/sec = `(iterations × gates_per_iteration) / (time_ms / 1000)`. For naive_loop the iteration count is `iterations/10`; the displayed speedup is scaled by ×10 to make it comparable.

**Memory:** qlib's own formats (`qlib_packed`, `qlib_tiled`) report the exact theoretical allocation size via `bench_packed_size()` / `bench_tiled_size()` — these are exact because the allocator is known. `blas_dense` and `naive_loop` report `bench_dense_size()` (exact). `blas_dense` and `naive_loop` pre-allocate scratch buffers (`tmp`, `out`) once outside the hot loop; these scratch buffers are **not** included in the reported memory (they are implementation overhead, not state storage). QuEST reports RSS delta (before/after `createDensityQureg` + `initDebugState` to force page faults) measured inside a forked child. Qulacs reports theoretical matrix size (`dim² × 16 bytes`) because RSS delta after `fork()` is unreliable on macOS.

## External Dependencies

| Dependency | Location | Build condition |
|------------|----------|----------------|
| QuEST v4 | `extern/QuEST/` | `QuEST_FOUND` (pre-built; `libQuEST` + `quest/include/quest.h`) |
| Qulacs | `extern/qulacs/` | `QULACS_FOUND` (built as `csim_static` from `src/csim/*.cpp`; needs Eigen) |

Setup instructions (submodule init, macOS Homebrew deps, build commands): see `benchmark/README.md`.

Custom install paths: pass `-DQuEST_DIR=/path` to cmake. No `QULACS_DIR` hint is supported; Qulacs must reside at `extern/qulacs/`.

## Memory Layout Background

The benchmark compares two production storage formats for density matrices (full specification in `docs/STATE_MODULE.md`):

- **Packed** (`MIXED_PACKED`): LAPACK packed lower-triangular, column-major. ~N(N+1)/2 elements. Favors BLAS-1 but has variable-stride access for high target qubits.
- **Tiled** (`MIXED_TILED`): Lower triangle partitioned into `TILE_DIM × TILE_DIM` tiles. Intra-tile layout is row-major. Eliminates the variable-stride bottleneck via contiguous cache-friendly access.

## Status

**In Progress.** Core benchmark (`bench_mixed`) is functional for single-qubit gates (X, H, Z) and two-qubit gates (CX, SWAP). External framework adapters are complete.

Known gaps / TODO:
- No separate timing of gate matrix construction vs. gate application.
- Measurement (non-unitary) operations not benchmarked.
- No thread-count sweep (OpenMP uses all available threads by default; controlled externally via `OMP_NUM_THREADS`).
- Scientific rigor checklist (environment, build config, statistical caveats) documented in `benchmark/README.md` — not yet enforced by the benchmark binary itself.

## Open Questions

1. **`--step` default intent**: The CLI default is `--step 1` (every integer qubit count), confirmed in source. Whether the default should be 2 for thesis-quality runs is a design decision, not a code question. Left open for author decision.
2. ~~**`--pgfplots` output destination**~~: **Resolved.** Output goes to stdout via `printf()`. The user pipes to a file (`> results.dat`). One timing table per gate followed by one memory table, all in a single stdout stream.
3. ~~**CSV header format**~~: **Resolved.** Header is `qubits,dim,gate,method,time_ms,ops_per_sec,memory_bytes,iterations`. Row format: `%u,%ld,%s,%s,%.6f,%.0f,%zu,%d`. Stable.
4. ~~**`bench_mixed.c` enum names**~~: **Resolved.** Enum at line 721 of `bench_mixed.c`: `M_QLIB=0, M_QLIB_TILED=1, M_BLAS=2, M_NAIVE=3, M_QUEST=4, M_QULACS=5`. The first value is `M_QLIB`, not `M_QLIB_PACKED`. The method *string* reported in output is `"qlib_packed"`.
5. ~~**`blas_dense` upper qubit cutoff**~~: **Resolved.** Hardcoded as `qubits <= 8` in `bench_mixed.c` (`run_dense` flag). Not derived from `--max-qubits`, not configurable. There is no separate size guard; `qlib_packed` and `qlib_tiled` run unrestricted up to `--max-qubits`.
6. ~~**QuEST version**~~: **Resolved.** CMakeLists.txt line 17 comment explicitly states "QuEST v4 requires the root directory as include path." v4 confirmed.
7. ~~**`bench` target dependency list**~~: **Resolved.** The `bench` custom target depends on exactly `test_state` and `test_gate` (CMakeLists.txt line 163). This is a hardcoded list; new test targets require manual update.
8. **Qulacs `--step` interaction with OpenMP threshold**: Qulacs enables OpenMP at dim ≥ 1024 (10 qubits). With `--step 1`, intermediate qubit counts (9, 11) are also benchmarked; their behavior relative to the OpenMP threshold is undocumented. Left open.
