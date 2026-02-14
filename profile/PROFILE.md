# Profiling Framework

Micro-profiling framework for measuring individual gate function performance.

## Directory Layout

```
profile/
  CMakeLists.txt          # Build definitions; one target per profiled gate
  include/profile.h       # Header-only framework (timing, stats, runner, output)
  src/profile_x_tiled.c   # Profiling harness for x_tiled() (exemplar)
  results_x_tiled.csv     # Raw CSV output from a profiling run
  REPORT_x_tiled.md       # Analysis report written from those results
```

Included from the top-level `CMakeLists.txt` via `add_subdirectory(profile)`.
Binaries go to `<build>/profile/`.

## The Framework (`profile.h`)

A single header providing everything needed to profile a gate function:

| Component | What it does |
|-----------|-------------|
| `prof_time_ns()` | `clock_gettime(CLOCK_MONOTONIC)` in nanoseconds |
| `prof_stats_t` / `prof_compute_stats()` | Sorts samples, computes mean, median, min, max, stddev, p5, p95 |
| `prof_tiled_bytes()` / `prof_packed_bytes()` | Estimates total bytes read+written per gate call (for bandwidth calculation) |
| `prof_result_t` | Collects function name, qubit count, target, code path, stats, and derived metrics (BW in GB/s, elements/ns) |
| `prof_print_header()` / `prof_print_result()` | Human-readable table output |
| `prof_print_csv_header()` / `prof_print_csv_row()` | CSV output |
| `prof_run_1q()` | Generic profiling loop for any `void fn(state_t*, qubit_t)` gate: allocates state to \|+\>^n, runs warmup iterations, measures, computes stats |

### Defaults

| Parameter | Value |
|-----------|-------|
| Warmup iterations | 200 |
| Measurement iterations | 2000 |
| Min qubits | 2 |
| Max qubits | 12 |

## Writing a New Profiling Harness

Each gate gets its own `src/profile_<name>.c`. Follow the existing `profile_x_tiled.c` pattern:

1. Declare the backend function(s) under test via `extern` (they are not in the public `gate.h`).
2. Parse CLI options (the `options_t` struct and `getopt_long` boilerplate can be copied).
3. Loop over qubit counts and target qubits, calling `prof_run_1q()` for each.
4. Print results with `prof_print_result()`; optionally write CSV with `prof_print_csv_row()`.
5. Register the executable in `CMakeLists.txt`:
   ```cmake
   add_executable(profile_<name> "${PROF_SRC_DIR}/profile_<name>.c")
   set_target_properties(profile_<name> PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/profile")
   target_include_directories(profile_<name> PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
   target_link_libraries(profile_<name> PRIVATE q)
   ```

### CLI Options (convention from `profile_x_tiled`)

```
--min-qubits N    --max-qubits N    --warmup N    --iters N
--csv FILE        --target N        --compare     -v/--verbose
```

## Running

```bash
cd cmake-build-debug && cmake --build .
./profile/profile_x_tiled                              # defaults
./profile/profile_x_tiled --csv results.csv            # with CSV output
./profile/profile_x_tiled --min-qubits 10 --iters 5000 # custom range
./profile/profile_x_tiled --compare                    # tiled vs packed
```

## Output Artifacts

- **Console**: aligned table with median/mean/min/p5/p95/stddev/BW/elem-per-ns.
- **CSV** (`--csv`): one row per (function, qubits, target) combination. Columns: `function, qubits, target, path, median_ns, mean_ns, min_ns, p5_ns, p95_ns, stddev_ns, max_ns, data_bytes, touch_bytes, bw_gbs, elem_per_ns, samples`.
- **Report**: manually written `REPORT_<name>.md` interpreting the results.

## Key Design Decisions

- **Release builds** (`-O3 -march=native`) for production-representative timings.
- **State initialized to \|+\>^n** so all matrix elements are nonzero (no short-circuit branches).
- **Median** used for derived metrics (BW, elem/ns) because it is robust to outliers.
- **Within-tile vs cross-tile** distinction tracked via `path` field (target < `LOG_TILE_DIM` is within-tile).
- **`--compare`** flag enables side-by-side profiling of tiled vs packed implementations of the same gate.
