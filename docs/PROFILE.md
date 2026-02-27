# Profiling Framework

Micro-profiling framework for measuring individual gate function performance.

## Status

In Progress. Two harnesses exist (`x_tiled`, `h_tiled`). The framework API is stable. Additional gate harnesses and a `prof_run_2q()` runner (for two-qubit gates) are not yet implemented.

## Directory Layout

```
profile/
  CMakeLists.txt                  # Build definitions; one target per profiled gate
  include/profile.h               # Header-only framework (timing, stats, runner, output)
  src/profile_x_tiled.c           # Profiling harness for x_tiled() (exemplar)
  src/profile_h_tiled.c           # Profiling harness for h_tiled()
  results_x_tiled.csv             # Raw CSV output from a profile_x_tiled run
  results_h_tiled.csv             # Raw CSV output from a profile_h_tiled run
  REPORT_x_tiled.md               # Analysis report interpreting x_tiled results
  PLAN_h_tiled_optimization.md    # Planning notes for h_tiled optimization work
```

## Dependencies (Inputs from Other Modules)

- **`q` library** (`src/`, `include/`): all harnesses link against `libq` via `target_link_libraries(... PRIVATE q)`, inheriting compiler flags, BLAS, OpenMP, and public include dirs transitively.
- **Gate internals**: each harness declares backend functions (e.g. `x_tiled`, `h_tiled`) via `extern`; these are not exposed in the public `gate.h`.
- **`state_t` / `qubit_t`**: types from `include/state.h` used by `prof_run_1q()`.

## The Framework (`profile.h`)

A single header providing everything needed to profile a gate function:

| Component | What it does |
|-----------|-------------|
| `prof_time_ns()` | `clock_gettime(CLOCK_MONOTONIC)` in nanoseconds |
| `prof_stats_t` / `prof_compute_stats()` | Sorts samples in-place, computes mean, median, min, max, stddev, p5, p95 |
| `prof_tiled_bytes()` / `prof_packed_bytes()` | Estimates total bytes read+written per gate call (for bandwidth calculation) |
| `gate_1q_fn` | Typedef `void (*gate_1q_fn)(state_t *, qubit_t)` — function pointer type accepted by `prof_run_1q()` |
| `prof_result_t` | Collects function name, qubit count, target, code path, stats, data\_bytes, touch\_bytes, and derived metrics (BW in GB/s, elements/ns) |
| `prof_print_header()` / `prof_print_result()` | Human-readable table output |
| `prof_print_csv_header()` / `prof_print_csv_row()` | CSV output |
| `prof_run_1q()` | Generic profiling loop: initializes state to \|+\>^n on the stack, runs warmup iterations, measures, computes stats, frees state |

### `prof_run_1q()` Signature

```c
prof_result_t prof_run_1q(
    gate_1q_fn   fn,
    const char  *func_name,
    state_type_t stype,       /* MIXED_TILED or MIXED_PACKED */
    qubit_t      qubits,
    qubit_t      target,
    int          warmup,
    int          iters);
```

`stype` controls storage layout. When `MIXED_TILED`, the `path` field is set to `"within-tile"` or `"cross-tile"` based on whether `target < LOG_TILE_DIM`. When `MIXED_PACKED`, `path` is set to `"packed"`.

### Defaults

| Parameter | Framework constant | Default value |
|-----------|-------------------|---------------|
| Warmup iterations | `PROF_DEFAULT_WARMUP` | 200 |
| Measurement iterations | `PROF_DEFAULT_ITERS` | 2000 |
| Min qubits (harness default) | — | `LOG_TILE_DIM` (3 in Debug, 5 in Release); floored to 2 |
| Max qubits | `PROF_DEFAULT_MAX_QUBITS` | 12 |

Note: the effective default for `--min-qubits` in existing harnesses is `LOG_TILE_DIM`, not 2. The value 2 is only used when `LOG_TILE_DIM < 2`.

## Build Integration

Included from the top-level `CMakeLists.txt` via `add_subdirectory(profile)`. Binaries go to `<build>/profile/`.

Current registered targets:

| Target | Source |
|--------|--------|
| `profile_x_tiled` | `src/profile_x_tiled.c` |
| `profile_h_tiled` | `src/profile_h_tiled.c` |

Both targets follow the same pattern:
```cmake
add_executable(profile_<name> "${PROF_SRC_DIR}/profile_<name>.c")
set_target_properties(profile_<name> PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/profile")
target_include_directories(profile_<name> PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(profile_<name> PRIVATE q)
```

Profiling is typically run against Release builds (`-O3 -march=native`) for production-representative timings, but Debug builds are also supported to exercise small-tile (8×8) code paths.

## Running

```bash
# Build (from repo root)
cmake --preset debug && cmake --build --preset debug
# Or for production timings:
cmake --preset release && cmake --build --preset release

# Run with defaults
./profile/profile_x_tiled
./profile/profile_h_tiled

# Common options
./profile/profile_x_tiled --csv results.csv
./profile/profile_x_tiled --min-qubits 10 --iters 5000
./profile/profile_x_tiled --compare       # tiled vs packed side-by-side
```

### CLI Options (convention from existing harnesses)

```
--min-qubits N    --max-qubits N    --warmup N    --iters N
--csv FILE        --target N        --compare     -v/--verbose
```

## Output Artifacts

Each harness produces two output phases:

1. **Main sweep**: all (qubits, target) combinations, printed as an aligned table and optionally written to CSV.
2. **Scaling analysis**: a second table fixing target=0 (within-tile) and target=`LOG_TILE_DIM` (cross-tile) while varying qubit count, to isolate scaling behaviour from target variation. This section is always printed; CSV rows are written if `--csv` was specified.

After the sweep, `profile_x_tiled` prints suggested follow-up profiling commands (Apple Instruments, llvm-mca, compiler vectorization reports, `sample`). `profile_h_tiled` does not include this block.

- **Console**: aligned table with median/mean/min/p5/p95/stddev/BW/elem-per-ns.
- **CSV** (`--csv`): one row per (function, qubits, target) combination. Columns: `function, qubits, target, path, median_ns, mean_ns, min_ns, p5_ns, p95_ns, stddev_ns, max_ns, data_bytes, touch_bytes, bw_gbs, elem_per_ns, samples`.
- **Report**: manually written `REPORT_<name>.md` interpreting the results.
- **Plan**: `PLAN_<name>.md` for optimization planning (optional, not produced by the binary).

## Key Design Decisions

- **State initialized to \|+\>^n** so all matrix elements are nonzero (no short-circuit branches).
- **Median** used for derived metrics (BW, elem/ns) because it is robust to outliers.
- **Within-tile vs cross-tile** distinction tracked via `path` field (target < `LOG_TILE_DIM` is within-tile). This threshold is compile-time, so Debug (LOG\_TILE\_DIM=3) and Release (LOG\_TILE\_DIM=5) classify the same target qubit differently. No tagging convention exists to distinguish profiles collected under the two builds.
- **`--compare`** is a per-harness opt-in, not a framework feature. Each harness hardcodes a second `prof_run_1q()` call with the packed variant and `MIXED_PACKED`. It is only meaningful for gates that have both a tiled and a packed implementation.
- **Backend functions accessed via `extern`** rather than the public API, so the harness can measure internal functions not exported in `gate.h`.

## Writing a New Profiling Harness

Each gate gets its own `src/profile_<name>.c`. Follow the existing `profile_x_tiled.c` pattern:

1. Declare the backend function(s) under test via `extern`.
2. Copy the `options_t` struct and `getopt_long` boilerplate from an existing harness.
3. Call `prof_run_1q(fn, name, MIXED_TILED, n, t, warmup, iters)` for each (qubit count, target) pair. Pass `MIXED_PACKED` and the packed variant for `--compare` support.
4. Print results with `prof_print_result()`; optionally write CSV with `prof_print_csv_row()`.
5. Add a scaling analysis section (fixed target=0 and target=`LOG_TILE_DIM`, varying qubits).
6. Register the executable in `CMakeLists.txt` using the pattern shown in Build Integration above.

## TODO

- [ ] Add `prof_run_2q()` runner to `profile.h` for two-qubit gate harnesses.
- [ ] Add harnesses for remaining gate backends (e.g. `y_tiled`, `z_tiled`, `h_packed`, `cnot_tiled`).
- [ ] Decide whether `PLAN_<name>.md` files should live under `profile/` or `docs/`.
- [ ] Add a `REPORT_h_tiled.md` once h_tiled optimization work is complete.
- [ ] Add profiling tool suggestion block to `profile_h_tiled.c` (currently only in `profile_x_tiled.c`).

## Open Questions

1. **CSV file ownership**: `results_x_tiled.csv` and `results_h_tiled.csv` are committed to the repository. Is this intentional, or should result files be `.gitignore`d as generated artifacts?
2. **`prof_run_2q()` signature**: What should the function signature look like for two-qubit gates? The `qubit_t` pair (control, target) and whether it handles both tiled and packed variants is unspecified.
3. ~~**`--compare` semantics**~~: **Resolved.** `--compare` is a per-harness opt-in. Each harness calls `prof_run_1q()` a second time with the packed variant (`x_packed`, `h_packed`) and `MIXED_PACKED`. The flag is only valid for gates that have both implementations; the framework does not enforce this.
4. **`LOG_TILE_DIM` interaction**: The `path` classification threshold is compile-time. Debug (LOG\_TILE\_DIM=3) and Release (LOG\_TILE\_DIM=5) will classify the same target qubit differently, and no tagging or naming convention currently exists to distinguish profiles collected under the two builds. Profiles from different build types are not directly comparable for within-tile/cross-tile breakdown.
5. **`PLAN_h_tiled_optimization.md` format**: No convention is defined for what a `PLAN_<name>.md` file should contain or whether it is part of the profiling module's documented interface.
6. ~~**Max qubit ceiling**~~: **Resolved.** The default max of 12 (`PROF_DEFAULT_MAX_QUBITS`) is a soft default, overridable with `--max-qubits N`. No hard memory limit is enforced in the code.
