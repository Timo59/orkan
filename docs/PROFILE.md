# Profiling Framework

Micro-profiling framework for measuring individual gate function performance. Each gate gets its own harness; the framework provides timing, statistics, bandwidth estimation, and CSV/console output.

---

## Directory Layout

```
profile/
  CMakeLists.txt                  Build definitions; one target per profiled gate
  include/profile.h               Header-only framework (timing, stats, runner, output)
  src/profile_x_tiled.c           Profiling harness for x_tiled() (exemplar)
  src/profile_h_tiled.c           Profiling harness for h_tiled()
```

---

## Public API (header-only, `profile/include/profile.h`)

| Component | Purpose |
|---|---|
| `prof_time_ns()` | `clock_gettime(CLOCK_MONOTONIC)` in nanoseconds |
| `prof_stats_t` / `prof_compute_stats()` | Sorts samples in-place; computes mean, median, min, max, stddev, p5, p95 |
| `prof_tiled_bytes()` / `prof_packed_bytes()` | Estimates total bytes read+written per gate call (for bandwidth) |
| `gate_1q_fn` | `typedef void (*gate_1q_fn)(state_t *, qubit_t)` — accepted by `prof_run_1q()` |
| `prof_result_t` | Holds function name, qubit count, target, code path, stats, data\_bytes, touch\_bytes, derived metrics (BW in GB/s, elements/ns) |
| `prof_print_header()` / `prof_print_result()` | Aligned-column console output |
| `prof_print_csv_header()` / `prof_print_csv_row()` | CSV output |
| `prof_run_1q()` | Generic loop: initialises state to \|+>^n, runs warmup, measures, computes stats, frees state |

### `prof_run_1q()` signature

```c
prof_result_t prof_run_1q(
    gate_1q_fn   fn,
    const char  *func_name,
    state_type_t stype,       // MIXED_TILED or MIXED_PACKED
    qubit_t      qubits,
    qubit_t      target,
    int          warmup,
    int          iters);
```

`stype` controls storage layout. When `MIXED_TILED`, the `path` field is set to `"within-tile"` or `"cross-tile"` based on whether `target < LOG_TILE_DIM`. When `MIXED_PACKED`, `path` is `"packed"`.

### Defaults

| Parameter | Constant | Default |
|---|---|---|
| Warmup iterations | `PROF_DEFAULT_WARMUP` | 200 |
| Measurement iterations | `PROF_DEFAULT_ITERS` | 2000 |
| Min qubits (harness default) | — | `LOG_TILE_DIM` (3 in Debug, 5 in Release); floored to 2 |
| Max qubits | `PROF_DEFAULT_MAX_QUBITS` | 12 |

---

## Architecture

### Key Design Decisions

- **State initialised to \|+>^n** so all matrix elements are nonzero (no short-circuit branches).
- **Median** for derived metrics (BW, elem/ns) — robust to outliers.
- **Within-tile vs cross-tile** distinction tracked via `path` field. `LOG_TILE_DIM` is compile-time, so Debug (3) and Release (5) classify the same target qubit differently. No tagging convention distinguishes profiles collected under the two builds — profiles from different build types are not directly comparable for within-tile/cross-tile breakdown.
- **`--compare` is a per-harness opt-in**, not a framework feature. Each harness hardcodes a second `prof_run_1q()` call with the packed variant (`x_packed`, `h_packed`) and `MIXED_PACKED`. Only meaningful for gates with both tiled and packed implementations.
- **Backend functions accessed via `extern`** rather than the public API, so harnesses can measure internal functions not exported in `gate.h`.

### Output Phases (per harness)

1. **Main sweep**: all (qubits, target) combinations, printed as an aligned table and optionally written to CSV.
2. **Scaling analysis**: a second table fixing `target=0` (within-tile) and `target=LOG_TILE_DIM` (cross-tile) while varying qubit count, isolating scaling behaviour from target variation. Always printed; CSV rows written if `--csv` was given.

---

## Build Integration

Included from the top-level `CMakeLists.txt` via `add_subdirectory(profile)`. Binaries go to `<build>/profile/`.

| Target | Source |
|---|---|
| `profile_x_tiled` | `src/profile_x_tiled.c` |
| `profile_h_tiled` | `src/profile_h_tiled.c` |

Both targets follow the same pattern:

```cmake
add_executable(profile_<name> "${PROF_SRC_DIR}/profile_<name>.c")
set_target_properties(profile_<name> PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/profile")
target_include_directories(profile_<name> PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(profile_<name> PRIVATE ${PROJECT_NAME_LOWER})
```

Profiling is typically run against Release builds (`-O3 -march=native`) for production-representative timings; Debug builds are also supported to exercise small-tile (8x8) code paths.

```bash
cmake --preset release && cmake --build --preset release
./cmake-build-release/profile/profile_x_tiled
./cmake-build-release/profile/profile_h_tiled
```

### CLI Options (convention from existing harnesses)

```
--min-qubits N    --max-qubits N    --warmup N    --iters N
--csv FILE        --target N        --compare     -v/--verbose
```

CSV row schema: `function, qubits, target, path, median_ns, mean_ns, min_ns, p5_ns, p95_ns, stddev_ns, max_ns, data_bytes, touch_bytes, bw_gbs, elem_per_ns, samples`.

---

## Writing a New Harness

1. Declare the backend function(s) under test via `extern`.
2. Copy the `options_t` struct and `getopt_long` boilerplate from `profile_x_tiled.c`.
3. Call `prof_run_1q(fn, name, MIXED_TILED, n, t, warmup, iters)` for each (qubit count, target). Pass `MIXED_PACKED` and the packed variant for `--compare`.
4. Print results with `prof_print_result()`; optionally `prof_print_csv_row()`.
5. Add a scaling-analysis section (fixed `target=0` and `target=LOG_TILE_DIM`, varying qubits).
6. Register the executable in `CMakeLists.txt` using the pattern above.
