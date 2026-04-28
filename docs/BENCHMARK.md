# Benchmark Module — Technical Specification

## 1. Purpose

Benchmark single-, two-, and three-qubit gate execution on density matrices across
six backends, producing per-gate-application timing statistics and effective memory
bandwidth.  The benchmark sweeps qubit counts and (optionally) qubit positions,
outputting results in console, CSV, or PGFplots format.

---

## 2. Backends

| # | Backend        | Language | Storage Format                        | Notes                        |
|---|----------------|----------|---------------------------------------|------------------------------|
| 1 | orkan\_tiled    | C        | Tiled Hermitian (lower-triangle tiles)| Always available             |
| 2 | orkan\_packed   | C        | Packed Hermitian (column-major lower) | Always available             |
| 3 | BLAS baseline  | C        | Full 2^n x 2^n dense matrix           | Max 10 qubits; `nan` above  |
| 4 | QuEST          | C        | Flat vector of 4^n `complex<double>`  | Built if found by CMake      |
| 5 | Qulacs         | C++      | Flat vector of 4^n `complex<double>`  | Built if found by CMake      |
| 6 | Qiskit Aer     | C++      | Flat vector of 4^n `complex<double>`  | Built if found by CMake      |

CMake determines at build time which competitors are available.  Unavailable backends
print `nan` in all output columns.  BLAS prints `nan` for qubit counts above 8.

---

## 3. Supported Gates

| Symbol(s)         | Qubits | Parametrised | Parameter     |
|-------------------|--------|--------------|---------------|
| `X`               | 1      | No           | —             |
| `Y`               | 1      | No           | —             |
| `Z`               | 1      | No           | —             |
| `H`               | 1      | No           | —             |
| `S`               | 1      | No           | —             |
| `T`               | 1      | No           | —             |
| `RX`              | 1      | Yes          | angle (rad)   |
| `RY`              | 1      | Yes          | angle (rad)   |
| `RZ`              | 1      | Yes          | angle (rad)   |
| `CX` / `CNOT`     | 2      | No           | —             |
| `CY`              | 2      | No           | —             |
| `CZ`              | 2      | No           | —             |
| `CCX` / `TOFFOLI` | 3      | No           | —             |

The gate parser is **case-insensitive**: `x`, `X`, `cnot`, `CNOT`, `Toffoli` all resolve
correctly.  Aliases (`CX`/`CNOT`, `CCX`/`TOFFOLI`) map to the same benchmark.

---

## 4. Command-Line Interface

```
bench_gate <GATE> [OPTIONS]
```

### 4.1 Positional Argument

| Argument | Required | Description                                |
|----------|----------|--------------------------------------------|
| `GATE`   | Yes      | Gate symbol from the table in Section 3    |

### 4.2 Options

| Flag             | Default   | Range              | Description                                             |
|------------------|-----------|--------------------|---------------------------------------------------------|
| `--par <float>`  | —         | —                  | Gate parameter (radians). **Mandatory** for RX/RY/RZ.   |
| `--min-qubits`   | `2`       | 0–30, ≥ gate qubits | Minimum number of qubits in sweep                     |
| `--max-qubits`   | `8`       | 0–30, ≥ min-qubits | Maximum number of qubits in sweep                      |
| `--samples`      | `10`      | 1–1024             | Number of timed blocks (each yields one sample)          |
| `--iterations`   | `100`     | ≥ 1                | Number of repetitions within each timed block            |
| `--warm-up`      | `3`       | ≥ 0                | Number of warm-up iterations (discarded)                 |
| `--per-qubit`    | off       | —                  | Report per-qubit-position statistics                     |
| `--output`       | `console` | —                  | Output format: `console`, `csv`, `pgfplots`              |

#### Iteration Tuning Guide

The `--iterations` parameter exists to push the timed block above clock resolution
noise (`CLOCK_MONOTONIC` has 1 us resolution on macOS).  At high qubit counts a
single gate application already takes milliseconds, making large iteration counts
wasteful.  The table below shows recommended minimum iterations per qubit count,
targeting a timed-block duration of at least 1 ms (1000x clock resolution):

| Qubits | Iterations | Rationale |
|--------|-----------|-----------|
| 2-4 | 5000 | Fastest backends need many reps to exceed clock noise |
| 5 | 1000 | |
| 6 | 200 | |
| 7 | 50 | |
| 8 | 10 | |
| 9 | 4 | |
| 10 | 2 | |
| 11+ | 1 | Single gate application already takes 4+ ms |

Using this schedule instead of a flat `--iterations 100` reduces the total wall-clock
time for a full 2-14 qubit H-gate sweep from ~6-22 hours to ~6 minutes while
improving measurement accuracy at low qubit counts.

### 4.3 Exit Codes

| Code | Condition                                        |
|------|--------------------------------------------------|
| 0    | Success                                          |
| 1    | Unrecognised gate symbol                         |
| 1    | Parametrised gate invoked without `--par`        |
| 1    | Unknown flag or invalid argument                 |

Unavailable backends (not built) are **not** errors — they silently produce `nan`.

---

## 5. Measurement Methodology

### 5.1 State Initialisation

Each backend receives a randomly initialised **Hermitian** matrix of dimension
2^n x 2^n in its respective memory layout.  Hermiticity is the only constraint — it
is required by orkan's symmetric storage formats.  No positive semi-definiteness or
trace constraint is imposed; all competitors (QuEST, Qulacs, Qiskit Aer) accept
arbitrary complex matrices and perform pure `ρ → UρU†` linear algebra without
validating the input.

Initialisation is **not** timed.  The state is **not** re-initialised between
iterations; gates are applied to an evolving state throughout warm-up and timed runs.

### 5.2 Timing

**Clock**: `clock_gettime(CLOCK_MONOTONIC)` — nanosecond resolution, immune to NTP
adjustments, ~20 ns overhead per call.

### 5.3 Aggregate Mode (default)

One **sweep** = apply the gate sequentially to every valid qubit position (or qubit
combination for multi-qubit gates).  For a 1-qubit gate on n qubits, one sweep
applies the gate n times.  For a 2-qubit gate, one sweep covers all n(n-1)
permutations (both (0,1) and (1,0), since control vs target matters).  For a 3-qubit
gate, all n(n-1)(n-2) ordered triples.

Timing uses a **two-level structure** to reduce clock overhead noise.  The inner
level (`iterations`) repeats the sweep many times within a single timed block.  The
outer level (`samples`) produces independent measurements for statistics.

```
for i in 1..warm_up:
    execute sweep                          # discarded

P = positions_in_sweep
for s in 1..samples:
    t0 = clock_gettime(MONOTONIC)
    for i in 1..iterations:
        execute sweep
    t1 = clock_gettime(MONOTONIC)
    sample[s] = (t1 - t0) / (iterations * P)
```

Each sample is the **mean time per gate application** within that block.  Statistics
(mean, median, min, 95% CI) are computed over the `samples` values.

### 5.4 Per-Qubit Mode (`--per-qubit`)

Each qubit position (or combination) is timed independently, using the same
two-level structure:

```
for each position p:
    for i in 1..warm_up:
        apply gate at position p               # discarded

    for s in 1..samples:
        t0 = clock_gettime(MONOTONIC)
        for i in 1..iterations:
            apply gate at position p
        t1 = clock_gettime(MONOTONIC)
        sample[p][s] = (t1 - t0) / iterations
```

Each sample is the **mean time per gate application** at position p within that
block.  Statistics are computed per position over the `samples` values.

### 5.5 Qubit Positions

| Gate qubits | Positions per n-qubit system | Notation           |
|-------------|------------------------------|--------------------|
| 1           | n                            | target t           |
| 2           | n(n-1)                       | (control, target)  |
| 3           | n(n-1)(n-2)                  | (c1, c2, target)   |

---

## 6. Reported Quantities

### 6.1 Timing Statistics

Computed over the `samples` measurements for each (backend, qubit-count) pair
(aggregate mode) or each (backend, qubit-count, position) triple (per-qubit mode).

| Quantity   | Definition                                                |
|------------|-----------------------------------------------------------|
| Mean       | Arithmetic mean of samples (microseconds)                 |
| Median     | 50th percentile of sorted samples (microseconds)          |
| Min        | Minimum sample (microseconds)                             |
| 95% CI lo  | 2.5th percentile of sorted samples (microseconds)         |
| 95% CI hi  | 97.5th percentile of sorted samples (microseconds)        |

The 95% confidence interval is **percentile-based** (non-parametric), not parametric.
All times are reported in **microseconds (us)**.

### 6.2 Effective Bandwidth

Bandwidth measures how efficiently each backend utilises the memory hierarchy.

```
bw_gbs = touch_bytes / (median_ns)
```

where `median_ns` is the median sample converted to nanoseconds, and `touch_bytes`
depends on the backend's storage format:

| Backend        | touch\_bytes formula                                     |
|----------------|----------------------------------------------------------|
| orkan\_tiled    | `2 * n_tile_pairs * TILE_SIZE * sizeof(cplx_t)` (see below) |
| orkan\_packed   | `2 * dim * (dim + 1) / 2 * sizeof(cplx_t)`               |
| BLAS           | `2 * 4^n * sizeof(cplx_t)`                               |
| QuEST          | `2 * 4^n * sizeof(cplx_t)`                               |
| Qulacs         | `2 * 4^n * sizeof(cplx_t)`                               |
| Qiskit Aer     | `2 * 4^n * sizeof(cplx_t)`                               |

where `sizeof(cplx_t) = 16` bytes, `dim = 2^n`, `TILE_SIZE = TILE_DIM^2`,
`n_tile_pairs = n_tiles * (n_tiles + 1) / 2`, and the factor of 2 accounts for
one read and one write per element.

**orkan\_tiled small-state special case:** When `qubits < LOG_TILE_DIM`, the entire
density matrix fits inside a single diagonal tile.  The tile storage is padded to
`TILE_DIM x TILE_DIM`, but the kernel only reads and writes the `dim x dim` upper-left
submatrix (the rest is zero padding).  In this case `touch_bytes = 2 * dim * dim *
sizeof(cplx_t)` instead of the tile-pair formula.

The orkan backends report **lower** `touch_bytes` than competitors because they
exploit Hermitian symmetry.  This is intentional: timing columns provide the direct
performance comparison; bandwidth shows memory-hierarchy utilisation efficiency.

---

## 7. BLAS Baseline

The BLAS baseline computes density-matrix gate conjugation via standard matrix
multiplication:

```
rho' = G * rho * G^dagger
```

implemented as two `zgemm` calls.

**Pre-computation** (not timed): For every qubit position in the sweep, the full
2^n x 2^n gate matrix is pre-computed by Kronecker products:

```
G_full = I_{2^a} (x) G (x) I_{2^b}
```

All such matrices (and their conjugate transposes) are generated before the timing
loop begins.

**Cap**: BLAS is benchmarked only for qubit counts <= 10.  Above 10, it is skipped
and `nan` is printed for all its columns.

---

## 8. Output Formats

All output is written to **stdout**.  The `--output` flag controls the format, not
the destination.  To save results, redirect: `bench_gate X --output csv > results.csv`.

### 8.1 Console — Aggregate Mode

One table per qubit count:

```
Gate: X   Qubits: 4   Samples: 10   Iterations: 100   Warm-up: 3

Backend          Mean (us)   Median (us)   Min (us)   95% CI (us)       BW (GB/s)
----------------------------------------------------------------------------------
orkan_tiled          12.3        11.8        10.5     [10.2, 14.4]         82.3
orkan_packed         10.1         9.7         8.9     [ 8.5, 11.7]         91.1
BLAS                45.2        44.8        43.1     [42.0, 48.4]         36.2
QuEST               18.7        18.2        16.9     [15.8, 21.6]         55.4
Qulacs              15.4        15.1        14.2     [13.5, 17.3]         66.8
Qiskit Aer          22.1        21.5        20.3     [19.1, 25.1]         46.9
```

### 8.2 Console — Per-Qubit Mode

One table per qubit count.  Rows = backend x qubit position:

```
Gate: X   Qubits: 4   Samples: 10   Iterations: 100   Warm-up: 3

Backend        Target   Mean (us)   Median (us)   Min (us)   95% CI (us)       BW (GB/s)
-----------------------------------------------------------------------------------------
orkan_tiled        0       11.2        10.9        10.1     [ 9.8, 12.6]         84.1
orkan_tiled        1       12.5        12.1        11.0     [10.5, 14.5]         75.8
...
Qiskit Aer        3       23.4        22.8        21.1     [20.0, 26.8]         44.2
```

For 2-qubit gates, the columns are `Control  Target`.
For 3-qubit gates: `C1  C2  Target`.

### 8.3 CSV

Header row followed by data.  Columns:

**Aggregate mode:**
```
qubits,backend,mean_us,median_us,min_us,ci_lo_us,ci_hi_us,bw_gbs
```

**Per-qubit mode (1-qubit gate):**
```
qubits,backend,target,mean_us,median_us,min_us,ci_lo_us,ci_hi_us,bw_gbs
```

**Per-qubit mode (2-qubit gate):**
```
qubits,backend,control,target,mean_us,median_us,min_us,ci_lo_us,ci_hi_us,bw_gbs
```

**Per-qubit mode (3-qubit gate):**
```
qubits,backend,c1,c2,target,mean_us,median_us,min_us,ci_lo_us,ci_hi_us,bw_gbs
```

Unavailable backends and BLAS above 8 qubits emit `nan` for all numeric fields.

### 8.4 PGFplots

**Aggregate mode:** One table per statistical quantity.  Number of qubits in rows,
backends in columns.

```
% Gate: X   Samples: 10   Iterations: 100   Warm-up: 3
% Quantity: mean_us
qubits  orkan_tiled  orkan_packed  BLAS    QuEST   Qulacs  Aer
2       1.2         1.0          3.5     2.1     1.8     2.5
3       2.8         2.4          8.1     4.9     4.0     5.8
...
```

Tables emitted: `mean_us`, `median_us`, `min_us`, `ci_lo_us`, `ci_hi_us`, `bw_gbs`.

**Per-qubit mode:** One table per (qubit count, statistical quantity).  Qubit
position(s) in rows, backends in columns.

For a 1-qubit gate:
```
% Gate: X   Qubits: 4   Samples: 10   Iterations: 100   Warm-up: 3
% Quantity: mean_us
target  orkan_tiled  orkan_packed  BLAS    QuEST   Qulacs  Aer
0       11.2        9.8          40.1    17.3    14.1    20.5
1       12.5        10.3         41.8    18.1    15.0    21.8
2       11.8        9.9          39.5    17.0    14.3    20.1
3       13.1        10.8         42.3    18.9    15.5    22.4
```

For 2-qubit gates, rows use `control  target`.
For 3-qubit gates, rows use `c1  c2  target`.

---

## 9. Kraus Channel Benchmark (`bench_kraus`)

A separate executable benchmarking quantum channel (Kraus map) application on
density matrices.  Measures the full operation `ρ → Σ_k A_k ρ A_k†` for
predefined noise channels across the same competitor backends.

### 9.1 Channels

| Channel       | Qubits | Operators | Parameter | Kraus Operators |
|---------------|--------|-----------|-----------|-----------------|
| `Depolarize`  | 1      | 4         | p (prob)  | `{√(1-p) I, √(p/3) X, √(p/3) Y, √(p/3) Z}` |

The channel parser is case-insensitive.  Default parameter: `p = 0.1`.

### 9.2 Backends

All six gate-benchmark backends participate.

| # | Backend        | Strategy | Cap | Notes |
|---|----------------|----------|-----|-------|
| 1 | orkan\_tiled    | Pre-built superoperator via `kraus_to_superop`, applied by `channel_1q` on tiled Hermitian storage | None | Superop built once at init |
| 2 | orkan\_packed   | Pre-built superoperator via `kraus_to_superop`, applied by `channel_1q` on packed Hermitian storage | None | Superop built once at init |
| 3 | BLAS baseline  | Kronecker-expand each operator to full N×N, then 2K zgemm + accumulate | 10 qubits | Pre-computes Kronecker products per position |
| 4 | QuEST          | Pre-built superoperator via `createKrausMap`/`setKrausMap`, applied by `mixKrausMap` | None | Superop built once at init |
| 5 | Qulacs         | Manual iterate-and-sum: copy ρ, apply `dm_single_qubit_dense_matrix_gate`, accumulate via `dm_state_add` | 10 qubits | No native Kraus API at csim level |
| 6 | Qiskit Aer     | Pre-built superoperator via `kraus_superop` + `vectorize_matrix`, applied by `apply_superop_matrix` | None | Superop built once at init |

**Fairness**: Each backend does what it would do per-application in production.
orkan, QuEST and Aer pre-build superoperators at init (not timed); BLAS pre-builds
Kronecker-expanded matrices at init (not timed); Qulacs copies small operator
matrices at init.  Only the `apply` call is timed.

### 9.3 Command-Line Interface

```
bench_kraus [OPTIONS]
```

| Flag             | Default       | Description                                             |
|------------------|---------------|---------------------------------------------------------|
| `--channel`      | `depolarize`  | Channel name from the channel table                     |
| `--par <float>`  | `0.1`         | Channel parameter                                       |
| `--min-qubits`   | `2`           | Minimum system size                                     |
| `--max-qubits`   | `10`          | Maximum system size                                     |
| `--samples`      | `10`          | Timed blocks                                            |
| `--iterations`   | auto-tuned    | Repetitions per block                                   |
| `--warm-up`      | `3`           | Warm-up iterations                                      |
| `--per-qubit`    | off           | Report per-position statistics                          |
| `--output`       | `console`     | Output format: `console`, `csv`, `pgfplots`             |

### 9.4 Measurement Methodology

Identical to the gate benchmark (Section 5): same two-level timing, same position
sweep, same statistics computation.  Position generation uses the channel's target
qubit count (1 for depolarizing → n positions).

### 9.5 Kraus-Specific Files

| File | Responsibility | Language |
|------|----------------|----------|
| `benchmark/include/bench_kraus.h` | Channel table, Kraus backend vtable, run declarations | C |
| `benchmark/src/bench_kraus_main.c` | CLI parsing, backend table, orchestration | C |
| `benchmark/src/bench_kraus_chan.c` | Analytic Kraus operator construction | C |
| `benchmark/src/bench_kraus_run.c` | Timing loop + position generation | C |
| `benchmark/src/bench_kraus_orkan.c` | orkan\_tiled and orkan\_packed Kraus backends | C |
| `benchmark/src/bench_kraus_blas.c` | BLAS Kraus backend | C |
| `benchmark/src/bench_kraus_quest.c` | QuEST Kraus backend | C |
| `benchmark/src/bench_kraus_qulacs.cpp` | Qulacs Kraus backend | C++ |
| `benchmark/src/bench_kraus_aer.cpp` | Qiskit Aer Kraus backend | C++ |

The Kraus benchmark reuses the gate benchmark's formatter source files
(`bench_console.c`, `bench_csv.c`, `bench_pgfplots.c`), result computation
(`bench_result.c`), and utilities (`bench_util.c`) by compiling them into the
`bench_kraus` target.  No gate benchmark source files are modified.

### 9.6 Build Integration

The `bench_kraus` target is defined in `benchmark/CMakeLists.txt` alongside
`bench_gate`.  It is `EXCLUDE_FROM_ALL` and must be built explicitly:

```
cmake --build --preset release --target bench_kraus
```

Same conditional compilation pattern as `bench_gate`: `HAVE_QUEST`, `HAVE_QULACS`,
`HAVE_AER` guard the competitor backends.

---

## 10. File Organisation

No source file exceeds 500 lines of well-documented code.

### Gate Benchmark (`bench_gate`)

| File                           | Responsibility                              | Language |
|--------------------------------|---------------------------------------------|----------|
| `benchmark/src/bench_main.c`   | CLI parsing, gate dispatch, orchestration   | C        |
| `benchmark/src/bench_run.c`    | Timing loop (aggregate & per-qubit modes)   | C        |
| `benchmark/src/bench_result.c` | Mean, median, min, percentile CI, bandwidth | C        |
| `benchmark/src/bench_util.c`   | PRNG, Hermitian init helpers, huge-page alloc | C        |
| `benchmark/src/bench_orkan.c`   | orkan\_tiled and orkan\_packed wrappers       | C        |
| `benchmark/src/bench_blas.c`   | BLAS baseline (Kronecker + zgemm)           | C        |
| `benchmark/src/bench_quest.c`  | QuEST density-matrix wrapper                | C        |
| `benchmark/src/bench_qulacs.cpp` | Qulacs density-matrix wrapper             | C++      |
| `benchmark/src/bench_aer.cpp`  | Qiskit Aer density-matrix wrapper           | C++      |
| `benchmark/src/bench_console.c`| Console output formatter                    | C        |
| `benchmark/src/bench_csv.c`    | CSV output formatter                        | C        |
| `benchmark/src/bench_pgfplots.c`| PGFplots output formatter                  | C        |
| `benchmark/include/bench.h`    | Shared types, constants, function decls     | C        |
| `benchmark/include/bench_result.h` | Result type and statistics API           | C        |

### Kraus Benchmark (`bench_kraus`)

| File                                  | Responsibility                        | Language |
|---------------------------------------|---------------------------------------|----------|
| `benchmark/include/bench_kraus.h`     | Channel table, Kraus vtable, run decls | C       |
| `benchmark/src/bench_kraus_main.c`    | CLI parsing, orchestration            | C        |
| `benchmark/src/bench_kraus_chan.c`     | Analytic Kraus operator construction  | C        |
| `benchmark/src/bench_kraus_run.c`     | Timing loop + position generation     | C        |
| `benchmark/src/bench_kraus_orkan.c`    | orkan\_tiled and orkan\_packed backends | C        |
| `benchmark/src/bench_kraus_blas.c`    | BLAS Kraus backend                    | C        |
| `benchmark/src/bench_kraus_quest.c`   | QuEST Kraus backend                   | C        |
| `benchmark/src/bench_kraus_qulacs.cpp`| Qulacs Kraus backend                  | C++      |
| `benchmark/src/bench_kraus_aer.cpp`   | Qiskit Aer Kraus backend              | C++      |

Shared files compiled into both targets: `bench_result.c`, `bench_util.c`,
`bench_console.c`, `bench_csv.c`, `bench_pgfplots.c`.

### 10.1 Build Integration

The root `CMakeLists.txt` defines `option(BUILD_BENCHMARKS … OFF)`.  When enabled,
it adds the `benchmark/` subdirectory.  The benchmark's own `CMakeLists.txt` requires
a Release build (`CMAKE_BUILD_TYPE STREQUAL "Release"`) and returns early otherwise.
Both `bench_gate` and `bench_kraus` targets are marked `EXCLUDE_FROM_ALL` — they must
be built explicitly (e.g. `cmake --build --preset release --target bench_gate`).

Backend wrappers are guarded by `#ifdef HAVE_QUEST` / `HAVE_QULACS` / `HAVE_AER`.
When a backend is not available, its wrapper returns `nan` for all measurements.

#### 10.1.1 QuEST

Detected via `find_path` (for `quest/include/quest.h`) and `find_library`.  Searched
in `extern/QuEST`, `$QuEST_DIR`, and standard system paths.  Sets `HAVE_QUEST` and
links the found library directly.

#### 10.1.2 Qulacs

Built from source as a `csim_static` static library from the submodule at
`extern/qulacs/src/csim/`.  Requires Eigen3 headers.  Sets `HAVE_QULACS`.  Links
OpenMP when available.

#### 10.1.3 Qiskit Aer

Header-only density-matrix kernel vendored to `extern/aer-dm/`.  Requires
`nlohmann/json.hpp` (searched via Homebrew and standard paths).  Sets `HAVE_AER`.
Links OpenMP when available.

---

## 11. Fairness Invariants

1. **Same initial state**: All backends receive the same randomly initialised density
   matrix for each qubit count (seeded identically where possible).
2. **Same gate sequence**: All backends apply the same gate to the same qubit
   positions in the same order.
3. **No re-initialisation**: The state evolves across iterations — no backend gets a
   fresh matrix between timed runs.
4. **Pre-computation excluded**: BLAS Kronecker products are computed outside the
   timing loop.  Library initialisation costs are absorbed by warm-up iterations.
5. **Warm-up**: All backends execute `--warm-up` iterations before timing begins,
   ensuring caches are warm and lazy initialisations complete.
6. **Uniform timing**: The same `clock_gettime(CLOCK_MONOTONIC)` call brackets every
   backend's execution — no backend-specific timing.
