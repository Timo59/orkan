# Per-Qubit Benchmark Extension — Implementation Plan

Target binary: `bench_gate` (built from `benchmark/src/`)
New flag: `--per-qubit`

This plan is written for parallel execution. Steps with no shared dependencies can
be handed to independent agents or developers simultaneously. Each step is sized at
roughly 1–3 hours of focused implementation work.

---

## Dependency Graph

```
         ┌─────────────────────────────────────┐
         │           Step A: bench.h            │
         │  (types, constants, declarations)    │
         └──┬────────┬────────┬────────┬────────┘
            │        │        │        │
            ▼        ▼        ▼        ▼
         Step B   Step C   Step D   Step E   ← all depend only on A; run in parallel
       bench_    bench_   bench_   bench_
       mixed.c  baselines quest.c  qulacs.cpp
            │        │        │        │
            └────────┴────┬───┴────────┘
                          │
         Step F ──────────┤   ← depends only on A; can run in parallel with B–E
       (CLI parsing)      │
                          │
         Step G ──────────┤   ← depends only on A; can run in parallel with B–F
       (output fns)       │
                          │
                          ▼
               ┌─────────────────────┐
               │   Step H: main()    │  ← depends on B, C, D, E, F, G (all prior)
               │  orchestration loop │
               └──────────┬──────────┘
                           │
               ┌───────────┴───────────┐
               ▼                       ▼
            Step I                  Step J
         docs/BENCHMARK.md      thesis notes
         (depends on A,F,G,H)   (depends on H,I)
```

**Parallel batches:**

- Batch 0 (prerequisite): Step A
- Batch 1 (all independent, all depend only on A): Steps B, C, D, E, F, G
- Batch 2 (integration): Step H (depends on all of Batch 1)
- Batch 3 (docs): Step I (depends on A, F, G, H)
- Batch 4 (thesis): Step J (depends on H, I)

---

## Step A — `bench.h` additions

**Depends on:** nothing
**File:** `benchmark/include/bench.h`
**Estimated time:** 1 h

### What to add

**1. Per-qubit default constants** — insert after the existing `BENCH_MAX_RUNS` line:

```c
#define BENCH_PQ_DEFAULT_RUNS        5    /**< Default --runs when --per-qubit is active */
#define BENCH_PQ_DEFAULT_ITERATIONS  100  /**< Default --iterations when --per-qubit is active */
```

**2. Sentinel constant** — insert near the top of the Types section, before the structs:

```c
/** @brief Sentinel value for bench_result_perq_t.target2 when the gate is single-qubit. */
#define QUBIT_NONE ((qubit_t)0xFF)
```

Note: `qubit_t` is `uint8_t` (confirmed via `q_types.h`). `0xFF` is safe as a sentinel
because no valid qubit index reaches 255; the maximum supported qubit count is 30.

**3. `bench_result_perq_t` struct** — insert after `bench_result_t`:

```c
/**
 * @brief Per-qubit benchmark result: timing for a single (gate, qubits, method, target) tuple.
 *
 * Parallel to bench_result_t but records a fixed target qubit rather than sweeping.
 * For 2Q gates, target = q1 and target2 = q2 (q1 < q2). For 1Q gates, target2 = QUBIT_NONE.
 */
typedef struct bench_result_perq {
    qubit_t    qubits;           /**< Number of qubits in the state */
    dim_t      dim;              /**< Hilbert space dimension (2^qubits) */
    const char *gate_name;       /**< Name of the gate tested */
    const char *method;          /**< Implementation method name */
    qubit_t    target;           /**< Target qubit (q1 for 2Q gates) */
    qubit_t    target2;          /**< Second qubit for 2Q gates; QUBIT_NONE for 1Q */
    int        is_2q;            /**< 1 if this is a 2Q measurement, 0 otherwise */
    double     time_ms;          /**< Mean total time per run (ms) */
    double     time_ms_std;      /**< Sample std dev across runs, Bessel-corrected (ms) */
    double     time_ms_min;      /**< Minimum run time (ms) */
    double     time_ms_median;   /**< Median run time (ms) */
    double     time_ms_cv;       /**< Coefficient of variation = std/mean x 100 (%) */
    double     ops_per_sec;      /**< Gate applications per second, from mean time */
    double     time_per_gate_ms; /**< Mean time per gate application: time_ms / iterations */
    size_t     memory_bytes;     /**< Memory used for state representation */
    int        iterations;       /**< Gate calls per timed run */
    int        runs;             /**< Number of independent timing runs */
} bench_result_perq_t;
```

Note: `sweep_size` is omitted from `bench_result_perq_t` because per-qubit results always
have exactly 1 gate application per iteration (sweep_size is implicitly 1). `time_per_gate_ms`
is therefore simply `time_ms / iterations`.

**4. `per_qubit` field in `bench_options_t`** — add to the existing struct:

```c
    int per_qubit;   /**< 1 if --per-qubit flag was given */
```

**5. `_at` adapter declarations** — add a new section after the existing 2Q benchmark
function declarations, before the Output functions section:

```c
/*
 * =====================================================================================================================
 * Per-qubit adapter functions — fixed target, single measurement (bench_mixed.c, bench_baselines.c)
 * =====================================================================================================================
 */

/* Single-qubit, fixed target */
bench_result_perq_t bench_qlib_packed_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t),
    qubit_t target, int iterations, int warmup, int runs);

bench_result_perq_t bench_qlib_tiled_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t),
    qubit_t target, int iterations, int warmup, int runs);

bench_result_perq_t bench_blas_dense_at(qubit_t qubits, const char *gate_name,
    const cplx_t gate_mat[4],
    qubit_t target, int iterations, int warmup, int runs);

/* Two-qubit, fixed pair */
bench_result_perq_t bench_qlib_packed_2q_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t, qubit_t),
    qubit_t q1, qubit_t q2, int iterations, int warmup, int runs);

bench_result_perq_t bench_qlib_tiled_2q_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t, qubit_t),
    qubit_t q1, qubit_t q2, int iterations, int warmup, int runs);

bench_result_perq_t bench_blas_dense_2q_at(qubit_t qubits, const char *gate_name,
    const cplx_t gate_mat[16],
    qubit_t q1, qubit_t q2, int iterations, int warmup, int runs);

#ifdef WITH_QUEST
bench_result_perq_t bench_quest_at(qubit_t qubits, const char *gate_name,
    qubit_t target, int iterations, int warmup, int runs);

bench_result_perq_t bench_quest_2q_at(qubit_t qubits, const char *gate_name,
    qubit_t q1, qubit_t q2, int iterations, int warmup, int runs);
#endif

#ifdef WITH_QULACS
bench_result_perq_t bench_qulacs_at(qubit_t qubits, const char *gate_name,
    qubit_t target, int iterations, int warmup, int runs);

bench_result_perq_t bench_qulacs_2q_at(qubit_t qubits, const char *gate_name,
    qubit_t q1, qubit_t q2, int iterations, int warmup, int runs);
#endif
```

### Constraints

- Do NOT modify any `.c` or `.cpp` files in this step.
- Do NOT change any existing declarations or definitions.
- The `bench_options_t` zero-initializer in `bench_parse_options()` (bench_main.c) will
  pick up `per_qubit = 0` automatically because it uses designated initializers with all
  other fields named; no change to bench_main.c is needed in this step.

---

## Step B — `bench_mixed.c` additions

**Depends on:** Step A
**File:** `benchmark/src/bench_mixed.c`
**Estimated time:** 2 h

### Context

The existing `run_1q_timing()` and `run_2q_timing()` internal functions sweep over all
qubits/pairs per iteration. The new `_at` variants fix the target and time only that
single application per iteration. State is re-initialized once before the warmup; it is
NOT re-initialized per run (same pattern as the sweep functions).

### What to add

Append to the end of `bench_mixed.c`, after the existing 2Q implementations.

**1. `run_1q_timing_at()` internal function:**

```c
/**
 * @brief Shared timing loop for single-qubit per-qubit benchmarks.
 *
 * Times gate_fn applied to a single fixed target qubit for `iterations` calls
 * per run, across `runs` independent runs. State is initialised once before
 * warmup and shared across all runs (same as run_1q_timing).
 */
static bench_result_perq_t run_1q_timing_at(qubit_t qubits, const char *gate_name,
                                             void (*gate_fn)(state_t*, qubit_t),
                                             qubit_t target,
                                             int iterations, int warmup, int runs,
                                             const qlib_bench_desc_t *desc) {
    bench_result_perq_t result = {0};
    result.qubits       = qubits;
    result.dim          = (dim_t)1 << qubits;
    result.gate_name    = gate_name;
    result.method       = desc->method;
    result.target       = target;
    result.target2      = QUBIT_NONE;
    result.is_2q        = 0;
    result.iterations   = iterations;
    result.memory_bytes = desc->memory_bytes;

    state_t state = {0};
    init_random_mixed_state(&state, qubits, desc->state_type);
    if (!state.data) return result;

    for (int i = 0; i < warmup; ++i)
        gate_fn(&state, target);

    double *run_times = malloc(runs * sizeof(double));
    if (!run_times) { state_free(&state); return result; }
    for (int r = 0; r < runs; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i)
            gate_fn(&state, target);
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }
    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    free(run_times);

    result.time_ms          = stats.mean;
    result.time_ms_std      = stats.std_dev;
    result.time_ms_min      = stats.min;
    result.time_ms_median   = stats.median;
    result.time_ms_cv       = stats.cv;
    result.runs             = runs;
    result.ops_per_sec      = (stats.mean > 0.0)
        ? (double)iterations / (stats.mean / 1000.0) : 0.0;
    result.time_per_gate_ms = (iterations > 0)
        ? result.time_ms / (double)iterations : 0.0;

    state_free(&state);
    return result;
}
```

**2. `run_2q_timing_at()` internal function:**

```c
/**
 * @brief Shared timing loop for two-qubit per-qubit benchmarks.
 *
 * Times gate_fn applied to a single fixed pair (q1, q2) for `iterations` calls
 * per run, across `runs` independent runs.
 */
static bench_result_perq_t run_2q_timing_at(qubit_t qubits, const char *gate_name,
                                             void (*gate_fn)(state_t*, qubit_t, qubit_t),
                                             qubit_t q1, qubit_t q2,
                                             int iterations, int warmup, int runs,
                                             const qlib_bench_desc_t *desc) {
    bench_result_perq_t result = {0};
    result.qubits       = qubits;
    result.dim          = (dim_t)1 << qubits;
    result.gate_name    = gate_name;
    result.method       = desc->method;
    result.target       = q1;
    result.target2      = q2;
    result.is_2q        = 1;
    result.iterations   = iterations;
    result.memory_bytes = desc->memory_bytes;

    state_t state = {0};
    init_random_mixed_state(&state, qubits, desc->state_type);
    if (!state.data) return result;

    for (int i = 0; i < warmup; ++i)
        gate_fn(&state, q1, q2);

    double *run_times = malloc(runs * sizeof(double));
    if (!run_times) { state_free(&state); return result; }
    for (int r = 0; r < runs; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i)
            gate_fn(&state, q1, q2);
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }
    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    free(run_times);

    result.time_ms          = stats.mean;
    result.time_ms_std      = stats.std_dev;
    result.time_ms_min      = stats.min;
    result.time_ms_median   = stats.median;
    result.time_ms_cv       = stats.cv;
    result.runs             = runs;
    result.ops_per_sec      = (stats.mean > 0.0)
        ? (double)iterations / (stats.mean / 1000.0) : 0.0;
    result.time_per_gate_ms = (iterations > 0)
        ? result.time_ms / (double)iterations : 0.0;

    state_free(&state);
    return result;
}
```

**3. Four public wrapper functions** (thin, mirroring the existing packed/tiled wrappers):

```c
bench_result_perq_t bench_qlib_packed_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t),
    qubit_t target, int iterations, int warmup, int runs)
{
    const qlib_bench_desc_t desc = {
        .state_type   = MIXED_PACKED,
        .memory_bytes = bench_packed_size(qubits),
        .method       = "qlib_packed",
    };
    return run_1q_timing_at(qubits, gate_name, gate_fn, target,
                            iterations, warmup, runs, &desc);
}

bench_result_perq_t bench_qlib_tiled_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t),
    qubit_t target, int iterations, int warmup, int runs)
{
    const qlib_bench_desc_t desc = {
        .state_type   = MIXED_TILED,
        .memory_bytes = bench_tiled_size(qubits),
        .method       = "qlib_tiled",
    };
    return run_1q_timing_at(qubits, gate_name, gate_fn, target,
                            iterations, warmup, runs, &desc);
}

bench_result_perq_t bench_qlib_packed_2q_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t, qubit_t),
    qubit_t q1, qubit_t q2, int iterations, int warmup, int runs)
{
    const qlib_bench_desc_t desc = {
        .state_type   = MIXED_PACKED,
        .memory_bytes = bench_packed_size(qubits),
        .method       = "qlib_packed",
    };
    return run_2q_timing_at(qubits, gate_name, gate_fn, q1, q2,
                            iterations, warmup, runs, &desc);
}

bench_result_perq_t bench_qlib_tiled_2q_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t, qubit_t),
    qubit_t q1, qubit_t q2, int iterations, int warmup, int runs)
{
    const qlib_bench_desc_t desc = {
        .state_type   = MIXED_TILED,
        .memory_bytes = bench_tiled_size(qubits),
        .method       = "qlib_tiled",
    };
    return run_2q_timing_at(qubits, gate_name, gate_fn, q1, q2,
                            iterations, warmup, runs, &desc);
}
```

### Constraints

- Do NOT modify `run_1q_timing()` or `run_2q_timing()` or any existing public functions.
- The `qlib_bench_desc_t` type is already defined in this file (file-local struct); the
  new static functions can use it directly.

---

## Step C — `bench_baselines.c` additions

**Depends on:** Step A
**File:** `benchmark/src/bench_baselines.c`
**Estimated time:** 2 h

### Context

The existing `bench_blas_dense()` builds one full gate matrix per qubit (array of
`qubits` matrices) and sweeps all targets per iteration. The `_at` variant builds
exactly one matrix for the fixed target and applies only that matrix per iteration.
The naive baseline is intentionally excluded from the per-qubit path: it is already
excluded from pgfplots output and is O(N³); the BLAS baseline suffices as the dense
reference.

### What to add

Append to the end of `bench_baselines.c`.

**1. `bench_blas_dense_at()`:**

```c
bench_result_perq_t bench_blas_dense_at(qubit_t qubits, const char *gate_name,
                                         const cplx_t gate_mat[4],
                                         qubit_t target,
                                         int iterations, int warmup, int runs) {
    bench_result_perq_t result = {0};
    result.qubits     = qubits;
    result.dim        = (dim_t)1 << qubits;
    result.gate_name  = gate_name;
    result.method     = "blas_dense";
    result.target     = target;
    result.target2    = QUBIT_NONE;
    result.is_2q      = 0;
    result.iterations = iterations;

    dim_t dim = result.dim;

    cplx_t *rho = malloc(dim * dim * sizeof(cplx_t));
    if (!rho) return result;
    for (dim_t i = 0; i < dim * dim; ++i)
        rho[i] = (double)rand() / RAND_MAX + ((double)rand() / RAND_MAX) * I;

    /* Build gate matrix only for the fixed target — outside the timed loop */
    cplx_t *U = build_full_gate_matrix(qubits, gate_mat, target);
    if (!U) { free(rho); return result; }

    cplx_t *tmp = malloc(dim * dim * sizeof(cplx_t));
    if (!tmp) { free(U); free(rho); return result; }

    result.memory_bytes = bench_dense_size(qubits);

    for (int i = 0; i < warmup; ++i)
        apply_uruh_blas(dim, U, rho, tmp);

    double *run_times = malloc(runs * sizeof(double));
    if (!run_times) { free(tmp); free(U); free(rho); return result; }
    for (int r = 0; r < runs; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i)
            apply_uruh_blas(dim, U, rho, tmp);
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }
    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    free(run_times);

    result.time_ms          = stats.mean;
    result.time_ms_std      = stats.std_dev;
    result.time_ms_min      = stats.min;
    result.time_ms_median   = stats.median;
    result.time_ms_cv       = stats.cv;
    result.runs             = runs;
    result.ops_per_sec      = (stats.mean > 0.0)
        ? (double)iterations / (stats.mean / 1000.0) : 0.0;
    result.time_per_gate_ms = (iterations > 0)
        ? result.time_ms / (double)iterations : 0.0;

    free(tmp); free(U); free(rho);
    return result;
}
```

**2. `bench_blas_dense_2q_at()`:**

```c
bench_result_perq_t bench_blas_dense_2q_at(qubit_t qubits, const char *gate_name,
                                             const cplx_t gate_mat[16],
                                             qubit_t q1, qubit_t q2,
                                             int iterations, int warmup, int runs) {
    bench_result_perq_t result = {0};
    result.qubits     = qubits;
    result.dim        = (dim_t)1 << qubits;
    result.gate_name  = gate_name;
    result.method     = "blas_dense";
    result.target     = q1;
    result.target2    = q2;
    result.is_2q      = 1;
    result.iterations = iterations;

    dim_t dim = result.dim;

    cplx_t *rho = malloc(dim * dim * sizeof(cplx_t));
    if (!rho) return result;
    for (dim_t i = 0; i < dim * dim; ++i)
        rho[i] = (double)rand() / RAND_MAX + ((double)rand() / RAND_MAX) * I;

    /* Build gate matrix only for the fixed pair — outside the timed loop */
    cplx_t *U = build_full_gate_matrix_2q(qubits, gate_mat, q1, q2);
    if (!U) { free(rho); return result; }

    cplx_t *tmp = malloc(dim * dim * sizeof(cplx_t));
    if (!tmp) { free(U); free(rho); return result; }

    result.memory_bytes = bench_dense_size(qubits);

    for (int i = 0; i < warmup; ++i)
        apply_uruh_blas(dim, U, rho, tmp);

    double *run_times = malloc(runs * sizeof(double));
    if (!run_times) { free(tmp); free(U); free(rho); return result; }
    for (int r = 0; r < runs; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i)
            apply_uruh_blas(dim, U, rho, tmp);
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }
    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    free(run_times);

    result.time_ms          = stats.mean;
    result.time_ms_std      = stats.std_dev;
    result.time_ms_min      = stats.min;
    result.time_ms_median   = stats.median;
    result.time_ms_cv       = stats.cv;
    result.runs             = runs;
    result.ops_per_sec      = (stats.mean > 0.0)
        ? (double)iterations / (stats.mean / 1000.0) : 0.0;
    result.time_per_gate_ms = (iterations > 0)
        ? result.time_ms / (double)iterations : 0.0;

    free(tmp); free(U); free(rho);
    return result;
}
```

### Constraints

- `build_full_gate_matrix()` and `build_full_gate_matrix_2q()` are already defined
  earlier in `bench_baselines.c` as file-local `static` functions; they are directly
  callable from the new functions appended to the same file.
- `apply_uruh_blas()` is likewise file-local static; same rule applies.
- Do NOT expose `build_full_gate_matrix*` or `apply_uruh_blas` in `bench.h`; they remain
  private to this translation unit.
- The dense baseline only applies up to 8 qubits (enforced in bench_main.c with
  `run_dense = (qubits <= 8)`). The `_at` functions themselves impose no qubit cap;
  the orchestration loop in Step H applies the cap.

---

## Step D — `bench_quest.c` additions

**Depends on:** Step A
**File:** `benchmark/src/bench_quest.c`
**Estimated time:** 1.5 h

### Context

The existing `bench_quest()` function dispatches on `gate_name` to select the correct
QuEST function pointer, then runs the sweep loop. The `_at` variants take the same
approach but fix the target and call only once per iteration.

Both new functions must be inside the `#ifdef WITH_QUEST` guard.

### What to add

Append inside the `#ifdef WITH_QUEST` block, after the existing `bench_quest()` function:

**1. `bench_quest_at()` — single-qubit, fixed target:**

```c
bench_result_perq_t bench_quest_at(qubit_t qubits, const char *gate_name,
                                    qubit_t target, int iterations, int warmup, int runs) {
    bench_result_perq_t result = {0};
    result.qubits     = qubits;
    result.dim        = (dim_t)1 << qubits;
    result.gate_name  = gate_name;
    result.method     = "quest";
    result.target     = target;
    result.target2    = QUBIT_NONE;
    result.is_2q      = 0;
    result.iterations = iterations;
    result.runs       = runs;

    Qureg rho = createDensityQureg((int)qubits);
    initDebugState(rho);

    dim_t dim = (dim_t)1 << qubits;
    result.memory_bytes = (size_t)dim * (size_t)dim * sizeof(double) * 2;

    void (*gate_fn_1q)(Qureg, int) = NULL;
    if      (strcmp(gate_name, "X") == 0) gate_fn_1q = applyPauliX;
    else if (strcmp(gate_name, "H") == 0) gate_fn_1q = applyHadamard;
    else if (strcmp(gate_name, "Z") == 0) gate_fn_1q = applyPauliZ;

    if (!gate_fn_1q) { destroyQureg(rho); return result; }

    int cnt = (runs < 200) ? runs : 200;
    double run_times[200];

    for (int i = 0; i < warmup; ++i)
        gate_fn_1q(rho, (int)target);

    for (int r = 0; r < cnt; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i)
            gate_fn_1q(rho, (int)target);
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }

    bench_run_stats_t stats = bench_compute_stats(run_times, cnt);
    result.time_ms          = stats.mean;
    result.time_ms_std      = stats.std_dev;
    result.time_ms_min      = stats.min;
    result.time_ms_median   = stats.median;
    result.time_ms_cv       = stats.cv;
    result.ops_per_sec      = (stats.mean > 0.0)
        ? (double)iterations / (stats.mean / 1000.0) : 0.0;
    result.time_per_gate_ms = (iterations > 0)
        ? result.time_ms / (double)iterations : 0.0;

    destroyQureg(rho);
    return result;
}
```

**2. `bench_quest_2q_at()` — two-qubit, fixed pair:**

```c
bench_result_perq_t bench_quest_2q_at(qubit_t qubits, const char *gate_name,
                                       qubit_t q1, qubit_t q2,
                                       int iterations, int warmup, int runs) {
    bench_result_perq_t result = {0};
    result.qubits     = qubits;
    result.dim        = (dim_t)1 << qubits;
    result.gate_name  = gate_name;
    result.method     = "quest";
    result.target     = q1;
    result.target2    = q2;
    result.is_2q      = 1;
    result.iterations = iterations;
    result.runs       = runs;

    Qureg rho = createDensityQureg((int)qubits);
    initDebugState(rho);

    dim_t dim = (dim_t)1 << qubits;
    result.memory_bytes = (size_t)dim * (size_t)dim * sizeof(double) * 2;

    void (*gate_fn_2q)(Qureg, int, int) = NULL;
    if      (strcmp(gate_name, "CX")   == 0) gate_fn_2q = applyControlledPauliX;
    else if (strcmp(gate_name, "SWAP") == 0) gate_fn_2q = applySwap;

    if (!gate_fn_2q) { destroyQureg(rho); return result; }

    int cnt = (runs < 200) ? runs : 200;
    double run_times[200];

    for (int i = 0; i < warmup; ++i)
        gate_fn_2q(rho, (int)q1, (int)q2);

    for (int r = 0; r < cnt; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i)
            gate_fn_2q(rho, (int)q1, (int)q2);
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }

    bench_run_stats_t stats = bench_compute_stats(run_times, cnt);
    result.time_ms          = stats.mean;
    result.time_ms_std      = stats.std_dev;
    result.time_ms_min      = stats.min;
    result.time_ms_median   = stats.median;
    result.time_ms_cv       = stats.cv;
    result.ops_per_sec      = (stats.mean > 0.0)
        ? (double)iterations / (stats.mean / 1000.0) : 0.0;
    result.time_per_gate_ms = (iterations > 0)
        ? result.time_ms / (double)iterations : 0.0;

    destroyQureg(rho);
    return result;
}
```

### Constraints

- Both functions must be added inside the `#ifdef WITH_QUEST` / `#endif /* WITH_QUEST */`
  guards that already wrap this translation unit.
- The init-once pattern (`bench_quest_init` / `bench_quest_cleanup`) is unchanged. The
  `_at` functions call `createDensityQureg` / `destroyQureg` per call, matching the
  existing `bench_quest()` pattern.
- Do NOT add `#ifdef WITH_QUEST` guards around individual functions — the file-level
  guard at line 24 already covers all code in the file.

---

## Step E — `bench_qulacs.cpp` additions

**Depends on:** Step A
**File:** `benchmark/src/bench_qulacs.cpp`
**Estimated time:** 1.5 h

### Context

The existing `bench_qulacs()` dispatches on `gate_name`, selects a Qulacs csim function
pointer (typed `void (*)(UINT, CTYPE*, ITYPE)` for 1Q or `void (*)(UINT, UINT, CTYPE*, ITYPE)`
for 2Q), and runs the sweep loop. The `_at` variants fix the target.

Both new functions must be inside the `#ifdef WITH_QULACS` guard and inside the
existing `extern "C"` block.

### What to add

Append inside the `extern "C" { ... }` block and the enclosing `#ifdef WITH_QULACS`
guard, after the existing `bench_qulacs()` function:

**1. `bench_qulacs_at()` — single-qubit, fixed target:**

```cpp
bench_result_perq_t bench_qulacs_at(qubit_t qubits, const char *gate_name,
                                     qubit_t target, int iterations, int warmup, int runs) {
    bench_result_perq_t result = {0};
    result.qubits     = qubits;
    result.dim        = (dim_t)1 << qubits;
    result.gate_name  = gate_name;
    result.method     = "qulacs";
    result.target     = target;
    result.target2    = QUBIT_NONE;
    result.is_2q      = 0;
    result.iterations = iterations;
    result.runs       = runs;

    ITYPE dim = static_cast<ITYPE>(1) << qubits;
    size_t matrix_size = static_cast<size_t>(dim) * static_cast<size_t>(dim) * sizeof(CTYPE);
    CTYPE *rho = static_cast<CTYPE*>(std::malloc(matrix_size));
    if (!rho) return result;

    size_t n_elems = static_cast<size_t>(dim) * static_cast<size_t>(dim);
    for (size_t i = 0; i < n_elems; ++i) {
        double re = (double)std::rand() / RAND_MAX - 0.5;
        double im = (double)std::rand() / RAND_MAX - 0.5;
        rho[i] = CTYPE(re, im);
    }

    result.memory_bytes = matrix_size;

    void (*gate_fn_1q)(UINT, CTYPE*, ITYPE) = nullptr;
    if      (std::strcmp(gate_name, "X") == 0) gate_fn_1q = dm_X_gate;
    else if (std::strcmp(gate_name, "H") == 0) gate_fn_1q = dm_H_gate;
    else if (std::strcmp(gate_name, "Z") == 0) gate_fn_1q = dm_Z_gate;

    if (!gate_fn_1q) { std::free(rho); return result; }

    int cnt = (runs < 200) ? runs : 200;
    double run_times[200];

    for (int i = 0; i < warmup; ++i)
        gate_fn_1q(static_cast<UINT>(target), rho, dim);

    for (int r = 0; r < cnt; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i)
            gate_fn_1q(static_cast<UINT>(target), rho, dim);
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }

    bench_run_stats_t stats = bench_compute_stats(run_times, cnt);
    result.time_ms          = stats.mean;
    result.time_ms_std      = stats.std_dev;
    result.time_ms_min      = stats.min;
    result.time_ms_median   = stats.median;
    result.time_ms_cv       = stats.cv;
    result.ops_per_sec      = (stats.mean > 0.0)
        ? static_cast<double>(iterations) / (stats.mean / 1000.0) : 0.0;
    result.time_per_gate_ms = (iterations > 0)
        ? result.time_ms / static_cast<double>(iterations) : 0.0;

    std::free(rho);
    return result;
}
```

**2. `bench_qulacs_2q_at()` — two-qubit, fixed pair:**

```cpp
bench_result_perq_t bench_qulacs_2q_at(qubit_t qubits, const char *gate_name,
                                         qubit_t q1, qubit_t q2,
                                         int iterations, int warmup, int runs) {
    bench_result_perq_t result = {0};
    result.qubits     = qubits;
    result.dim        = (dim_t)1 << qubits;
    result.gate_name  = gate_name;
    result.method     = "qulacs";
    result.target     = q1;
    result.target2    = q2;
    result.is_2q      = 1;
    result.iterations = iterations;
    result.runs       = runs;

    ITYPE dim = static_cast<ITYPE>(1) << qubits;
    size_t matrix_size = static_cast<size_t>(dim) * static_cast<size_t>(dim) * sizeof(CTYPE);
    CTYPE *rho = static_cast<CTYPE*>(std::malloc(matrix_size));
    if (!rho) return result;

    size_t n_elems = static_cast<size_t>(dim) * static_cast<size_t>(dim);
    for (size_t i = 0; i < n_elems; ++i) {
        double re = (double)std::rand() / RAND_MAX - 0.5;
        double im = (double)std::rand() / RAND_MAX - 0.5;
        rho[i] = CTYPE(re, im);
    }

    result.memory_bytes = matrix_size;

    void (*gate_fn_2q)(UINT, UINT, CTYPE*, ITYPE) = nullptr;
    if      (std::strcmp(gate_name, "CX")   == 0) gate_fn_2q = dm_CNOT_gate;
    else if (std::strcmp(gate_name, "SWAP") == 0) gate_fn_2q = dm_SWAP_gate;

    if (!gate_fn_2q) { std::free(rho); return result; }

    int cnt = (runs < 200) ? runs : 200;
    double run_times[200];

    for (int i = 0; i < warmup; ++i)
        gate_fn_2q(static_cast<UINT>(q1), static_cast<UINT>(q2), rho, dim);

    for (int r = 0; r < cnt; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i)
            gate_fn_2q(static_cast<UINT>(q1), static_cast<UINT>(q2), rho, dim);
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }

    bench_run_stats_t stats = bench_compute_stats(run_times, cnt);
    result.time_ms          = stats.mean;
    result.time_ms_std      = stats.std_dev;
    result.time_ms_min      = stats.min;
    result.time_ms_median   = stats.median;
    result.time_ms_cv       = stats.cv;
    result.ops_per_sec      = (stats.mean > 0.0)
        ? static_cast<double>(iterations) / (stats.mean / 1000.0) : 0.0;
    result.time_per_gate_ms = (iterations > 0)
        ? result.time_ms / static_cast<double>(iterations) : 0.0;

    std::free(rho);
    return result;
}
```

### Constraints

- Both functions must remain inside the `extern "C" { ... }` block (for C linkage) and
  inside the `#ifdef WITH_QULACS` / `#endif` guards.
- The `CTYPE` and `ITYPE` types are from `<csim/type.hpp>`, already included at the top.
- `dm_X_gate`, `dm_H_gate`, `dm_Z_gate`, `dm_CNOT_gate`, `dm_SWAP_gate` are from
  `<csim/update_ops_dm.hpp>`, already included at the top.
- The `UINT` cast pattern matches the existing code exactly.

---

## Step F — `bench_main.c`: CLI parsing

**Depends on:** Step A
**File:** `benchmark/src/bench_main.c`
**Estimated time:** 1 h

### What to change

**1. Initialization in `bench_parse_options()`** — add `per_qubit = 0` to the designated
initializer for `bench_options_t opts`:

```c
    bench_options_t opts = {
        ...
        .per_qubit       = 0,
    };
```

**2. Argument parsing loop** — add the `--per-qubit` branch in the `for` loop in
`bench_parse_options()`, alongside the existing flag branches:

```c
        } else if (strcmp(argv[i], "--per-qubit") == 0) {
            opts.per_qubit = 1;
```

**3. Post-parse default override** — after the existing validation block, add:

```c
    /*
     * --per-qubit uses lighter defaults to keep total runtime manageable.
     * Only override iterations/runs if the user did not specify them explicitly.
     * Detection: compare against BENCH_DEFAULT_* values. This is a heuristic
     * but is correct for the normal use case.
     */
    if (opts.per_qubit) {
        if (opts.iterations == BENCH_DEFAULT_ITERATIONS)
            opts.iterations = BENCH_PQ_DEFAULT_ITERATIONS;
        if (opts.runs == BENCH_DEFAULT_RUNS)
            opts.runs = BENCH_PQ_DEFAULT_RUNS;
    }
```

**4. Mutual exclusion check** — add after the `per_qubit` default override:

```c
    if (opts.per_qubit && opts.pgfplots_output) {
        fprintf(stderr, "error: --per-qubit and --pgfplots are mutually exclusive\n");
        exit(1);
    }
```

**5. `print_usage()` update** — add the `--per-qubit` line to the help text, after
`--pgfplots`:

```c
    printf("  --per-qubit      Time each target qubit independently (default runs=%d, iters=%d)\n",
           BENCH_PQ_DEFAULT_RUNS, BENCH_PQ_DEFAULT_ITERATIONS);
```

### Risk note

The "did the user explicitly set `--runs`?" detection using `== BENCH_DEFAULT_RUNS` is
a heuristic: if the user happens to pass `--runs 10` (the default) plus `--per-qubit`,
the override will incorrectly kick in and lower it to 5. This is the same approach used
by many CLI tools and is acceptable for a benchmark binary. An alternative is to add a
`runs_explicit` boolean field to `bench_options_t`, set it in the `--runs` parsing
branch, and check it instead of comparing against the default. Implementers may choose
either approach; document the chosen approach in a comment.

---

## Step G — `bench_main.c`: per-qubit output functions

**Depends on:** Step A
**File:** `benchmark/src/bench_main.c`
**Estimated time:** 2 h

### Context

Three new functions are needed. They only emit output; they do not call any adapter.
They should be placed in the Output functions section of `bench_main.c`, after the
existing `bench_print_csv()` function.

### Bandwidth column derivation

The `bw_gbs` (GB/s) column represents the memory bandwidth consumed by a single gate
application, estimated as:

```
bw_gbs = (2 * memory_bytes) / (time_per_gate_ms * 1e-3) / 1e9
```

Rationale: a gate application reads the full density matrix once (rho) and writes it
once (rho updated in-place), so the minimum memory traffic is 2 × `memory_bytes`.
`time_per_gate_ms * 1e-3` converts ms to seconds.

Combined: `bw_gbs = 2.0 * memory_bytes / (time_per_gate_ms * 1e6)`.

If `time_per_gate_ms <= 0`, emit 0.0.

### What to add

**1. CSV header function:**

```c
void bench_print_perq_csv_header(void) {
    printf("qubits,dim,gate,method,runs,iterations,"
           "target,target2,is_2q,"
           "time_ms_mean,time_ms_std,time_ms_min,time_ms_median,time_ms_cv,"
           "time_per_gate_ms,ops_per_sec,bw_gbs,memory_bytes\n");
}
```

**2. CSV row function:**

```c
void bench_print_perq_csv(const bench_result_perq_t *r) {
    double bw_gbs = (r->time_per_gate_ms > 0.0)
        ? 2.0 * (double)r->memory_bytes / (r->time_per_gate_ms * 1e6)
        : 0.0;
    const char *t2_str = (r->is_2q) ? "" : "";  /* target2 is always numeric below */
    if (r->is_2q) {
        printf("%u,%lld,%s,%s,%d,%d,%u,%u,%d,"
               "%.6f,%.6f,%.6f,%.6f,%.3f,%.9f,%.0f,%.3f,%zu\n",
               r->qubits, (long long)r->dim, r->gate_name, r->method,
               r->runs, r->iterations,
               r->target, r->target2, r->is_2q,
               r->time_ms, r->time_ms_std, r->time_ms_min,
               r->time_ms_median, r->time_ms_cv,
               r->time_per_gate_ms, r->ops_per_sec, bw_gbs, r->memory_bytes);
    } else {
        printf("%u,%lld,%s,%s,%d,%d,%u,,%d,"
               "%.6f,%.6f,%.6f,%.6f,%.3f,%.9f,%.0f,%.3f,%zu\n",
               r->qubits, (long long)r->dim, r->gate_name, r->method,
               r->runs, r->iterations,
               r->target, r->is_2q,
               r->time_ms, r->time_ms_std, r->time_ms_min,
               r->time_ms_median, r->time_ms_cv,
               r->time_per_gate_ms, r->ops_per_sec, bw_gbs, r->memory_bytes);
    }
}
```

Note: for 1Q gates, `target2` is left empty in the CSV (`,` with no value between
commas). This keeps the column count constant: target2 is always present as a column
but may be empty. Alternatively implementers may emit the string `"none"` — pick one
convention and document it.

**3. Console table function:**

```c
/**
 * @brief Print per-qubit results as a human-readable table.
 *
 * Groups results by (gate, qubits, method) and prints one row per target.
 * Callers pass a flat array of results for a single (gate, qubits, method) combination.
 *
 * @param results  Array of per-qubit results (all same gate/qubits/method)
 * @param n        Number of entries
 */
void bench_print_perq_console(const bench_result_perq_t *results, int n) {
    if (n == 0) return;
    const bench_result_perq_t *r0 = &results[0];

    printf("\n  Gate: %-6s  qubits: %u  method: %s\n",
           r0->gate_name, r0->qubits, r0->method);

    if (r0->is_2q) {
        printf("  %-6s  %-6s  %10s  %10s  %8s  %12s\n",
               "q1", "q2", "mean(ms)", "min(ms)", "CV(%)", "ops/sec");
        printf("  %-6s  %-6s  %10s  %10s  %8s  %12s\n",
               "──────", "──────", "──────────", "──────────", "────────", "────────────");
        for (int i = 0; i < n; ++i) {
            const bench_result_perq_t *r = &results[i];
            printf("  %-6u  %-6u  %10.4f  %10.4f  %8.2f  %12.0f\n",
                   r->target, r->target2,
                   r->time_ms, r->time_ms_min, r->time_ms_cv, r->ops_per_sec);
        }
    } else {
        printf("  %-6s  %10s  %10s  %8s  %12s\n",
               "target", "mean(ms)", "min(ms)", "CV(%)", "ops/sec");
        printf("  %-6s  %10s  %10s  %8s  %12s\n",
               "──────", "──────────", "──────────", "────────", "────────────");
        for (int i = 0; i < n; ++i) {
            const bench_result_perq_t *r = &results[i];
            printf("  %-6u  %10.4f  %10.4f  %8.2f  %12.0f\n",
                   r->target,
                   r->time_ms, r->time_ms_min, r->time_ms_cv, r->ops_per_sec);
        }
    }
}
```

### Declarations in bench.h (already handled by Step A)

Step A declares these three functions in `bench.h`. Step G only provides the
implementations.

Add to the Output functions section of `bench.h` (in Step A):

```c
void bench_print_perq_csv_header(void);
void bench_print_perq_csv(const bench_result_perq_t *r);
void bench_print_perq_console(const bench_result_perq_t *results, int n);
```

---

## Step H — `bench_main.c`: per-qubit orchestration loop

**Depends on:** Steps B, C, D, E, F, G (all must be complete)
**File:** `benchmark/src/bench_main.c`
**Estimated time:** 3 h

This is the most complex step. It wires together everything from the preceding steps
and must leave the existing non-per-qubit code path completely unchanged.

### Design decisions to make before implementing

**Array sizing:** The per-qubit result arrays must be sized to accommodate the maximum
number of targets. For 1Q gates, max targets = `max_qubits`. For 2Q gates, max targets
= `MAX_PAIRS` (128). Define:

```c
#define MAX_TARGETS  MAX_PAIRS   /* 128; sufficient for up to 16 qubits */
```

**Method count for per-qubit:** Per-qubit output does not need the `naive_loop` method
(excluded for the same reason as in pgfplots). The effective per-qubit method list is:
qlib_packed, qlib_tiled, blas_dense, quest, qulacs — matching the existing 5-method
set used in pgfplots output. Define `NUM_PQ_METHODS 5` if useful, or just use separate
arrays.

**Result storage:** Use a flat 2D array:
`bench_result_perq_t perq[MAX_TARGETS][NUM_METHODS]`
with the same `M_QLIB=0, M_QLIB_TILED=1, M_BLAS=2, M_QUEST=3, M_QULACS=4` enum that
already exists (note: in the existing code `M_QUEST=4, M_QULACS=5`; for per-qubit use
only 5 methods so re-number or keep existing values and skip index 3 which was
`M_NAIVE`).

Concrete recommendation: use the existing enum values unchanged:
- `M_QLIB = 0`, `M_QLIB_TILED = 1`, `M_BLAS = 2`, `M_QUEST = 4`, `M_QULACS = 5`
- For per-qubit, skip index 3 (M_NAIVE). Array size remains `NUM_METHODS = 6`.

### Where to insert the per-qubit branch

In `main()`, immediately after `bench_parse_options()` returns and the existing
pgfplots allocation block. The per-qubit path is a complete early return:

```c
    if (opts.per_qubit) {
        /* ... per-qubit implementation ... */
        /* returns 0 at the end; does not fall through to the existing loop */
        return perq_main(&opts);
    }
```

Define a static helper function `perq_main()` above `main()` to keep the function body
manageable. Prototype: `static int perq_main(const bench_options_t *opts)`.

### Structure of `perq_main()`

```
perq_main(opts):
    if csv_output: bench_print_perq_csv_header()
    else: print banner

    bench_quest_init()  [if WITH_QUEST]

    for qubits in [opts->min_qubits .. opts->max_qubits step opts->step]:
        dim = 1 << qubits
        run_dense = (qubits <= 8)

        if !csv_output: print "qubits: N, dim: D"

        /* ── 1Q gates ── */
        for each gate g in gates[]:
            n_targets = qubits
            bench_result_perq_t perq[MAX_TARGETS][NUM_METHODS] = {0}

            for target = 0 to qubits-1:
                perq[target][M_QLIB]       = bench_qlib_packed_at(...)
                perq[target][M_QLIB_TILED] = bench_qlib_tiled_at(...)
                if run_dense:
                    perq[target][M_BLAS]   = bench_blas_dense_at(...)
                #ifdef WITH_QUEST
                    perq[target][M_QUEST]  = bench_quest_at(...)
                #endif
                #ifdef WITH_QULACS
                    perq[target][M_QULACS] = bench_qulacs_at(...)
                #endif

            /* Output */
            if csv_output:
                for target in 0..qubits-1:
                    for each active method m:
                        bench_print_perq_csv(&perq[target][m])
            else:
                for each active method m:
                    collect perq[0..qubits-1][m] into a flat array of n_targets results
                    bench_print_perq_console(flat_array, n_targets)

        /* ── 2Q gates ── */
        qubit_t q1s[MAX_PAIRS], q2s[MAX_PAIRS]
        n_pairs = build_all_pairs(qubits, q1s, q2s)

        for each gate g in gates_2q[]:
            bench_result_perq_t perq[MAX_PAIRS][NUM_METHODS] = {0}

            for p = 0 to n_pairs-1:
                perq[p][M_QLIB]       = bench_qlib_packed_2q_at(... q1s[p], q2s[p] ...)
                if gates_2q[g].has_tiled:
                    perq[p][M_QLIB_TILED] = bench_qlib_tiled_2q_at(...)
                if run_dense:
                    perq[p][M_BLAS]   = bench_blas_dense_2q_at(...)
                #ifdef WITH_QUEST
                    perq[p][M_QUEST]  = bench_quest_2q_at(...)
                #endif
                #ifdef WITH_QULACS
                    perq[p][M_QULACS] = bench_qulacs_2q_at(...)
                #endif

            /* Output */
            if csv_output:
                for p in 0..n_pairs-1:
                    for each active method m:
                        bench_print_perq_csv(&perq[p][m])
            else:
                for each active method m:
                    collect perq[0..n_pairs-1][m] into flat array
                    bench_print_perq_console(flat_array, n_pairs)

    bench_quest_cleanup()  [if WITH_QUEST]
    return 0
```

### Stack allocation note

`bench_result_perq_t perq[MAX_PAIRS][NUM_METHODS]` = 128 × 6 × ~160 bytes ≈ 120 KB
per 2Q gate configuration. This is within reasonable stack limits on macOS/Linux
(typical 8 MB stack). However to be safe, declare these arrays at function scope in
`perq_main()` rather than inside the inner gate loop, and zero them with `memset` at
the start of each gate iteration.

Alternatively, heap-allocate once at the start of `perq_main()` and reuse:

```c
bench_result_perq_t (*perq)[NUM_METHODS] =
    calloc(MAX_TARGETS, sizeof(*perq));
if (!perq) { fprintf(stderr, "..."); return 1; }
/* ... loops ... */
free(perq);
```

This avoids stack concerns entirely. Either approach is acceptable; document the choice.

### Constraints

- The existing `for (qubit_t qubits = ...)` main loop, the `pgfplots_data_t` struct,
  and all existing output paths must remain completely untouched.
- `bench_quest_init()` / `bench_quest_cleanup()` calls in `perq_main()` are independent
  of those in the existing path; each path calls them exactly once.
- The `SET_SWEEP` macro is NOT used in the per-qubit path (per-qubit results have
  implicit sweep_size=1; `time_per_gate_ms` is computed in the `_at` functions directly).

---

## Step I — `docs/BENCHMARK.md` update

**Depends on:** Steps A, F, G, H (all changes final)
**File:** `docs/BENCHMARK.md`
**Estimated time:** 1 h

### What to add / update

**1. CLI Reference table** — add `--per-qubit` row:

| Flag | Argument | Default | Description |
|------|----------|---------|-------------|
| `--per-qubit` | — | off | Time each target qubit independently; sets `--runs=5 --iterations=100` if not overridden; mutually exclusive with `--pgfplots` |

**2. Configuration constants** — add to the constants table:

| Constant | Value | Description |
|----------|-------|-------------|
| `BENCH_PQ_DEFAULT_RUNS` | 5 | Default `--runs` when `--per-qubit` is active |
| `BENCH_PQ_DEFAULT_ITERATIONS` | 100 | Default `--iterations` when `--per-qubit` is active |
| `QUBIT_NONE` | 0xFF | Sentinel for `bench_result_perq_t.target2` on 1Q measurements |

**3. Key Data Structures section** — add `bench_result_perq_t` entry documenting all
fields, noting the difference from `bench_result_t` (no `sweep_size`; adds `target`,
`target2`, `is_2q`; `time_per_gate_ms` = time_ms/iterations rather than
time_ms/(iterations×sweep_size)).

**4. CSV Output section** — add a subsection for the per-qubit CSV format:

- Header: `qubits,dim,gate,method,runs,iterations,target,target2,is_2q,time_ms_mean,...,bw_gbs,memory_bytes`
- Document `bw_gbs` derivation: `2 × memory_bytes / (time_per_gate_ms × 10⁶)` (GB/s)
- Note: `target2` column is empty (not present) for 1Q gate rows.

**5. Adapter functions section** — list all `_at` signatures (or add them to the
existing function listing table).

**6. Status** — update module status from "In Progress" to reflect the new capability.

---

## Step J — Thesis notes

**Depends on:** Steps H, I
**File:** `~/Projects/thesis/2.Simulation/notes/BENCHMARK_THESIS_NOTES.md`
**Estimated time:** 1.5 h

### What to add

Add a new section (e.g. `## Per-Qubit Analysis`) covering:

1. **What it measures**: `--per-qubit` isolates the cost of applying a gate to a
   specific (target qubit, state representation) pair, removing the averaging effect
   of the sweep loop. This exposes whether gate cost depends on target position within
   the density matrix layout.

2. **Interpreting target-position differences**: In packed lower-triangular storage,
   lower-index qubits modulate inner indices of the data structure, while higher-index
   qubits modulate outer indices. A gate on qubit 0 touches many cache lines in a
   scattered pattern; a gate on qubit n-1 touches a contiguous block. The
   `--per-qubit` output shows whether this asymmetry is measurable in practice (it
   typically becomes visible at qubit counts where the density matrix no longer fits
   in L3 cache).

3. **Bandwidth column**: `bw_gbs` = effective memory bandwidth (GB/s) estimated as
   2 × (state size) / (time per gate application). This is a lower bound on actual
   bandwidth (assumes a single read-modify-write pass). Use this to compare how close
   each implementation comes to hardware memory bandwidth limits (measurable via
   `stream` benchmark). A high `bw_gbs` relative to peak bandwidth means the
   implementation is memory-bound; a low value suggests compute-bound or cache-resident.

4. **Citing per-qubit data in a thesis chapter**: The per-qubit CSV is suitable for
   pgfplots `\addplot table` with `x=target, y=time_per_gate_ms` to produce a
   target-position sweep plot for a fixed (gate, qubits, method) combination. For
   thesis narrative: describe whether packed, tiled, QuEST, and Qulacs show flat or
   position-dependent profiles. A flat profile supports the claim that the implementation
   is target-position agnostic. A sloped profile is a finding worth discussing.

5. **Runtime estimates** at the per-qubit defaults (5 runs × 100 iterations):
   compute approximate total measurement count and wall-clock estimate at the default
   qubit range (2–12), and compare to the ~13.5 h full-default estimate to justify
   the lighter defaults.

---

## Time Estimates

| Step | Description | Estimate |
|------|-------------|----------|
| A | bench.h additions | 1 h |
| B | bench_mixed.c additions | 2 h |
| C | bench_baselines.c additions | 2 h |
| D | bench_quest.c additions | 1.5 h |
| E | bench_qulacs.cpp additions | 1.5 h |
| F | bench_main.c CLI parsing | 1 h |
| G | bench_main.c output functions | 2 h |
| H | bench_main.c orchestration | 3 h |
| I | docs/BENCHMARK.md update | 1 h |
| J | Thesis notes | 1.5 h |
| **Total sequential** | | **16.5 h** |
| **Critical path** (A→{B,C,D,E,F,G}→H→I→J) | | **9.5 h** |
| **With max parallelism** (Batch 0 + 1 in parallel + H + I + J) | | **A(1h) + H(3h) + I(1h) + J(1.5h) = 6.5 h elapsed** |

The critical path is A → any-one-of(B,C,D,E,F,G) → H → I → J = 9.5 h.
With 6 agents running Batch 1 in parallel, the Batch 1 wall time collapses to the
longest of {B=2h, C=2h, D=1.5h, E=1.5h, F=1h, G=2h} = 2 h, giving a practical
elapsed time of approximately **7.5 h** for the full feature.

---

## Risks and Ambiguities

### R1 — Default override heuristic (Step F)
The check `if (opts.iterations == BENCH_DEFAULT_ITERATIONS)` to detect whether the
user explicitly set `--iterations` is a heuristic. A user who runs `--per-qubit
--iterations 1000` intends to keep the default sweep count, but the check will
correctly detect the explicit value and not override it. The edge case is `--per-qubit
--iterations 1000` where 1000 happens to equal `BENCH_DEFAULT_ITERATIONS`: the system
will (correctly) not override. The only ambiguous case is `--per-qubit --runs 10`:
`10 == BENCH_DEFAULT_RUNS` so the system will override to 5 even though the user
specified 10. Adding a `runs_explicit` boolean flag to `bench_options_t` eliminates
this entirely; the plan leaves this as a decision for the implementer.

### R2 — QuEST argument order for 2Q gates (Step D)
`applyControlledPauliX(Qureg, int control, int target)` — the argument convention in
QuEST v4 is `(control, target)`. In the existing `bench_quest()` sweep code, pairs are
passed as `gate_fn_2q(rho, q1s[p], q2s[p])` where `q1 < q2`. For CX this means q1 is
always the control and q2 the target. `bench_quest_2q_at()` must pass `(rho, q1, q2)`
in the same order. Implementers should verify that `applySwap`'s argument order is
also `(Qureg, int, int)` with no semantic asymmetry — SWAP is symmetric so the order
does not affect results, but CX is not.

### R3 — Qulacs argument order for 2Q gates (Step E)
`dm_CNOT_gate(UINT control, UINT target, CTYPE*, ITYPE)`. Verify the convention matches
what is used in `bench_qulacs()`. The existing code passes `(q1s[p], q2s[p], rho, dim)`
with no documentation of which is control vs target. The test result is correct only if
the convention is consistent across all methods (all treat q1 as control for CX).

### R4 — State sharing across runs in `_at` functions (Step B)
The `_at` functions initialize one state, run all warmup calls, then run all `runs`
timed rounds on that same state. This mirrors the existing `run_1q_timing()` pattern.
However, for a fixed target, the state is repeatedly transformed by the same gate,
potentially driving it toward an eigenstate (e.g., repeated X on a qubit in |0⟩ or
|1⟩). Since the state is initialized with random data (not a valid density matrix),
this is unlikely to be a problem in practice, but it is worth noting: results for
very small `iterations` (e.g., 1) with a fixed target may differ from results with a
sweep because the state is not re-randomized between runs.

### R5 — Stack size for `perq` array in Step H
`bench_result_perq_t[128][6]` is approximately 120 KB. This is safe on Linux (8 MB
default stack) but macOS sets the stack to 8 MB by default for the main thread and
512 KB for secondary threads. Since `bench_gate` is single-threaded at the
orchestration level, 120 KB is fine. If future changes introduce threading at the
orchestration level, switch to heap allocation.

### R6 — `dim_t` type width
`dim_t` is `int64_t` based on the existing `printf` format `%lld` in `bench_print_csv`.
The `bench_result_perq_t` struct uses `dim_t` for the `dim` field and `bench_print_perq_csv`
must use `%lld` (with the `(long long)` cast already used in `bench_print_csv`).
Verify this matches before Step G is merged.
