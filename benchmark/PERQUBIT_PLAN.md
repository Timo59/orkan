# Per-Qubit Benchmark Extension — Implementation Plan

Target binary: `bench_gate` (built from `benchmark/src/`)
New flag: `--per-qubit`

This plan is written for parallel execution. Steps with no shared dependencies can
be handed to independent agents or developers simultaneously. Each step is sized at
roughly 1–3 hours of focused implementation work.

---

## Measurement Semantics

**What `--per-qubit` measures**: steady-state single-target throughput.

Each `_at` function applies the same gate to the same fixed target qubit thousands of
times on a single state that is never re-initialized between runs. By the time the
timed loop starts (after warmup), the branch predictor has saturated, the hardware
prefetcher has learned the access stride, and all working-set pages are resident in
TLB. This is **not** circuit-representative throughput. It is a deliberate isolation
measurement: it removes target-sweep averaging to expose whether gate cost depends on
target position within the memory layout of each representation, and it measures the
steady-state ceiling each implementation can sustain.

Consequences:
- Results are **not** directly comparable to the sweep benchmark's `time_per_gate_ms`.
  The sweep benchmark is a better proxy for circuit-like workloads.
- Results **are** appropriate for comparing representations head-to-head at each
  target position, for identifying memory-layout asymmetries, and for bounding
  achievable throughput.
- Document this distinction whenever citing per-qubit numbers in a thesis or paper.

Note: runs within a single `_at` call are NOT statistically independent — the state is not re-initialized between runs. The reported standard deviation and CV are measures of timing variability (scheduler jitter, thermal effects), not statistical sampling error. Confidence intervals derived from these values are approximate.

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
         Step E2 ─────────┤   ← Aer-DM adapter; depends only on A; parallel with B–E
       bench_aer_dm.cpp   │
                          │
         Step F ──────────┤   ← depends only on A; can run in parallel with B–E2
       (CLI parsing)      │
                          │
         Step G ──────────┤   ← depends only on A; can run in parallel with B–F
       (output fns)       │
                          │
                          ▼
               ┌─────────────────────┐
               │   Step H: main()    │  ← depends on B, C, D, E, E2, F, G (all prior)
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
- Batch 1 (all independent, all depend only on A): Steps B, C, D, E, E2, F, G
- Batch 2 (integration): Step H (depends on all of Batch 1)
- Batch 3 (docs): Step I (depends on A, F, G, H)
- Batch 4 (thesis): Step J (depends on H, I)

---

## Step A — `bench.h` additions

**Depends on:** nothing
**File:** `benchmark/include/bench.h`
**Estimated time:** 1.5 h

### What to add

**1. Per-qubit default constants** — insert after the existing `BENCH_MAX_RUNS` line:

```c
#define BENCH_PQ_DEFAULT_RUNS        5    /**< Default --runs when --per-qubit is active.
                                           *   5 runs is sufficient for steady-state throughput
                                           *   measurements — variance reflects scheduler jitter,
                                           *   not sampling error. Increase with --runs for
                                           *   publication-quality data, but expect runtime to
                                           *   scale linearly. */
#define BENCH_PQ_DEFAULT_ITERATIONS  100  /**< Default --iterations when --per-qubit is active */
/**
 * Hard cap on --runs for all _at backends. Must equal BENCH_MAX_RUNS because
 * bench_compute_stats() sorts into a stack buffer of size 200. If more runs are
 * needed in the future, bench_compute_stats must be updated to heap-allocate its
 * sort buffer before this constant is raised.
 */
#define BENCH_PQ_MAX_RUNS                  BENCH_MAX_RUNS  /* = 200 */
#define BENCH_PQ_FALLBACK_ITERATIONS       10000 /**< Fallback when probe duration rounds to 0 ns */
/**
 * Hard cap on auto-calibrated iteration count. Without this cap, a nanosecond-scale
 * gate at 2 qubits yields probe_ms ≈ 1e-9, causing ceil(1.0/1e-9) = 1,000,000,000
 * iterations — a multi-hour run. Cap at 1,000,000 to keep any single _at call under
 * ~1 second at the smallest realistic gate time.
 */
#define BENCH_PQ_MAX_CALIBRATED_ITERATIONS 1000000
/**
 * Number of iterations in the auto-calibration probe run. Using 100 iterations
 * rather than 1 ensures the probe measures gate cost, not clock_gettime overhead
 * (~20-50 ns on Linux/macOS). At 2-4 qubits, a single gate call is 5-50 ns;
 * timing 1 call would measure the timer, not the gate.
 */
#define BENCH_PQ_PROBE_ITERS               100
/** Minimum timed run duration in ms for auto-calibration.
 *  All _at functions target at least this many ms per run when computing
 *  calibrated_iterations via ceil(MIN_MEASUREMENT_MS / probe_ms).
 *  Step H uses the same constant when calibrating outside the _at functions. */
#define MIN_MEASUREMENT_MS 1.0
/**
 * Maximum qubit count for build_all_ordered_pairs. At N=16, N*(N-1)=240 < 2*MAX_PAIRS=256.
 * At N=17, N*(N-1)=272 > 256, which would overflow the output arrays.
 * Enforced by a runtime guard inside build_all_ordered_pairs.
 */
#define MAX_QUBITS_FOR_ORDERED_PAIRS       16
```

Note: `BENCH_PQ_DEFAULT_RUNS` is 5. 5 runs is sufficient for steady-state throughput
measurements — variance reflects scheduler jitter, not sampling error. Increase with
--runs for publication-quality data, but expect runtime to scale linearly.

Note on `BENCH_PQ_MAX_RUNS`: `bench_compute_stats()` in `bench.h` contains
`double sorted[200]` with `cnt = (n < 200) ? n : 200`. If `BENCH_PQ_MAX_RUNS` exceeded
200, the median would be computed over a truncated 200-sample population while mean and
stddev would use all samples — a silent statistical inconsistency. `BENCH_PQ_MAX_RUNS`
is therefore capped at `BENCH_MAX_RUNS`. To support more runs, heap-allocate the sort
buffer in `bench_compute_stats` first.

**2. Compiler barrier macro** — insert after the constants block, before the Types section:

```c
/**
 * @brief Prevent the compiler from reordering or eliding the preceding store/load
 *        across this point. Used to ensure probe and warmup calls are not hoisted
 *        out of timing loops by the optimizer.
 *
 * On GCC/Clang this emits an empty inline assembly clobber. On other compilers it
 * is a no-op (relying on the volatile clock call to bound reordering in practice).
 */
#if defined(__GNUC__) || defined(__clang__)
#define BENCH_COMPILER_BARRIER() __asm__ volatile("" ::: "memory")
#else
#define BENCH_COMPILER_BARRIER() /* no-op on unsupported compilers */
#endif
```

**3. Sentinel constant** — insert near the top of the Types section, before the structs:

```c
/** @brief Sentinel value for bench_result_perq_t.target2 when the gate is single-qubit. */
#define QUBIT_NONE ((qubit_t)0xFF)

/* Compile-time guard: QUBIT_NONE must exceed the maximum supported qubit count (30).
 * Use static_assert (C11/C++11) rather than _Static_assert (C11-only) so this
 * compiles cleanly in bench_qulacs.cpp and bench_aer_dm.cpp as well. */
static_assert((qubit_t)0xFF > 30,
    "QUBIT_NONE sentinel must exceed max supported qubit count");
```

Note: `qubit_t` is `uint8_t` (confirmed via `q_types.h`). `0xFF` is safe as a sentinel
because no valid qubit index reaches 255; the maximum supported qubit count is 30. The
`static_assert` documents this invariant in the code without runtime cost.

**4. `bench_result_perq_t` struct** — insert after `bench_result_t`:

```c
/**
 * @brief Per-qubit benchmark result: timing for a single (gate, qubits, method, target) tuple.
 *
 * Parallel to bench_result_t but records a fixed target qubit rather than sweeping.
 * For 2Q gates, target = q1 and target2 = q2 (q1 < q2). For 1Q gates, target2 = QUBIT_NONE.
 *
 * MEASUREMENT SEMANTICS: This struct records steady-state single-target throughput.
 * The gate is applied to the same target on the same (non-reinitialized) state for
 * all runs. Branch predictor saturation and prefetcher warm-up are intentional; this
 * isolates memory-layout effects from circuit-representative costs.
 */
typedef struct bench_result_perq {
    qubit_t    qubits;           /**< Number of qubits in the state */
    dim_t      dim;              /**< Hilbert space dimension (2^qubits) */
    const char *gate_name;       /**< Name of the gate tested */
    const char *method;          /**< Implementation method name */
    qubit_t    target;           /**< Target qubit (q1 for 2Q gates) */
    qubit_t    target2;          /**< Second qubit for 2Q gates; QUBIT_NONE for 1Q */
    int        is_2q;            /**< 1 if this is a 2Q measurement, 0 otherwise */
    double     time_ms;               /**< Mean total time per run (ms) */
    double     time_ms_std;           /**< Sample std dev across runs, Bessel-corrected (ms) */
    double     time_ms_min;           /**< Minimum run time (ms) */
    double     time_ms_median;        /**< Median run time (ms) */
    double     time_ms_cv;            /**< Coefficient of variation = std/mean x 100 (%) */
    double     ops_per_sec;           /**< Gate applications per second, from mean time */
    double     time_per_gate_ms;      /**< Mean time per gate application: time_ms / iterations */
    double     time_per_gate_ms_median; /**< Median time per gate application: time_ms_median / iterations */
    double     time_per_gate_ms_min;    /**< Min time per gate application: time_ms_min / iterations */
    size_t     memory_bytes;     /**< Memory used for state representation */
    int        iterations;       /**< Gate calls per timed run */
    int        runs;             /**< Number of independent timing runs */
    /**
     * @brief Theoretical peak memory bandwidth (GB/s), assuming full-matrix read+write.
     *
     * Formula: 2 * memory_bytes / (time_per_gate_ms * 1e-3) / 1e9
     *        = 2 * memory_bytes / (time_per_gate_ms * 1e6)
     *
     * THIS IS AN UPPER-BOUND ESTIMATE, not an exact measurement. Caveats vary by
     * backend — see bw_notes field for a per-call annotation:
     * - "half_matrix": qlib_packed / qlib_tiled, 1Q gates — only ~half the matrix
     *   elements are touched per single-qubit gate. Formula OVERSTATES consumption by up to 2x.
     * - "three_quarter_matrix": qlib_packed / qlib_tiled, 2Q gates — approximately 3/4
     *   of matrix elements are touched (all elements where at least one of the two target
     *   bits differs). Formula OVERSTATES consumption by ~4/3.
     * - "full_matrix": an approximate full read+write pass per gate.
     * - "multi_pass_est": blas_dense — two zgemm calls + temp buffer ≈ 5x memory_bytes
     *   of traffic. Formula UNDERSTATES consumption.
     * - "unknown": quest / qulacs / aer_dm — internal access patterns are opaque.
     *
     * Use theoretical_bw_gbs to compare methods RELATIVE to each other and to the
     * hardware STREAM bandwidth ceiling. Do not treat it as a precise physical
     * measurement. Filter or scale by bw_notes in downstream analysis.
     */
    double     theoretical_bw_gbs;
    /**
     * @brief Bandwidth annotation for downstream filtering (see theoretical_bw_gbs).
     *
     * One of: "half_matrix", "three_quarter_matrix", "full_matrix", "multi_pass_est", "unknown".
     * Set by bench_fill_perq_stats via its bw_notes parameter.
     */
    const char *bw_notes;
} bench_result_perq_t;
```

Note: `sweep_size` is omitted from `bench_result_perq_t` because per-qubit results always
have exactly 1 gate application per iteration (sweep_size is implicitly 1). `time_per_gate_ms`
is therefore simply `time_ms / iterations`.

**5. `bench_options_t` additions** — add to the existing struct:

```c
    int per_qubit;         /**< 1 if --per-qubit flag was given */
    int runs_explicit;     /**< 1 if --runs was given on the command line */
    int iterations_explicit; /**< 1 if --iterations was given on the command line */
    const char *gate_filter; /**< Non-NULL if --gate <NAME> was given; NULL means all gates */
```

The `runs_explicit` and `iterations_explicit` flags are required. The heuristic
approach of comparing `opts.runs == BENCH_DEFAULT_RUNS` is incorrect when a user
passes `--per-qubit --runs 5` (5 happens to equal `BENCH_PQ_DEFAULT_RUNS`; the
override would fire). Explicit flags cost three lines and are unambiguous.

**6. `bench_fill_perq_stats` helper** — insert before the `_at` adapter declarations:

```c
/**
 * @brief Populate timing statistics fields of a bench_result_perq_t from a stats struct.
 *
 * Eliminates the 8-way copy-paste of the stats-to-result assignment block across
 * bench_mixed.c, bench_baselines.c, bench_quest.c, and bench_qulacs.cpp.
 *
 * @param r          Result struct to populate (other fields must already be set)
 * @param s          Computed statistics
 * @param iterations Gate calls per run (used to compute time_per_gate_ms, ops_per_sec,
 *                   and theoretical_bw_gbs)
 * @param runs       Number of timed runs (stored in r->runs so callers need not set it separately)
 * @param bw_notes   Bandwidth annotation string (e.g. "half_matrix", "multi_pass_est");
 *                   stored directly — must point to a string literal or static storage.
 */
static inline void bench_fill_perq_stats(bench_result_perq_t *r,
                                          const bench_run_stats_t *s,
                                          int iterations,
                                          int runs,
                                          const char *bw_notes) {
    r->time_ms          = s->mean;
    r->time_ms_std      = s->std_dev;
    r->time_ms_min      = s->min;
    r->time_ms_median   = s->median;
    r->time_ms_cv       = s->cv;
    r->iterations       = iterations;
    r->runs             = runs;
    r->ops_per_sec      = (s->mean > 0.0)
        ? (double)iterations / (s->mean / 1000.0) : 0.0;
    r->time_per_gate_ms = (iterations > 0)
        ? s->mean / (double)iterations : 0.0;
    r->time_per_gate_ms_median = (iterations > 0)
        ? s->median / (double)iterations : 0.0;
    r->time_per_gate_ms_min = (iterations > 0)
        ? s->min / (double)iterations : 0.0;
    r->theoretical_bw_gbs = (r->time_per_gate_ms > 0.0)
        ? 2.0 * (double)r->memory_bytes / (r->time_per_gate_ms * 1e6)
        : 0.0;
    r->bw_notes = bw_notes;
}
```

Note: `bench_run_stats_t` and `bench_compute_stats()` are already defined in `bench.h`
(or `bench_stats.c`). This helper is `static inline` so it is usable from both C and
C++ translation units without a separate compilation unit.

**7. `_at` adapter declarations** — add a new section after the existing 2Q benchmark
function declarations, before the Output functions section:

```c
/*
 * =====================================================================================================================
 * Per-qubit adapter functions — fixed target, single measurement
 * Files: bench_mixed.c, bench_baselines.c, bench_quest.c, bench_qulacs.cpp, bench_aer_dm.cpp
 *
 * All _at functions accept a pre-allocated state pointer (NULL = allocate internally).
 * Passing a pre-allocated, pre-initialized state avoids repeated malloc/free/init cycles
 * in the orchestration loop (see Step H). The state must match the backend's expected
 * layout; passing NULL is always safe but incurs per-call allocation cost.
 * =====================================================================================================================
 */

/* Single-qubit, fixed target */
bench_result_perq_t bench_qlib_packed_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t),
    qubit_t target, int iterations, int warmup, int runs,
    state_t *preinit_state);  /* NULL = allocate internally */

bench_result_perq_t bench_qlib_tiled_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t),
    qubit_t target, int iterations, int warmup, int runs,
    state_t *preinit_state);  /* NULL = allocate internally */

bench_result_perq_t bench_blas_dense_at(qubit_t qubits, const char *gate_name,
    const cplx_t gate_mat[4],
    qubit_t target, int iterations, int warmup, int runs,
    cplx_t *preinit_rho);  /* NULL = allocate internally */

/* Two-qubit, fixed pair */
bench_result_perq_t bench_qlib_packed_2q_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t, qubit_t),
    qubit_t q1, qubit_t q2, int iterations, int warmup, int runs,
    state_t *preinit_state);  /* NULL = allocate internally */

bench_result_perq_t bench_qlib_tiled_2q_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t, qubit_t),
    qubit_t q1, qubit_t q2, int iterations, int warmup, int runs,
    state_t *preinit_state);  /* NULL = allocate internally */

bench_result_perq_t bench_blas_dense_2q_at(qubit_t qubits, const char *gate_name,
    const cplx_t gate_mat[16],
    qubit_t q1, qubit_t q2, int iterations, int warmup, int runs,
    cplx_t *preinit_rho);  /* NULL = allocate internally */

#ifdef WITH_QUEST
bench_result_perq_t bench_quest_at(qubit_t qubits, const char *gate_name,
    qubit_t target, int iterations, int warmup, int runs,
    Qureg *preinit_qureg);  /* NULL = allocate internally; non-NULL = use directly, skip create/destroy */

bench_result_perq_t bench_quest_2q_at(qubit_t qubits, const char *gate_name,
    qubit_t q1, qubit_t q2, int iterations, int warmup, int runs,
    Qureg *preinit_qureg);  /* NULL = allocate internally */
#endif

#ifdef WITH_QULACS
bench_result_perq_t bench_qulacs_at(qubit_t qubits, const char *gate_name,
    qubit_t target, int iterations, int warmup, int runs,
    void *preinit_rho);  /* NULL = allocate internally; CTYPE* cast internally */

bench_result_perq_t bench_qulacs_2q_at(qubit_t qubits, const char *gate_name,
    qubit_t q1, qubit_t q2, int iterations, int warmup, int runs,
    void *preinit_rho);  /* NULL = allocate internally */
#endif

#ifdef WITH_AER_DM
bench_result_perq_t bench_aer_dm_at(qubit_t qubits, const char *gate_name,
    qubit_t target, int iterations, int warmup, int runs,
    void *preinit_dm);  /* NULL = allocate internally; AER::QV::DensityMatrix<double>* cast internally */

bench_result_perq_t bench_aer_dm_2q_at(qubit_t qubits, const char *gate_name,
    qubit_t q1, qubit_t q2, int iterations, int warmup, int runs,
    void *preinit_dm);  /* NULL = allocate internally */
#endif
```

**8. `build_all_ordered_pairs` declaration** — add to `bench.h` immediately after the
`build_all_pairs` declaration:

```c
/**
 * @brief Build all ordered qubit pairs (q1 != q2, both orderings) for an n-qubit system.
 *
 * Writes N*(N-1) pairs into q1_out/q2_out. Unlike build_all_pairs (which only writes
 * pairs with q1 < q2), this function writes both (q1, q2) and (q2, q1), exposing
 * memory-layout asymmetry for gates that are not symmetric (e.g. CX).
 *
 * Callers must ensure q1_out/q2_out have room for at least 2*MAX_PAIRS entries.
 * Returns the number of pairs written (always qubits*(qubits-1) for valid input).
 */
int build_all_ordered_pairs(qubit_t qubits, qubit_t *q1_out, qubit_t *q2_out);
```

**9. Per-qubit output function declarations** — add to the Output functions section:

```c
void bench_print_perq_csv_header(void);
void bench_print_perq_csv(const bench_result_perq_t *r);
void bench_print_perq_console(const bench_result_perq_t *results, int n);
```

### Constraints

- Do NOT modify any `.c` or `.cpp` files in this step.
- Do NOT change any existing declarations or definitions.
- The `bench_options_t` zero-initializer in `bench_parse_options()` (bench_main.c) will
  pick up the new fields as 0/NULL automatically if designated initializers are used
  with all fields named; verify this in Step F.

### Iteration auto-calibration

At small qubit counts (2–6 qubits), a single gate application takes nanoseconds and
`clock_gettime` overhead (20–50 ns on Linux/macOS) is a non-trivial fraction of the
measured time. `BENCH_PQ_DEFAULT_ITERATIONS = 100` provides reasonable coverage at
medium qubit counts but may be insufficient at small counts.

**Sentinel convention**: Pass `iterations=0` to request auto-calibration inside an `_at`
function. Pass a positive value to use that exact count and skip internal calibration.
This eliminates the double-calibration conflict: Step H calibrates once per
(qubits, method) on target 0 and passes the result to all `_at` calls for that method
(for all backends including quest/qulacs/aer_dm). The auto-calibration probe inside
each `_at` function **only fires when `iterations == 0` is passed**.

The `_at` functions in Steps B–E2 **must** auto-scale `iterations` using a probe run
before warmup. The probe pattern and zero-safe formula are:

```c
/* --- Auto-calibration probe (before warmup) ---
 * Convention: only runs when iterations == 0 (sentinel meaning "please auto-calibrate").
 * When iterations > 0, that exact count is used directly and the probe is skipped.
 * Step H calibrates once per (qubits, method) for all backends and passes a positive
 * calibrated_iterations to every _at call. Pass iterations=0 only if calling an _at
 * function directly without an external calibration step.
 * MIN_MEASUREMENT_MS is defined in the constants block above (bench.h). */

if (iterations == 0) {
    /* Probe: run BENCH_PQ_PROBE_ITERS iterations to get a reliable estimate.
     * At 2-4 qubits, a single gate call is 5-50 ns while clock_gettime costs
     * 20-50 ns — timing 1 call measures the timer, not the gate. */
    uint64_t probe_t0 = bench_time_ns();
    for (int pi = 0; pi < BENCH_PQ_PROBE_ITERS; ++pi) {
        gate_fn(&state, target);
        BENCH_COMPILER_BARRIER();
    }
    uint64_t probe_t1 = bench_time_ns();
    double probe_ms = (double)(probe_t1 - probe_t0) / 1e6 / BENCH_PQ_PROBE_ITERS;

    if (probe_ms <= 0.0) {
        /* clock_gettime rounded to 0 ns; use fallback — avoids 1.0/0.0 = +inf and
         * the subsequent (int)+inf which is undefined behaviour in C */
        iterations = BENCH_PQ_FALLBACK_ITERATIONS;
    } else {
        iterations = (int)ceil(MIN_MEASUREMENT_MS / probe_ms);
    }
    /* Cap: a nanosecond-scale gate at 2 qubits would otherwise produce ~1e9 iterations. */
    if (iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
        iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
}
result.iterations = iterations;   /* update result so callers see actual value */
```

Document the final `iterations` value in the result struct; the caller must be able to
see what was actually used.

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

**State initialization**: use `init_random_mixed_state(&state, qubits, desc->state_type)`
with a fixed seed (`srand(42)` before allocation) for within-backend reproducibility.
Cross-backend state content differs — see Known Limitations L1.

### What to add

**Note on `build_all_ordered_pairs` placement**: `build_all_pairs` is defined in
`bench_mixed.c` (line 43) and declared in `bench.h`. Add `build_all_ordered_pairs`
immediately after `build_all_pairs` in `bench_mixed.c`:

```c
/**
 * @brief Build all ordered qubit pairs (q1 != q2, both orderings).
 *
 * Writes N*(N-1) pairs into q1_out/q2_out (both (i,j) and (j,i) for i != j).
 * Callers must ensure q1_out/q2_out have room for at least 2*MAX_PAIRS entries.
 * Returns the number of pairs written, or 0 on error (qubits exceeds
 * MAX_QUBITS_FOR_ORDERED_PAIRS).
 *
 * The runtime guard replaces a previous assert() which was compiled away with
 * -DNDEBUG and could silently overflow the output arrays in release builds.
 */
int build_all_ordered_pairs(qubit_t qubits, qubit_t *q1_out, qubit_t *q2_out) {
    if (qubits > MAX_QUBITS_FOR_ORDERED_PAIRS) {
        fprintf(stderr, "build_all_ordered_pairs: qubits=%d exceeds max (%d)\n",
                (int)qubits, MAX_QUBITS_FOR_ORDERED_PAIRS);
        return 0;
    }
    int n = 0;
    for (qubit_t i = 0; i < qubits; ++i) {
        for (qubit_t j = 0; j < qubits; ++j) {
            if (i == j) continue;
            q1_out[n] = i;
            q2_out[n] = j;
            ++n;
        }
    }
    return n;
}
```

The declaration in `bench.h` (Step A item 8) makes it visible to all callers.

Append to the end of `bench_mixed.c`, after the existing 2Q implementations.

**1. `run_1q_timing_at()` internal function:**

```c
/**
 * @brief Shared timing loop for single-qubit per-qubit benchmarks.
 *
 * Times gate_fn applied to a single fixed target qubit for `iterations` calls
 * per run, across `runs` independent runs. State is initialised once before
 * warmup and shared across all runs (same as run_1q_timing).
 *
 * MEASUREMENT SEMANTICS: steady-state single-target throughput; see plan §Measurement Semantics.
 */
static bench_result_perq_t run_1q_timing_at(qubit_t qubits, const char *gate_name,
                                             void (*gate_fn)(state_t*, qubit_t),
                                             qubit_t target,
                                             int iterations, int warmup, int runs,
                                             const qlib_bench_desc_t *desc,
                                             state_t *preinit_state) {
    bench_result_perq_t result = {0};
    if (gate_fn == NULL) {
        fprintf(stderr, "run_1q_timing_at: NULL gate_fn for gate '%s'\n", gate_name);
        return result;
    }
    if (target >= qubits) {
        fprintf(stderr, "run_1q_timing_at: target=%u >= qubits=%u\n",
                (unsigned)target, (unsigned)qubits);
        return result;
    }
    if (runs <= 0) {
        fprintf(stderr, "run_1q_timing_at: runs=%d must be > 0\n", runs);
        return result;
    }
    if (runs > BENCH_PQ_MAX_RUNS) {
        fprintf(stderr, "run_1q_timing_at: runs=%d exceeds BENCH_PQ_MAX_RUNS=%d\n",
                runs, BENCH_PQ_MAX_RUNS);
        return result;
    }
    result.qubits       = qubits;
    result.dim          = (dim_t)1 << qubits;
    result.gate_name    = gate_name;
    result.method       = desc->method;
    result.target       = target;
    result.target2      = QUBIT_NONE;
    result.is_2q        = 0;
    result.memory_bytes = desc->memory_bytes;
    result.runs         = runs;

    int owns_state = (preinit_state == NULL);
    state_t local_state = {0};
    state_t *state = preinit_state ? preinit_state : &local_state;

    if (owns_state) {
        srand(42);  /* fixed seed for within-backend reproducibility only;
                     * cross-backend state content differs (see Known Limitations L1) */
        init_random_mixed_state(state, qubits, desc->state_type);
        if (!state->data) return result;
    }

    /* Auto-calibration probe — only when iterations == 0 (sentinel).
     * Pass iterations=0 to auto-calibrate; pass a positive value to use it directly.
     * See Step A §Iteration auto-calibration. */
    if (iterations == 0) {
        uint64_t probe_t0 = bench_time_ns();
        for (int pi = 0; pi < BENCH_PQ_PROBE_ITERS; ++pi) {
            gate_fn(state, target);
            BENCH_COMPILER_BARRIER();
        }
        uint64_t probe_t1 = bench_time_ns();
        double probe_ms = (double)(probe_t1 - probe_t0) / 1e6 / BENCH_PQ_PROBE_ITERS;
        if (probe_ms <= 0.0) {
            iterations = BENCH_PQ_FALLBACK_ITERATIONS;
        } else {
            iterations = (int)ceil(MIN_MEASUREMENT_MS / probe_ms);
        }
        if (iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
            iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
    }
    result.iterations = iterations;

    for (int i = 0; i < warmup; ++i) {
        gate_fn(state, target);
        BENCH_COMPILER_BARRIER();
    }

    double *run_times = malloc((size_t)runs * sizeof(double));
    if (!run_times) {
        if (owns_state) state_free(state);
        return result;
    }
    for (int r = 0; r < runs; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i) {
            gate_fn(state, target);
        }
        BENCH_COMPILER_BARRIER();  /* once per run: prevents loop hoisting without per-call overhead */
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }
    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    free(run_times);

    bench_fill_perq_stats(&result, &stats, iterations, runs, "half_matrix");
    if (owns_state) state_free(state);
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
 *
 * MEASUREMENT SEMANTICS: steady-state single-target throughput; see plan §Measurement Semantics.
 */
static bench_result_perq_t run_2q_timing_at(qubit_t qubits, const char *gate_name,
                                             void (*gate_fn)(state_t*, qubit_t, qubit_t),
                                             qubit_t q1, qubit_t q2,
                                             int iterations, int warmup, int runs,
                                             const qlib_bench_desc_t *desc,
                                             state_t *preinit_state) {
    bench_result_perq_t result = {0};
    if (gate_fn == NULL) {
        fprintf(stderr, "run_2q_timing_at: NULL gate_fn for gate '%s'\n", gate_name);
        return result;
    }
    if (q1 >= qubits || q2 >= qubits) {
        fprintf(stderr, "run_2q_timing_at: q1=%u or q2=%u >= qubits=%u\n",
                (unsigned)q1, (unsigned)q2, (unsigned)qubits);
        return result;
    }
    if (q1 == q2) {
        fprintf(stderr, "run_2q_timing_at: q1 == q2 (%u)\n", (unsigned)q1);
        return result;
    }
    if (runs <= 0) {
        fprintf(stderr, "run_2q_timing_at: runs=%d must be > 0\n", runs);
        return result;
    }
    if (runs > BENCH_PQ_MAX_RUNS) {
        fprintf(stderr, "run_2q_timing_at: runs=%d exceeds BENCH_PQ_MAX_RUNS=%d\n",
                runs, BENCH_PQ_MAX_RUNS);
        return result;
    }
    result.qubits       = qubits;
    result.dim          = (dim_t)1 << qubits;
    result.gate_name    = gate_name;
    result.method       = desc->method;
    result.target       = q1;
    result.target2      = q2;
    result.is_2q        = 1;
    result.memory_bytes = desc->memory_bytes;
    result.runs         = runs;

    int owns_state = (preinit_state == NULL);
    state_t local_state = {0};
    state_t *state = preinit_state ? preinit_state : &local_state;

    if (owns_state) {
        srand(42);  /* fixed seed for within-backend reproducibility only;
                     * cross-backend state content differs (see Known Limitations L1) */
        init_random_mixed_state(state, qubits, desc->state_type);
        if (!state->data) return result;
    }

    /* Auto-calibration probe — only when iterations == 0 (sentinel).
     * Pass iterations=0 to auto-calibrate; pass a positive value to use it directly. */
    if (iterations == 0) {
        uint64_t probe_t0 = bench_time_ns();
        for (int pi = 0; pi < BENCH_PQ_PROBE_ITERS; ++pi) {
            gate_fn(state, q1, q2);
            BENCH_COMPILER_BARRIER();
        }
        uint64_t probe_t1 = bench_time_ns();
        double probe_ms = (double)(probe_t1 - probe_t0) / 1e6 / BENCH_PQ_PROBE_ITERS;
        if (probe_ms <= 0.0) {
            iterations = BENCH_PQ_FALLBACK_ITERATIONS;
        } else {
            iterations = (int)ceil(MIN_MEASUREMENT_MS / probe_ms);
        }
        if (iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
            iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
    }
    result.iterations = iterations;

    for (int i = 0; i < warmup; ++i) {
        gate_fn(state, q1, q2);
        BENCH_COMPILER_BARRIER();
    }

    double *run_times = malloc((size_t)runs * sizeof(double));
    if (!run_times) {
        if (owns_state) state_free(state);
        return result;
    }
    for (int r = 0; r < runs; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i) {
            gate_fn(state, q1, q2);
        }
        BENCH_COMPILER_BARRIER();  /* once per run: prevents loop hoisting without per-call overhead */
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }
    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    free(run_times);

    bench_fill_perq_stats(&result, &stats, iterations, runs, "three_quarter_matrix");
    if (owns_state) state_free(state);
    return result;
}
```

**3. Four public wrapper functions** (thin, mirroring the existing packed/tiled wrappers):

```c
bench_result_perq_t bench_qlib_packed_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t),
    qubit_t target, int iterations, int warmup, int runs,
    state_t *preinit_state)
{
    const qlib_bench_desc_t desc = {
        .state_type   = MIXED_PACKED,
        .memory_bytes = bench_packed_size(qubits),
        .method       = "qlib_packed",
    };
    return run_1q_timing_at(qubits, gate_name, gate_fn, target,
                            iterations, warmup, runs, &desc, preinit_state);
}

bench_result_perq_t bench_qlib_tiled_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t),
    qubit_t target, int iterations, int warmup, int runs,
    state_t *preinit_state)
{
    const qlib_bench_desc_t desc = {
        .state_type   = MIXED_TILED,
        .memory_bytes = bench_tiled_size(qubits),
        .method       = "qlib_tiled",
    };
    return run_1q_timing_at(qubits, gate_name, gate_fn, target,
                            iterations, warmup, runs, &desc, preinit_state);
}

bench_result_perq_t bench_qlib_packed_2q_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t, qubit_t),
    qubit_t q1, qubit_t q2, int iterations, int warmup, int runs,
    state_t *preinit_state)
{
    const qlib_bench_desc_t desc = {
        .state_type   = MIXED_PACKED,
        .memory_bytes = bench_packed_size(qubits),
        .method       = "qlib_packed",
    };
    return run_2q_timing_at(qubits, gate_name, gate_fn, q1, q2,
                            iterations, warmup, runs, &desc, preinit_state);
}

bench_result_perq_t bench_qlib_tiled_2q_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t, qubit_t),
    qubit_t q1, qubit_t q2, int iterations, int warmup, int runs,
    state_t *preinit_state)
{
    const qlib_bench_desc_t desc = {
        .state_type   = MIXED_TILED,
        .memory_bytes = bench_tiled_size(qubits),
        .method       = "qlib_tiled",
    };
    return run_2q_timing_at(qubits, gate_name, gate_fn, q1, q2,
                            iterations, warmup, runs, &desc, preinit_state);
}
```

### Constraints

- Do NOT modify `run_1q_timing()` or `run_2q_timing()` or any existing public functions.
- The `qlib_bench_desc_t` type is already defined in this file (file-local struct); the
  new static functions can use it directly.
- Use `bench_fill_perq_stats()` (defined in `bench.h` as `static inline`) for all
  stats-to-result assignment, passing `runs` as the fourth argument; do NOT hand-write
  the block again.
- `bw_notes` for qlib backends is `"half_matrix"` for 1Q gates (roughly half the
  matrix elements touched) and `"three_quarter_matrix"` for 2Q gates (approximately
  3/4 of elements touched — all where at least one target bit differs).
- **Inlining and LTO**: `gate_fn` is called via a function pointer. It must NOT be
  marked `__attribute__((always_inline))`, and LTO must not devirtualize across the
  pointer call site. If the compiler can see through the pointer and inline the gate
  body into the timing loop, it may perform loop-wide optimizations (CSE, hoisting,
  vectorization across iterations) that eliminate realistic gate cost, producing
  artificially low timings. Ensure `gate_fn` remains an opaque function pointer at
  the call site; use `BENCH_COMPILER_BARRIER()` between calls if additional protection
  is needed.

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

**BLAS bandwidth note**: `apply_uruh_blas` computes `U * rho * U^dagger` via two
`zgemm` calls with an intermediate temporary buffer. The minimum memory traffic is
approximately 5 × `memory_bytes` (read U, read rho, write tmp, read U^dagger, read tmp,
write rho). The `theoretical_bw_gbs` field in the result (computed by
`bench_fill_perq_stats` as `2 * memory_bytes / time`) severely understates actual
bandwidth consumption for this backend. The `bw_notes` field is set to `"multi_pass_est"`
to flag this. The CSV column is labeled `theoretical_bw_gbs` specifically to indicate
that the formula does not apply uniformly across all backends.

**State initialization**: use `srand(42)` before initializing `rho` for within-backend
reproducibility. Cross-backend state content differs — see Known Limitations L1.
Comment: `/* fixed seed for within-backend reproducibility only; cross-backend state content differs (see L1) */`.

### What to add

Append to the end of `bench_baselines.c`.

**1. `bench_blas_dense_at()`:**

```c
bench_result_perq_t bench_blas_dense_at(qubit_t qubits, const char *gate_name,
                                         const cplx_t gate_mat[4],
                                         qubit_t target,
                                         int iterations, int warmup, int runs,
                                         cplx_t *preinit_rho) {
    bench_result_perq_t result = {0};
    if (target >= qubits) {
        fprintf(stderr, "bench_blas_dense_at: target=%u >= qubits=%u\n",
                (unsigned)target, (unsigned)qubits);
        return result;
    }
    if (runs <= 0) {
        fprintf(stderr, "bench_blas_dense_at: runs=%d must be > 0\n", runs);
        return result;
    }
    if (runs > BENCH_PQ_MAX_RUNS) {
        fprintf(stderr, "bench_blas_dense_at: runs=%d exceeds BENCH_PQ_MAX_RUNS=%d\n",
                runs, BENCH_PQ_MAX_RUNS);
        return result;
    }
    result.qubits       = qubits;
    result.dim          = (dim_t)1 << qubits;
    result.gate_name    = gate_name;
    result.method       = "blas_dense";
    result.target       = target;
    result.target2      = QUBIT_NONE;
    result.is_2q        = 0;
    result.memory_bytes = bench_dense_size(qubits);
    result.runs         = runs;

    dim_t dim = result.dim;

    int owns_rho = (preinit_rho == NULL);
    cplx_t *rho = preinit_rho;
    if (owns_rho) {
        rho = malloc((size_t)dim * (size_t)dim * sizeof(cplx_t));
        if (!rho) return result;
        srand(42);  /* fixed seed for within-backend reproducibility only;
                     * cross-backend state content differs (see Known Limitations L1) */
        for (dim_t i = 0; i < dim * dim; ++i)
            rho[i] = (double)rand() / RAND_MAX + ((double)rand() / RAND_MAX) * I;
    }

    /* Build gate matrix only for the fixed target — outside the timed loop */
    cplx_t *U = build_full_gate_matrix(qubits, gate_mat, target);
    if (!U) { if (owns_rho) free(rho); return result; }

    cplx_t *tmp = malloc((size_t)dim * (size_t)dim * sizeof(cplx_t));
    if (!tmp) { free(U); if (owns_rho) free(rho); return result; }

    /* Auto-calibration probe — only when iterations == 0 (sentinel).
     * Pass iterations=0 to auto-calibrate; pass a positive value to use it directly. */
    if (iterations == 0) {
        uint64_t probe_t0 = bench_time_ns();
        for (int pi = 0; pi < BENCH_PQ_PROBE_ITERS; ++pi) {
            apply_uruh_blas(dim, U, rho, tmp);
            BENCH_COMPILER_BARRIER();
        }
        uint64_t probe_t1 = bench_time_ns();
        double probe_ms = (double)(probe_t1 - probe_t0) / 1e6 / BENCH_PQ_PROBE_ITERS;
        if (probe_ms <= 0.0) {
            iterations = BENCH_PQ_FALLBACK_ITERATIONS;
        } else {
            iterations = (int)ceil(MIN_MEASUREMENT_MS / probe_ms);
        }
        if (iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
            iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
    }
    result.iterations = iterations;

    for (int i = 0; i < warmup; ++i) {
        apply_uruh_blas(dim, U, rho, tmp);
        BENCH_COMPILER_BARRIER();
    }

    double *run_times = malloc((size_t)runs * sizeof(double));
    if (!run_times) { free(tmp); free(U); if (owns_rho) free(rho); return result; }
    for (int r = 0; r < runs; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i) {
            apply_uruh_blas(dim, U, rho, tmp);
        }
        BENCH_COMPILER_BARRIER();  /* once per run: prevents loop hoisting without per-call overhead */
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }
    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    free(run_times);

    bench_fill_perq_stats(&result, &stats, iterations, runs, "multi_pass_est");
    free(tmp); free(U); if (owns_rho) free(rho);
    return result;
}
```

**2. `bench_blas_dense_2q_at()`:**

```c
bench_result_perq_t bench_blas_dense_2q_at(qubit_t qubits, const char *gate_name,
                                             const cplx_t gate_mat[16],
                                             qubit_t q1, qubit_t q2,
                                             int iterations, int warmup, int runs,
                                             cplx_t *preinit_rho) {
    bench_result_perq_t result = {0};
    if (q1 >= qubits || q2 >= qubits) {
        fprintf(stderr, "bench_blas_dense_2q_at: q1=%u or q2=%u >= qubits=%u\n",
                (unsigned)q1, (unsigned)q2, (unsigned)qubits);
        return result;
    }
    if (q1 == q2) {
        fprintf(stderr, "bench_blas_dense_2q_at: q1 == q2 (%u)\n", (unsigned)q1);
        return result;
    }
    if (runs <= 0) {
        fprintf(stderr, "bench_blas_dense_2q_at: runs=%d must be > 0\n", runs);
        return result;
    }
    if (runs > BENCH_PQ_MAX_RUNS) {
        fprintf(stderr, "bench_blas_dense_2q_at: runs=%d exceeds BENCH_PQ_MAX_RUNS=%d\n",
                runs, BENCH_PQ_MAX_RUNS);
        return result;
    }
    result.qubits       = qubits;
    result.dim          = (dim_t)1 << qubits;
    result.gate_name    = gate_name;
    result.method       = "blas_dense";
    result.target       = q1;
    result.target2      = q2;
    result.is_2q        = 1;
    result.memory_bytes = bench_dense_size(qubits);
    result.runs         = runs;

    dim_t dim = result.dim;

    int owns_rho = (preinit_rho == NULL);
    cplx_t *rho = preinit_rho;
    if (owns_rho) {
        rho = malloc((size_t)dim * (size_t)dim * sizeof(cplx_t));
        if (!rho) return result;
        srand(42);  /* fixed seed for within-backend reproducibility only;
                     * cross-backend state content differs (see Known Limitations L1) */
        for (dim_t i = 0; i < dim * dim; ++i)
            rho[i] = (double)rand() / RAND_MAX + ((double)rand() / RAND_MAX) * I;
    }

    /* Build gate matrix only for the fixed pair — outside the timed loop */
    cplx_t *U = build_full_gate_matrix_2q(qubits, gate_mat, q1, q2);
    if (!U) { if (owns_rho) free(rho); return result; }

    cplx_t *tmp = malloc((size_t)dim * (size_t)dim * sizeof(cplx_t));
    if (!tmp) { free(U); if (owns_rho) free(rho); return result; }

    /* Auto-calibration probe — only when iterations == 0 (sentinel).
     * Pass iterations=0 to auto-calibrate; pass a positive value to use it directly. */
    if (iterations == 0) {
        uint64_t probe_t0 = bench_time_ns();
        for (int pi = 0; pi < BENCH_PQ_PROBE_ITERS; ++pi) {
            apply_uruh_blas(dim, U, rho, tmp);
            BENCH_COMPILER_BARRIER();
        }
        uint64_t probe_t1 = bench_time_ns();
        double probe_ms = (double)(probe_t1 - probe_t0) / 1e6 / BENCH_PQ_PROBE_ITERS;
        if (probe_ms <= 0.0) {
            iterations = BENCH_PQ_FALLBACK_ITERATIONS;
        } else {
            iterations = (int)ceil(MIN_MEASUREMENT_MS / probe_ms);
        }
        if (iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
            iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
    }
    result.iterations = iterations;

    for (int i = 0; i < warmup; ++i) {
        apply_uruh_blas(dim, U, rho, tmp);
        BENCH_COMPILER_BARRIER();
    }

    double *run_times = malloc((size_t)runs * sizeof(double));
    if (!run_times) { free(tmp); free(U); if (owns_rho) free(rho); return result; }
    for (int r = 0; r < runs; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i) {
            apply_uruh_blas(dim, U, rho, tmp);
        }
        BENCH_COMPILER_BARRIER();  /* once per run: prevents loop hoisting without per-call overhead */
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }
    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    free(run_times);

    bench_fill_perq_stats(&result, &stats, iterations, runs, "multi_pass_est");
    free(tmp); free(U); if (owns_rho) free(rho);
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
- Use `bench_fill_perq_stats()` with `runs` and `bw_notes = "multi_pass_est"` for stats assignment.

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

**Stack allocation fix**: the original plan used `double run_times[200]` with a
`cnt = min(runs, 200)` cap — silently truncating to 200 if `runs > 200`. This is
incorrect: if qlib and blas use 500 runs but QuEST caps at 200, results are not
comparable and no warning is issued. Use heap allocation for all backends. A hard
error if `runs > BENCH_PQ_MAX_RUNS` (200) is appropriate.

**State initialization**: use `initDebugState(rho)` (existing QuEST pattern). QuEST's
`initDebugState` produces a fixed deterministic state and ignores `srand(42)` entirely.
This is an accepted limitation for 1.0 — see Known Limitations L1. Document it with
a comment.

### What to add

Append inside the `#ifdef WITH_QUEST` block, after the existing `bench_quest()` function:

**1. `bench_quest_at()` — single-qubit, fixed target:**

```c
bench_result_perq_t bench_quest_at(qubit_t qubits, const char *gate_name,
                                    qubit_t target, int iterations, int warmup, int runs,
                                    Qureg *preinit_qureg) {
    bench_result_perq_t result = {0};
    if (target >= qubits) {
        fprintf(stderr, "bench_quest_at: target=%u >= qubits=%u\n",
                (unsigned)target, (unsigned)qubits);
        return result;
    }
    if (runs <= 0) {
        fprintf(stderr, "bench_quest_at: runs=%d must be > 0\n", runs);
        return result;
    }
    if (runs > BENCH_PQ_MAX_RUNS) {
        fprintf(stderr, "bench_quest_at: runs=%d exceeds BENCH_PQ_MAX_RUNS=%d\n",
                runs, BENCH_PQ_MAX_RUNS);
        return result;
    }
    result.qubits       = qubits;
    result.dim          = (dim_t)1 << qubits;
    result.gate_name    = gate_name;
    result.method       = "quest";
    result.target       = target;
    result.target2      = QUBIT_NONE;
    result.is_2q        = 0;
    result.runs         = runs;

    dim_t dim = (dim_t)1 << qubits;
    result.memory_bytes = (size_t)dim * (size_t)dim * sizeof(double) * 2;

    int owns_qureg = (preinit_qureg == NULL);
    Qureg local_qureg;
    Qureg rho;
    if (owns_qureg) {
        local_qureg = createDensityQureg((int)qubits);
        /* initDebugState produces a fixed deterministic state and ignores srand(42).
         * Cross-backend state content therefore differs from qlib/BLAS/Qulacs.
         * This is an accepted limitation — see Known Limitations L1. */
        initDebugState(local_qureg);
        rho = local_qureg;
    } else {
        rho = *preinit_qureg;
    }

    void (*gate_fn_1q)(Qureg, int) = NULL;
    if      (strcmp(gate_name, "X") == 0) gate_fn_1q = applyPauliX;
    else if (strcmp(gate_name, "H") == 0) gate_fn_1q = applyHadamard;
    else if (strcmp(gate_name, "Z") == 0) gate_fn_1q = applyPauliZ;

    if (!gate_fn_1q) {
        if (owns_qureg) destroyQureg(rho);
        return result;
    }

    /* Auto-calibration probe — only when iterations == 0 (sentinel).
     * Pass iterations=0 to auto-calibrate; pass a positive value to use it directly. */
    if (iterations == 0) {
        uint64_t probe_t0 = bench_time_ns();
        for (int pi = 0; pi < BENCH_PQ_PROBE_ITERS; ++pi) {
            gate_fn_1q(rho, (int)target);
            BENCH_COMPILER_BARRIER();
        }
        uint64_t probe_t1 = bench_time_ns();
        double probe_ms = (double)(probe_t1 - probe_t0) / 1e6 / BENCH_PQ_PROBE_ITERS;
        if (probe_ms <= 0.0) {
            iterations = BENCH_PQ_FALLBACK_ITERATIONS;
        } else {
            iterations = (int)ceil(MIN_MEASUREMENT_MS / probe_ms);
        }
        if (iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
            iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
    }
    result.iterations = iterations;

    double *run_times = malloc((size_t)runs * sizeof(double));
    if (!run_times) {
        if (owns_qureg) destroyQureg(rho);
        return result;
    }

    for (int i = 0; i < warmup; ++i) {
        gate_fn_1q(rho, (int)target);
        BENCH_COMPILER_BARRIER();
    }

    for (int r = 0; r < runs; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i) {
            gate_fn_1q(rho, (int)target);
        }
        BENCH_COMPILER_BARRIER();  /* once per run: prevents loop hoisting without per-call overhead */
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }

    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    free(run_times);

    bench_fill_perq_stats(&result, &stats, iterations, runs, "unknown");
    if (owns_qureg) destroyQureg(rho);
    return result;
}
```

**2. `bench_quest_2q_at()` — two-qubit, fixed pair:**

```c
bench_result_perq_t bench_quest_2q_at(qubit_t qubits, const char *gate_name,
                                       qubit_t q1, qubit_t q2,
                                       int iterations, int warmup, int runs,
                                       Qureg *preinit_qureg) {
    bench_result_perq_t result = {0};
    if (q1 >= qubits || q2 >= qubits) {
        fprintf(stderr, "bench_quest_2q_at: q1=%u or q2=%u >= qubits=%u\n",
                (unsigned)q1, (unsigned)q2, (unsigned)qubits);
        return result;
    }
    if (q1 == q2) {
        fprintf(stderr, "bench_quest_2q_at: q1 == q2 (%u)\n", (unsigned)q1);
        return result;
    }
    if (runs <= 0) {
        fprintf(stderr, "bench_quest_2q_at: runs=%d must be > 0\n", runs);
        return result;
    }
    if (runs > BENCH_PQ_MAX_RUNS) {
        fprintf(stderr, "bench_quest_2q_at: runs=%d exceeds BENCH_PQ_MAX_RUNS=%d\n",
                runs, BENCH_PQ_MAX_RUNS);
        return result;
    }
    result.qubits       = qubits;
    result.dim          = (dim_t)1 << qubits;
    result.gate_name    = gate_name;
    result.method       = "quest";
    result.target       = q1;
    result.target2      = q2;
    result.is_2q        = 1;
    result.runs         = runs;

    dim_t dim = (dim_t)1 << qubits;
    result.memory_bytes = (size_t)dim * (size_t)dim * sizeof(double) * 2;

    int owns_qureg = (preinit_qureg == NULL);
    Qureg local_qureg;
    Qureg rho;
    if (owns_qureg) {
        local_qureg = createDensityQureg((int)qubits);
        /* initDebugState ignores srand(42); state content differs from other backends.
         * Accepted limitation — see Known Limitations L1. */
        initDebugState(local_qureg);
        rho = local_qureg;
    } else {
        rho = *preinit_qureg;
    }

    void (*gate_fn_2q)(Qureg, int, int) = NULL;
    if      (strcmp(gate_name, "CX")   == 0) gate_fn_2q = applyControlledPauliX;
    else if (strcmp(gate_name, "SWAP") == 0) gate_fn_2q = applySwap;

    if (!gate_fn_2q) {
        if (owns_qureg) destroyQureg(rho);
        return result;
    }

    /* Auto-calibration probe — only when iterations == 0 (sentinel).
     * Pass iterations=0 to auto-calibrate; pass a positive value to use it directly. */
    if (iterations == 0) {
        uint64_t probe_t0 = bench_time_ns();
        for (int pi = 0; pi < BENCH_PQ_PROBE_ITERS; ++pi) {
            gate_fn_2q(rho, (int)q1, (int)q2);
            BENCH_COMPILER_BARRIER();
        }
        uint64_t probe_t1 = bench_time_ns();
        double probe_ms = (double)(probe_t1 - probe_t0) / 1e6 / BENCH_PQ_PROBE_ITERS;
        if (probe_ms <= 0.0) {
            iterations = BENCH_PQ_FALLBACK_ITERATIONS;
        } else {
            iterations = (int)ceil(MIN_MEASUREMENT_MS / probe_ms);
        }
        if (iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
            iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
    }
    result.iterations = iterations;

    double *run_times = malloc((size_t)runs * sizeof(double));
    if (!run_times) {
        if (owns_qureg) destroyQureg(rho);
        return result;
    }

    for (int i = 0; i < warmup; ++i) {
        gate_fn_2q(rho, (int)q1, (int)q2);
        BENCH_COMPILER_BARRIER();
    }

    for (int r = 0; r < runs; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i) {
            gate_fn_2q(rho, (int)q1, (int)q2);
        }
        BENCH_COMPILER_BARRIER();  /* once per run: prevents loop hoisting without per-call overhead */
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }

    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    free(run_times);

    bench_fill_perq_stats(&result, &stats, iterations, runs, "unknown");
    if (owns_qureg) destroyQureg(rho);
    return result;
}
```

### Constraints

- Both functions must be added inside the `#ifdef WITH_QUEST` / `#endif /* WITH_QUEST */`
  guards that already wrap this translation unit.
- The init-once pattern (`bench_quest_init` / `bench_quest_cleanup`) is unchanged.
- When `preinit_qureg` is NULL the function allocates with `createDensityQureg` /
  `initDebugState` and destroys with `destroyQureg` on return (same as the old pattern).
  When `preinit_qureg` is non-NULL the caller owns the Qureg lifecycle; the function
  neither creates nor destroys it.
- Do NOT add `#ifdef WITH_QUEST` guards around individual functions — the file-level
  guard at line 24 already covers all code in the file.
- Use heap allocation for `run_times` — never stack-allocated `double[200]` with a silent cap.
- Use `bench_fill_perq_stats()` with `runs` and `bw_notes = "unknown"` for stats assignment.
- QuEST's `initDebugState` ignores `srand(42)`; this is an accepted limitation for 1.0.

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

**Stack allocation fix**: same as Step D — use heap allocation for `run_times`.
Hard error if `runs > BENCH_PQ_MAX_RUNS`.

**State initialization**: `std::srand(42)` before initialization. Qulacs uses
`rand()` in `[-0.5, 0.5]`, which differs from both qlib (`init_random_mixed_state`)
and BLAS (`rand()` in `[0,1]`). This is a known limitation (see L1). The seed
ensures within-backend reproducibility; cross-backend state content differs.
Comment: `/* fixed seed for within-backend reproducibility only; cross-backend state content differs (see L1) */`.

### What to add

Append inside the `extern "C" { ... }` block and the enclosing `#ifdef WITH_QULACS`
guard, after the existing `bench_qulacs()` function:

**1. `bench_qulacs_at()` — single-qubit, fixed target:**

```cpp
bench_result_perq_t bench_qulacs_at(qubit_t qubits, const char *gate_name,
                                     qubit_t target, int iterations, int warmup, int runs,
                                     void *preinit_rho) {
    bench_result_perq_t result = {0};
    if (target >= qubits) {
        std::fprintf(stderr, "bench_qulacs_at: target=%u >= qubits=%u\n",
                     (unsigned)target, (unsigned)qubits);
        return result;
    }
    if (runs <= 0) {
        std::fprintf(stderr, "bench_qulacs_at: runs=%d must be > 0\n", runs);
        return result;
    }
    if (runs > BENCH_PQ_MAX_RUNS) {
        std::fprintf(stderr, "bench_qulacs_at: runs=%d exceeds BENCH_PQ_MAX_RUNS=%d\n",
                     runs, BENCH_PQ_MAX_RUNS);
        return result;
    }
    result.qubits       = qubits;
    result.dim          = (dim_t)1 << qubits;
    result.gate_name    = gate_name;
    result.method       = "qulacs";
    result.target       = target;
    result.target2      = QUBIT_NONE;
    result.is_2q        = 0;
    result.runs         = runs;

    ITYPE dim = static_cast<ITYPE>(1) << qubits;
    size_t matrix_size = static_cast<size_t>(dim) * static_cast<size_t>(dim) * sizeof(CTYPE);

    bool owns_rho = (preinit_rho == nullptr);
    CTYPE *rho = preinit_rho ? static_cast<CTYPE*>(preinit_rho)
                             : static_cast<CTYPE*>(std::malloc(matrix_size));
    if (!rho) return result;

    if (owns_rho) {
        std::srand(42);  /* fixed seed for within-backend reproducibility only;
                          * cross-backend state content differs (see Known Limitations L1) */
        size_t n_elems = static_cast<size_t>(dim) * static_cast<size_t>(dim);
        for (size_t i = 0; i < n_elems; ++i) {
            /* NOTE: distribution [-0.5,0.5] differs from qlib's init_random_mixed_state */
            double re = (double)std::rand() / RAND_MAX - 0.5;
            double im = (double)std::rand() / RAND_MAX - 0.5;
            rho[i] = CTYPE(re, im);
        }
    }

    result.memory_bytes = matrix_size;

    void (*gate_fn_1q)(UINT, CTYPE*, ITYPE) = nullptr;
    if      (std::strcmp(gate_name, "X") == 0) gate_fn_1q = dm_X_gate;
    else if (std::strcmp(gate_name, "H") == 0) gate_fn_1q = dm_H_gate;
    else if (std::strcmp(gate_name, "Z") == 0) gate_fn_1q = dm_Z_gate;

    if (!gate_fn_1q) { if (owns_rho) std::free(rho); return result; }

    /* Auto-calibration probe — only when iterations == 0 (sentinel).
     * Pass iterations=0 to auto-calibrate; pass a positive value to use it directly. */
    if (iterations == 0) {
        uint64_t probe_t0 = bench_time_ns();
        for (int pi = 0; pi < BENCH_PQ_PROBE_ITERS; ++pi) {
            gate_fn_1q(static_cast<UINT>(target), rho, dim);
            BENCH_COMPILER_BARRIER();
        }
        uint64_t probe_t1 = bench_time_ns();
        double probe_ms = (double)(probe_t1 - probe_t0) / 1e6 / BENCH_PQ_PROBE_ITERS;
        if (probe_ms <= 0.0) {
            iterations = BENCH_PQ_FALLBACK_ITERATIONS;
        } else {
            iterations = (int)std::ceil(MIN_MEASUREMENT_MS / probe_ms);
        }
        if (iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
            iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
    }
    result.iterations = iterations;

    double *run_times = static_cast<double*>(std::malloc(static_cast<size_t>(runs) * sizeof(double)));
    if (!run_times) { if (owns_rho) std::free(rho); return result; }

    for (int i = 0; i < warmup; ++i) {
        gate_fn_1q(static_cast<UINT>(target), rho, dim);
        BENCH_COMPILER_BARRIER();
    }

    for (int r = 0; r < runs; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i) {
            gate_fn_1q(static_cast<UINT>(target), rho, dim);
        }
        BENCH_COMPILER_BARRIER();  /* once per run: prevents loop hoisting without per-call overhead */
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }

    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    std::free(run_times);

    bench_fill_perq_stats(&result, &stats, iterations, runs, "unknown");
    if (owns_rho) std::free(rho);
    return result;
}
```

**2. `bench_qulacs_2q_at()` — two-qubit, fixed pair:**

```cpp
bench_result_perq_t bench_qulacs_2q_at(qubit_t qubits, const char *gate_name,
                                         qubit_t q1, qubit_t q2,
                                         int iterations, int warmup, int runs,
                                         void *preinit_rho) {
    bench_result_perq_t result = {0};
    if (q1 >= qubits || q2 >= qubits) {
        std::fprintf(stderr, "bench_qulacs_2q_at: q1=%u or q2=%u >= qubits=%u\n",
                     (unsigned)q1, (unsigned)q2, (unsigned)qubits);
        return result;
    }
    if (q1 == q2) {
        std::fprintf(stderr, "bench_qulacs_2q_at: q1 == q2 (%u)\n", (unsigned)q1);
        return result;
    }
    if (runs <= 0) {
        std::fprintf(stderr, "bench_qulacs_2q_at: runs=%d must be > 0\n", runs);
        return result;
    }
    if (runs > BENCH_PQ_MAX_RUNS) {
        std::fprintf(stderr, "bench_qulacs_2q_at: runs=%d exceeds BENCH_PQ_MAX_RUNS=%d\n",
                     runs, BENCH_PQ_MAX_RUNS);
        return result;
    }
    result.qubits       = qubits;
    result.dim          = (dim_t)1 << qubits;
    result.gate_name    = gate_name;
    result.method       = "qulacs";
    result.target       = q1;
    result.target2      = q2;
    result.is_2q        = 1;
    result.runs         = runs;

    ITYPE dim = static_cast<ITYPE>(1) << qubits;
    size_t matrix_size = static_cast<size_t>(dim) * static_cast<size_t>(dim) * sizeof(CTYPE);

    bool owns_rho = (preinit_rho == nullptr);
    CTYPE *rho = preinit_rho ? static_cast<CTYPE*>(preinit_rho)
                             : static_cast<CTYPE*>(std::malloc(matrix_size));
    if (!rho) return result;

    if (owns_rho) {
        std::srand(42);  /* fixed seed for within-backend reproducibility only;
                          * cross-backend state content differs (see Known Limitations L1) */
        size_t n_elems = static_cast<size_t>(dim) * static_cast<size_t>(dim);
        for (size_t i = 0; i < n_elems; ++i) {
            double re = (double)std::rand() / RAND_MAX - 0.5;
            double im = (double)std::rand() / RAND_MAX - 0.5;
            rho[i] = CTYPE(re, im);
        }
    }

    result.memory_bytes = matrix_size;

    void (*gate_fn_2q)(UINT, UINT, CTYPE*, ITYPE) = nullptr;
    if      (std::strcmp(gate_name, "CX")   == 0) gate_fn_2q = dm_CNOT_gate;
    else if (std::strcmp(gate_name, "SWAP") == 0) gate_fn_2q = dm_SWAP_gate;

    if (!gate_fn_2q) { if (owns_rho) std::free(rho); return result; }

    /* Auto-calibration probe — only when iterations == 0 (sentinel).
     * Pass iterations=0 to auto-calibrate; pass a positive value to use it directly. */
    if (iterations == 0) {
        uint64_t probe_t0 = bench_time_ns();
        for (int pi = 0; pi < BENCH_PQ_PROBE_ITERS; ++pi) {
            gate_fn_2q(static_cast<UINT>(q1), static_cast<UINT>(q2), rho, dim);
            BENCH_COMPILER_BARRIER();
        }
        uint64_t probe_t1 = bench_time_ns();
        double probe_ms = (double)(probe_t1 - probe_t0) / 1e6 / BENCH_PQ_PROBE_ITERS;
        if (probe_ms <= 0.0) {
            iterations = BENCH_PQ_FALLBACK_ITERATIONS;
        } else {
            iterations = (int)std::ceil(MIN_MEASUREMENT_MS / probe_ms);
        }
        if (iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
            iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
    }
    result.iterations = iterations;

    double *run_times = static_cast<double*>(std::malloc(static_cast<size_t>(runs) * sizeof(double)));
    if (!run_times) { if (owns_rho) std::free(rho); return result; }

    for (int i = 0; i < warmup; ++i) {
        gate_fn_2q(static_cast<UINT>(q1), static_cast<UINT>(q2), rho, dim);
        BENCH_COMPILER_BARRIER();
    }

    for (int r = 0; r < runs; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i) {
            gate_fn_2q(static_cast<UINT>(q1), static_cast<UINT>(q2), rho, dim);
        }
        BENCH_COMPILER_BARRIER();  /* once per run: prevents loop hoisting without per-call overhead */
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }

    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    std::free(run_times);

    bench_fill_perq_stats(&result, &stats, iterations, runs, "unknown");
    if (owns_rho) std::free(rho);
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
- Use heap allocation for `run_times` — never the stack-allocated `double[200]` pattern.
- Use `bench_fill_perq_stats()` with `runs` and `bw_notes = "unknown"` for stats assignment.
- When `preinit_rho` is non-NULL the function casts it to `CTYPE*` and uses it
  directly (skipping allocation and free). When NULL the function allocates internally
  as before. Step H passes a non-NULL pre-allocated buffer for state reuse.

---

## Step E2 — `bench_aer_dm.cpp` additions

**Depends on:** Step A
**File:** `benchmark/src/bench_aer_dm.cpp`
**Estimated time:** 1.5 h

### Context

Commit `67a7239` added the Aer-DM backend (`M_AER_DM = 6`) to the sweep benchmark.
This step adds `bench_aer_dm_at()` and `bench_aer_dm_2q_at()` to complete the per-qubit
coverage.

The method string is `"qiskit_aer_dm"` to match the existing `bench_aer_dm()` function
(which sets `result.method = "qiskit_aer_dm"`). This is critical for downstream CSV
analysis scripts that key on the method string.

**Real API**: the existing `bench_aer_dm.cpp` uses:
- `AER::QV::DensityMatrix<double>` — the density matrix class (from
  `simulators/density_matrix/densitymatrix.hpp`)
- `dm.initialize()` — resets the state to |0⟩⟨0|
- `dm.data()` — returns a `std::complex<double>*` buffer for direct element access
- `dm.apply_x(t)` — apply X gate at qubit `t`
- `dm.apply_unitary_matrix(reg, vec)` — apply a vectorized unitary; `reg` is
  `AER::reg_t{(AER::uint_t)t}`, `vec` is produced by `AER::Utils::vectorize_matrix(mat)`
- `dm.apply_cnot(q1, q2)` — CNOT gate
- `dm.apply_swap(q1, q2)` — SWAP gate
- Gate kind resolved once before warmup via a local `gate_kind_t` enum (BK_X, BK_H,
  BK_Z, BK_CX, BK_SWAP, BK_UNKNOWN) using if/else on `std::strcmp`; returns if
  `BK_UNKNOWN`.
- `AER::cvector_t` — the vectorized matrix type; H and Z are precomputed before warmup.
- The entire function body is wrapped in `try { ... } catch (const std::exception &e)`
  and `catch (...)` with `fprintf` to stderr and an early return.
- `BENCH_MAX_RUNS` (not `BENCH_PQ_MAX_RUNS`) is used in `bench_aer_dm()` to clamp
  `runs`; the new `_at` functions use `BENCH_PQ_MAX_RUNS` (which equals `BENCH_MAX_RUNS`).

**Random init pattern** (from `bench_aer_dm()`): after `dm.initialize()`, use
`std::mt19937_64` seeded with a function of `qubits`, then
`std::uniform_real_distribution<double>(-0.5, 0.5)` to fill all elements of `dm.data()`.
The seed must be fixed so the state is reproducible; use the same expression as the
existing function to keep within-backend consistency.

**Note**: the existing `bench_aer_dm()` sweeps all qubits per iteration (1Q) or all
pairs per iteration (2Q). The `_at` variants fix the target to a single qubit/pair and
call only that per iteration.

### What to add

Append inside the existing `#ifdef WITH_AER_DM` guard and `extern "C"` block, after
the existing `bench_aer_dm()` function:

**1. `bench_aer_dm_at()` — single-qubit, fixed target:**

```cpp
bench_result_perq_t bench_aer_dm_at(qubit_t qubits, const char *gate_name,
                                      qubit_t target, int iterations, int warmup, int runs,
                                      void *preinit_dm) {
    bench_result_perq_t result = {0};
    result.qubits    = qubits;
    result.dim       = (dim_t)1 << qubits;
    result.gate_name = gate_name;
    result.method    = "qiskit_aer_dm";
    result.target    = target;
    result.target2   = QUBIT_NONE;
    result.is_2q     = 0;
    result.runs      = runs;

    try {
        if (qubits == 0) return result;

        if (target >= qubits) {
            std::fprintf(stderr, "bench_aer_dm_at: target=%u >= qubits=%u\n",
                         (unsigned)target, (unsigned)qubits);
            return result;
        }
        if (runs <= 0) {
            std::fprintf(stderr, "bench_aer_dm_at: runs=%d must be > 0\n", runs);
            return result;
        }
        if (runs > BENCH_PQ_MAX_RUNS) {
            std::fprintf(stderr, "bench_aer_dm_at: runs=%d exceeds BENCH_PQ_MAX_RUNS=%d\n",
                         runs, BENCH_PQ_MAX_RUNS);
            return result;
        }

        size_t n_elems = (size_t)result.dim * (size_t)result.dim;
        result.memory_bytes = n_elems * sizeof(std::complex<double>);

        bool owns_dm = (preinit_dm == nullptr);
        AER::QV::DensityMatrix<double> local_dm;
        AER::QV::DensityMatrix<double> *dm_ptr;
        if (owns_dm) {
            local_dm = AER::QV::DensityMatrix<double>(qubits);
            local_dm.initialize();
            /* Random init — same pattern as bench_aer_dm(); fixed seed for within-backend
             * reproducibility; cross-backend state content differs (see Known Limitations L1) */
            std::mt19937_64 rng(static_cast<uint64_t>(qubits) * 6364136223846793005ULL
                                + 1442695040888963407ULL);
            std::uniform_real_distribution<double> dist(-0.5, 0.5);
            auto *buf = local_dm.data();
            for (size_t i = 0; i < n_elems; ++i)
                buf[i] = std::complex<double>(dist(rng), dist(rng));
            dm_ptr = &local_dm;
        } else {
            dm_ptr = static_cast<AER::QV::DensityMatrix<double>*>(preinit_dm);
        }
        AER::QV::DensityMatrix<double> &dm = *dm_ptr;

        /* Resolve gate kind once — no strcmp in hot path */
        enum gate_kind_t { BK_X, BK_H, BK_Z, BK_UNKNOWN };
        gate_kind_t gate_kind;
        if      (std::strcmp(gate_name, "X") == 0) gate_kind = BK_X;
        else if (std::strcmp(gate_name, "H") == 0) gate_kind = BK_H;
        else if (std::strcmp(gate_name, "Z") == 0) gate_kind = BK_Z;
        else                                        gate_kind = BK_UNKNOWN;

        if (gate_kind == BK_UNKNOWN) return result;

        /* Precompute vectorized matrices for H and Z */
        AER::cvector_t h_vec, z_vec;
        if (gate_kind == BK_H)
            h_vec = AER::Utils::vectorize_matrix(AER::Linalg::Matrix::H);
        if (gate_kind == BK_Z) {
            AER::cmatrix_t z_mat(2, 2);
            z_mat(0, 0) = {1, 0}; z_mat(0, 1) = {0, 0};
            z_mat(1, 0) = {0, 0}; z_mat(1, 1) = {-1, 0};
            z_vec = AER::Utils::vectorize_matrix(z_mat);
        }

        AER::reg_t qubit_reg{static_cast<AER::uint_t>(target)};

        /* Auto-calibration probe — only when iterations == 0 (sentinel).
         * Pass iterations=0 to auto-calibrate; pass a positive value to use it directly. */
        if (iterations == 0) {
            uint64_t probe_t0 = bench_time_ns();
            for (int pi = 0; pi < BENCH_PQ_PROBE_ITERS; ++pi) {
                switch (gate_kind) {
                case BK_X: dm.apply_x(target); break;
                case BK_H: dm.apply_unitary_matrix(qubit_reg, h_vec); break;
                case BK_Z: dm.apply_unitary_matrix(qubit_reg, z_vec); break;
                default:   break;
                }
                BENCH_COMPILER_BARRIER();
            }
            uint64_t probe_t1 = bench_time_ns();
            double probe_ms = (double)(probe_t1 - probe_t0) / 1e6 / BENCH_PQ_PROBE_ITERS;
            if (probe_ms <= 0.0) {
                iterations = BENCH_PQ_FALLBACK_ITERATIONS;
            } else {
                iterations = (int)std::ceil(MIN_MEASUREMENT_MS / probe_ms);
            }
            if (iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
                iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
        }
        result.iterations = iterations;

        /* Warmup */
        for (int i = 0; i < warmup; ++i) {
            switch (gate_kind) {
            case BK_X: dm.apply_x(target); break;
            case BK_H: dm.apply_unitary_matrix(qubit_reg, h_vec); break;
            case BK_Z: dm.apply_unitary_matrix(qubit_reg, z_vec); break;
            default:   break;
            }
            BENCH_COMPILER_BARRIER();
        }

        double *run_times = static_cast<double*>(
            std::malloc(static_cast<size_t>(runs) * sizeof(double)));
        if (!run_times) return result;

        for (int r = 0; r < runs; ++r) {
            uint64_t start = bench_time_ns();
            for (int i = 0; i < iterations; ++i) {
                switch (gate_kind) {
                case BK_X: dm.apply_x(target); break;
                case BK_H: dm.apply_unitary_matrix(qubit_reg, h_vec); break;
                case BK_Z: dm.apply_unitary_matrix(qubit_reg, z_vec); break;
                default:   break;
                }
            }
            BENCH_COMPILER_BARRIER();  /* once per run: prevents loop hoisting without per-call overhead */
            run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
        }

        bench_run_stats_t stats = bench_compute_stats(run_times, runs);
        std::free(run_times);

        bench_fill_perq_stats(&result, &stats, iterations, runs, "unknown");
    } catch (const std::exception &e) {
        std::fprintf(stderr, "bench_aer_dm_at: exception for %u qubits, gate %s: %s\n",
                     (unsigned)qubits, gate_name, e.what());
        return result;
    } catch (...) {
        std::fprintf(stderr, "bench_aer_dm_at: unknown exception for %u qubits, gate %s\n",
                     (unsigned)qubits, gate_name);
        return result;
    }
    return result;
}
```

**2. `bench_aer_dm_2q_at()` — two-qubit, fixed pair:**

```cpp
bench_result_perq_t bench_aer_dm_2q_at(qubit_t qubits, const char *gate_name,
                                         qubit_t q1, qubit_t q2,
                                         int iterations, int warmup, int runs,
                                         void *preinit_dm) {
    bench_result_perq_t result = {0};
    result.qubits    = qubits;
    result.dim       = (dim_t)1 << qubits;
    result.gate_name = gate_name;
    result.method    = "qiskit_aer_dm";
    result.target    = q1;
    result.target2   = q2;
    result.is_2q     = 1;
    result.runs      = runs;

    try {
        if (qubits == 0) return result;

        if (q1 >= qubits || q2 >= qubits) {
            std::fprintf(stderr, "bench_aer_dm_2q_at: q1=%u or q2=%u >= qubits=%u\n",
                         (unsigned)q1, (unsigned)q2, (unsigned)qubits);
            return result;
        }
        if (q1 == q2) {
            std::fprintf(stderr, "bench_aer_dm_2q_at: q1 == q2 (%u)\n", (unsigned)q1);
            return result;
        }
        if (runs <= 0) {
            std::fprintf(stderr, "bench_aer_dm_2q_at: runs=%d must be > 0\n", runs);
            return result;
        }
        if (runs > BENCH_PQ_MAX_RUNS) {
            std::fprintf(stderr, "bench_aer_dm_2q_at: runs=%d exceeds BENCH_PQ_MAX_RUNS=%d\n",
                         runs, BENCH_PQ_MAX_RUNS);
            return result;
        }

        size_t n_elems = (size_t)result.dim * (size_t)result.dim;
        result.memory_bytes = n_elems * sizeof(std::complex<double>);

        bool owns_dm = (preinit_dm == nullptr);
        AER::QV::DensityMatrix<double> local_dm;
        AER::QV::DensityMatrix<double> *dm_ptr;
        if (owns_dm) {
            local_dm = AER::QV::DensityMatrix<double>(qubits);
            local_dm.initialize();
            /* Random init — same pattern as bench_aer_dm() */
            std::mt19937_64 rng(static_cast<uint64_t>(qubits) * 6364136223846793005ULL
                                + 1442695040888963407ULL);
            std::uniform_real_distribution<double> dist(-0.5, 0.5);
            auto *buf = local_dm.data();
            for (size_t i = 0; i < n_elems; ++i)
                buf[i] = std::complex<double>(dist(rng), dist(rng));
            dm_ptr = &local_dm;
        } else {
            dm_ptr = static_cast<AER::QV::DensityMatrix<double>*>(preinit_dm);
        }
        AER::QV::DensityMatrix<double> &dm = *dm_ptr;

        /* Resolve gate kind once — no strcmp in hot path */
        enum gate_kind_t { BK_CX, BK_SWAP, BK_UNKNOWN };
        gate_kind_t gate_kind;
        if      (std::strcmp(gate_name, "CX")   == 0) gate_kind = BK_CX;
        else if (std::strcmp(gate_name, "SWAP") == 0) gate_kind = BK_SWAP;
        else                                           gate_kind = BK_UNKNOWN;

        if (gate_kind == BK_UNKNOWN) return result;

        /* Auto-calibration probe — only when iterations == 0 (sentinel).
         * Pass iterations=0 to auto-calibrate; pass a positive value to use it directly. */
        if (iterations == 0) {
            uint64_t probe_t0 = bench_time_ns();
            for (int pi = 0; pi < BENCH_PQ_PROBE_ITERS; ++pi) {
                switch (gate_kind) {
                case BK_CX:   dm.apply_cnot(q1, q2); break;
                case BK_SWAP: dm.apply_swap(q1, q2); break;
                default:      break;
                }
                BENCH_COMPILER_BARRIER();
            }
            uint64_t probe_t1 = bench_time_ns();
            double probe_ms = (double)(probe_t1 - probe_t0) / 1e6 / BENCH_PQ_PROBE_ITERS;
            if (probe_ms <= 0.0) {
                iterations = BENCH_PQ_FALLBACK_ITERATIONS;
            } else {
                iterations = (int)std::ceil(MIN_MEASUREMENT_MS / probe_ms);
            }
            if (iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
                iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
        }
        result.iterations = iterations;

        /* Warmup */
        for (int i = 0; i < warmup; ++i) {
            switch (gate_kind) {
            case BK_CX:   dm.apply_cnot(q1, q2); break;
            case BK_SWAP: dm.apply_swap(q1, q2); break;
            default:      break;
            }
            BENCH_COMPILER_BARRIER();
        }

        double *run_times = static_cast<double*>(
            std::malloc(static_cast<size_t>(runs) * sizeof(double)));
        if (!run_times) return result;

        for (int r = 0; r < runs; ++r) {
            uint64_t start = bench_time_ns();
            for (int i = 0; i < iterations; ++i) {
                switch (gate_kind) {
                case BK_CX:   dm.apply_cnot(q1, q2); break;
                case BK_SWAP: dm.apply_swap(q1, q2); break;
                default:      break;
                }
            }
            BENCH_COMPILER_BARRIER();  /* once per run: prevents loop hoisting without per-call overhead */
            run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
        }

        bench_run_stats_t stats = bench_compute_stats(run_times, runs);
        std::free(run_times);

        bench_fill_perq_stats(&result, &stats, iterations, runs, "unknown");
    } catch (const std::exception &e) {
        std::fprintf(stderr, "bench_aer_dm_2q_at: exception for %u qubits, gate %s: %s\n",
                     (unsigned)qubits, gate_name, e.what());
        return result;
    } catch (...) {
        std::fprintf(stderr, "bench_aer_dm_2q_at: unknown exception for %u qubits, gate %s\n",
                     (unsigned)qubits, gate_name);
        return result;
    }
    return result;
}
```

### Constraints

- Both functions must be inside the `#ifdef WITH_AER_DM` / `#endif` guard and inside
  the existing `extern "C"` block.
- Use the real class/method names from `bench_aer_dm.cpp` as shown above: no placeholder
  names like `aer_dm_alloc`, `AerDMState`, etc.
- The entire body of each function (excluding the result struct initialization before
  `try`) must be wrapped in `try { ... } catch (const std::exception &e)` and
  `catch (...)` — Aer throws on invalid arguments and allocation failure.
- When `preinit_dm` is non-NULL it is cast to `AER::QV::DensityMatrix<double>*` and
  used directly; no allocation or destruction occurs. When NULL, a local instance is
  constructed and initialized with random data as before. The caller owns the lifecycle
  of any passed object.
- Use switch/if-else on `gate_kind_t` (same as existing code), not function pointer
  dispatch.
- Use heap allocation for `run_times`.
- Use `bench_fill_perq_stats()` with `runs` and `bw_notes = "unknown"` for stats assignment.
- Method string must be `"qiskit_aer_dm"` to match the sweep benchmark.

---

## Step F — `bench_main.c`: CLI parsing

**Depends on:** Step A
**File:** `benchmark/src/bench_main.c`
**Estimated time:** 1 h

### What to change

**1. Initialization in `bench_parse_options()`** — add new fields to the designated
initializer for `bench_options_t opts`:

```c
    bench_options_t opts = {
        ...
        .per_qubit           = 0,
        .runs_explicit       = 0,
        .iterations_explicit = 0,
        .gate_filter         = NULL,
    };
```

**2. Mark explicit flags in `--runs` and `--iterations` parsing** — in the existing
argument parsing branches:

```c
        } else if (strcmp(argv[i], "--runs") == 0 && i + 1 < argc) {
            opts.runs = atoi(argv[++i]);
            opts.runs_explicit = 1;  /* ADD THIS LINE */
            if (opts.runs <= 0) {
                fprintf(stderr, "error: --runs requires a positive integer (got %d)\n",
                        opts.runs);
                exit(EXIT_FAILURE);
            }
        } else if (strcmp(argv[i], "--iterations") == 0 && i + 1 < argc) {
            opts.iterations = atoi(argv[++i]);
            opts.iterations_explicit = 1;  /* ADD THIS LINE */
            if (opts.iterations <= 0) {
                fprintf(stderr, "error: --iterations requires a positive integer (got %d)\n",
                        opts.iterations);
                exit(EXIT_FAILURE);
            }
```

Note: rejecting `--iterations 0` at parse time prevents the confusing case where a zero
value accidentally triggers auto-calibration. Any caller that wants auto-calibration
must rely on the default (i.e., omit `--iterations` entirely when `--per-qubit` is active).

**3. `--per-qubit` and `--gate` flags** — add branches in the parsing loop:

```c
        } else if (strcmp(argv[i], "--per-qubit") == 0) {
            opts.per_qubit = 1;
        } else if (strcmp(argv[i], "--gate") == 0 && i + 1 < argc) {
            /* Case-insensitive match is performed at call sites using strcasecmp;
             * store the raw string here for error reporting. */
            opts.gate_filter = argv[++i];
```

Note: `--gate` matching is case-insensitive. All `gate_filter` comparisons at call
sites should use `strcasecmp` (POSIX) rather than `strcmp`. This allows `--gate hd`
and `--gate HD` to produce the same result.

Note: `--gate` without `--per-qubit` silently has no effect in the current sweep path
(the sweep benchmark does not implement gate filtering). Consider printing a warning
when `opts.gate_filter != NULL && !opts.per_qubit` to alert the user their filter will
be ignored.

**4. Post-parse default override** — after the existing validation block, add:

```c
    /*
     * --per-qubit uses lighter defaults to keep total runtime manageable.
     * Only override iterations/runs if the user did NOT explicitly set them.
     * We track this with runs_explicit / iterations_explicit rather than
     * comparing against BENCH_DEFAULT_* values (which is ambiguous when the
     * user passes the same value as the default).
     */
    if (opts.per_qubit) {
        if (!opts.iterations_explicit)
            opts.iterations = BENCH_PQ_DEFAULT_ITERATIONS;
        if (!opts.runs_explicit)
            opts.runs = BENCH_PQ_DEFAULT_RUNS;
    }
```

**5. Mutual exclusion check** — add after the `per_qubit` default override:

```c
    if (opts.per_qubit && opts.pgfplots_output) {
        fprintf(stderr, "error: --per-qubit and --pgfplots are mutually exclusive\n");
        exit(1);
    }
```

**6. `print_usage()` update** — add the `--per-qubit` and `--gate` lines to the help text:

```c
    printf("  --per-qubit      Time each target qubit independently (default runs=%d, iters=%d)\n",
           BENCH_PQ_DEFAULT_RUNS, BENCH_PQ_DEFAULT_ITERATIONS);
    printf("  --gate NAME      Benchmark only the named gate (case-insensitive; e.g. --gate X, --gate CX)\n");
```

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

The `theoretical_bw_gbs` column represents an upper-bound estimate of memory bandwidth.
For in-place representations (qlib_packed, qlib_tiled):

```
theoretical_bw_gbs = (2 * memory_bytes) / (time_per_gate_ms * 1e-3) / 1e9
                   = 2.0 * memory_bytes / (time_per_gate_ms * 1e6)
```

This assumes a full-matrix read+write pass per gate. For single-qubit gates, the actual
bytes touched are approximately `memory_bytes` (only elements whose index differs in bit
`target`), so the formula overstates consumption by up to 2x for qlib backends.
For blas_dense, the formula understates consumption (5+ passes, not 2).

The column is named `theoretical_bw_gbs` in the CSV header — not `bw_gbs` — to make
the caveats visible to downstream consumers. The companion `bw_notes` column provides
a per-row annotation for downstream filtering.

If `time_per_gate_ms <= 0`, emit 0.0.

### What to add

**1. CSV header function:**

```c
void bench_print_perq_csv_header(void) {
    printf("qubits,dim,gate,method,runs,iterations,"
           "target,target2,is_2q,"
           "time_ms_mean,time_ms_std,time_ms_min,time_ms_median,time_ms_cv,"
           "time_per_gate_ms,time_per_gate_ms_median,time_per_gate_ms_min,"
           "ops_per_sec,theoretical_bw_gbs,bw_notes,memory_bytes\n");
}
```

**2. CSV row function:**

```c
void bench_print_perq_csv(const bench_result_perq_t *r) {
    /* skip failed measurements (zero-initialized on error) */
    if (r->time_ms <= 0.0) return;
    /* target2: always emit an integer. Use QUBIT_NONE (255) for 1Q gates.
     * Never leave the column empty — downstream tools must not need to handle
     * both numeric and empty fields in the same column. */
    unsigned t2 = r->is_2q ? (unsigned)r->target2 : (unsigned)QUBIT_NONE;
    printf("%u,%lld,%s,%s,%d,%d,%u,%u,%d,"
           "%.6f,%.6f,%.6f,%.6f,%.3f,%.9f,%.9f,%.9f,%.0f,%.3f,%s,%zu\n",
           r->qubits, (long long)r->dim, r->gate_name, r->method,
           r->runs, r->iterations,
           (unsigned)r->target, t2, r->is_2q,
           r->time_ms, r->time_ms_std, r->time_ms_min,
           r->time_ms_median, r->time_ms_cv,
           r->time_per_gate_ms, r->time_per_gate_ms_median, r->time_per_gate_ms_min,
           r->ops_per_sec,
           r->theoretical_bw_gbs,
           r->bw_notes ? r->bw_notes : "unknown",
           r->memory_bytes);
}
```

Note: `target2` is always emitted as an integer. For 1Q gates the value is `QUBIT_NONE`
(255). This keeps the column uniformly numeric; downstream tools never encounter an
empty field.

Note: `time_per_gate_ms_median` and `time_per_gate_ms_min` appear immediately after
`time_per_gate_ms` in the column order so that the three per-gate-time variants are
contiguous for easy column-range selection in downstream analysis.

**3. Console table function:**

```c
/**
 * @brief Print per-qubit results as a human-readable table.
 *
 * For 1Q gates: prints one row per target (full table).
 * For 2Q gates: prints a summary (min/mean/max across all pairs) rather than
 * all pairs (N*(N-1)/2 for symmetric gates, N*(N-1) for asymmetric ordered pairs),
 * which becomes unreadable at 12 qubits (66–132 pairs × methods).
 * Full per-target 2Q data is always available in the CSV.
 *
 * @param results  Array of per-qubit results (all same gate/qubits/method)
 * @param n        Number of entries
 */
void bench_print_perq_console(const bench_result_perq_t *results, int n) {
    if (n == 0) return;
    const bench_result_perq_t *r0 = &results[0];
    /* Note: individual result rows with time_ms <= 0.0 are skipped inside
     * the per-row printf calls below (failed measurements, zero-initialized on error). */

    printf("\n  Gate: %-6s  qubits: %u  method: %s\n",
           r0->gate_name, r0->qubits, r0->method);

    if (r0->is_2q) {
        /* Summary mode: compute min/mean/max across all pairs using time_per_gate_ms.
         * time_ms is NOT used here: auto-calibration may produce different iteration
         * counts per pair, making total-run time_ms values non-comparable. Only
         * time_per_gate_ms (= time_ms / iterations) is comparable across pairs. */
        /* Find first valid result to initialize min/max */
        int first_valid = -1;
        for (int i = 0; i < n; ++i) {
            if (results[i].time_ms > 0.0) { first_valid = i; break; }
        }
        if (first_valid < 0) {
            /* All measurements failed -- nothing to print */
            return;
        }
        double t_min = results[first_valid].time_per_gate_ms;
        double t_max = results[first_valid].time_per_gate_ms;
        int min_idx = first_valid, max_idx = first_valid;
        double t_sum = 0.0;
        int valid_count = 0;
        for (int i = 0; i < n; ++i) {
            /* skip failed measurements (zero-initialized on error) */
            if (results[i].time_ms <= 0.0) continue;
            double t = results[i].time_per_gate_ms;
            if (t < t_min) { t_min = t; min_idx = i; }
            if (t > t_max) { t_max = t; max_idx = i; }
            t_sum += t;
            ++valid_count;
        }
        double t_mean = (valid_count > 0) ? t_sum / (double)valid_count : 0.0;
        printf("  2Q summary across %d pairs (per-gate times):\n", valid_count);
        printf("  min: %.4f us/gate @ (q1=%u, q2=%u)   max: %.4f us/gate @ (q1=%u, q2=%u)"
               "   mean: %.4f us/gate   n_pairs=%d\n",
               t_min * 1000.0,
               (unsigned)results[min_idx].target, (unsigned)results[min_idx].target2,
               t_max * 1000.0,
               (unsigned)results[max_idx].target, (unsigned)results[max_idx].target2,
               t_mean * 1000.0, valid_count);
        printf("  (full per-pair data in CSV)\n");
    } else {
        printf("  %-6s  %12s  %14s  %12s  %8s  %12s\n",
               "target", "mean(us/gate)", "median(us/gate)", "min(us/gate)", "CV(%)", "ops/sec");
        printf("  %-6s  %12s  %14s  %12s  %8s  %12s\n",
               "------", "------------", "---------------", "------------", "--------", "------------");
        for (int i = 0; i < n; ++i) {
            const bench_result_perq_t *r = &results[i];
            /* skip failed measurements (zero-initialized on error) */
            if (r->time_ms <= 0.0) continue;
            /* Noise warning: CV > 50% indicates the measurement is dominated by
             * scheduler noise or thermal events; results should not be cited. */
            const char *noise_flag = (r->time_ms_cv > 50.0) ? " (noisy)" : "";
            printf("  %-6u  %12.4f  %14.4f  %12.4f  %8.2f  %12.0f%s\n",
                   r->target,
                   r->time_per_gate_ms * 1000.0,
                   r->time_per_gate_ms_median * 1000.0,
                   r->time_per_gate_ms_min * 1000.0,
                   r->time_ms_cv, r->ops_per_sec, noise_flag);
        }
    }
}
```

Note: the console 1Q table displays per-gate times in microseconds (`time_per_gate_ms * 1000`)
since per-gate times are typically sub-millisecond. Columns are `mean(us/gate)`,
`median(us/gate)`, `min(us/gate)` — all using `time_per_gate_ms{,_median,_min}` so rows
are directly comparable regardless of auto-calibrated `iterations`. A `(noisy)` flag is
appended when CV exceeds 50%. The 2Q console output also uses microseconds (`* 1000.0`)
and labels them `us/gate` to match the 1Q table; it shows `min/mean/max` with qubit-pair
coordinates rather than a full per-pair listing.

### Declarations in bench.h (already handled by Step A)

Step A declares these three functions in `bench.h`. Step G only provides the
implementations.

---

## Step H — `bench_main.c`: per-qubit orchestration loop

**Depends on:** Steps B, C, D, E, E2, F, G (all must be complete)
**File:** `benchmark/src/bench_main.c`
**Estimated time:** 3 h

This is the most complex step. It wires together everything from the preceding steps
and must leave the existing non-per-qubit code path completely unchanged.

### Design decisions to make before implementing

**Array sizing:** The per-qubit result arrays must be sized to accommodate the maximum
number of targets. For 1Q gates, max targets = `max_qubits`. For 2Q gates with ordered
pairs (asymmetric gates like CX), max targets = N*(N-1) — e.g. 12*11=132 pairs at 12
qubits — which exceeds `MAX_PAIRS` (128). Define:

```c
#define MAX_TARGETS  (2 * MAX_PAIRS)  /* 256; accommodates ordered pairs: N*(N-1) <= 2*MAX_PAIRS */
static_assert(MAX_TARGETS >= 2 * MAX_PAIRS,
    "MAX_TARGETS must accommodate ordered pairs");
```

Note: `2 * MAX_PAIRS * NUM_METHODS * sizeof(bench_result_perq_t)` is at most
256 × 7 × ~120 bytes ≈ 215 KB. The array is heap-allocated (see below), so stack size
is not a concern.

**`gate_def_2q_perq_t` struct for per-qubit mode**: Step H needs `is_symmetric` to choose
between `build_all_pairs` (symmetric gates) and `build_all_ordered_pairs` (asymmetric gates).
Do NOT add `is_symmetric` to the existing `gate_def_2q_t` struct — that would require
modifying existing struct initializers in the sweep path, which contradicts the constraint
"do not modify existing public functions." Instead, define a SEPARATE struct used only
within `perq_main()`:

```c
/**
 * @brief 2Q gate descriptor for per-qubit mode only.
 *
 * Defined as a separate type from gate_def_2q_t to avoid modifying existing sweep-path
 * initializers. Used exclusively within perq_main(); the sweep path is untouched.
 */
typedef struct {
    const char *name;
    void (*fn)(state_t*, qubit_t, qubit_t);
    const cplx_t *mat;
    int has_tiled;
    int is_symmetric;  /**< 1 if gate is symmetric (e.g. SWAP); 0 if asymmetric (e.g. CX) */
} gate_def_2q_perq_t;

static const gate_def_2q_perq_t gates_2q_perq[] = {
    {"CX",   cx,        CXMAT,   1, 0},   /* asymmetric: use ordered pairs */
    {"SWAP", swap_gate, SWAPMAT, 1, 1},   /* symmetric:  use unordered pairs */
    {NULL, NULL, NULL, 0, 0}
};
```

Use `gates_2q_perq` and `gate_def_2q_perq_t` throughout `perq_main()` for the 2Q gate
loop. The existing `gate_def_2q_t` and `gates_2q` table remain completely unmodified.

**Method count for per-qubit:** The effective per-qubit method list is:
qlib_packed, qlib_tiled, blas_dense, quest, qulacs, aer_dm — 6 methods.
Use the existing enum values (`M_QLIB=0, M_QLIB_TILED=1, M_BLAS=2, M_QUEST=4,
M_QULACS=5, M_AER_DM=6`), skipping index 3 (M_NAIVE). Array size remains `NUM_METHODS`.

**Result storage:** Use a flat 2D array:
`bench_result_perq_t perq[MAX_TARGETS][NUM_METHODS]`
heap-allocated once per `perq_main()` call to avoid stack concerns:

```c
bench_result_perq_t (*perq)[NUM_METHODS] =
    calloc(MAX_TARGETS, sizeof(*perq));
if (!perq) { fprintf(stderr, "per-qubit: allocation failed\n"); return 1; }
/* ... loops ... */
free(perq);
```

**State reuse across targets (MANDATORY):** State must be allocated once per
(qubits, method), not once per `_at` call. Without reuse, a 2Q run at 12 qubits with
blas_dense incurs 66 cycles of allocating and freeing a 256 MB matrix — adding
minutes of pure allocation overhead.

Use method-first, target-second loop ordering. Pass the pre-allocated state to each
`_at` call via the `preinit_state` / `preinit_rho` parameter (see Step A declarations).
The `_at` functions skip allocation when a non-NULL pre-allocated pointer is passed.

```c
/* MANDATORY: method-first, target-second for state reuse */
for each method m:
    allocate + init state once for (qubits, method m)
    for each target t:
        result = backend_at(..., preinit_state=state)
    free state
```

The target-first ordering (simpler but allocates per call) is NOT acceptable for
production use. Comment the chosen ordering to make the intent clear.

**Pre-allocated state types by method:**
- `M_QLIB` (packed): `state_t` with `MIXED_PACKED` layout
- `M_QLIB_TILED`: `state_t` with `MIXED_TILED` layout
- `M_BLAS`: raw `cplx_t*` of size `dim*dim`
- `M_QUEST`: `Qureg` via `createDensityQureg` / `initDebugState`; pass `&qureg`
- `M_QULACS`: raw `CTYPE*` of size `dim*dim`, random-initialized; pass `(void*)rho`
- `M_AER_DM`: `AER::QV::DensityMatrix<double>` stack object, random-initialized;
  pass `(void*)&dm`

### Where to insert the per-qubit branch

In `main()`, immediately after `bench_parse_options()` returns and the existing
pgfplots allocation block. The per-qubit path is a complete early return:

```c
    if (opts.per_qubit) {
        return perq_main(&opts);
    }
```

Define a static helper function `perq_main()` above `main()` to keep the function body
manageable. Prototype: `static int perq_main(const bench_options_t *opts)`.

### Structure of `perq_main()`

```
perq_main(opts):
    bench_result_perq_t (*perq)[NUM_METHODS] = calloc(MAX_TARGETS, sizeof(*perq))
    if !perq: return 1

    if csv_output: bench_print_perq_csv_header()
    else: print banner

    bench_quest_init()  [if WITH_QUEST]

    fprintf(stderr, "per-qubit: starting qubits=%d..%d\n",
            opts->min_qubits, opts->max_qubits)  /* progress reporting */

    results_count = 0  /* for --gate filter empty-match detection */

    for qubits in [opts->min_qubits .. opts->max_qubits step opts->step]:
        dim = 1 << qubits
        run_dense = (qubits <= 8)

        fprintf(stderr, "per-qubit: qubits=%d ...\n", qubits)  /* progress */

        if !csv_output: print "qubits: N, dim: D"

        /* ── 1Q gates ── */
        for each gate g in gates[]:
            gate_matches = (opts->gate_filter == NULL ||
                            strcasecmp(opts->gate_filter, g.name) == 0)
            if !gate_matches: continue
            results_count++

            n_targets = qubits
            memset(perq, 0, MAX_TARGETS * sizeof(*perq))

            /* Method-first loop for state reuse.
             *
             * Calibrate once per (qubits, method) on target 0, then pass that
             * calibrated_iterations to all _at calls for this method+qubit combination.
             * This ensures all targets use the same iterations count, making
             * time_per_gate_ms directly comparable across targets without needing
             * the L6 caveat about variable iteration counts.
             *
             * Sentinel convention: pass calibrated_iterations (a positive value) to
             * skip internal auto-calibration in the _at function. All backends —
             * including quest/qulacs/aer_dm — now receive a pre-allocated state and a
             * pre-calibrated iteration count; none use the iterations=0 sentinel during
             * normal Step H execution.
             *
             * Note: the initial probe on target 0 runs on the pre-allocated state,
             * which for the first gate of a method is the freshly initialised random
             * state (no prior operations). For subsequent gates within the same method
             * block the state may have been mutated; this is acceptable since calibration
             * only needs an order-of-magnitude estimate of gate cost.
             *
             * Warmup: all _at calls use `opts->warmup` (= BENCH_DEFAULT_WARMUP = 100).
             * The per-qubit mode uses the same warmup value as the sweep benchmark.
             * For low-qubit-count measurement where state fits in L1/L2 cache, even a
             * few warmup iterations suffice; `opts->warmup` is used as-is. */

            /* qlib_packed */
            {
                int calibrated_iterations = opts->iterations;
                state_t state = alloc_and_init(MIXED_PACKED, qubits, seed=42)
                /* Calibrate once for this (qubits, method) pair using target 0 */
                if (!opts->iterations_explicit) {
                    bench_result_perq_t cal = bench_qlib_packed_at(qubits, g.gate_name,
                        g.gate_fn, /*target=*/0, /*iterations=*/BENCH_PQ_PROBE_ITERS,
                        /*warmup=*/0, /*runs=*/1, &state);
                    double probe_ms = (cal.iterations > 0)
                        ? cal.time_ms / (double)cal.iterations : 0.0;
                    calibrated_iterations = (probe_ms <= 0.0)
                        ? BENCH_PQ_FALLBACK_ITERATIONS
                        : (int)ceil(MIN_MEASUREMENT_MS / probe_ms);
                    if (calibrated_iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
                        calibrated_iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
                }
                for target = 0 to qubits-1:
                    perq[target][M_QLIB] = bench_qlib_packed_at(...,
                        /*iterations=*/calibrated_iterations, /*warmup=*/opts->warmup,
                        /*runs=*/opts->runs, preinit_state=&state)
                state_free(&state)
            }
            /* qlib_tiled */
            {
                int calibrated_iterations = opts->iterations;
                state_t state = alloc_and_init(MIXED_TILED, qubits, seed=42)
                if (!opts->iterations_explicit) {
                    bench_result_perq_t cal = bench_qlib_tiled_at(qubits, g.gate_name,
                        g.gate_fn, /*target=*/0, /*iterations=*/BENCH_PQ_PROBE_ITERS,
                        /*warmup=*/0, /*runs=*/1, &state);
                    double probe_ms = (cal.iterations > 0)
                        ? cal.time_ms / (double)cal.iterations : 0.0;
                    calibrated_iterations = (probe_ms <= 0.0)
                        ? BENCH_PQ_FALLBACK_ITERATIONS
                        : (int)ceil(MIN_MEASUREMENT_MS / probe_ms);
                    if (calibrated_iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
                        calibrated_iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
                }
                for target = 0 to qubits-1:
                    perq[target][M_QLIB_TILED] = bench_qlib_tiled_at(...,
                        /*iterations=*/calibrated_iterations, /*warmup=*/opts->warmup,
                        /*runs=*/opts->runs, preinit_state=&state)
                state_free(&state)
            }
            /* blas_dense (capped at 8 qubits) */
            if run_dense:
            {
                int calibrated_iterations = opts->iterations;
                cplx_t *rho = alloc_and_init_dense(dim, seed=42)
                if (!opts->iterations_explicit) {
                    bench_result_perq_t cal = bench_blas_dense_at(qubits, g.gate_name,
                        g.gate_mat, /*target=*/0, /*iterations=*/BENCH_PQ_PROBE_ITERS,
                        /*warmup=*/0, /*runs=*/1, rho);
                    double probe_ms = (cal.iterations > 0)
                        ? cal.time_ms / (double)cal.iterations : 0.0;
                    calibrated_iterations = (probe_ms <= 0.0)
                        ? BENCH_PQ_FALLBACK_ITERATIONS
                        : (int)ceil(MIN_MEASUREMENT_MS / probe_ms);
                    if (calibrated_iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
                        calibrated_iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
                }
                for target = 0 to qubits-1:
                    perq[target][M_BLAS] = bench_blas_dense_at(...,
                        /*iterations=*/calibrated_iterations, /*warmup=*/opts->warmup,
                        /*runs=*/opts->runs, preinit_rho=rho)
                free(rho)
            }
            /* quest, qulacs, aer_dm: allocate state once per (qubits, method) and pass
             * it to every _at call. This collapses N alloc/init/free cycles (one per
             * target) into a single cycle per method — the dominant runtime saving at
             * 13+ qubits where each allocation is ~1 GB.
             *
             * Calibrate once per (qubits, method) on target 0 (same pattern as
             * qlib/blas above). Pass iterations=0 if opts->iterations_explicit is
             * false so the first target calibrates internally; reuse that count for all
             * subsequent targets (opts->iterations_explicit ? opts->iterations : 0 for
             * the first call, then the returned result.iterations for the rest).
             * Simpler: calibrate externally once with a probe call and pass the result. */
            #ifdef WITH_QUEST
            {
                int calibrated_iterations = opts->iterations;
                Qureg qureg = createDensityQureg((int)qubits);
                /* initDebugState for fixed deterministic state (ignores srand(42)) */
                initDebugState(qureg);
                if (!opts->iterations_explicit) {
                    /* Probe on target 0 to calibrate — reuse qureg */
                    bench_result_perq_t cal = bench_quest_at(qubits, g.gate_name,
                        /*target=*/0, /*iterations=*/BENCH_PQ_PROBE_ITERS,
                        /*warmup=*/0, /*runs=*/1, &qureg);
                    double probe_ms = (cal.iterations > 0)
                        ? cal.time_ms / (double)cal.iterations : 0.0;
                    calibrated_iterations = (probe_ms <= 0.0)
                        ? BENCH_PQ_FALLBACK_ITERATIONS
                        : (int)ceil(MIN_MEASUREMENT_MS / probe_ms);
                    if (calibrated_iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
                        calibrated_iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
                }
                for target = 0 to qubits-1:
                    perq[target][M_QUEST] = bench_quest_at(qubits, g.gate_name,
                        target, /*iterations=*/calibrated_iterations,
                        /*warmup=*/opts->warmup, /*runs=*/opts->runs, &qureg)
                destroyQureg(qureg);
            }
            #endif
            #ifdef WITH_QULACS
            {
                ITYPE dim_ql = (ITYPE)1 << qubits;
                size_t matrix_size = (size_t)dim_ql * (size_t)dim_ql * sizeof(CTYPE);
                CTYPE *ql_rho = (CTYPE*)malloc(matrix_size);
                if (ql_rho) {
                    /* Random init matching bench_qulacs_at() internal pattern */
                    std::srand(42);
                    size_t n_ql = (size_t)dim_ql * (size_t)dim_ql;
                    for (size_t i = 0; i < n_ql; ++i) {
                        double re = (double)std::rand() / RAND_MAX - 0.5;
                        double im = (double)std::rand() / RAND_MAX - 0.5;
                        ql_rho[i] = CTYPE(re, im);
                    }
                    int calibrated_iterations = opts->iterations;
                    if (!opts->iterations_explicit) {
                        bench_result_perq_t cal = bench_qulacs_at(qubits, g.gate_name,
                            /*target=*/0, /*iterations=*/BENCH_PQ_PROBE_ITERS,
                            /*warmup=*/0, /*runs=*/1, (void*)ql_rho);
                        double probe_ms = (cal.iterations > 0)
                            ? cal.time_ms / (double)cal.iterations : 0.0;
                        calibrated_iterations = (probe_ms <= 0.0)
                            ? BENCH_PQ_FALLBACK_ITERATIONS
                            : (int)ceil(MIN_MEASUREMENT_MS / probe_ms);
                        if (calibrated_iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
                            calibrated_iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
                    }
                    for target = 0 to qubits-1:
                        perq[target][M_QULACS] = bench_qulacs_at(qubits, g.gate_name,
                            target, /*iterations=*/calibrated_iterations,
                            /*warmup=*/opts->warmup, /*runs=*/opts->runs, (void*)ql_rho)
                    free(ql_rho);
                }
            }
            #endif
            #ifdef WITH_AER_DM
            {
                AER::QV::DensityMatrix<double> aer_dm(qubits);
                aer_dm.initialize();
                /* Random init matching bench_aer_dm_at() internal pattern */
                std::mt19937_64 rng((uint64_t)qubits * 6364136223846793005ULL
                                    + 1442695040888963407ULL);
                std::uniform_real_distribution<double> dist(-0.5, 0.5);
                size_t n_aer = (size_t)dim * (size_t)dim;
                auto *buf = aer_dm.data();
                for (size_t i = 0; i < n_aer; ++i)
                    buf[i] = std::complex<double>(dist(rng), dist(rng));
                int calibrated_iterations = opts->iterations;
                if (!opts->iterations_explicit) {
                    bench_result_perq_t cal = bench_aer_dm_at(qubits, g.gate_name,
                        /*target=*/0, /*iterations=*/BENCH_PQ_PROBE_ITERS,
                        /*warmup=*/0, /*runs=*/1, (void*)&aer_dm);
                    double probe_ms = (cal.iterations > 0)
                        ? cal.time_ms / (double)cal.iterations : 0.0;
                    calibrated_iterations = (probe_ms <= 0.0)
                        ? BENCH_PQ_FALLBACK_ITERATIONS
                        : (int)ceil(MIN_MEASUREMENT_MS / probe_ms);
                    if (calibrated_iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
                        calibrated_iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
                }
                for target = 0 to qubits-1:
                    perq[target][M_AER_DM] = bench_aer_dm_at(qubits, g.gate_name,
                        target, /*iterations=*/calibrated_iterations,
                        /*warmup=*/opts->warmup, /*runs=*/opts->runs, (void*)&aer_dm)
                /* aer_dm destructor frees internal buffer */
            }
            #endif

            /* Output */
            if csv_output:
                for target in 0..qubits-1:
                    for each active method m:
                        /* bench_print_perq_csv skips rows where time_ms <= 0.0
                         * (zero-initialized on error from _at functions) */
                        bench_print_perq_csv(&perq[target][m])
            else:
                for each active method m:
                    /* IMPORTANT: perq[][m] is not contiguous -- must copy into flat array.
                     * Taking &perq[0][m] as a pointer to a contiguous array of
                     * bench_result_perq_t is WRONG: perq is declared as
                     * bench_result_perq_t perq[MAX_TARGETS][NUM_METHODS], so
                     * perq[t][m] elements for fixed m are strided by NUM_METHODS
                     * elements apart in memory. */
                    bench_result_perq_t flat[MAX_TARGETS];
                    for (int t = 0; t < n_targets; ++t)
                        flat[t] = perq[t][m];
                    /* bench_print_perq_console skips individual rows where time_ms <= 0.0 */
                    bench_print_perq_console(flat, n_targets);

        /* ── 2Q gates ── */
        /* Pair arrays: sized for ordered pairs (N*(N-1) entries, at most 2*MAX_PAIRS).
         * build_all_pairs     generates N*(N-1)/2 unordered pairs — used for symmetric gates (SWAP).
         * build_all_ordered_pairs generates N*(N-1) ordered pairs  — used for asymmetric gates (CX).
         * Asymmetric gates require both (q_control, q_target) and (q_target, q_control) because
         * memory-layout cost can differ between the two orderings (packed/tiled layout asymmetry). */
        qubit_t q1s[2 * MAX_PAIRS], q2s[2 * MAX_PAIRS]

        for each gate g in gates_2q_perq[]:
            gate_matches = (opts->gate_filter == NULL ||
                            strcasecmp(opts->gate_filter, g.name) == 0)
            if !gate_matches: continue
            results_count++

            /* Choose pair generator: symmetric gates (SWAP) use unordered pairs;
             * asymmetric gates (CX) use ordered pairs to expose layout asymmetry. */
            int n_pairs;
            if (g.is_symmetric) {
                n_pairs = build_all_pairs(qubits, q1s, q2s);       /* N*(N-1)/2 pairs */
            } else {
                n_pairs = build_all_ordered_pairs(qubits, q1s, q2s); /* N*(N-1) pairs */
            }
            /* Guard applies to BOTH pair generators:
             *   - build_all_ordered_pairs returns 0 on error (qubits > MAX_QUBITS_FOR_ORDERED_PAIRS).
             *   - build_all_pairs uses assert(n < MAX_PAIRS) which is compiled away with
             *     -DNDEBUG; at 17 qubits, N*(N-1)/2 = 136 > MAX_PAIRS = 128, causing
             *     a silent buffer overflow without this runtime check.
             * This guard must be unconditional and must appear after both branches. */
            if (n_pairs <= 0 || n_pairs > MAX_TARGETS) {
                fprintf(stderr, "per-qubit: skipping gate %s at %d qubits: "
                        "n_pairs=%d out of range (max %d)\n",
                        g.name, (int)qubits, n_pairs, MAX_TARGETS);
                continue;
            }

            memset(perq, 0, MAX_TARGETS * sizeof(*perq))

            /* Method-first loop for state reuse — same calibrate-once pattern as 1Q.
             * Calibrate on pair (q1s[0], q2s[0]) to get iterations for all pairs. */
            {
                int calibrated_iterations = opts->iterations;
                state_t state = alloc_and_init(MIXED_PACKED, qubits, seed=42)
                if (!opts->iterations_explicit) {
                    bench_result_perq_t cal = bench_qlib_packed_2q_at(qubits, g.gate_name,
                        g.gate_fn, /*q1=*/q1s[0], /*q2=*/q2s[0],
                        /*iterations=*/BENCH_PQ_PROBE_ITERS, /*warmup=*/0, /*runs=*/1, &state);
                    double probe_ms = (cal.iterations > 0)
                        ? cal.time_ms / (double)cal.iterations : 0.0;
                    calibrated_iterations = (probe_ms <= 0.0)
                        ? BENCH_PQ_FALLBACK_ITERATIONS
                        : (int)ceil(MIN_MEASUREMENT_MS / probe_ms);
                    if (calibrated_iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
                        calibrated_iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
                }
                for p = 0 to n_pairs-1:
                    perq[p][M_QLIB] = bench_qlib_packed_2q_at(... q1s[p], q2s[p] ...,
                        /*iterations=*/calibrated_iterations, /*warmup=*/opts->warmup,
                        /*runs=*/opts->runs, preinit_state=&state)
                state_free(&state)
            }
            if g.has_tiled:
            {
                int calibrated_iterations = opts->iterations;
                state_t state = alloc_and_init(MIXED_TILED, qubits, seed=42)
                if (!opts->iterations_explicit) {
                    bench_result_perq_t cal = bench_qlib_tiled_2q_at(qubits, g.gate_name,
                        g.gate_fn, /*q1=*/q1s[0], /*q2=*/q2s[0],
                        /*iterations=*/BENCH_PQ_PROBE_ITERS, /*warmup=*/0, /*runs=*/1, &state);
                    double probe_ms = (cal.iterations > 0)
                        ? cal.time_ms / (double)cal.iterations : 0.0;
                    calibrated_iterations = (probe_ms <= 0.0)
                        ? BENCH_PQ_FALLBACK_ITERATIONS
                        : (int)ceil(MIN_MEASUREMENT_MS / probe_ms);
                    if (calibrated_iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
                        calibrated_iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
                }
                for p = 0 to n_pairs-1:
                    perq[p][M_QLIB_TILED] = bench_qlib_tiled_2q_at(...,
                        /*iterations=*/calibrated_iterations, /*warmup=*/opts->warmup,
                        /*runs=*/opts->runs, preinit_state=&state)
                state_free(&state)
            }
            if run_dense:
            {
                int calibrated_iterations = opts->iterations;
                cplx_t *rho = alloc_and_init_dense(dim, seed=42)
                if (!opts->iterations_explicit) {
                    bench_result_perq_t cal = bench_blas_dense_2q_at(qubits, g.gate_name,
                        g.gate_mat, /*q1=*/q1s[0], /*q2=*/q2s[0],
                        /*iterations=*/BENCH_PQ_PROBE_ITERS, /*warmup=*/0, /*runs=*/1, rho);
                    double probe_ms = (cal.iterations > 0)
                        ? cal.time_ms / (double)cal.iterations : 0.0;
                    calibrated_iterations = (probe_ms <= 0.0)
                        ? BENCH_PQ_FALLBACK_ITERATIONS
                        : (int)ceil(MIN_MEASUREMENT_MS / probe_ms);
                    if (calibrated_iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
                        calibrated_iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
                }
                for p = 0 to n_pairs-1:
                    perq[p][M_BLAS] = bench_blas_dense_2q_at(... q1s[p], q2s[p] ...,
                        /*iterations=*/calibrated_iterations, /*warmup=*/opts->warmup,
                        /*runs=*/opts->runs, preinit_rho=rho)
                free(rho)
            }
            /* quest, qulacs, aer_dm: allocate state once per (qubits, method) and
             * pass to every _at call — same state-reuse pattern as qlib/blas above.
             * At 13 qubits with 156 ordered CX pairs this eliminates 155 redundant
             * ~1 GB alloc/init/free cycles per backend. */
            #ifdef WITH_QUEST
            {
                int calibrated_iterations = opts->iterations;
                Qureg qureg = createDensityQureg((int)qubits);
                initDebugState(qureg);
                if (!opts->iterations_explicit) {
                    bench_result_perq_t cal = bench_quest_2q_at(qubits, g.gate_name,
                        /*q1=*/q1s[0], /*q2=*/q2s[0],
                        /*iterations=*/BENCH_PQ_PROBE_ITERS, /*warmup=*/0, /*runs=*/1, &qureg);
                    double probe_ms = (cal.iterations > 0)
                        ? cal.time_ms / (double)cal.iterations : 0.0;
                    calibrated_iterations = (probe_ms <= 0.0)
                        ? BENCH_PQ_FALLBACK_ITERATIONS
                        : (int)ceil(MIN_MEASUREMENT_MS / probe_ms);
                    if (calibrated_iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
                        calibrated_iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
                }
                for p = 0 to n_pairs-1:
                    perq[p][M_QUEST] = bench_quest_2q_at(qubits, g.gate_name,
                        q1s[p], q2s[p], /*iterations=*/calibrated_iterations,
                        /*warmup=*/opts->warmup, /*runs=*/opts->runs, &qureg)
                destroyQureg(qureg);
            }
            #endif
            #ifdef WITH_QULACS
            {
                ITYPE dim_ql = (ITYPE)1 << qubits;
                size_t matrix_size = (size_t)dim_ql * (size_t)dim_ql * sizeof(CTYPE);
                CTYPE *ql_rho = (CTYPE*)malloc(matrix_size);
                if (ql_rho) {
                    std::srand(42);
                    size_t n_ql = (size_t)dim_ql * (size_t)dim_ql;
                    for (size_t i = 0; i < n_ql; ++i) {
                        double re = (double)std::rand() / RAND_MAX - 0.5;
                        double im = (double)std::rand() / RAND_MAX - 0.5;
                        ql_rho[i] = CTYPE(re, im);
                    }
                    int calibrated_iterations = opts->iterations;
                    if (!opts->iterations_explicit) {
                        bench_result_perq_t cal = bench_qulacs_2q_at(qubits, g.gate_name,
                            /*q1=*/q1s[0], /*q2=*/q2s[0],
                            /*iterations=*/BENCH_PQ_PROBE_ITERS, /*warmup=*/0, /*runs=*/1,
                            (void*)ql_rho);
                        double probe_ms = (cal.iterations > 0)
                            ? cal.time_ms / (double)cal.iterations : 0.0;
                        calibrated_iterations = (probe_ms <= 0.0)
                            ? BENCH_PQ_FALLBACK_ITERATIONS
                            : (int)ceil(MIN_MEASUREMENT_MS / probe_ms);
                        if (calibrated_iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
                            calibrated_iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
                    }
                    for p = 0 to n_pairs-1:
                        perq[p][M_QULACS] = bench_qulacs_2q_at(qubits, g.gate_name,
                            q1s[p], q2s[p], /*iterations=*/calibrated_iterations,
                            /*warmup=*/opts->warmup, /*runs=*/opts->runs, (void*)ql_rho)
                    free(ql_rho);
                }
            }
            #endif
            #ifdef WITH_AER_DM
            {
                AER::QV::DensityMatrix<double> aer_dm(qubits);
                aer_dm.initialize();
                std::mt19937_64 rng((uint64_t)qubits * 6364136223846793005ULL
                                    + 1442695040888963407ULL);
                std::uniform_real_distribution<double> dist(-0.5, 0.5);
                size_t n_aer = (size_t)dim * (size_t)dim;
                auto *buf = aer_dm.data();
                for (size_t i = 0; i < n_aer; ++i)
                    buf[i] = std::complex<double>(dist(rng), dist(rng));
                int calibrated_iterations = opts->iterations;
                if (!opts->iterations_explicit) {
                    bench_result_perq_t cal = bench_aer_dm_2q_at(qubits, g.gate_name,
                        /*q1=*/q1s[0], /*q2=*/q2s[0],
                        /*iterations=*/BENCH_PQ_PROBE_ITERS, /*warmup=*/0, /*runs=*/1,
                        (void*)&aer_dm);
                    double probe_ms = (cal.iterations > 0)
                        ? cal.time_ms / (double)cal.iterations : 0.0;
                    calibrated_iterations = (probe_ms <= 0.0)
                        ? BENCH_PQ_FALLBACK_ITERATIONS
                        : (int)ceil(MIN_MEASUREMENT_MS / probe_ms);
                    if (calibrated_iterations > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
                        calibrated_iterations = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
                }
                for p = 0 to n_pairs-1:
                    perq[p][M_AER_DM] = bench_aer_dm_2q_at(qubits, g.gate_name,
                        q1s[p], q2s[p], /*iterations=*/calibrated_iterations,
                        /*warmup=*/opts->warmup, /*runs=*/opts->runs, (void*)&aer_dm)
                /* aer_dm destructor frees internal buffer */
            }
            #endif

            /* Output */
            /* CSV note: for asymmetric gates (CX), both (q1,q2) and (q2,q1) appear
             * as separate rows. The `target` column holds q1 and `target2` holds q2;
             * control/target semantics follow the same convention as the sweep benchmark
             * (q1 = control for CX). Downstream analysis should group by (target, target2)
             * or use the ordered pair to detect layout asymmetry. */
            if csv_output:
                for p in 0..n_pairs-1:
                    for each active method m:
                        /* bench_print_perq_csv skips rows where time_ms <= 0.0
                         * (zero-initialized on error from _at functions) */
                        bench_print_perq_csv(&perq[p][m])
            else:
                for each active method m:
                    /* IMPORTANT: perq[][m] is not contiguous -- must copy into flat array.
                     * Taking &perq[0][m] as a pointer to a contiguous array of
                     * bench_result_perq_t is WRONG: perq[p][m] elements for fixed m are
                     * strided by NUM_METHODS elements apart in memory. */
                    bench_result_perq_t flat[MAX_TARGETS];
                    for (int p = 0; p < n_pairs; ++p)
                        flat[p] = perq[p][m];
                    /* bench_print_perq_console prints 2Q summary (min/mean/max), not all pairs */
                    /* bench_print_perq_console skips individual rows where time_ms <= 0.0 */
                    bench_print_perq_console(flat, n_pairs);

    bench_quest_cleanup()  [if WITH_QUEST]

    /* --gate filter empty-match detection */
    if results_count == 0 && opts->gate_filter != NULL:
        fprintf(stderr, "warning: --gate '%s' matched no gates\n", opts->gate_filter)
        free(perq)
        return EXIT_FAILURE

    free(perq)
    return 0
```

### Constraints

- The existing `for (qubit_t qubits = ...)` main loop, the `pgfplots_data_t` struct,
  and all existing output paths must remain completely untouched.
- `bench_quest_init()` / `bench_quest_cleanup()` calls in `perq_main()` are independent
  of those in the existing path; each path calls them exactly once.
- The `SET_SWEEP` macro is NOT used in the per-qubit path (per-qubit results have
  implicit sweep_size=1; `time_per_gate_ms` is computed in the `_at` functions directly).
- State reuse (method-first loop ordering) is MANDATORY for all backends: qlib, blas,
  quest, qulacs, and aer_dm. Each allocates once per (qubits, method) before the target
  loop and frees/destroys after it.
- All `gate_filter` comparisons use `strcasecmp`, not `strcmp`.

---

## Step I — `docs/BENCHMARK.md` update

**Depends on:** Steps A, F, G, H (all changes final)
**File:** `docs/BENCHMARK.md`
**Estimated time:** 1 h

### What to add / update

**1. CLI Reference table** — add `--per-qubit` and `--gate` rows:

| Flag | Argument | Default | Description |
|------|----------|---------|-------------|
| `--per-qubit` | — | off | Time each target qubit independently; sets `--runs=5 --iterations=100` if not overridden; mutually exclusive with `--pgfplots` |
| `--gate` | NAME | (all) | Benchmark only the named gate (case-insensitive; e.g. `X`, `CX`); applies to per-qubit mode; unmatched names produce a warning and `EXIT_FAILURE` |

**2. Configuration constants** — add to the constants table:

| Constant | Value | Description |
|----------|-------|-------------|
| `BENCH_PQ_DEFAULT_RUNS` | 5 | Default `--runs` when `--per-qubit` is active |
| `BENCH_PQ_DEFAULT_ITERATIONS` | 100 | Default `--iterations` when `--per-qubit` is active |
| `BENCH_PQ_MAX_RUNS` | `BENCH_MAX_RUNS` (200) | Hard cap on `--runs` for all backends; tied to `bench_compute_stats` sort buffer size |
| `BENCH_PQ_FALLBACK_ITERATIONS` | 10000 | Iterations used when probe duration rounds to 0 ns |
| `QUBIT_NONE` | 0xFF | Sentinel for `bench_result_perq_t.target2` on 1Q measurements |

**3. Key Data Structures section** — add `bench_result_perq_t` entry documenting all
fields, noting the difference from `bench_result_t` (no `sweep_size`; adds `target`,
`target2`, `is_2q`, `theoretical_bw_gbs`, `bw_notes`; `time_per_gate_ms` =
time_ms/iterations rather than time_ms/(iterations×sweep_size)).

**4. CSV Output section** — add a subsection for the per-qubit CSV format:

- Header: `qubits,dim,gate,method,runs,iterations,target,target2,is_2q,time_ms_mean,time_ms_std,time_ms_min,time_ms_median,time_ms_cv,time_per_gate_ms,time_per_gate_ms_median,time_per_gate_ms_min,ops_per_sec,theoretical_bw_gbs,bw_notes,memory_bytes`
- Document `theoretical_bw_gbs` formula and its caveats per backend (see struct field doc)
- Document `bw_notes` values: `"half_matrix"` (qlib 1Q), `"three_quarter_matrix"` (qlib 2Q),
  `"multi_pass_est"` (blas), `"unknown"` (quest/qulacs/aer_dm)
- Note: `target2` column is always an integer; 1Q gate rows use value 255 (`QUBIT_NONE`)
- Note: measurement semantics (steady-state single-target throughput, not circuit throughput)

**5. Adapter functions section** — list all `_at` signatures (or add them to the
existing function listing table), including `bench_aer_dm_at` and `bench_aer_dm_2q_at`.

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
   the density matrix layout. This is steady-state single-target throughput — the
   branch predictor has saturated and the prefetcher has learned the access pattern
   before the timed loop starts. This is NOT circuit-representative throughput; the
   sweep benchmark is a better proxy for that.

2. **Interpreting target-position differences**: In packed lower-triangular storage,
   lower-index qubits modulate inner indices of the data structure, while higher-index
   qubits modulate outer indices. A gate on qubit 0 touches many cache lines in a
   scattered pattern; a gate on qubit n-1 touches a contiguous block. The
   `--per-qubit` output shows whether this asymmetry is measurable in practice (it
   typically becomes visible at qubit counts where the density matrix no longer fits
   in L3 cache).

3. **Bandwidth column**: `theoretical_bw_gbs` = estimated memory bandwidth (GB/s) as
   2 × (state size) / (time per gate application). This is an upper bound for in-place
   representations and severely understates consumption for the BLAS baseline (which
   uses ~5 passes). Use this to compare representations relative to each other and to
   hardware memory bandwidth limits (measurable via `stream` benchmark), not for
   absolute physical claims. Always note backend-specific caveats when citing this
   column; filter by the companion `bw_notes` CSV column.

4. **Citing per-qubit data in a thesis chapter**: The per-qubit CSV is suitable for
   pgfplots `\addplot table` with `x=target, y=time_per_gate_ms` to produce a
   target-position sweep plot for a fixed (gate, qubits, method) combination. For
   thesis narrative: describe whether packed, tiled, QuEST, Qulacs, and Aer-DM show
   flat or position-dependent profiles. A flat profile supports the claim that the
   implementation is target-position agnostic. A sloped profile is a finding worth
   discussing.

5. **Runtime estimates** at the per-qubit defaults (20 runs × 100 iterations):
   compute approximate total measurement count and wall-clock estimate at the default
   qubit range (2–12), and compare to the full-default sweep benchmark estimate to
   justify the lighter per-qubit defaults.

---

## Time Estimates

### Implementation time

| Step | Description | Estimate |
|------|-------------|----------|
| A | bench.h additions | 1.5 h |
| B | bench_mixed.c additions | 2 h |
| C | bench_baselines.c additions | 2 h |
| D | bench_quest.c additions | 1.5 h |
| E | bench_qulacs.cpp additions | 1.5 h |
| E2 | bench_aer_dm.cpp additions | 1.5 h |
| F | bench_main.c CLI parsing | 1 h |
| G | bench_main.c output functions | 2 h |
| H | bench_main.c orchestration | 3 h |
| I | docs/BENCHMARK.md update | 1 h |
| J | Thesis notes | 1.5 h |
| **Total sequential** | | **18 h** |
| **Critical path** (A→{B,C,D,E,E2,F,G}→H→I→J) | | **9.5 h** |
| **With max parallelism** | | **A(1.5h) + H(3h) + I(1h) + J(1.5h) ≈ 7 h elapsed** |

### Benchmark execution time (wall clock)

Rough estimates for a representative run (`--per-qubit --min-qubits 2 --max-qubits 12`,
all 6 methods enabled, default 5 runs × 100 iterations, 3 gates: X, H, CX).

Key costs per `_at` call:
- 12-qubit state size: 4 GB (packed/tiled), 256 MB (dense, capped at 8Q),
  ~64 GB (QuEST DM at 13Q, but capped at 12Q here), ~1 GB (Qulacs/Aer-DM at 12Q)
- Time per call: dominated by gate application time; at 12Q, ~1–10 ms per run × 100 iters
- Allocation per call: eliminated for ALL backends via method-first state-reuse ordering

Note: for asymmetric gates (CX), `build_all_ordered_pairs` produces N*(N-1) pairs instead
of N*(N-1)/2. At 12 qubits this doubles the CX pair count from 66 to 132. The estimates
below reflect this for the CX gate.

| Qubit range | 1Q gates | 2Q pairs (SWAP) | 2Q pairs (CX, ordered) | Rough total |
|-------------|----------|-----------------|------------------------|-------------|
| 2–8 | ~4 targets avg | ~3 pairs avg | ~6 pairs avg | < 5 min |
| 9–12 | ~11 targets avg | ~45 pairs avg | ~90 pairs avg | 25–50 min |
| **Total** | | | | **~30–55 min** |

These estimates assume qlib backends dominate runtime (QuEST/Qulacs/Aer-DM are similar
per-gate). With state reuse now implemented for all backends, the N alloc/free cycles
per 2Q gate at 12 qubits (N = 132 for ordered CX pairs) collapse to 1 alloc per method.
At 13 qubits, QuEST/Qulacs/Aer-DM each manage ~1 GB states; without reuse that would
be 156 × 1 GB alloc/init/free cycles per backend per gate type. The improvement to the
allocation portion alone is roughly 100–500x depending on the allocator and OS page
fault cost.

---

## Known Limitations

### L1 — State initialization differs between backends

qlib uses `init_random_mixed_state()` (seeded with `srand(42)`), BLAS uses `rand()` in
`[0,1]` (seeded with `srand(42)`), Qulacs uses `rand()` in `[-0.5,0.5]` (seeded with
`srand(42)`). QuEST's `initDebugState` produces a fixed deterministic state and ignores
`srand(42)` entirely. All backends that honour `srand` use seed 42 for within-backend
reproducibility, but cross-backend state content differs.

The `srand(42)` comment in all `_at` functions reads:
`/* fixed seed for within-backend reproducibility only; cross-backend state content differs (see L1) */`

At low qubit counts where the matrix fits in L1/L2 cache, different data values can
affect SIMD fast-paths in Qulacs/QuEST internals. This is an accepted limitation for a
cross-framework comparison. Document it when citing cross-backend comparisons.

Ideal fix (out of scope for this plan): extract a single
`bench_fill_random_dm(double *rho, size_t n, uint32_t seed)` C helper and call it from
all backends via an `extern "C"` shim.

### L2 — Measurement semantics vs. circuit workloads

See §Measurement Semantics at the top of this plan. Per-qubit numbers must not be
directly compared to sweep-benchmark `time_per_gate_ms` values without qualification.

### L3 — Bandwidth column validity varies by backend

`theoretical_bw_gbs` is computed as `2 * memory_bytes / time_per_gate_ms`. This is a
reasonable proxy for qlib in-place backends for 1Q gates (`bw_notes = "half_matrix"`)
and an overestimate by ~4/3 for 2Q gates (`bw_notes = "three_quarter_matrix"`). For
BLAS it understates consumption by ~2.5x (`bw_notes = "multi_pass_est"`). For
QuEST/Qulacs/Aer-DM, internal memory access patterns are opaque (`bw_notes = "unknown"`).
Use the `bw_notes` CSV column to filter or scale in downstream analysis.

### L4 — NUMA and thread affinity

On a multi-socket machine, different target qubits may touch different NUMA nodes
depending on how the operating system pages the state array. This can produce a
2–4× latency difference that reflects NUMA topology, not representation cost. When
running per-qubit benchmarks for publication, use `numactl --interleave=all` (Linux)
to distribute pages uniformly across sockets, or `OMP_PLACES=cores OMP_PROC_BIND=close`
when OpenMP threading is active. On macOS with Unified Memory (M-series), NUMA effects
are absent. Document the hardware configuration and memory policy in any published results.

### L5 — IQR outlier filtering not implemented

With 20 runs, IQR-based outlier detection is feasible but not currently implemented in
`bench_compute_stats()`. A stability warning is issued when CV exceeds 50% in console
output. Future work: add IQR trimming (discard runs outside [Q1 - 1.5*IQR, Q3 + 1.5*IQR])
and report trimmed mean alongside the raw mean. The `bench_run_stats_t` struct should
gain `trimmed_mean` and `n_outliers_discarded` fields.

### L6 — Variable iteration counts

Auto-calibration produces different iteration counts for different backends at the same
qubit count. `time_ms` is NOT comparable across backends or targets. Use
`time_per_gate_ms` exclusively for all cross-method and cross-target comparisons.

For all backends (qlib, blas, quest, qulacs, aer_dm), Step H calibrates once per
(qubits, method) pair on target 0 and passes the same `calibrated_iterations` to all
targets. `time_per_gate_ms` values within a method are therefore directly comparable
across targets for every backend.

---

## Risks and Ambiguities

### R1 — `runs_explicit` / `iterations_explicit` flags (Step F)

Resolved: Step F now requires `runs_explicit` and `iterations_explicit` bool fields in
`bench_options_t`. The comparison-against-default heuristic is not used. If `--runs`
appears on the command line, `runs_explicit` is set to 1 regardless of the value, and
the per-qubit default override is skipped.

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
For a fixed target, the state is repeatedly transformed by the same gate, potentially
driving it toward an eigenstate (e.g., repeated X on a qubit in |0⟩ or |1⟩). Since
the state is initialized with random data (not a valid density matrix), this is
unlikely to produce pathological behavior in practice. This is by design: see
§Measurement Semantics. The timed loop measures steady-state performance, not
transient or convergence behavior.

### R5 — Stack size for `perq` array in Step H

Resolved: `perq` is heap-allocated via `calloc` in `perq_main()`. No stack allocation
concerns.

### R6 — `dim_t` type width

`dim_t` is `int64_t` based on the existing `printf` format `%lld` in `bench_print_csv`.
The `bench_result_perq_t` struct uses `dim_t` for the `dim` field and `bench_print_perq_csv`
must use `%lld` (with the `(long long)` cast already used in `bench_print_csv`).
Verify this matches before Step G is merged.

### R7 — Aer-DM API names (Step E2)

Resolved: Step E2 now uses the real `AER::QV::DensityMatrix<double>` API matching
`bench_aer_dm.cpp` exactly — `dm.apply_x`, `dm.apply_unitary_matrix`,
`dm.apply_cnot`, `dm.apply_swap`, `dm.initialize()`, `dm.data()`, `AER::reg_t`,
`AER::cvector_t`, `AER::Utils::vectorize_matrix`. The placeholder names
(`aer_dm_alloc`, `AerDMState`, etc.) have been removed. No further reconciliation
needed at implementation time.
