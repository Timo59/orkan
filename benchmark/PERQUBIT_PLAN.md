# Per-Qubit Benchmark Extension — Implementation Plan (v3)

Target binary: `bench_gate` (built from `benchmark/src/`)
New flag: `--per-qubit`

---

## Measurement Semantics

**What `--per-qubit` measures**: steady-state single-target throughput.

Each `_at` function applies the same gate to the same fixed target qubit many times on
a single state that is **re-initialized before each target's warmup phase**. After warmup,
the branch predictor and prefetcher are saturated. This is a deliberate isolation
measurement: it exposes whether gate cost depends on target position in the memory
layout of each representation. It is **not** circuit-representative throughput.

Consequences:
- Results are **not** directly comparable to the sweep benchmark's `time_per_gate_ms`.
  The sweep benchmark is a better proxy for circuit-like workloads.
- Results **are** appropriate for comparing representations head-to-head at each target
  position and for identifying memory-layout asymmetries.
- Document this distinction whenever citing per-qubit numbers in a thesis or paper.

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
       mixed.c  baseline quest.c  qulacs.cpp
            │        │        │        │
            └────────┴────┬───┴────────┘
                          │
         Step E2 ─────────┤   ← Aer-DM adapter; depends only on A
       bench_aer_dm.cpp   │
                          │
         Step F ──────────┤   ← CLI parsing; depends only on A
         Step G ──────────┤   ← output functions; depends only on A
                          ▼
               ┌─────────────────────┐
               │   Step H: main()    │  ← depends on B, C, D, E, E2, F, G
               └──────────┬──────────┘
                           │
               ┌───────────┴───────────┐
               ▼                       ▼
            Step I                  Step J
         docs/BENCHMARK.md      thesis notes
```

**Parallel batches:**
- Batch 0: Step A
- Batch 1: Steps B, C, D, E, E2, F, G (all depend only on A)
- Batch 2: Step H
- Batch 3: Step I
- Batch 4: Step J

---

## Platform and Performance Requirements

### CPU Frequency Scaling

On Linux with `schedutil` or `powersave` governor, early benchmark runs execute at
reduced frequency; sustained load causes turbo boost to decay after ~100ms. Either:
- (preferred) Set the performance governor before running: `cpupower frequency-set -g
  performance` (requires root). Document the requirement in the README.
- Or use time-based warmup (see §Warmup Strategy below), which automatically absorbs
  frequency ramp-up by running until at least 100ms of wall time has elapsed.

On macOS with Apple Silicon (M-series), the OS manages frequency automatically and
this is not a concern. On Intel macOS, `sudo powermetrics -s cpu_power` can confirm
clock frequency during the run.

**Document the CPU governor and clock frequency in any published results.**

### Huge Pages (10+ Qubits)

At 12 qubits, a density matrix is 4 GB containing ~1M 4K pages, far exceeding L1 TLB
capacity (~64 entries). TLB miss cost dominates gate application time at 10+ qubits.

Mitigations:
- **Linux**: `madvise(MADV_HUGEPAGE)` immediately after `malloc` for state buffers ≥ 1 GB.
  Alternatively, allocate with `mmap(MAP_ANONYMOUS | MAP_HUGETLB)` (requires `hugepages`
  kernel support). Add `bench_alloc_huge(size_t n)` in `bench.h` that wraps this.
- **macOS**: The VM subsystem handles large allocations with 2 MB pages automatically
  for `mmap(MAP_ANONYMOUS)` mappings above a threshold (~4 MB). No explicit action needed.

**Log the actual huge-page mapping in verbose mode** by checking
`/proc/PID/smaps` (Linux) or `vm_region_info` (macOS) after state allocation. This
documents whether the mitigation is active for a given run.

### Timing Serialization

`clock_gettime(CLOCK_MONOTONIC)` reads the hardware time counter, but on ARM64 this
reads `CNTVCT_EL0` via vDSO which is **not serializing** — the CPU can reorder it
relative to adjacent instructions.

Add a platform-specific serialization barrier to `bench.h`:

```c
static inline void bench_timing_barrier(void) {
#if defined(__aarch64__)
    __asm__ volatile("isb" ::: "memory");
#elif defined(__x86_64__)
    __asm__ volatile("lfence" ::: "memory");
#else
    __asm__ volatile("" ::: "memory");  /* compiler barrier only */
#endif
}
```

Use `bench_timing_barrier()` immediately before each `bench_time_ns()` call in
timing loops. Note on ARM: `bench_timing_barrier()` after `t0 = bench_time_ns()`
prevents the CPU from speculatively executing the first gate call before the timestamp
is committed. The volatile re-read of `h->call` provides a data dependency, but the
explicit barrier makes this guarantee architectural rather than implementation-dependent.

---

## Single Timing Harness

**The old plan copy-pasted the timing loop and calibration probe across 12+ functions.
This is replaced by a single parameterized harness.**

### Callback Interface

Define in `bench.h`:

```c
/**
 * @brief Callback invoked by bench_run_timed() to apply one gate operation.
 *
 * The callback must be stored in a volatile function pointer to prevent the
 * compiler from seeing through the indirection and hoisting loop body across
 * iterations. See bench_run_timed() for usage.
 *
 * @param ctx  Backend-specific context (state pointer, qubit indices, etc.)
 */
typedef void (*bench_gate_cb_t)(void *ctx);

/**
 * @brief Context for the timing harness, parameterizing one (gate, target) combination.
 *
 * All fields must be fully populated before calling bench_run_timed().
 */
typedef struct {
    volatile bench_gate_cb_t call;  /**< Opaque callback — volatile prevents devirtualization */
    void    *ctx;                   /**< Backend context passed to call() */
    int      iterations;            /**< Gate applications per timed run (> 0) */
    int      runs;                  /**< Number of independent timed runs (1..BENCH_MAX_RUNS) */
} bench_harness_t;
```

### `bench_run_timed()`

Define in `bench.h` as `static inline`:

```c
/**
 * @brief Run the timing harness: warmup, then runs × iterations gate applications.
 *
 * The callback is stored as a volatile function pointer in bench_harness_t, which
 * prevents the compiler from devirtualizing the call even with LTO. bench_timing_barrier()
 * prevents both compiler and CPU reordering of the time reads relative to the gate
 * application.
 *
 * Warmup: runs until at least BENCH_PQ_WARMUP_MS wall time has elapsed (time-based,
 * not count-based). This absorbs CPU frequency ramp-up and TLB cold-start.
 *
 * Warmup is batched to avoid clock_gettime overhead at small qubit counts: at 2 qubits
 * gates take ~10ns each, so calling clock_gettime every iteration would cost more than
 * the gate itself. Batch size scales with qubit count.
 *
 * @param h        Populated harness descriptor
 * @param out      Array of `h->runs` doubles to receive per-run times in ms
 * @param qubits   Qubit count, used to select warmup batch size
 * @return         0 on success
 */
static inline int bench_run_timed(const bench_harness_t *h, double *out, int qubits) {
    /* Batch size for warmup: 128 calls per clock check at small qubits, scaling down
     * to 1 at large qubits where a single gate already takes many milliseconds. */
    int warmup_batch = (qubits <= 4) ? 128 : (qubits <= 7) ? 16 : (qubits <= 10) ? 4 : 1;

    uint64_t wt0 = bench_time_ns();
    do {
        for (int wb = 0; wb < warmup_batch; ++wb)
            h->call(h->ctx);
    } while (bench_ns_to_ms(bench_time_ns() - wt0) < BENCH_PQ_WARMUP_MS);

    for (int r = 0; r < h->runs; ++r) {
        bench_timing_barrier();
        uint64_t t0 = bench_time_ns();
        bench_timing_barrier();   /* prevent first gate from executing before t0 committed */
        for (int i = 0; i < h->iterations; ++i) {
            h->call(h->ctx);  /* volatile fn ptr — not inlined, not hoisted */
        }
        bench_timing_barrier();
        out[r] = bench_ns_to_ms(bench_time_ns() - t0);
    }
    return 0;
}
```

**Why volatile function pointer prevents hoisting**: The compiler cannot prove that
`h->call` is the same callable across iterations (it is declared volatile, so re-reading
it each iteration is required by the standard). This makes each iteration a black-box
call from the optimizer's perspective, preventing CSE and loop-body hoisting. Note:
`volatile` on a function pointer is a pragmatic (not standards-guaranteed) devirtualization
barrier; `__attribute__((noinline))` on callbacks provides an additional defense layer.

### Calibration

**There is one calibration mechanism: an external probe in Step H.** The `_at`
functions do not auto-calibrate; they always receive a positive `iterations` count.
The `iterations=0` sentinel is eliminated.

The probe in Step H:
1. Runs the harness with `iterations = BENCH_PQ_PROBE_ITERS` and `runs = BENCH_PQ_PROBE_RUNS`.
2. Takes the **minimum** run time across probe runs (not a single cold-cache run).
3. Computes `calibrated_iterations = ceil(MIN_MEASUREMENT_MS / (probe_min_ms / BENCH_PQ_PROBE_ITERS))`.
4. Clamps to `[1, BENCH_PQ_MAX_CALIBRATED_ITERATIONS]`.

Taking the minimum instead of a single run eliminates systematic cold-cache bias
(a single probe run can be inflated 10x by a context switch).

**Adaptive probe cap**: At 12 qubits a single gate application takes ~100ms. Running
`BENCH_PQ_PROBE_ITERS = 100` iterations × `BENCH_PQ_PROBE_RUNS = 5` runs would take
~50 seconds before any measurement. If a single probe iteration exceeds
`MIN_MEASUREMENT_MS`, cap `BENCH_PQ_PROBE_ITERS` to 1 for that configuration:

```c
/* Before probe: quick single-shot to detect very slow gates */
bench_timing_barrier();
uint64_t probe_single_t0 = bench_time_ns();
bench_timing_barrier();   /* prevent first gate from executing before t0 committed */
h->call(h->ctx);
bench_timing_barrier();
double single_iter_ms = bench_ns_to_ms(bench_time_ns() - probe_single_t0);
int probe_iters = (single_iter_ms >= MIN_MEASUREMENT_MS) ? 1 : BENCH_PQ_PROBE_ITERS;
probe_harness.iterations = probe_iters;
```

```c
/* Calibration constants */
#define BENCH_PQ_PROBE_ITERS              100   /**< iterations in each probe run */
#define BENCH_PQ_PROBE_RUNS               5     /**< runs for calibration probe; take minimum */
#define BENCH_PQ_WARMUP_MS               100.0  /**< minimum warmup duration in ms */
#define MIN_MEASUREMENT_MS                1.0   /**< minimum timed run duration to target */
#define BENCH_PQ_FALLBACK_ITERATIONS      10000 /**< when probe rounds to 0 ns */
#define BENCH_PQ_MAX_CALIBRATED_ITERATIONS 1000000
```

---

## Data Structures

### `bench_result_perq_t`

Per-qubit result struct. **No bandwidth column**: the `theoretical_bw_gbs` and
`bw_notes` fields from the old plan are dropped. The formula is wrong by different
factors for every backend, and a companion annotation string to interpret it adds
complexity for minimal value. Report only directly measured quantities.

```c
typedef struct bench_result_perq {
    qubit_t    qubits;
    dim_t      dim;
    const char *gate_name;
    const char *method;
    qubit_t    target;           /**< Target qubit (q1 for 2Q gates) */
    qubit_t    target2;          /**< Second qubit for 2Q gates; QUBIT_NONE for 1Q */
    int        is_2q;
    double     time_ms_mean;
    double     time_ms_std;      /**< Sample stddev (Bessel-corrected); honest n=15 */
    double     time_ms_min;
    double     time_ms_median;
    double     time_ms_cv;
    double     ops_per_sec;
    double     time_per_gate_us;        /**< mean time per gate application (microseconds) */
    double     time_per_gate_us_median;
    double     time_per_gate_us_min;
    size_t     memory_bytes;
    int        iterations;
    int        runs;
} bench_result_perq_t;
```

Units: all per-gate times in **microseconds** to avoid sub-millisecond precision loss
in the CSV. `time_ms_*` fields remain in milliseconds (total run time).

**Statistical note**: with the default of 15 runs, Bessel-corrected stddev has a 95%
CI of approximately [0.74σ, 1.83σ]. Report min/median/mean and note n in any citation.
Do not derive confidence intervals from CV without acknowledging this limitation.

### Array Layout

```c
bench_result_perq_t perq[NUM_METHODS][MAX_TARGETS]
```

Method-major layout so that inner loops over targets for a single method are contiguous.
This is the inverse of the old plan's `perq[MAX_TARGETS][NUM_METHODS]`.

When calling `bench_print_perq_console()`, pass `perq[m]` directly — it is a contiguous
array of `MAX_TARGETS` elements. No copy to a flat array needed.

### Constants

```c
#define BENCH_PQ_DEFAULT_RUNS          15    /**< 15 gives honest CI; old plan had 5 */
#define BENCH_PQ_DEFAULT_ITERATIONS   100    /**< overridden by calibration */
#define MAX_TARGETS                   (2 * MAX_PAIRS)  /* 256 */
#define QUBIT_NONE                    ((qubit_t)-1)    /**< max value of qubit_t regardless of width */
```

Note: `BENCH_MAX_RUNS` (defined elsewhere in bench.h) is used directly wherever a
maximum runs bound is needed. No separate alias is introduced.

### Random Number Generation

**Do not use `srand`/`rand`**: they are global state, not thread-safe, and not portable
across libcs.

For state initialization in C files, use a fast PRNG seeded from the OS entropy source.
The fill does **not** need cryptographic quality — it only needs to avoid pathological
structure (e.g., all-zeros, which can trigger SIMD fast-paths that skew results).

```c
/* Seed xoshiro256** from 8 bytes of OS entropy, then fill buf with pseudo-random bytes.
 * On macOS/BSD: arc4random_buf provides the seed bytes directly.
 * On Linux: getrandom(2) (glibc >= 2.25, #include <sys/random.h>) provides them.
 *
 * Operates on raw uint8_t* bytes only — never touches cplx_t directly — so it is
 * safe to call from both C and C++ translation units even though _Complex is a C99 type.
 */
static inline void bench_fill_random(void *buf, size_t n) {
    /* Seed one uint64_t from OS entropy */
    uint64_t seed;
#if defined(__APPLE__) || defined(__FreeBSD__)
    arc4random_buf(&seed, sizeof(seed));
#else
    /* Linux: requires #include <sys/random.h> (glibc >= 2.25) */
    getrandom(&seed, sizeof(seed), 0);
#endif

    /* Expand seed via splitmix64 into 4 x uint64 xoshiro256** state */
    uint64_t z;
    z = (seed += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; uint64_t s0 = z ^ (z >> 31);
    z = (seed += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; uint64_t s1 = z ^ (z >> 31);
    z = (seed += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; uint64_t s2 = z ^ (z >> 31);
    z = (seed += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; uint64_t s3 = z ^ (z >> 31);

    /* xoshiro256** bulk fill — 8 bytes per step */
    uint8_t *p = (uint8_t *)buf;
    size_t remaining = n;
    while (remaining >= 8) {
        /* xoshiro256** output word */
        uint64_t r = s1 * 5; r = ((r << 7) | (r >> 57)) * 9;
        /* write 8 bytes directly via memcpy to avoid strict-aliasing UB */
        memcpy(p, &r, 8); p += 8; remaining -= 8;
        /* xoshiro256** advance */
        uint64_t t = s1 << 17;
        s2 ^= s0; s3 ^= s1; s1 ^= s2; s0 ^= s3; s2 ^= t;
        s3 = (s3 << 45) | (s3 >> 19);
    }
    if (remaining > 0) {
        uint64_t r = s1 * 5; r = ((r << 7) | (r >> 57)) * 9;
        memcpy(p, &r, remaining);  /* tail bytes */
    }
}
```

The key requirement: seed from `arc4random_buf` (macOS/BSD) or `getrandom(2)` (Linux).
Do **not** use `fread` from `/dev/urandom` for large buffers: at 4 GB (12-qubit density
matrix) this would take minutes. The PRNG fill runs in microseconds.

- **C++ files** (bench_qulacs.cpp, bench_aer_dm.cpp): use `std::mt19937_64` seeded
  from `std::random_device{}()`. The existing bench_aer_dm.cpp already does this
  correctly with a fixed seed derived from `qubits`. Keep that pattern.

The old plan's `srand(42)` is removed from all C backends. Each `_at` call re-seeds
by filling the state buffer with `bench_fill_random`. Within-backend reproducibility
is no longer required; cross-backend state comparison was already invalid.

---

## Step A — `bench.h` additions

**File:** `benchmark/include/bench.h`

Add in order after the existing constants block:
1. `bench_timing_barrier()` inline function (platform-specific, see §Timing Serialization)
2. `bench_fill_random()` static inline helper (see §Random Number Generation). Note:
   `bench_fill_random` operates on raw `uint8_t*` bytes only and never on `cplx_t`,
   making it safe to call from C++ translation units even though `_Complex` is a C99 type.
3. Constants block: `BENCH_PQ_DEFAULT_RUNS=15`, `BENCH_PQ_DEFAULT_ITERATIONS=100`,
   `BENCH_PQ_PROBE_ITERS=100`, `BENCH_PQ_PROBE_RUNS=5`,
   `BENCH_PQ_WARMUP_MS=100.0`, `MIN_MEASUREMENT_MS=1.0`,
   `BENCH_PQ_FALLBACK_ITERATIONS=10000`, `BENCH_PQ_MAX_CALIBRATED_ITERATIONS=1000000`,
   `MAX_TARGETS=(2*MAX_PAIRS)`, `QUBIT_NONE=((qubit_t)-1)`,
   `MAX_QUBITS_FOR_ORDERED_PAIRS=16`,
   `BENCH_HUGEPAGE_THRESHOLD=((size_t)1 << 30)`
4. `bench_gate_cb_t` typedef and `bench_harness_t` struct
5. `bench_run_timed()` static inline function
6. `bench_fill_perq_stats()` static inline helper
7. `bench_result_perq_t` struct
8. `bench_options_t` additions: `per_qubit`, `runs_explicit`, `iterations_explicit`,
   `gate_filter`
9. `_at` adapter declarations (see §Adapter Function Signatures)
10. `build_all_ordered_pairs` declaration
11. Per-qubit output function declarations

Do NOT modify any existing declarations or definitions.

### Adapter Function Signatures

All `_at` functions receive:
- A pre-allocated, pre-initialized state pointer specific to the backend.
- A positive `iterations` count (calibrated by Step H; no auto-calibration inside `_at`).
- `bench_run_timed()` performs time-based warmup internally; no `warmup` parameter.

```c
/* qlib backends (bench_mixed.c) */
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

/* BLAS baseline (bench_baselines.c) */
bench_result_perq_t bench_blas_dense_at(qubit_t qubits, const char *gate_name,
    const cplx_t gate_mat[4],
    qubit_t target, int iterations, int runs, cplx_t *rho);

bench_result_perq_t bench_blas_dense_2q_at(qubit_t qubits, const char *gate_name,
    const cplx_t gate_mat[16],
    qubit_t q1, qubit_t q2, int iterations, int runs, cplx_t *rho);

#ifdef WITH_QUEST
bench_result_perq_t bench_quest_at(qubit_t qubits, const char *gate_name,
    qubit_t target, int iterations, int runs, Qureg *qureg);

bench_result_perq_t bench_quest_2q_at(qubit_t qubits, const char *gate_name,
    qubit_t q1, qubit_t q2, int iterations, int runs, Qureg *qureg);
#endif

#ifdef WITH_QULACS
bench_result_perq_t bench_qulacs_at(qubit_t qubits, const char *gate_name,
    qubit_t target, int iterations, int runs, void *rho);

bench_result_perq_t bench_qulacs_2q_at(qubit_t qubits, const char *gate_name,
    qubit_t q1, qubit_t q2, int iterations, int runs, void *rho);
#endif

#ifdef WITH_AER_DM
bench_result_perq_t bench_aer_dm_at(qubit_t qubits, const char *gate_name,
    qubit_t target, int iterations, int runs, void *dm);

bench_result_perq_t bench_aer_dm_2q_at(qubit_t qubits, const char *gate_name,
    qubit_t q1, qubit_t q2, int iterations, int runs, void *dm);
#endif
```

### `bench_fill_perq_stats` Helper

Inline helper that populates stats fields of `bench_result_perq_t` from a
`bench_run_stats_t`. Convert `time_ms_{mean,min,median}` → `time_per_gate_us` by
dividing by `iterations` and multiplying by 1000. Compute `ops_per_sec` from mean.
This eliminates copy-paste across all backends.

---

## Step B — `bench_mixed.c` additions

**File:** `benchmark/src/bench_mixed.c`

### State Management

The `_at` functions receive a pre-allocated `state_t*`. They **re-initialize the state**
before each target's warmup (call `state_reinit_random()` described below), then call
`bench_run_timed()`.

Add a file-local helper:
```c
/* Re-fill state->data with random values without re-allocating.
 * Uses bench_fill_random() from bench.h. */
static void state_reinit_random(state_t *state);
```

This is the key difference from the sweep functions: the state is re-initialized per
target, not once globally. After thousands of gate applications to the same target,
the state can be driven toward eigenstates, changing measured cost for backends with
data-dependent SIMD fast-paths. Re-initializing prevents this.

### Callback Structs

Each `_at` function wraps its gate call in a minimal context struct and callback:

```c
typedef struct { state_t *state; void (*fn)(state_t*, qubit_t); qubit_t target; }
    cb1q_ctx_t;
static void cb1q(void *ctx) {
    cb1q_ctx_t *c = ctx; c->fn(c->state, c->target);
}
```

The callback address is stored in `bench_harness_t.call` as a volatile function
pointer. Two-qubit variants extend the struct with `q1`/`q2`.

### `_at` Function Structure

Each `_at` function:
1. Validates inputs; fills fixed `bench_result_perq_t` fields.
2. Calls `state_reinit_random(state)`.
3. Builds a `cb1q_ctx_t` (or `cb2q_ctx_t`) and a `bench_harness_t`.
4. Declares `double run_times[BENCH_MAX_RUNS]`; asserts `runs <= BENCH_MAX_RUNS`;
   calls `bench_run_timed()`.
5. Calls `bench_compute_stats()` then `bench_fill_perq_stats()`.
6. Returns result.

Four public functions: `bench_qlib_packed_at`, `bench_qlib_tiled_at`,
`bench_qlib_packed_2q_at`, `bench_qlib_tiled_2q_at`. Each is a thin wrapper filling
`qlib_bench_desc_t` and calling a shared `run_1q_at()` or `run_2q_at()` internal.

Also add `build_all_ordered_pairs()` immediately after `build_all_pairs()`. Same
structure: nested loop over all `i != j` pairs. Returns 0 if `qubits >
MAX_QUBITS_FOR_ORDERED_PAIRS`. Declared in `bench.h`.

### Constraints

- Do NOT modify `run_1q_timing()`, `run_2q_timing()`, or any existing public functions.
- `bw_notes` field removed; no bandwidth computation.
- No `srand(42)`: use `bench_fill_random()` from bench.h.

---

## Step C — `bench_baselines.c` additions

**File:** `benchmark/src/bench_baselines.c`

Same pattern as Step B. Two `_at` functions: `bench_blas_dense_at` and
`bench_blas_dense_2q_at`.

Key difference from qlib: the callback invokes `apply_uruh_blas(dim, U, rho, tmp)`.
`U` and `tmp` are allocated **before** the harness call and freed after. `rho` is
re-initialized before each target's warmup via a call to `bench_fill_random()`.

`U` is built once per `_at` call (before warmup) via `build_full_gate_matrix()` /
`build_full_gate_matrix_2q()`. It is NOT rebuilt per run or per iteration.

The callback context holds `(dim_t dim, cplx_t *U, cplx_t *rho, cplx_t *tmp)`.

`run_times` is a stack array `double run_times[BENCH_MAX_RUNS]` with
`assert(runs <= BENCH_MAX_RUNS)`.

The dense baseline is capped at 8 qubits in the Step H orchestration loop, not in
the `_at` functions themselves.

---

## Step D — `bench_quest.c` additions

**File:** `benchmark/src/bench_quest.c`

Both new functions inside `#ifdef WITH_QUEST`.

**State re-init**: QuEST's `initDebugState(qureg)` before each target's warmup. This
produces a fixed deterministic state; it ignores randomness entirely. This is an
accepted limitation for cross-backend state consistency (it was already a known
limitation; re-initializing before each target at least makes it consistent within the
QuEST backend across targets).

**Gate dispatch**: resolve `gate_name` to a `void (*)(Qureg, int)` or
`void (*)(Qureg, int, int)` function pointer once before warmup. Store in callback
context. The callback invokes it through a struct field (opaque to the harness).

**Important**: `applyControlledPauliX(Qureg, int control, int target)` — verify
argument order in Step D matches the sweep path (q1 = control, q2 = target for CX).

`run_times` is a stack array `double run_times[BENCH_MAX_RUNS]` with
`assert(runs <= BENCH_MAX_RUNS)`.

---

## Step E — `bench_qulacs.cpp` additions

**File:** `benchmark/src/bench_qulacs.cpp`

Both new functions inside `#ifdef WITH_QULACS` and `extern "C"`.

**C/C++ boundary**: `bench_qulacs_at` and `bench_qulacs_2q_at` are C-callable wrappers
(defined in this `.cpp` file, declared in `bench.h` inside `extern "C"`). All
Qulacs-specific types (`CTYPE`, `ITYPE`, `UINT`, csim function pointers) are confined
to this file. No Qulacs types appear in `bench.h`.

**State re-init**: call `bench_fill_random()` on the `CTYPE*` buffer before each
target's warmup. Replace `std::srand(42)` / `std::rand()` with `std::mt19937_64`
seeded from `std::random_device{}()` each time `state_reinit` is needed, or use a
simple fill via `bench_fill_random` (C helper) since the buffer is plain bytes.

**Gate dispatch**: resolve `gate_name` to `void (*)(UINT, CTYPE*, ITYPE)` or
`void (*)(UINT, UINT, CTYPE*, ITYPE)` once before warmup. The callback context holds
the function pointer, `rho`, `dim`, and target qubit(s).

The callback type wraps the Qulacs csim signature. `bench_run_timed()` sees only
`void (*)(void*)`.

---

## Step E2 — `bench_aer_dm.cpp` additions

**File:** `benchmark/src/bench_aer_dm.cpp`

Both new functions inside `#ifdef WITH_AER_DM` and `extern "C"`.

**C/C++ boundary**: same pattern as Step E. `AER::QV::DensityMatrix<double>` and all
Aer types are invisible outside this file.

**State re-init**: call `dm.initialize()` and then refill `dm.data()` with random
values (using `std::mt19937_64` as in the existing `bench_aer_dm()` function) before
each target's warmup.

**Gate dispatch**: resolve to `gate_kind_t` enum (BK_X, BK_H, BK_Z, BK_CX,
BK_SWAP, BK_UNKNOWN) once before warmup. The callback context holds the
`DensityMatrix*`, `gate_kind`, precomputed `h_vec`/`z_vec`, and qubit indices.
The switch inside the callback dispatches to the correct Aer method.

Wrap the entire function body in `try { ... } catch (...)` as in `bench_aer_dm()`.

Method string: `"qiskit_aer_dm"` to match the sweep benchmark.

---

## Step F — `bench_main.c`: CLI parsing

**File:** `benchmark/src/bench_main.c`

### Changes to `bench_options_t` (Step A declares; Step F implements)

New fields: `per_qubit`, `runs_explicit`, `iterations_explicit`, `gate_filter`.

### `bench_parse_options()` changes

1. Zero-initialize new fields in the designated initializer.
2. On `--runs`: set `opts.runs_explicit = 1`.
3. On `--iterations`: set `opts.iterations_explicit = 1`. Reject `--iterations 0`
   at parse time (prevents accidental trigger of any future auto-calibrate code).
4. Add `--per-qubit` branch: sets `opts.per_qubit = 1`.
5. Add `--gate NAME` branch: stores `argv[++i]` in `opts.gate_filter`.
6. Post-parse defaults: if `opts.per_qubit` and not explicit, set `opts.runs =
   BENCH_PQ_DEFAULT_RUNS` and `opts.iterations = BENCH_PQ_DEFAULT_ITERATIONS`.
7. Mutual exclusion: if `opts.per_qubit && opts.pgfplots_output`, print error and exit.
8. Update `print_usage()` with `--per-qubit` and `--gate` entries.

**Note**: `--gate` without `--per-qubit` has no effect in the sweep path. Print a
warning when `opts.gate_filter != NULL && !opts.per_qubit`.

**`--gate` matching**: use `strcasecmp` (POSIX) at all call sites, not `strcmp`.

---

## Step G — `bench_main.c`: per-qubit output functions

**File:** `benchmark/src/bench_main.c`

Three new functions, declared in `bench.h` (Step A), implemented here.

### `bench_print_perq_csv_header()`

Column order:
```
qubits,dim,gate,method,runs,iterations,target,target2,is_2q,
time_ms_mean,time_ms_std,time_ms_min,time_ms_median,time_ms_cv,
time_per_gate_us,time_per_gate_us_median,time_per_gate_us_min,
ops_per_sec,memory_bytes
```

No bandwidth column. `time_per_gate_us` in microseconds (not ms) for numeric precision.

### `bench_print_perq_csv()`

Skip rows where `time_ms_mean <= 0.0` (failed measurements). `target2`: print `-1`
when `target2 == QUBIT_NONE` (1Q gates), not the raw sentinel value `255` — downstream
analysis scripts must not mistake the sentinel for a valid qubit index. Use a helper:
```c
/* Emit -1 for 1Q gates; actual index for 2Q gates */
if (r->target2 == QUBIT_NONE) fprintf(f, ",-1");
else                           fprintf(f, ",%u", (unsigned)r->target2);
```
Use `%lld` with `(long long)` cast for `dim`.

### `bench_print_perq_console()`

**1Q gates**: one row per target. Columns: `target`, `mean(us/gate)`,
`median(us/gate)`, `min(us/gate)`, `CV(%)`, `ops/sec`. Append `(noisy)` when
`cv > 50%`. Skip rows where `time_ms_mean <= 0.0`.

**2Q gates**: summary line showing min/mean/max `time_per_gate_us` across pairs with
the pair coordinates at min and max. Full per-pair data is in the CSV. Since
`perq[m]` is now a contiguous array (method-major layout), pass `perq[m]` directly;
no copy needed.

---

## Step H — `bench_main.c`: per-qubit orchestration

**File:** `benchmark/src/bench_main.c`

Implement `static int perq_main(const bench_options_t *opts)`. Call from `main()` as
an early return when `opts.per_qubit` is set.

### State Re-initialization Policy

**Per-target re-init is mandatory.** Before each target's `_at` call, the state must
be re-initialized to a fresh random state. This prevents gates from driving the state
toward eigenstates across targets.

The `_at` functions handle this internally (see Steps B–E2). Step H does not need to
call any separate reinit function; it simply calls `_at` for each target in sequence,
and each `_at` call re-initializes before its own warmup.

**Do not** share a single state across targets within the same method loop. The
method-first loop order is used solely to avoid repeated allocations (allocate once
per method, pass the same buffer to each `_at` call), not to share state content
across targets.

### Result Storage

```c
bench_result_perq_t (*perq)[MAX_TARGETS] = calloc(NUM_METHODS, sizeof(*perq));
```

Layout: `perq[method][target]` — method-major. Access pattern for output:
```c
bench_print_perq_console(perq[m], n_targets);  /* contiguous, no copy */
```

Heap-allocated once per `perq_main()` call. Total size:
`NUM_METHODS × MAX_TARGETS × sizeof(bench_result_perq_t)` ≈ 7 × 256 × 120 bytes ≈ 215 KB.

### Calibration (external probe)

For each (gate, qubits, method) before the target loop:

```
/* Quick single-shot to detect very slow gates (e.g. 12-qubit DM) */
single_iter_ms = time one call to h->call(h->ctx)
probe_iters = (single_iter_ms >= MIN_MEASUREMENT_MS) ? 1 : BENCH_PQ_PROBE_ITERS

probe_harness.iterations = probe_iters
probe_harness.runs       = BENCH_PQ_PROBE_RUNS
run bench_run_timed(&probe_harness, probe_times, qubits)
probe_min_ms = min(probe_times[0..BENCH_PQ_PROBE_RUNS-1])
probe_per_iter_ms = probe_min_ms / probe_iters
if probe_per_iter_ms <= 0:
    calibrated_iterations = BENCH_PQ_FALLBACK_ITERATIONS
else:
    calibrated_iterations = ceil(MIN_MEASUREMENT_MS / probe_per_iter_ms)
    clamp to [1, BENCH_PQ_MAX_CALIBRATED_ITERATIONS]
```

If `opts->iterations_explicit`, skip calibration and use `opts->iterations` directly.

The probe runs on the initial target 0 with a freshly initialized state. The probe
itself counts as warmup; it does not corrupt the calibrated measurement since each
`_at` call re-initializes before its own warmup.

### `perq_main()` Structure

```
allocate perq[NUM_METHODS][MAX_TARGETS] via calloc
emit CSV header or banner

#ifdef WITH_QUEST: bench_quest_init()

results_count = 0

for qubits in [min_qubits .. max_qubits step step]:
    dim = 1 << qubits
    run_dense = (qubits <= 8)

    /* 1Q gates */
    for each gate g in gates[]:
        if gate_filter && strcasecmp(gate_filter, g.name) != 0: continue
        results_count++
        n_targets = qubits
        memset(perq, 0, NUM_METHODS * sizeof(*perq))

        /* method-first loop: one allocation per (qubits, method) */
        for each active method m:
            allocate state buffer for method m (use bench_alloc_huge for >= BENCH_HUGEPAGE_THRESHOLD)
            initialize state buffer with bench_fill_random / appropriate backend init
            calibrate to get calibrated_iterations (probe on target 0)
            for target = 0 to qubits-1:
                perq[m][target] = backend_at(..., calibrated_iterations, runs, state)
                /* _at re-inits state before warmup; state buffer ownership stays here */
            free/destroy state buffer

        output(perq, n_targets)

    /* 2Q gates */
    build q1s/q2s arrays (ordered or unordered depending on gate symmetry)
    for each gate g in gates_2q_perq[]:
        ... same method-first pattern, with n_pairs instead of n_targets ...

#ifdef WITH_QUEST: bench_quest_cleanup()

if results_count == 0 && gate_filter != NULL: warn and return EXIT_FAILURE
free(perq)
return 0
```

### `gate_def_2q_perq_t` — Separate Struct for Per-Qubit Mode

Define a separate struct for 2Q gates in per-qubit mode to avoid modifying existing
sweep-path `gate_def_2q_t` struct and its initializers:

```c
typedef struct {
    const char *name;
    void (*fn)(state_t*, qubit_t, qubit_t);
    const cplx_t *mat;
    int has_tiled;
    int is_symmetric;  /* 1 = SWAP (unordered pairs); 0 = CX (ordered pairs) */
} gate_def_2q_perq_t;
```

Use `build_all_pairs` for symmetric gates; `build_all_ordered_pairs` for asymmetric
gates. Check `n_pairs <= 0 || n_pairs > MAX_TARGETS` after either builder and skip
with a warning if out of range.

### Huge Page Allocation

Add a helper `bench_alloc_huge(size_t n)` to `bench.h`:
- If `n >= BENCH_HUGEPAGE_THRESHOLD` on Linux: allocate with
  `mmap(MAP_ANONYMOUS|MAP_PRIVATE)` then `madvise(MADV_HUGEPAGE)`.
- Otherwise: plain `malloc`.
- On macOS: `mmap(MAP_ANONYMOUS|MAP_PRIVATE)` (kernel handles huge pages).
- Always pair with a corresponding `bench_free_huge(void*, size_t)`.

`BENCH_HUGEPAGE_THRESHOLD` is `(size_t)1 << 30` (1 GB). Huge pages are only
worthwhile for very large allocations; the per-page overhead of 2 MB pages is not
worth it below ~100 MB.

### C/C++ Compilation Boundary

`perq_main()` is in `bench_main.c` (a `.c` file). It calls:
- `bench_qulacs_at`, `bench_qulacs_2q_at` — defined in `bench_qulacs.cpp`, C-callable
  via `extern "C"` in that file.
- `bench_aer_dm_at`, `bench_aer_dm_2q_at` — same pattern in `bench_aer_dm.cpp`.

**No C++ types or headers appear in `bench_main.c`.** The Qulacs state buffer is
typed as `void*` in the C side; the `.cpp` file casts it to `CTYPE*` internally.
The Aer-DM object is allocated inside `bench_aer_dm.cpp` via a C-callable allocator:

```c
/* In bench_aer_dm.cpp, inside extern "C": */
void *bench_aer_dm_alloc(qubit_t qubits);      /* returns AER::QV::DensityMatrix<double>* */
void  bench_aer_dm_reinit(void *dm, qubit_t qubits); /* initialize() + random fill */
void  bench_aer_dm_free(void *dm);
```

These three functions are declared in `bench.h` inside `#ifdef WITH_AER_DM` and
implemented in `bench_aer_dm.cpp`. `perq_main()` calls them to manage Aer-DM state
without any C++ types crossing the translation boundary.

Qulacs state buffers are plain `cplx_t*` arrays (or `void*` aliased to `CTYPE*`);
allocation can happen in `bench_main.c` using `bench_alloc_huge`. Re-initialization
is handled inside the `_at` functions.

---

## Step I — `docs/BENCHMARK.md` update

Add or update:
1. CLI reference: `--per-qubit` and `--gate` flags with defaults and mutual exclusions.
2. Constants table: new `BENCH_PQ_*` constants.
3. Key Data Structures: `bench_result_perq_t` field documentation. Note absence of
   bandwidth column and reason (accuracy varies too much per backend).
4. CSV format: per-qubit CSV column order, units (microseconds for per-gate times).
5. Warmup strategy: document time-based warmup (`BENCH_PQ_WARMUP_MS = 100ms`),
   batched warmup loop, and CPU governor requirement.
6. Huge pages: document threshold (`BENCH_HUGEPAGE_THRESHOLD = 1 GB`) and
   platform-specific behavior.
7. Timing serialization: document `bench_timing_barrier()` and platform differences.
8. Adapter function listing: all `_at` signatures including Aer-DM.
9. Statistical note: 15 runs, report min/median/mean; CI not claimed.
10. Known limitations (see §Known Limitations below).

---

## Step J — Thesis Notes

**File:** `~/Projects/thesis/2.Simulation/notes/BENCHMARK_THESIS_NOTES.md`

Add `## Per-Qubit Analysis` section covering:
1. What it measures and what it does not (steady-state, not circuit throughput).
2. Interpreting target-position differences in packed vs. tiled layouts.
3. State re-initialization rationale (eigenstate bias, data-dependent SIMD fast-paths).
4. Platform requirements for valid results (CPU governor, huge pages).
5. Citing per-qubit data in pgfplots: `x=target, y=time_per_gate_us` for
   target-position sweep plots.
6. Runtime estimates at default settings (15 runs, qubits 2–12, all methods).

---

## Known Limitations

**L1 — State initialization differs between backends.** qlib and BLAS use
`bench_fill_random()`; Qulacs and Aer-DM use `std::mt19937_64`. QuEST's
`initDebugState` produces a deterministic state unrelated to the others. At low qubit
counts where the matrix fits in L1/L2 cache, different data values can affect
SIMD fast-paths. Accept this for cross-framework comparison; document it when citing.

**L2 — Measurement semantics vs. circuit workloads.** See §Measurement Semantics.
Do not compare per-qubit `time_per_gate_us` directly to sweep `time_per_gate_ms`
without qualification.

**L3 — Bandwidth not reported.** The `theoretical_bw_gbs` field is removed because
the formula `2 * memory_bytes / time_per_gate` is wrong by different factors (2×, 4/3×,
2.5× underestimate) for each backend. If bandwidth analysis is needed, measure it
with STREAM on the same hardware and compare independently.

**L4 — NUMA and thread affinity.** On multi-socket machines, target-qubit cost
differences may reflect NUMA topology, not representation cost. Use `numactl
--interleave=all` on Linux for publication results. Not applicable on Apple Silicon.

**L5 — No IQR outlier filtering.** With 15 runs, IQR trimming is feasible but not
implemented. A `(noisy)` flag at CV > 50% is the only outlier indicator.

**L6 — Variable iteration counts across methods.** Step H calibrates once per
(gate, qubits, method). Different methods will have different `iterations` counts.
Use `time_per_gate_us` (not total run time) for all cross-method comparisons.

**L7 — Calibration-on-target-0 bias.** The probe runs on target 0 only. For packed
storage representations, gates on higher-index targets touch more of the state buffer
(the stride pattern means higher qubits access a wider spread of memory), so the
calibrated iteration count may be too high for high targets. Runs for targets near
`qubits-1` may exceed `MIN_MEASUREMENT_MS` by 2–5× in packed storage due to the
increased working set. This inflates wall time but does not affect `time_per_gate_us`
accuracy. Accept this as a known limitation; do not re-calibrate per target.

---

## Time Estimates

| Step | Description | Estimate |
|------|-------------|----------|
| A | bench.h additions | 2 h |
| B | bench_mixed.c additions | 2 h |
| C | bench_baselines.c additions | 1.5 h |
| D | bench_quest.c additions | 1.5 h |
| E | bench_qulacs.cpp additions | 1.5 h |
| E2 | bench_aer_dm.cpp additions | 2 h |
| F | bench_main.c CLI parsing | 1 h |
| G | bench_main.c output functions | 1.5 h |
| H | bench_main.c orchestration | 3 h |
| I | docs/BENCHMARK.md update | 1 h |
| J | Thesis notes | 1 h |
| **Total sequential** | | **18 h** |
| **Critical path** (A→Batch1→H→I→J) | | **9 h** |

### Runtime Estimates (wall clock)

Default: `--per-qubit --min-qubits 2 --max-qubits 12`, 15 runs × calibrated iterations,
3 gates (X, H, CX), 6 methods.

| Qubit range | Notes | Rough total |
|-------------|-------|-------------|
| 2–8 | States fit in cache/RAM; fast | < 10 min |
| 9–12 | 4 GB states; TLB-miss dominated | 20–45 min |
| **Total** | | **~30–55 min** |

Time-based warmup (100 ms per target) adds approximately
`n_targets × n_methods × 100ms` in warmup overhead:
- At 12 qubits: 12 targets × 6 methods × 0.1s ≈ 7s per gate type.
- For 3 gate types × 11 qubit counts: ≈ 4 minutes additional warmup total.
This is small relative to gate application time at 12 qubits.
