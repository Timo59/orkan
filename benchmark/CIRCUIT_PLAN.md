# bench_circuit Implementation Plan

Implementation plan for `bench_circuit`: a circuit-level benchmark comparing qlib_packed,
qlib_tiled, blas_dense, QuEST, and Qulacs on three application-motivated quantum circuits in
density matrix mode. This is a separate binary from `bench_gate`; it shares `bench.h`
infrastructure (types, timing helpers, stats) but introduces new circuit abstraction files.

---

## Selected Circuits

### 1. QAOA p=1 on cycle graph (mandatory, VQA thesis)

n-qubit MaxCut QAOA with one layer. Circuit structure:
1. H on all n qubits.
2. For each edge (i, i+1 mod n) of the cycle: RZZ(γ) decomposed as CX · Rz(γ) · CX.
3. Rx(β) on all n qubits.

Total gate count: n + 3n + n = 5n gates. Fixed parameters: γ = π/4, β = π/4 (compile-time
defaults; overridable via `--gamma` and `--beta` CLI flags).

Gate dependencies from `gate.h`: `h`, `rx`, `rz`, `cx` — all implemented.

### 2. Hardware-Efficient VQE Ansatz

d = ⌊n/2⌋ layers of: Ry(π/4) on all n qubits + CX ladder (qubit 0→1, 1→2, ..., n−2→n−1).
Total gate count: d × (2n − 1). Angle θ = π/4 fixed at compile time for reproducibility.

Gate dependencies: `ry`, `cx` — all implemented.

### 3. Quantum Volume (QV) circuit

n layers, each: random permutation of n/2 qubit pairs, each pair receives a random SU(4)
unitary decomposed via KAK: 3 CX + 6 single-qubit Ry/Rz rotations per pair.
Total gate count: approximately 9n²/2 gates. Angles and permutations generated from a
deterministic LCG seeded at 42; stored as static compile-time tables in
`bench_circuit_ops.c` for fixed seed = 42. Only even qubit counts are valid for QV.

Gate dependencies: `ry`, `rz`, `cx` — all implemented.

All circuits use **density matrix mode** for all competitors.

---

## File Map

New files to be created:

```
benchmark/
├── include/
│   └── bench_circuit.h              # circuit_op_t, gate enum, circuit_t, constructor declarations
└── src/
    ├── bench_circuit_ops.c          # circuit_t helpers + all three circuit constructors
    ├── bench_circuit_qlib.c         # qlib packed/tiled circuit runner
    ├── bench_circuit_baselines.c    # blas_dense circuit runner
    ├── bench_circuit_quest.c        # QuEST circuit runner (WITH_QUEST guard)
    ├── bench_circuit_qulacs.cpp     # Qulacs circuit runner (WITH_QULACS guard)
    └── bench_circuit_main.c        # main(), CLI parsing, output formatting, orchestration
```

Modified files:

```
benchmark/CMakeLists.txt             # Add bench_circuit target (Step H)
docs/BENCHMARK.md                    # Update spec (Step I)
~/Projects/thesis/2.Simulation/notes/BENCHMARK_THESIS_NOTES.md   # Thesis notes (Step J)
```

---

## Dependency Graph

```
Step A (bench_circuit.h + bench_circuit_ops.c skeleton)
  |
  +-----+-----+-----+-----+
  |     |     |     |     |
  v     v     v     v     v
Step B  Step C  Step D  Step E  Step F
(circuit  (qlib   (QuEST  (Qulacs  (blas_dense
constructors) runner) runner) runner)  runner)
  |     |     |     |     |
  +-----+-----+-----+-----+
                   |
                   v
              Step G (bench_circuit_main.c)
                   |
                   v
              Step H (CMakeLists.txt)
                   |
                   v
              Step I (docs/BENCHMARK.md)
                   |
                   v
              Step J (BENCHMARK_THESIS_NOTES.md)
```

Steps B, C, D, E, F can all be assigned to independent agents simultaneously after Step A
completes. Steps G, H, I, J form a strict serial chain from G onward because each depends on
the previous step's artifacts being finalized.

---

## Step-by-Step Specification

---

### Step A — Define shared circuit infrastructure

**Files:** `benchmark/include/bench_circuit.h` + `benchmark/src/bench_circuit_ops.c`
**Depends on:** nothing (can start immediately)
**Estimated time:** 1.5 hours

#### bench_circuit.h

```c
#ifndef BENCH_CIRCUIT_H
#define BENCH_CIRCUIT_H

#ifdef __cplusplus
extern "C" {
#endif

/* Gate type identifiers used in circuit_op_t */
typedef enum {
    GATE_X  = 0,
    GATE_Y  = 1,
    GATE_Z  = 2,
    GATE_H  = 3,
    GATE_S  = 4,
    GATE_T  = 5,
    GATE_RX = 6,
    GATE_RY = 7,
    GATE_RZ = 8,
    GATE_CX = 9,
    GATE_CZ = 10,
    GATE_SWAP = 11,
} circuit_gate_t;

/*
 * A single gate operation in a circuit.
 * One-qubit gates use q0; two-qubit gates use q0 (control) and q1 (target).
 * q2 is reserved for future three-qubit gates. angle is used by rotation gates.
 */
typedef struct {
    circuit_gate_t gate_type;
    int q0, q1, q2;
    double angle;
} circuit_op_t;

/* A flat, dynamically-allocated gate sequence */
typedef struct {
    circuit_op_t *ops;       /* Heap-allocated array of gate operations */
    int           n_ops;     /* Number of gates currently in the circuit */
    int           cap;       /* Allocated capacity (internal; do not read) */
    int           n_qubits;  /* Width of the circuit */
    const char   *name;      /* Human-readable circuit name (static string) */
} circuit_t;

/* Lifecycle */
circuit_t circuit_alloc(int n_qubits, int initial_cap, const char *name);
void      circuit_free(circuit_t *c);

/* Append helpers — grow backing array as needed */
void circuit_append_1q(circuit_t *c, circuit_gate_t g, int q0);
void circuit_append_1q_angle(circuit_t *c, circuit_gate_t g, int q0, double angle);
void circuit_append_2q(circuit_t *c, circuit_gate_t g, int q0, int q1);

/* RZZ composite: appends CX(q0,q1) · Rz(q1, gamma) · CX(q0,q1) */
void circuit_rzz(circuit_t *c, int q0, int q1, double gamma);

/* Circuit constructors */
circuit_t circuit_qaoa_p1(int n, double gamma, double beta);
circuit_t circuit_vqe_hea(int n, int depth);     /* depth = -1 → use floor(n/2) */
circuit_t circuit_qv(int n, int seed);           /* n must be even; seed = 42 recommended */

#ifdef __cplusplus
}
#endif

#endif /* BENCH_CIRCUIT_H */
```

#### bench_circuit_ops.c skeleton

Implement `circuit_alloc`, `circuit_free`, `circuit_append_1q`, `circuit_append_1q_angle`,
`circuit_append_2q`, and `circuit_rzz` here. Leave constructor bodies (`circuit_qaoa_p1`,
`circuit_vqe_hea`, `circuit_qv`) as stubs returning an empty circuit — these are filled in
Step B.

Implementation notes:
- `circuit_alloc`: `malloc(initial_cap * sizeof(circuit_op_t))`; set `n_ops = 0`, `cap = initial_cap`.
- `circuit_append_*`: `realloc` when `n_ops == cap` (double capacity).
- `circuit_rzz`: three calls to `circuit_append_*` in order: CX(q0,q1), Rz(q1, gamma), CX(q0,q1).

---

### Step B — Implement circuit constructors

**Files:** `benchmark/src/bench_circuit_ops.c` (fill in the three constructor stubs)
**Depends on:** Step A
**Estimated time:** 2.5 hours

#### `circuit_qaoa_p1(int n, double gamma, double beta)`

```
initial_cap = 5*n
circuit = circuit_alloc(n, initial_cap, "qaoa_p1")
for q in 0..n-1:
    circuit_append_1q(c, GATE_H, q)
for i in 0..n-1:
    q0 = i
    q1 = (i + 1) % n
    circuit_rzz(c, q0, q1, gamma)   /* CX · Rz(gamma) · CX */
for q in 0..n-1:
    circuit_append_1q_angle(c, GATE_RX, q, beta)
return circuit
```

Total ops: n + 3n + n = 5n. Verify: `assert(c.n_ops == 5*n)`.

#### `circuit_vqe_hea(int n, int depth)`

```
if depth < 0: depth = n / 2   /* floor division */
initial_cap = depth * (2*n - 1)
circuit = circuit_alloc(n, initial_cap, "vqe_hea")
for d in 0..depth-1:
    for q in 0..n-1:
        circuit_append_1q_angle(c, GATE_RY, q, M_PI / 4.0)
    for q in 0..n-2:
        circuit_append_2q(c, GATE_CX, q, q+1)
return circuit
```

Total ops: `depth * (n + (n-1)) = depth * (2n-1)`. Verify: `assert(c.n_ops == depth * (2*n-1))`.

#### `circuit_qv(int n, int seed)`

Requirements:
- n must be even (assert or return empty circuit with a warning if odd).
- Use a 32-bit LCG: `state = state * 1664525u + 1013904223u` (Numerical Recipes constants).
- Seed the LCG with the given `seed` value.
- n_pairs = n / 2 per layer; n_layers = n.
- For each layer:
  - Generate a random permutation of qubits 0..n-1 using Fisher-Yates with the LCG.
  - Pair adjacent elements of the permutation: (perm[0], perm[1]), (perm[2], perm[3]), ...
  - For each pair (qa, qb):
    - Append 3 Ry rotations on qa with angles drawn from the LCG.
    - Append 3 Ry rotations on qb with angles drawn from the LCG.
    - Append CX(qa, qb), CX(qb, qa), CX(qa, qb)  [standard KAK CX structure].
    - Append 3 Rz rotations on qa with angles from the LCG.
    - Append 3 Rz rotations on qb with angles from the LCG.
    - Note: exact KAK angle parameterization is flexible; what matters is that all 5
      competitors receive the identical `circuit_op_t` sequence with identical angles.

LCG angle generation: map the raw uint32 to [0, 2π) via `angle = (lcg_next(&state) / 4294967296.0) * 2.0 * M_PI`.

Static compile-time table note: for `seed = 42` and `n` up to 20, the angle table can be
pre-generated and stored as a static array to avoid LCG recomputation per benchmark run.
However, the LCG approach is simpler and fast enough; prefer it unless profiling shows
otherwise. Document whichever approach is chosen.

Implementation risk: the KAK decomposition used here is a structural approximation (CX-Ry-Rz
pattern). A referee may ask whether the 9 single-qubit + 3 CX gates per pair are a valid SU(4)
decomposition. The answer: they are sufficient for benchmarking purposes — all competitors
execute the same sequence; correctness of the SU(4) representation is not required for
performance measurement. Document this in `bench_circuit_main.c`'s header comment.

---

### Step C — qlib circuit runner

**File:** `benchmark/src/bench_circuit_qlib.c`
**Depends on:** Step A
**Estimated time:** 2 hours

#### Function signature

```c
#include "bench.h"
#include "bench_circuit.h"

bench_result_t bench_circuit_qlib(const circuit_t *c, state_type_t type,
                                   int iterations, int warmup, int runs);
```

`type` must be `MIXED_PACKED` or `MIXED_TILED`. Assert or early-return on `PURE`.

#### Structure

```
Allocate state once with state_init(type, n_qubits) — outside timed loop.
Memory: bench_packed_size(n_qubits) or bench_tiled_size(n_qubits).

Build static dispatch table:
  apply_op(state_t *s, const circuit_op_t *op) → switch(op->gate_type):
    GATE_H   → h(s, op->q0)
    GATE_RX  → rx(s, op->q0, op->angle)
    GATE_RY  → ry(s, op->q0, op->angle)
    GATE_RZ  → rz(s, op->q0, op->angle)
    GATE_CX  → cx(s, op->q0, op->q1)
    GATE_X   → x(s, op->q0)
    GATE_Z   → z(s, op->q0)
    default  → (no-op; warn in debug builds only)

Warmup loop (iterations = warmup):
  Re-initialize state to maximally mixed (I/dim) before warmup.
  Apply all c->n_ops operations once per warmup iteration.

Timed loop:
  for r in 0..runs-1:
    Re-initialize state to maximally mixed (I/dim).
    start = bench_time_ns()
    for i in 0..iterations-1:
      for op in 0..c->n_ops-1:
        apply_op(&state, &c->ops[op])
    run_times[r] = bench_ns_to_ms(bench_time_ns() - start)

state_free(&state)
stats = bench_compute_stats(run_times, runs)
```

Populate `bench_result_t`:
- `gate_name = c->name`
- `method = (type == MIXED_PACKED) ? "qlib_packed" : "qlib_tiled"`
- `iterations = iterations`
- `sweep_size = c->n_ops`  ← total gate count is the natural sweep unit for circuit benchmarks
- `ops_per_sec = (iterations * c->n_ops) / (time_ms / 1000.0)`
- `time_per_gate_ms = time_ms / ((double)iterations * c->n_ops)`
- `memory_bytes` = packed or tiled size as above

Maximally mixed state initialization: set all diagonal elements to `1.0/dim` and all
off-diagonals to 0. Use `state_set(state, i, i, 1.0/dim)` in a loop over i = 0..dim-1, or
call `memset(state.data, 0, ...)` then set diagonals. Do NOT use `state_plus` — that sets
all elements to `1/dim`, not only the diagonal.

Note: The plan spec says "re-initialize to maximally mixed state (all diagonal entries =
1/dim)". Confirm that this is the intended initial condition (it is the maximally mixed
density matrix ρ = I/dim). Cross-check against QuEST's `initMixedState` or equivalent.

---

### Step D — QuEST circuit runner

**File:** `benchmark/src/bench_circuit_quest.c`
**Depends on:** Step A; requires `WITH_QUEST` at compile time
**Estimated time:** 2 hours

#### Function signature (inside `#ifdef WITH_QUEST`)

```c
bench_result_t bench_circuit_quest(const circuit_t *c, int iterations, int warmup, int runs);
```

#### Structure

Pattern mirrors `bench_quest.c` exactly:
- Requires `bench_quest_init()` / `bench_quest_cleanup()` called by `bench_circuit_main.c`.
- Creates `Qureg rho = createDensityQureg(c->n_qubits)` once before the loop.
- Re-initializes to maximally mixed state before each run.

QuEST maximally mixed initialization: call `initMixedState(rho)` if available in QuEST v4,
otherwise set elements manually via `setDensityAmps`. Investigate the QuEST v4 API — check
`extern/QuEST/quest/include/quest.h` for the correct function name.

Gate dispatch function (inside static helper):

```c
static void apply_circuit_op_quest(Qureg rho, const circuit_op_t *op) {
    switch (op->gate_type) {
        case GATE_H:   applyHadamard(rho, op->q0); break;
        case GATE_RX:  applyRx(rho, op->q0, op->angle); break;
        case GATE_RY:  applyRy(rho, op->q0, op->angle); break;
        case GATE_RZ:  applyRz(rho, op->q0, op->angle); break;
        case GATE_CX:  applyControlledPauliX(rho, op->q0, op->q1); break;
        case GATE_X:   applyPauliX(rho, op->q0); break;
        case GATE_Z:   applyPauliZ(rho, op->q0); break;
        default: break;
    }
}
```

Verify QuEST v4 API names for `applyRx`, `applyRy`, `applyRz` — these may differ from v3.
Check `extern/QuEST/quest/include/quest.h` before implementation.

Memory: `(size_t)dim * dim * sizeof(double) * 2` (matches `bench_quest.c` formula).
`gate_name = c->name`, `method = "quest"`.

---

### Step E — Qulacs circuit runner

**File:** `benchmark/src/bench_circuit_qulacs.cpp`
**Depends on:** Step A; requires `WITH_QULACS` at compile time
**Estimated time:** 2 hours

#### Function signature (inside `#ifdef WITH_QULACS`, `extern "C"` block)

```cpp
bench_result_t bench_circuit_qulacs(const circuit_t *c, int iterations, int warmup, int runs);
```

#### Structure

Pattern mirrors `bench_qulacs.cpp` exactly:
- Allocate `CTYPE *rho` once (full dim×dim density matrix).
- Re-initialize to maximally mixed state before each run: set diagonal elements to
  `CTYPE(1.0/dim, 0.0)` and all off-diagonals to `CTYPE(0.0, 0.0)`.

Gate dispatch:

```cpp
static void apply_circuit_op_qulacs(const circuit_op_t *op, CTYPE *rho, ITYPE dim) {
    switch (op->gate_type) {
        case GATE_H:   dm_H_gate(op->q0, rho, dim); break;
        case GATE_RX:  dm_RX_gate(op->q0, op->angle, rho, dim); break;
        case GATE_RY:  dm_RY_gate(op->q0, op->angle, rho, dim); break;
        case GATE_RZ:  dm_RZ_gate(op->q0, op->angle, rho, dim); break;
        case GATE_CX:  dm_CNOT_gate(op->q0, op->q1, rho, dim); break;
        case GATE_X:   dm_X_gate(op->q0, rho, dim); break;
        case GATE_Z:   dm_Z_gate(op->q0, rho, dim); break;
        default: break;
    }
}
```

Verify Qulacs csim API names for `dm_RX_gate`, `dm_RY_gate`, `dm_RZ_gate` — check
`extern/qulacs/src/csim/update_ops_dm.hpp` for exact signatures. The rotation gates likely
accept `(UINT target, double angle, CTYPE *state, ITYPE dim)` but confirm parameter order.

Memory: `static_cast<size_t>(dim) * dim * sizeof(CTYPE)` (matches `bench_qulacs.cpp`).
`gate_name = c->name`, `method = "qulacs"`.

---

### Step F — blas_dense circuit runner

**File:** `benchmark/src/bench_circuit_baselines.c`
**Depends on:** Step A
**Estimated time:** 2.5 hours

#### Function signature

```c
bench_result_t bench_circuit_blas_dense(const circuit_t *c, int iterations,
                                         int warmup, int runs);
```

#### Structure

This runner applies each gate as a full N×N unitary acting on a full N×N density matrix via
UρU†, using two `zgemm` calls per gate operation.

Outside the timed loop:
1. Allocate `rho` (dim×dim complex double, full density matrix).
2. For each distinct gate in `c->ops`, pre-compute the full N×N gate matrix using the
   tensor-product expansion (embed the 1Q or 2Q gate matrix into the N×N space).
   - Store as an array of pointers indexed by op index OR keyed on (gate_type, q0, q1, angle).
   - Simpler: pre-compute one N×N matrix per unique `circuit_op_t` by iterating ops and
     building a map. For circuits with repeating ops (QAOA, VQE), this deduplicates work.
   - Alternative (simpler to implement but slower): compute each gate matrix on the fly
     inside the timed loop. This is acceptable if the gate-matrix construction time is
     amortized over `iterations`.
   - **Recommended approach:** pre-compute ALL op matrices (one per position in `c->ops`)
     outside the timed loop. Memory cost: `c->n_ops × dim² × sizeof(cplx_t)`. For n=10,
     n_ops ≈ 50, dim = 1024: 50 × 1M × 16B ≈ 800 MB. This is prohibitive. Use the
     deduplication strategy: build a table keyed on (gate_type, q0, q1, angle_quantized),
     or fall back to on-the-fly construction inside the timed loop for correctness.
   - **Practical recommendation:** pre-compute one matrix per unique (gate_type, qubit tuple)
     ignoring angle (use a canonical angle for the key), or pre-compute and re-use by
     constructing the N×N matrix once per op at init time but storing only the 2×2 or 4×4
     kernel, then expanding in-place just before `zgemm`. This keeps memory bounded.
   - The simplest correct implementation: build the N×N unitary for each op position once
     before the timing loop and store in `op_mats[n_ops]`. Cap the dense runner at a qubit
     limit (e.g. n ≤ 8) matching the per-gate baseline behaviour.

Inside the timed loop:
```
for i in 0..iterations-1:
    for k in 0..c->n_ops-1:
        /* Apply rho <- U_k * rho * U_k† via two zgemm calls */
        zgemm("N", "N", dim, dim, dim, 1, op_mats[k], dim, rho, dim, 0, tmp, dim)  /* tmp = U*rho */
        zgemm("N", "C", dim, dim, dim, 1, tmp, dim, op_mats[k], dim, 0, rho, dim)  /* rho = tmp*U† */
```

Re-initialize `rho` to maximally mixed (diagonal = 1/dim, off-diagonal = 0) before each run.

Memory reported: `dim² × sizeof(cplx_t)` (state only; auxiliary buffers excluded).

Qubit cap: enforce `n <= CIRCUIT_DENSE_MAX_QUBITS` (define as 8 in `bench_circuit.h` or
`bench_circuit_baselines.c`). Return a zero `bench_result_t` with a warning to stderr if
exceeded.

Risks: N×N matrix storage for `c->n_ops` distinct gates may exceed available RAM at n > 8.
The `n <= 8` cap prevents this. At n = 8: `op_mats` = 40 ops × 256² × 16B = 40 MB,
acceptable.

---

### Step G — Write `benchmark/src/bench_circuit_main.c`

**File:** `benchmark/src/bench_circuit_main.c`
**Depends on:** Steps A, B, C, D, E, F (all runners and constructors complete)
**Estimated time:** 3 hours

This is the most architecturally complex step because it owns the orchestration logic,
output formatting, and CLI interface for a multi-circuit, multi-method benchmark.

#### Extended options struct

Define a new `circuit_options_t` in this file (or extend `bench_options_t` with a wrapper):

```c
typedef struct {
    bench_options_t base;    /* min_qubits, max_qubits, step, iterations, warmup, runs, csv, pgfplots, verbose */
    int  run_qaoa;           /* 1 if "all" or "qaoa" selected */
    int  run_vqe;            /* 1 if "all" or "vqe" selected */
    int  run_qv;             /* 1 if "all" or "qv" selected */
    int  seed;               /* --seed N, default 42 */
    double gamma;            /* --gamma FLOAT, default M_PI/4 */
    double beta;             /* --beta FLOAT, default M_PI/4 */
} circuit_options_t;
```

#### CLI flags (superset of bench_gate flags)

| Flag | Default | Description |
|------|---------|-------------|
| `--min-qubits N` | 2 | Minimum qubit count |
| `--max-qubits N` | 12 | Maximum qubit count |
| `--step N` | 1 | Qubit count increment |
| `--iterations N` | 100 | Circuit runs per timed window (lower default than bench_gate because circuits are larger) |
| `--warmup N` | 10 | Untimed warm-up iterations |
| `--runs N` | 10 | Independent timing runs |
| `--csv` | off | CSV output to stdout |
| `--pgfplots` | off | pgfplots-compatible .dat output |
| `--verbose` | off | Extended per-result info |
| `--circuit NAME` | all | Select circuit: `qaoa`, `vqe`, `qv`, or `all` |
| `--seed N` | 42 | RNG seed for QV circuit |
| `--gamma FLOAT` | π/4 | QAOA phase separation parameter |
| `--beta FLOAT` | π/4 | QAOA mixing parameter |
| `--help` | — | Print usage and exit |

Validation:
- `--circuit` must be one of `{qaoa, vqe, qv, all}`; otherwise exit(1).
- `--seed` must be ≥ 0.
- `--gamma`, `--beta`: any finite double.
- QV requires even qubit counts; skip odd qubit counts with a warning when `run_qv = 1`.

#### Output formats

CSV: add a `circuit_name` column as the first field after `qubits`:

```
circuit,qubits,dim,method,runs,iterations,n_ops,time_ms_mean,time_ms_std,time_ms_min,
time_ms_median,time_ms_cv,time_per_gate_ms,ops_per_sec,memory_bytes
```

Console: for each `(circuit, qubits)` combination, print a block:

```
Circuit: qaoa_p1  qubits: 6  n_ops: 30
  qlib_packed    12.345 ± 0.123 ms (CV  1.2%)   1234567 ops/sec
  qlib_tiled      9.876 ± 0.089 ms (CV  0.9%)   1543210 ops/sec
  blas_dense    456.789 ± 5.678 ms (CV  1.4%)     23456 ops/sec
  quest          34.567 ± 0.456 ms (CV  1.3%)    123456 ops/sec
  qulacs         18.234 ± 0.234 ms (CV  1.2%)    234567 ops/sec
  Speedup: 1.25x tiled vs packed, 37.0x vs blas_dense ...
  Memory: packed=XX KB, tiled=XX KB, dense=XX MB, quest=XX MB, qulacs=XX MB
─────────────────────────────────────────────────────────────
```

pgfplots: emit one table per circuit per statistic, with columns for each method. Follow the
same structure as `bench_gate`'s `print_pgfplots_output()`: five tables per circuit (mean,
ci95, min, median, cv) plus one memory table. Units: µs per gate (time_ms × 1000 / n_ops).

Note: `sweep_size` in `bench_result_t` is repurposed as `n_ops` for circuit benchmarks —
it represents the number of gate applications per iteration. This is consistent with the
field's semantics in `bench_gate` (sweep_size = gates applied per iteration).

#### main() outline

```c
int main(int argc, char *argv[]) {
    circuit_options_t opts = circuit_parse_options(argc, argv);
    /* Emit header */
#ifdef WITH_QUEST
    bench_quest_init();
#endif
    for qubits in opts.base.min_qubits .. opts.base.max_qubits step opts.base.step:
        if opts.run_qaoa:
            circuit_t c = circuit_qaoa_p1(qubits, opts.gamma, opts.beta)
            /* run all methods; collect results; emit output */
            circuit_free(&c)
        if opts.run_vqe:
            circuit_t c = circuit_vqe_hea(qubits, -1)
            /* run all methods; collect results; emit output */
            circuit_free(&c)
        if opts.run_qv:
            if qubits % 2 != 0: continue (warn)
            circuit_t c = circuit_qv(qubits, opts.seed)
            /* run all methods; collect results; emit output */
            circuit_free(&c)
#ifdef WITH_QUEST
    bench_quest_cleanup();
#endif
    return 0;
}
```

Reuse `bench_quest_init` and `bench_quest_cleanup` from `bench_quest.c` — these are already
declared in `bench.h` under `#ifdef WITH_QUEST`. No duplication needed.

---

### Step H — Add `bench_circuit` target to `benchmark/CMakeLists.txt`

**File:** `benchmark/CMakeLists.txt`
**Depends on:** Step G (all source files must exist before CMake target is written)
**Estimated time:** 1 hour

Append after the existing `bench_gate` target block (after the status summary):

```cmake
#=======================================================================================================================
# bench_circuit - Circuit-level benchmark
#=======================================================================================================================

set(BENCH_CIRCUIT_SRC
        "${BENCH_SRC_DIR}/bench_circuit_ops.c"
        "${BENCH_SRC_DIR}/bench_circuit_qlib.c"
        "${BENCH_SRC_DIR}/bench_circuit_baselines.c"
        "${BENCH_SRC_DIR}/bench_circuit_main.c"
)

if(QuEST_FOUND)
    list(APPEND BENCH_CIRCUIT_SRC "${BENCH_SRC_DIR}/bench_circuit_quest.c")
endif()

if(QULACS_FOUND)
    list(APPEND BENCH_CIRCUIT_SRC "${BENCH_SRC_DIR}/bench_circuit_qulacs.cpp")
endif()

add_executable(bench_circuit ${BENCH_CIRCUIT_SRC})

set_target_properties(bench_circuit PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/benchmark"
    EXCLUDE_FROM_ALL TRUE
)

target_include_directories(bench_circuit
        PRIVATE
            "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_link_libraries(bench_circuit
        PRIVATE
            q
)

if(QuEST_FOUND)
    target_compile_definitions(bench_circuit PRIVATE WITH_QUEST)
    target_include_directories(bench_circuit PRIVATE
        ${QuEST_INCLUDE_DIR}
        ${EXTERN_DIR}/QuEST/build
    )
    target_link_libraries(bench_circuit PRIVATE ${QuEST_LIBRARY})
endif()

if(QULACS_FOUND)
    target_compile_definitions(bench_circuit PRIVATE WITH_QULACS)
    target_include_directories(bench_circuit PRIVATE
        ${QULACS_INCLUDE_DIR}
    )
    target_link_libraries(bench_circuit PRIVATE csim_static)
    set_property(TARGET bench_circuit PROPERTY CXX_STANDARD 17)
    set_property(TARGET bench_circuit PROPERTY CXX_STANDARD_REQUIRED ON)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(bench_circuit PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()
```

Update the status summary block to add a circuit benchmark line:

```cmake
message(STATUS "  bench_circuit: enabled (circuit-level benchmarks)")
```

Note: `bench_circuit` links `bench_mixed.c` (for `build_all_pairs`) indirectly through `q`
— wait, `build_all_pairs` is declared in `bench.h` and defined in `bench_mixed.c`, which is
**not** part of `bench_circuit`'s source list. Two options:

1. Move `build_all_pairs` to a shared `bench_utils.c` compiled into both targets.
2. Duplicate the definition in `bench_circuit_ops.c` (simplest for now).
3. Add `bench_mixed.c` to `BENCH_CIRCUIT_SRC` (causes duplicate `main` — impossible).

**Recommended:** add a new file `benchmark/src/bench_utils.c` containing `build_all_pairs`,
and add it to both `BENCH_MIXED_SRC` and `BENCH_CIRCUIT_SRC`, removing it from
`bench_mixed.c`. Alternatively, define it `static` in each file that needs it. For
`bench_circuit`, `build_all_pairs` may not be needed at all (circuits use fixed qubit pairs
from the circuit_op_t sequence). Confirm whether any circuit runner needs it before deciding.

---

### Step I — Update `docs/BENCHMARK.md`

**File:** `docs/BENCHMARK.md`
**Depends on:** Steps G, H (implementation complete; CLI and output format finalized)
**Estimated time:** 1 hour

Changes required:

1. **File Layout section**: add `bench_circuit_ops.c`, `bench_circuit_qlib.c`,
   `bench_circuit_baselines.c`, `bench_circuit_main.c`, `bench_circuit_quest.c`,
   `bench_circuit_qulacs.cpp`, and `bench_circuit.h` to the tree.

2. **CMake targets table**: add `bench_circuit | Executable | Always (Release only)`.

3. **CLI Reference**: add a new subsection `### bench_circuit` with the full flag table
   (all bench_gate flags plus `--circuit`, `--seed`, `--gamma`, `--beta`), default values,
   valid ranges, and validation behaviour.

4. **Circuits Benchmarked**: new section describing the three circuits, their gate counts,
   parameter conventions, and qubit-range recommendations (see runtime estimates below).

5. **Output Formats**: add `circuit_name` column to CSV description; document the
   circuit-name column position (first column) and that `sweep_size` = `n_ops`.

6. **Status section**: update to reflect `bench_circuit` is implemented.

7. **Intro paragraph**: remove "Circuit-level benchmarking is out of scope; a future
   `bench_circuit` executable is reserved for that purpose." Replace with a sentence
   pointing to `bench_circuit`.

---

### Step J — Update thesis notes

**File:** `~/Projects/thesis/2.Simulation/notes/BENCHMARK_THESIS_NOTES.md`
**Depends on:** Steps G, H, I (full implementation and spec finalized)
**Estimated time:** 1.5 hours

Add a new section (e.g., `## Circuit-Level Benchmarks (bench_circuit)`) covering:

1. **Circuit selection and scientific motivation**
   - QAOA p=1: primary VQA thesis result; demonstrates end-to-end algorithm performance on
     the exact circuit structure used in experiments.
   - VQE hardware-efficient ansatz: continuous-parameter variational circuit; tests
     sustained rotation-gate throughput rather than Clifford-only workloads.
   - Quantum Volume: generalisation beyond VQA; tests random deep circuits with SU(4)
     unitaries; the QV metric is widely used in industry benchmarking.

2. **Parameter reproducibility**
   - QAOA: γ = β = π/4 fixed compile-time defaults. Values are not optimised; they are
     chosen to ensure all rotation gates are non-trivial (angle ≠ 0, π/2, π).
     The thesis should note that these are benchmark angles, not variational optima.
   - QV: seed = 42; angles and permutations are deterministic across all competitors and
     all runs. The LCG formula and constants must be stated in the thesis for full
     reproducibility (LCG: state = state × 1664525 + 1013904223, mod 2³²).

3. **Cross-competitor fairness**
   - All competitors receive the identical `circuit_op_t` sequence (same gate types, same
     qubit indices, same angles).
   - All competitors initialize to the same maximally mixed state ρ = I/dim before each run.
   - All competitors use density matrix mode.
   - Timing covers only the gate application loop; state allocation and re-initialization
     are outside the timed region.

4. **Thesis claims supported by each circuit**
   - QAOA → VQA performance: absolute throughput and speedup vs. QuEST/Qulacs on the exact
     circuit structure used in VQA experiments supports the claim that qlib is competitive
     for VQA simulation workloads.
   - VQE → continuous-parameter variational: demonstrates that qlib's rotation-gate
     implementation (rx, ry, rz with arbitrary angle) is efficient at scale.
   - QV → generalisation: a strong QV result shows that performance advantages are not
     specific to VQA circuit structure but extend to arbitrary unitary circuits.

5. **Runtime estimates and recommended qubit ranges**

   | Circuit | Recommended range | Notes |
   |---------|------------------|-------|
   | QAOA    | 4–16 qubits      | 5n gates; fast even at n=16 (80 gates). Extend to 20 if runtimes permit. |
   | VQE HEA | 4–16 qubits      | d(2n−1) gates; d=⌊n/2⌋. At n=16: d=8, 248 gates. |
   | QV      | 4–14 qubits (even only) | ~9n²/2 gates. At n=14: ~882 gates. Dense baseline capped at n≤8. |

   blas_dense runner is capped at n ≤ 8 for all circuits due to O(dim²) memory per gate
   matrix. For thesis plots, run blas_dense only in the 4–8 qubit range and qlib/QuEST/
   Qulacs from 4–16 (QAOA, VQE) or 4–14 (QV).

---

## Dependency Summary Table

| Step | Description | Depends on | Parallelizable with |
|------|-------------|------------|---------------------|
| A | Header + ops skeleton | — | (start point) |
| B | Circuit constructors | A | C, D, E, F |
| C | qlib runner | A | B, D, E, F |
| D | QuEST runner | A | B, C, E, F |
| E | Qulacs runner | A | B, C, D, F |
| F | blas_dense runner | A | B, C, D, E |
| G | main + CLI + output | A, B, C, D, E, F | — |
| H | CMakeLists.txt update | G | — |
| I | docs/BENCHMARK.md update | G, H | — |
| J | Thesis notes update | G, H, I | — |

Maximum parallelism: 5 agents can work simultaneously on Steps B through F after Step A
completes. The total critical path is: A → (B|C|D|E|F) → G → H → I → J.

---

## Time Estimates

| Step | Time |
|------|------|
| A | 1.5 h |
| B | 2.5 h |
| C | 2.0 h |
| D | 2.0 h |
| E | 2.0 h |
| F | 2.5 h |
| G | 3.0 h |
| H | 1.0 h |
| I | 1.0 h |
| J | 1.5 h |
| **Total (serial)** | **19.0 h** |
| **Critical path** | **A + max(B..F) + G + H + I + J = 1.5 + 2.5 + 3.0 + 1.0 + 1.0 + 1.5 = 10.5 h** |

With 5 agents in parallel on Steps B–F, wall-clock time is approximately 10.5 hours.

---

## Risks and Open Questions

1. **QuEST v4 API for rotation gates and maximally mixed initialization.** The names
   `applyRx`, `applyRy`, `applyRz`, and `initMixedState` are assumed but must be confirmed
   against `extern/QuEST/quest/include/quest.h` before starting Step D. QuEST v4 changed
   several API names relative to v3.

2. **Qulacs csim rotation gate signatures.** `dm_RX_gate`, `dm_RY_gate`, `dm_RZ_gate` must
   be verified against `extern/qulacs/src/csim/update_ops_dm.hpp`. Parameter order
   (especially whether angle comes before or after the qubit index) must be confirmed.

3. **blas_dense N×N gate matrix construction.** The tensor-product embedding of a 2×2 or
   4×4 gate matrix into an N×N unitary is non-trivial to implement correctly and efficiently.
   Confirm that the existing `bench_baselines.c` code can be referenced as a model. The
   `bench_blas_dense` function there already constructs full N×N matrices for 1Q gates.

4. **`build_all_pairs` ownership.** Currently defined in `bench_mixed.c` and declared in
   `bench.h`. Step H's CMake target must not link `bench_mixed.c`. Resolve by adding
   `bench_utils.c` (recommended) or by confirming that `bench_circuit`'s runners do not
   need `build_all_pairs` at all (likely true, since circuit ops already specify qubit indices).

5. **QV qubit-count restriction.** QV requires even qubit counts. The main loop must skip
   odd qubit counts when `run_qv = 1` and emit a warning. This is a user experience concern
   but must be implemented correctly in Step G.

6. **KAK decomposition validity.** The 9 single-qubit + 3 CX gates per pair used in the QV
   constructor are structurally correct for SU(4) coverage but are not a rigorous KAK
   decomposition with computed interaction angles. For benchmarking, this is acceptable.
   Document explicitly in the source file header comment to pre-empt reviewer questions.

7. **Memory pressure at large qubit counts.** For blas_dense at n = 8: one full 256×256
   density matrix = 512 KB; `n_ops` pre-computed gate matrices at n=8 for QV ≈ 882 ops ×
   512 KB ≈ 440 MB. Pre-computing all gate matrices for QV will exceed available RAM even
   at the n ≤ 8 cap. Implement on-the-fly gate matrix construction inside the timed loop
   for the blas_dense runner, or pre-compute only a small set (deduplicated by gate type
   and qubit indices). The time cost of matrix construction inside the timed loop will
   inflate the blas_dense numbers, but this affects all qubit counts equally and the runner
   is only a lower-bound reference anyway.

8. **Default iteration counts.** `bench_gate` defaults to `--iterations 1000`. For circuit
   benchmarks, a full circuit application is 5n to ~9n² gate calls. Recommend defaulting
   `bench_circuit` to `--iterations 100` (not 1000) to keep wall-clock time manageable at
   large qubit counts. Document this difference.
