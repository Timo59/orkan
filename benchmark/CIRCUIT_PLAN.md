# bench_circuit Implementation Plan

Implementation plan for `bench_circuit`: a circuit-level benchmark comparing qlib_packed,
qlib_tiled, QuEST, Qulacs, and Aer-DM on three application-motivated quantum circuits in
density matrix mode. This is a separate binary from `bench_gate`; it shares
`bench.h` infrastructure (types, timing helpers, stats) but introduces new circuit abstraction
files.

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
**Note:** The VQE angle (π/4) is intentionally not CLI-overridable. It is hardcoded for
reproducibility — a fixed, non-trivial angle ensures all rotation gates are non-trivial while
keeping results deterministic across runs without requiring any seed. This is not an oversight.

Gate dependencies: `ry`, `cx` — all implemented.

### 3. Quantum Volume (QV) circuit

n layers, each: random permutation of n/2 qubit pairs, each pair receives a random
parameterized two-qubit circuit decomposed as: 3 Ry rotations per qubit + CX(a,b)·CX(b,a)·CX(a,b)
+ 3 Rz rotations per qubit. This is a standard gate template exercising the same gate *mix*
(9 single-qubit rotations + 3 two-qubit gates per pair) as a KAK decomposition. **Do not
call it a "KAK decomposition"** — the interaction angles are random, not KAK-computed, so the
resulting unitary is not a general SU(4) from KAK theory. Call it the "random parameterized
two-qubit circuit template" in all source comments and documentation.

Total gate count: n layers × n/2 pairs × (6 Ry + 3 CX + 6 Rz) = 15n²/2 gates (approximately).
Concretely: n=8 has 8 layers × 4 pairs × 15 ops/pair = 480 gate ops. Angles and permutations
generated from a deterministic 32-bit LCG seeded at 42. Only even qubit counts are valid for QV.

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
    ├── bench_circuit_quest.c        # QuEST circuit runner (WITH_QUEST guard)
    ├── bench_circuit_qulacs.cpp     # Qulacs circuit runner (WITH_QULACS guard)
    ├── bench_circuit_aer_dm.cpp     # Aer-DM circuit runner (WITH_AER_DM guard)
    └── bench_circuit_main.c         # main(), CLI parsing, output formatting, orchestration
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
  +-----+-----+-----+-----+-----+
  |     |     |     |     |     |
  v     v     v     v     v     v
Step B  Step C  Step D  Step E  Step F
(circuit  (qlib   (QuEST  (Qulacs  (Aer-DM
constructors) runner) runner) runner)  runner)
  |     |     |     |     |     |
  +-----+-----+-----+-----+-----+
                   |
                   v
              Step H (bench_circuit_main.c)
                   |
                   v
              Step I (CMakeLists.txt)
                   |
                   v
              Step J (docs/BENCHMARK.md)
                   |
                   v
              Step K (BENCHMARK_THESIS_NOTES.md)
```

Steps B–F can all be assigned to independent agents simultaneously after Step A completes.
Steps H, I, J, K form a strict serial chain from H onward.

---

## Implementation Strategy: Parallel Agent Execution

This section defines how a coordinator can maximally parallelize implementation using multiple
agents working on independent files simultaneously.

### Dependency Graph (agent view)

```
Wave 1 (sequential):  Agent 1 → Step A (bench_circuit.h + bench_circuit_ops.c skeleton)
                                              |
                         +--------------------+--------------------+
                         |          |          |          |        |
Wave 2 (parallel):   Agent 2    Agent 3    Agent 4    Agent 5   Agent 6
                     Step B     Step C     Step D     Step E    Step F
                   (circuit   (qlib      (QuEST     (Qulacs   (Aer-DM
                  constructors) runner)   runner)   runner)   runner)
                         |          |          |          |        |
                         +--------------------+--------------------+
                                              |
Wave 3 (sequential):              Agent → Step H (bench_circuit_main.c)
                                              |
                                        Step I (CMakeLists.txt) — can be
                                        drafted in parallel with Wave 2,
                                        finalized after Step H
                                              |
                                        Steps J, K (docs)
```

Steps B–F depend only on Step A and are fully independent of each other — each touches a
distinct file and has no shared mutable state.

### Wave 1: Step A (unblocks everything)

**Assigned to:** Agent 1
**Output:** `benchmark/include/bench_circuit.h` + `benchmark/src/bench_circuit_ops.c` skeleton
**Unblocks:** Steps B, C, D, E, F, and the CMake draft

Agent 1 must complete and commit (or hand off via a shared plan) the `circuit_t`,
`circuit_op_t`, `circuit_gate_t` struct/enum definitions and the lifecycle/append function
implementations before any Wave 2 agent starts. The `bench_circuit.h` header is the shared
contract: all Wave 2 agents read it but do not modify it.

### Wave 2: Five agents in parallel (Steps B–F)

Launch all five agents simultaneously once Step A's `bench_circuit.h` is finalized.

Each agent receives in its prompt:
1. The complete `bench_circuit.h` content (copy verbatim from Step A output) — the agent must
   not re-read or modify it.
2. The `runner_fn` typedef from Step H:
   ```c
   typedef void (*runner_fn)(const circuit_t *c, void *ctx, int iterations);
   ```
3. The handoff contract (see below): the exact function signatures the agent must export.
4. The backend-specific API context (QuEST headers path, Qulacs header path, Aer-DM header
   path, or qlib gate function signatures) so the agent can work without reading unrelated
   files.
5. The relevant step specification from this plan (Step B, C, D, E, or F respectively).

**Agent 2 — Step B (circuit constructors):** needs the three circuit specifications from this
plan (QAOA, VQE, QV) and the `circuit_append_*` / `circuit_rzz` function signatures from
`bench_circuit.h`. No external backend knowledge required.

**Agent 3 — Step C (qlib packed/tiled runner):** needs the qlib gate function signatures
(`h`, `rx`, `ry`, `rz`, `cx`, `state_init`, `state_set`), `bench_packed_size`,
`bench_tiled_size`, `bench_timing_barrier`, `bench_time_ns`, `bench_ns_to_ms`,
`BENCH_PQ_WARMUP_MS` from `bench.h`.

**Agent 4 — Step D (QuEST runner):** needs the QuEST v4 API header path
(`extern/QuEST/quest/include/quest.h`) and the list of function names to verify
(`applyHadamard`, `applyPauliX`, `applyPauliZ`, `applyRx`, `applyRy`, `applyRz`,
`applyControlledPauliX`, `createDensityQureg`, `destroyQureg`, `initMixedState` or
equivalent).

**Agent 5 — Step E (Qulacs runner):** needs the Qulacs csim header path
(`extern/qulacs/src/csim/update_ops_dm.hpp`) and the function signatures to verify
(`dm_H_gate`, `dm_X_gate`, `dm_Z_gate`, `dm_RX_gate`, `dm_RY_gate`, `dm_RZ_gate`,
`dm_CNOT_gate` with their argument orders).

**Agent 6 — Step F (Aer-DM runner):** needs the Aer-DM class path
(`vendor/aer-dm/src/simulators/density_matrix/...`) and the existing `bench_aer_dm.cpp`
patterns for `AER::QV::DensityMatrix<double>`, `AER::Utils::vectorize_matrix`,
`apply_unitary_matrix`, `AER::reg_t`.

### Wave 3: Step H (main orchestration)

**Assigned to:** single agent after all five Wave 2 files exist.

**What it needs:** the handoff contracts below — the exact function signatures exported by each
Wave 2 agent — so it can write `bench_circuit_main.c` without reading the runner source files.
Step I (CMakeLists.txt) is finalized in this wave after Step H is complete; it may be drafted
speculatively during Wave 2.

### Handoff Contracts

Every Wave 2 agent must export functions matching these exact signatures. The coordinator must
verify conformance before starting Wave 3.

**Step B (circuit constructors) — already declared in bench_circuit.h; no additional export.**

**Step C (qlib runner):**
```c
/* Exported from bench_circuit_qlib.c */
bench_result_t bench_circuit_qlib_packed(const circuit_t *c, int iterations,
                                          int warmup, int runs);
bench_result_t bench_circuit_qlib_tiled(const circuit_t *c, int iterations,
                                         int warmup, int runs);
```

**Step D (QuEST runner — guarded by WITH_QUEST):**
```c
/* Exported from bench_circuit_quest.c, inside #ifdef WITH_QUEST */
bench_result_t bench_circuit_quest(const circuit_t *c, int iterations,
                                    int warmup, int runs);
```

**Step E (Qulacs runner — guarded by WITH_QULACS, extern "C"):**
```cpp
/* Exported from bench_circuit_qulacs.cpp, inside #ifdef WITH_QULACS + extern "C" */
bench_result_t bench_circuit_qulacs(const circuit_t *c, int iterations,
                                     int warmup, int runs);
```

**Step F (Aer-DM runner — guarded by WITH_AER_DM, extern "C"):**
```cpp
/* Exported from bench_circuit_aer_dm.cpp, inside #ifdef WITH_AER_DM + extern "C" */
bench_result_t bench_circuit_aer_dm(const circuit_t *c, int iterations,
                                     int warmup, int runs);
```

**Context allocation (for calibrate_iterations compatibility):**

Each runner file must also export a pair of context alloc/free functions so that
`calibrate_iterations` in `bench_circuit_main.c` can call any runner through the `runner_fn`
function pointer:

```c
/* qlib (in bench_circuit_qlib.c) */
typedef struct { state_t *s; state_type_t type; } qlib_circuit_ctx_t;
qlib_circuit_ctx_t *qlib_circuit_ctx_alloc(int n_qubits, state_type_t type);
void                qlib_circuit_ctx_free(qlib_circuit_ctx_t *ctx);

/* QuEST (in bench_circuit_quest.c, #ifdef WITH_QUEST) */
typedef struct { Qureg rho; } quest_circuit_ctx_t;
quest_circuit_ctx_t *quest_circuit_ctx_alloc(int n_qubits);
void                 quest_circuit_ctx_free(quest_circuit_ctx_t *ctx);

/* Qulacs (in bench_circuit_qulacs.cpp, extern "C", #ifdef WITH_QULACS) */
typedef struct { CTYPE *rho; ITYPE dim; } qulacs_circuit_ctx_t;
qulacs_circuit_ctx_t *qulacs_circuit_ctx_alloc(int n_qubits);
void                  qulacs_circuit_ctx_free(qulacs_circuit_ctx_t *ctx);

/* Aer-DM (in bench_circuit_aer_dm.cpp, extern "C", #ifdef WITH_AER_DM) */
typedef struct { void *dm; int n_qubits; } aer_dm_circuit_ctx_t;
aer_dm_circuit_ctx_t *aer_dm_circuit_ctx_alloc(int n_qubits);
void                  aer_dm_circuit_ctx_free(aer_dm_circuit_ctx_t *ctx);
```

The `runner_fn` typedef in `bench_circuit_main.c` matches:
```c
typedef void (*runner_fn)(const circuit_t *c, void *ctx, int iterations);
```

Each runner wraps its internal loop in a static function with this exact signature. `ctx` is
cast to the appropriate `*_circuit_ctx_t` type inside the wrapper.

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

/*
 * Gate type identifiers used in circuit_op_t.
 * Only gates actually used by the three circuit constructors are defined.
 * Dispatch switches must handle all listed cases; add a default branch that
 * prints to stderr in debug builds for unimplemented gate types.
 */
typedef enum {
    GATE_H  = 0,
    GATE_X  = 1,
    GATE_Z  = 2,
    GATE_RX = 3,
    GATE_RY = 4,
    GATE_RZ = 5,
    GATE_CX = 6,
} circuit_gate_t;

/*
 * A single gate operation in a circuit.
 * One-qubit gates use q0; two-qubit gates use q0 (control/first) and q1 (target/second).
 * angle is used by rotation gates; ignored for non-parametric gates.
 */
typedef struct {
    circuit_gate_t gate_type;
    qubit_t        q0;
    qubit_t        q1;      /* used only by 2Q gates; 0 otherwise */
    double         angle;   /* used only by rotation gates; 0.0 otherwise */
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

/*
 * Append helpers.
 * All return 0 on success, -1 on allocation failure.
 * On failure, the circuit is left unchanged and an error is printed to stderr.
 */
int circuit_append_1q(circuit_t *c, circuit_gate_t g, qubit_t q0);
int circuit_append_1q_angle(circuit_t *c, circuit_gate_t g, qubit_t q0, double angle);
int circuit_append_2q(circuit_t *c, circuit_gate_t g, qubit_t q0, qubit_t q1);

/* RZZ composite: appends CX(q0,q1) · Rz(q1, gamma) · CX(q0,q1). Returns 0 or -1. */
int circuit_rzz(circuit_t *c, qubit_t q0, qubit_t q1, double gamma);

/* Circuit constructors. Return an empty circuit (n_ops == 0) on allocation failure. */
circuit_t circuit_qaoa_p1(int n, double gamma, double beta);
circuit_t circuit_vqe_hea(int n, int depth);     /* depth < 0 → use floor(n/2) */
circuit_t circuit_qv(int n, int seed);            /* n must be even */

#ifdef __cplusplus
}
#endif

#endif /* BENCH_CIRCUIT_H */
```

**Design decisions locked in:**
- `q2` field removed: no three-qubit gates are used; the extra 4 bytes per op wastes memory
  and creates false expectations. Remove it entirely.
- `qubit_t` used for `q0`/`q1` to match the rest of the codebase (not `int`).
- Enum contains only the 7 gates actually used by the three circuits. No dead entries.
- Append functions return `int` for explicit failure signaling.

#### bench_circuit_ops.c skeleton

Implement `circuit_alloc`, `circuit_free`, `circuit_append_1q`, `circuit_append_1q_angle`,
`circuit_append_2q`, and `circuit_rzz` here. Leave constructor bodies (`circuit_qaoa_p1`,
`circuit_vqe_hea`, `circuit_qv`) as stubs returning `circuit_alloc(n_qubits, 0, name)` — filled
in Step B.

Implementation notes:
- `circuit_alloc`: call `malloc(initial_cap * sizeof(circuit_op_t))`; check return; set
  `n_ops = 0`, `cap = initial_cap`. If `initial_cap == 0`, set `ops = NULL`.
- `circuit_append_*`: when `n_ops == cap`, `realloc` to `max(1, cap * 2)`. If realloc fails,
  print to stderr and return -1 without modifying the circuit.
- `circuit_rzz`: three calls in order: `circuit_append_2q(c, GATE_CX, q0, q1)`,
  `circuit_append_1q_angle(c, GATE_RZ, q1, gamma)`, `circuit_append_2q(c, GATE_CX, q0, q1)`.
  Return -1 if any append fails.
- `circuit_free`: free `c->ops`; zero out the struct.

---

### Step B — Implement circuit constructors

**Files:** `benchmark/src/bench_circuit_ops.c` (fill in the three constructor stubs)
**Depends on:** Step A
**Estimated time:** 2 hours

#### `circuit_qaoa_p1(int n, double gamma, double beta)`

```
initial_cap = 5*n
c = circuit_alloc(n, initial_cap, "qaoa_p1")
for q in 0..n-1:
    circuit_append_1q(&c, GATE_H, q)
for i in 0..n-1:
    q0 = i
    q1 = (i + 1) % n
    circuit_rzz(&c, q0, q1, gamma)   /* CX · Rz(gamma) · CX */
for q in 0..n-1:
    circuit_append_1q_angle(&c, GATE_RX, q, beta)
return c
```

Total ops: n + 3n + n = 5n. Add: `assert(c.n_ops == 5*n || c.n_ops == 0 /* alloc failed */)`.

#### `circuit_vqe_hea(int n, int depth)`

```
if depth < 0: depth = n / 2   /* floor division */
initial_cap = depth * (2*n - 1)
c = circuit_alloc(n, initial_cap, "vqe_hea")
for d in 0..depth-1:
    for q in 0..n-1:
        circuit_append_1q_angle(&c, GATE_RY, q, M_PI / 4.0)
    for q in 0..n-2:
        circuit_append_2q(&c, GATE_CX, q, q+1)
return c
```

Total ops: `depth * (n + (n-1)) = depth * (2n-1)`.

#### `circuit_qv(int n, int seed)`

Requirements:
- Assert `n % 2 == 0`; if odd, print to stderr and return `circuit_alloc(n, 0, "qv")`.
- Use a 32-bit LCG: `state = state * 1664525u + 1013904223u` (Numerical Recipes constants).
- Seed with the given `seed` value cast to `uint32_t`.
- n_pairs = n / 2 per layer; n_layers = n.
- Approximate total ops: n_layers × n_pairs × 15 = 15n²/2. Use this for initial_cap.

For each layer:
1. Generate a random permutation of qubits 0..n-1 using Fisher-Yates with the LCG
   (swap index j with a random index in [j, n-1]).
2. Pair adjacent elements: (perm[0], perm[1]), (perm[2], perm[3]), ...
3. For each pair (qa, qb):
   - 3 Ry on qa, angles from LCG.
   - 3 Ry on qb, angles from LCG.
   - CX(qa, qb), CX(qb, qa), CX(qa, qb).
   - 3 Rz on qa, angles from LCG.
   - 3 Rz on qb, angles from LCG.

LCG angle generation: `angle = (lcg_next(&state) / 4294967296.0) * 2.0 * M_PI`.

The resulting circuit is NOT a rigorous KAK decomposition; the interaction angles are random
draws, not computed from KAK theory. This is intentional: all competitors receive the same
gate sequence, so the benchmark measures throughput on a circuit of the same *structure* as
KAK-decomposed random unitaries. Add a prominent comment in the source file:

```c
/*
 * QV circuit template: random parameterized two-qubit circuit template.
 *
 * Each qubit pair receives: Ry(a0)⊗Ry(a1) · CX·CX†·CX · Rz(a2)⊗Rz(a3) where
 * a_i are random angles drawn from a uniform distribution over [0, 2π).
 * This exercises the same gate MIX as a KAK decomposition (9 single-qubit rotations
 * and 3 CNOT gates per pair) but the angles are not KAK-computed interaction parameters.
 * Do not describe this as a KAK decomposition in any external-facing text.
 */
```

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

`type` must be `MIXED_PACKED` or `MIXED_TILED`. On `PURE` or any unexpected value, fill
`result.method` with the appropriate name and return a zeroed result with an error to stderr.

#### Structure

```
Allocate state once with state_init(type, n_qubits) — outside timed loop.
Memory: bench_packed_size(n_qubits) or bench_tiled_size(n_qubits).

Build static dispatch:
  apply_op(state_t *s, const circuit_op_t *op) → switch(op->gate_type):
    GATE_H   → h(s, op->q0)
    GATE_X   → x(s, op->q0)
    GATE_Z   → z(s, op->q0)
    GATE_RX  → rx(s, op->q0, op->angle)
    GATE_RY  → ry(s, op->q0, op->angle)
    GATE_RZ  → rz(s, op->q0, op->angle)
    GATE_CX  → cx(s, op->q0, op->q1)
    default  → fprintf(stderr, "bench_circuit_qlib: unimplemented gate type %d\n", op->gate_type)

Maximally mixed state initialization helper:
  void init_maximally_mixed(state_t *s, qubit_t qubits, state_type_t type, dim_t dim):
    size_t nbytes = (type == MIXED_PACKED) ? bench_packed_size(qubits) : bench_tiled_size(qubits)
    memset(s->data, 0, nbytes)
    for i in 0..dim-1:
        state_set(s, i, i, 1.0 / (double)dim)
  (Do NOT use state_plus — that sets ALL elements to 1/dim, not only the diagonal.
   Note: state_byte_size does not exist in bench.h; use bench_packed_size / bench_tiled_size.)

Warmup (time-based, matching bench_run_timed semantics):
  init_maximally_mixed(&state, dim)
  t0 = bench_time_ns()
  do:
    for op in 0..c->n_ops-1: apply_op(&state, &c->ops[op])
  while bench_ns_to_ms(bench_time_ns() - t0) < BENCH_PQ_WARMUP_MS

Timed loop:
  for r in 0..runs-1:
    init_maximally_mixed(&state, dim)   /* re-initialize before each timed run */
    bench_timing_barrier()
    t0 = bench_time_ns()
    bench_timing_barrier()
    for i in 0..iterations-1:
      for op in 0..c->n_ops-1:
        apply_op(&state, &c->ops[op])
    bench_timing_barrier()
    run_times[r] = bench_ns_to_ms(bench_time_ns() - t0)
```

**State re-initialization policy (explicit decision):** State is re-initialized to maximally
mixed (ρ = I/dim) before each of the `runs` timed runs. Within a single run, state evolves
across `iterations` circuit applications — this is intentional. The benchmark measures
**sustained throughput** (the steady-state gate-application rate) rather than single-shot
latency from a fixed initial state. This is the correct metric for comparing simulator
performance in a VQO inner-loop context where the same circuit is applied many times.
Document this in `bench_circuit_main.c`'s file header comment.

Populate `bench_result_t`:
- `gate_name = c->name`
- `method = (type == MIXED_PACKED) ? "qlib_packed" : "qlib_tiled"`
- `iterations = iterations`
- `sweep_size = c->n_ops`  (total gate count is the natural sweep unit)
- `ops_per_sec = (iterations * c->n_ops) / (time_ms / 1000.0)` (from mean time)
- `time_per_gate_ms = time_ms / ((double)iterations * c->n_ops)` (from mean time)
- `memory_bytes` = packed or tiled size

---

### Step D — QuEST circuit runner

**File:** `benchmark/src/bench_circuit_quest.c`
**Depends on:** Step A; requires `WITH_QUEST` at compile time
**Estimated time:** 2.5 hours

#### Function signature (inside `#ifdef WITH_QUEST`)

```c
bench_result_t bench_circuit_quest(const circuit_t *c, int iterations, int warmup, int runs);
```

#### Structure

Pattern mirrors `bench_quest.c` exactly:
- Requires `bench_quest_init()` / `bench_quest_cleanup()` called by `bench_circuit_main.c`.
- Creates `Qureg rho = createDensityQureg(c->n_qubits)` once before the loop.
- Re-initializes to maximally mixed state before each timed run.

QuEST maximally mixed initialization: check `extern/QuEST/quest/include/quest.h` for the
correct function. Candidate names (QuEST v4): `initMixedState(rho)` or `setQuregToRenormMixed`.
If no direct function exists, zero all amplitudes then set diagonals to `1.0/dim` via
`setDensityAmps` or equivalent. Confirm before implementation — this is the primary API
uncertainty for Step D.

Gate dispatch:

```c
static void apply_circuit_op_quest(Qureg rho, const circuit_op_t *op) {
    switch (op->gate_type) {
        case GATE_H:   applyHadamard(rho, op->q0); break;
        case GATE_X:   applyPauliX(rho, op->q0); break;
        case GATE_Z:   applyPauliZ(rho, op->q0); break;
        case GATE_RX:  applyRx(rho, op->q0, op->angle); break;
        case GATE_RY:  applyRy(rho, op->q0, op->angle); break;
        case GATE_RZ:  applyRz(rho, op->q0, op->angle); break;
        case GATE_CX:  applyControlledPauliX(rho, op->q0, op->q1); break;
        default:
            fprintf(stderr, "bench_circuit_quest: unimplemented gate type %d\n", op->gate_type);
    }
}
```

Verify all function names against `extern/QuEST/quest/include/quest.h` before implementation.
QuEST v4 changed several API names from v3; `applyRx`/`applyRy`/`applyRz` in particular may
differ. This verification step is mandatory before writing any code for Step D.

Use `bench_timing_barrier()` around the timed region (matching the pattern in `bench_run_timed`):

```c
bench_timing_barrier();
t0 = bench_time_ns();
bench_timing_barrier();
for (int i = 0; i < iterations; ++i) { /* apply all ops */ }
bench_timing_barrier();
run_times[r] = bench_ns_to_ms(bench_time_ns() - t0);
```

Memory: `(size_t)dim * dim * sizeof(double) * 2` (matches `bench_quest.c` formula).

**Estimated time revised:** 2.5 hours (up from 2.0) to account for mandatory API verification.

---

### Step E — Qulacs circuit runner

**File:** `benchmark/src/bench_circuit_qulacs.cpp`
**Depends on:** Step A; requires `WITH_QULACS` at compile time
**Estimated time:** 2.5 hours

#### Function signature (inside `#ifdef WITH_QULACS`, `extern "C"` block)

```cpp
bench_result_t bench_circuit_qulacs(const circuit_t *c, int iterations, int warmup, int runs);
```

#### Structure

Pattern mirrors `bench_qulacs.cpp` exactly:
- Allocate `CTYPE *rho = (CTYPE *)malloc((size_t)dim * dim * sizeof(CTYPE))`.
- Re-initialize to maximally mixed state before each timed run: set diagonal elements to
  `CTYPE(1.0/dim, 0.0)` via `memset(rho, 0, ...)` followed by a loop over diagonals.

Gate dispatch:

```cpp
static void apply_circuit_op_qulacs(const circuit_op_t *op, CTYPE *rho, ITYPE dim) {
    switch (op->gate_type) {
        case GATE_H:   dm_H_gate(op->q0, rho, dim); break;
        case GATE_X:   dm_X_gate(op->q0, rho, dim); break;
        case GATE_Z:   dm_Z_gate(op->q0, rho, dim); break;
        case GATE_RX:  dm_RX_gate(op->q0, op->angle, rho, dim); break;
        case GATE_RY:  dm_RY_gate(op->q0, op->angle, rho, dim); break;
        case GATE_RZ:  dm_RZ_gate(op->q0, op->angle, rho, dim); break;
        case GATE_CX:  dm_CNOT_gate(op->q0, op->q1, rho, dim); break;
        default:
            fprintf(stderr, "bench_circuit_qulacs: unimplemented gate type %d\n", op->gate_type);
    }
}
```

Verify exact signatures against `extern/qulacs/src/csim/update_ops_dm.hpp` before
implementation. In particular:
- Confirm that `dm_RX_gate`, `dm_RY_gate`, `dm_RZ_gate` accept `(UINT target, double angle,
  CTYPE *state, ITYPE dim)` in that order.
- Confirm that `dm_CNOT_gate` accepts `(UINT control, UINT target, CTYPE *state, ITYPE dim)`.

This verification is mandatory before writing any dispatch code for this step.

Use `bench_timing_barrier()` around the timed region (same pattern as Step D).

Memory: `static_cast<size_t>(dim) * dim * sizeof(CTYPE)`.

**Estimated time revised:** 2.5 hours (up from 2.0) to account for mandatory API verification.

---

### Step F — Aer-DM circuit runner

**File:** `benchmark/src/bench_circuit_aer_dm.cpp`
**Depends on:** Step A; requires `WITH_AER_DM` at compile time
**Estimated time:** 2.5 hours

#### Function signature (inside `#ifdef WITH_AER_DM`, `extern "C"` block)

```cpp
bench_result_t bench_circuit_aer_dm(const circuit_t *c, int iterations, int warmup, int runs);
```

#### Structure

Pattern mirrors `bench_aer_dm.cpp` exactly:
- Use the `bench_aer_dm_alloc` / `bench_aer_dm_reinit` / `bench_aer_dm_free` API from `bench.h`
  if those functions support density matrix initialization to maximally mixed state. If not,
  replicate the pattern from `bench_aer_dm.cpp` directly.
- Re-initialize to maximally mixed state before each timed run.

#### Gate dispatch

The existing `bench_aer_dm.cpp` defines `gate_kind_t` as `{ BK_X, BK_H, BK_Z, BK_CX,
BK_SWAP, BK_UNKNOWN }`. It has **no Rx, Ry, or Rz cases** — only X, H, Z, CX, SWAP are
handled, all dispatched via `apply_unitary_matrix`. All three benchmark circuits (QAOA, VQE,
QV) use rotation gates heavily; they cannot be dispatched via the existing gate_kind mechanism.

For each rotation gate op in the circuit loop, build the 2×2 rotation matrix explicitly and
vectorize it with `AER::Utils::vectorize_matrix`, then dispatch via `apply_unitary_matrix`.
The matrices are:

- **Rx(θ)**: `[[cos(θ/2), -i·sin(θ/2)], [-i·sin(θ/2), cos(θ/2)]]`
- **Ry(θ)**: `[[cos(θ/2), -sin(θ/2)], [sin(θ/2), cos(θ/2)]]`
- **Rz(θ)**: `[[exp(-iθ/2), 0], [0, exp(iθ/2)]]`

This construction must happen per-op inside the circuit loop because QV angles vary per gate.
H, X, Z can reuse precomputed vectorized matrices (as `bench_aer_dm.cpp` does for H and Z).
CX must use the two-qubit `apply_unitary_matrix` overload matching `bench_aer_dm.cpp`.

**Per-op overhead note:** Building and vectorizing a 2×2 matrix per rotation gate adds
measurable per-op overhead for the Aer-DM runner that is absent from qlib's direct
`rx`/`ry`/`rz` calls. This overhead is inherent to the Aer-DM API and must be documented
as part of the Aer-DM runner's measurement in the thesis (it is not a benchmark artifact).

Wrap all dispatch in an `apply_circuit_op_aer_dm` helper:

```cpp
static void apply_circuit_op_aer_dm(AER::QV::DensityMatrix<double> &dm,
                                     const circuit_op_t *op) {
    const AER::reg_t q0_reg = {(uint_t)op->q0};
    switch (op->gate_type) {
        case GATE_H:  dm.apply_unitary_matrix(q0_reg, /* precomputed h_vec */); break;
        case GATE_X:  dm.apply_unitary_matrix(q0_reg, /* precomputed x_vec */); break;
        case GATE_Z:  dm.apply_unitary_matrix(q0_reg, /* precomputed z_vec */); break;
        case GATE_RX: {
            double c = cos(op->angle / 2.0), s = sin(op->angle / 2.0);
            AER::cmatrix_t m = {{c, {0,-s}}, {{0,-s}, c}};
            dm.apply_unitary_matrix(q0_reg, AER::Utils::vectorize_matrix(m));
            break;
        }
        case GATE_RY: {
            double c = cos(op->angle / 2.0), s = sin(op->angle / 2.0);
            AER::cmatrix_t m = {{c, -s}, {s, c}};
            dm.apply_unitary_matrix(q0_reg, AER::Utils::vectorize_matrix(m));
            break;
        }
        case GATE_RZ: {
            using namespace std::complex_literals;
            double h = op->angle / 2.0;
            AER::cmatrix_t m = {{std::exp(-1i*h), 0}, {0, std::exp(1i*h)}};
            dm.apply_unitary_matrix(q0_reg, AER::Utils::vectorize_matrix(m));
            break;
        }
        case GATE_CX: {
            const AER::reg_t cx_reg = {(uint_t)op->q0, (uint_t)op->q1};
            dm.apply_unitary_matrix(cx_reg, /* precomputed cx_vec */);
            break;
        }
        default:
            fprintf(stderr, "bench_circuit_aer_dm: unimplemented gate type %d\n",
                    op->gate_type);
    }
}
```

Precomputed vectorized matrices for H, X, Z, CX must be built once before the warmup loop
(matching the pattern of `bench_aer_dm.cpp`) and passed in or captured via lambda.

Add a `default` case to any other switch that prints to stderr.

Use `bench_timing_barrier()` around the timed region.

Memory: same formula as `bench_aer_dm.cpp` uses for `memory_bytes`:
`((size_t)1 << qubits) * ((size_t)1 << qubits) * sizeof(std::complex<double>)`.

---

### Step H — Write `benchmark/src/bench_circuit_main.c`

**File:** `benchmark/src/bench_circuit_main.c`
**Depends on:** Steps A, B, C, D, E (Qulacs), F (Aer-DM) — all runners and constructors complete
**Estimated time:** 3 hours

This step owns orchestration, output formatting, CLI interface, and correctness checks.

#### File header comment (mandatory)

```c
/*
 * bench_circuit_main.c — Circuit-level benchmark orchestration
 *
 * Measurement methodology:
 *   - All competitors are initialized to the maximally mixed state rho = I/dim before
 *     each of the `runs` independent timed runs.
 *   - Within a single timed run, the state evolves across `iterations` circuit applications.
 *     State is NOT re-initialized between iterations within a run.
 *   - This measures SUSTAINED THROUGHPUT — the steady-state gate-application rate —
 *     rather than single-shot latency from a fixed initial state.
 *   - Timing covers only the gate application loop. State allocation and re-initialization
 *     are outside the timed region.
 *   - All competitors receive the identical circuit_op_t sequence: same gate types,
 *     same qubit indices, same angles.
 *   - All competitors run in density matrix mode.
 *   - OpenMP: benchmarks run with OpenMP enabled (matching qlib's production configuration).
 *     Thread count is controlled by OMP_NUM_THREADS. For thesis results, record the
 *     thread count alongside benchmark output.
 */
```

#### Iteration calibration

Do NOT hardcode `--iterations 100`. Instead, implement calibration matching the existing
infrastructure.

Define the `runner_fn` callback type at file scope in `bench_circuit_main.c`:

```c
/* runner_fn: apply the circuit c for `iterations` repetitions using runner-specific
 * context `ctx`. Called by calibrate_iterations once per method. */
typedef void (*runner_fn)(const circuit_t *c, void *ctx, int iterations);
```

`ctx` is a runner-specific context struct allocated by the caller. For qlib, it contains
`state_t *s` and `state_type_t type`; for Aer-DM it contains the `DensityMatrix<double>*`
pointer; for QuEST it contains the `Qureg`; etc. `run_all_methods` is responsible for
allocating each runner's context struct before calling `calibrate_iterations`.

```c
static int calibrate_iterations(const circuit_t *c, runner_fn apply_circuit, void *ctx) {
    /* Run BENCH_PQ_PROBE_ITERS iterations, BENCH_PQ_PROBE_RUNS times, take minimum */
    /* Scale to target MIN_MEASUREMENT_MS total timed run duration */
    /* Clamp to [1, BENCH_PQ_MAX_CALIBRATED_ITERATIONS] */
    /* Fall back to BENCH_PQ_FALLBACK_ITERATIONS if probe rounds to 0 ns */
}
```

The calibration probes each (circuit, method) combination independently. When `--iterations N`
is specified explicitly on the CLI (detect via `iterations_explicit` flag matching
`bench_options_t`), skip calibration and use N directly.

Default for `--iterations`: use calibration (not a fixed constant).

#### Extended options struct

Define a `circuit_options_t` in this file:

```c
typedef struct {
    bench_options_t base;    /* reuses min_qubits, max_qubits, step, iterations,
                                warmup, runs, csv_output, pgfplots_output, verbose,
                                iterations_explicit, runs_explicit */
    int    run_qaoa;         /* 1 if "all" or "qaoa" selected */
    int    run_vqe;          /* 1 if "all" or "vqe" selected */
    int    run_qv;           /* 1 if "all" or "qv" selected */
    int    seed;             /* --seed N, default 42 */
    double gamma;            /* --gamma FLOAT, default M_PI/4 */
    double beta;             /* --beta FLOAT, default M_PI/4 */
} circuit_options_t;
```

#### CLI flags

| Flag | Default | Description |
|------|---------|-------------|
| `--min-qubits N` | 2 | Minimum qubit count |
| `--max-qubits N` | 12 | Maximum qubit count |
| `--step N` | 1 | Qubit count increment |
| `--iterations N` | calibrated | Circuit applications per timed window; disables calibration when specified |
| `--warmup` | — | Warmup runs until `BENCH_PQ_WARMUP_MS` ms elapsed (time-based, not count-based) |
| `--runs N` | 15 | Independent timing runs (default matches `BENCH_PQ_DEFAULT_RUNS`) |
| `--csv` | off | CSV output to stdout |
| `--pgfplots` | off | CSV to stdout in pgfplots-compatible format (see Output Formats) |
| `--verbose` | off | Extended per-result info (min, median, memory) |
| `--circuit NAME` | all | Select circuit: `qaoa`, `vqe`, `qv`, or `all` |
| `--seed N` | 42 | RNG seed for QV circuit |
| `--gamma FLOAT` | π/4 | QAOA phase separation parameter |
| `--beta FLOAT` | π/4 | QAOA mixing parameter |
| `--help` | — | Print usage and exit |

Validation:
- `--circuit` must be one of `{qaoa, vqe, qv, all}`; print usage and exit(1) otherwise.
- `--seed` must be ≥ 0; exit(1) otherwise.
- `--gamma`, `--beta`: any finite double; exit(1) on non-finite value.
- QV requires even qubit counts; skip odd qubit counts silently with a per-qubit warning.
- `--runs` clamped to `[1, BENCH_MAX_RUNS]`.

#### Correctness verification (debug builds only)

For n ≤ 4, after collecting all timed results, run a verification pass:
- Apply the circuit once from the same initial state (ρ = I/dim) using each available method.
- Compare all output density matrices pairwise using Frobenius norm: `||ρ_a - ρ_b||_F < 1e-6`.
- Print `[PASS]` or `[FAIL: method_a vs method_b: err = X.XXe-YY]` to stderr.
- This check runs ONLY in debug builds (`#ifndef NDEBUG`) and is NOT part of the timed region.

**Aer-DM correctness init:** `bench_aer_dm_reinit` fills the density matrix with random bytes
(it is designed to defeat fast-path optimizations for special states, not to set ρ = I/dim).
For the correctness verification pass, the Aer-DM backend must use a dedicated initialization:
call `dm->initialize()` (which sets ρ = |0⟩⟨0|) and then set each diagonal element to
`1.0/dim` while zeroing off-diagonals — or use whichever `AER::QV::DensityMatrix` method sets
a mixed state directly. This dedicated init is used **only** in the correctness pass; the
timed runs continue to use `bench_aer_dm_reinit` (random bytes). Document this distinction
with a comment at the correctness-check call site.

#### Output formats

**CSV output** (`--csv`): emit one row per (circuit, qubits, method) combination:

```
circuit,qubits,dim,method,runs,iterations,n_ops,time_ms_mean,time_ms_std,time_ms_min,
time_ms_median,time_ms_cv,time_per_gate_ms,ops_per_sec,memory_bytes
```

The `circuit` column (circuit name) is the first column. `n_ops` replaces `sweep_size` in the
header to be self-documenting. The `bench_result_t.sweep_size` field holds `c->n_ops`.

**Console output**: for each `(circuit, qubits)` combination, print a block:

```
Circuit: qaoa_p1  qubits: 6  n_ops: 30  [OpenMP threads: N]
  qlib_packed    12.345 ± 0.123 ms (CV  1.2%)   1234567 ops/sec
  qlib_tiled      9.876 ± 0.089 ms (CV  0.9%)   1543210 ops/sec
  quest          34.567 ± 0.456 ms (CV  1.3%)    123456 ops/sec
  qulacs         18.234 ± 0.234 ms (CV  1.2%)    234567 ops/sec
  qiskit_aer_dm  21.000 ± 0.300 ms (CV  1.4%)    200000 ops/sec
  Speedup: 1.25x tiled vs packed
─────────────────────────────────────────────────────────────
```

The OpenMP thread count line is printed once per binary invocation, not once per block.

**pgfplots output** (`--pgfplots`): emit CSV to stdout in a format suitable for direct
`\pgfplotstableread` consumption. Do NOT embed LaTeX formatting. Produce one CSV file per
circuit (written to stdout with a `# CIRCUIT: name` comment header to allow splitting by a
post-processing script). Columns: `qubits method time_us_per_gate_mean time_us_per_gate_min
time_us_per_gate_median time_us_per_gate_cv_pct memory_bytes`.
Use microseconds per gate (`time_ms × 1000 / n_ops`) as the primary unit — consistent with
`bench_gate`'s per-gate presentation and suitable for pgfplots axis labels.

Generating actual `.dat` files or LaTeX table code is the responsibility of a post-processing
script, not this binary.

#### main() outline

```c
int main(int argc, char *argv[]) {
    circuit_options_t opts = circuit_parse_options(argc, argv);
    /* Print OpenMP thread count once */
#ifdef WITH_QUEST
    bench_quest_init();
#endif

    for (int qubits = opts.base.min_qubits;
         qubits <= opts.base.max_qubits;
         qubits += opts.base.step) {

        if (opts.run_qaoa) {
            circuit_t c = circuit_qaoa_p1(qubits, opts.gamma, opts.beta);
            run_all_methods(&c, &opts);
            circuit_free(&c);
        }
        if (opts.run_vqe) {
            circuit_t c = circuit_vqe_hea(qubits, -1);
            run_all_methods(&c, &opts);
            circuit_free(&c);
        }
        if (opts.run_qv && qubits % 2 == 0) {
            circuit_t c = circuit_qv(qubits, opts.seed);
            run_all_methods(&c, &opts);
            circuit_free(&c);
        } else if (opts.run_qv) {
            fprintf(stderr, "warn: QV skipped for odd qubit count %d\n", qubits);
        }
    }

#ifdef WITH_QUEST
    bench_quest_cleanup();
#endif
    return 0;
}
```

`run_all_methods` calls each runner, collects results into a local array, then calls the
appropriate output function.

---

### Step I — Add `bench_circuit` target to `benchmark/CMakeLists.txt`

**File:** `benchmark/CMakeLists.txt`
**Depends on:** Step H (all source files must exist)
**Estimated time:** 1 hour

Append after the existing `bench_gate` target block:

```cmake
#=======================================================================================================================
# bench_circuit - Circuit-level benchmark
#=======================================================================================================================

set(BENCH_CIRCUIT_SRC
        "${BENCH_SRC_DIR}/bench_circuit_ops.c"
        "${BENCH_SRC_DIR}/bench_circuit_qlib.c"
        "${BENCH_SRC_DIR}/bench_circuit_main.c"
)

if(QuEST_FOUND)
    list(APPEND BENCH_CIRCUIT_SRC "${BENCH_SRC_DIR}/bench_circuit_quest.c")
endif()

if(QULACS_FOUND)
    list(APPEND BENCH_CIRCUIT_SRC "${BENCH_SRC_DIR}/bench_circuit_qulacs.cpp")
endif()

if(AER_DM_FOUND)
    list(APPEND BENCH_CIRCUIT_SRC "${BENCH_SRC_DIR}/bench_circuit_aer_dm.cpp")
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

# q target provides qlib gates, state, and transitive BLAS/OpenMP dependencies.
# Verify BLAS linkage scope in the root CMakeLists.txt: if BLAS is PRIVATE to q,
# add an explicit BLAS link here (see note in Step I specification).
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

if(AER_DM_FOUND)
    target_compile_definitions(bench_circuit PRIVATE WITH_AER_DM)
    target_link_libraries(bench_circuit PRIVATE aer_dm_static)
    set_property(TARGET bench_circuit PROPERTY CXX_STANDARD 17)
    set_property(TARGET bench_circuit PROPERTY CXX_STANDARD_REQUIRED ON)
endif()

message(STATUS "  bench_circuit: enabled (circuit-level benchmarks)")
```

**`build_all_pairs` resolution (definitive decision):** `bench_circuit`'s runners do not
need `build_all_pairs`. Circuit operations already carry explicit qubit indices in `q0`/`q1`.
Do NOT create `bench_utils.c`. Do NOT add `bench_mixed.c` to `BENCH_CIRCUIT_SRC` — it would
pull in unnecessary gate-benchmark infrastructure (sweep loops, per-qubit pair tables) that
conflicts with `bench_circuit`'s design and introduces symbol conflicts with functions already
defined or redefined in the circuit runner files. (`main` is in `bench_main.c`, not
`bench_mixed.c`.) The symbol remains defined in `bench_mixed.c` and linked only by
`bench_gate`. No action required for `bench_circuit`.

---

### Step J — Update `docs/BENCHMARK.md`

**File:** `docs/BENCHMARK.md`
**Depends on:** Steps H, I (implementation complete; CLI and output format finalized)
**Estimated time:** 1 hour

Changes required:

1. **Intro paragraph**: replace "Circuit-level benchmarking is out of scope" with a sentence
   describing `bench_circuit` and pointing to its subsection.

2. **File Layout section**: add all new `bench_circuit_*` files to the tree.

3. **CMake targets table**: add `bench_circuit | Executable (EXCLUDE_FROM_ALL) | Always (Release only)`.

4. **CLI Reference**: add `### bench_circuit` subsection with the full flag table from Step H,
   default values, and validation behavior. Note that `--iterations` uses calibration by
   default (unlike `bench_gate`'s fixed 1000).

5. **New section — Circuits Benchmarked**: describe the three circuits, gate counts, parameter
   conventions, and recommended qubit ranges:
   - QAOA: 4–16 qubits (5n gates; extend to 20 if runtimes permit).
   - VQE HEA: 4–16 qubits (d(2n-1) gates; d=⌊n/2⌋).
   - QV: 4–14 qubits, even only (~15n²/2 gates).
   - Note: QV circuit uses "random parameterized two-qubit circuit template" — not KAK.

6. **Output Formats**: document the `circuit` column, `n_ops` field, pgfplots CSV format, and
   the requirement to use a post-processing script for `.dat` file generation.

7. **Measurement methodology note**: document the sustained-throughput vs. single-shot
   distinction.

8. **OpenMP note**: state that benchmarks run with OpenMP enabled; thread count should be
   recorded alongside results.

---

### Step K — Update thesis notes

**File:** `~/Projects/thesis/2.Simulation/notes/BENCHMARK_THESIS_NOTES.md`
**Depends on:** Steps H, I, J (full implementation and spec finalized)
**Estimated time:** 1.5 hours

Add a new section `## Circuit-Level Benchmarks (bench_circuit)` covering:

1. **Circuit selection and scientific motivation**
   - QAOA p=1: primary VQA thesis result; demonstrates end-to-end algorithm performance.
   - VQE HEA: continuous-parameter variational circuit; tests sustained rotation-gate throughput.
   - QV: generalization beyond VQA; tests random deep circuits with the gate mix of random
     two-qubit unitaries; QV metric widely used in industry benchmarking.

2. **QV circuit description for thesis**

   State explicitly: "The QV circuit uses a random parameterized two-qubit circuit template
   in which each qubit pair receives three Ry rotations per qubit, three CX gates
   (CX(a,b)·CX(b,a)·CX(a,b)), and three Rz rotations per qubit, with angles drawn uniformly
   from [0, 2π). This exercises the same gate mix as a KAK decomposition but is not a
   rigorous KAK decomposition — the angles are random draws, not KAK-computed interaction
   parameters. All competitors receive the identical gate sequence."

3. **Parameter reproducibility**
   - QAOA: γ = β = π/4 (not optimized; chosen so all rotation gates are non-trivial).
   - QV: seed = 42. LCG formula: `state ← state × 1664525 + 1013904223 (mod 2³²)`.
     State these constants explicitly in the thesis for full reproducibility.

4. **Cross-competitor fairness**
   - Identical `circuit_op_t` sequence for all competitors.
   - All initialized to ρ = I/dim before each timed run.
   - All use density matrix mode.

5. **Sustained throughput vs. single-shot latency**

   State that the benchmark measures sustained throughput (state evolves across `iterations`
   circuit applications per timed run). This matches the workload of a VQO inner loop where
   the same circuit is applied many times with the same parameters. Single-shot latency (first
   application from a fixed initial state) is a distinct metric not measured here.

6. **Cache behavior note**

   At n=12, the state density matrix is 256 MB (packed: less; full: 256² × 16 B ≈ 256 MB),
   exceeding typical L3 cache. State re-initialization between timed runs flushes cache,
   making the first gates of each run cold. QAOA's sequential qubit pair pattern is
   cache-friendly; QV's random permutations are worst-case. Document this in the thesis
   methodology section when presenting large-n results.

7. **OpenMP thread configuration**: record the value of `OMP_NUM_THREADS` and the hardware
   thread count alongside all benchmark results.

8. **Recommended qubit ranges for thesis plots**

   | Circuit | qlib/QuEST/Qulacs/Aer | Notes |
   |---------|----------------------|-------|
   | QAOA    | 4–16                 | Extend to 20 if time permits |
   | VQE HEA | 4–16                | |
   | QV      | 4–14 (even)         | ~15n²/2 gates; can be slow at n=14 |

---

## Dependency Summary Table

| Step | Description | Depends on | Parallelizable with |
|------|-------------|------------|---------------------|
| A | Header + ops skeleton | — | (start point) |
| B | Circuit constructors | A | C, D, E, F |
| C | qlib runner | A | B, D, E, F |
| D | QuEST runner | A | B, C, E, F |
| E | Qulacs runner | A | B, C, D, F |
| F | Aer-DM runner | A | B, C, D, E |
| H | main + CLI + output | A, B, C, D, E, F | — |
| I | CMakeLists.txt update | H | — |
| J | docs/BENCHMARK.md update | H, I | — |
| K | Thesis notes update | H, I, J | — |

Maximum parallelism: 5 agents can work simultaneously on Steps B–F after Step A completes.
Critical path: A → max(B, C, D, E, F) → H → I → J → K.

---

## Time Estimates

| Step | Time |
|------|------|
| A | 1.5 h |
| B | 2.0 h |
| C | 2.0 h |
| D (QuEST) | 2.5 h |
| E (Qulacs) | 2.5 h |
| F (Aer-DM) | 2.5 h |
| H | 3.0 h |
| I | 1.0 h |
| J | 1.0 h |
| K | 1.5 h |
| **Total (serial)** | **19.5 h** |
| **Critical path** | **A + max(B..F) + H + I + J + K = 1.5 + 2.5 + 3.0 + 1.0 + 1.0 + 1.5 = 10.5 h** |

With 5 agents in parallel on Steps B–F, wall-clock time is approximately 10.5 hours.

---

## Resolved Issues (reference for reviewers)

The following issues raised in review were resolved; decisions are final.

**State re-initialization:** re-initialized to ρ = I/dim before each timed run; NOT between
iterations within a run. This measures sustained throughput, which is the correct metric for
VQO inner-loop simulation. Documented in the file header comment of `bench_circuit_main.c`.

**QV naming:** renamed from "KAK decomposition" to "random parameterized two-qubit circuit
template." All source comments, documentation, and thesis notes must use this terminology.

**Aer-DM runner:** added as Step F with `bench_circuit_aer_dm.cpp` guarded by
`WITH_AER_DM`, matching the pattern of the existing `bench_aer_dm.cpp`.

**`build_all_pairs` ownership:** no action required for `bench_circuit`. Circuit ops carry
explicit qubit indices; `build_all_pairs` is not needed. No `bench_utils.c` will be created.

**`q2` field:** removed from `circuit_op_t`. No three-qubit gates are used.

**`int` vs `qubit_t`:** `q0`/`q1` use `qubit_t` throughout, matching the codebase.

**Dead enum entries:** removed. Only the 7 gate types used by the three circuits are defined.

**`default` branch in dispatch switches:** all switch statements in all runners print to
stderr rather than silently discarding unknown gate types.

**Warmup strategy:** time-based warmup (run until `BENCH_PQ_WARMUP_MS` elapsed), matching
`bench_run_timed`. Count-based warmup is not used.

**Iteration calibration:** implemented in `bench_circuit_main.c` using the same calibration
constants (`BENCH_PQ_PROBE_ITERS`, `BENCH_PQ_PROBE_RUNS`, `MIN_MEASUREMENT_MS`) as the
existing per-qubit infrastructure. Fixed `--iterations 100` default is rejected.

**`bench_timing_barrier()`:** all runners bracket the timed region with
`bench_timing_barrier()` calls, matching the `bench_run_timed` pattern in `bench.h`.

**pgfplots output:** emits CSV to stdout with a comment header for script-based splitting.
No LaTeX formatting baked into the binary. Post-processing produces `.dat` files.

**OpenMP:** benchmarks run with OpenMP enabled. Thread count printed once per invocation in
console output; must be recorded with thesis results.

**Correctness verification:** debug-only (guarded by `#ifndef NDEBUG`), runs for n ≤ 4,
compares all competitors pairwise using Frobenius norm tolerance 1e-6, outside the timed
region.

**Memory calculation for QV n=8:** n=8 has 8 layers × 4 pairs × 15 ops/pair = 480 ops
(not 40 as previously stated).

**Allocation failures:** `circuit_append_*` returns -1 on failure; constructors return an
empty circuit; runners check for NULL state/rho and return zeroed results with stderr messages.

**blas_dense removal:** the blas_dense runner was cut from this plan. For per-gate benchmarks
it provided a useful O(dim³) floor comparison; for full-circuit benchmarks it adds no insight
and its O(dim³) per-gate cost makes it impractical at any qubit count where the other
competitors are interesting. The competitor set is: qlib_packed, qlib_tiled, QuEST, Qulacs,
Aer-DM.
