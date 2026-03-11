# Qiskit-Aer Density Matrix Adapter: Implementation Plan

**Date:** 2026-03-11
**Author:** Systems integration plan for adding `qiskit_aer_dm` as a 7th competitor to `bench_gate`.
**Approach confirmed analogous to Qulacs.** Differences noted where they exist.

---

## Overview

The plan adds Qiskit-Aer's internal C++ density matrix kernel (`AER::QV::DensityMatrix<double>`)
as a 7th benchmark competitor (enum index 6, `M_AER_DM`), following the same vendored-header /
static-adapter pattern used for Qulacs. Because existing benchmark data for other methods already
took up to 24 hours to collect, an additional `bench_aer_only` executable is included to collect
only the missing Qiskit-Aer column without re-running any other method.

**`bench_aer_only` implementation strategy:** `bench_aer_only` is **not** a standalone
reimplementation. It is compiled from the *identical* source list as `bench_gate` (same
`bench_main.c`, same `bench.h`, same output functions, same pgfplots logic) with one additional
compile definition `-DBENCH_AER_ONLY_MODE=1`. `bench_main.c` uses this symbol to skip every
method call except `M_AER_DM`. This gives `bench_aer_only` bit-for-bit identical CLI flags,
defaults, output modes (console, CSV, pgfplots), statistical options, verbose mode, and
`pgfplots_data_t` structure — with zero code duplication.

**Naming note:** `bench_gate` has 7 method slots in the enum (`M_QLIB`, `M_QLIB_TILED`,
`M_BLAS`, `M_NAIVE`, `M_QUEST`, `M_QULACS`, `M_AER_DM`). The pgfplots output function
`print_pgfplots_output()` emits 6 columns because `M_NAIVE` is excluded from pgfplots tables
(naive runs at reduced iterations and is only useful in console/CSV mode). This is consistent
with the current code where `method_names[]` has 5 entries omitting naive; adding Aer-DM brings
the pgfplots column count to 6.

---

## Dependency Graph

```
Step 1 (vendor headers)
    │
    ├──────────────────────────────────────────────────┐
    ▼                                                  ▼
Step 2 (write bench_aer_dm.cpp)             Step 3a (CMake detection + INTERFACE target)
    │                                                  │
    └──────────────────┬───────────────────────────────┘
                       │  (both 2 and 3a must be done)
                       ▼
          Step 3b (CMake source list + target wiring)
                       │
                       ▼
          Step 4 (update bench.h)
                       │
                ┌──────┴──────┐
                ▼             ▼
            Step 5          Step 8
        (bench_main.c)  (docs/BENCHMARK.md)
                │         (can run in parallel
          ┌─────┴─────┐    with Step 5)
          ▼           ▼
      Step 6       Step 7
  (thesis notes) (bench_aer_only CMake only)
  (depends on     (CMake additions are independent
   Step 5;         of Step 5 source edits, but the
   does not        binary is not correct until Step 5
   block 7/8)      guards are in place)
```

**Parallel batches:**

| Batch | Steps | Prerequisite |
|-------|-------|--------------|
| A | Step 1 | None |
| B | Steps 2, 3a | Step 1 |
| C | Step 3b | Steps 2 and 3a |
| D | Step 4 | Step 3b |
| E | Steps 5, 8 | Step 4 |
| F | Steps 6, 7 | Step 5 |

Step 3 is split into two sub-steps: **3a** (CMake detection block and INTERFACE target) is
independent of Step 2 and can run in parallel with it. **3b** (appending `bench_aer_dm.cpp` to
`BENCH_MIXED_SRC` and wiring the `bench_gate` target) depends on Step 2, because
`bench_aer_dm.cpp` must exist before CMake references it in the source list. Steps 6 and 7 both
depend on Step 5 and can run concurrently with each other. Step 6 depends on Step 5 because the
thesis notes describe the adapter design (method names, API surface) which is determined in Step
2 and validated in Step 5. Step 7 depends on Step 5 because the `#ifndef BENCH_AER_ONLY_MODE`
guards added in Step 5 must exist before `bench_aer_only` produces correct output (without them,
the binary compiles but runs all methods, defeating its purpose). Step 8 has no dependency on
Steps 5, 6, or 7 and can run concurrently with Step 5.

---

## Step 1 — Vendor Qiskit-Aer Headers

**Estimated effort:** 3–5 hours
**Depends on:** nothing
**Blocks:** Steps 2, 3a

### What to do

Clone the Qiskit-Aer repository (tag `0.14.x` or latest stable) and copy the required headers
into `extern/aer-dm/`. The goal is to vendor exactly the headers needed by
`AER::QV::DensityMatrix<double>` and nothing else.

**Target directory structure (indicative, not exhaustive):**

The transitive include closure through `DensityMatrix`'s base class hierarchy, SIMD dispatch
layers, and AER matrix types is likely 15–30 headers in practice. The listing below is a
starting point; the compiler will identify every missing header during the probe step. Add only
what is required.

```
extern/aer-dm/
└── include/
    ├── simulators/
    │   ├── density_matrix/
    │   │   └── densitymatrix.hpp          # Primary: AER::QV::DensityMatrix<double>
    │   └── statevector/
    │       └── qubitvector.hpp            # Base class: DensityMatrix inherits from QubitVector
    └── framework/
        ├── types.hpp                      # AER complex types, AERcomplex typedef
        ├── matrix.hpp                     # AER::matrix<T> container
        ├── utils.hpp                      # Utility functions used by DM kernel
        ├── avx2_types.hpp                 # SIMD type aliases (include if present)
        └── json.hpp                       # nlohmann/json — required hard dependency (see note below)
```

**`nlohmann/json.hpp` is a required hard dependency, not optional.** `framework/types.hpp`
unconditionally includes `nlohmann/json.hpp` (~25,000 lines of header-only template code).
Always copy it; omitting it will cause a compile error. Be aware that this header significantly
increases compilation time for `bench_aer_dm.cpp` (expect 10–30 seconds per TU on typical
hardware). If build times become intolerable, investigate whether a stub `json.hpp` that defines
only the types actually used by `types.hpp` is feasible — but this requires careful inspection
and is not part of the standard vendoring path.

Note: `DensityMatrix<double>` inherits from `QubitVector` (in
`simulators/statevector/qubitvector.hpp`), not from `UnitaryMatrix`. Do not copy
`simulators/unitary/unitarymatrix.hpp` unless the compiler explicitly demands it.

**Shell commands to perform the copy:**

```bash
# Assumes Qiskit-Aer source is checked out to /tmp/qiskit-aer
QISKIT_REPO=/tmp/qiskit-aer
DEST=/path/to/qlib/extern/aer-dm/include

mkdir -p ${DEST}/simulators/density_matrix
mkdir -p ${DEST}/simulators/statevector
mkdir -p ${DEST}/framework

cp ${QISKIT_REPO}/src/simulators/density_matrix/densitymatrix.hpp \
   ${DEST}/simulators/density_matrix/

cp ${QISKIT_REPO}/src/simulators/statevector/qubitvector.hpp \
   ${DEST}/simulators/statevector/

cp ${QISKIT_REPO}/src/framework/types.hpp \
   ${QISKIT_REPO}/src/framework/matrix.hpp \
   ${QISKIT_REPO}/src/framework/utils.hpp \
   ${DEST}/framework/
# Always copy json.hpp — it is unconditionally included by types.hpp:
cp ${QISKIT_REPO}/src/framework/json.hpp ${DEST}/framework/
# Add avx2_types.hpp if the compiler requires it
```

**Version pinning — record the exact commit used:**

After copying headers, record the Qiskit-Aer git commit hash so that the vendor snapshot is
reproducible:

```bash
git -C /tmp/qiskit-aer rev-parse HEAD > /path/to/qlib/extern/aer-dm/VERSION
# Also add a comment in the CMake detection block: "vendored from <hash>"
```

The `extern/aer-dm/VERSION` file must be committed alongside the headers.

**Verification — compile a minimal probe program:**

The probe serves three purposes: (1) confirm which headers are transitively required, (2)
verify the actual API surface of `AER::QV::DensityMatrix<double>` before Step 2 is written,
and (3) verify the object model (heap vs. stack allocation of internal buffers).
**Step 2 must not be written until this probe compiles successfully and the API is confirmed.**

The probe tests six things: constructor convention, data accessor, sizeof, and all five gate
methods used in the adapter.

```bash
cat > /tmp/aer_probe.cpp << 'EOF'
#include "simulators/density_matrix/densitymatrix.hpp"
#include <cassert>
#include <complex>
#include <cstdio>

int main() {
    // (A) Constructor convention: confirm whether the argument is num_qubits or num_states.
    //     If the constructor takes num_qubits, dm(4) produces a 16x16 density matrix.
    //     If it takes num_states (= 2^n), dm(4) produces a 4-qubit state only if num_states=16.
    //     Check the resulting state dimension to disambiguate.
    AER::QV::DensityMatrix<double> dm(4);

    // (B) sizeof check: confirm the DensityMatrix object itself is small (internal buffer
    //     is almost certainly heap-allocated by the constructor via std::vector or similar).
    //     If this prints a value larger than a few KB, the object model is unusual.
    //     The object itself should be small (pointer + metadata); all state is on the heap.
    printf("sizeof(DensityMatrix<double>) = %zu\n", sizeof(dm));
    // Expected: a small value (< 1 KB). The 4-qubit density matrix buffer is
    //   (2^4)^2 * 16 bytes = 4096 bytes and must be on the heap, not inside the object.

    // (C) Data accessor: verify dm.data() returns std::complex<double>*
    //     or find the actual accessor exposed by the QubitVector base class.
    auto *p = dm.data();
    p[0] = {1.0, 0.0};   // write a value; confirms raw-pointer access

    // (D) Gate methods: test all five gates used in bench_aer_dm.cpp.
    //     These method names are NOT guaranteed to exist — the probe confirms them.
    //     If any call fails to compile, comment it out and note the correct alternative
    //     (e.g., apply_unitary_matrix with an explicit gate matrix).
    dm.apply_x(0);         // X gate on qubit 0
    dm.apply_h(0);         // H gate on qubit 0 (may be apply_hadamard)
    dm.apply_z(0);         // Z gate on qubit 0
    dm.apply_cnot(0, 1);   // CX gate, control=0 target=1
    dm.apply_swap(0, 1);   // SWAP gate on qubits 0 and 1

    return 0;
}
EOF
g++ -std=c++17 -I extern/aer-dm/include -O2 /tmp/aer_probe.cpp -o /dev/null
```

Iterate until the probe compiles without errors. Each missing include will name the required
header explicitly. Add only what the compiler demands.

**Resolving API mismatches:** If any gate method does not exist on
`AER::QV::DensityMatrix<double>`, use `apply_unitary_matrix(qubit_set, mat)` where `mat` is the
2×2 (or 4×4 for 2Q) gate matrix as an `AER::matrix<std::complex<double>>`. Document the actual
method names found and carry them into Step 2.

**Resolving data accessor mismatches:** If `dm.data()` does not exist or does not return a raw
pointer, inspect the `QubitVector` base class for an alternative (e.g., `data_pointer()`,
`begin()`, or an array subscript operator). If no raw pointer is available, initialize the state
by applying a sequence of random unitary gates to force non-trivial values before the timed
region.

**Pybind11 check — confirm no Python dependency leaks in:**

```bash
grep -r "pybind11" extern/aer-dm/include/
# Must produce no output.
```

If `pybind11` appears in any vendored header, that header pulls in the Python wrapper layer and
must NOT be included. The density matrix simulation core (`densitymatrix.hpp`,
`qubitvector.hpp`, `framework/`) should be entirely free of pybind11 references.

**Do NOT copy:**
- `src/controllers/` — Python controller layer
- `src/transpile/` — circuit transformation layer
- Any file matching `*pybind*.hpp` or `*python*.hpp`

---

## Step 2 — Write `benchmark/src/bench_aer_dm.cpp`

**Estimated effort:** 2–3 hours
**Depends on:** Step 1 (probe must have compiled and API surface confirmed)
**Blocks:** Step 3b, Step 4

> **WARNING: ALL CODE SNIPPETS IN THIS STEP ARE TEMPLATES AND PLACEHOLDERS.**
>
> The gate method names (`apply_x`, `apply_h`, `apply_z`, `apply_cnot`, `apply_swap`), the
> constructor convention (num_qubits vs. num_states), and the data accessor (`dm.data()`) shown
> below are **assumed names** that have NOT been verified. They MUST be replaced with the
> names confirmed by the Step 1 probe before this file is written. Do not write
> `bench_aer_dm.cpp` until the probe compiles successfully and the actual API is known.
> Every occurrence of a method name in the snippets below is a placeholder.

### What to do

Create `benchmark/src/bench_aer_dm.cpp` modelled directly on `benchmark/src/bench_qulacs.cpp`.

**Key differences from Qulacs:**

| Aspect | Qulacs | Aer DM |
|--------|--------|--------|
| Include | `<csim/update_ops_dm.hpp>` | `"simulators/density_matrix/densitymatrix.hpp"` |
| State type | `CTYPE*` (raw pointer, manually heap-allocated) | `AER::QV::DensityMatrix<double>` (local variable; internal buffer is heap-allocated by constructor) |
| Allocation | `std::malloc` / `std::free` | Constructor / destructor |
| 1Q gate API | `dm_X_gate(qubit, rho, dim)` | confirmed in Step 1 probe (placeholder: `dm.apply_x(qubit)`) |
| 1Q H gate | `dm_H_gate(qubit, rho, dim)` | confirmed in Step 1 probe (placeholder: `dm.apply_h(qubit)`) |
| 1Q Z gate | `dm_Z_gate(qubit, rho, dim)` | confirmed in Step 1 probe (placeholder: `dm.apply_z(qubit)`) |
| 2Q CX gate | `dm_CNOT_gate(q1, q2, rho, dim)` | confirmed in Step 1 probe (placeholder: `dm.apply_cnot(ctrl, tgt)`) |
| 2Q SWAP gate | `dm_SWAP_gate(q1, q2, rho, dim)` | confirmed in Step 1 probe (placeholder: `dm.apply_swap(q0, q1)`) |
| Global state | None | None — constructor/destructor only |
| Random init | Manual loop over `CTYPE*` | Fill underlying data buffer via accessor confirmed in Step 1 |

**Note on object model:** `AER::QV::DensityMatrix<double>` is declared as a local variable in
`bench_aer_dm()`. The object itself is small (pointer + metadata); the internal state buffer
(which is `dim² × sizeof(std::complex<double>)` bytes) is almost certainly heap-allocated by the
constructor (likely via `std::vector` or similar). Verify with the `sizeof` check in the Step 1
probe. Do NOT rely on the state buffer being on the stack — at 12 qubits the density matrix is
256 MB, which would overflow any reasonable stack limit (default 8 MB).

**OpenMP thread contention warning:** `AER::QV::DensityMatrix<double>` inherits from
`QubitVector`, which uses OpenMP internally. If `bench_gate` or `bench_aer_only` is run with
`OMP_NUM_THREADS > 1`, both qlib (through `libq`) and Aer's kernel may compete for the same
thread pool, measuring contended behavior rather than peak single-method throughput. For
controlled comparisons:

```bash
OMP_NUM_THREADS=1 OMP_MAX_ACTIVE_LEVELS=1 ./benchmark/bench_gate
OMP_NUM_THREADS=1 OMP_MAX_ACTIVE_LEVELS=1 ./benchmark/bench_aer_only
```

`OMP_MAX_ACTIVE_LEVELS=1` disables nested parallelism. See also the Notes section at the end
of this plan.

**Exact function signatures (C linkage, matching bench_qulacs.cpp pattern):**

```cpp
#ifdef WITH_AER_DM

#include "simulators/density_matrix/densitymatrix.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <complex>
#include <exception>
#include <random>
#include "bench.h"

extern "C" {

bench_result_t bench_aer_dm(qubit_t qubits, const char *gate_name,
                             int iterations, int warmup, int runs);

} /* extern "C" */

#endif /* WITH_AER_DM */
```

**Exception handling — mandatory:** `AER::QV::DensityMatrix<double>` is a C++ class whose
constructor, gate methods, and destructor can throw (at minimum `std::bad_alloc`). If an
exception propagates through the `extern "C"` boundary into `bench_main.c`, behavior is
undefined. The entire body of `bench_aer_dm()` must be wrapped in a try/catch block:

```cpp
bench_result_t bench_aer_dm(qubit_t qubits, const char *gate_name,
                             int iterations, int warmup, int runs) {
    bench_result_t result = {0};
    result.qubits     = qubits;
    result.dim        = (dim_t)1 << qubits;
    result.gate_name  = gate_name;
    result.method     = "qiskit_aer_dm";
    result.iterations = iterations;
    result.runs       = runs;

    try {
        /* ... all Aer DM logic here ... */
    } catch (const std::exception &e) {
        fprintf(stderr, "bench_aer_dm: exception for %u qubits, gate %s: %s\n",
                (unsigned)qubits, gate_name, e.what());
        return result;
    } catch (...) {
        fprintf(stderr, "bench_aer_dm: unknown exception for %u qubits, gate %s\n",
                (unsigned)qubits, gate_name);
        return result;
    }

    return result;
}
```

**Allocation and initialisation pattern:**

```cpp
// Constructor convention confirmed in Step 1 probe.
// If the argument is num_qubits: AER::QV::DensityMatrix<double> dm(qubits);
// If the argument is num_states: AER::QV::DensityMatrix<double> dm(1 << qubits);
// Use whichever form the probe confirmed.
AER::QV::DensityMatrix<double> dm(qubits);  // replace with confirmed form

// Fill with random values to force page faults and avoid trivial gate no-ops.
// Use a local mt19937 PRNG with a deterministic seed derived from qubits.
// Do NOT use std::rand(): it is not thread-safe and the global state is shared.
std::mt19937_64 rng(static_cast<uint64_t>(qubits) * 6364136223846793005ULL + 1442695040888963407ULL); /* Knuth LCG constants — deterministic seed varies with qubit count */
std::uniform_real_distribution<double> dist(-0.5, 0.5);

// Use the data accessor confirmed in Step 1. If dm.data() returns std::complex<double>*:
size_t dim_val = (size_t)1 << qubits;
size_t n_elems = dim_val * dim_val;
auto *buf = dm.data();  // replace with confirmed accessor; verify it returns std::complex<double>*
for (size_t i = 0; i < n_elems; ++i) {
    buf[i] = std::complex<double>(dist(rng), dist(rng));
}
// If no raw pointer accessor is available, apply random unitaries before the timed region instead.

result.memory_bytes = n_elems * sizeof(std::complex<double>);
```

**Dispatch pattern (1Q):**

```cpp
// Function pointer approach does not work with member functions.
// Use a string comparison + direct call, identical in structure to Qulacs.
// Replace method names below with the confirmed API from Step 1.
double run_times[BENCH_MAX_RUNS];

if (std::strcmp(gate_name, "X") == 0) {
    // warmup
    for (int i = 0; i < warmup; ++i)
        dm.apply_x(i % qubits);  // replace with confirmed method name
    // timed runs
    for (int r = 0; r < runs; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i)
            for (qubit_t t = 0; t < qubits; ++t)
                dm.apply_x(t);  // replace with confirmed method name
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }
} else if (std::strcmp(gate_name, "H") == 0) {
    // use confirmed H gate method or apply_unitary_matrix as determined in Step 1
    ...
} else if (std::strcmp(gate_name, "Z") == 0) {
    ...
}
```

Note: the `int cnt = (runs < 200) ? runs : 200;` clamp from `bench_qulacs.cpp` is not needed
here. `bench_parse_options` already enforces `runs <= BENCH_MAX_RUNS (200)`, so the clamp is
dead code. Use `runs` directly.

Because `AER::QV::DensityMatrix<double>` is a value type (not a raw pointer), it cannot be
passed to a C function pointer. The if/else dispatch structure is identical to the Qulacs
approach in terms of logic; only the call syntax changes from `gate_fn(t, rho, dim)` to
a member call on `dm`.

**Method string:**

```cpp
result.method = "qiskit_aer_dm";
```

**Stats and ops_per_sec computation:** copy verbatim from `bench_qulacs.cpp`. The logic is
identical (bench_compute_stats, sweep over qubits for 1Q, build_all_pairs for 2Q).

**`build_all_pairs` and `bench_compute_stats`:** both are already declared in `bench.h` with C
linkage and will link correctly from a `.cpp` translation unit compiled with `extern "C" {}`.

---

## Step 3a — Update `benchmark/CMakeLists.txt`: Detection + INTERFACE Target

**Estimated effort:** 30–45 minutes
**Depends on:** Step 1
**Blocks:** Step 3b
**Parallel with:** Step 2

### What to do

Add Qiskit-Aer detection and the `aer_dm_static` INTERFACE library target immediately after the
Qulacs block (around line 90).

**Detection block (add after the `endif()` that closes the Qulacs block):**

```cmake
# Qiskit-Aer DM — header-only kernel vendored to extern/aer-dm/
# Vendored from commit: <record the hash from extern/aer-dm/VERSION here>
set(AER_DM_INCLUDE_DIR "${EXTERN_DIR}/aer-dm/include")
if(EXISTS "${AER_DM_INCLUDE_DIR}/simulators/density_matrix/densitymatrix.hpp")
    set(AER_DM_FOUND TRUE)
endif()
```

There is no `.cpp` source to compile for Aer DM (the kernel is header-only). The
`aer_dm_static` target is therefore an `INTERFACE` library that carries the include path and
compile flags. If the kernel proves to require compiled translation units (unlikely; confirm by
inspecting the headers), convert to a `STATIC` library collecting those `.cpp` files.

**Interface target:**

```cmake
if(AER_DM_FOUND)
    add_library(aer_dm_static INTERFACE)
    # Use SYSTEM includes so the compiler treats vendored headers as system headers,
    # suppressing all warnings from them without leaking -w to consuming targets.
    target_include_directories(aer_dm_static SYSTEM INTERFACE "${AER_DM_INCLUDE_DIR}")
    # -march=native is intentionally omitted: bench_gate already receives it transitively
    # through libq's PUBLIC compiler flags in Release builds.
    # Note: CXX_STANDARD cannot be set on an INTERFACE target.
    # Set it on consuming targets (bench_gate, bench_aer_only) instead.
endif()
```

Note: `-w` is intentionally absent. Suppressing warnings from vendored headers is handled by the
`SYSTEM` keyword on `target_include_directories`, which is the correct mechanism and does not
leak warning suppression to consumer translation units (`bench_main.c`, `bench_aer_dm.cpp`, etc.).

---

## Step 3b — Update `benchmark/CMakeLists.txt`: Source List + Target Wiring

**Estimated effort:** 30–45 minutes
**Depends on:** Steps 2 and 3a (file `bench_aer_dm.cpp` must exist; `aer_dm_static` target must exist)
**Blocks:** Step 4

### What to do

**Append to `BENCH_MIXED_SRC` (add after the `if(QULACS_FOUND)` source-append block):**

```cmake
if(AER_DM_FOUND)
    list(APPEND BENCH_MIXED_SRC "${BENCH_SRC_DIR}/bench_aer_dm.cpp")
endif()
```

**In the `bench_gate` target wiring section (add after the `if(QULACS_FOUND)` wiring block):**

```cmake
if(AER_DM_FOUND)
    target_compile_definitions(bench_gate PRIVATE WITH_AER_DM)
    target_link_libraries(bench_gate PRIVATE aer_dm_static)
    # Redundant if Qulacs is also found (Qulacs block already sets CXX_STANDARD 17 on bench_gate);
    # required if only Aer DM is present and Qulacs is absent.
    set_property(TARGET bench_gate PROPERTY CXX_STANDARD 17)
    set_property(TARGET bench_gate PROPERTY CXX_STANDARD_REQUIRED ON)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(bench_gate PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()
```

**Seamless integration confirmed:** When `extern/aer-dm/include/simulators/density_matrix/densitymatrix.hpp`
is present at configure time, `AER_DM_FOUND` is set, `WITH_AER_DM` is defined on `bench_gate`,
and the Aer-DM column appears automatically in all output modes (console, CSV, pgfplots) with no
additional flags or manual steps required. If the directory is absent (e.g., on a machine where
headers have not been copied), the adapter silently drops out and `bench_gate` still builds.

**Update the configuration summary printout (at the bottom of the file):**

```cmake
message(STATUS "Benchmark configuration:")
message(STATUS "  QuEST:     ${QuEST_FOUND}")
message(STATUS "  Qulacs:    ${QULACS_FOUND}")
message(STATUS "  Aer DM:    ${AER_DM_FOUND}")
message(STATUS "  OpenMP:    ${OpenMP_CXX_FOUND}")
```

---

## Step 4 — Update `benchmark/include/bench.h`

**Estimated effort:** 30 minutes
**Depends on:** Steps 2 and 3b (must know function signature and compile guard name)
**Blocks:** Steps 5, 6, 7, 8

### What to do

Add the `bench_aer_dm` declaration in the "External framework benchmarks" section of `bench.h`,
following the exact pattern of `bench_qulacs`.

**Add immediately after the `#endif /* WITH_QULACS */` block (around line 300):**

```c
#ifdef WITH_AER_DM
/**
 * @brief Run Qiskit-Aer density matrix benchmark in-process.
 *
 * AER::QV::DensityMatrix<double> has no global state; this is a direct
 * construct/compute/destruct call. Performs `runs` independent timed rounds
 * and returns the full statistical summary. Exceptions from C++ code are
 * caught internally; a zero time_ms result indicates failure.
 */
bench_result_t bench_aer_dm(qubit_t qubits, const char *gate_name,
                            int iterations, int warmup, int runs);
#endif
```

No new constants are needed. The existing `MAX_PAIRS`, `BENCH_MAX_RUNS`, `BENCH_DEFAULT_*`
constants all apply unchanged.

---

## Step 5 — Update `benchmark/src/bench_main.c`

**Estimated effort:** 2–3 hours
**Depends on:** Step 4
**Blocks:** Step 6

### What to do

Seven distinct edits are required in `bench_main.c`. They are listed in top-to-bottom file
order and labelled Edit A through Edit G to avoid confusion with the top-level Step numbering.

**Edit A — `bench_summary_t` struct (around line 57):** add `aer_dm` result slot:

```c
typedef struct {
    bench_result_t packed;
    bench_result_t tiled;   int has_tiled;
    bench_result_t dense;   int has_dense;
    bench_result_t naive;   int has_naive;
    double naive_scale;
    bench_result_t quest;   int has_quest;
    bench_result_t qulacs;  int has_qulacs;
    bench_result_t aer_dm;  int has_aer_dm;   /* NEW */
} bench_summary_t;
```

**Edit B — `print_speedup_and_memory()` (around line 67):** add `aer_dm` to both the speedup
and memory sections, following the `has_qulacs` pattern exactly.

In the speedup section, note that `print_speedup_and_memory` divides by `s->packed.time_ms`. In
`BENCH_AER_ONLY_MODE`, `s->packed` is zero-initialized (`time_ms == 0`), causing division by
zero. Guard the entire `print_speedup_and_memory` call site in Edit F and Edit G with
`#ifndef BENCH_AER_ONLY_MODE` (see those edits below). The function body itself does not need
guarding; it is unreachable in that mode.

```c
// In speedup section (after the has_qulacs block):
if (s->has_aer_dm && s->aer_dm.time_ms > 0) {
    printf("%s%.1fx vs Aer-DM", has ? ", " : " ", s->aer_dm.time_ms / s->packed.time_ms);
    has = 1;
}

// In max_mem scan (after the has_qulacs line):
if (s->has_aer_dm && s->aer_dm.memory_bytes > max_mem) max_mem = s->aer_dm.memory_bytes;

// In memory print (after the has_qulacs block):
if (s->has_aer_dm && s->aer_dm.memory_bytes > 0)
    printf(", Aer-DM=%.1f %s", (double)s->aer_dm.memory_bytes / scale, unit);
```

**Edit C — `NUM_METHODS` macro and method enum (around line 259–282):**

Change `NUM_METHODS` from 6 to 7 and add `M_AER_DM` to the enum:

```c
#define NUM_METHODS     7   /* qlib_packed, qlib_tiled, blas, naive, quest, qulacs, aer_dm */

enum { M_QLIB = 0, M_QLIB_TILED = 1, M_BLAS = 2, M_NAIVE = 3,
       M_QUEST = 4, M_QULACS = 5, M_AER_DM = 6 };
```

**`pgfplots_data_t` struct — automatic growth:** The struct's timing and memory arrays are
declared as `[NUM_GATES][MAX_QUBIT_CONFIGS][NUM_METHODS]`. Because `NUM_METHODS` is a
compile-time constant used directly in the array dimensions, changing `NUM_METHODS` from 6 to 7
automatically expands every array in `pgfplots_data_t` — no struct definition changes are needed
beyond the `#define` update. The struct itself requires no separate edit.

Note: `method_names[]`, `method_indices[]`, and `num_methods` in `print_pgfplots_output()` are
local arrays NOT driven by `NUM_METHODS` — Edit D updating them is a hard requirement, not
optional cleanup.

**Edit D — `print_pgfplots_output()` method table (around line 364):** add `qiskit_aer_dm` as
the 6th pgfplots column. Note: `M_NAIVE` is intentionally excluded from pgfplots output (naive
runs at reduced iterations and is only meaningful in console/CSV mode). The pgfplots method list
therefore has 6 entries (not 7):

```c
const char *method_names[]   = {"qlib_packed", "qlib_tiled", "blas", "quest", "qulacs", "qiskit_aer_dm"};
int         method_indices[] = {M_QLIB, M_QLIB_TILED, M_BLAS, M_QUEST, M_QULACS, M_AER_DM};
int         num_methods      = 6;
```

This change applies to the single `print_pgfplots_output()` function call site. Because
`method_names` and `method_indices` are iterated with `num_methods`, adding the 6th entry
automatically propagates to all 5 timing tables (mean, ci95, min, median, cv) and the memory
table — no per-table edits are needed.

**Edit E — header banner `printf` (around line 510):** add Aer-DM to the "Comparing:" line
and the enabled-framework status lines. Gate the non-Aer components of the banner with
`#ifndef BENCH_AER_ONLY_MODE` so that in `BENCH_AER_ONLY_MODE`, the banner accurately reflects
only the Aer-DM method running:

```c
#ifndef BENCH_AER_ONLY_MODE
    printf("Comparing: qlib packed");
    printf(" vs qlib tiled");
    printf(" vs BLAS dense");
    printf(" vs naive loop");
  #ifdef WITH_QUEST
    printf(" vs QuEST");
  #endif
  #ifdef WITH_QULACS
    printf(" vs Qulacs");
  #endif
  #ifdef WITH_AER_DM
    printf(" vs Qiskit-Aer-DM");
  #endif
#else /* BENCH_AER_ONLY_MODE */
  #ifdef WITH_AER_DM
    printf("Comparing: Qiskit-Aer-DM only (bench_aer_only mode)\n");
  #endif
#endif /* !BENCH_AER_ONLY_MODE */

// Status lines — similarly gated:
#ifndef BENCH_AER_ONLY_MODE
  #ifdef WITH_QUEST
    printf("QuEST: enabled\n");
  #endif
  #ifdef WITH_QULACS
    printf("Qulacs: enabled\n");
  #endif
#endif
#ifdef WITH_AER_DM
    printf("Qiskit-Aer-DM: enabled\n");
#endif
```

**Edit F — 1Q gate loop (around line 585):** four sub-edits.

(a) **Method call sites** — wrap all non-Aer calls in `#ifndef BENCH_AER_ONLY_MODE`, then add
the Aer call unconditionally (within `#ifdef WITH_AER_DM`):

```c
#ifndef BENCH_AER_ONLY_MODE
    s.packed = bench_qlib_packed(...);
    s.tiled  = bench_qlib_tiled(...);
    s.has_tiled = ...;
    if (run_dense) {
        s.dense = bench_blas_dense(...);
        s.has_dense = 1;
        if (!opts.pgfplots_output) {
            s.naive = bench_naive_loop(...);
            s.has_naive = 1;
        }
    }
  #ifdef WITH_QUEST
    s.quest     = bench_quest(...);
    s.has_quest = ...;
  #endif
  #ifdef WITH_QULACS
    s.qulacs     = bench_qulacs(...);
    s.has_qulacs = ...;
  #endif
#endif /* !BENCH_AER_ONLY_MODE */

#ifdef WITH_AER_DM
    s.aer_dm     = bench_aer_dm(qubits, gates[g].name, opts.iterations, opts.warmup, opts.runs);
    s.has_aer_dm = (s.aer_dm.time_ms > 0);
#endif
```

(b) **pgf storage** (inside `if (pgf && qi < MAX_QUBIT_CONFIGS)`, after the `#ifdef WITH_QULACS`
pgf block):

```c
#ifdef WITH_AER_DM
    pgf->time_ms[g][qi][M_AER_DM]        = s.aer_dm.time_ms;
    pgf->time_ms_std[g][qi][M_AER_DM]    = s.aer_dm.time_ms_std;
    pgf->time_ms_min[g][qi][M_AER_DM]    = s.aer_dm.time_ms_min;
    pgf->time_ms_median[g][qi][M_AER_DM] = s.aer_dm.time_ms_median;
    pgf->time_ms_cv[g][qi][M_AER_DM]     = s.aer_dm.time_ms_cv;
    pgf->memory[g][qi][M_AER_DM]         = s.aer_dm.memory_bytes;
#endif
```

(c) **CSV emit** (inside `if (opts.csv_output)`, after the `#ifdef WITH_QULACS` CSV block):

```c
#ifdef WITH_AER_DM
    if (s.has_aer_dm) { SET_SWEEP(s.aer_dm, sw1q); bench_print_csv(&s.aer_dm); }
#endif
```

(d) **Console emit** (inside `else if (!opts.pgfplots_output)`, after the `#ifdef WITH_QULACS`
console block):

```c
#ifndef BENCH_AER_ONLY_MODE
    bench_print_result(&s.packed, opts.verbose);
    /* ... other non-Aer results ... */
#endif /* !BENCH_AER_ONLY_MODE */

#ifdef WITH_AER_DM
    if (s.has_aer_dm) bench_print_result(&s.aer_dm, opts.verbose);
#endif

/* Speedup and memory summary: divide by s.packed.time_ms, which is zero in
 * BENCH_AER_ONLY_MODE. Skip entirely in that mode to avoid division by zero
 * and a meaningless zero-valued line. */
#ifndef BENCH_AER_ONLY_MODE
    print_speedup_and_memory(&s);
#endif
```

**Edit G — 2Q gate loop (around line 685):** same four sub-edits as Edit F, with two
substitutions throughout: `pg` (= `NUM_1Q_GATES + g`) replaces `g` in all pgf array indices,
and `sw2q` replaces `sw1q` in `SET_SWEEP`. Note: the 2Q loop does not include `naive_loop`,
consistent with the existing benchmark design. The `#ifndef BENCH_AER_ONLY_MODE` and
`#ifndef BENCH_AER_ONLY_MODE` guard on `bench_print_result(&s.packed, ...)` and
`print_speedup_and_memory(&s)` apply identically in the 2Q loop.

(a) **`bench_aer_dm()` call site** — the `#ifndef BENCH_AER_ONLY_MODE` guard covers ALL
non-Aer calls: `bench_qlib_packed_2q`, `bench_qlib_tiled_2q`, `bench_blas_dense_2q`, and the
Quest/Qulacs calls. The Aer call sits OUTSIDE the guard (inside its own `#ifdef WITH_AER_DM`).
The complete guard scope is shown below — do not start the guard at the Quest/Qulacs block; it
must begin at the first qlib call:

```c
#ifndef BENCH_AER_ONLY_MODE
    s.packed = bench_qlib_packed_2q(qubits, gates_2q[g].name, gates_2q[g].fn,
                                    opts.iterations, opts.warmup, opts.runs);

    if (gates_2q[g].has_tiled) {
        s.tiled = bench_qlib_tiled_2q(qubits, gates_2q[g].name, gates_2q[g].fn,
                                      opts.iterations, opts.warmup, opts.runs);
        s.has_tiled = (s.tiled.time_ms > 0);
    }

    if (run_dense) {
        s.dense = bench_blas_dense_2q(qubits, gates_2q[g].name, gates_2q[g].mat,
                                      opts.iterations, opts.warmup, opts.runs);
        s.has_dense = 1;
    }

  #ifdef WITH_QUEST
    s.quest     = bench_quest(qubits, gates_2q[g].name, opts.iterations, opts.warmup, opts.runs);
    s.has_quest = (s.quest.time_ms > 0);
  #endif
  #ifdef WITH_QULACS
    s.qulacs     = bench_qulacs(qubits, gates_2q[g].name, opts.iterations, opts.warmup, opts.runs);
    s.has_qulacs = (s.qulacs.time_ms > 0);
  #endif
#endif /* !BENCH_AER_ONLY_MODE */

#ifdef WITH_AER_DM
    s.aer_dm     = bench_aer_dm(qubits, gates_2q[g].name, opts.iterations, opts.warmup, opts.runs);
    s.has_aer_dm = (s.aer_dm.time_ms > 0);
#endif
```

(b) **pgf storage** (inside `if (pgf && qi < MAX_QUBIT_CONFIGS)`, after the `#ifdef WITH_QULACS`
pgf block; uses `pg` = `NUM_1Q_GATES + g`):

```c
#ifdef WITH_AER_DM
    pgf->time_ms[pg][qi][M_AER_DM]        = s.aer_dm.time_ms;
    pgf->time_ms_std[pg][qi][M_AER_DM]    = s.aer_dm.time_ms_std;
    pgf->time_ms_min[pg][qi][M_AER_DM]    = s.aer_dm.time_ms_min;
    pgf->time_ms_median[pg][qi][M_AER_DM] = s.aer_dm.time_ms_median;
    pgf->time_ms_cv[pg][qi][M_AER_DM]     = s.aer_dm.time_ms_cv;
    pgf->memory[pg][qi][M_AER_DM]         = s.aer_dm.memory_bytes;
#endif
```

(c) **CSV emit** (inside `if (opts.csv_output)`, after the `#ifdef WITH_QULACS` CSV block;
uses `sw2q`):

```c
#ifdef WITH_AER_DM
    if (s.has_aer_dm) { SET_SWEEP(s.aer_dm, sw2q); bench_print_csv(&s.aer_dm); }
#endif
```

(d) **Console emit** (inside `else if (!opts.pgfplots_output)`, after the `#ifdef WITH_QULACS`
console block):

```c
#ifndef BENCH_AER_ONLY_MODE
    bench_print_result(&s.packed, opts.verbose);
    /* ... other non-Aer results ... */
#endif /* !BENCH_AER_ONLY_MODE */

#ifdef WITH_AER_DM
    if (s.has_aer_dm) bench_print_result(&s.aer_dm, opts.verbose);
#endif

#ifndef BENCH_AER_ONLY_MODE
    print_speedup_and_memory(&s);
#endif
```

---

## Step 6 — Update Thesis Notes

**Estimated effort:** 1–2 hours
**Depends on:** Step 5 (adapter design — method names, API surface — is confirmed and validated)
**Blocks:** nothing (leaf step)
**File:** `~/Projects/thesis/2.Simulation/notes/BENCHMARK_THESIS_NOTES.md`

### What to add

Append a new section to the thesis notes describing the Qiskit-Aer adapter. The section must
address:

1. **What is extracted:** `AER::QV::DensityMatrix<double>` from the `src/simulators/density_matrix/`
   and `src/framework/` directories of the Qiskit-Aer source tree. Approximately 15–30 C++ headers
   are vendored to `extern/aer-dm/include/` (exact count determined during Step 1). No Python,
   no pybind11, no Qiskit runtime.

2. **Why this is fair:** The extracted kernel is the actual computation core that the Qiskit
   runtime calls when executing circuits on the `AerSimulator` with `method='density_matrix'`.
   Benchmarking it directly removes Python dispatch overhead, making the comparison with qlib's
   C kernel and Qulacs's csim kernel a like-for-like measurement of the underlying
   density-matrix arithmetic.

3. **What the method name means:** The column labelled `qiskit_aer_dm` in all output modes
   (console, CSV, pgfplots) reports throughput of Aer's internal C++ density matrix kernel,
   compiled and invoked directly without Python, Qiskit-Terra, or any other Qiskit infrastructure.

4. **Comparison with Qulacs and QuEST adapters:**
   - Qulacs (`qulacs`): vendored C++ sources from `extern/qulacs/src/csim/`, compiled as
     `csim_static`, called via C-style function pointers.
   - QuEST (`quest`): pre-built `.dylib`/`.so` at `extern/QuEST/`, requires `initQuESTEnv` /
     `finalizeQuESTEnv` around the session (init-once pattern).
   - Aer DM (`qiskit_aer_dm`): vendored header-only kernel, no global init/finalize, constructed
     and destroyed per `bench_aer_dm()` call (same no-global-state property as Qulacs).

5. **Thesis caveat:** State that results for qiskit_aer_dm were collected using the standalone
   `bench_aer_only` binary (see Step 7) because the full `bench_gate` run for other competitors
   had already been completed; merging the results is straightforward since both binaries emit
   the same CSV schema.

---

## Step 7 — Add `bench_aer_only` CMake Target

**Estimated effort:** 1 hour
**Depends on:** Step 5 (the `#ifndef BENCH_AER_ONLY_MODE` guards in `bench_main.c` must exist before `bench_aer_only` produces correct output)
**Blocks:** nothing (leaf step)

> **Lifecycle note:** This target is a one-shot data-collection tool. Once Qiskit-Aer data has
> been merged with the existing `bench_gate` CSV output, the `bench_aer_only` target should be
> removed from `CMakeLists.txt`. Do not commit the binary.

### Rationale

`bench_gate` runs all competitors (qlib, BLAS, naive, QuEST, Qulacs, Aer-DM) sequentially.
At 12 qubits with default parameters, a full run takes up to 24 hours. Since data for qlib,
BLAS, naive, QuEST, and Qulacs is already collected, re-running those methods wastes wall-clock
time. `bench_aer_only` runs only the Aer-DM method against the same gates and qubit range,
produces the same output in all modes, and can be merged with existing data in post-processing.

### Implementation approach — compile-definition gating (no new source file)

`bench_aer_only` shares **all** source files with `bench_gate`. No new `.c` or `.cpp` file is
created. The CMake target compiles the identical source list but adds
`-DBENCH_AER_ONLY_MODE=1`, which causes `bench_main.c` to skip all method calls except
`bench_aer_dm()`.

The `#ifndef BENCH_AER_ONLY_MODE` guards in `bench_main.c` are documented in Step 5 (Edits F
and G) and must be in place before `bench_aer_only` produces correct output. The CMake
additions in this step can be written in parallel with Step 5 (the CMake file does not reference
the source guards), but the resulting binary is not testable until Step 5 is complete.

### `bench_mixed.c` and `bench_baselines.c` link dependency

`bench_aer_only` compiles `bench_mixed.c` and `bench_baselines.c` as part of `BENCH_MIXED_SRC`.
These files contain function bodies that reference qlib and BLAS symbols unconditionally (the
function bodies are not gated by `BENCH_AER_ONLY_MODE`). Two strategies exist:

**Option (a) — accept the `libq` link dependency (chosen approach):** Link `bench_aer_only`
against `q` (and transitively BLAS). This is the simpler path: the link works correctly, the
binary is slightly larger, and the existing source files require no modification beyond the call
site guards in Step 5. This is the approach used in the CMake wiring below.

**Option (b) — guard function bodies with `#ifndef BENCH_AER_ONLY_MODE` (not chosen):**
Wrap the bodies of functions in `bench_mixed.c` and `bench_baselines.c` that reference qlib/BLAS
symbols with `#ifndef BENCH_AER_ONLY_MODE` stubs. This would allow dropping the `libq` and BLAS
link dependency from `bench_aer_only`. The trade-off is significant additional source
modifications that increase maintenance burden. Option (a) is preferred unless the `libq`
dependency becomes unavailable or causes build problems on the target machine.

### CMake wiring for `bench_aer_only` (add at the bottom of `benchmark/CMakeLists.txt`)

```cmake
#=======================================================================================================================
# bench_aer_only — Qiskit-Aer data-collection tool (temporary; delete after data collection)
#=======================================================================================================================
if(AER_DM_FOUND)
    # Same source list as bench_gate — no new source file needed.
    # BENCH_AER_ONLY_MODE=1 gates out all non-Aer method calls in bench_main.c.
    add_executable(bench_aer_only ${BENCH_MIXED_SRC})
    set_target_properties(bench_aer_only PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/benchmark"
        EXCLUDE_FROM_ALL TRUE
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
    )
    target_include_directories(bench_aer_only PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
    )
    # WITH_AER_DM enables bench_aer_dm() calls; BENCH_AER_ONLY_MODE skips all other methods.
    target_compile_definitions(bench_aer_only PRIVATE WITH_AER_DM BENCH_AER_ONLY_MODE=1)
    # q is required: bench_mixed.c and bench_baselines.c function bodies reference qlib/BLAS
    # symbols unconditionally. Option (a) accepted — see Step 7 rationale.
    target_link_libraries(bench_aer_only PRIVATE q aer_dm_static)
    # Note: QuEST and Qulacs link dependencies are intentionally omitted here.
    # Without WITH_QUEST or WITH_QULACS defined, bench_quest.c and bench_qulacs.cpp compile
    # to empty translation units with no external symbol references, so no link dependency is
    # needed. Adding those libraries would be dead weight.
    # Note: compiling bench_quest.c and bench_qulacs.cpp without their WITH_QUEST/WITH_QULACS
    # guards may produce an 'empty translation unit' compiler warning; this is cosmetic and
    # does not affect the build.
    if(OpenMP_CXX_FOUND)
        target_link_libraries(bench_aer_only PRIVATE OpenMP::OpenMP_CXX)
    endif()
    message(STATUS "  bench_aer_only: enabled (temporary data-collection tool)")
endif()
```

### Build commands

```bash
cmake -B cmake-build-release -DCMAKE_BUILD_TYPE=Release -DBUILD_BENCHMARKS=ON
cmake --build cmake-build-release --target bench_aer_only

# Console mode (default):
OMP_NUM_THREADS=1 OMP_MAX_ACTIVE_LEVELS=1 \
    ./cmake-build-release/benchmark/bench_aer_only

# CSV mode:
OMP_NUM_THREADS=1 OMP_MAX_ACTIVE_LEVELS=1 \
    ./cmake-build-release/benchmark/bench_aer_only --csv > aer_dm_results.csv

# pgfplots mode:
OMP_NUM_THREADS=1 OMP_MAX_ACTIVE_LEVELS=1 \
    ./cmake-build-release/benchmark/bench_aer_only --pgfplots > aer_dm_pgf.dat

# Verbose console mode:
OMP_NUM_THREADS=1 OMP_MAX_ACTIVE_LEVELS=1 \
    ./cmake-build-release/benchmark/bench_aer_only --verbose

# Custom range:
OMP_NUM_THREADS=1 OMP_MAX_ACTIVE_LEVELS=1 \
    ./cmake-build-release/benchmark/bench_aer_only \
    --min-qubits 8 --max-qubits 12 --step 1 \
    --runs 20 --warmup 200 --iterations 5000 --csv
```

All flags listed above are identical to `bench_gate` — they come from the same
`bench_parse_options()` call. Run with `OMP_NUM_THREADS=1 OMP_MAX_ACTIVE_LEVELS=1` for
controlled single-threaded measurements that avoid OpenMP thread contention between qlib and
Aer's internal parallelism.

### Merging CSV data

CSV is the canonical merge format. The output CSV of `bench_aer_only --csv` has the identical
header and row format as `bench_gate --csv`. Concatenate after removing the header from the
second file:

```bash
# Assuming existing data is in bench_gate_results.csv:
tail -n +2 aer_dm_results.csv >> bench_gate_results.csv
```

This is sufficient for all downstream analysis. If pgfplots figures are generated from CSV
(rather than the raw pgfplots dat files), no further merge step is needed.

### Merging pgfplots data

pgfplots output from `bench_aer_only --pgfplots` emits `nan` for all non-Aer columns
(`qlib_packed`, `qlib_tiled`, `blas`, `quest`, `qulacs`). This output is **not** directly
concatenatable with `bench_gate --pgfplots` output; the Aer column must be extracted and joined
column-wise with the existing pgfplots data.

**Recommended approach:** use Python/pandas to load both files, join on the `(gate, qubits)` key
columns, and write the merged result. This is robust to comment lines, variable whitespace, and
the multi-table structure emitted by `print_pgfplots_output()`.

If pgfplots merging is not required (e.g., thesis figures are generated from CSV), skip this
step entirely and use the CSV merge procedure above.

---

## Step 8 — Update `docs/BENCHMARK.md`

**Estimated effort:** 1 hour
**Depends on:** Step 4
**Blocks:** nothing (leaf step)

### What to change

**File layout table:** add `bench_aer_dm.cpp` to the `src/` listing:

```
│   ├── bench_aer_dm.cpp        # Qiskit-Aer DM adapter (compiled only when WITH_AER_DM is defined)
```

Also add `extern/aer-dm/` to the repo layout description if relevant.

**Dependencies table:** add a row for Aer-DM:

| Dependency | Location | Detection |
|---|---|---|
| Qiskit-Aer DM kernel | `extern/aer-dm/include/` | Presence of `simulators/density_matrix/densitymatrix.hpp` |

**CMake targets table:** add:

| Target | Type | Condition |
|---|---|---|
| `aer_dm_static` | Interface library | `AER_DM_FOUND` — carries SYSTEM include path only (`-march=native` is omitted; consumers receive it transitively through `libq`) |
| `bench_aer_only` | Executable | `AER_DM_FOUND` (Release only; temporary — delete after data collection) |

**Compile definitions:** add:

> `WITH_AER_DM` is defined on `bench_gate` and `bench_aer_only` when the Aer DM headers are
> found at `extern/aer-dm/include/`. Set automatically; no manual flag required. When headers are
> present, the `qiskit_aer_dm` column appears in all output modes (console, CSV, pgfplots)
> without any additional configuration.
>
> `BENCH_AER_ONLY_MODE` is defined only on `bench_aer_only`. It gates out all non-Aer method
> call sites in `bench_main.c` while preserving the full CLI, output formatting, and statistical
> infrastructure.

**Configuration summary section:** update the example output:

```
-- Benchmark configuration:
--   QuEST:     TRUE/FALSE
--   Qulacs:    TRUE/FALSE
--   Aer DM:    TRUE/FALSE
--   OpenMP:    TRUE/FALSE
```

**Implementations Benchmarked table:** add row 7 (enum index 6) and update the `NUM_METHODS`
note to 7:

| Index | Enum | Method string | Description | Source |
|---|---|---|---|---|
| 6 | `M_AER_DM` | `qiskit_aer_dm` | Qiskit-Aer internal C++ density matrix kernel | `bench_aer_dm.cpp` |

Update the introductory sentence: "Seven method slots are defined (indices 0–6; enum in
`bench_main.c`). Six appear in pgfplots output (naive is excluded); all seven appear in CSV and
console output.".

**pgfplots output section:** update the method column list. All five timing tables (mean,
ci95, min, median, cv) and the memory table now emit 6 columns:
`qlib_packed`, `qlib_tiled`, `blas`, `quest`, `qulacs`, `qiskit_aer_dm`.
The `naive` method is excluded from pgfplots output (as it has always been).

**External framework adapters section:** add a paragraph for Aer DM following the Qulacs
paragraph:

> **Qiskit-Aer DM** runs directly in-process. `bench_aer_dm()` constructs an
> `AER::QV::DensityMatrix<double>` (whose internal state buffer is heap-allocated by the
> constructor), fills its internal buffer with non-trivial values to fault all pages before the
> timed region, runs warmup + timed rounds, then destroys the object. No global state; no
> init/finalize. The kernel is the same C++ code executed by the Qiskit `AerSimulator` for
> `method='density_matrix'` workloads, extracted and compiled without Python infrastructure.
> Exceptions from C++ are caught at the `extern "C"` boundary; a zero `time_ms` result indicates
> a failure (e.g., `std::bad_alloc` at large qubit counts).

**Memory reporting table:** add row:

| Method | Measurement |
|---|---|
| `qiskit_aer_dm` | `dim² × sizeof(std::complex<double>)` |

**Status section:** update "Known gaps" — remove Qiskit-Aer from any future-work mention if
present. Change status note to reflect that three external competitors (QuEST, Qulacs, Aer-DM)
are complete.

---

## Notes for All Steps

- The Aer DM adapter is purely a **density-matrix** competitor. Qiskit-Aer also has statevector
  and matrix-product state backends; these are out of scope for this plan.
- All gate operations on `AER::QV::DensityMatrix<double>` apply the gate in the Schrödinger
  picture, consistent with the other density-matrix competitors (QuEST, Qulacs). The mathematical
  operation is `ρ → UρU†` for 1Q gates and the analogous 2Q form for CX and SWAP.
- `AER_DM_FOUND` is controlled entirely by the presence of the vendored headers in
  `extern/aer-dm/include/`. If the directory is absent (e.g., on a machine where headers have
  not been copied), the adapter silently drops out and the binary still builds.
- The `bench_aer_only` target is marked `EXCLUDE_FROM_ALL TRUE`, so it is not built by default.
  It must be requested explicitly: `cmake --build ... --target bench_aer_only`.
- `bench_aer_only` and `bench_gate` produce identical output for the `qiskit_aer_dm` rows.
  pgfplots output from `bench_aer_only` will print `nan` for all non-Aer columns, which is
  correct and expected. Use CSV as the canonical merge format (see Step 7).
- **OpenMP thread contention:** `AER::QV::DensityMatrix<double>` (via its `QubitVector` base
  class) uses OpenMP internally. Running `bench_gate` or `bench_aer_only` with `OMP_NUM_THREADS
  > 1` may cause thread pool contention between qlib's OpenMP parallelism and Aer's internal
  parallelism, producing measurements that reflect contended behavior rather than peak
  single-method throughput. Always run with `OMP_NUM_THREADS=1 OMP_MAX_ACTIVE_LEVELS=1` for
  controlled, reproducible comparisons. `OMP_MAX_ACTIVE_LEVELS=1` disables nested parallelism.
