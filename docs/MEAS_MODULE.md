# MEAS Module: Expectation Value, Gradient, and Moment Matrix

**Module:** Measurement and Gradient Computation
**Status:** Design complete — implementation pending
**Last Updated:** 2026-02-23

---

## Scope

This module provides three computational deliverables, all restricted to
**diagonal cost Hamiltonians** (combinatorial optimisation problems):

| Deliverable       | Ideal circuit | Noisy circuit | Requirement              |
|-------------------|:-------------:|:-------------:|--------------------------|
| Expectation value | ✓             | ✓             | —                        |
| Gradient          | ✓             | ✗             | Adjoint differentiation  |
| Moment matrix     | ✓             | ✗             | Pure initial state only  |

**Gradient computation is not provided for noisy circuits.** Noise is
non-invertible; there is no backward pass. If an external optimiser requires
gradients for a noisy circuit it must call `expectation_value` repeatedly using
the parameter-shift rule or finite differences.

State preparation — including LCU-based state preparation — is the caller's
responsibility. This module receives a `state_t *` initial state and applies
no knowledge of how it was produced.

---

## Preliminaries

The following work in other modules must be completed before any part of this
module can be implemented. Items are ordered by dependency and criticality.

### 1. Gate Module — Diagonal Phase Multiply *(critical path)*

A new operation must be added to the gate module across all four backends
(pure/packed, pure/tiled, mixed/packed, mixed/tiled):

```c
void gate_diag_phase(state_t *state, const double *h, double gamma);
```

Semantics:

- **Pure:** `ψ(x) ← ψ(x) · exp(−iγ h_x)` for all x in 0 … 2^n − 1
- **Mixed:** `ρ_{xy} ← ρ_{xy} · exp(−iγ (h_x − h_y))` for all x, y

This is the only unresolved gate dependency. Without it the cost channel
(`qaoa.c`), the QAOA circuit constructor, and the adjoint gradient cannot be
implemented. `expectation_value` and `gradient_adjoint` with a mixer-only
circuit could be written first but would be of limited utility.

### 2. Shared Type — `unitary_fn_t`

The type used for unitaries in `lcu_layer_t` (see Data Structures) must be
defined in a shared header (`q_types.h` or `gate.h`) before the meas data
structures can be declared:

```c
// Applies a single unitary U in-place to *state.
// All parameters are baked in at construction time by the caller.
typedef void (*unitary_fn_t)(state_t *state);
```

### 3. Pauli Module — Pauli Exponential Unitary *(blocks moment matrix)*

`lcu_layer_t` supports unitaries of the form exp(−iθ P) for a Pauli string P.
This is delivered by a separate module not yet fully specified. The delivery
contract required by this module is:

- A constructor that returns a `unitary_fn_t` with angle θ and Pauli string P
  baked in.
- The returned function applies exp(−iθ P) to the input state in place.

`moment_matrix` implementation must be deferred until this interface is
published, even if the Pauli module itself is incomplete.

### 4. Build System — `ENABLE_MPI` Option *(blocks moment matrix)*

The `ENABLE_MPI` CMake option referenced in the Build Integration section does
not yet exist in `CMakeLists.txt` or `CMakePresets.json`. It must be added
before `moment_matrix.c` can be compiled with MPI support. The single-process
(OpenMP-only) path can be compiled without this option.

### 5. Architecture Decision — MPI Write Protocol *(blocks moment matrix)*

The `moment_matrix` write path requires a decision on distribution and
serialisation strategy before `moment_matrix.c` is written. The three options
are:

| Option | Notes |
|--------|-------|
| OpenMP-only, single process | Lowest risk; defers MPI entirely |
| MPI gather to rank 0, then sequential write | Simple MPI; rank 0 must hold full lower triangle |
| MPI-IO (`MPI_File_write_at`) | Each rank writes its own rows; no gather needed |

**Decision must be recorded here before implementation begins.**

---

## Cost Hamiltonian

H_C is restricted to **diagonal Hamiltonians in the Z basis**:

    H_C = Σ_x h_x |x⟩⟨x|

`h` is stored as a real vector of length 2^n (a `double` array). All other
Hamiltonian representations are out of scope.

Expectation values follow directly from the diagonal:

    Pure state:  ⟨H_C⟩ = Σ_x h_x |ψ(x)|²      (weighted sum over populations)
    Mixed state: ⟨H_C⟩ = Σ_x h_x ρ_{xx}         (weighted sum over diagonal of ρ)

**Co-state initialisation** — H_C|ψ_f⟩ is element-wise multiplication by h_x,
costing O(2^n) with no matrix multiply.

**Mixed-state backward pass** — At the start of Pass 2 (see below), H_f is
initialised from h_x and immediately expanded into the same MIXED_PACKED or
MIXED_TILED memory layout as the density matrix, with all off-diagonal elements
set to zero. Subsequent unitary gates fill in the off-diagonal entries in place.

---

## Expectation Value

A single forward pass through the circuit, updating one state buffer in place.

**Pure state:**
```
tmp = copy(initial)
for k = 1 to N_ch:
    tmp = apply(channel[k], tmp, params[k])
return Σ_x h[x] * |tmp[x]|²
```

**Mixed state:** identical, with `Σ_x h[x] * tmp[x][x]` as the final sum.

For noisy circuits `apply` includes the Kraus step; the algorithm is otherwise
identical. Memory: one state-sized buffer throughout.

---

## Gradient: Adjoint Differentiation (Ideal Circuits Only)

### Mathematical Foundation

For a parametrised gate U_k(θ_k) = exp(−iθ_k G_k) at position k in a circuit
U_f · U_k · V:

    ∂⟨H_C⟩/∂θ_k = −2 Im Tr(H_f^(k) · G_k · ρ_k)

where:
- ρ_k = V ρ_0 V†  is the state immediately before gate k  (Schrödinger picture)
- H_f^(k) = U_f† H_C U_f  is H_C pulled back to position k (Heisenberg picture)

For **pure states** both ρ_k and H_f^(k) are represented as 2^n-dimensional
complex vectors; the trace reduces to a vector inner product. H_f^(k) is never
expanded into a matrix. For **mixed states** both are 2^n × 2^n Hermitian
operators stored in the density matrix memory layout.

### Parametrised Channel

The fundamental unit is a **parametrised channel**: all gates sharing one scalar
parameter θ, grouped together. Each channel exposes three operations:

```c
typedef struct param_channel {
    void   (*apply)     (state_t *state, double theta, const void *ctx);
    double (*grad_inner)(const state_t *lambda, const state_t *tmp,
                         double theta, const void *ctx);
    const void *ctx;        // channel-specific data (qubit indices, Hamiltonian, …)
    int    has_gradient;    // 0 for fixed non-parametrised channels
} param_channel_t;
```

- **`apply(state, θ)`** — Forward unitary evolution: state → U(θ) state U(θ)†.
  For noisy channels includes the Kraus step after the unitary. The unitary
  adjoint needed by Pass 2 is obtained by calling `apply(state, −θ)`, which is
  equivalent for all channel types appearing in ideal circuits. A separate
  adjoint callback is not needed: the Kraus adjoint ε†(A) = Σ_k K†_k A K_k is
  structurally different and cannot be obtained by negating a parameter, but
  Pass 2 is never invoked for noisy circuits.
- **`grad_inner(λ, tmp, θ)`** — Returns −2 Im ⟨λ|G_k|tmp⟩ (pure) or
  −2 Im Tr(H_f · G_k · σ) (mixed), the gradient contribution of this channel.
  Never called for noisy circuits.

Parameters are **not** shared across channels. Each channel owns exactly one
scalar parameter.

### Circuit

```c
typedef struct {
    param_channel_t **channels;
    int               n_channels;
    double           *params;     // length n_channels, one parameter per channel
    int               is_noisy;   // 0: ideal (gradient supported), 1: noisy (forward only)
} pqc_t;
```

A circuit is **either entirely ideal or entirely noisy**. Mixed circuits are not
supported. For ideal circuits, `gradient_adjoint` uses the 3-pass algorithm
below. For noisy circuits, only `expectation_value` is available.

### 3-Pass Algorithm

Exploits unitarity U_k U_k† = I to carry both the state and co-state through the
circuit without storing any intermediate states. Total cost: 3 × N_ch channel
applications, two state-sized buffers.

**Pass 1 — forward, compute final state and expectation value:**
```
tmp = copy(initial)
for k = 1 to N_ch:
    apply(channel[k], tmp, +params[k])

E = Σ_x h[x] * |tmp[x]|²     (pure)
E = Σ_x h[x] * tmp[x][x]     (mixed)
```

**Initialise co-state from final state:**
```
mu = h[x] * tmp[x]            (pure:  element-wise multiply, gives vector)
mu = expand h[x] into DM layout, then propagate as H_f    (mixed)
```

**Pass 2 — backward, pull co-state to circuit entrance:**
```
for k = N_ch, N_ch-1, ..., 1:
    apply(channel[k], mu, -params[k])    # unitary adjoint = apply with negated theta
mu_0 = mu
```

**Pass 3 — forward, accumulate gradients:**
```
tmp = copy(initial)
λ   = mu_0
for k = 1 to N_ch:
    apply(channel[k], tmp, +params[k])
    apply(channel[k], λ,   +params[k])
    if channel[k].has_gradient:
        grad[k] = grad_inner(channel[k], λ, tmp, params[k])
```

At the point of the `grad_inner` call, `tmp` holds the post-gate forward state
ρ_k (Schrödinger picture) and `λ` holds the co-state H_f^(k)|ψ_k⟩ (pure) or
H_f^(k) (mixed) at position k.

---

## QAOA: First Implementation Target

### Circuit Structure

Depth-p QAOA on n qubits: 2p channels, 2p parameters.

```
Channel 1:    cost layer,  param γ_1
Channel 2:    mixer layer, param β_1
  ...
Channel 2p-1: cost layer,  param γ_p
Channel 2p:   mixer layer, param β_p
```

### Cost Channel

Implements U_C(γ) = exp(−iγ H_C).

Since H_C is diagonal, this is **element-wise phase multiplication** in O(2^n):

    Pure:  ψ(x)   ←  ψ(x) · exp(−iγ h_x)
    Mixed: ρ_{xy} ←  ρ_{xy} · exp(−iγ (h_x − h_y))

Generator G_C = H_C. Gradient inner product:

    −2 Im ⟨μ|H_C|tmp⟩  =  −2 Im Σ_x h_x μ*(x) tmp(x)          (pure)
    −2 Im Tr(H_f H_C σ) =  −2 Im Σ_x h_x [H_f · σ]_{xx}        (mixed)

No RZZ decomposition required. No Trotter approximation.

### Mixer Channel

Implements U_B(β) = Π_j exp(−iβ X_j / 2), applied as n sequential Rx(β) calls
(all commute: [X_i, X_j] = 0 for i ≠ j).

Generator G_B = Σ_j X_j / 2. Gradient inner product:

    −2 Im ⟨μ|G_B|tmp⟩  =  Σ_j (−Im ⟨μ|X_j|tmp⟩)

Computed with a single temporary state copy: for each qubit j, apply X_j to the
copy of `tmp`, accumulate −Im ⟨μ | copy_j⟩ into the running sum, restore copy.
This costs one temporary buffer and n sequential X applications regardless of
state type.

For mixed states the same structure applies with the Hilbert-Schmidt inner
product: Tr(H_f · X_j · σ) computed by applying X_j to a copy of σ and
contracting with H_f.

### Initial State

|ψ_0⟩ = |+⟩^⊗n (uniform superposition) via `state_plus()` — already implemented.

### Required Gates

| Gate / Operation                      | Gate module status     | Needed for    |
|---------------------------------------|------------------------|---------------|
| Rx(θ)                                 | ✓ All backends         | Mixer channel |
| Diagonal phase multiply exp(−iγ h_x)  | ✗ New operation needed | Cost channel  |

The cost channel must be implemented as a direct diagonal phase multiply, not as
decomposed RZZ gates. The new gate operation must be added to the gate module
before the cost channel can be implemented.

---

## Moment Matrix

### Definition

Given L layers each offering K_k unitaries {U_{k;1}, …, U_{k;K_k}}, a
multi-index I = (i_1, …, i_L) selects one unitary per layer and defines the
composed unitary U_I = U_{L;i_L} · … · U_{1;i_1}. The moment matrix is:

    M_{IJ} = Tr(H_C · U_I ρ_0 U_J†)

M is N × N (N = ΠK_k), Hermitian (M_{IJ} = M_{JI}*), and computed once for
a given (ρ_0, H_C, circuit). It is stored to file and reused by the external
σ-optimisation library (SDP solver or otherwise) without recomputation.

### Pure Initial State

For |ψ_0⟩ pure, define |φ_I⟩ = U_I|ψ_0⟩. Then:

    M_{IJ} = ⟨φ_J|H_C|φ_I⟩ = Σ_x h_x φ_I(x) φ_J*(x)

The diagonal M_{II} = Σ_x h_x |φ_I(x)|² is real. Mixed initial states are out
of scope (see Deferred).

### Unitaries

Each U_{k;i} is either a gate sequence built from the gate module or an
exponential of a Pauli string from the Pauli module.

### Computation Algorithm (MPI)

The matrix is computed column by column. For each column J:

```
1. Root rank computes |φ_J⟩ = U_J|ψ_0⟩ and broadcasts it to all ranks.
2. Each rank iterates its assigned rows I ≥ J:
       compute |φ_I⟩ = U_I|ψ_0⟩
       M[I][J] = Σ_x h[x] * φ_I[x] * conj(φ_J[x])
       M[J][I] = conj(M[I][J])                      # Hermitian symmetry
3. Ranks write their portion of the lower triangle to the output file.
```

Memory per rank: two state vectors (|φ_I⟩ and |φ_J⟩) plus the partial M
region. For shared-memory systems OpenMP parallelism over rows within one rank
is sufficient.

### File Format

The output file stores the lower triangle of M as raw binary with a fixed
header. All multi-byte values are little-endian.

**Header (fixed layout):**

| Field         | Type       | Description                                  |
|---------------|------------|----------------------------------------------|
| `magic`       | `uint32_t` | `0x4D4D4154` ("MMAT")                        |
| `version`     | `uint32_t` | Format version, currently 1                  |
| `n_qubits`    | `int32_t`  | Number of qubits n                           |
| `L`           | `int32_t`  | Number of layers                             |
| `K[0..L-1]`  | `int32_t`  | Unitaries per layer                          |
| `N`           | `int64_t`  | Total multi-indices N = ΠK_k                 |
| `timestamp`   | `int64_t`  | Unix timestamp of computation                |
| `h[0..2^n-1]`| `double`   | Hamiltonian diagonal, length 2^n             |

**Data:** Lower triangle of M in row-major order (I ≥ J), stored as
`double complex` pairs (real, imag). Length: N*(N+1)/2 elements.

**Writing:** Header fields must be written individually (not via `fwrite` of a
struct) to avoid compiler-inserted padding between fields of different sizes.
Data elements are written as contiguous `double` pairs (real, imag) in a
single `fwrite` call per row segment.

The caller is responsible for converting to the format required by the SDP
solver (e.g., sparse (i, j, value) triples for MOSEK's C API).

---

## Data Structures

```c
// Diagonal cost Hamiltonian: h[x] for x = 0, …, 2^n_qubits − 1
typedef struct {
    double  *h;
    int      n_qubits;
} diag_hamiltonian_t;

// Abstract parametrised channel (grad_inner never called for noisy circuits)
typedef struct param_channel {
    void   (*apply)     (state_t *state, double theta, const void *ctx);
    double (*grad_inner)(const state_t *lambda, const state_t *tmp,
                         double theta, const void *ctx);
    const void *ctx;
    int    has_gradient;    // 0 for fixed non-parametrised channels
} param_channel_t;

// Flat parametrised circuit.
// params is caller-owned; pqc_free does NOT free it.
typedef struct {
    param_channel_t **channels;
    int               n_channels;
    double           *params;      // length n_channels, caller-owned
    int               is_noisy;    // 0: ideal, 1: noisy (expectation value only)
} pqc_t;

// Function pointer type for a single unitary in a moment-matrix circuit.
// Applies U in-place to *state. All parameters are baked in at construction
// time by the caller. Defined in q_types.h (see Preliminary 2).
typedef void (*unitary_fn_t)(state_t *state);

// One layer of a moment-matrix circuit: K_k unitaries.
// Each unitary is a gate sequence (gate module) or Pauli exponential
// (Pauli module, see Preliminary 3), exposed as a unitary_fn_t.
typedef struct {
    unitary_fn_t *unitaries;   // array of K function pointers
    int           K;
} lcu_layer_t;

// Full moment-matrix circuit
typedef struct {
    lcu_layer_t *layers;
    int          L;
    int          N;    // ΠK_k, total number of multi-indices
} lcu_circuit_t;
```

---

## API

**Error handling convention** (consistent with the rest of the library):
programming errors (NULL pointers, mismatched `n_qubits`, `gradient_adjoint`
called with `is_noisy == 1`) trigger `assert` with a descriptive message.
Resource failures (allocation, file I/O) return `NULL` or write nothing and
set `errno`.

```c
// Forward pass only. Works for ideal and noisy circuits.
double expectation_value(
    const pqc_t              *circuit,
    const state_t            *initial,
    const diag_hamiltonian_t *H_C
);

// 3-pass adjoint differentiation. Ideal circuits only (circuit->is_noisy == 0).
// Returns ⟨H_C⟩. Writes ∂⟨H_C⟩/∂params[c] into grad_out[c] for each channel c.
double gradient_adjoint(
    const pqc_t              *circuit,
    const state_t            *initial,
    const diag_hamiltonian_t *H_C,
    double                   *grad_out   // length circuit->n_channels
);

// QAOA circuit constructors
pqc_t *qaoa_circuit_new(
    int                       p,          // number of layers
    const diag_hamiltonian_t *H_C,
    int                       n_qubits
);
void pqc_free(pqc_t *circuit);

// Compute and write the moment matrix (pure initial state, MPI-parallel).
// Writes lower triangle + header to output_path in the binary format above.
void moment_matrix(
    const lcu_circuit_t      *circuit,
    const state_t            *initial,    // must be PURE
    const diag_hamiltonian_t *H_C,
    const char               *output_path
);

// Load a previously written moment matrix from file.
// Allocates *M_out (caller must free). Sets *N_out = circuit N.
void moment_matrix_load(
    const char   *path,
    int          *N_out,
    double complex **M_out    // lower triangle, length N*(N+1)/2
);
```

---

## Build Integration

The measurement module is a static library linked against the state and gate
modules. When moment matrix computation is enabled, it additionally requires
MPI.

Planned `CMakeLists.txt` additions (inside `QSim/`):

```cmake
add_library(meas STATIC
    src/meas/expectation.c
    src/meas/gradient.c
    src/meas/moment_matrix.c
    src/meas/qaoa.c
)
target_include_directories(meas PUBLIC include)
target_link_libraries(meas PRIVATE state gate)

if(ENABLE_MPI)
    find_package(MPI REQUIRED)
    target_link_libraries(meas PRIVATE MPI::MPI_C)
    target_compile_definitions(meas PRIVATE MEAS_MPI)
endif()
```

`ENABLE_MPI` defaults to `OFF` and is set independently of `ENABLE_OPENMP`.
When `ENABLE_MPI` is off, `moment_matrix` runs single-process with OpenMP
parallelism over rows. When on, MPI is used for cross-node distribution with
OpenMP within each rank.

---

## Deferred — Out of Scope for First Implementation

- **Gradient for noisy circuits** — caller uses parameter-shift or finite
  differences via repeated `expectation_value` calls
- **Mixed initial state moment matrix** — requires Hilbert-Schmidt (Liouville)
  space; deferred indefinitely
- **Noise in moment matrix** — moot given pure-only initial state
- **Hessian / second-order gradients**
- **Natural gradient** (quantum geometric tensor / Fisher information metric)
- **General non-diagonal observables**
- **Multi-parameter channels** (more than one independent parameter per channel)
- **Trotter-based Hamiltonian simulation**
- **General Pauli-string Hamiltonian evolution** (non-QAOA; no decision on
  Trotterisation vs. direct exponentiation)
- **RSWAP and other non-standard parametrised gates**
