# Design Specification: Quantum State Simulation Library

## 1. Purpose and Scope

### 1.1 Purpose

A C library for studying the **performance and scalability of quantum heuristics**, particularly QAOA and related variational algorithms. This is a research tool (emulator), not a hardware simulator — it directly computes quantities that would be difficult on actual quantum hardware.

### 1.2 Primary Use Cases

- Classical simulation cost analysis (memory, time scaling)
- Quantum algorithm behavior studies (solution quality vs circuit depth, parameter optimization)
- Comparison of algorithmic strategies (e.g., QAOA vs LCU-enhanced QAOA)

### 1.3 In Scope

| Feature | Description |
|---------|-------------|
| State representations | State vector (pure), density matrix lower-triangle (mixed) |
| Evolution picture | Schrödinger (evolve states/density matrices forward) |
| Gates/Channels | Hardcoded optimized gates, arbitrary Kraus channels, LCU |
| Parameterized circuits | Product circuits (QAOA/VQE) and sum circuits (LCU) via function pointers |
| Observables | Diagonal Hamiltonians, sparse Pauli sums |
| Measurements | Expectation values, gradients (adjoint method), shot-based sampling |
| Noise | Runtime-toggled noise on gates |
| Parallelism | Shared-memory (OpenMP, optional), cache optimization, SIMD |
| Utilities | State copy, save/load to disk |

### 1.4 Out of Scope

- Stabilizer/Clifford simulation
- Visualization / plotting
- Predefined error-correction circuits (custom LCU protocols supported)

### 1.5 Scalability Targets

| Representation | Max qubits | Memory required |
|----------------|------------|-----------------|
| Pure state | n = 32 | 64 GB RAM |
| Mixed state | n = 16 | 32 GB (packed storage) |

**Performance targets:**
- Runtime: <48h for typical optimization campaigns
- Test suite: <2 minutes for all gate tests

**Design constraint:** RAM-limited rather than compute-limited. Packed storage and direct sparse operations avoid constructing full gate matrices via Kronecker products.

### 1.6 Implementation Order

1. Pure state simulation (noise-free)
2. Mixed state simulation (noise-free)
3. Full QAOA circuit support
4. Noise channels (runtime toggle)

### 1.7 Outlook (Future Extensions)

**Performance optimizations:**
- GPU acceleration
- Distributed memory (MPI)
- Split complex representation for enhanced SIMD

**Heisenberg picture (observable evolution):**

Evolve observables backward through circuits: O → U†OU. This enables efficient computation of moment matrices for advanced protocols like nested LCU parameter optimization.

*Use case:* For layered LCU circuits (LCU_1 → ... → LCU_p), updating parameters of layer k requires:
```
M_ij = tr[ U_i · ρ_k · U_j† · C_k ]

where:
  ρ_k = LCU_{k-1}(...LCU_1(ρ))      [forward evolved state]
  C_k = LCU_p†(...LCU_{k+1}†(C))    [backward evolved observable]
```

This dual evolution approach is O(p·m + m²) vs O(p·m³) for naive branching, providing 10-1000× speedup for typical m (number of LCU terms).

*Required additions:*
- `hermitian_t` type with lower-triangle storage for evolved observables
- `apply_gate_heisenberg(hermitian_t* obs, ...)` for backward evolution
- `trace_product(hermitian_t* A, hermitian_t* B)` for moment matrix entries

---

## 2. Architecture

### 2.1 Module Structure

```
┌─────────────────────────────────────────────────────────────────┐
│                        Public API                                │
├─────────────────────────────────────────────────────────────────┤
│  Module A          Module B           Module C        Module D   │
│  ──────────        ──────────         ──────────      ────────── │
│  State             Gates/Channels     Measurement     Utilities  │
│  - Allocation      - Unitary gates    - Expectation   - Memory   │
│  - Initialization  - Parameterized    - Gradients     - I/O      │
│  - Copy            - Noise channels   - Sampling      - Parallel │
│  - Serialization   - LCU operations   - Pauli sums    - RNG      │
└─────────────────────────────────────────────────────────────────┘
```

### 2.2 Design Principles

1. **Modularity**: Each module is independently testable. Gates receive states and operate on them without knowledge of initialization. Measurements receive states without knowledge of circuit history.

2. **Type dispatch**: A quantum state can be pure or mixed. Gates dispatch to the appropriate implementation (`gate_pure()` or `gate_mixed()`) based on state type.

3. **Separation of concerns**: States evolve through circuits (Schrödinger picture). Observables are used only at measurement time and are not evolved. (Heisenberg picture is a future extension — see Outlook.)

4. **Simulation type fixed at init**: Pure or mixed is chosen at allocation time and remains constant.

---

## 3. Data Structures

### 3.1 Quantum State

See `docs/STATE_MODULE.md` Section 2 for complete type definitions and invariants.

```c
typedef enum state_type {
    PURE,           // State vector |ψ⟩
    MIXED_PACKED,   // Density matrix ρ, Hermitian packed lower triangle
    MIXED_TILED     // Density matrix ρ, blocked layout for cache optimization
} state_type_t;

typedef struct state {
    state_type_t type;   // Storage format (set before init)
    cplx_t *data;        // Heap-allocated array
    qubit_t qubits;      // Number of qubits n
} state_t;
```

### 3.2 Memory Layout

See `docs/STATE_MODULE.md` Section 4 for index formulas and detailed specifications.

**PURE:** Contiguous array of size N = 2ⁿ with little-endian qubit convention.

**MIXED_PACKED:** Lower triangle, column-major (LAPACK `uplo='L'`). Size: N(N+1)/2.

```
Matrix (N=4):                  Packed storage:
┌──────────────────────┐
│ ρ00                  │       [0]: ρ00  ┐
│ ρ10  ρ11             │       [1]: ρ10  │ Col 0
│ ρ20  ρ21  ρ22        │       [2]: ρ20  │
│ ρ30  ρ31  ρ32  ρ33   │       [3]: ρ30  ┘
└──────────────────────┘       [4]: ρ11  ┐ Col 1
                               ...
```

**MIXED_TILED:** Lower-triangular tiles, each TILE_DIM×TILE_DIM (default 32×32 = 16KB). Tiles stored row-major at tile level, elements row-major within tiles. Enables L1/L2 cache locality for gate operations.

```
Tile grid (N=64, TILE_DIM=32):        Linear storage:
┌─────────┬─────────┐
│  T(0,0) │         │                 [0..1023]:    T(0,0) elements
├─────────┼─────────┤                 [1024..2047]: T(1,0) elements
│  T(1,0) │  T(1,1) │                 [2048..3071]: T(1,1) elements
└─────────┴─────────┘
(store lower triangle of tiles)
```

### 3.3 Observables

**Diagonal Hamiltonian** (for combinatorial optimization):
```c
// Array of 2^n real values: H_diag[i] = ⟨i|H|i⟩
double* H_diag;
```

**Sparse Pauli Sum** (for general observables):
```c
typedef enum { PAULI_I = 0, PAULI_X = 1, PAULI_Y = 2, PAULI_Z = 3 } pauli_t;

typedef struct {
    double coeff;           // Real coefficient (observables are Hermitian)
    pauli_t* ops;           // Array of n_qubits Pauli operators
} pauli_term_t;

typedef struct {
    int n_terms;
    int n_qubits;
    pauli_term_t* terms;
} pauli_sum_t;
```

### 3.4 Parameterized Circuits

```c
// Function pointer type for parameterized circuits
// Works for both product circuits (QAOA/VQE) and sum circuits (LCU)
typedef void (*circuit_fn)(
    state_t* s,
    const double* params,
    int n_params,
    void* ctx              // User context (e.g., graph structure)
);
```

---

## 4. Public API

### 4.1 State Management

See `docs/STATE_MODULE.md` Section 5 for complete API reference.

```c
// Initialization (set state->type before calling)
void state_init(state_t *state, qubit_t qubits, cplx_t **data);  // Zero-init or transfer ownership
void state_plus(state_t *state, qubit_t qubits);                 // |+⟩^⊗n or |+⟩⟨+|^⊗n
void state_free(state_t *state);

// Copy
state_t state_cp(const state_t *state);  // Deep copy

// Element access (handles Hermitian symmetry for mixed types)
cplx_t state_get(const state_t *state, dim_t row, dim_t col);
void   state_set(state_t *state, dim_t row, dim_t col, cplx_t val);
```

### 4.2 Gates

```c
// Single-qubit gates
qs_error_t apply_X(state_t* s, int q);
qs_error_t apply_Y(state_t* s, int q);
qs_error_t apply_Z(state_t* s, int q);
qs_error_t apply_H(state_t* s, int q);
qs_error_t apply_S(state_t* s, int q);
qs_error_t apply_T(state_t* s, int q);

// Parameterized single-qubit gates
qs_error_t apply_Rx(state_t* s, int q, double theta);
qs_error_t apply_Ry(state_t* s, int q, double theta);
qs_error_t apply_Rz(state_t* s, int q, double theta);

// Two-qubit gates
qs_error_t apply_CX(state_t* s, int ctrl, int tgt);   // CNOT
qs_error_t apply_CY(state_t* s, int ctrl, int tgt);
qs_error_t apply_CZ(state_t* s, int ctrl, int tgt);
qs_error_t apply_SWAP(state_t* s, int q1, int q2);

// Parameterized two-qubit gates
qs_error_t apply_CRz(state_t* s, int ctrl, int tgt, double theta);
qs_error_t apply_RZZ(state_t* s, int q1, int q2, double theta);  // For QAOA

// General channel (Kraus operators)
qs_error_t apply_channel(state_t* s, int n_kraus,
                         const double complex** kraus_ops, int n_qubits_affected,
                         const int* qubits);
```

### 4.3 Runtime Configuration

```c
// Noise toggle
qs_error_t set_noise_enabled(int enabled);     // 0 = off, 1 = on
qs_error_t get_noise_enabled(int* enabled);
qs_error_t set_noise_param(double p);          // Noise strength parameter

// Thread safety mode
qs_error_t set_thread_safe_mode(int enabled);  // 0 = single-threaded (faster), 1 = thread-safe
qs_error_t get_thread_safe_mode(int* enabled);
```

### 4.4 Measurement

```c
// Expectation values (result via pointer)
qs_error_t expect_diag(const state_t* s, const double* H_diag, double* result);
qs_error_t expect_pauli(const state_t* s, const pauli_sum_t* H, double* result);

// Shot-based measurement
qs_error_t expect_diag_shots(const state_t* s, const double* H_diag,
                             int n_shots, double* result);
qs_error_t expect_pauli_shots(const state_t* s, const pauli_sum_t* H,
                              int n_shots, double* result);

// Gradients (method auto-selected by library)
qs_error_t grad_expect_diag(
    circuit_fn circuit,
    const double* params,
    int n_params,
    void* ctx,
    state_t* initial_state,
    const double* H_diag,
    double* grad_out           // Output: array of n_params gradients
);

qs_error_t grad_expect_pauli(
    circuit_fn circuit,
    const double* params,
    int n_params,
    void* ctx,
    state_t* initial_state,
    const pauli_sum_t* H,
    double* grad_out           // Output: array of n_params gradients
);
```

### 4.5 LCU (Parameterized Sum Circuits)

LCU circuits implement linear combinations of unitaries: Σᵢ cᵢ Uᵢ. Users build LCU circuits using state arithmetic primitives:

```c
// State arithmetic for LCU construction
qs_error_t state_scale(state_t* s, double complex c);           // s → c·s
qs_error_t state_add(state_t* dest, const state_t* src);        // dest → dest + src
qs_error_t state_add_scaled(state_t* dest, const state_t* src,
                            double complex c);                   // dest → dest + c·src

// Workspace management for LCU
state_t* state_alloc_like(const state_t* s);                    // Allocate with same type/size
```

LCU circuits use the same `circuit_fn` interface as product circuits. Example pattern:
```c
void my_lcu_circuit(state_t* s, const double* params, int n_params, void* ctx) {
    state_t* temp = state_alloc_like(s);
    state_t* accum = state_alloc_like(s);
    state_init_zero(accum);  // Zero accumulator

    for (int i = 0; i < n_terms; i++) {
        state_copy_into(temp, s);           // temp = input state
        apply_unitary_i(temp, params, ctx); // temp = Uᵢ|ψ⟩
        state_add_scaled(accum, temp, coeffs[i]); // accum += cᵢ·Uᵢ|ψ⟩
    }
    state_copy_into(s, accum);  // s = Σᵢ cᵢ Uᵢ|ψ⟩
    state_free(temp);
    state_free(accum);
}
```

---

## 5. Internal Design

### 5.1 Gate Implementation Pattern

```c
// Public dispatcher
qs_error_t apply_X(state_t* s, int q) {
    if (s == NULL) return QS_ERR_NULL;
    if (q < 0 || q >= s->n_qubits) return QS_ERR_QUBIT;

    qs_error_t err;
    if (s->type == STATE_PURE) {
        err = apply_X_pure(s, q);
    } else {
        err = apply_X_mixed(s, q);
    }
    if (err != QS_OK) return err;

    int noise_on;
    err = get_noise_enabled(&noise_on);
    if (err != QS_OK) return err;
    if (noise_on && s->type == STATE_MIXED) {
        err = apply_noise_1q(s, q);
    }
    return err;
}

// Internal implementations (static)
static qs_error_t apply_X_pure(state_t* s, int q);
static qs_error_t apply_X_mixed(state_t* s, int q);
```

**Noise and pure states:** Noise channels (Kraus operators) produce mixed states and are only applied when `s->type == STATE_MIXED`. For pure state simulations, enable noise by allocating a mixed state and initializing it as a pure density matrix (e.g., via `state_init_zero()`).

### 5.2 Stride-Based Amplitude Access

For qubit q, amplitudes pair with stride 2^q:
```c
size_t stride = (size_t)1 << q;
size_t dim = (size_t)1 << n_qubits;
for (size_t i = 0; i < dim; i += 2 * stride) {
    for (size_t j = 0; j < stride; j++) {
        size_t idx0 = i + j;           // Qubit q = 0
        size_t idx1 = i + j + stride;  // Qubit q = 1
        // Apply 2x2 gate matrix to (amp[idx0], amp[idx1])
    }
}
```

### 5.3 Gradient Computation

Two methods supported, auto-selected based on memory and parameter count.

**Notation:** Let *d* denote circuit depth (total number of gate layers), *p* denote the number of parameters, and *n* denote the number of qubits.

**Method 1: Adjoint (backpropagation)**

Best when memory is available and many parameters exist.

| Variant | Memory | Time |
|---------|--------|------|
| Full caching | O(d × 2ⁿ) | O(d × 2ⁿ) |
| Checkpointing | O(√d × 2ⁿ) | O(d × 2ⁿ) |

*Full caching:*
1. Forward pass: apply gates, cache all intermediate states
2. Backward pass: propagate gradients through cached states

*Checkpointing:*
1. Forward pass: apply gates, cache only √d checkpoints
2. Backward pass: recompute intermediates between checkpoints

**Method 2: Parameter-shift rule**

Best when memory is severely constrained (large n).

| Memory | Time |
|--------|------|
| O(2ⁿ) | O(p × d × 2ⁿ) |

For each parameter θ in gate e^{-iθG/2}:
```
∂⟨H⟩/∂θ = ½(⟨H⟩|_{θ+π/2} - ⟨H⟩|_{θ-π/2})
```

Requires 2p circuit evaluations but only one state in memory at a time.

**Auto-selection criteria:**
- Available memory vs √d × 2ⁿ requirement
- Number of parameters p vs √d tradeoff
- Library always selects optimal method automatically

### 5.4 Parallelization Strategy

**Branch-level vs simulation-level parallelism:**
- n ≤ ~14: Prefer branch-level (run different states in parallel)
- n > ~14: Prefer simulation-level (parallelize 2ⁿ operations per gate)
- Auto-selected based on n and detected core count

**Within-gate parallelism (OpenMP):**

OpenMP parallelization is enabled by default (`-DENABLE_OPENMP=ON`) and can be disabled at build time.
Gates use conditional parallelization with workload-appropriate thresholds:
- **Pure states**: `OMP_THRESHOLD = 2048` (n ≥ 11 qubits) — O(dim) work per gate
- **Mixed states**: `OMP_THRESHOLD = 64` (n ≥ 6 qubits) — O(dim²) work per gate

```c
// Pure states: parallelize over independent blocks
#pragma omp parallel for if(dim >= OMP_THRESHOLD)
for (dim_t i = 0; i < dim; i += step) {
    cblas_zswap(stride, state->data + i, 1, state->data + i + stride, 1);
}

// Leftmost qubit edge case: parallelize element-wise when outer loop has single iteration
if (step >= dim && stride >= OMP_THRESHOLD) {
    #pragma omp parallel for
    for (dim_t j = 0; j < stride; ++j) {
        // Process element pairs directly
    }
}
```

**Mixed state parallelization:**

The `TRAVERSE_MIXED_1Q` macro parallelizes the outer `col_block` loop when multiple independent
column blocks exist. For the Hadamard gate, which processes independent 2×2 blocks, dynamic
scheduling balances load across threads:

```c
#pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
for (dim_t c = 0; c < dim; ++c) { /* process 2×2 blocks */ }
```

**Note:** When targeting the leftmost qubit (target = n-1), the outer loop has only one iteration.
Pure state gates handle this by switching to element-wise parallelization. Mixed state gates
using `TRAVERSE_MIXED_1Q` do not parallelize in this case due to cross-column dependencies in
the packed storage format.

### 5.5 Thread Safety

Controlled via `set_thread_safe_mode()` (see Section 4.3):
- **Single-threaded mode** (enabled=0): Global RNG, shared workspace (faster)
- **Thread-safe mode** (enabled=1): Thread-local RNG, no shared mutable state

---

## 6. Serialization Format

### 6.1 Binary Format (Primary)

```
Header (32 bytes):
  [0-3]   Magic: "QSIM"
  [4-7]   Version: uint32
  [8-11]  n_qubits: uint32
  [12-15] type: uint32 (0=pure, 1=mixed)
  [16-19] endianness: uint32 (0x01020304 for detection)
  [20-31] Reserved

Body:
  Amplitudes as double complex array
  - Pure:  2^n entries
  - Mixed: 2^n * (2^n + 1) / 2 entries (lower triangle)
```

### 6.2 Text Format (Debug)

Optional human-readable export for debugging small states.

---

## 7. Error Handling

### 7.1 State Module

See `docs/STATE_MODULE.md` Section 6 for complete specification.

| Error Type | Handling | Rationale |
|------------|----------|-----------|
| Programmer error | `assert()` | Bugs should crash immediately for easy debugging |
| Out-of-memory | Set `data = NULL`, print to stderr | Runtime condition, caller can handle |

**Caller responsibility:** Check `state.data == NULL` after `state_init()` or `state_plus()`.

### 7.2 Other Modules

Functions return `qs_error_t` status codes. Caller is responsible for checking return codes.

```c
typedef enum {
    QS_OK           =  0,   // Success
    QS_ERR_NULL     = -1,   // Null pointer argument
    QS_ERR_OOM      = -2,   // Out of memory
    QS_ERR_QUBIT    = -3,   // Invalid qubit index
    QS_ERR_TYPE     = -4,   // Invalid state type for operation
    QS_ERR_FILE     = -5,   // File I/O error
    QS_ERR_FORMAT   = -6,   // Invalid file format
    QS_ERR_PARAM    = -7,   // Invalid parameter value
} qs_error_t;
```

---

## 8. Dependencies

| Dependency | Purpose | Required |
|------------|---------|----------|
| CBLAS | Linear algebra primitives | Yes |
| LAPACK | Matrix operations (for unit tests) | Yes |
| complex.h | C99 complex numbers | Yes |
| OpenMP | Shared-memory parallelism | Optional (enabled by default) |
| Unity | Unit testing framework | Testing only |
| CMake | Build system | Yes |

**CMake options:**
```bash
cmake -DENABLE_OPENMP=ON ..   # Enable OpenMP parallelization (default)
cmake -DENABLE_OPENMP=OFF ..  # Disable OpenMP (single-threaded)
```

---

## 9. Build System

### 9.1 CMake Structure

```
qlib/
├── CMakeLists.txt
├── src/
│   ├── state.c
│   ├── gate.c
│   ├── measure.c
│   └── util.c
├── include/
│   └── qlib.h          # Public API header
├── tests/
│   ├── test_state.c
│   ├── test_gate.c
│   └── test_measure.c
└── examples/
    └── qaoa_maxcut.c
```

### 9.2 Build Outputs

- `libqlib.a` — Static library
- `libqlib.so` / `libqlib.dylib` — Shared library (optional)
- Test executables

---

## 10. Testing Strategy

### 10.1 Unit Tests (Priority 1)

| Module | Tests |
|--------|-------|
| State | Allocation, initialization, index calculations, copy |
| Gates | Each gate against known matrix multiplication |
| Measurement | Expectation values against analytical results |

### 10.2 Integration Tests (Priority 2)

- Small circuits with known outputs
- Pure vs mixed giving same ⟨O⟩ for pure states lifted to density matrices
- LCU circuits producing correct weighted superpositions

### 10.3 Gradient Tests (Priority 3)

- Compare adjoint-method gradients to finite differences
- Tolerance: relative error < 1e-6

### 10.4 Numerical Tolerance

- Standard: 1e-12 for double precision
- Runtime checks: optional norm/trace verification after operations

---

## 11. Appendix: Qubit Convention

**Little-endian convention:**

For n qubits, basis state |bₙ₋₁ bₙ₋₂ ... b₁ b₀⟩ has array index:

```
index = Σᵢ bᵢ × 2ⁱ = b₀ + 2b₁ + 4b₂ + ... + 2^(n-1) × bₙ₋₁
```

Qubit 0 is the **rightmost** (least significant) bit.

**Example (3 qubits):**
```
|000⟩ → index 0
|001⟩ → index 1    (qubit 0 = 1)
|010⟩ → index 2    (qubit 1 = 1)
|011⟩ → index 3
|100⟩ → index 4    (qubit 2 = 1)
|101⟩ → index 5
|110⟩ → index 6
|111⟩ → index 7
```

**Gate stride:** A gate on qubit q accesses pairs with stride 2^q.

---

## 12. Glossary

| Term | Definition |
|------|------------|
| Adjoint method | Gradient computation via backward propagation through circuit |
| Heisenberg picture | Observables evolve (O → U†OU), states fixed. *Future extension.* |
| Kraus operators | Matrices {Kᵢ} defining a quantum channel: ρ → Σᵢ Kᵢ ρ Kᵢ† |
| LCU | Linear Combination of Unitaries: Σᵢ cᵢ Uᵢ |
| Pauli string | Tensor product of Pauli matrices: σᵢ₁ ⊗ σᵢ₂ ⊗ ... ⊗ σᵢₙ |
| Product circuit | Sequential application of parameterized unitaries |
| Schrödinger picture | States evolve (|ψ⟩ → U|ψ⟩), observables fixed. *Current implementation.* |
| Shot noise | Statistical fluctuations from finite measurement samples |
| Sum circuit | Weighted sum of unitaries (LCU) |

---

*Document generated via Socratic design process.*
*Last updated: 2026-01-23 (OpenMP parallelization for gate operations)*
