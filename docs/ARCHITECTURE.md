# qSim Architecture Overview

**Document Version**: 1.0
**Last Updated**: 2025-12-22
**Status**: Draft
**Related**: [PRD.md](../PRD.md), [TECHNICAL_SPEC.md](./TECHNICAL_SPEC.md), [API.md](./API.md)

---

## Table of Contents

1. [System Overview](#1-system-overview)
2. [Design Principles](#2-design-principles)
3. [Layer Architecture](#3-layer-architecture)
4. [Module Design](#4-module-design)
5. [Data Flow](#5-data-flow)
6. [Concurrency Model](#6-concurrency-model)
7. [Platform Abstraction](#7-platform-abstraction)
8. [Build System Architecture](#8-build-system-architecture)
9. [Testing Architecture](#9-testing-architecture)
10. [Extension Points](#10-extension-points)

---

## 1. System Overview

### 1.1 Purpose

qSim is a high-performance quantum computing simulation library designed for:
- Resource-constrained environments
- C/C++ application integration
- Workstation-scale quantum simulations (up to 30 qubits)
- High memory efficiency and execution speed

### 1.2 Key Characteristics

- **Language**: Pure C17 implementation
- **Dependencies**: Minimal (BLAS/LAPACK only)
- **Paradigm**: Imperative, procedural
- **Memory Management**: Explicit (manual allocation/deallocation)
- **Parallelism**: OpenMP (optional, controlled by user)
- **Portability**: Cross-platform (macOS, Linux, x86, ARM)

### 1.3 High-Level Architecture

```
┌────────────────────────────────────────────────┐
│         User Application (C/C++)               │
└────────────────┬───────────────────────────────┘
                 │
      ┌──────────▼──────────┐
      │   Public API        │  qlib.h
      │   (qlib.h)          │
      └──────────┬──────────┘
                 │
      ┌──────────▼──────────┐
      │  Circuit Layer      │  circuit.h (v1.0)
      │  (circuit.h)        │
      └──────────┬──────────┘
                 │
      ┌──────────▼──────────┐
      │  Gate Operations    │  gates.h (v1.0)
      │  (gates.h)          │
      └──────────┬──────────┘
                 │
      ┌──────────▼──────────┐
      │  State Management   │  state.h
      │  (state.h)          │
      └──────────┬──────────┘
                 │
      ┌──────────▼──────────┐
      │  Core Types         │  q_types.h
      │  (q_types.h)        │
      └──────────┬──────────┘
                 │
      ┌──────────▼──────────┐
      │  BLAS/LAPACK        │  Platform-specific
      │  (Accelerate/       │
      │   OpenBLAS)         │
      └─────────────────────┘
```

---

## 2. Design Principles

### 2.1 Core Principles

1. **Simplicity**: Prefer simple, readable code over clever optimizations
2. **Explicitness**: No hidden state, clear function contracts
3. **Performance**: Optimize hot paths, profile-driven optimization
4. **Portability**: Avoid platform-specific code outside abstraction layers
5. **Reliability**: Comprehensive testing, defensive programming
6. **Modularity**: Clear separation of concerns, minimal coupling

### 2.2 Code Organization Philosophy

```
Guideline: Each module has a single, well-defined responsibility

- state.h/c:   State representation and lifecycle management
- gates.h/c:   Gate definitions and application logic
- circuit.h/c: Circuit construction and execution
- measure.h/c: Measurement operations and expectation values
```

### 2.3 API Design Philosophy

**Principles**:
- **Orthogonality**: Functions do one thing well
- **Consistency**: Similar operations have similar signatures
- **Safety**: Validate inputs, return error codes
- **Clarity**: Self-documenting names, clear parameter order

**Example**:
```c
// Good: Clear, consistent naming
state_t* state_init(qubit_t num_qubits, state_type_t type);
void state_free(state_t* state);
state_t* state_cp(const state_t* src);

// Bad: Inconsistent naming
state_t* create_state(qubit_t n, state_type_t t);
void destroy(state_t* s);
state_t* copy_state(const state_t* src);
```

---

## 3. Layer Architecture

### 3.1 Layer Responsibilities

#### Layer 1: Core Types (q_types.h)
**Responsibility**: Fundamental type definitions

**Contains**:
- `cplx_t`: Complex number type
- `qubit_t`: Qubit identifier
- `dim_t`: Dimension type
- Mathematical macros and constants

**Dependencies**: Standard C library only

---

#### Layer 2: State Management (state.h)
**Responsibility**: Quantum state lifecycle and representation

**Contains**:
- `state_t`: State structure
- `state_init()`, `state_free()`: Lifecycle
- `state_cp()`: Deep copy
- `state_plus()`: Special state initialization

**Dependencies**: Core Types, BLAS (for allocation alignment)

---

#### Layer 3: Gate Operations (gates.h)
**Responsibility**: Quantum gate application to states

**Contains**:
- Single-qubit gates (X, Y, Z, H, S, T, Rx, Ry, Rz)
- Two-qubit gates (CNOT, SWAP)
- Multi-qubit gates (Toffoli)
- Controlled gate variants

**Dependencies**: State Management, BLAS (for matrix operations)

---

#### Layer 4: Circuit Management (circuit.h)
**Responsibility**: Circuit construction and execution

**Contains**:
- `circuit_t`: Circuit structure
- `circuit_create()`, `circuit_free()`: Lifecycle
- `circuit_add_gate_*()`: Gate addition
- `circuit_execute()`: Execution engine

**Dependencies**: Gate Operations, State Management

---

#### Layer 5: Measurement (measure.h)
**Responsibility**: Measurement and observable calculations

**Contains**:
- `measure_all()`, `measure_qubit()`: Basis measurement
- `expectation_pauli()`: Pauli expectation values
- `expectation_hermitian()`: General observables
- Observable moment matrices

**Dependencies**: State Management, BLAS

---

#### Layer 6: Public API (qlib.h)
**Responsibility**: Unified API header

**Contains**: Includes all lower layers

**Purpose**: Single include for users

---

### 3.2 Dependency Graph

```
    qlib.h
      │
      ├─── measure.h ──┐
      │                │
      ├─── circuit.h ──┼─→ gates.h ──┐
      │                │              │
      │                └──────────────┼─→ state.h ──→ q_types.h
      │                               │
      └───────────────────────────────┘

All layers depend on:
  - BLAS/LAPACK (platform-specific)
  - OpenMP (optional)
```

---

## 4. Module Design

### 4.1 State Module (state.h/c)

**Purpose**: Manage quantum state representation and lifecycle

**Key Structures**:
```c
typedef struct {
    state_type_t type;    // PURE or MIXED
    cplx_t* data;         // State data (aligned)
    qubit_t num_qubits;   // Number of qubits
    dim_t dim;            // Hilbert space dimension
} state_t;
```

**Internal Functions** (not exposed in API):
```c
// Memory allocation with alignment
static cplx_t* allocate_aligned(size_t num_elements);

// Validate state invariants
static bool validate_state(const state_t* state);
```

**Public API**:
- Initialization: `state_init()`, `state_init_basis()`, `state_plus()`
- Lifecycle: `state_free()`
- Operations: `state_cp()`
- Inspection: `state_len()`, `state_get_amplitude()`

---

### 4.2 Gates Module (gates.h/c)

**Purpose**: Implement quantum gate operations

**Design Pattern**: Strategy pattern for gate application

**Internal Architecture**:
```c
// Gate application engine
static void apply_single_qubit_gate(state_t* state, qubit_t target,
                                    const cplx_t U[4]);

static void apply_two_qubit_gate(state_t* state, qubit_t q1, qubit_t q2,
                                 const cplx_t U[16]);
```

**Gate Definitions**:
```c
// Pauli-X matrix (static constant)
static const cplx_t GATE_X[4] = {
    {0.0, 0.0}, {1.0, 0.0},
    {1.0, 0.0}, {0.0, 0.0}
};
```

**Optimization Strategy**:
- Single-qubit: Custom vectorized loops (faster than BLAS for this pattern)
- Two-qubit: Custom loops with cache-friendly access
- Multi-qubit: BLAS matrix-vector operations (only if beneficial)

---

### 4.3 Circuit Module (circuit.h/c)

**Purpose**: Circuit construction and execution

**Design Decision**: Immediate vs. Deferred Execution

**Option A: Immediate Execution**
```c
circuit_t* circuit = circuit_create(3);
circuit_execute_gate_h(circuit, state, 0);  // Apply immediately
```

**Option B: Deferred Execution** (Recommended for v1.0)
```c
circuit_t* circuit = circuit_create(3);
circuit_add_gate_h(circuit, 0);    // Add to circuit
circuit_add_gate_cnot(circuit, 0, 1);
circuit_execute(circuit, state);    // Execute all at once
```

**Internal Representation** (for Option B):
```c
typedef struct {
    gate_type_t type;  // GATE_H, GATE_CNOT, etc.
    qubit_t targets[3]; // Target qubits (max 3 for Toffoli)
    double param;       // For rotation gates
} gate_op_t;

typedef struct {
    qubit_t num_qubits;
    gate_op_t* gates;   // Dynamic array of gate operations
    size_t num_gates;
    size_t capacity;
} circuit_t;
```

---

### 4.4 Measurement Module (measure.h/c)

**Purpose**: Quantum measurement and expectation values

**Key Algorithms**:

**Computational Basis Measurement**:
```c
uint64_t measure_all(state_t* state) {
    // 1. Compute probabilities: P(i) = |ψ_i|^2
    double* probs = compute_probabilities(state);

    // 2. Sample from distribution
    uint64_t outcome = sample_categorical(probs, state->dim);

    // 3. Collapse state
    collapse_to_basis(state, outcome);

    return outcome;
}
```

**Expectation Value**:
```c
double expectation_hermitian(const state_t* state, const cplx_t* H) {
    // ⟨ψ|H|ψ⟩
    // 1. Compute H|ψ⟩ using BLAS
    cplx_t* H_psi = apply_matrix(H, state->data, state->dim);

    // 2. Compute ⟨ψ|H|ψ⟩ = ψ† · (H|ψ⟩)
    cplx_t result = dot_product(state->data, H_psi, state->dim);

    return result.real;  // Expectation is real for Hermitian operators
}
```

---

## 5. Data Flow

### 5.1 Typical Simulation Flow

```
┌─────────────────┐
│ 1. Create State │  state = state_init(n, STATE_PURE)
└────────┬────────┘
         │
┌────────▼────────┐
│ 2. Create       │  circuit = circuit_create(n)
│    Circuit      │
└────────┬────────┘
         │
┌────────▼────────┐
│ 3. Add Gates    │  circuit_add_gate_h(circuit, 0)
│                 │  circuit_add_gate_cnot(circuit, 0, 1)
└────────┬────────┘
         │
┌────────▼────────┐
│ 4. Execute      │  circuit_execute(circuit, state)
│    Circuit      │
└────────┬────────┘
         │
┌────────▼────────┐
│ 5. Measure or   │  outcome = measure_all(state)
│    Compute      │  OR
│    Expectation  │  exp_val = expectation_pauli(state, "ZZ")
└────────┬────────┘
         │
┌────────▼────────┐
│ 6. Cleanup      │  state_free(state)
│                 │  circuit_free(circuit)
└─────────────────┘
```

### 5.2 Memory Flow

```
Allocation:
  User → state_init() → allocate_aligned() → [Aligned Memory]
                                                     │
                                                     ▼
                                              [state_t structure]
Usage:
  [state_t] → gates/circuit operations → [Modified state data]

Deallocation:
  User → state_free() → free(data) → free(state) → NULL
```

---

## 6. Concurrency Model

### 6.1 Thread Safety

**Policy**: Library is **NOT thread-safe by default** in v1.0

**Rationale**:
- Simplifies implementation
- Avoids synchronization overhead
- Users can manage concurrency at application level

**Documentation Requirement**: Clearly document thread-safety guarantees

### 6.2 OpenMP Parallelization

**Usage**: Internal parallelization within gate operations

**Pattern**:
```c
void apply_single_qubit_gate(state_t* state, qubit_t target,
                             const cplx_t U[4]) {
    #pragma omp parallel for if(state->num_qubits >= 20)
    for (uint64_t i = 0; i < state->dim; i += 2) {
        // Apply gate to 2D subspace
    }
}
```

**Threshold**: Only parallelize for large state vectors (n ≥ 20 qubits)

**Environment Variable**: User controls thread count via `OMP_NUM_THREADS`

### 6.3 Future Considerations (Post-v1.0)

- **Thread-safe API option**: Mutex-protected operations
- **Concurrent circuit execution**: Execute multiple circuits in parallel
- **GPU offloading**: Transparent GPU acceleration for large simulations

---

## 7. Platform Abstraction

### 7.1 Platform-Specific Code Organization

```
qsim/
├── src/
│   ├── state.c           # Platform-independent
│   ├── gates.c           # Platform-independent
│   └── platform/         # Platform-specific code
│       ├── blas_wrapper.h    # BLAS abstraction
│       ├── blas_macos.c      # macOS Accelerate
│       └── blas_linux.c      # Linux OpenBLAS
```

### 7.2 BLAS Abstraction Layer

**Purpose**: Provide unified interface to platform-specific BLAS

**Header** (blas_wrapper.h):
```c
// Abstract BLAS interface
typedef int64_t blasint;  // ILP64 integer type

// Wrapper functions
void qsim_zdotc(blasint n, const cplx_t* x, blasint incx,
                const cplx_t* y, blasint incy, cplx_t* result);

void qsim_zgemv(/* ... */);  // Matrix-vector product
```

**Implementation**:
```c
// blas_macos.c
#include <Accelerate/Accelerate.h>

void qsim_zdotc(blasint n, const cplx_t* x, blasint incx,
                const cplx_t* y, blasint incy, cplx_t* result) {
    cblas_zdotc_sub(n, x, incx, y, incy, result);
}
```

### 7.3 Build-Time Platform Detection

**CMake**:
```cmake
if(APPLE)
    set(BLAS_IMPL "Accelerate")
    target_link_libraries(qsim PRIVATE "-framework Accelerate")
elseif(UNIX)
    set(BLAS_IMPL "OpenBLAS")
    # Build or find OpenBLAS
endif()
```

---

## 8. Build System Architecture

### 8.1 CMake Structure

```
qsim/
├── CMakeLists.txt              # Root configuration
├── cmake/
│   ├── PlatformConfig.cmake    # Platform detection
│   ├── CompilerFlags.cmake     # Compiler options
│   └── Dependencies.cmake      # External dependencies
├── src/
│   └── CMakeLists.txt          # Library build rules
├── test/
│   └── CMakeLists.txt          # Test build rules
└── extern/
    └── Unity/                  # Testing framework
```

### 8.2 Build Targets

**Primary Targets**:
- `qsim`: Main library (static or shared)
- `qsim_tests`: Test executable
- `install`: Install headers and library

**Custom Targets**:
- `verify_ilp64`: BLAS verification tool
- `benchmark`: Performance benchmarks (future)

### 8.3 Dependency Management

**Unity** (Testing):
- Git submodule
- Automatically initialized via CMake

**OpenBLAS** (Linux only):
- ExternalProject_Add in CMake
- Downloaded, configured, and built automatically
- Cached for subsequent builds

---

## 9. Testing Architecture

### 9.1 Test Organization

```
test/
├── test_state.c        # State management tests
├── test_gates.c        # Gate operation tests (v1.0)
├── test_circuit.c      # Circuit tests (v1.0)
├── test_measure.c      # Measurement tests (v1.0)
├── test_qhipster.c     # Comparative tests vs. QHipster
├── test.h              # Test utilities and helpers
└── CMakeLists.txt      # Test build configuration
```

### 9.2 Test Categories

**Unit Tests**:
- Test individual functions in isolation
- Mock dependencies where needed
- Fast execution (< 1 second total)

**Integration Tests**:
- Test complete workflows (circuit creation → execution → measurement)
- Verify module interactions
- Moderate execution time

**Validation Tests**:
- Compare against analytical results
- Cross-validate with reference implementations (Qiskit, Cirq)
- Ensure numerical correctness

**Performance Tests**:
- Benchmark critical operations
- Track performance regressions
- Generate performance reports

### 9.3 Test Utilities (test.h)

```c
// Assertion helpers
#define ASSERT_COMPLEX_EQUAL(expected, actual, tolerance) \
    TEST_ASSERT_DOUBLE_WITHIN(tolerance, (expected).real, (actual).real); \
    TEST_ASSERT_DOUBLE_WITHIN(tolerance, (expected).imag, (actual).imag)

// State validation
bool states_equal(const state_t* s1, const state_t* s2, double tol);

// Test state generators
state_t* create_bell_state(void);
state_t* create_ghz_state(qubit_t n);
```

---

## 10. Extension Points

### 10.1 Future Architecture Enhancements

**Plugin System** (v2.0):
- Allow custom gate implementations
- Support alternative state representations (tensor networks, stabilizers)

**Backend Abstraction** (v2.0):
- CPU backend (current)
- GPU backend (CUDA/OpenCL)
- Distributed backend (MPI)

**Noise Models** (v1.1):
```c
typedef struct {
    noise_model_t* noise_model;
    // ... existing fields
} circuit_t;
```

### 10.2 API Versioning Strategy

**Versioning Scheme**: Semantic versioning (MAJOR.MINOR.PATCH)

**Compatibility Promise**:
- MAJOR: Breaking API changes
- MINOR: Backward-compatible features
- PATCH: Bug fixes only

**Deprecation Policy**:
- Deprecated functions marked with compiler warnings
- Minimum 1 minor version before removal
- Migration guide provided

---

## 11. Design Decisions Log

| Decision | Date | Rationale |
|----------|------|-----------|
| Use C17 instead of C++ | 2025-12-22 | Minimize dependencies, maximize portability |
| Manual memory management | 2025-12-22 | Explicit control, predictable performance |
| Bundle OpenBLAS for Linux | 2025-12-22 | Ensure ILP64 support, simplify installation |
| Deferred circuit execution | TBD | Enable future optimization passes |
| Little-endian qubit ordering | TBD | Match common quantum computing conventions |

---

## 12. References

- [PRD.md](../PRD.md) - Product requirements
- [TECHNICAL_SPEC.md](./TECHNICAL_SPEC.md) - Implementation details
- [API.md](./API.md) - Public API reference

---

**Document Maintenance**: Update as architectural decisions are made during implementation
