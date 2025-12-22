# qSim Technical Specification v1.0

**Document Version**: 1.0
**Last Updated**: 2025-12-22
**Status**: Draft
**Related**: [PRD.md](../PRD.md), [API.md](./API.md), [ARCHITECTURE.md](./ARCHITECTURE.md)

---

## Table of Contents

1. [System Overview](#1-system-overview)
2. [Data Structures](#2-data-structures)
3. [Quantum State Representation](#3-quantum-state-representation)
4. [Gate Implementation](#4-gate-implementation)
5. [Circuit Execution](#5-circuit-execution)
6. [Measurement Operations](#6-measurement-operations)
7. [Memory Management](#7-memory-management)
8. [Performance Optimization](#8-performance-optimization)
9. [Numerical Precision](#9-numerical-precision)
10. [Platform-Specific Details](#10-platform-specific-details)

---

## 1. System Overview

### 1.1 Architecture Layers

```
┌─────────────────────────────────────┐
│      Public API (qlib.h)           │
├─────────────────────────────────────┤
│   Circuit Layer (circuit.h)        │
├─────────────────────────────────────┤
│   Gate Operations (gates.h)        │
├─────────────────────────────────────┤
│   State Management (state.h)       │
├─────────────────────────────────────┤
│   Linear Algebra (BLAS/LAPACK)     │
└─────────────────────────────────────┘
```

### 1.2 Module Dependencies

[To be completed with actual implementation]

---

## 2. Data Structures

### 2.1 Core Types (q_types.h)

#### 2.1.1 Complex Numbers

```c
// Complex number representation
typedef struct {
    double real;
    double imag;
} cplx_t;
```

**Rationale**: Using custom complex type instead of C99 `_Complex` for:
- Explicit memory layout control
- Better BLAS integration
- Portability across compilers

#### 2.1.2 Dimension and Qubit Types

```c
typedef uint32_t qubit_t;      // Qubit identifier (0 to n-1)
typedef uint64_t dim_t;        // Dimension type (2^n for n qubits)
```

**Constraints**:
- Maximum qubits: 30 (dim_t can represent 2^30 = 1,073,741,824)
- Qubit indexing: 0-based, little-endian convention

### 2.2 Quantum State (state.h)

#### 2.2.1 State Type Enumeration

```c
typedef enum {
    STATE_PURE,    // Statevector representation
    STATE_MIXED    // Density matrix representation
} state_type_t;
```

#### 2.2.2 State Structure

```c
typedef struct {
    state_type_t type;    // PURE or MIXED
    cplx_t* data;         // State data (statevector or density matrix)
    qubit_t num_qubits;   // Number of qubits
    dim_t dim;            // Hilbert space dimension (2^num_qubits)
} state_t;
```

**Memory Layout**:
- **Pure state**: Array of `2^n` complex numbers (statevector)
- **Mixed state**: Lower-triangular density matrix (compact storage)

**Invariants**:
- `dim == (1ULL << num_qubits)`
- `data != NULL` for valid states
- For pure states: Sum of |amplitude|^2 = 1.0 (normalized)
- For mixed states: Trace(ρ) = 1.0, ρ = ρ†, ρ ≥ 0

### 2.3 Circuit Structure (circuit.h)

[To be defined in implementation phase]

```c
// Placeholder structure
typedef struct {
    qubit_t num_qubits;
    // Gate list representation TBD
    // Optimization metadata TBD
} circuit_t;
```

---

## 3. Quantum State Representation

### 3.1 Pure States (Statevectors)

#### 3.1.1 Mathematical Representation

For n qubits, the state |ψ⟩ is represented as:

```
|ψ⟩ = Σ(i=0 to 2^n-1) α_i |i⟩
```

Where `α_i` are complex amplitudes stored in `data[]` array.

#### 3.1.2 Memory Layout

```
Index  Binary  State     Amplitude
0      000     |000⟩     data[0]
1      001     |001⟩     data[1]
2      010     |010⟩     data[2]
...
2^n-1  111...  |111...⟩  data[2^n-1]
```

**Storage**: `2^n * sizeof(cplx_t)` = `2^n * 16` bytes

**Example**: 30 qubits requires `2^30 * 16` = 16 GiB

#### 3.1.3 Qubit Ordering Convention

**Little-endian**: Qubit 0 is the least significant bit

```
|q₂q₁q₀⟩ → Computational basis state
```

### 3.2 Mixed States (Density Matrices)

#### 3.2.1 Mathematical Representation

```
ρ = Σ(i,j) ρ_ij |i⟩⟨j|
```

Where ρ is a 2^n × 2^n Hermitian matrix.

#### 3.2.2 Compact Storage (Lower Triangle)

Only store lower triangle due to Hermitian property: ρ_ij = ρ*_ji

**Storage**: `(2^n * (2^n + 1) / 2) * sizeof(cplx_t)` bytes

**Example**: 15 qubits requires ~8.6 GiB (vs. ~17.2 GiB for full matrix)

---

## 4. Gate Implementation

### 4.1 Gate Representation

#### 4.1.1 Single-Qubit Gates

Represented as 2×2 unitary matrices:

```
U = [u00  u01]
    [u10  u11]
```

**Pauli-X Gate**:
```c
cplx_t X[4] = {
    {0.0, 0.0}, {1.0, 0.0},  // [0 1]
    {1.0, 0.0}, {0.0, 0.0}   // [1 0]
};
```

#### 4.1.2 Two-Qubit Gates

Represented as 4×4 unitary matrices (16 complex numbers).

**CNOT Gate** (control=q₁, target=q₀):
```
|00⟩ → |00⟩
|01⟩ → |01⟩
|10⟩ → |11⟩
|11⟩ → |10⟩
```

### 4.2 Gate Application Algorithm

#### 4.2.1 Single-Qubit Gate on Pure State

**Algorithm**: Apply 2×2 matrix to target qubit subspace

```
For each basis state |...b...⟩:
    Extract 2D subspace for target qubit
    Apply matrix multiplication
    Write back results
```

**Complexity**: O(2^n) operations

**BLAS Usage**: Custom vectorized loops (more efficient than BLAS for this pattern)

#### 4.2.2 Two-Qubit Gate on Pure State

**Algorithm**: Apply 4×4 matrix to target qubit pair subspace

**Complexity**: O(2^n) operations

**Optimization**: Cache-friendly access patterns, vectorization

### 4.3 Rotation Gates

#### 4.3.1 Parametric Rotation

```c
// Rx(θ) = exp(-iθX/2)
void apply_rx(state_t* state, qubit_t target, double theta);

// Matrix form:
// [cos(θ/2)      -i*sin(θ/2)]
// [-i*sin(θ/2)    cos(θ/2)  ]
```

Similar for Ry, Rz rotations.

---

## 5. Circuit Execution

### 5.1 Execution Model

[To be defined]

**Options**:
1. **Immediate execution**: Apply each gate immediately to state
2. **Deferred execution**: Build gate list, optimize, then apply
3. **Hybrid**: Optimization passes followed by sequential application

**v1.0 Decision**: [To be determined based on performance testing]

### 5.2 Optimization Opportunities

[Out of scope for v1.0, but architecture should allow future optimization]

- Gate fusion (combine adjacent single-qubit gates)
- Commutation analysis
- Dead code elimination

---

## 6. Measurement Operations

### 6.1 Computational Basis Measurement

#### 6.1.1 Full Measurement

**Algorithm**:
```
1. Compute probabilities: P(i) = |α_i|²
2. Generate random number r ∈ [0,1)
3. Select outcome based on cumulative probability
4. Collapse state to |outcome⟩
```

**Implementation**:
```c
// Pseudocode
uint64_t measure_all(state_t* state) {
    double* probs = compute_probabilities(state);
    uint64_t outcome = sample_distribution(probs, state->dim);
    collapse_to_basis_state(state, outcome);
    return outcome;
}
```

#### 6.1.2 Partial Measurement

Measure subset of qubits, preserve others.

**Complexity**: More involved, requires partial trace operation.

### 6.2 Expectation Values

#### 6.2.1 Pauli Observable

For observable O and state |ψ⟩:

```
⟨O⟩ = ⟨ψ|O|ψ⟩
```

**Implementation** (pure state):
```c
double expectation_pauli(state_t* state, pauli_op_t* op) {
    cplx_t* O_psi = apply_operator(op, state);  // |O|ψ⟩
    cplx_t result = dot_product(state->data, O_psi, state->dim);
    return result.real;  // Expectation value is real
}
```

**BLAS**: Use `cblas_zdotc` for dot product

#### 6.2.2 General Hermitian Observable

Similar algorithm, but operator is general Hermitian matrix.

### 6.3 Observable Moment Matrices

[To be specified - advanced feature]

Calculate matrix elements ⟨i|U†OU|j⟩ for custom basis defined by unitary U.

---

## 7. Memory Management

### 7.1 Allocation Strategy

#### 7.1.1 State Allocation

```c
state_t* state_init(qubit_t num_qubits, state_type_t type);
```

**Implementation**:
1. Allocate `state_t` structure
2. Calculate required data size
3. Allocate aligned memory for data array (64-byte alignment for SIMD)
4. Initialize state (default: |00...0⟩)

**Alignment**: Use `aligned_alloc()` or `posix_memalign()` for cache-line alignment

#### 7.1.2 Deallocation

```c
void state_free(state_t* state);
```

**Requirements**:
- Free `data` array
- Free `state_t` structure
- Set pointers to NULL (defensive programming)
- Idempotent: safe to call multiple times

### 7.2 Memory Pooling

**v1.0**: Not implemented (future optimization)

**Rationale**: Complex to implement correctly, marginal benefit for typical use cases

### 7.3 Memory Safety

**Tools**:
- AddressSanitizer for leak detection
- Valgrind for memory errors
- Manual testing with small qubit counts

**Best Practices**:
- Clear ownership semantics
- Document which functions allocate/free memory
- Provide examples of correct usage

---

## 8. Performance Optimization

### 8.1 BLAS Integration

#### 8.1.1 ILP64 Interface

**Requirement**: 64-bit integer BLAS interface

**Rationale**: State vectors for n > 16 qubits exceed 2^31 elements (int32 limit)

**Configuration**:
- macOS: Accelerate framework (native ILP64)
- Linux: OpenBLAS with `INTERFACE64=1`

#### 8.1.2 Key BLAS Operations

| Operation | BLAS Function | Usage |
|-----------|---------------|-------|
| Vector dot product | `cblas_zdotc` | Expectation values |
| Matrix-vector product | `cblas_zgemv` | Gate application (density matrices) |
| Matrix-matrix product | `cblas_zgemm` | Density matrix evolution |
| Vector scaling | `cblas_zdscal` | Normalization |

### 8.2 OpenMP Parallelization

#### 8.2.1 Parallel Regions

**Candidates**:
- Single-qubit gate application (parallelize over independent 2D subspaces)
- Probability computation (parallelize over state vector)
- Multi-qubit gate application (parallelize outer loops)

**Example**:
```c
#pragma omp parallel for schedule(static)
for (uint64_t i = 0; i < dim; i += stride) {
    // Apply gate to subspace
}
```

#### 8.2.2 Performance Considerations

- **Overhead**: OpenMP thread creation cost
- **Threshold**: Only parallelize for n ≥ 20 qubits (empirically determined)
- **Scalability**: Test on 4, 8, 16 core systems

### 8.3 Vectorization

#### 8.3.1 Compiler Hints

```c
#pragma GCC ivdep  // Ignore vector dependencies
#pragma clang loop vectorize(enable)
```

#### 8.3.2 Data Alignment

Ensure 64-byte alignment for AVX-512, 32-byte for AVX2.

### 8.4 Cache Optimization

**Strategy**: Access memory in contiguous blocks

**Example**: Process state vector in chunks that fit in L2 cache

---

## 9. Numerical Precision

### 9.1 Floating-Point Representation

**Type**: IEEE 754 double precision (64-bit)
- Significand: 53 bits (~15.9 decimal digits)
- Machine epsilon: ~2.22 × 10^-16

### 9.2 Error Accumulation

**Sources**:
1. Rounding errors in gate application
2. Normalization errors
3. BLAS library numerical errors

**Mitigation**:
- Use Kahan summation for critical accumulations
- Periodic renormalization (when |⟨ψ|ψ⟩ - 1| > ε)
- Validate unitary gate matrices (U†U = I)

### 9.3 Validation Tests

**Test Suite**:
- Identity gate preserves state (to machine precision)
- Norm conservation after gate application
- Hermiticity of density matrices
- Expectation values of projection operators ∈ [0, 1]

---

## 10. Platform-Specific Details

### 10.1 macOS (Accelerate Framework)

**BLAS**: vecLib (part of Accelerate)

**Linking**:
```cmake
target_link_libraries(qsim PRIVATE "-framework Accelerate")
```

**ILP64**: Native support via `cblas_*` functions

### 10.2 Linux (OpenBLAS)

**BLAS**: OpenBLAS 0.3.27+ (bundled)

**Build Options**:
```bash
INTERFACE64=1 BINARY=64 USE_OPENMP=1
```

**First Build**: ~15-30 minutes (compiles OpenBLAS)

**Subsequent Builds**: Fast (cached)

### 10.3 ARM Architecture

**Status**: Supported, but not optimized for v1.0

**Considerations**:
- NEON SIMD instructions (ARM64)
- Different cache hierarchy
- May require tuning for Raspberry Pi

### 10.4 Compiler-Specific Features

**GCC**:
- Use `-march=native` for automatic vectorization
- OpenMP support via `-fopenmp`

**Clang**:
- Similar optimizations as GCC
- Better diagnostic messages

---

## 11. Testing Strategy

### 11.1 Unit Tests

**Framework**: Unity

**Coverage**:
- State initialization (all basis states)
- Gate application (single, two-qubit)
- Measurement (deterministic cases)
- Memory management (allocation, free)

### 11.2 Integration Tests

**Scenarios**:
- Complete quantum circuits (Grover, VQE)
- Cross-platform consistency
- Performance benchmarks

### 11.3 Numerical Validation

**Reference**:
- Compare against known analytical results
- Cross-validate with Qiskit, Cirq (where applicable)

---

## 12. Build Configuration

### 12.1 CMake Options

```cmake
option(USE_SYSTEM_OPENBLAS "Use system OpenBLAS" OFF)
option(ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(BUILD_SHARED_LIBS "Build shared library" OFF)
```

### 12.2 Compiler Flags

**Debug**:
```cmake
-g -O0 -Wall -Wextra -Wpedantic
```

**Release**:
```cmake
-O3 -march=native -DNDEBUG -flto
```

**AddressSanitizer**:
```cmake
-fsanitize=address -fno-omit-frame-pointer
```

---

## 13. Open Questions & Decisions

[To be resolved during implementation]

1. **Circuit representation**: Immediate vs. deferred execution?
2. **Error handling**: Error codes vs. exceptions (in future C++ wrapper)?
3. **Thread safety**: Make API thread-safe or document as non-thread-safe?
4. **Random number generation**: Mersenne Twister vs. system RNG?

---

## Appendix A: Performance Benchmarks

[To be populated with actual benchmark results]

**Target Metrics**:
- Time to apply Hadamard gate to n qubits
- Time to apply CNOT gate to n qubits
- Time to measure n-qubit state
- Memory usage vs. theoretical minimum

---

## Appendix B: References

1. Nielsen & Chuang, "Quantum Computation and Quantum Information"
2. BLAS Quick Reference: http://www.netlib.org/blas/
3. OpenBLAS Documentation: https://github.com/OpenMathLib/OpenBLAS/wiki
4. Qiskit Aer Source (reference implementation): https://github.com/Qiskit/qiskit-aer

---

**Document Maintenance**: Update this spec as implementation decisions are made
