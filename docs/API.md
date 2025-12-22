# qSim API Reference v1.0

**Document Version**: 1.0
**Last Updated**: 2025-12-22
**Status**: Draft
**Related**: [PRD.md](../PRD.md), [TECHNICAL_SPEC.md](./TECHNICAL_SPEC.md)

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Getting Started](#2-getting-started)
3. [Core Types](#3-core-types)
4. [State Management](#4-state-management)
5. [Quantum Gates](#5-quantum-gates)
6. [Circuit Operations](#6-circuit-operations)
7. [Measurement](#7-measurement)
8. [Utilities](#8-utilities)
9. [Error Handling](#9-error-handling)
10. [Examples](#10-examples)

---

## 1. Introduction

This document provides a complete API reference for the qSim quantum computing simulation library.

### 1.1 Header Files

```c
#include <qlib.h>  // Main header (includes all public APIs)
```

Alternatively, include specific headers:

```c
#include <q_types.h>  // Type definitions
#include <state.h>    // State management
#include <gates.h>    // Gate operations (v1.0)
#include <circuit.h>  // Circuit operations (v1.0)
```

### 1.2 Linking

**CMake**:
```cmake
find_package(qsim REQUIRED)
target_link_libraries(your_target PRIVATE qsim::qsim)
```

**Manual**:
```bash
gcc -o program program.c -lqsim -lm
```

---

## 2. Getting Started

### 2.1 Minimal Example

```c
#include <qlib.h>
#include <stdio.h>

int main() {
    // Create a 3-qubit state in |000⟩
    state_t* state = state_init(3, STATE_PURE);
    if (!state) {
        fprintf(stderr, "Failed to initialize state\n");
        return 1;
    }

    // Apply Hadamard gate to qubit 0
    // gate_h(state, 0);  // v1.0 API TBD

    // Measure all qubits
    // uint64_t result = measure_all(state);  // v1.0 API TBD

    // Clean up
    state_free(state);

    return 0;
}
```

### 2.2 Build and Run

```bash
gcc -o example example.c -lqsim
./example
```

---

## 3. Core Types

### 3.1 Complex Numbers

```c
typedef struct {
    double real;
    double imag;
} cplx_t;
```

**Functions**:

```c
// Create complex number
cplx_t cplx_make(double real, double imag);

// Complex arithmetic (utilities)
cplx_t cplx_add(cplx_t a, cplx_t b);
cplx_t cplx_sub(cplx_t a, cplx_t b);
cplx_t cplx_mul(cplx_t a, cplx_t b);
cplx_t cplx_conj(cplx_t a);
double cplx_abs(cplx_t a);
double cplx_abs2(cplx_t a);  // |a|^2 (faster)
```

**Macros**:

```c
#define CPLX_ZERO  ((cplx_t){0.0, 0.0})
#define CPLX_ONE   ((cplx_t){1.0, 0.0})
#define CPLX_I     ((cplx_t){0.0, 1.0})
```

### 3.2 Qubit and Dimension Types

```c
typedef uint32_t qubit_t;  // Qubit identifier (0 to n-1)
typedef uint64_t dim_t;    // Hilbert space dimension (2^n)
```

### 3.3 State Type

```c
typedef enum {
    STATE_PURE,    // Pure state (statevector)
    STATE_MIXED    // Mixed state (density matrix)
} state_type_t;
```

---

## 4. State Management

### 4.1 state_t Structure

```c
typedef struct {
    state_type_t type;    // PURE or MIXED
    cplx_t* data;         // State data
    qubit_t num_qubits;   // Number of qubits
    dim_t dim;            // Hilbert space dimension
} state_t;
```

**Note**: This structure is opaque in v1.0. Do not access fields directly.

### 4.2 State Initialization

#### state_init

```c
state_t* state_init(qubit_t num_qubits, state_type_t type);
```

**Description**: Creates a new quantum state initialized to |00...0⟩.

**Parameters**:
- `num_qubits`: Number of qubits (0 < num_qubits ≤ 30 for pure states)
- `type`: `STATE_PURE` or `STATE_MIXED`

**Returns**: Pointer to newly allocated state, or `NULL` on failure

**Example**:
```c
state_t* state = state_init(5, STATE_PURE);  // 5-qubit state
```

**Error Conditions**:
- Returns `NULL` if memory allocation fails
- Returns `NULL` if `num_qubits` is 0 or exceeds limits

---

#### state_init_basis

```c
state_t* state_init_basis(qubit_t num_qubits, uint64_t basis_index);
```

**Description**: Creates a state initialized to a specific computational basis state.

**Parameters**:
- `num_qubits`: Number of qubits
- `basis_index`: Index of basis state (0 ≤ basis_index < 2^num_qubits)

**Returns**: Pointer to state in |basis_index⟩, or `NULL` on failure

**Example**:
```c
// Create 3-qubit state |101⟩ (binary 101 = 5)
state_t* state = state_init_basis(3, 5);
```

---

#### state_plus

```c
void state_plus(state_t* state);
```

**Description**: Reinitializes state to uniform superposition |+⟩^⊗n.

**Parameters**:
- `state`: State to modify (must be pure state)

**Example**:
```c
state_t* state = state_init(3, STATE_PURE);
state_plus(state);  // Now in |+++⟩ = (|000⟩ + |001⟩ + ... + |111⟩)/√8
```

**Notes**:
- Only works for pure states
- Equivalent to applying Hadamard to all qubits of |00...0⟩

---

### 4.3 State Copying

#### state_cp

```c
state_t* state_cp(const state_t* src);
```

**Description**: Creates a deep copy of a quantum state.

**Parameters**:
- `src`: Source state to copy

**Returns**: New state identical to `src`, or `NULL` on failure

**Example**:
```c
state_t* original = state_init(5, STATE_PURE);
state_t* copy = state_cp(original);
```

**Notes**:
- Allocates new memory for copy
- Caller is responsible for freeing both states

---

### 4.4 State Deallocation

#### state_free

```c
void state_free(state_t* state);
```

**Description**: Frees all memory associated with a quantum state.

**Parameters**:
- `state`: State to free (may be `NULL`)

**Example**:
```c
state_t* state = state_init(3, STATE_PURE);
// ... use state ...
state_free(state);
state = NULL;  // Good practice
```

**Notes**:
- Safe to call with `NULL` pointer (no-op)
- After calling, pointer is invalid

---

### 4.5 State Inspection

#### state_len

```c
dim_t state_len(const state_t* state);
```

**Description**: Returns the number of complex elements in state data array.

**Parameters**:
- `state`: State to query

**Returns**:
- Pure state: 2^n
- Mixed state: 2^n * (2^n + 1) / 2

**Example**:
```c
state_t* state = state_init(5, STATE_PURE);
dim_t len = state_len(state);  // Returns 32 (2^5)
```

---

#### state_get_amplitude

```c
cplx_t state_get_amplitude(const state_t* state, uint64_t index);
```

**Description**: Retrieves amplitude of a specific basis state (pure states only).

**Parameters**:
- `state`: Quantum state (must be pure)
- `index`: Basis state index (0 ≤ index < 2^n)

**Returns**: Complex amplitude α_index

**Example**:
```c
state_t* state = state_init_basis(3, 5);  // |101⟩
cplx_t amp = state_get_amplitude(state, 5);  // Returns (1.0, 0.0)
cplx_t amp0 = state_get_amplitude(state, 0); // Returns (0.0, 0.0)
```

---

## 5. Quantum Gates

[To be completed in v1.0 implementation]

### 5.1 Single-Qubit Gates

#### gate_x, gate_y, gate_z (Pauli Gates)

```c
int gate_x(state_t* state, qubit_t target);
int gate_y(state_t* state, qubit_t target);
int gate_z(state_t* state, qubit_t target);
```

**Description**: Apply Pauli X, Y, or Z gate to target qubit.

**Parameters**:
- `state`: Quantum state to modify
- `target`: Target qubit index (0 ≤ target < num_qubits)

**Returns**: 0 on success, error code on failure

---

#### gate_h (Hadamard Gate)

```c
int gate_h(state_t* state, qubit_t target);
```

**Description**: Apply Hadamard gate to target qubit.

**Matrix**:
```
H = (1/√2) [1   1]
           [1  -1]
```

---

#### gate_s, gate_t (Phase Gates)

```c
int gate_s(state_t* state, qubit_t target);
int gate_t(state_t* state, qubit_t target);
int gate_sdg(state_t* state, qubit_t target);  // S†
int gate_tdg(state_t* state, qubit_t target);  // T†
```

**Description**: Apply phase gates S, T, and their adjoints.

---

#### gate_rx, gate_ry, gate_rz (Rotation Gates)

```c
int gate_rx(state_t* state, qubit_t target, double theta);
int gate_ry(state_t* state, qubit_t target, double theta);
int gate_rz(state_t* state, qubit_t target, double theta);
```

**Description**: Apply rotation gates around X, Y, or Z axis.

**Parameters**:
- `state`: Quantum state
- `target`: Target qubit
- `theta`: Rotation angle in radians

**Example**:
```c
gate_rx(state, 0, M_PI/4);  // Rotate qubit 0 by π/4 around X
```

---

### 5.2 Two-Qubit Gates

#### gate_cnot (Controlled-NOT)

```c
int gate_cnot(state_t* state, qubit_t control, qubit_t target);
```

**Description**: Apply CNOT gate with specified control and target qubits.

**Parameters**:
- `state`: Quantum state
- `control`: Control qubit index
- `target`: Target qubit index

**Example**:
```c
gate_cnot(state, 1, 0);  // Qubit 1 controls, qubit 0 is target
```

---

#### gate_swap

```c
int gate_swap(state_t* state, qubit_t q1, qubit_t q2);
```

**Description**: Swap the states of two qubits.

---

### 5.3 Multi-Qubit Gates

#### gate_toffoli (CCNOT)

```c
int gate_toffoli(state_t* state, qubit_t ctrl1, qubit_t ctrl2, qubit_t target);
```

**Description**: Apply Toffoli (CCNOT) gate with two control qubits.

---

### 5.4 Controlled Gates

```c
int gate_controlled(state_t* state, qubit_t control, qubit_t target,
                    const cplx_t U[4]);
```

**Description**: Apply arbitrary controlled single-qubit gate.

**Parameters**:
- `U`: 2×2 unitary matrix (row-major order)

---

## 6. Circuit Operations

[To be completed in v1.0 implementation]

### 6.1 Circuit Creation

```c
circuit_t* circuit_create(qubit_t num_qubits);
void circuit_free(circuit_t* circuit);
```

### 6.2 Circuit Construction

```c
int circuit_add_gate_x(circuit_t* circuit, qubit_t target);
int circuit_add_gate_h(circuit_t* circuit, qubit_t target);
int circuit_add_gate_cnot(circuit_t* circuit, qubit_t control, qubit_t target);
// ... similar for all gates
```

### 6.3 Circuit Execution

```c
int circuit_execute(circuit_t* circuit, state_t* state);
```

**Description**: Execute entire circuit on quantum state.

---

## 7. Measurement

[To be completed in v1.0 implementation]

### 7.1 Computational Basis Measurement

#### measure_all

```c
uint64_t measure_all(state_t* state);
```

**Description**: Measure all qubits in computational basis, collapse state.

**Returns**: Measurement outcome as bit string (little-endian)

**Example**:
```c
uint64_t outcome = measure_all(state);
// outcome = 0b101 means qubits measured as |101⟩
```

---

#### measure_qubit

```c
int measure_qubit(state_t* state, qubit_t target);
```

**Description**: Measure single qubit, collapse appropriately.

**Returns**: 0 or 1 (measurement outcome)

---

### 7.2 Expectation Values

#### expectation_pauli

```c
double expectation_pauli(const state_t* state, const char* pauli_string);
```

**Description**: Compute expectation value of Pauli observable.

**Parameters**:
- `state`: Quantum state (pure or mixed)
- `pauli_string`: String like "XYZ" (qubit 0=X, 1=Y, 2=Z)

**Returns**: Real expectation value

**Example**:
```c
double exp_val = expectation_pauli(state, "ZZI");  // ⟨Z₀Z₁I₂⟩
```

---

#### expectation_hermitian

```c
double expectation_hermitian(const state_t* state, const cplx_t* H, dim_t dim);
```

**Description**: Compute expectation value of arbitrary Hermitian operator.

**Parameters**:
- `H`: Hermitian matrix (2^n × 2^n)
- `dim`: Matrix dimension (must equal state dimension)

---

### 7.3 Observable Moment Matrices

[Advanced feature - API TBD]

```c
int observable_moments(const state_t* state, const cplx_t* O,
                       const cplx_t* U, cplx_t* result);
```

---

## 8. Utilities

### 8.1 State Validation

```c
bool state_is_normalized(const state_t* state, double tolerance);
bool state_is_valid(const state_t* state);
```

### 8.2 State Printing (Debug)

```c
void state_print(const state_t* state, FILE* fp);
```

**Description**: Print state amplitudes to file (for debugging).

**Example**:
```c
state_print(state, stdout);
```

**Output** (for 2-qubit state):
```
|00⟩: (0.500000, 0.000000)
|01⟩: (0.500000, 0.000000)
|10⟩: (0.500000, 0.000000)
|11⟩: (0.500000, 0.000000)
```

---

## 9. Error Handling

### 9.1 Error Codes

```c
typedef enum {
    QSIM_OK = 0,
    QSIM_ERR_ALLOC,          // Memory allocation failed
    QSIM_ERR_INVALID_ARG,    // Invalid argument
    QSIM_ERR_OUT_OF_RANGE,   // Index out of range
    QSIM_ERR_NOT_NORMALIZED, // State not normalized
    QSIM_ERR_NOT_UNITARY,    // Matrix not unitary
    QSIM_ERR_UNKNOWN         // Unknown error
} qsim_error_t;
```

### 9.2 Error Messages

```c
const char* qsim_error_string(qsim_error_t error);
```

**Returns**: Human-readable error message

### 9.3 Error Handling Pattern

```c
int result = gate_h(state, 0);
if (result != QSIM_OK) {
    fprintf(stderr, "Error: %s\n", qsim_error_string(result));
    return result;
}
```

---

## 10. Examples

### 10.1 Bell State Preparation

```c
#include <qlib.h>

int main() {
    // Create 2-qubit state |00⟩
    state_t* state = state_init(2, STATE_PURE);

    // Create Bell state |Φ+⟩ = (|00⟩ + |11⟩)/√2
    gate_h(state, 0);        // Apply H to qubit 0
    gate_cnot(state, 0, 1);  // CNOT with control=0, target=1

    // Measure
    uint64_t outcome = measure_all(state);
    printf("Measured: %llu\n", outcome);  // 0 (|00⟩) or 3 (|11⟩)

    state_free(state);
    return 0;
}
```

### 10.2 Grover's Algorithm (Sketch)

```c
#include <qlib.h>

int main() {
    const qubit_t n = 3;
    state_t* state = state_init(n, STATE_PURE);

    // Initialize to uniform superposition
    for (qubit_t i = 0; i < n; i++) {
        gate_h(state, i);
    }

    // Grover iteration
    const int iterations = 2;  // ~π/4 * √(2^n) optimal
    for (int iter = 0; iter < iterations; iter++) {
        // Oracle (marks target state)
        oracle_mark(state, target);

        // Diffusion operator
        grover_diffusion(state, n);
    }

    // Measure
    uint64_t result = measure_all(state);
    printf("Found: %llu\n", result);

    state_free(state);
    return 0;
}
```

### 10.3 Expectation Value Calculation

```c
#include <qlib.h>

int main() {
    state_t* state = state_init(2, STATE_PURE);

    // Prepare state |+⟩|0⟩
    gate_h(state, 0);

    // Compute ⟨X₀ ⊗ Z₁⟩
    double exp_val = expectation_pauli(state, "XZ");
    printf("Expectation value: %f\n", exp_val);  // Should be 1.0

    state_free(state);
    return 0;
}
```

---

## 11. Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2025-12-22 | Initial API specification |

---

## 12. See Also

- [TECHNICAL_SPEC.md](./TECHNICAL_SPEC.md) - Implementation details
- [ARCHITECTURE.md](./ARCHITECTURE.md) - System architecture
- [PRD.md](../PRD.md) - Product requirements

---

**Note**: This API is subject to change during v1.0 development. Final API will be documented upon release.
