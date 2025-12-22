# Memory Layout: State Representation

**Module:** State Representation
**Last Updated:** 2025-12-22

---

## Overview

Quantum states are stored as complex-valued arrays (`cplx_t*`) with different layouts depending on whether the state is pure or mixed.

- **Pure states:** Linear array (state vector)
- **Mixed states:** Packed lower-triangular array (density matrix)

---

## Pure States (State Vectors)

### Memory Layout

A pure state with `n` qubits is represented as a state vector `|ψ⟩` in a Hilbert space of dimension `d = 2^n`.

**Array structure:**
```
data[0], data[1], data[2], ..., data[d-1]
```

**Interpretation:**
```
|ψ⟩ = data[0]|0⟩ + data[1]|1⟩ + data[2]|2⟩ + ... + data[d-1]|d-1⟩
```

Where `|i⟩` denotes the i-th computational basis state.

### Example: 1-Qubit Pure State

**Qubits:** 1
**Dimension:** `d = 2`
**Array length:** 2

```
Index | Basis State | Coefficient
------|-------------|-------------
  0   |     |0⟩     |   data[0]
  1   |     |1⟩     |   data[1]
```

**Example values:**
```c
// |0⟩ state
data[0] = 1.0 + 0.0*I;
data[1] = 0.0 + 0.0*I;

// |+⟩ = (|0⟩ + |1⟩)/√2 state
data[0] = 0.70710678 + 0.0*I;
data[1] = 0.70710678 + 0.0*I;

// |i⟩ = (|0⟩ + i|1⟩)/√2 state
data[0] = 0.70710678 + 0.0*I;
data[1] = 0.0 + 0.70710678*I;
```

### Example: 2-Qubit Pure State

**Qubits:** 2
**Dimension:** `d = 4`
**Array length:** 4

```
Index | Basis State | Coefficient
------|-------------|-------------
  0   |    |00⟩     |   data[0]
  1   |    |01⟩     |   data[1]
  2   |    |10⟩     |   data[2]
  3   |    |11⟩     |   data[3]
```

**Basis state encoding:** Index `i` represents qubit values in binary (little-endian):
- `i = 0`: `|00⟩` (qubit 0 = 0, qubit 1 = 0)
- `i = 1`: `|01⟩` (qubit 0 = 1, qubit 1 = 0)
- `i = 2`: `|10⟩` (qubit 0 = 0, qubit 1 = 1)
- `i = 3`: `|11⟩` (qubit 0 = 1, qubit 1 = 1)

**Example values:**
```c
// |00⟩ state
data = [1.0, 0.0, 0.0, 0.0]

// Bell state |Φ+⟩ = (|00⟩ + |11⟩)/√2
data = [0.70710678, 0.0, 0.0, 0.70710678]

// Uniform superposition |++⟩ = (|00⟩ + |01⟩ + |10⟩ + |11⟩)/2
data = [0.5, 0.5, 0.5, 0.5]
```

### Accessing Elements

**Direct indexing:** `data[i]` gives the amplitude of basis state `|i⟩`

**Getting amplitude of a specific basis state:**
```c
// Get amplitude of |101⟩ for 3-qubit state
qubit_t qubits = 3;
unsigned basis_state = 0b101;  // Binary: qubit0=1, qubit1=0, qubit2=1
cplx_t amplitude = state.data[basis_state];  // data[5]
```

**Setting amplitude:**
```c
state.data[5] = 0.5 + 0.3*I;
```

---

## Mixed States (Density Matrices)

### Memory Layout

A mixed state with `n` qubits is represented as a density matrix `ρ` of dimension `d × d` where `d = 2^n`.

**Storage format:** Lower triangle in **column-major packed format** (LAPACK convention)

**Why only lower triangle?**
- Density matrices are Hermitian: `ρ† = ρ`
- This means `ρ[i,j] = conj(ρ[j,i])`
- Only need to store `i ≥ j` (lower triangle); upper triangle is redundant

**Array length:** `d * (d + 1) / 2`

### Packed Storage Indexing

For a `d × d` Hermitian matrix, the lower triangle is stored column-by-column:

```
Column 0: elements (0,0)
Column 1: elements (1,0), (1,1)
Column 2: elements (2,0), (2,1), (2,2)
...
Column j: elements (j,0), (j,1), ..., (j,j)
```

**Index formula:** Element `ρ[i,j]` where `i ≥ j` is stored at:
```
packed_index = j * d - j * (j - 1) / 2 + (i - j)
             = j * (2*d - j + 1) / 2 + i - j
```

Simplified for `i ≥ j`:
```c
dim_t packed_index(dim_t i, dim_t j, dim_t d) {
    if (i < j) {
        // Upper triangle: use Hermitian symmetry
        dim_t temp = i;
        i = j;
        j = temp;
    }
    return j * (2 * d - j - 1) / 2 + i;
}
```

**Accessing upper triangle:** Use Hermitian symmetry:
```c
// For i < j (upper triangle)
cplx_t element = conj(data[packed_index(j, i, d)]);
```

### Example: 1-Qubit Mixed State

**Qubits:** 1
**Dimension:** `d = 2`
**Array length:** `2 * 3 / 2 = 3`

**Matrix representation:**
```
ρ = [ ρ00  ρ01 ]
    [ ρ10  ρ11 ]
```

**Packed storage:**
```
Index | Matrix Element | Description
------|----------------|-------------
  0   |     ρ[0,0]     | Diagonal (population of |0⟩)
  1   |     ρ[1,0]     | Off-diagonal (coherence)
  2   |     ρ[1,1]     | Diagonal (population of |1⟩)
```

**Note:** `ρ[0,1] = conj(ρ[1,0])` (not stored, use Hermitian property)

**Example values:**
```c
// Pure state |0⟩ (ρ = |0⟩⟨0|)
data[0] = 1.0 + 0.0*I;  // ρ[0,0]
data[1] = 0.0 + 0.0*I;  // ρ[1,0]
data[2] = 0.0 + 0.0*I;  // ρ[1,1]

// Pure state |+⟩ (ρ = |+⟩⟨+| = [[0.5, 0.5], [0.5, 0.5]])
data[0] = 0.5 + 0.0*I;  // ρ[0,0]
data[1] = 0.5 + 0.0*I;  // ρ[1,0]
data[2] = 0.5 + 0.0*I;  // ρ[1,1]

// Maximally mixed state (ρ = I/2 = [[0.5, 0], [0, 0.5]])
data[0] = 0.5 + 0.0*I;  // ρ[0,0]
data[1] = 0.0 + 0.0*I;  // ρ[1,0]
data[2] = 0.5 + 0.0*I;  // ρ[1,1]
```

### Example: 2-Qubit Mixed State

**Qubits:** 2
**Dimension:** `d = 4`
**Array length:** `4 * 5 / 2 = 10`

**Matrix representation:**
```
ρ = [ ρ00  ρ01  ρ02  ρ03 ]
    [ ρ10  ρ11  ρ12  ρ13 ]
    [ ρ20  ρ21  ρ22  ρ23 ]
    [ ρ30  ρ31  ρ32  ρ33 ]
```

**Packed storage (column-by-column):**
```
Index | Matrix Element | Position in Matrix
------|----------------|--------------------
  0   |     ρ[0,0]     | (row 0, col 0) - diagonal
  1   |     ρ[1,0]     | (row 1, col 0)
  2   |     ρ[2,0]     | (row 2, col 0)
  3   |     ρ[3,0]     | (row 3, col 0)
  4   |     ρ[1,1]     | (row 1, col 1) - diagonal
  5   |     ρ[2,1]     | (row 2, col 1)
  6   |     ρ[3,1]     | (row 3, col 1)
  7   |     ρ[2,2]     | (row 2, col 2) - diagonal
  8   |     ρ[3,2]     | (row 3, col 2)
  9   |     ρ[3,3]     | (row 3, col 3) - diagonal
```

**Diagonal indices:** 0, 4, 7, 9
**Pattern:** For element `(i,i)`, index is `i * (2*d - i + 1) / 2`

**Example: Maximally mixed state** (ρ = I/4):
```c
// All off-diagonals are zero, diagonals are 0.25
data[0] = 0.25;  // ρ[0,0]
data[1] = 0.0;   // ρ[1,0]
data[2] = 0.0;   // ρ[2,0]
data[3] = 0.0;   // ρ[3,0]
data[4] = 0.25;  // ρ[1,1]
data[5] = 0.0;   // ρ[2,1]
data[6] = 0.0;   // ρ[3,1]
data[7] = 0.25;  // ρ[2,2]
data[8] = 0.0;   // ρ[3,2]
data[9] = 0.25;  // ρ[3,3]
```

### Example: 3-Qubit Mixed State

**Qubits:** 3
**Dimension:** `d = 8`
**Array length:** `8 * 9 / 2 = 36`

**Diagonal element indices:**
```
ρ[0,0] → index 0
ρ[1,1] → index 8
ρ[2,2] → index 15
ρ[3,3] → index 21
ρ[4,4] → index 26
ρ[5,5] → index 30
ρ[6,6] → index 33
ρ[7,7] → index 35
```

---

## Visual Layout Comparison

### 2-Qubit Pure State (4 elements)
```
Memory:  [  0  ][  1  ][  2  ][  3  ]
         |  c0  | c1   | c2   | c3  |
         |______|______|______|_____|

Interpretation: |ψ⟩ = c0|00⟩ + c1|01⟩ + c2|10⟩ + c3|11⟩
```

### 2-Qubit Mixed State (10 elements)
```
Matrix Layout:               Packed Storage:
   col 0  col 1  col 2  col 3      [0][1][2][3] [4][5][6] [7][8] [9]
r0 [ 0    ·     ·     ·    ]         |___col0__| |_col1__| |col2| |c3|
r1 [ 1    4    ·     ·    ]
r2 [ 2    5    7    ·    ]       · = upper triangle (not stored)
r3 [ 3    6    8    9    ]       Use Hermitian symmetry to access
```

---

## Calculating Trace (for validation)

**Pure state:** Trace is not defined (state vectors don't have a trace)

**Mixed state:** `Tr(ρ) = sum of diagonal elements = 1` (for valid density matrix)

```c
double calculate_trace(const state_t *state) {
    if (state->type != MIXED) return 0.0;  // Not applicable

    dim_t d = POW2(state->qubits, dim_t);
    double trace = 0.0;

    for (dim_t i = 0; i < d; i++) {
        // Diagonal element ρ[i,i] is at packed index:
        dim_t idx = i * (2 * d - i + 1) / 2;
        trace += creal(state->data[idx]);
    }

    return trace;
}
```

**Note:** For a valid density matrix, `trace ≈ 1.0` (within numerical precision).

---

## Element Access Helper Functions (Recommended)

### For Pure States
```c
// Get amplitude of basis state |i⟩
cplx_t state_get_amplitude(const state_t *state, dim_t i) {
    return state->data[i];
}

// Set amplitude
void state_set_amplitude(state_t *state, dim_t i, cplx_t value) {
    state->data[i] = value;
}
```

### For Mixed States
```c
// Get matrix element ρ[i,j]
cplx_t state_get_element(const state_t *state, dim_t i, dim_t j) {
    dim_t d = POW2(state->qubits, dim_t);

    if (i >= j) {
        // Lower triangle (stored)
        dim_t idx = j * (2 * d - j - 1) / 2 + i;
        return state->data[idx];
    } else {
        // Upper triangle (use Hermitian symmetry)
        dim_t idx = i * (2 * d - i - 1) / 2 + j;
        return conj(state->data[idx]);
    }
}

// Set matrix element ρ[i,j] (also sets ρ[j,i] = conj(value) implicitly)
void state_set_element(state_t *state, dim_t i, dim_t j, cplx_t value) {
    dim_t d = POW2(state->qubits, dim_t);

    if (i >= j) {
        dim_t idx = j * (2 * d - j - 1) / 2 + i;
        state->data[idx] = value;
    } else {
        // Store in lower triangle with conjugate
        dim_t idx = i * (2 * d - i - 1) / 2 + j;
        state->data[idx] = conj(value);
    }
}
```

**Note:** These functions are NOT currently implemented but would be useful additions.

---

## BLAS/LAPACK Compatibility

The packed storage format is compatible with LAPACK Hermitian packed (HP) routines:

| Operation | LAPACK Routine | Description |
|-----------|----------------|-------------|
| Matrix-vector multiply | `zhpmv()` | y = α*A*x + β*y |
| Rank-1 update | `zhpr()` | A = A + α*x*x† |
| Eigenvalues | `zhpev()` | Eigendecomposition |
| Matrix-matrix multiply | `zhpmm()` | C = α*A*B + β*C |

**Storage parameter:** Use `uplo = 'L'` (lower triangle) when calling LAPACK routines.

---

## Memory Footprint Comparison

| Qubits | Dimension | Pure State (elements) | Mixed State (elements) | Savings |
|--------|-----------|----------------------|------------------------|---------|
| 1      | 2         | 2                    | 3                      | -50%    |
| 2      | 4         | 4                    | 10                     | -150%   |
| 3      | 8         | 8                    | 36                     | -350%   |
| 4      | 16        | 16                   | 136                    | -750%   |
| 10     | 1024      | 1,024                | 524,800                | 51,200% |

**Key observation:** Mixed states scale as O(4^n), pure states as O(2^n)

**Memory usage:**
- Each `cplx_t` is 16 bytes (2 × 8-byte doubles)
- 10-qubit pure state: 16 KB
- 10-qubit mixed state: 8.2 MB
- Packed storage saves 50% compared to full matrix storage

---

## Summary

| State Type | Array Length | Element Access | Complexity |
|------------|--------------|----------------|------------|
| **PURE**   | `2^n`        | `data[i]`      | Simple     |
| **MIXED**  | `2^n(2^n+1)/2` | Formula or helper | Complex |

**Recommendation:** For mixed states, implement helper functions for element access to hide the packed indexing complexity.
