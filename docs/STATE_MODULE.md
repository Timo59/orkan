# State Module Technical Specification

**Module:** Quantum State Representation
**Header:** `include/state.h`
**Implementation:** `src/state.c`
**Dependencies:** `q_types.h`, `utils.h`, BLAS (`cblas_zcopy`)

---

## 1. Overview

The state module provides data structures and functions for representing quantum states of n-qubit systems. It supports three storage formats:

| Type | Mathematical Object | Storage Format |
|------|---------------------|----------------|
| `PURE` | State vector ∈ ℂ^(2^n) | Dense array |
| `MIXED_PACKED` | Density matrix ∈ ℂ^(2^n × 2^n) | Hermitian packed (lower triangle, column-major) |
| `MIXED_TILED` | Density matrix ∈ ℂ^(2^n × 2^n) | Blocked layout (tiles of lower triangle) |

**Design principles:**
- Runtime dispatch via type enum
- No derived values stored: `dim`, `n_tiles` computed on demand
- Assert on programmer errors, NULL on out-of-memory

---

## 2. Type Definitions

### 2.1 State Type Enum

```c
typedef enum state_type {
    PURE,           // State vector |ψ⟩
    MIXED_PACKED,   // Density matrix ρ, lower triangle column-major
    MIXED_TILED     // Density matrix ρ, blocked/tiled layout
} state_type_t;
```

### 2.2 State Structure

```c
typedef struct state {
    state_type_t type;   // Storage format (must be set before init)
    cplx_t *data;        // Heap-allocated array (owned by struct)
    qubit_t qubits;      // Number of qubits n
} state_t;
```

**Invariants:**
- `type` must be set before calling `state_init()` or `state_len()`
- `data != NULL` implies valid, 64-byte-aligned allocation of `state_len()` elements
- `data == NULL` implies uninitialized, freed, or allocation failure

### 2.3 Supporting Types (from q_types.h)

| Type | Definition | Size | Purpose |
|------|------------|------|---------|
| `cplx_t` | Platform-specific double complex | 16 bytes | Complex amplitudes |
| `dim_t` | BLAS integer, **must be 64-bit** (ILP64) | 8 bytes | Array dimensions |
| `qubit_t` | `unsigned char` | 1 byte | Qubit count (max 255) |

`dim_t` is 64-bit on both platforms: `__LAPACK_int` on macOS (vecLib) and `blasint` on Linux (OpenBLAS compiled with ILP64). This guarantees all index and length formulas are overflow-safe for any memory-feasible qubit count.

---

## 3. Constants

```c
#ifndef LOG_TILE_DIM
#define LOG_TILE_DIM 5
#endif
#define TILE_DIM (1 << LOG_TILE_DIM)
```

**Location:** `state.h`

**Rationale for default (LOG_TILE_DIM=5 → TILE_DIM=32):**
- 32 × 32 × 16 bytes = 16 KB per tile
- 1 tile fits in L1 cache (~32-64 KB), enabling fast within-tile operations
- 4 tiles fit in L2 cache (~128 KB), enabling cross-tile operations without thrashing

Override at compile time with `-DLOG_TILE_DIM=6` for 64×64 tiles.

---

## 4. Memory Layouts

### 4.1 PURE: State Vector

For n qubits, dimension d = 2^n:

```
|ψ⟩ = Σᵢ cᵢ|i⟩    where i ∈ {0, 1, ..., d-1}
```

**Storage:** `data[i]` = amplitude cᵢ of basis state |i⟩

**Basis encoding:** Little-endian binary
- Index i = 5 (binary 101) → |101⟩ → qubit₀ = 1, qubit₁ = 0, qubit₂ = 1

**Example (2 qubits):**
```
Index:  0     1     2     3
State: |00⟩  |01⟩  |10⟩  |11⟩
```

**Array length:** 2^n

### 4.2 MIXED_PACKED: Hermitian Packed Lower Triangle

For n qubits, dimension d = 2^n, the density matrix ρ is Hermitian (ρ† = ρ). Store only the lower triangle in column-major order:

```
ρ = [ρ₀₀  ρ₀₁* ρ₀₂* ρ₀₃*]     Packed storage (column-major lower):
    [ρ₁₀  ρ₁₁  ρ₁₂* ρ₁₃*]     [ρ₀₀, ρ₁₀, ρ₂₀, ρ₃₀, ρ₁₁, ρ₂₁, ρ₃₁, ρ₂₂, ρ₃₂, ρ₃₃]
    [ρ₂₀  ρ₂₁  ρ₂₂  ρ₂₃*]      idx:  0    1    2    3    4    5    6    7    8    9
    [ρ₃₀  ρ₃₁  ρ₃₂  ρ₃₃ ]
```

**Index formula for ρ[row, col] where row ≥ col:**
```c
packed_idx = col * d - col * (col + 1) / 2 + row
```

**Upper triangle access:** `ρ[row, col] = conj(ρ[col, row])` for row < col

**Array length:** d(d + 1)/2

**LAPACK compatibility:** Use `uplo = 'L'` for `zhpmv()`, `zhpr()`, `zhpevd()`

### 4.3 MIXED_TILED: Blocked Layout

Tiles of size TILE_DIM × TILE_DIM, stored in row-major tile order, lower-triangular at tile level.

```
Tile grid (n_tiles × n_tiles where n_tiles = ceil(d / TILE_DIM)):
- Store tiles T(tr, tc) where tr ≥ tc (lower triangular)
- Number of stored tiles: n_tiles × (n_tiles + 1) / 2
- Within each tile: full TILE_DIM × TILE_DIM, row-major
```

**Derived quantities (computed, not stored):**
```c
dim_t dim = 1 << qubits;
dim_t n_tiles = (dim + TILE_DIM - 1) / TILE_DIM;  // ceil(dim / TILE_DIM)
```

**Index formula:**
```c
// Tile containing element (row, col)
tile_row = row / TILE_DIM;
tile_col = col / TILE_DIM;

// Tile offset in storage (lower-triangular, row-major tile ordering)
tile_offset = tile_row * (tile_row + 1) / 2 + tile_col;

// Position within tile (row-major)
local_row = row % TILE_DIM;
local_col = col % TILE_DIM;
elem_offset = local_row * TILE_DIM + local_col;

// Full index
index = tile_offset * TILE_DIM * TILE_DIM + elem_offset;
```

**Hermitian access:** For row < col at tile level, access via conjugate of (col, row).

**Intra-tile layout is row-major.** Do not pass tile pointers to LAPACK routines (`zhpmv`, `zhpr`, etc.) which expect column-major packed storage.

**Padding:** When `dim < TILE_DIM` (i.e., `qubits < LOG_TILE_DIM`), a single full tile is allocated but only the `dim × dim` corner is physical. Padding elements must remain zero. When `dim` is a multiple of `TILE_DIM` (i.e., `qubits ≥ LOG_TILE_DIM`), no padding exists.

**Array length:** n_tiles × (n_tiles + 1) / 2 × TILE_DIM²

---

## 5. API Reference

### 5.1 state_len

```c
dim_t state_len(const state_t *state);
```

Returns array length for the state's data buffer (see Section 4 for formulas).

---

### 5.2 state_init

```c
void state_init(state_t *state, qubit_t qubits, cplx_t **data);
```

Initialize state with given qubit count. Frees previous allocation if `state->data` was non-NULL.

- **`data == NULL`:** Allocates a zero-initialized, 64-byte-aligned buffer via `aligned_alloc`.
- **`data != NULL`:** Transfers ownership (`*data` set to NULL). If `*data` is already 64-byte-aligned, takes it directly; otherwise copies into aligned storage and frees the original.

---

### 5.3 state_free

```c
void state_free(state_t *state);
```

Deallocate data buffer, set `data = NULL` and `qubits = 0`. Safe to call on NULL or already-freed state.

---

### 5.4 state_plus

```c
void state_plus(state_t *state, qubit_t qubits);
```

Initialize to uniform superposition |+⟩^⊗n = (1/√d) Σᵢ |i⟩ for PURE, or ρ = |+⟩⟨+|^⊗n = (1/d) Σᵢⱼ |i⟩⟨j| for mixed types.

---

### 5.5 state_cp

```c
state_t state_cp(const state_t *state);
```

Return deep copy with independent data allocation. Uses `cblas_zcopy()`.

---

### 5.6 state_print

```c
void state_print(const state_t *state);
```

Print state to stdout. Uses `vprint()` for PURE, `mprint_packed()` for MIXED_PACKED, `mprint_tiled()` for MIXED_TILED.

---

### 5.7 state_get

```c
cplx_t state_get(const state_t *state, dim_t row, dim_t col);
```

Get element at (row, col). For PURE, asserts `col == 0` and returns `data[row]`. For mixed types, handles Hermitian symmetry: returns `conj(data[index(col, row)])` when row < col.

---

### 5.8 state_set

```c
void state_set(state_t *state, dim_t row, dim_t col, cplx_t val);
```

Set element at (row, col). For PURE, asserts `col == 0`. For mixed types, handles Hermitian symmetry: stores `conj(val)` at (col, row) when row < col.

---

## 6. Error Handling

### 6.1 Strategy

| Error Type | Handling | Rationale |
|------------|----------|-----------|
| Programmer error | `assert()` | Bugs should crash immediately for easy debugging |
| Out-of-memory | Set `data = NULL`, print to stderr | Runtime condition, caller can handle |

### 6.2 Programmer Errors (Asserted)

- `state == NULL`
- `state->type` not in `{PURE, MIXED_PACKED, MIXED_TILED}`
- `state->data == NULL` when data access required (e.g., `state_cp`)

### 6.3 Runtime Errors (NULL + stderr)

- `malloc`/`calloc` returns NULL

### 6.4 Caller Responsibility

After `state_init()` or `state_plus()`, check:
```c
if (state.data == NULL) {
    // Handle allocation failure
}
```

---

## 7. Thread Safety

No internal shared state. Concurrent reads of same state are safe; concurrent modifications require caller synchronization. Gate module achieves thread safety through disjoint memory regions in butterfly decomposition.

---

## 8. Memory Ownership

### 8.1 Rules

1. **`state_t` owns its data:** After initialization, the struct owns the `data` pointer
2. **Ownership transfer:** `state_init(..., &data)` transfers ownership; `data` becomes NULL
3. **No shared ownership:** Each buffer has exactly one owner
4. **Cleanup:** `state_free()` releases ownership and deallocates

### 8.2 Ownership Transfer Protocol

```c
cplx_t *buffer = malloc(len * sizeof(cplx_t));
// ... populate buffer ...

state_init(&state, n, &buffer);
assert(buffer == NULL);  // Ownership transferred

state_free(&state);  // Proper cleanup
```

### 8.3 Reinitialization

Calling `state_init()` on an initialized state frees the previous allocation:
```c
state_init(&state, 2, NULL);  // First init
state_init(&state, 3, NULL);  // Frees previous, allocates new
```

---

## 9. Memory Footprint

| Qubits | Dim | PURE | MIXED_PACKED | MIXED_TILED (TILE=32) |
|--------|-----|------|--------------|----------------------|
| 1 | 2 | 32 B | 48 B | 16 KB |
| 2 | 4 | 64 B | 160 B | 16 KB |
| 5 | 32 | 512 B | 8.4 KB | 16 KB |
| 8 | 256 | 4 KB | 526 KB | 576 KB |
| 10 | 1024 | 16 KB | 8.4 MB | 8.6 MB |
| 12 | 4096 | 64 KB | 134 MB | 136 MB |

Each `cplx_t` is 16 bytes. MIXED_TILED has overhead for small systems but converges to ~same as packed for n ≥ 8.

---

## 10. Platform Considerations

### 10.1 macOS (vecLib/Accelerate)

```c
#include <vecLib/cblas_new.h>
typedef __LAPACK_double_complex cplx_t;
typedef __LAPACK_int dim_t;
```

### 10.2 Linux (OpenBLAS)

```c
#include <cblas.h>
typedef openblas_complex_double cplx_t;
typedef blasint dim_t;
```

First build compiles OpenBLAS with ILP64 (~10-30 min).

### 10.3 Complex Number Operations

Use C99 functions: `creal()`, `cimag()`, `cabs()`, `conj()`

---

## 11. Utility Functions

**Location:** `include/utils.h`, `src/utils.c`

| Function | Purpose |
|----------|---------|
| `vprint(v, n)` | Print complex vector |
| `mprint_packed(m, n)` | Print Hermitian matrix from packed lower triangle |
| `mprint_tiled(m, n)` | Print Hermitian matrix from tiled layout |
| `packed_to_full(packed, full, dim)` | Convert packed to full matrix |
| `full_to_packed(full, packed, dim)` | Convert full matrix to packed |
| `tiled_to_full(tiled, full, dim)` | Convert tiled to full matrix |
| `full_to_tiled(full, tiled, dim)` | Convert full matrix to tiled |

---

## 12. Test Strategy

### 12.1 Unit Tests (test/src/test_state.c)

| Function | Test Coverage |
|----------|---------------|
| `state_len()` | All three types, various qubit counts |
| `state_init()` | NULL data (zero init), ownership transfer, reinitialization |
| `state_free()` | Deallocation, double-free safety |
| `state_plus()` | Normalization (pure), trace = 1 (mixed) |
| `state_cp()` | Deep copy independence |
| Edge cases | Single qubit, large qubit counts |

---

## 13. Usage Examples

### 13.1 Pure State |0⟩

```c
state_t psi = {.type = PURE};
state_init(&psi, 3, NULL);      // 3-qubit |000...0⟩
if (!psi.data) { /* handle OOM */ }

psi.data[0] = 1.0;              // Set to |000⟩
state_free(&psi);
```

### 13.2 Mixed State (Maximally Mixed)

```c
state_t rho = {.type = MIXED_PACKED};
state_init(&rho, 2, NULL);

dim_t d = 1 << rho.qubits;
for (dim_t i = 0; i < d; i++) {
    dim_t diag_idx = i * d - i * (i + 1) / 2 + i;
    rho.data[diag_idx] = 1.0 / d;  // ρ = I/d
}
state_free(&rho);
```

### 13.3 Uniform Superposition

```c
state_t plus = {.type = PURE};
state_plus(&plus, 4);           // |+⟩^⊗4
// All amplitudes = 1/√16 = 0.25

state_free(&plus);
```

### 13.4 Deep Copy

```c
state_t original = {.type = MIXED_TILED};
state_plus(&original, 6);

state_t copy = state_cp(&original);
// Modify original without affecting copy

state_free(&original);
state_free(&copy);
```

---

