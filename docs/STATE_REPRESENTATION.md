# State Representation - Technical Requirements

**Module:** State Representation
**Status:** In Progress
**Last Updated:** 2025-12-22

## Overview

This document specifies the technical requirements for quantum state representation in qSim, covering both pure states (state vectors) and mixed states (density matrices).

---

## Implementation Status

### ✅ Completed

#### Core Data Structures (`include/state.h`)
- [x] `state_type_t` enum: Distinguishes PURE vs MIXED states
- [x] `state_t` struct: Contains type, data pointer, and qubit count
- [x] Platform-agnostic complex type `cplx_t` (q_types.h:37-42)
- [x] Platform-agnostic dimension type `dim_t` for 64-bit indexing (q_types.h:37-42)

#### Memory Management Functions (`src/state.c`)
- [x] `state_free()`: Deallocates state data and resets metadata (state.c:23-28)
- [x] `state_len()`: Computes array size based on state type (state.c:31-40)
  - Pure states: `2^n` elements
  - Mixed states: `2^n * (2^n + 1) / 2` elements (lower triangle)
- [x] `state_init()`: Initializes state with optional data transfer (state.c:43-63)
  - Handles memory ownership transfer via double pointer
  - Zero-initializes if no data provided

#### State Initialization Functions (`src/state.c`)
- [x] `state_plus()`: Creates uniform superposition (state.c:66-89)
  - Pure: amplitudes = `1/sqrt(2^n)`
  - Mixed: entries = `1/2^n`
- [x] `state_cp()`: Deep copy using BLAS `cblas_zcopy()` (state.c:92-106)

#### Test Infrastructure (`test/include/test.h`, `test/src/test.c`)
- [x] Unity framework integration
- [x] Custom assertion macros for complex arrays with tolerance
- [x] Test state generators:
  - [x] Computational basis states (Z eigenstates)
  - [x] Hadamard basis states (X eigenstates)
  - [x] Circular basis states (Y eigenstates)
  - [x] Combined test state generator (3 × 2^n states)

---

## 🔲 Missing Requirements

### 1. **State Validation**

#### 1.1 Pure State Validation
**Priority:** High
**Functions needed:**
- `bool state_is_normalized(const state_t *state, double tol)`
  - Check: `|<ψ|ψ> - 1| < tol`
  - Implementation: `cblas_dznrm2()` for L2-norm

#### 1.2 Density Matrix Validation
**Priority:** High
**Functions needed:**
- `bool state_is_valid_density_matrix(const state_t *state, double tol)`
  - Check trace = 1: `Tr(ρ) = 1`
  - Check Hermiticity (implicit via lower-triangle storage)
  - Check positive semi-definite: eigenvalues ≥ 0
  - Implementation: LAPACK eigenvalue routines

### 2. **State Properties**

#### 2.1 Purity Calculation
**Priority:** Medium
**Functions needed:**
- `double state_purity(const state_t *state)`
  - Pure states: Return 1.0 directly
  - Mixed states: Compute `Tr(ρ²)` using BLAS matrix multiplication

#### 2.2 Reduced Density Matrices
**Priority:** Low (future feature)
**Functions needed:**
- `state_t state_partial_trace(const state_t *state, const qubit_t *qubits, qubit_t n_trace_out)`
  - Trace out specified qubits
  - Returns reduced density matrix

### 3. **Conversion Between Representations**

#### 3.1 Pure to Mixed Conversion
**Priority:** High (needed for mixed-state simulations)
**Functions needed:**
- `state_t state_pure_to_mixed(const state_t *pure_state)`
  - Compute `ρ = |ψ⟩⟨ψ|` (outer product)
  - Store only lower triangle
  - Implementation: `cblas_zher()` (Hermitian rank-1 update) or custom loop

#### 3.2 Mixed to Pure Projection
**Priority:** Low
**Functions needed:**
- `state_t state_mixed_to_pure(const state_t *mixed_state, double purity_threshold)`
  - Check if state is pure enough (purity > threshold)
  - Extract dominant eigenvector if pure
  - Implementation: LAPACK `zheevd()` for eigendecomposition

### 4. **Memory Optimization**

#### 4.1 In-Place Operations
**Priority:** Medium
**Enhancement:** Add in-place variants where applicable:
- `void state_normalize_inplace(state_t *state)`
- `void state_scale_inplace(state_t *state, cplx_t scalar)`

#### 4.2 Packed Storage Utilities
**Priority:** Medium
**Functions needed:**
- `dim_t packed_index(dim_t i, dim_t j, dim_t dim)`
  - Convert (row, col) to packed lower-triangle index
  - Formula: `i >= j ? i*(i+1)/2 + j : j*(j+1)/2 + i`
- `cplx_t state_get_element(const state_t *state, dim_t i, dim_t j)`
  - Safe element access for mixed states
- `void state_set_element(state_t *state, dim_t i, dim_t j, cplx_t value)`
  - Safe element modification

### 5. **Documentation**

#### 5.1 API Documentation
**Priority:** High
**Needed:**
- Doxygen-style comments for all public functions
- Memory ownership documentation (who owns what after function calls)
- Example usage patterns

#### 5.2 Memory Layout Documentation
**Priority:** High
**Needed:**
- Document lower-triangle storage scheme with examples
- Illustrate packed indexing for 2-qubit and 3-qubit systems
- BLAS/LAPACK function mapping (which BLAS routines operate on packed storage)

### 6. **Testing**

#### 6.1 Unit Tests for state.c
**Priority:** High
**File:** `test/src/test_state.c` (currently nearly empty)
**Tests needed:**
- `state_init()` with NULL and non-NULL data
- `state_free()` doesn't leak memory
- `state_len()` correct for both types, various qubit counts
- `state_plus()` normalization and uniformity
- `state_cp()` deep copy verification (modify original, check copy unchanged)
- Memory ownership transfer in `state_init()`

#### 6.2 Edge Case Tests
**Priority:** Medium
**Tests needed:**
- Zero qubits (should fail gracefully or be documented as unsupported)
- Maximum qubit count (test dimension overflow with `dim_t`)
- Allocation failure handling (requires mock allocator)

#### 6.3 Mixed State Tests
**Priority:** High (once pure→mixed conversion implemented)
**Tests needed:**
- Lower-triangle storage correctness
- Hermiticity preservation
- Trace preservation after operations

---

## Technical Decisions to Document

### 1. **Memory Ownership**
**Current behavior:** `state_init()` with non-NULL `data` parameter transfers ownership via double pointer.
**Question:** Should we also provide a non-transfer initialization (copy data instead)?

### 2. **Error Handling**
**Current behavior:** Print to stderr and return/set invalid state.
**Question:** Should we provide error codes or errno-style error reporting?

### 3. **Thread Safety**
**Current status:** Not thread-safe (modifies state_t in place).
**Question:** Document thread-safety guarantees? Require external synchronization?

### 4. **Qubit Limits**
**Current:** `qubit_t` is `unsigned char` → max 255 qubits.
**Question:** Is this sufficient? (2^255 is astronomically large, likely not a practical limit)

### 5. **Normalization Policy**
**Current:** `state_plus()` creates normalized states, but no automatic normalization after operations.
**Question:** Should gate operations auto-normalize, or is this the user's responsibility?

---

## Implementation Priority

### Phase 1: Core Completion (Immediate)
1. Complete `test/src/test_state.c` with comprehensive unit tests
2. Add state validation functions (normalization check, density matrix validity)
3. Implement `state_pure_to_mixed()` conversion
4. Document memory layout and ownership

### Phase 2: Utilities (After gates module)
5. Add packed storage helper functions
6. Implement purity calculation
7. Add in-place operation variants

### Phase 3: Advanced Features (Future)
8. Partial trace for subsystem analysis
9. Mixed-to-pure projection for near-pure states
10. Eigendecomposition-based state analysis

---

## BLAS/LAPACK Functions Relevant to State Representation

| Operation | BLAS/LAPACK Function | Notes |
|-----------|---------------------|-------|
| Vector norm | `cblas_dznrm2()` | L2 norm for normalization check |
| Vector copy | `cblas_zcopy()` | Used in `state_cp()` |
| Vector scaling | `cblas_zdscal()` | For normalization |
| Outer product | `cblas_zher()` | Hermitian rank-1 update (ρ = |ψ⟩⟨ψ|) |
| Matrix-matrix multiply | `cblas_zhemm()` | For Tr(ρ²) calculation |
| Eigendecomposition | `zheevd()` (LAPACK) | Purity, validity checks |

**Note:** Packed storage (lower triangle) requires specific BLAS variants:
- `cblas_zhpmv()`: Hermitian packed matrix-vector multiplication
- `cblas_zhpr()`: Hermitian packed rank-1 update
- `zhpevd()`: Eigenvalues of Hermitian packed matrix
