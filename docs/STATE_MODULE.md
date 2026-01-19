# State Module: Implementation Status

**Module:** State Representation
**Header:** `include/state.h`
**Implementation:** `src/state.c`
**Last Updated:** 2026-01-19

---

## Overview

Pure states: State vectors with `2^n` complex elements
Mixed states: Density matrices with `2^n(2^n+1)/2` elements (lower-triangle packed storage)
Storage format: LAPACK Hermitian packed (column-major), compatible with `zhpmv()`, `zhpr()`, `zhpevd()`

**Memory Layout:**
- Pure: `data[i]` = amplitude of basis state `|i⟩`
- Mixed: `data[j*(2d-j-1)/2 + i]` = element ρ[i,j] for i≥j; use conjugate for upper triangle

---

## Data Structures

| Type | Status | Description |
|------|--------|-------------|
| `state_type_t` | ✅ Complete | Enum: PURE, MIXED |
| `state_t` | ✅ Complete | Struct: type, data, qubits |
| `cplx_t` | ✅ Complete | Platform-agnostic complex type (q_types.h) |
| `dim_t` | ✅ Complete | 64-bit dimension type (q_types.h) |

---

## Implemented Functions

| Function | Status | Location | Description |
|----------|--------|----------|-------------|
| `state_free()` | ✅ Complete | state.c:23-28 | Deallocate and reset |
| `state_len()` | ✅ Complete | state.c:31-40 | Calculate array size |
| `state_init()` | ✅ Complete | state.c:43-63 | Initialize with ownership transfer |
| `state_plus()` | ✅ Complete | state.c:66-89 | Create uniform superposition |
| `state_cp()` | ✅ Complete | state.c:92-106 | Deep copy using cblas_zcopy() |

**Key:** Ownership transfer via `**data`; errors → stderr + data=NULL

---

## Unit Tests

| Test Function | Status | Target | Coverage |
|---------------|--------|--------|----------|
| `test_state_len_pure()` | ✅ Complete | state_len() | Verify 2^n for n=1,2,3,4 |
| `test_state_len_mixed()` | ✅ Complete | state_len() | Verify d(d+1)/2 for n=1,2,3 |
| `test_state_init_pure_null_data()` | ✅ Complete | state_init() | Zero initialization |
| `test_state_init_pure_with_data()` | ✅ Complete | state_init() | Ownership transfer verification |
| `test_state_init_mixed_null_data()` | ✅ Complete | state_init() | Zero initialization |
| `test_state_init_mixed_with_data()` | ✅ Complete | state_init() | Ownership transfer verification |
| `test_state_init_reinitialization()` | ✅ Complete | state_init() | Memory leak check on reinit |
| `test_state_free_pure()` | ✅ Complete | state_free() | Deallocation and reset |
| `test_state_free_mixed()` | ✅ Complete | state_free() | Double-free safety |
| `test_state_plus_pure_normalization()` | ✅ Complete | state_plus() | Verify 1/√(2^n), sum\|a\|²=1 |
| `test_state_plus_mixed_trace()` | ✅ Complete | state_plus() | Verify Tr(ρ)=1, all entries=1/2^n |
| `test_state_cp_pure_deep_copy()` | ✅ Complete | state_cp() | Independent allocations |
| `test_state_cp_mixed_deep_copy()` | ✅ Complete | state_cp() | Independent allocations |
| `test_state_single_qubit()` | ✅ Complete | Edge cases | 1-qubit correctness |
| `test_state_many_qubits()` | ✅ Complete | Edge cases | Large n, overflow check |

**Total:** 15 tests — ALL PASSING ✅
**File:** `test/src/test_state.c`

---

## Validation Functions (Missing)

| Function | Priority | Implementation | Status |
|----------|----------|----------------|--------|
| `state_is_normalized()` | High | cblas_dznrm2() | 🔲 TODO |
| `state_is_valid_density_matrix()` | High | Check Tr(ρ)=1, eigenvalues≥0 | 🔲 TODO |
| `state_purity()` | Medium | Pure: return 1.0; Mixed: Tr(ρ²) | 🔲 TODO |
| `state_partial_trace()` | Low | Trace out qubits → reduced ρ | 🔲 Future |

---

## Helper Functions (Missing)

| Function | Priority | Description | Status |
|----------|----------|-------------|--------|
| `packed_index()` | Medium | Convert (i,j) to packed index | 🔲 TODO |
| `state_get_element()` | Medium | Safe access for mixed states | 🔲 TODO |
| `state_set_element()` | Medium | Safe modification for mixed states | 🔲 TODO |
| `state_normalize_inplace()` | Medium | In-place normalization | 🔲 TODO |
| `state_scale_inplace()` | Medium | In-place scaling | 🔲 TODO |

**Packed index formula:**
`i >= j ? j*(2*d - j - 1)/2 + i : conj(data[i*(2*d - i - 1)/2 + j])`

---

## Documentation Tasks

| Task | Priority | Status |
|------|----------|--------|
| Doxygen comments for public API | High | 🔲 TODO |
| Memory ownership documentation | High | 🔲 TODO |
| Usage pattern examples | Medium | 🔲 TODO |

---

## BLAS/LAPACK Functions

| Operation | Function | Usage |
|-----------|----------|-------|
| Vector norm | `cblas_dznrm2()` | Normalization check |
| Vector copy | `cblas_zcopy()` | Used in state_cp() |
| Vector scaling | `cblas_zdscal()` | Normalization |
| Outer product | `cblas_zher()` | ρ = \|ψ⟩⟨ψ\| (not needed per user) |
| Matrix multiply | `cblas_zhemm()` | Tr(ρ²) calculation |
| Eigendecomposition | `zheevd()` | Purity, validity checks |
| HP matrix-vector | `cblas_zhpmv()` | Hermitian packed operations |
| HP rank-1 update | `cblas_zhpr()` | Hermitian packed operations |
| HP eigenvalues | `zhpevd()` | Hermitian packed eigenvalues |

**Note:** Use `uplo = 'L'` for lower triangle in LAPACK calls.

---

## Design Decisions

| Topic | Current | Question | Status |
|-------|---------|----------|--------|
| Memory ownership | Transfer via `**data` | Add non-transfer init? | 🔲 Undecided |
| Error handling | `qs_error_t` return codes | — | ✅ Implemented |
| Thread safety | Not thread-safe | Document guarantees? | 🔲 Undecided |
| Qubit limits | 255 (unsigned char) | Sufficient? | 🔲 Undecided |
| Normalization | Manual | Auto-normalize gates? | 🔲 Undecided |

**Error Codes** (defined in `q_types.h`):
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

## Test Infrastructure (Complete)

| Component | Status | Location |
|-----------|--------|----------|
| Unity framework | ✅ Complete | test/include/test.h |
| Complex assertions | ✅ Complete | TEST_ASSERT_COMPLEX_WITHIN |
| Computational basis generator | ✅ Complete | test_cb_pure() |
| Hadamard basis generator | ✅ Complete | test_xb_pure() |
| Circular basis generator | ✅ Complete | test_yb_pure() |

---

## Implementation Priorities

**Phase 1:** Unit tests, validation functions, Doxygen
**Phase 2:** Helper functions, purity calculation, in-place ops
**Phase 3:** Partial trace, eigendecomposition

---

## Key Memory Layout Facts

**Pure states (2-qubit):** 4 elements
`|ψ⟩ = c₀|00⟩ + c₁|01⟩ + c₂|10⟩ + c₃|11⟩`

**Mixed states (2-qubit):** 10 elements (lower triangle, column-major)
Diagonal indices: 0, 4, 7, 9
Formula: `diag_idx(i) = i*(2*d - i + 1)/2`

**Trace calculation:**
```c
for (i = 0; i < d; i++) {
    trace += creal(data[i*(2*d - i + 1)/2]);
}
```

**Memory footprint:**
- Each cplx_t: 16 bytes
- 10-qubit pure: 16 KB
- 10-qubit mixed: 8.2 MB
- Packed storage saves 50% vs full matrix

**Basis encoding:** Little-endian binary
`i=5` → `|101⟩` → q0=1, q1=0, q2=1

---

## Common Pitfalls

1. **Shallow copy:** Use `state_cp()`, not `state2 = state1`
2. **Ownership transfer:** After `state_init(&s, n, &data)`, `data` is NULL
3. **Uninitialized type:** Always set `.type` before `state_init()`
4. **Thread safety:** Requires external synchronization for concurrent writes
5. **Error checking:** Check `state.data != NULL` after initialization
