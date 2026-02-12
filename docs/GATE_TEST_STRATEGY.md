# Gate Module: Comprehensive Unit Test Strategy

**Version:** 1.0
**Date:** 2026-02-06
**Status:** Analysis Complete

---

## Executive Summary

This document provides a comprehensive test strategy for the quantum gate module, covering all gate types (single-qubit, two-qubit, and three-qubit) across three storage formats (PURE, MIXED_PACKED, MIXED_TILED). The analysis identifies current coverage, gaps, and provides a concrete implementation plan with test counts.

**Key Findings:**
- **Current Coverage**: 11/31 gates implemented and tested (PURE + MIXED_PACKED only)
- **Missing Coverage**: 20 gates, entire MIXED_TILED format
- **Total Test Count Target**: ~1M tests across 31 gates × 3 formats
- **Priority**: Complete MIXED_TILED infrastructure, then implement remaining two-qubit and three-qubit gates

---

## 1. Current Test Infrastructure Status

### 1.1 What Exists and Works ✅

#### Test Harnesses (Complete)
| Harness | Storage Type | Gate Type | Status | Location |
|---------|--------------|-----------|--------|----------|
| `testSingleQubitGate()` | PURE | Single-qubit | ✅ | test/src/test_qhipster.c |
| `testSingleQubitGateMixed()` | MIXED_PACKED | Single-qubit | ✅ | test/src/test_mhipster.c |
| `testRotationGate()` | PURE | Parameterized | ✅ | test/src/test_qhipster.c |
| `testRotationGateMixed()` | MIXED_PACKED | Parameterized | ✅ | test/src/test_mhipster.c |
| `testTwoQubitGate()` | PURE | Two-qubit | ✅ | test/src/test_qhipster.c |
| `testTwoQubitGateMixed()` | MIXED_PACKED | Two-qubit | ✅ | test/src/test_mhipster.c |

#### Test State Generators (Complete)
| Generator | Type | States Generated | Count |
|-----------|------|------------------|-------|
| `test_cb_pure()` | PURE | Computational basis |0⟩, |1⟩ | 2 |
| `test_xb_pure()` | PURE | Hadamard basis |+⟩, |-⟩ | 2 |
| `test_yb_pure()` | PURE | Circular basis |+i⟩, |-i⟩ | 2 |
| `test_bell_pure()` | PURE | Bell states (n=2) | 4 |
| `test_ghz_pure()` | PURE | GHZ state | 1 |
| `test_w_pure()` | PURE | W state | 1 |
| `test_*_mixed()` | MIXED | All above as density matrices | Same |
| `test_maximally_mixed()` | MIXED | I/2^n | 1 |
| `test_random_mixture()` | MIXED | Random statistical mixtures | 10-20 |

#### Reference Implementation Utilities (Complete)
| Function | Purpose | Status |
|----------|---------|--------|
| `zmv()` | Matrix-vector product (pure reference) | ✅ |
| `zumu()` | U·ρ·U† computation (mixed reference) | ✅ |
| `zkron()` | Kronecker product for matrix building | ✅ |
| `mat_id()` | Identity matrix builder | ✅ |
| `mat_single_qubit_gate()` | Build I⊗...⊗U⊗...⊗I | ✅ |
| `mat_two_qubit_gate()` | Build full two-qubit gate matrix | ✅ |
| `mat_rx/ry/rz()` | Rotation gate matrix builders | ✅ |
| `density_pack()` | Full → packed storage | ✅ |
| `density_unpack()` | Packed → full storage | ✅ |

#### Implemented Gates (11/31 complete)
| Gate | Pure | Mixed_Packed | Mixed_Tiled | Tests Passing |
|------|------|--------------|-------------|---------------|
| X | ✅ | ✅ | ❌ | Yes |
| Y | ✅ | ✅ | ❌ | Yes |
| Z | ✅ | ✅ | ❌ | Yes |
| H | ✅ | ✅ | ❌ | Yes |
| S | ✅ | ✅ | ❌ | Yes |
| Sdg | ✅ | ✅ | ❌ | Yes |
| T | ✅ | ✅ | ❌ | Yes |
| Tdg | ✅ | ✅ | ❌ | Yes |
| Rx | ✅ | ✅ | ❌ | Yes |
| Ry | ✅ | ✅ | ❌ | Yes |
| Rz | ✅ | ✅ | ❌ | Yes |
| CX | ✅ | ✅ | ❌ | Yes |

### 1.2 Coverage Analysis by Gate Type

#### Single-Qubit Gates (11/13 complete for PURE+MIXED_PACKED)

**Test Coverage per Gate:**
- **System sizes**: n = 1, 2, 3, 4 qubits
- **Target positions**: All positions (0 to n-1)
- **Test states**:
  - n=1: 6 states (cb, xb, yb)
  - n=2: 10 states (6 + 4 Bell)
  - n=3: 12 states (6 + 4 Bell + GHZ + W)
  - n=4: 12 states (same as n=3, Bell states are n=2 only)
- **Mixed states**: Same pure states as density matrices + maximally mixed + 10 random mixtures

**Test count formula (non-rotation gates):**
```
n=1: 1 position × 6 states = 6 tests
n=2: 2 positions × 10 states = 20 tests
n=3: 3 positions × 12 states = 36 tests
n=4: 4 positions × 12 states = 48 tests
Total per gate (pure): 110 tests

Mixed adds:
n=2: 2 positions × (10 + 1 + 10) = 42 tests
n=3: 3 positions × (12 + 1 + 10) = 69 tests
n=4: 4 positions × (12 + 1 + 10) = 92 tests
Total per gate (mixed): 203 tests

Grand total per gate: ~313 tests
```

**Note:** Actual counts from GATES_MODULE.md show ~5,900 pure + ~5,900 mixed = ~11,800 per gate. This suggests more test states or larger n values are being used than documented.

**Rotation gates** (Rx, Ry, Rz): 6 angles × base tests = ~1,878 tests per gate

#### Two-Qubit Gates (1/4 complete)

**Test Coverage per Gate:**
- **System sizes**: n = 2, 3, 4 qubits
- **Qubit pairs**: All (q1, q2) where q1 ≠ q2
  - n=2: 2 qubits → 2×1 = 2 pairs
  - n=3: 3 qubits → 3×2 = 6 pairs
  - n=4: 4 qubits → 4×3 = 12 pairs
- **Test states**: Same as single-qubit (entangled states crucial here)

**Test count formula:**
```
n=2: 2 pairs × 10 states = 20 tests
n=3: 6 pairs × 12 states = 72 tests
n=4: 12 pairs × 12 states = 144 tests
Total per gate (pure): 236 tests

Mixed (similar ratio):
Total per gate (mixed): ~400 tests

Grand total per gate: ~636 tests
```

**Note:** GATES_MODULE.md estimates ~18,800 pure + ~18,900 mixed = ~37,700 per two-qubit gate, again suggesting more coverage.

#### Three-Qubit Gates (0/1 complete)

**Toffoli (CCX) only**

**Test Coverage:**
- **System sizes**: n = 3, 4 qubits (minimum 3 for three-qubit gate)
- **Qubit triples**: All (control1, control2, target) where all different
  - n=3: 3 qubits → 3×2×1 = 6 triples
  - n=4: 4 qubits → 4×3×2 = 24 triples
- **Test states**: Same as above

**Test count estimate:**
```
n=3: 6 triples × 12 states = 72 tests
n=4: 24 triples × 12 states = 288 tests
Total per gate (pure): 360 tests
Total per gate (mixed): ~600 tests
Grand total: ~960 tests
```

**Note:** GATES_MODULE.md estimates ~20,000 pure + ~20,000 mixed = ~40,000 per three-qubit gate. Needs verification.

---

## 2. Missing Test Coverage

### 2.1 Missing Gates (20 gates)

#### Single-Qubit Gates (2 gates)
| Gate | Description | Priority | Notes |
|------|-------------|----------|-------|
| `hy()` | Hadamard-Y | Medium | Clifford gate, less common |
| `p(θ)` | Phase rotation P(θ) | High | Parameterized, used in many algorithms |

#### Two-Qubit Gates (17 gates)
| Gate | Description | Priority | Notes |
|------|-------------|----------|-------|
| `cy()` | Controlled-Y | High | Pauli gate variant |
| `cz()` | Controlled-Z | High | Pauli gate variant, symmetric |
| `swap()` | SWAP | High | Fundamental two-qubit operation |
| `cs()` | Controlled-S | Medium | Phase gate variant |
| `csdg()` | Controlled-S† | Medium | Phase gate variant |
| `ch()` | Controlled-H | Medium | Clifford variant |
| `chy()` | Controlled-Hy | Low | Less common |
| `ct()` | Controlled-T | Medium | Non-Clifford variant |
| `ctdg()` | Controlled-T† | Medium | Non-Clifford variant |
| `cp(θ)` | Controlled-P(θ) | High | Parameterized phase |
| `cpdg(θ)` | Controlled-P†(θ) | Medium | Less common |
| `rswap()` | Rotated SWAP | Low | Specialized gate |

**Controlled gate pattern:** Most controlled gates follow CX pattern: control=1 applies gate to target

#### Three-Qubit Gates (1 gate)
| Gate | Description | Priority | Notes |
|------|-------------|----------|-------|
| `ccx()` | Toffoli (CCNOT) | High | Fundamental for reversible computing |

### 2.2 Missing Tiled Storage Tests (ALL gates × MIXED_TILED)

**Status:** BLOCKED_MEMORY_LAYOUT.md spec exists, but no implementation or tests

**Required Infrastructure:**
1. State module for tiled format:
   - `state_tiled_init()`, `state_tiled_free()`, etc.
   - `state_tiled_get()`, `state_tiled_set()` (handles Hermitian conjugation)
   - `state_tiled_from_packed()`, `state_tiled_to_packed()` (for test conversion)

2. Test harnesses for tiled format:
   - `testSingleQubitGateTiled()`
   - `testRotationGateTiled()`
   - `testTwoQubitGateTiled()`
   - `testThreeQubitGateTiled()`

3. Tiled pack/unpack utilities:
   - `density_unpack_tiled()` - tiled → full matrix
   - `density_pack_tiled()` - full matrix → tiled

4. Gate implementations:
   - All 11 existing gates × tiled format
   - All 20 missing gates × tiled format

**Test Strategy for Tiled:** Independent duplicate of packed tests
- Use same test state generators
- Convert packed → tiled for input
- Apply tiled gate
- Extract tiled → full matrix
- Compare with same BLAS reference (zumu)

---

## 3. Edge Cases and Important Test Scenarios

### 3.1 Qubit Position Edge Cases

#### Adjacent vs Non-Adjacent Qubits (Two-Qubit Gates)

**Why it matters:** Memory access patterns differ significantly

| Scenario | Example | Memory Access Pattern |
|----------|---------|----------------------|
| Adjacent qubits | q1=0, q2=1 on n=3 | Consecutive blocks in memory |
| Non-adjacent | q1=0, q2=2 on n=3 | Strided access with gap |
| Maximally separated | q1=0, q2=n-1 | Maximum stride |

**Test coverage:** ALL pairs are tested by current harness

#### Control vs Target Ordering

**Question:** Does control > target vs control < target exercise different code paths?

**Analysis:**
- **Physical symmetry:** For CZ (symmetric gate), order shouldn't matter
- **Implementation:** qHiPSTER uses bit masks - different orderings MAY affect:
  - Bit insertion logic
  - Memory access patterns
  - Cache behavior

**Current coverage:** All (q1, q2) pairs tested ✅ (includes both orderings)

**CX-specific:** CX(control, target) is NOT symmetric
- CX(0, 1) ≠ CX(1, 0)
- Current test covers both directions ✅

#### Boundary Qubits (Leftmost and Rightmost)

**Why it matters:**
- **Rightmost (qubit 0):** Smallest stride, most cache-friendly
- **Leftmost (qubit n-1):** Largest stride, different OpenMP strategy

**Edge cases to verify:**
1. Single-qubit gate on qubit 0 (rightmost)
2. Single-qubit gate on qubit n-1 (leftmost)
3. Two-qubit gate with q1=0 or q2=0
4. Two-qubit gate with q1=n-1 or q2=n-1
5. Three-qubit gate with any qubit at boundary

**Current coverage:** All positions tested ✅

### 3.2 Small System Edge Cases

#### Minimum System Sizes

| Gate Type | Minimum n | Why |
|-----------|-----------|-----|
| Single-qubit | n=1 | Tests gate in isolation |
| Two-qubit | n=2 | Tests gate with no spectator qubits |
| Three-qubit | n=3 | Tests gate with no spectator qubits |

**Current coverage:**
- Single-qubit: n=1,2,3,4 ✅
- Two-qubit: n=2,3,4 ✅
- Three-qubit: n=3,4 planned ✅

#### Tiled Format Edge Case (n < LOG_TILE_DIM)

**Issue:** When dim < TILE_DIM (e.g., n < 5 with TILE_DIM=32):
- Single tile allocated
- Only dim×dim corner is physical
- Padding must remain zero

**Critical tests:**
1. n=1,2,3,4 with tiled format (all fit in single tile)
2. Verify padding stays zero after gate operations
3. n=5 transition case (dim=32 = TILE_DIM exactly)
4. n=6 first case requiring multiple tiles

**Status:** Not implemented ⚠️

### 3.3 State-Specific Edge Cases

#### Pure State Edge Cases

| State | Why Test | What It Catches |
|-------|----------|-----------------|
| \|0⟩^⊗n | All-zero computational | Amplitude errors at zero |
| \|1⟩^⊗n | All-one computational | Bit flip handling |
| \|+⟩^⊗n | Uniform superposition | Phase errors (all equal real) |
| Bell states | Entanglement propagation | Two-qubit gate correctness |
| GHZ state | Multi-partite entanglement | Global phase handling |

**Current coverage:** All included ✅

#### Mixed State Edge Cases

| State | Why Test | What It Catches |
|-------|----------|-----------------|
| Pure as mixed | Rank-1 density matrix | Basic correctness |
| Maximally mixed | I/2^n (rank=2^n) | Trace preservation |
| Random mixtures | Arbitrary eigenvalues | Full conjugation UρU† |

**Current coverage:** All included ✅

### 3.4 Numerical Edge Cases

#### Rotation Gates: Special Angles

**Current angles tested:**
- θ = 0 (identity)
- θ = π/4 (small rotation)
- θ = π/2 (quarter turn)
- θ = π (half turn)
- θ = 3π/2 (three-quarter turn)
- θ = -π/3 (negative angle)

**What these catch:**
- θ=0: Gate becomes identity
- θ=π: Maximum rotation
- Negative θ: Sign handling
- Irrational multiples of π: Rounding errors

**Recommendation:** Current coverage is EXCELLENT ✅

#### Precision Tolerance

**Current:** ε = 1e-12

**Rationale:**
- Double precision: ~15-16 decimal digits
- Accumulated rounding over O(2^n) operations
- Sparse operations have better stability than dense

**Edge case tests needed:**
1. Norm preservation: |⟨ψ'|ψ'⟩ - 1| < ε
2. Trace preservation: |Tr(ρ') - Tr(ρ)| < ε
3. Hermiticity: |ρ[i,j] - conj(ρ[j,i])| < ε

**Status:** Norm/trace tested implicitly, Hermiticity not explicitly tested ⚠️

---

## 4. Test Matrix: Gate × State Type × Coverage

### 4.1 Single-Qubit Gates

| Gate | PURE | MIXED_PACKED | MIXED_TILED | Test Count (per format) | Priority |
|------|------|--------------|-------------|------------------------|----------|
| X | ✅ | ✅ | ❌ | ~5,900 | DONE |
| Y | ✅ | ✅ | ❌ | ~5,900 | DONE |
| Z | ✅ | ✅ | ❌ | ~5,900 | DONE |
| H | ✅ | ✅ | ❌ | ~5,900 | DONE |
| S | ✅ | ✅ | ❌ | ~5,900 | DONE |
| Sdg | ✅ | ✅ | ❌ | ~5,900 | DONE |
| T | ✅ | ✅ | ❌ | ~5,900 | DONE |
| Tdg | ✅ | ✅ | ❌ | ~5,900 | DONE |
| Hy | ❌ | ❌ | ❌ | ~5,900 | P2 |
| P(θ) | ❌ | ❌ | ❌ | ~35,000 (6 angles) | P2 |
| Rx(θ) | ✅ | ✅ | ❌ | ~35,000 | DONE |
| Ry(θ) | ✅ | ✅ | ❌ | ~35,000 | DONE |
| Rz(θ) | ✅ | ✅ | ❌ | ~35,000 | DONE |

**Subtotal:** 13 gates × 3 formats = 39 implementations needed
- **Complete:** 11 gates × 2 formats = 22 implementations (56%)
- **Missing:** 17 implementations (44%)

### 4.2 Two-Qubit Gates

| Gate | PURE | MIXED_PACKED | MIXED_TILED | Test Count (per format) | Priority |
|------|------|--------------|-------------|------------------------|----------|
| CX | ✅ | ✅ | ❌ | ~18,800 | DONE |
| CY | ❌ | ❌ | ❌ | ~18,800 | P3 |
| CZ | ❌ | ❌ | ❌ | ~18,800 | P3 |
| SWAP | ❌ | ❌ | ❌ | ~18,800 | P3 |
| CS | ❌ | ❌ | ❌ | ~18,800 | P4 |
| CSdg | ❌ | ❌ | ❌ | ~18,800 | P4 |
| CH | ❌ | ❌ | ❌ | ~18,800 | P4 |
| CHy | ❌ | ❌ | ❌ | ~18,800 | P4 |
| CT | ❌ | ❌ | ❌ | ~18,800 | P4 |
| CTdg | ❌ | ❌ | ❌ | ~18,800 | P4 |
| CP(θ) | ❌ | ❌ | ❌ | ~113,000 (6 angles) | P3 |
| CPdg(θ) | ❌ | ❌ | ❌ | ~113,000 | P4 |
| RSWAP | ❌ | ❌ | ❌ | ~18,800 | P4 |

**Subtotal:** 13 gates × 3 formats = 39 implementations needed
- **Complete:** 1 gate × 2 formats = 2 implementations (5%)
- **Missing:** 37 implementations (95%)

### 4.3 Three-Qubit Gates

| Gate | PURE | MIXED_PACKED | MIXED_TILED | Test Count (per format) | Priority |
|------|------|--------------|-------------|------------------------|----------|
| CCX (Toffoli) | ❌ | ❌ | ❌ | ~20,000 | P5 |

**Subtotal:** 1 gate × 3 formats = 3 implementations needed
- **Complete:** 0 implementations (0%)
- **Missing:** 3 implementations (100%)

### 4.4 Grand Total

**Total implementations needed:** 31 gates × 3 formats = 93 implementations
**Complete:** 24 implementations (26%)
**Missing:** 69 implementations (74%)

**Test count projection:**
- Single-qubit: 13 gates × ~5,900 tests × 3 formats = ~230,000 tests
- Single-qubit rotation: 4 gates × ~35,000 tests × 3 formats = ~420,000 tests
- Two-qubit: 11 gates × ~18,800 tests × 3 formats = ~620,000 tests
- Two-qubit parameterized: 2 gates × ~113,000 tests × 3 formats = ~678,000 tests
- Three-qubit: 1 gate × ~20,000 tests × 3 formats = ~60,000 tests

**Grand total: ~2M tests** (revised from 1M due to tiled format inclusion)

---

## 5. Infrastructure Needs

### 5.1 Missing Test Infrastructure

#### Tiled Storage Infrastructure (HIGH PRIORITY)

**State Module Extensions:**
```c
// File: include/state.h, src/state/state_tiled.c
void state_tiled_init(state_tiled_t *state, qubit_t qubits, cplx_t **data);
void state_tiled_free(state_tiled_t *state);
dim_t state_tiled_len(qubit_t qubits);
cplx_t state_tiled_get(const state_tiled_t *state, dim_t row, dim_t col);
void state_tiled_set(state_tiled_t *state, dim_t row, dim_t col, cplx_t val);
void state_tiled_print(const state_tiled_t *state);
void state_tiled_plus(state_tiled_t *state, qubit_t qubits);
state_tiled_t state_tiled_cp(const state_tiled_t *state);
```

**Conversion Utilities:**
```c
// File: test/src/test_mhipster_tiled.c (static helpers)
void state_tiled_from_packed(state_tiled_t *tiled, const cplx_t *packed, qubit_t qubits);
void state_tiled_to_packed(const state_tiled_t *tiled, cplx_t *packed);
cplx_t* density_unpack_tiled(const state_tiled_t *tiled);
cplx_t* density_pack_tiled(const cplx_t *full, qubit_t qubits);
```

**Test Harnesses:**
```c
// File: test/src/test_mhipster_tiled.c
void testSingleQubitGateTiled(single_qubit_gate_tiled gate, const cplx_t *mat);
void testRotationGateTiled(rotation_gate_tiled gate, void (*mat_fn)(double, cplx_t*));
void testTwoQubitGateTiled(two_qubit_gate_tiled gate, const cplx_t *mat);
void testThreeQubitGateTiled(three_qubit_gate_tiled gate, const cplx_t *mat);
```

#### Three-Qubit Gate Infrastructure (MEDIUM PRIORITY)

**Reference Matrix Builder:**
```c
// File: test/src/gatemat.c
cplx_t* mat_three_qubit_gate(unsigned nqubits, const cplx_t *gate,
                              unsigned q1, unsigned q2, unsigned q3);
```

**Test Harnesses:**
```c
// File: test/src/test_qhipster.c
void testThreeQubitGate(three_qubit_gate gate, const cplx_t *mat);

// File: test/src/test_mhipster.c
void testThreeQubitGateMixed(three_qubit_gate gate, const cplx_t *mat);
```

**Typedef:**
```c
// File: test/include/test_gate.h
typedef qs_error_t (*three_qubit_gate)(state_t *state, qubit_t q1, qubit_t q2, qubit_t q3);
```

#### Missing Gate Matrices

**File: test/include/test_gate.h**

```c
// Controlled-Y: Control=1 applies Y to target
static const cplx_t CYMAT[16] = {
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 0, -I,
    0, 0, I, 0
};

// Controlled-Z: Symmetric (control and target interchangeable)
static const cplx_t CZMAT[16] = {
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, -1
};

// SWAP: Swaps states of two qubits
// Already defined as SWAPMAT ✅

// Toffoli (CCX): 8×8 matrix, flips target if both controls are 1
static const cplx_t CCXMAT[64] = {
    // Basis: |000>, |001>, |010>, |011>, |100>, |101>, |110>, |111>
    // Only |110> and |111> are swapped
    1, 0, 0, 0, 0, 0, 0, 0,  // |000> -> |000>
    0, 1, 0, 0, 0, 0, 0, 0,  // |001> -> |001>
    0, 0, 1, 0, 0, 0, 0, 0,  // |010> -> |010>
    0, 0, 0, 1, 0, 0, 0, 0,  // |011> -> |011>
    0, 0, 0, 0, 1, 0, 0, 0,  // |100> -> |100>
    0, 0, 0, 0, 0, 1, 0, 0,  // |101> -> |101>
    0, 0, 0, 0, 0, 0, 0, 1,  // |110> -> |111>
    0, 0, 0, 0, 0, 0, 1, 0   // |111> -> |110>
};
```

**Matrix builders for parameterized gates:**
```c
// File: test/src/test_qhipster.c
void mat_p(double theta, cplx_t *mat);    // Phase gate P(θ) = diag(1, e^(iθ))
void mat_cp(double theta, cplx_t *mat);   // Controlled-P: 4×4 matrix
void mat_hy(cplx_t *mat);                 // Hadamard-Y (2×2)
```

### 5.2 Test Utilities Status

| Utility | Status | Priority to Add |
|---------|--------|-----------------|
| `mat_id()` | ✅ Complete | N/A |
| `mat_single_qubit_gate()` | ✅ Complete | N/A |
| `mat_two_qubit_gate()` | ✅ Complete | N/A |
| `mat_three_qubit_gate()` | ❌ Missing | HIGH |
| `mat_rx/ry/rz()` | ✅ Complete | N/A |
| `mat_p()` | ❌ Missing | MEDIUM |
| `mat_cp()` | ❌ Missing | MEDIUM |
| `mat_hy()` | ❌ Missing | MEDIUM |
| `density_pack()` | ✅ Complete (packed) | N/A |
| `density_unpack()` | ✅ Complete (packed) | N/A |
| `density_pack_tiled()` | ❌ Missing | HIGH |
| `density_unpack_tiled()` | ❌ Missing | HIGH |
| `state_tiled_*` functions | ❌ Missing | HIGH |

---

## 6. Prioritized Implementation Plan

### Phase 0: Tiled Storage Infrastructure (CRITICAL PATH)

**Goal:** Enable testing all gates on MIXED_TILED format

**Tasks:**
1. **State module for tiled format** [HIGH]
   - Implement `state_tiled.c` and `state_tiled.h`
   - Test-first: Write `test/src/test_state_tiled.c`
   - Verify: init, free, get, set, copy, print

2. **Tiled conversion utilities** [HIGH]
   - Implement pack/unpack in `test/src/test_mhipster_tiled.c`
   - Test against known matrices (identity, single-element)

3. **Tiled test harnesses** [HIGH]
   - Duplicate harnesses from `test_mhipster.c`
   - Verify with one simple gate (Z gate - diagonal only)

**Exit Criteria:**
- State module tests pass
- Z gate passes on tiled format
- Infrastructure ready for all gates

**Estimated Test Count:** ~15 state module tests + ~6,000 Z gate tests = ~6,015 tests

### Phase 1: Complete Two-Qubit Core Gates

**Goal:** Implement CY, CZ, SWAP for all three formats

**Tasks per gate:**
1. Add gate matrix constant to `test_gate.h`
2. Add test functions to `test_gate.c`
3. Implement `cy_pure()`, `cy_mixed()`, `cy_tiled()` in gate modules
4. Add dispatcher `cy()` in `gate.c`
5. Verify all tests pass

**Gates:**
- CY (Controlled-Y)
- CZ (Controlled-Z)
- SWAP

**Exit Criteria:**
- All three gates pass ~37,700 tests × 3 formats = ~113,100 tests each
- Total: ~339,300 new tests passing

### Phase 2: Complete Remaining Single-Qubit Gates

**Goal:** Implement Hy and P(θ) for all three formats

**Tasks:**
1. Hy gate (non-parameterized)
   - Add HYMAT constant
   - Implement for all formats
   - Test count: ~11,800 × 3 = ~35,400 tests

2. P(θ) gate (parameterized)
   - Add `mat_p()` builder
   - Implement for all formats
   - Test 6 angles
   - Test count: ~70,000 × 3 = ~210,000 tests

**Exit Criteria:**
- Hy: ~35,400 tests passing
- P(θ): ~210,000 tests passing
- Total: ~245,400 new tests

### Phase 3: Parameterized Two-Qubit Gates

**Goal:** Implement CP(θ) for all three formats

**Tasks:**
1. Add `mat_cp()` matrix builder
2. Implement `cp_pure()`, `cp_mixed()`, `cp_tiled()`
3. Test 6 angles × all qubit pairs

**Exit Criteria:**
- CP(θ): ~113,000 × 3 formats = ~339,000 tests passing

### Phase 4: Remaining Controlled Gates

**Goal:** Complete controlled gate library

**Gates:**
- CS, CSdg (Controlled phase gates)
- CH, CHy (Controlled Hadamard variants)
- CT, CTdg (Controlled T gates)
- CPdg(θ) (Controlled phase dagger)
- RSWAP (Rotated swap)

**Exit Criteria:**
- 8 gates × ~37,700 tests × 3 formats = ~905,000 tests passing

### Phase 5: Three-Qubit Gates

**Goal:** Complete Toffoli gate

**Tasks:**
1. **Infrastructure:**
   - Implement `mat_three_qubit_gate()` in gatemat.c
   - Implement `testThreeQubitGate()` harnesses (pure, mixed, tiled)
   - Add `three_qubit_gate` typedef

2. **Toffoli Implementation:**
   - Add CCXMAT (8×8 matrix)
   - Implement `ccx_pure()`, `ccx_mixed()`, `ccx_tiled()`
   - Add dispatcher `ccx()`

**Exit Criteria:**
- Infrastructure tests pass
- Toffoli: ~40,000 × 3 formats = ~120,000 tests passing

### Phase 6: Final Validation

**Goal:** Complete test coverage verification

**Tasks:**
1. Run full test suite (~2M tests)
2. Verify runtime < 2 minutes
3. Coverage analysis (line, branch, condition)
4. Document any known limitations

**Exit Criteria:**
- All 93 implementations complete
- All ~2M tests passing
- Documentation updated

---

## 7. Test Count Summary

### By Gate Category

| Category | Gates | Tests per Gate (×3 formats) | Total Tests |
|----------|-------|----------------------------|-------------|
| Single-qubit (fixed) | 9 | ~35,400 | ~318,600 |
| Single-qubit (rotation) | 4 | ~105,000 | ~420,000 |
| Two-qubit (fixed) | 11 | ~56,400 | ~620,400 |
| Two-qubit (parameterized) | 2 | ~339,000 | ~678,000 |
| Three-qubit | 1 | ~60,000 | ~60,000 |
| **GRAND TOTAL** | **27** | - | **~2,097,000** |

### By Implementation Status

| Status | Implementations | Tests | Percentage |
|--------|----------------|-------|------------|
| Complete | 24/93 | ~450,000/2.1M | 26%/21% |
| Missing | 69/93 | ~1.65M/2.1M | 74%/79% |

### By Storage Format

| Format | Complete | Missing | Total Needed |
|--------|----------|---------|--------------|
| PURE | 12/31 | 19/31 | 31 |
| MIXED_PACKED | 12/31 | 19/31 | 31 |
| MIXED_TILED | 0/31 | 31/31 | 31 |
| **TOTAL** | **24/93** | **69/93** | **93** |

---

## 8. Quality Assurance Checklist

### Per-Gate Validation

For each gate implementation, verify:

- [ ] **Correctness:** All test states match reference within ε=1e-12
- [ ] **Norm preservation (pure):** |⟨ψ'|ψ'⟩ - 1| < ε for all test vectors
- [ ] **Trace preservation (mixed):** |Tr(ρ') - Tr(ρ)| < ε for all test matrices
- [ ] **Hermiticity (mixed):** ρ'[i,j] = conj(ρ'[j,i]) for all i,j
- [ ] **Error handling:** NULL checks, range checks return appropriate error codes
- [ ] **All positions:** Gate tested at all valid qubit positions
- [ ] **All pairs/triples:** Multi-qubit gates tested for all valid combinations
- [ ] **Entangled states:** Bell/GHZ/W states tested (multi-qubit gates)
- [ ] **Edge angles:** Rotation gates tested at θ=0, π/2, π, -π/3

### Global Validation

After completing all gates:

- [ ] **Build time:** Full test suite runs in < 2 minutes
- [ ] **Memory leaks:** Valgrind clean on representative subset
- [ ] **Thread safety:** No data races under ThreadSanitizer
- [ ] **Numerical stability:** No NaN or Inf in any test
- [ ] **Documentation:** All functions have @brief, @param, @note tags
- [ ] **Code review:** All implementations reviewed for correctness
- [ ] **Regression testing:** No existing tests broken

---

## 9. Risk Assessment

### High Risk Items

1. **Tiled format complexity**
   - **Risk:** Tile indexing bugs difficult to debug
   - **Mitigation:** Test-first approach, start with single tile (n < 5)
   - **Validation:** Compare tiled vs packed output for all gates

2. **Test count explosion**
   - **Risk:** 2M tests may exceed 2-minute runtime budget
   - **Mitigation:** Profile and optimize test harnesses, parallelize test execution
   - **Fallback:** Reduce MAXQUBITS from 4 to 3 (cuts tests by ~40%)

3. **Three-qubit gate correctness**
   - **Risk:** 8×8 matrix more complex, harder to verify
   - **Mitigation:** Cross-validate with external simulator (Qiskit/Cirq)
   - **Validation:** Hand-verify Toffoli truth table for n=3 case

### Medium Risk Items

1. **Controlled gate symmetry assumptions**
   - **Risk:** Assuming CZ is symmetric may hide bugs
   - **Mitigation:** Test both orderings explicitly even for symmetric gates

2. **Parameterized gate angle coverage**
   - **Risk:** 6 test angles may miss edge cases
   - **Mitigation:** Current angles well-chosen, covers identity/π/negative

3. **Random mixture reproducibility**
   - **Risk:** Non-deterministic test failures
   - **Mitigation:** Fixed random seed in test infrastructure ✅

### Low Risk Items

1. **Numerical precision edge cases**
   - **Risk:** ε=1e-12 may be too tight for large n
   - **Mitigation:** Empirically validated for n≤4, can adjust if needed

2. **Memory overhead for large tests**
   - **Risk:** n=4 density matrices are 16×16 (4KB each)
   - **Mitigation:** Acceptable for test environment, not production

---

## 10. Recommendations

### Immediate Actions

1. **Implement tiled storage infrastructure (Phase 0)**
   - Highest priority: blocks all tiled gate testing
   - Estimated effort: 2-3 days
   - Deliverable: Z gate passing on tiled format

2. **Complete CY, CZ, SWAP gates (Phase 1)**
   - High priority: fundamental two-qubit operations
   - Estimated effort: 1 day per gate × 3 gates = 3 days
   - Deliverable: Core two-qubit library complete

3. **Add three-qubit infrastructure (Phase 5 prep)**
   - Medium priority: needed for Toffoli
   - Estimated effort: 1 day
   - Deliverable: `testThreeQubitGate()` harness ready

### Long-Term Strategy

1. **Parallelize test execution**
   - Use CTest parallel execution: `ctest -j8`
   - Expected speedup: 4-8× on multi-core systems
   - Keeps runtime under 2-minute budget

2. **Benchmark-driven optimization**
   - Profile gate implementations after correctness established
   - Focus on H gate (full mixing, most expensive)
   - Compare tiled vs packed performance on n=8-14 systems

3. **Extended validation**
   - Cross-validate with Qiskit/Cirq for n=5,6 systems
   - Test against QuEST library (external/QuEST)
   - Publish verification results

### Optional Enhancements

1. **Property-based testing**
   - Verify gate identities: X²=I, H²=I, CNOT²=I
   - Verify commutation relations: XZ=-ZX
   - Verify composition: HXH=Z

2. **Performance regression tests**
   - Track gate throughput (gates/second) over time
   - Alert on >10% performance degradation
   - Compare tiled vs packed benchmarks

3. **Extended system sizes**
   - Add n=5,6 tests for critical gates
   - Requires larger reference matrices (32×32, 64×64)
   - Validates stride indexing at larger scales

---

## 11. Conclusion

The quantum gate module has solid test infrastructure for PURE and MIXED_PACKED formats, with 11/31 gates fully tested. The critical path forward is:

1. **Phase 0:** Implement tiled storage infrastructure (~6K tests)
2. **Phase 1:** Complete core two-qubit gates (~339K tests)
3. **Phase 2-5:** Systematically implement remaining gates (~1.3M tests)

Total test coverage will reach **~2.1M tests across 93 gate implementations** (31 gates × 3 formats), providing comprehensive validation of the gate library. With parallel test execution, the 2-minute runtime budget is achievable.

**Key Success Metrics:**
- ✅ All 93 implementations complete
- ✅ All ~2.1M tests passing
- ✅ Runtime < 2 minutes
- ✅ Zero memory leaks
- ✅ Documentation complete

The test strategy ensures correctness through exhaustive coverage while remaining practical for continuous integration.
