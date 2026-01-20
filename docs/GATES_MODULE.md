# Gates Module: Technical Specification

**Module:** Quantum Gate Operations
**Header:** `include/gate.h`
**Implementation:** `src/qhipster.c`, `src/mhipster.c`, `src/gate.c`
**Last Updated:** 2026-01-20

---

## Overview

This module implements quantum gate operations for both pure states (state vectors) and mixed states
(density matrices). Gates exploit sparse, local structure of quantum operations—operating directly on
state arrays via position-dependent stride indexing rather than constructing full gate matrices via
Kronecker products.

**Key design characteristics:**
- **Position-dependent implementation**: Stride between elements depends on target qubit index—different
  positions exercise different code paths
- **Unified API**: Single `x(state, target)` function dispatches to pure or mixed implementation
  based on `state->type` enum
- **Direct packed operations**: Mixed state gates operate directly on lower-triangular packed storage
  (no unpacking required)
- **Macro-based traversal** (mixed states): Common loop structure factored into `TRAVERSE_MIXED_1Q`
  macro with gate-specific operations passed as macro parameters—enables adding new gates with minimal
  code while ensuring consistent traversal logic
- **Test-driven validation**: Reference implementations using full BLAS matrix operations verify
  correctness

---

## Unit Tests

### Verification Approach

Gate functions validated against reference implementations:
- **Pure states**: $ \lvert \psi^{\prime} \rangle = U \lvert \psi \rangle $ via matrix-vector
  multiplication (`zgemv()`)
- **Mixed states**: $ \rho^{\prime} = U \rho U^{\dagger} $ via matrix-matrix multiplication (`zgemm()`)

**Reference vs Production:**
- **Reference**: Full matrix representation with hardcoded Kronecker products (I⊗...⊗U⊗...⊗I)
  - Simple, verifiable correctness
  - Limited to n≤4 qubits (16×16 matrices, ~4KB memory)
  - Used only in test code

- **Production**: Sparse stride-based indexing on packed arrays
  - Exploits locality: single-qubit gate affects 2^(n-1) amplitude pairs
  - Scales to n=32 pure / n=16 mixed
  - Position-dependent logic validated by testing all qubit positions

**Separation of concerns:**
- Tests verify correctness via simple, trusted BLAS operations
- Production code optimized for performance and scalability
- Pack/unpack operations isolated to test infrastructure

### Test State Basis

**Single-qubit informationally complete set (6 states)**:
- Computational: $ \lvert 0 \rangle, \ \lvert 1 \rangle $
- Hadamard: $ \lvert + \rangle, \ \lvert - \rangle $
- Circular: $ \lvert +i \rangle, \ \lvert -i \rangle $

**Rationale**: Detects amplitude errors (computational), real phase errors (Hadamard), complex phase
errors (circular).

**Multi-qubit pure states**:
- Tensor products of single-qubit bases
- Bell states (n=2): $ \frac{1}{\sqrt{2}} (\lvert 00 \rangle \pm \lvert 11 \rangle) $,
  $ \frac{1}{\sqrt{2}} (\lvert 01 \rangle \pm \lvert 10 \rangle) $
- GHZ states (n≥3): $ \frac{1}{\sqrt{2}} (\lvert 0 \dots 0 \rangle + \lvert 1 \dots 1 \rangle) $
- W states (n≥3): $ \frac{1}{\sqrt{n}}(\lvert 100\dots0 \rangle + \lvert 010\dots0 \rangle + \dots +
  \lvert 000\dots1\rangle) $

**Rationale**: Bell/GHZ/W test entanglement propagation and multi-partite correlations critical for
two-qubit gate validation.

**Mixed states**:
- **Pure-as-mixed**: All pure states as density matrices $ \lvert \psi \rangle \langle \psi \rvert $
- **Maximally mixed**: $ I/2^{n} $
- **Random statistical mixtures**: $ \rho = \sum_{i} p_i \lvert \psi_i \rangle \langle \psi_i \rvert $
  - Random probabilities $ p_i \geq 0 $, $ \sum_i p_i = 1 $
  - States $ \lvert \psi_i \rangle $ randomly selected from pure test states (computational, Hadamard,
    circular, Bell, GHZ, W)
  - Variable rank (typically 2-10 components)
  - 10-20 random mixtures per system size
  - Deterministic random seed for reproducibility

**Rationale**: Rank>1 states verify conjugation $ U \rho U^{\dagger} $ operates on full density matrix,
not just extracted pure components. Random mixtures ensure broad coverage of density matrix space
including arbitrary eigenvalue distributions.

### Test Coverage Strategy

**Why test all positions:**
qHiPSTER-based implementation uses position-dependent stride indexing:
- Different target qubits → different memory access patterns
- Each position exercises distinct code paths in sparse operations
- Testing all positions validates stride calculation logic

**Coverage per gate type:**
- **Single-qubit gates**: All target positions (0 to n-1) on n=1,2,3,4 qubit systems
- **Two-qubit gates**: All (target, control) pairs on n=2,3,4 qubit systems
- **Three-qubit gates**: All (target, control1, control2) triples on n=3,4 qubit systems

**Why n≤4:**
- Reference implementations limited by Kronecker product memory (16×16 matrices manageable)
- Validates position-dependent logic across representative system sizes
- Trust that correct stride indexing scales to n=32

### Test Count Estimate

| Gate Type | Pure States | Mixed States | Total per Gate |
|-----------|-------------|--------------|----------------|
| Single-qubit | ~5,900 | ~5,900 | ~11,800 |
| Two-qubit | ~18,800 | ~18,900 | ~37,700 |
| Three-qubit | ~20,000 | ~20,000 | ~40,000 |

**Library totals**:
- ~30 single-qubit gates: 355k tests
- ~15 two-qubit gates: 566k tests
- ~2 three-qubit gates: 80k tests
- **Grand total: ~1M tests**

**Runtime requirement**: Complete all gate tests in <2 minutes (build-time regression testing)

**Empirical performance**: Tests for n=1,...,6 with full basis coverage execute in milliseconds per
gate.

---

## Test Infrastructure

### Current Status

**Phase 0: Testing Infrastructure - COMPLETE ✅**

All infrastructure for testing gates on both pure and mixed states is complete:
- **Pure state generators**: Computational, Hadamard, circular bases
- **Mixed state generators**: Pure-as-mixed, maximally mixed, random mixtures
- **Entangled state generators**: Bell, GHZ, W states (both pure and mixed)
- **Pure state test harness**: `testSingleQubitGate()` in `test/src/test_qhipster.c`
- **Mixed state test harness**: `testSingleQubitGateMixed()` in `test/src/test_mhipster.c`
- **Packed storage utilities**: `density_pack()` and `density_unpack()` (static in test_mhipster.c)
- **BLAS helper**: `zumu()` computes U*M*U† for reference comparisons
- **Unified test file**: `test/src/test_gate.c` with single main() for all gate tests

**Test harness features:**
- Function pointer + 2×2 matrix → automated testing
- Tests all positions on n=1,2,3,4 qubit systems
- Builds full unitary via Kronecker products
- Validates outputs against BLAS reference implementations

**Phase 1: First Single-Qubit Gate + Error Handling - COMPLETE ✅**

Pauli-X gate fully implemented and tested for both pure and mixed states:
- ✅ `x_pure()` in `src/qhipster.c` - Pure state implementation
- ✅ `x_mixed()` in `src/mhipster.c` - Mixed state implementation with packed storage
- ✅ `x()` dispatcher in `src/gate.c` with input validation (NULL checks, range checks)
- ✅ Function documentation with @brief, @param, and @note tags
- ✅ Comprehensive inline documentation of packed storage algorithm
- ✅ All tests passing (5/5 tests in test_gate)

**Error Handling Infrastructure - COMPLETE ✅**
- ✅ `qs_error_t` enum added to `q_types.h` with standard error codes
- ✅ Gate functions return `qs_error_t` instead of `void`
- ✅ Error handling tests: `test_error_null_state`, `test_error_null_data`, `test_error_qubit_out_of_range`
- ✅ Test harnesses updated to verify `QS_OK` return codes

**Next: Phase 2 - Implement remaining single-qubit gates (Y, Z, H, S, T, rotations)**

### State Generators

**Pure state generators** (location: `test/src/test_pure_states.c`):

| Function | Generates | Status |
|----------|-----------|--------|
| `test_cb_pure()` | Computational basis states | ✅ Complete |
| `test_xb_pure()` | Hadamard basis states | ✅ Complete |
| `test_yb_pure()` | Circular basis states | ✅ Complete |
| `test_bell_pure()` | 4 Bell states (n=2) | ✅ Complete |
| `test_ghz_pure()` | n-qubit GHZ state (n≥2) | ✅ Complete |
| `test_w_pure()` | n-qubit W state (n≥3) | ✅ Complete |
| `test_mk_states_pure()` | All pure test states | ✅ Complete |
| `test_rm_states_pure()` | Cleanup: free pure states | ✅ Complete |

**Mixed state generators** (location: `test/src/test_mixed_states.c`):

| Function | Generates | Status |
|----------|-----------|--------|
| `test_cb_mixed()` | Computational basis as density matrices | ✅ Complete |
| `test_xb_mixed()` | Hadamard basis as density matrices | ✅ Complete |
| `test_yb_mixed()` | Circular basis as density matrices | ✅ Complete |
| `test_bell_mixed()` | 4 Bell states as density matrices (n=2) | ✅ Complete |
| `test_ghz_mixed()` | GHZ state as density matrix (n≥2) | ✅ Complete |
| `test_w_mixed()` | W state as density matrix (n≥3) | ✅ Complete |
| `test_maximally_mixed()` | Maximally mixed state I/2ⁿ | ✅ Complete |
| `test_random_mixture()` | Random mixture ρ=Σpᵢ\|ψᵢ⟩⟨ψᵢ\| | ✅ Complete |
| `test_mk_states_mixed()` | All mixed test states | ✅ Complete |
| `test_rm_states_mixed()` | Cleanup: free mixed states | ✅ Complete |

### Test Harness Functions

| Function | Purpose | Location | Status |
|----------|---------|----------|--------|
| `testSingleQubitGate()` | Validates single-qubit gates on pure states | `test/src/test_qhipster.c` | ✅ Complete |
| `testSingleQubitGateMixed()` | Validates single-qubit gates on mixed states | `test/src/test_mhipster.c` | ✅ Complete |
| `testTwoQubitGate()` | Validates two-qubit gates on pure states | TBD | 🔲 Phase 4 |
| `testTwoQubitGateMixed()` | Validates two-qubit gates on mixed states | TBD | 🔲 Phase 4 |
| `testThreeQubitGate()` | Validates three-qubit gates on pure states | TBD | 🔲 Phase 5 |
| `testThreeQubitGateMixed()` | Validates three-qubit gates on mixed states | TBD | 🔲 Phase 5 |

### Test Helper Functions

| Function | Purpose | Location             | Status |
|----------|---------|----------------------|--------|
| `zmv()` | Matrix-vector multiplication M*v via `zgemv()` | `test/src/linalg.c` | ✅ Complete |
| `zumu()` | Similarity transformation U*M*U† via `zgemm()` | `test/src/linalg.c` | ✅ Complete |
| `zkron()` | Kronecker product A⊗B | `test/src/linalg.c`  | ✅ Complete |
| `mat_id()` | Identity matrix I(2ⁿ×2ⁿ) | `test/src/gatemat.c` | ✅ Complete |
| `mat_single_qubit_gate()` | Full matrix I⊗...⊗U⊗...⊗I for single-qubit gate | `test/src/gatemat.c` | ✅ Complete |
| `density_pack()` | Full Hermitian matrix → packed storage (static) | `test/src/test_mhipster.c` | ✅ Complete |
| `density_unpack()` | Packed lower-triangle → full Hermitian matrix (static) | `test/src/test_mhipster.c` | ✅ Complete |
| `mat_two_qubit_gate()` | Full matrix for two-qubit gate | `test/src/gatemat.c` | 🔲 Phase 4 |
| `mat_three_qubit_gate()` | Full matrix for three-qubit gate | `test/src/gatemat.c` | 🔲 Phase 5 |

### Test File Organization

**Unified header approach** minimizes file overhead while maintaining clear separation of concerns:

**Headers:**
- `test/include/test_gate.h` - Single unified header containing:
  - `single_qubit_gate` typedef: `qs_error_t (*)(state_t*, qubit_t)` (returns error code)
  - Gate matrix constants (XMAT, YMAT, ZMAT, HMAT, SMAT, TMAT, P0MAT, P1MAT, SWAPMAT)
  - Declarations for `testSingleQubitGate()` and `testSingleQubitGateMixed()`

**Source files:**
- `test/src/test_qhipster.c` - Pure state test harness implementation only
  - Contains `testSingleQubitGate()` function
  - No test functions or main()
- `test/src/test_mhipster.c` - Mixed state test harness implementation only
  - Contains `testSingleQubitGateMixed()` function
  - Contains static `density_pack()` and `density_unpack()` helpers
  - No test functions or main()
- `test/src/test_gates.c` - **All gate test definitions**
  - Contains `setUp()`, `tearDown()`, and `main()`
  - Contains all test functions: `test_X_pure()`, `test_X_mixed()`, etc.
  - Single executable for all gate tests (pure + mixed)

**Benefits:**
- Single source of truth for typedef and gate matrices
- Infrastructure (harnesses) separated from test definitions
- Easy to add new gates: just add test functions to test_gates.c
- Minimal compilation overhead: one test executable

---

## Production Code Implementation Status

**Location**: `src/qhipster.c`, `src/mhipster.c`, `src/gate.c`, `include/gate.h`

### Single-Qubit Gate Implementations

| Gate     | Description         | Pure State | Mixed State | Phase |
|----------|---------------------|------------|-------------|-------|
| `x()`    | Pauli-X             | ✅ Complete | ✅ Complete | 1 |
| `y()`    | Pauli-Y             | ✅ Complete | ✅ Complete | 2 |
| `z()`    | Pauli-Z             | ✅ Complete | ✅ Complete | 2 |
| `h()`    | Hadamard            | 🔲 TODO | 🔲 TODO | 2 |
| `s()`    | Phase gate          | 🔲 TODO | 🔲 TODO | 2 |
| `sdg()`  | S†                  | 🔲 TODO | 🔲 TODO | 2 |
| `t()`    | π/8 gate            | 🔲 TODO | 🔲 TODO | 2 |
| `tdg()`  | T†                  | 🔲 TODO | 🔲 TODO | 2 |
| `hy()`   | Hadamard-Y          | 🔲 TODO | 🔲 TODO | 2 |
| `p()`    | Phase rotation P(θ) | 🔲 TODO | 🔲 TODO | 2 |
| `pdg()`  | P†(θ)               | 🔲 TODO | 🔲 TODO | 2 |
| `rx()`   | X-rotation RX(θ)    | 🔲 TODO | 🔲 TODO | 2 |
| `rxdg()` | RX†(θ)              | 🔲 TODO | 🔲 TODO | 2 |
| `ry()`   | Y-rotation RY(θ)    | 🔲 TODO | 🔲 TODO | 2 |
| `rydg()` | RY†(θ)              | 🔲 TODO | 🔲 TODO | 2 |
| `rz()`   | Z-rotation RZ(θ)    | 🔲 TODO | 🔲 TODO | 2 |
| `rzdg()` | RZ†(θ)              | 🔲 TODO | 🔲 TODO | 2 |

### Two-Qubit Gate Implementations

| Gate | Description | Pure State | Mixed State | Phase |
|------|-------------|------------|-------------|-------|
| `cx()` | CNOT | 🔲 TODO | 🔲 TODO | 4 |
| `cy()` | Controlled-Y | 🔲 TODO | 🔲 TODO | 4 |
| `cz()` | Controlled-Z | 🔲 TODO | 🔲 TODO | 4 |
| `cs()` | Controlled-S | 🔲 TODO | 🔲 TODO | 4 |
| `csdg()` | Controlled-S† | 🔲 TODO | 🔲 TODO | 4 |
| `ch()` | Controlled-H | 🔲 TODO | 🔲 TODO | 4 |
| `chy()` | Controlled-Hy | 🔲 TODO | 🔲 TODO | 4 |
| `ct()` | Controlled-T | 🔲 TODO | 🔲 TODO | 4 |
| `ctdg()` | Controlled-T† | 🔲 TODO | 🔲 TODO | 4 |
| `cp()` | Controlled-P(θ) | 🔲 TODO | 🔲 TODO | 4 |
| `cpdg()` | Controlled-P†(θ) | 🔲 TODO | 🔲 TODO | 4 |
| `swap()` | SWAP | 🔲 TODO | 🔲 TODO | 4 |
| `rswap()` | Rotated SWAP | 🔲 TODO | 🔲 TODO | 4 |

### Three-Qubit Gate Implementations

| Gate | Description | Pure State | Mixed State | Phase |
|------|-------------|------------|-------------|-------|
| `ccx()` | Toffoli (CCNOT) | 🔲 TODO | 🔲 TODO | 5 |

---

## Design Notes

### Position-Dependent Stride Logic

qHiPSTER uses stride-based indexing where stride = 2^target. Different target qubits exercise different memory access patterns and code paths, requiring validation of all qubit positions.

### Mixed State Traversal Macro (`TRAVERSE_MIXED_1Q`)

Single-qubit gates on density matrices share identical traversal structure but differ in element operations.
The `TRAVERSE_MIXED_1Q` macro in `src/mhipster.c` factors out the common loop structure:

```c
TRAVERSE_MIXED_1Q(state, target, DIAG_OP, CONJ_OP, OFFDIAG_OP)
```

**Operation points** (gate-specific macros called at each point):
- `DIAG_OP(data, n, stride)`: Process (0,0)↔(1,1) quadrant pairs (diagonal blocks)
- `CONJ_OP(d10, count, start_k, dim, col_block)`: Process (1,0)↔(0,1) pairs where (0,1) is in
  upper triangle (accessed via conjugate)
- `OFFDIAG_OP(d10, d01, n)`: Process (1,0)↔(0,1) pairs where both are in lower triangle

**Example gate definition** (complete implementation):
```c
#define X_DIAG_OP(data, n, stride)      op_swap_stride(data, n, stride)
#define X_CONJ_OP(d10, count, ...)      op_swap_conj(d10, count, __VA_ARGS__)
#define X_OFFDIAG_OP(d10, d01, n)       op_swap(d10, d01, n)

void x_mixed(state_t *state, const qubit_t target) {
    TRAVERSE_MIXED_1Q(state, target, X_DIAG_OP, X_CONJ_OP, X_OFFDIAG_OP);
}
```

**Benefits**:
- Adding new single-qubit gates requires only defining 3 operation macros + 1 function
- Bug fixes to traversal logic automatically apply to all gates
- `static inline` helper functions ensure no runtime overhead from abstraction

### Numerical Precision

**Test tolerance**: ε = 1e-12 accounts for double-precision arithmetic plus accumulated rounding errors over O(2^n) sparse operations. Sparse operations provide better numerical stability than dense matrix multiplication (O(2^3n) operations).

### Packed Storage for Mixed States

**Production**: Gates operate directly on lower-triangular packed format achieving 50% memory savings (critical for n=16: 32GB vs 64GB).

**Testing**: Reference implementation unpacks → applies UρU† via BLAS → repacks, isolating pack/unpack logic from gate logic.

### Unified State Representation

Single `x(state, target)` API dispatches to `x_pure()` or `x_mixed()` based on `state->type` enum. Type-safe at runtime with optimized per-type implementations.

### Error Handling

All gate functions return `qs_error_t` (defined in `q_types.h`):
- `QS_OK` (0): Success
- `QS_ERR_NULL` (-1): Null state or null data pointer
- `QS_ERR_QUBIT` (-3): Target qubit index out of range

Internal implementations (`x_pure`, `x_mixed`, etc.) trust validated inputs and do not perform additional checks.

### Thread Safety

State objects not thread-safe for concurrent mutation. Users can create separate state objects for parallel evaluation. Gate operations are memory-bandwidth limited (~3× speedup on 10 cores with OpenMP).

---

## Implementation Workflow

**Iterative development strategy**: Build one gate end-to-end (pure + mixed, tested) before moving to next gate. This validates infrastructure early and establishes clear patterns.

### Phase 1: Complete First Single-Qubit Gate (Pure + Mixed) - COMPLETE ✅
**Goal**: First fully validated gate demonstrating complete workflow

**Status**: COMPLETE ✅

**Deliverables**:
- ✅ `density_pack()` and `density_unpack()` helper functions (static in test_mhipster.c)
- ✅ `testSingleQubitGateMixed()` test harness
- ✅ `zumu()` BLAS helper for U*M*U† computation
- ✅ Pure state X gate implementation in `src/qhipster.c`
- ✅ Mixed state X gate implementation in `src/mhipster.c`
  - Direct operation on packed density matrix using qHiPSTER stride logic
  - Handles both pure and mixed dispatch via `state->type` enum in `src/gate.c`
- ✅ Input validation in `x()` dispatcher (NULL checks, range checks)
- ✅ Function-level and inline documentation
- ✅ All tests passing (test_X_pure and test_X_mixed)

**Exit criteria met**: Both pure and mixed X implementations pass all tests, establishing clear pattern for remaining gates

**Key achievements**:
- Established pattern for unified gate API with type-based dispatch
- Validated packed storage approach for mixed states
- Demonstrated correct error handling without corrupting caller's state
- Proved qHiPSTER stride indexing works correctly for packed lower-triangular storage

### Phase 2: Remaining Single-Qubit Gates - IN PROGRESS 🔄
**Goal**: Complete single-qubit gate library

**Status**: Pauli gates complete, Clifford/rotation gates pending
- ✅ Y gate: Complete (pure + mixed), all tests passing
- ✅ Z gate: Complete (pure + mixed), all tests passing
- ✅ Mixed state traversal refactored: `TRAVERSE_MIXED_1Q` macro extracts common loop structure
- 🔲 Remaining: H, S, S†, T, T†, Hy, rotation gates

**Deliverables**: Implement all remaining single-qubit gates for both pure and mixed states following Phase 1 pattern:
- Pauli gates: Y, Z
- Clifford gates: H (Hadamard), S, S†
- Non-Clifford: T, T†, Hy (Hadamard-Y)
- Rotation gates: RX(θ), RY(θ), RZ(θ), P(θ) and their daggers

**Implementation pattern per gate** (established from X gate):
1. Add tests in `test/src/test_gate.c` (TDD approach)
2. Add declaration in `include/gate.h` with documentation
3. Add dispatcher in `src/gate.c` with input validation (returns `qs_error_t`)
4. Implement `{gate}_pure()` in `src/qhipster.c`
5. Implement `{gate}_mixed()` in `src/mhipster.c` with packed storage
6. Verify all tests pass

**Exit criteria**: All single-qubit gates pass tests (pure + mixed) - ~30 gates × ~11,800 tests/gate = ~355k tests total

### Phase 3: Two-Qubit Gate Test Infrastructure
**Goal**: Enable testing of two-qubit gates

**Deliverables**:
- `mat_two_qubit_gate()`: Build full matrix via Kronecker products (4×4 → 2^n×2^n)
- `testTwoQubitGate()` and `testTwoQubitGateMixed()`: Test harnesses with (target, control) pair loops

**Exit criteria**: Infrastructure ready for two-qubit gate implementation

**Note**: Entangled state generators (Bell, GHZ, W) already complete from Phase 0

### Phase 4: Two-Qubit Gates
**Goal**: Complete two-qubit gate library

**Deliverables**: Implement all two-qubit gates (CNOT, CZ, SWAP, controlled operations) for both pure and mixed states

**Exit criteria**: All two-qubit gates pass tests (pure + mixed)

### Phase 5: Three-Qubit Gates
**Goal**: Complete Toffoli gate

**Deliverables**:
- `mat_three_qubit_gate()` and test harnesses
- Toffoli implementation (pure + mixed)

**Exit criteria**: Toffoli passes all tests - complete gate library with comprehensive test coverage

---

## Quality Assurance

### Build-Time Testing
- **Runtime**: All tests complete in <2 minutes
- **Coverage**: All gates tested across all positions, system sizes (n=1-4), state types, and basis states
- **CI/CD**: Tests run on every build, abort on failure

### Validation Checks
- Output matches reference within ε = 1e-12
- Norm preservation (pure): |⟨ψ'|ψ'⟩ - 1| < ε
- Trace preservation (mixed): |Tr(ρ') - Tr(ρ)| < ε
- Memory leak testing via valgrind/Instruments during development

---

## Future Considerations

**GPU acceleration**: Gate operations are memory-bandwidth limited—potential for GPU speedup in
streaming operations

**Extended gate set**: Additional gates based on algorithm requirements (e.g., multi-controlled gates,
custom parameterized gates)

**Performance optimization**: Profile-guided optimization after correctness established

**Noise models**: Post-v1.0 extension for realistic hardware simulation
