# Gates Module: Technical Specification

**Module:** Quantum Gate Operations
**Header:** `include/gate.h`
**Implementation:** `src/qhipster.c`
**Last Updated:** 2025-12-24

---

## Overview

This module implements quantum gate operations for both pure states (state vectors) and mixed states
(density matrices). Gates exploit sparse, local structure of quantum operations—operating directly on
state arrays via position-dependent stride indexing rather than constructing full gate matrices via
Kronecker products.

**Key design characteristics:**
- **Position-dependent implementation**: Stride between elements depends on target qubit index—different
  positions exercise different code paths
- **Unified API**: Single `applyX(state, target)` function dispatches to pure or mixed implementation
  based on `state->type` enum
- **Direct packed operations**: Mixed state gates operate directly on lower-triangular packed storage
  (no unpacking required)
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

**Single-qubit pure state testing: COMPLETE ✅**

The `testSingleQubitGate()` harness provides complete infrastructure:
- Function pointer + 2×2 matrix → automated testing
- Loops through n=1,...,4 qubits
- Tests all target positions per system size
- Builds full unitary via Kronecker products
- Generates test vectors (passed to qHiPSTER function)
- Generates reference vectors (multiplied by full matrix)
- Compares outputs with numerical tolerance

**Example usage:**
```c
testSingleQubitGate(applyX, XMAT, num_qubits, target);
```

This architecture is clean, lean, and proven. Mixed state and multi-qubit testing will follow the
same pattern.

### State Generators

| Function | Generates | Location | Status |
|----------|-----------|----------|--------|
| `test_cb_pure()` | Computational basis states | `test.c` | ✅ Complete |
| `test_xb_pure()` | Hadamard basis states | `test.c` | ✅ Complete |
| `test_yb_pure()` | Circular basis states | `test.c` | ✅ Complete |
| `test_gen_states_pure()` | Combined: all three bases | `test.c` | ✅ Complete |
| `test_rm_states_pure()` | Cleanup: free test states | `test.c` | ✅ Complete |
| `test_pure_as_mixed()` | Convert pure states → density matrices | `test.c` | 🔲 Phase 0 |
| `test_maximally_mixed()` | I/2ⁿ | `test.c` | 🔲 Phase 0 |
| `test_random_mixture()` | Random mixtures ρ=Σpᵢ\|ψᵢ⟩⟨ψᵢ\| | `test.c` | 🔲 Phase 0 |
| `test_bell_states()` | 4 Bell states | `test.c` | 🔲 Phase 3 |
| `test_ghz_state()` | n-qubit GHZ | `test.c` | 🔲 Phase 3 |
| `test_w_state()` | n-qubit W | `test.c` | 🔲 Phase 3 |

### Test Harness Functions

| Function | Purpose | Location | Status |
|----------|---------|----------|--------|
| `testSingleQubitGate()` | Validates single-qubit gates on pure states | `test_qhipster.c` | ✅ Complete |
| `testSingleQubitGateMixed()` | Validates single-qubit gates on mixed states | `test_qhipster.c` | 🔲 Phase 0 |
| `testTwoQubitGate()` | Validates two-qubit gates on pure states | `test_qhipster.c` | 🔲 Phase 4 |
| `testTwoQubitGateMixed()` | Validates two-qubit gates on mixed states | `test_qhipster.c` | 🔲 Phase 4 |
| `testThreeQubitGate()` | Validates three-qubit gates on pure states | `test_qhipster.c` | 🔲 Phase 5 |
| `testThreeQubitGateMixed()` | Validates three-qubit gates on mixed states | `test_qhipster.c` | 🔲 Phase 5 |

### Test Helper Functions

| Function | Purpose | Location | Status |
|----------|---------|----------|--------|
| `mv()` | Matrix-vector multiplication via `zgemv()` | `test_qhipster.c` | ✅ Complete |
| `kron()` | Kronecker product A⊗B | `test_qhipster.c` | ✅ Complete |
| `mat_id()` | Identity matrix I(2ⁿ×2ⁿ) | `test_qhipster.c` | ✅ Complete |
| `mat_single_qubit_gate()` | Full matrix I⊗...⊗U⊗...⊗I for single-qubit gate | `test_qhipster.c` | ✅ Complete |
| `density_pack()` | Full Hermitian matrix → packed storage | `test_qhipster.c` | 🔲 Phase 0 |
| `density_unpack()` | Packed lower-triangle → full Hermitian matrix | `test_qhipster.c` | 🔲 Phase 0 |
| `mat_two_qubit_gate()` | Full matrix for two-qubit gate | `test_qhipster.c` | 🔲 Phase 4 |
| `mat_three_qubit_gate()` | Full matrix for three-qubit gate | `test_qhipster.c` | 🔲 Phase 5 |

**Note**: Test files include predefined gate matrices (XMAT, YMAT, ZMAT, HMAT, SMAT, TMAT, SWAPMAT,
etc.) in `test/include/test_qhipster.h`.

---

## Production Code Implementation Status

**Location**: `src/qhipster.c`, `include/gate.h`

### Single-Qubit Gate Implementations

| Gate | Description | Pure State | Mixed State | Phase |
|------|-------------|------------|-------------|-------|
| `applyX()` | Pauli-X | ✅ Complete | 🔲 TODO | 1 |
| `applyY()` | Pauli-Y | ✅ Complete | 🔲 TODO | 1 |
| `applyZ()` | Pauli-Z | ✅ Complete | 🔲 TODO | 1 |
| `applyH()` | Hadamard | 🔲 TODO | 🔲 TODO | 2 |
| `applyS()` | Phase gate | 🔲 TODO | 🔲 TODO | 2 |
| `applySdagger()` | S† | 🔲 TODO | 🔲 TODO | 2 |
| `applyT()` | π/8 gate | 🔲 TODO | 🔲 TODO | 2 |
| `applyTdagger()` | T† | 🔲 TODO | 🔲 TODO | 2 |
| `applyHy()` | Hadamard-Y | 🔲 TODO | 🔲 TODO | 2 |
| `applyP()` | Phase rotation P(θ) | 🔲 TODO | 🔲 TODO | 2 |
| `applyPdagger()` | P†(θ) | 🔲 TODO | 🔲 TODO | 2 |
| `applyRX()` | X-rotation RX(θ) | 🔲 TODO | 🔲 TODO | 2 |
| `applyRXdagger()` | RX†(θ) | 🔲 TODO | 🔲 TODO | 2 |
| `applyRY()` | Y-rotation RY(θ) | 🔲 TODO | 🔲 TODO | 2 |
| `applyRYdagger()` | RY†(θ) | 🔲 TODO | 🔲 TODO | 2 |
| `applyRZ()` | Z-rotation RZ(θ) | 🔲 TODO | 🔲 TODO | 2 |
| `applyRZdagger()` | RZ†(θ) | 🔲 TODO | 🔲 TODO | 2 |

### Two-Qubit Gate Implementations

| Gate | Description | Pure State | Mixed State | Phase |
|------|-------------|------------|-------------|-------|
| `applyCX()` | CNOT | 🔲 TODO | 🔲 TODO | 4 |
| `applyCY()` | Controlled-Y | 🔲 TODO | 🔲 TODO | 4 |
| `applyCZ()` | Controlled-Z | 🔲 TODO | 🔲 TODO | 4 |
| `applyCS()` | Controlled-S | 🔲 TODO | 🔲 TODO | 4 |
| `applyCSdagger()` | Controlled-S† | 🔲 TODO | 🔲 TODO | 4 |
| `applyCH()` | Controlled-H | 🔲 TODO | 🔲 TODO | 4 |
| `applyCHy()` | Controlled-Hy | 🔲 TODO | 🔲 TODO | 4 |
| `applyCT()` | Controlled-T | 🔲 TODO | 🔲 TODO | 4 |
| `applyCTdagger()` | Controlled-T† | 🔲 TODO | 🔲 TODO | 4 |
| `applyCP()` | Controlled-P(θ) | 🔲 TODO | 🔲 TODO | 4 |
| `applyCPdagger()` | Controlled-P†(θ) | 🔲 TODO | 🔲 TODO | 4 |
| `applySWAP()` | SWAP | 🔲 TODO | 🔲 TODO | 4 |
| `applyRSWAP()` | Rotated SWAP | 🔲 TODO | 🔲 TODO | 4 |

### Three-Qubit Gate Implementations

| Gate | Description | Pure State | Mixed State | Phase |
|------|-------------|------------|-------------|-------|
| `applyToffoli()` | Toffoli (CCNOT) | 🔲 TODO | 🔲 TODO | 5 |

---

## Design Notes

### Position-Dependent Stride Logic

**Why different positions matter:**
qHiPSTER implementation uses stride-based indexing where stride depends on target qubit index:
```c
// Simplified example for single-qubit gate
stride = 1ULL << target;  // 2^target
for (i = 0; i < (1ULL << n); i++) {
    if ((i & stride) == 0) {
        // Operate on amplitude pair: state[i] and state[i + stride]
    }
}
```

Different target qubits → different strides → different memory access patterns → genuinely different
code execution paths.

**Testing implication**: Must validate all qubit positions to ensure stride calculations are correct.

### Numerical Precision

**Test tolerance**: ε = 1e-12 for comparing production implementations against BLAS reference matrix
multiplication.

**Rationale**: Accounts for double-precision arithmetic (machine epsilon ~1e-16) plus accumulated
rounding errors over O(2^n) sparse operations.

**Numerical stability advantage**: Sparse operations accumulate errors through O(2^n) operations per
gate vs O(2^3n) for dense matrix multiplication—better stability at scale.

### Packed Storage for Mixed States

**Production**: Gates operate directly on lower-triangular packed format (no unpacking)
- Generalized qHiPSTER algorithms for packed density matrices
- 50% memory savings critical for n=16 mixed states (32GB vs 64GB)

**Testing**: Unpack → full matrix BLAS operations (UρU†) → repack
- Simple, verifiable reference implementation
- Isolates pack/unpack logic from gate logic
- Systematic failures → packing bugs; gate-specific failures → gate bugs

### State Representation

**Unified `state_t` structure:**
```c
typedef struct {
    type_t type;  // PURE or MIXED
    int n;        // number of qubits
    cplx_t* data; // amplitudes (pure) or packed density matrix (mixed)
} state_t;
```

**API design**: Single function per gate dispatches based on type
```c
int applyX(state_t* state, int target) {
    return (state->type == PURE)
        ? applyX_pure(state, target)
        : applyX_mixed(state, target);
}
```

**Benefits**: Type-safe at runtime, single API for users, internal implementations optimized per type.

### Thread Safety and Parallelization

**State objects**: Not thread-safe for concurrent mutation of same object. Users can create separate
state objects for parallel circuit evaluation.

**Empirical parallelization results** (OpenMP):
- Gate-level: ~3× speedup on 10 cores (memory bandwidth limited)
- Circuit-level: ~5× speedup
- Algorithm-level: Best scaling but RAM-constrained at target qubit counts

**Library scope**: Provides thread-safe object creation. Algorithm-level parallelization implemented
by users external to library.

---

## Implementation Workflow

### Overview

**Iterative development strategy**: Build one gate end-to-end (pure + mixed, tested) before moving to
next gate. Discover infrastructure needs organically rather than big-design-up-front.

**Rationale**:
- Early validation of test infrastructure
- Discover design issues sooner
- Clear "done" criteria per phase
- Working gates immediately available
- Common patterns emerge through iteration

### Phase 0: Mixed State Testing Infrastructure
**Goal**: Enable testing of single-qubit gates on mixed states

**Deliverables**:
- `density_pack()`: Full Hermitian matrix → packed lower-triangular storage
- `density_unpack()`: Packed storage → full Hermitian matrix
- `testSingleQubitGateMixed()`: Test harness analogous to `testSingleQubitGate()`
  - Function pointer + 2×2 matrix → automated testing
  - Unpack density matrix → apply UρU† via BLAS → repack
  - Compare with production implementation output
- Mixed state generators:
  - `test_pure_as_mixed()`: Convert existing pure states to density matrices
  - `test_maximally_mixed()`: Generate I/2^n
  - `test_random_mixture()`: Generate ρ = Σpᵢ|ψᵢ⟩⟨ψᵢ| with random weights

**Exit criteria**: Can test arbitrary single-qubit gate on mixed states using same pattern as pure
state testing.

### Phase 1: Complete One Single-Qubit Gate End-to-End
**Goal**: First fully validated gate (pure + mixed) demonstrating complete workflow

**Suggested gate**: X or H (simple, well-understood)

**Deliverables**:
- Implement mixed state version of chosen gate
  - Direct operation on packed density matrix
  - Use generalized qHiPSTER stride logic
- Validate against test harness (uses Phase 0 infrastructure)
- Document any issues or patterns discovered

**Exit criteria**:
- Both pure and mixed implementations pass all tests
- Clear implementation pattern established for remaining gates

**Learning outcomes**:
- Validate Phase 0 infrastructure works as intended
- Identify any missing test cases or helper functions
- Establish code patterns for subsequent gates

### Phase 2: Remaining Single-Qubit Gates
**Goal**: Complete single-qubit gate library

**Deliverables**:
- Implement all remaining single-qubit gates (H, S, T, rotations, daggers)
- Both pure and mixed versions
- Follow pattern established in Phase 1
- Mostly copy-paste-modify with different 2×2 matrices

**Exit criteria**: All single-qubit gates pass tests (pure + mixed)

**Expected discoveries**:
- Common helper functions needed
- Opportunities for code reuse
- Edge cases in rotation gates (small angles, boundary values)

### Phase 3: Entangled State Generators
**Goal**: Enable testing of two-qubit gates

**Deliverables**:
- `test_bell_states()`: Generate 4 Bell states
- `test_ghz_state()`: Generate n-qubit GHZ states
- `test_w_state()`: Generate n-qubit W states

**Rationale**: Two-qubit gates (especially CNOT) require entangled test states to validate entanglement
propagation.

**Exit criteria**: Can generate entangled states for n=2,3,4 qubit systems.

### Phase 4: Two-Qubit Gates
**Goal**: Complete two-qubit gate library

**Deliverables**:
- `mat_two_qubit_gate()`: Build full matrix for two-qubit gates via Kronecker products
- `testTwoQubitGate()`: Pure state test harness for two-qubit gates
  - Loop through (target, control) pairs
  - Generate full matrix, compare outputs
- `testTwoQubitGateMixed()`: Mixed state test harness
- Implement all two-qubit gates (CNOT, CZ, SWAP, controlled operations)
  - Pure and mixed versions
  - Iterate one gate at a time

**Exit criteria**: All two-qubit gates pass tests (pure + mixed)

**Challenges**:
- More complex Kronecker product construction (4×4 → 2^n×2^n)
- Two index parameters (target, control)
- Control qubit logic validation

### Phase 5: Three-Qubit Gates
**Goal**: Complete Toffoli gate

**Deliverables**:
- `mat_three_qubit_gate()`: Build full matrix for three-qubit gates
- `testThreeQubitGate()` and `testThreeQubitGateMixed()`: Test harnesses
- Implement Toffoli (pure + mixed)

**Exit criteria**: Toffoli passes all tests

**Final milestone**: Complete gate library with comprehensive test coverage

---

## Validation and Quality Assurance

### Build-Time Testing

**Requirement**: All tests complete in <2 minutes

**CI/CD integration**: Tests run on every build, abort on failure

**Coverage verification**: Automated check ensures all gates tested across:
- All qubit positions
- All system sizes (n=1,2,3,4)
- All state types (pure, mixed)
- All basis states and entangled states

### Memory Leak Testing

**Strategy**: Regular checks during development

**Tools**: Platform-specific (valgrind on Linux, Instruments on macOS)

**Discipline**: Manual verification during iterative development

### Numerical Validation

**Tolerance**: ε = 1e-12

**Checks per test**:
- Output state matches reference within tolerance
- Norm preservation (pure states): |⟨ψ'|ψ'⟩ - 1| < ε
- Trace preservation (mixed states): |Tr(ρ') - Tr(ρ)| < ε
- Hermiticity (mixed states): verified by storage format

**Failure analysis**: Distinguish between:
- Systematic failures (infrastructure bugs)
- Gate-specific failures (implementation bugs)
- Numerical precision issues (tolerance too tight)

---

## Future Considerations

**GPU acceleration**: Gate operations are memory-bandwidth limited—potential for GPU speedup in
streaming operations

**Extended gate set**: Additional gates based on algorithm requirements (e.g., multi-controlled gates,
custom parameterized gates)

**Performance optimization**: Profile-guided optimization after correctness established

**Noise models**: Post-v1.0 extension for realistic hardware simulation
