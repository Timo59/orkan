# Gates Module: Technical Specification

**Module:** Quantum Gate Operations
**Header:** `include/qhipster.h`
**Implementation:** `src/qhipster.c`
**Last Updated:** 2025-12-23

---

## Unit Tests

### Verification Approach

Gate functions validated against reference implementations:
- **Pure states**: $ \lvert \psi^{\prime} \rangle = U \lvert \psi \rangle $ via matrix-vector multiplication (`zgemv()`)
- **Mixed states**: $ \rho^{\prime} = U \rho U^{\ast} $ via matrix-matrix multiplication (`zgemm()`)

Reference implementations use full matrix representation with hardcoded Kronecker products (I⊗...⊗U⊗...⊗I). Production implementations exploit sparse local structure and operate directly on packed arrays.

### Test State Basis

**Single-qubit informationally complete set (6 states)**:
- Computational: $ \lvert 0 \rangle, \ \lvert 1 \rangle $
- Hadamard: $ \lvert + \rangle, \ \lvert - \rangle $
- Circular: $ \lvert +i \rangle, \ \lvert -i \rangle $

**Rationale**: Detects amplitude errors (computational), real phase errors (Hadamard), complex phase errors (circular).

**Multi-qubit pure states**:
- Tensor products of single-qubit bases
- Bell states (n=2): $ \frac{1}{\sqrt{2}} (\lvert 00 \rangle + \lvert 11 \rangle) $,
$ \frac{1}{\sqrt{2}} (\lvert 00 \rangle - \lvert 11 \rangle) $,
$ \frac{1}{\sqrt{2}} (\lvert 01 \rangle + \lvert 10 \rangle) $,
$ \frac{1}{\sqrt{2}} (\lvert 01 \rangle - \lvert 10 \rangle) $
- GHZ states (n≥3): $ \frac{1}{\sqrt{2}} (\lvert 0 \dots 0 \rangle + \lvert 1 \dots 1 \rangle) $
- W states (n≥3): $\frac{1}{\sqrt{n}}(\lvert 100\dots0 \rangle + \lvert 010\dots0 \rangle + \lvert 000\dots1\rangle)$

**Rationale**: Bell/GHZ/W test entanglement propagation and multi-partite correlations.

**Mixed states**:
- Pure-as-mixed: All pure states as density matrices $ \lvert \psi \rangle \langle \psi \rvert $
- Maximally mixed: $ I/2^{n} $
- Rank > 1 examples:
  - Thermal: $ (3 \lvert 0 \rangle \langle 0 \rvert + \lvert 1 \rangle \langle 1 \rvert )/4 $
    - Werner (n=2): $ p \lvert \phi{+} \rangle \langle \phi^{+} \rvert + \frac{(1-p)}{4} I $
  - Separable mixtures

**Rationale**: Rank>1 states verify conjugation operates on full density matrix, not just extracted pure components.

### Test Coverage

**Single-qubit gates**: Test all target positions (0 to n-1) on n=1,2,3,4 qubit systems
**Two-qubit gates**: Test all (target, control) pairs on n=2,3,4 qubit systems
**Three-qubit gates**: Test all (target, control1, control2) triples on n=3,4 qubit systems

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

**Estimated runtime**: 10-30 minutes (one-time build verification)

---

## Test Infrastructure Status

### State Generators

| Function | Generates | Location | Status |
|----------|-----------|----------|--------|
| `test_cb_pure()` | Computational basis states | `test.c` | ✅ Complete |
| `test_xb_pure()` | Hadamard basis states | `test.c` | ✅ Complete |
| `test_yb_pure()` | Circular basis states | `test.c` | ✅ Complete |
| `test_gen_states_pure()` | Combined: all three bases | `test.c` | ✅ Complete |
| `test_rm_states_pure()` | Cleanup: free test states | `test.c` | ✅ Complete |
| `test_bell_states()` | 4 Bell states | `test.c` | 🔲 TODO |
| `test_ghz_state()` | n-qubit GHZ | `test.c` | 🔲 TODO |
| `test_w_state()` | n-qubit W | `test.c` | 🔲 TODO |
| `test_maximally_mixed()` | I/2ⁿ | `test.c` | 🔲 TODO |
| `test_thermal_mixed()` | Thermal state (bias) | `test.c` | 🔲 TODO |
| `test_werner_state()` | Werner state | `test.c` | 🔲 TODO |

### Test Harness Functions

| Function | Purpose | Location | Status |
|----------|---------|----------|--------|
| `testSingleQubitGate()` | Validates single-qubit gates on pure states | `test_qhipster.c` | ✅ Complete |
| `testSingleQubitGateMixed()` | Validates single-qubit gates on mixed states | `test_qhipster.c` | 🔲 TODO |
| `testTwoQubitGate()` | Validates two-qubit gates on pure states | `test_qhipster.c` | 🔲 TODO |
| `testTwoQubitGateMixed()` | Validates two-qubit gates on mixed states | `test_qhipster.c` | 🔲 TODO |
| `testThreeQubitGate()` | Validates three-qubit gates on pure states | `test_qhipster.c` | 🔲 TODO |
| `testThreeQubitGateMixed()` | Validates three-qubit gates on mixed states | `test_qhipster.c` | 🔲 TODO |

### Test Helper Functions

| Function | Purpose | Location | Status |
|----------|---------|----------|--------|
| `mv()` | Matrix-vector multiplication via `zgemv()` | `test_qhipster.c` | ✅ Complete |
| `kron()` | Kronecker product A⊗B | `test_qhipster.c` | ✅ Complete |
| `mat_id()` | Identity matrix I(2ⁿ×2ⁿ) | `test_qhipster.c` | ✅ Complete |
| `mat_single_qubit_gate()` | Full matrix I⊗...⊗U⊗...⊗I for single-qubit gate | `test_qhipster.c` | ✅ Complete |
| `mat_two_qubit_gate()` | Full matrix for two-qubit gate | `test_qhipster.c` | 🔲 TODO |
| `mat_three_qubit_gate()` | Full matrix for three-qubit gate | `test_qhipster.c` | 🔲 TODO |
| `density_unpack()` | Packed lower-triangle → full Hermitian matrix | `test_qhipster.c` | 🔲 TODO |
| `density_pack()` | Full Hermitian matrix → packed storage | `test_qhipster.c` | 🔲 TODO |

**Note**: Test files include predefined gate matrices (XMAT, YMAT, ZMAT, HMAT, SMAT, TMAT, SWAPMAT, etc.) in `test/include/test_qhipster.h`.

---

## Production Code Implementation Status

**Location**: `src/qhipster.c`, `include/qhipster.h`

### Single-Qubit Gate Implementations

| Gate | Description | Location | Status |
|------|-------------|----------|--------|
| `applyX()` | Pauli-X | `qhipster.c` | ✅ Complete |
| `applyY()` | Pauli-Y | `qhipster.c` | ✅ Complete |
| `applyZ()` | Pauli-Z | `qhipster.c` | ✅ Complete |
| `applyH()` | Hadamard | `qhipster.c` | 🔲 TODO |
| `applyS()` | Phase gate | `qhipster.c` | 🔲 TODO |
| `applySdagger()` | S† | `qhipster.c` | 🔲 TODO |
| `applyT()` | π/8 gate | `qhipster.c` | 🔲 TODO |
| `applyTdagger()` | T† | `qhipster.c` | 🔲 TODO |
| `applyHy()` | Hadamard-Y | `qhipster.c` | 🔲 TODO |
| `applyP()` | Phase rotation P(θ) | `qhipster.c` | 🔲 TODO |
| `applyPdagger()` | P†(θ) | `qhipster.c` | 🔲 TODO |
| `applyRX()` | X-rotation RX(θ) | `qhipster.c` | 🔲 TODO |
| `applyRXdagger()` | RX†(θ) | `qhipster.c` | 🔲 TODO |
| `applyRY()` | Y-rotation RY(θ) | `qhipster.c` | 🔲 TODO |
| `applyRYdagger()` | RY†(θ) | `qhipster.c` | 🔲 TODO |
| `applyRZ()` | Z-rotation RZ(θ) | `qhipster.c` | 🔲 TODO |
| `applyRZdagger()` | RZ†(θ) | `qhipster.c` | 🔲 TODO |

### Two-Qubit Gate Implementations

| Gate | Description | Location | Status |
|------|-------------|----------|--------|
| `applyCX()` | CNOT | `qhipster.c` | 🔲 TODO |
| `applyCY()` | Controlled-Y | `qhipster.c` | 🔲 TODO |
| `applyCZ()` | Controlled-Z | `qhipster.c` | 🔲 TODO |
| `applyCS()` | Controlled-S | `qhipster.c` | 🔲 TODO |
| `applyCSdagger()` | Controlled-S† | `qhipster.c` | 🔲 TODO |
| `applyCH()` | Controlled-H | `qhipster.c` | 🔲 TODO |
| `applyCHy()` | Controlled-Hy | `qhipster.c` | 🔲 TODO |
| `applyCT()` | Controlled-T | `qhipster.c` | 🔲 TODO |
| `applyCTdagger()` | Controlled-T† | `qhipster.c` | 🔲 TODO |
| `applyCP()` | Controlled-P(θ) | `qhipster.c` | 🔲 TODO |
| `applyCPdagger()` | Controlled-P†(θ) | `qhipster.c` | 🔲 TODO |
| `applySWAP()` | SWAP | `qhipster.c` | 🔲 TODO |
| `applyRSWAP()` | Rotated SWAP | `qhipster.c` | 🔲 TODO |

### Three-Qubit Gate Implementations

| Gate | Description | Location | Status |
|------|-------------|----------|--------|
| `applyToffoli()` | Toffoli (CCNOT) | `qhipster.c` | 🔲 TODO |

### Mixed State Gate Implementations

All gates above require mixed state (density matrix) support:
- Pure state implementation (operate on statevector)
- Mixed state implementation (operate on packed density matrix via UρU†)
- **Status**: X, Y, Z pure-only; all mixed implementations TODO

---

## Design Notes

**Numerical tolerance**: BLAS precision (machine epsilon ~1e-16, relaxed to ~1e-12 for accumulated errors)

**Packed storage for mixed states**: Production gates operate on lower-triangular packed format; tests unpack → conjugate → repack for verification

**Optimization**: All gates exploit sparse local structure (stride/step indexing) rather than full matrix multiplication

**Thread safety**: Not thread-safe (same as state module)

---

## Implementation Priorities

**Phase 1**: Helper functions (pack/unpack, reference implementations, Kronecker products)
**Phase 2**: State generators (Bell, GHZ, W, mixed states)
**Phase 3**: Complete single-qubit gates (Clifford, rotations)
**Phase 4**: Two-qubit gates (controlled operations, SWAP)
**Phase 5**: Mixed state support for all gates
**Phase 6**: Full test suite implementation
