# Gates Module: Technical Specification

**Module:** Quantum Gate Operations
**Header:** `include/qhipster.h`
**Implementation:** `src/qhipster.c`
**Last Updated:** 2025-12-23

---

## Test Strategy

### Verification Approach

Gate functions validated against reference implementations:
- **Pure states**: `|ψ'⟩ = U|ψ⟩` via `zgemv()`
- **Mixed states**: `ρ' = UρU†` via `zgemm()`

Reference implementations use full matrix representation with hardcoded Kronecker products (I⊗...⊗U⊗...⊗I). Production implementations exploit sparse local structure and operate directly on packed arrays.

### Test State Basis

**Single-qubit informationally complete set (6 states)**:
- Computational: |0⟩, |1⟩
- Hadamard: |+⟩, |-⟩
- Circular: |+i⟩, |-i⟩

**Rationale**: Detects amplitude errors (computational), real phase errors (Hadamard), complex phase errors (circular).

**Multi-qubit pure states**:
- Tensor products of single-qubit bases
- Bell states (n=2): |Φ⁺⟩, |Φ⁻⟩, |Ψ⁺⟩, |Ψ⁻⟩
- GHZ states (n≥3): (|0...0⟩ + |1...1⟩)/√2
- W states (n≥3): Equal superposition of single-excitation states

**Rationale**: Bell/GHZ/W test entanglement propagation and multi-partite correlations.

**Mixed states**:
- Pure-as-mixed: All pure states as density matrices |ψ⟩⟨ψ|
- Maximally mixed: I/2ⁿ
- Rank>1 examples:
  - Thermal: (3|0⟩⟨0| + |1⟩⟨1|)/4
  - Werner (n=2): p|Φ⁺⟩⟨Φ⁺| + (1-p)I/4
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

## Implementation Status

### Helper Functions

| Function | Purpose | Status |
|----------|---------|--------|
| `density_unpack()` | Packed lower-triangle → full Hermitian matrix | 🔲 TODO |
| `density_pack()` | Full Hermitian matrix → packed storage | 🔲 TODO |
| `apply_full_matrix_pure()` | Reference: U\|ψ⟩ via zgemv() | 🔲 TODO |
| `apply_conjugation_mixed()` | Reference: UρU† via zgemm() | 🔲 TODO |
| `kronecker_product()` | Hardcoded I⊗U⊗I constructions | 🔲 TODO |

### State Generators (Test Utilities)

| Function | Generates | Status |
|----------|-----------|--------|
| `test_cb_pure()` | Computational basis states | ✅ Complete |
| `test_xb_pure()` | Hadamard basis states | ✅ Complete |
| `test_yb_pure()` | Circular basis states | ✅ Complete |
| `test_bell_states()` | 4 Bell states | 🔲 TODO |
| `test_ghz_state()` | n-qubit GHZ | 🔲 TODO |
| `test_w_state()` | n-qubit W | 🔲 TODO |
| `test_maximally_mixed()` | I/2ⁿ | 🔲 TODO |
| `test_thermal_mixed()` | Thermal state (bias) | 🔲 TODO |
| `test_werner_state()` | Werner state | 🔲 TODO |

### Single-Qubit Gate Implementations

| Gate | Description | Status |
|------|-------------|--------|
| `applyX()` | Pauli-X | ✅ Complete |
| `applyY()` | Pauli-Y | ✅ Complete |
| `applyZ()` | Pauli-Z | ✅ Complete |
| `applyH()` | Hadamard | 🔲 TODO |
| `applyS()` | Phase gate | 🔲 TODO |
| `applySdagger()` | S† | 🔲 TODO |
| `applyT()` | π/8 gate | 🔲 TODO |
| `applyTdagger()` | T† | 🔲 TODO |
| `applyHy()` | Hadamard-Y | 🔲 TODO |
| `applyP()` | Phase rotation P(θ) | 🔲 TODO |
| `applyPdagger()` | P†(θ) | 🔲 TODO |
| `applyRX()` | X-rotation RX(θ) | 🔲 TODO |
| `applyRXdagger()` | RX†(θ) | 🔲 TODO |
| `applyRY()` | Y-rotation RY(θ) | 🔲 TODO |
| `applyRYdagger()` | RY†(θ) | 🔲 TODO |
| `applyRZ()` | Z-rotation RZ(θ) | 🔲 TODO |
| `applyRZdagger()` | RZ†(θ) | 🔲 TODO |

### Two-Qubit Gate Implementations

| Gate | Description | Status |
|------|-------------|--------|
| `applyCX()` | CNOT | 🔲 TODO |
| `applyCY()` | Controlled-Y | 🔲 TODO |
| `applyCZ()` | Controlled-Z | 🔲 TODO |
| `applyCS()` | Controlled-S | 🔲 TODO |
| `applyCSdagger()` | Controlled-S† | 🔲 TODO |
| `applyCH()` | Controlled-H | 🔲 TODO |
| `applyCHy()` | Controlled-Hy | 🔲 TODO |
| `applyCT()` | Controlled-T | 🔲 TODO |
| `applyCTdagger()` | Controlled-T† | 🔲 TODO |
| `applyCP()` | Controlled-P(θ) | 🔲 TODO |
| `applyCPdagger()` | Controlled-P†(θ) | 🔲 TODO |
| `applySWAP()` | SWAP | 🔲 TODO |
| `applyRSWAP()` | Rotated SWAP | 🔲 TODO |

### Three-Qubit Gate Implementations

| Gate | Description | Status |
|------|-------------|--------|
| `applyToffoli()` | Toffoli (CCNOT) | 🔲 TODO |

### Mixed State Gate Implementations

All gates above require mixed state (density matrix) support:
- Pure state implementation (operate on statevector)
- Mixed state implementation (operate on packed density matrix via UρU†)
- **Status**: X, Y, Z pure-only; all mixed implementations TODO

### Unit Tests

| Test Suite | Target | Status |
|------------|--------|--------|
| `test_pauli_gates_pure()` | X, Y, Z on pure states | 🔲 TODO |
| `test_pauli_gates_mixed()` | X, Y, Z on mixed states | 🔲 TODO |
| `test_clifford_gates_pure()` | H, S, T on pure states | 🔲 TODO |
| `test_clifford_gates_mixed()` | H, S, T on mixed states | 🔲 TODO |
| `test_rotation_gates_pure()` | RX, RY, RZ on pure states | 🔲 TODO |
| `test_rotation_gates_mixed()` | RX, RY, RZ on mixed states | 🔲 TODO |
| `test_controlled_gates_pure()` | CX, CZ, etc. on pure states | 🔲 TODO |
| `test_controlled_gates_mixed()` | CX, CZ, etc. on mixed states | 🔲 TODO |
| `test_toffoli_pure()` | Toffoli on pure states | 🔲 TODO |
| `test_toffoli_mixed()` | Toffoli on mixed states | 🔲 TODO |
| `test_gate_invariants()` | Normalization, hermiticity | 🔲 TODO |

**Test file**: `test/src/test_gates.c` (not yet created)

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
