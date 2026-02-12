# Gates Module: Technical Specification

**Module:** Quantum Gate Operations
**Header:** `include/gate.h`
**Implementation:** `src/gate/gate_pure.c` (pure), `src/gate/packed/gate_packed_1q.c` + `src/gate/packed/gate_packed_<name>.c` (mixed packed, one file per multi-qubit gate), `src/gate/tiled/gate_tiled.c` (mixed tiled 1Q + rotation), `src/gate/tiled/gate_tiled_<name>.c` (mixed tiled, one file per multi-qubit gate), `src/gate/gate.c` (dispatchers)
**Dependencies:** `state.h`, `q_types.h`, `<complex.h>`, `<math.h>`, OpenMP (optional)

---

## 1. Overview

The gate module implements quantum gate operations on the `state_t` struct's data buffer. Gates exploit
sparse, local structure of quantum operations — operating directly on state arrays via position-dependent
stride indexing rather than constructing full gate matrices via Kronecker products.

**Design characteristics:**
- **Direct packed operations**: mixed state gates operate directly on lower-triangular packed or tiled
  storage (no unpacking required)
- **Macro-based traversal**: common loop structures factored into traversal macros with gate-specific
  operations passed as macro parameters
- **In-place**: all gates modify state data in-place with O(1) extra memory

---

## 2. Gate Inventory

### 2.1 Single-Qubit Gates

| Gate | Function | Matrix | Category | Pure | Packed | Tiled* |
|------|----------|--------|----------|------|--------|--------|
| Pauli-X | `x(state, target)` | `[[0,1],[1,0]]` | Pauli | ✓ | ✓ | ✓ |
| Pauli-Y | `y(state, target)` | `[[0,-i],[i,0]]` | Pauli | ✓ | ✓ | ✓ |
| Pauli-Z | `z(state, target)` | `[[1,0],[0,-1]]` | Pauli | ✓ | ✓ | ✓ |
| Hadamard | `h(state, target)` | `[[1,1],[1,-1]]/√2` | Clifford | ✓ | ✓ | ✓ |
| Phase S | `s(state, target)` | `[[1,0],[0,i]]` | Clifford | ✓ | ✓ | ✓ |
| S-dagger | `sdg(state, target)` | `[[1,0],[0,-i]]` | Clifford | ✓ | ✓ | ✓ |
| T | `t(state, target)` | `[[1,0],[0,e^(iπ/4)]]` | Non-Clifford | ✓ | ✓ | ✓ |
| T-dagger | `tdg(state, target)` | `[[1,0],[0,e^(-iπ/4)]]` | Non-Clifford | ✓ | ✓ | ✓ |
| Hadamard-Y | `hy(state, target)` | `[[1,-1],[1,1]]/√2` | Clifford | ✓ | ✓ | ✓ |
| P(θ) | `p(state, target, θ)` | `[[1,0],[0,e^(iθ)]]` | Rotation | ✓ | ✓ | ✓ |
| Rx(θ) | `rx(state, target, θ)` | `[[c,-is],[-is,c]]` | Rotation | ✓ | ✓ | ✓ |
| Ry(θ) | `ry(state, target, θ)` | `[[c,-s],[s,c]]` | Rotation | ✓ | ✓ | ✓ |
| Rz(θ) | `rz(state, target, θ)` | `[[e^(-iθ/2),0],[0,e^(iθ/2)]]` | Rotation | ✓ | ✓ | ✓ |

*Tiled implementations exist but dispatch is temporarily disabled (under construction).

### 2.2 Two-Qubit Gates

Two-qubit gates take `(state, control, target)` for controlled gates or `(state, q1, q2)` for symmetric
gates.

| Gate | Function | Category | Pure | Packed | Tiled* |
|------|----------|----------|------|--------|--------|
| CNOT | `cx(state, ctrl, tgt)` | Controlled-Pauli | ✓ | ✓ | ✓ |
| CY | `cy(state, ctrl, tgt)` | Controlled-Pauli | ✓ | ✓ | — |
| CZ | `cz(state, ctrl, tgt)` | Controlled-Pauli | ✓ | ✓ | — |
| SWAP | `swap(state, q1, q2)` | Permutation | ✓ | ✓ | ✓ |

**CZ symmetry:** CZ(a,b) = CZ(b,a).

**RZZ(θ) — future work:** RZZ(θ) = diag(e^(-iθ/2), e^(iθ/2), e^(iθ/2), e^(-iθ/2)) is relevant for
QAOA circuits. Since Z⊗Z is diagonal, RZZ applies only phases — no off-diagonal mixing. Its packed and
tiled implementations require (row, col)-dependent phase scaling rather than block mixing. Deferred to a
future revision.

### 2.3 Three-Qubit Gates

| Gate | Function | Description | Pure | Packed | Tiled* |
|------|----------|-------------|------|--------|--------|
| Toffoli | `ccx(state, ctrl1, ctrl2, tgt)` | Doubly-controlled X | ✓ | ✓ | — |

---

## 3. Dispatch Architecture

### 3.1 Three-Way Dispatch

Every public gate function validates inputs and dispatches to a type-specific implementation:

```c
void x(state_t *state, const qubit_t target) {
    GATE_VALIDATE(state && state->data, "x: null state or data pointer");
    GATE_VALIDATE(target < state->qubits, "x: target qubit out of range");

    switch (state->type) {
        case PURE:         x_pure(state, target);   break;
        case MIXED_PACKED: x_packed(state, target);  break;
        case MIXED_TILED:  /* temporarily disabled (under construction) */
    }
}
```

### 3.2 Error Handling

Gate functions use `GATE_VALIDATE` — a macro that prints a diagnostic to stderr and calls `exit(EXIT_FAILURE)`:

```c
#define GATE_VALIDATE(cond, msg) do {                           \
    if (!(cond)) {                                              \
        fprintf(stderr, "qlib: gate: %s\n", (msg));            \
        exit(EXIT_FAILURE);                                     \
    }                                                           \
} while(0)
```

**Validated conditions:**

| Condition | Message |
|-----------|---------|
| `state && state->data` | null state or data pointer |
| `target < state->qubits` | target qubit out of range |
| `control < state->qubits` | control qubit out of range |
| `control != target` | control and target must differ |

---

## 4. Pure State Gate Algorithms

### 4.1 Single-Qubit Gates: TRAVERSE_PURE_1Q

For qubit q in an n-qubit state, pairs (ψ[idx], ψ[idx + 2^q]) are processed independently.

```c
#define TRAVERSE_PURE_1Q(state, target, PAIR_OP)
    stride = 1 << target;
    for (i = 0; i < dim; i += 2*stride) {
        for (j = 0; j < stride; ++j) {
            a = data + i + j;           // target bit = 0
            b = data + i + j + stride;  // target bit = 1
            PAIR_OP(a, b);
        }
    }
```

OpenMP threshold: `dim >= 4096` (n ≥ 12 qubits).

### 4.2 Two-Qubit Gates: Flat Loop with Two-Bit Insertion

For a controlled-U gate with control c and target t, enumerate all base indices where both c and t bits
are 0 using `insertBits2_0()`. Define `lo = min(c, t)`, `hi = max(c, t)`:

```c
lo = min(control, target);
hi = max(control, target);
incr_ctrl = 1 << control;
incr_tgt  = 1 << target;
n_base = dim >> 2;  // dim / 4

#pragma omp parallel for if(dim >= OMP_THRESHOLD)
for (k = 0; k < n_base; ++k) {
    base = insertBits2_0(k, lo, hi);
    a = data + base + incr_ctrl;             // control=1, target=0
    b = data + base + incr_ctrl + incr_tgt;  // control=1, target=1
    PAIR_OP(a, b);
}
```

No traversal macro exists for two-qubit gates; each gate inlines the loop with gate-specific operations.

**SWAP:** Not a controlled-U gate. Uses the same flat loop but swaps a different pair:

```c
for (k = 0; k < n_base; ++k) {
    base = insertBits2_0(k, lo, hi);
    idx_01 = base | incr_lo;    // lo-bit=1, hi-bit=0
    idx_10 = base | incr_hi;    // lo-bit=0, hi-bit=1
    swap(data[idx_01], data[idx_10]);
}
```

### 4.3 Three-Qubit Gates: Flat Loop with Three-Bit Insertion

For Toffoli (CCX) with controls c1, c2 and target t, sort the three qubit positions: `lo < mid < hi`,
then enumerate base indices using chained `insertBit0()`:

```c
n_base = dim >> 3;  // dim / 8

#pragma omp parallel for if(dim >= OMP_THRESHOLD)
for (k = 0; k < n_base; ++k) {
    base = insertBit0(insertBit0(insertBit0(k, lo), mid), hi);
    idx_c11_t0 = base | incr_c1 | incr_c2;
    idx_c11_t1 = base | incr_c1 | incr_c2 | incr_tgt;
    swap(data[idx_c11_t0], data[idx_c11_t1]);
}
```

Higher-order multi-controlled gates belong to a separate decomposition module.

### 4.4 Per-Gate Optimizations (Pure State)

| Gate | Operation | Optimization |
|------|-----------|-------------|
| X | swap(a, b) | Zero arithmetic (memory only) |
| Y | a'=-ib, b'=ia | Zero multiplications (swap + negate components) |
| Z | b' = -b | Single negation |
| S | b' = ib | Zero multiplications (component swap) |
| Sdg | b' = -ib | Zero multiplications (component swap) |
| T | b' = e^{iπ/4}b | 2 adds + 2 muls via `(re-im, re+im)/√2` |
| Tdg | b' = e^{-iπ/4}b | 2 adds + 2 muls via `(re+im, im-re)/√2` |
| H | butterfly | 2 adds per pair, 2 muls (by 1/√2) |
| Rx | 2×2 unitary mixing | cos/sin precomputed |
| Ry | 2×2 real rotation | Real coefficients only |
| Rz | diagonal phases | a' = e^(-iθ/2)a, b' = e^(iθ/2)b |

---

## 5. Mixed State Gate Algorithms (Packed)

### 5.1 Single-Qubit Gates: TRAVERSE_PACKED_BLOCKS

For ρ → UρU†, the conjugation mixes 2×2 blocks of the density matrix indexed by the target qubit.

```c
#define TRAVERSE_PACKED_BLOCKS(state, target, BLOCK_OP)
    for (bc = 0; bc < half_dim; ++bc) {
        c0 = insertBit0(bc, target);
        c1 = c0 | incr;
        for (br = bc; br < half_dim; ++br) {
            r0 = insertBit0(br, target);
            r1 = r0 | incr;
            // Compute packed indices for 4 elements
            idx00 = pack_idx(dim, r0, c0);
            idx11 = pack_idx(dim, r1, c1);
            idx10 = pack_idx(dim, r1, c0);
            idx01 = lower_01 ? pack_idx(dim, r0, c1) : pack_idx(dim, c1, r0);
            BLOCK_OP(data, idx00, idx01, idx10, idx11, lower_01, diag_block);
        }
    }
```

OpenMP: `schedule(static)` on the outer `bc` loop. Threshold: `dim >= 64` (n ≥ 6 qubits).

### 5.2 Two-Qubit Gates: Two-Bit Insertion with 4×4 Block Pairs

Two-qubit gates use `insertBits2_0(idx, lo, hi)` to enumerate base indices where both qubit bits
are zero, producing 4×4 sub-blocks of the density matrix. Each block has 4 row indices
(r00, r01, r10, r11) and 4 column indices (c00, c01, c10, c11).

```c
quarter_dim = dim >> 2;
for (bc = 0; bc < quarter_dim; ++bc) {
    c00 = insertBits2_0(bc, lo, hi);
    c01 = c00 | incr_tgt;  c10 = c00 | incr_ctrl;  c11 = c10 | incr_tgt;
    for (br = bc; br < quarter_dim; ++br) {
        r00 = insertBits2_0(br, lo, hi);
        r01 = r00 | incr_tgt;  r10 = r00 | incr_ctrl;  r11 = r10 | incr_tgt;
        // Process element pairs within the 4×4 block
    }
}
```

The triangular iteration (br >= bc) visits each block once. Within each block, the gate's
action decomposes into swap/transform pairs. Triangle membership determines access mode:

- **Fast path** (r00 > c11): all 16 elements in lower triangle — direct access
- **Slow path**: some elements cross the diagonal; access via `conj(data[pack_idx(dim, c, r)])`
- **Diagonal blocks** (bc == br): Hermitian pairs share storage; one side must be skipped

OpenMP: `schedule(static, 1)` on the outer `bc` loop (round-robin balances the triangular
iteration). Threshold: `dim >= 64`.

**Per-gate operations:**

| Gate | Type | Operation |
|------|------|-----------|
| CX | Zero-flop permutation | 6 swap pairs; perm(x) = x ^ (1<<target) when control bit set |
| CY | Phased permutation | 6 swap pairs with per-pair phase factors (±i, negation) |
| CZ | Zero-flop diagonal | Negate elements where exactly one of (row, col) has both bits set |
| SWAP | Zero-flop permutation | 6 swap pairs; perm(x) = x with bits lo and hi exchanged |

**CX, CZ, and SWAP are zero-flop** — only memory moves and sign flips.

### 5.3 Three-Qubit Gates: Three-Bit Insertion with 8×8 Block Pairs

CCX (Toffoli) uses the same structural pattern as two-qubit packed gates but with three-bit insertion
via chained `insertBit0()` (same as the pure state algorithm in Section 4.3). The outer loop iterates
over `dim >> 3` base indices; each block produces 8 row × 8 column indices. OpenMP:
`schedule(static, 1)` on the outer loop. Threshold: `dim >= 64`.

---

## 6. Mixed State Gate Algorithms (Tiled)

### 6.1 Tile Structure

Tiles are TILE_DIM × TILE_DIM (default 32×32 = 16 KB). Full storage within each tile eliminates
per-element conjugation branches (unlike packed format). See `docs/STATE_MODULE.md` Section 4.3 for
tile layout details.

### 6.2 Within-Tile Operations (target < LOG_TILE_DIM)

When the target qubit index is less than LOG_TILE_DIM (default 5), all 4 butterfly elements lie within
the same tile.

```c
// Load rho00, rho01, rho10, rho11 from local offsets within tile
// Apply gate, write back directly — no lower_01/diag_block flags needed
```

### 6.3 Cross-Tile Operations (target ≥ LOG_TILE_DIM)

Butterfly elements span 4 different tiles. Process in tile groups:

```c
tile_bit = target - LOG_TILE_DIM;   // which bit in tile index
tile_incr = 1 << tile_bit;

for each tile group (base_tr, base_tc) where both have tile_bit = 0:
    T00 = tile(base_tr, base_tc)
    T01 = tile(base_tr, base_tc | tile_incr)
    T10 = tile(base_tr | tile_incr, base_tc)
    T11 = tile(base_tr | tile_incr, base_tc | tile_incr)
    // Some tiles may be in upper triangle → need conjugation logic at tile level
```

### 6.4 Controlled Gates on Tiled Storage

Same four-block decomposition as packed (Section 5.2). Within-tile: branch-free passes (no triangle
constraints). Cross-tile: triangle status determined once per tile group.

### 6.5 Parallelization (Tiled)

```c
// Within-tile: parallelize over tiles
#pragma omp parallel for schedule(dynamic)
for (dim_t t = 0; t < n_lower_tiles; ++t)

// Cross-tile: parallelize over tile groups
#pragma omp parallel for schedule(dynamic)
for (dim_t g = 0; g < n_tile_groups; ++g)
```

---

## 7. Unit Tests

### 7.1 Verification Approach

Gate functions are validated against reference implementations:
- **Pure states**: ψ' = Uψ via matrix-vector multiplication (`zgemv`)
- **Mixed states**: ρ' = UρU† via matrix-matrix multiplication (`zgemm`)

Reference implementations build full matrices via Kronecker products (I⊗...⊗U⊗...⊗I), limited to
n ≤ 4 qubits.

### 7.2 Test State Basis

**Single-qubit informationally complete set:** |0⟩, |1⟩, |+⟩, |-⟩, |+i⟩, |-i⟩.

**Multi-qubit states:** Tensor products of the above, plus Bell states (n=2), GHZ and W states (n≥3).

**Mixed states:** All pure states as density matrices, maximally mixed I/2ⁿ, random statistical
mixtures (10-20 per system size, deterministic seed).

### 7.3 Coverage Strategy

| Gate Type | System Sizes | Positions | Tests per Gate |
|-----------|-------------|-----------|----------------|
| Single-qubit | n=1,2,3,4 | All targets | ~11,800 |
| Two-qubit | n=2,3,4 | All (ctrl, tgt) pairs | ~37,700 |
| Three-qubit | n=3,4 | All (c1, c2, tgt) triples | ~40,000 |

**Numerical tolerance:** ε = 1e-12 (double precision with accumulated rounding over O(2ⁿ) operations).

**Runtime target:** Complete all gate tests in < 2 minutes.

### 7.4 Test Infrastructure

| Component | Location | Purpose |
|-----------|----------|---------|
| `testSingleQubitGate()` | `test/gate/test_gate_pure_1q.c` | Pure state 1Q harness |
| `testSingleQubitGateMixed()` | `test/gate/test_gate_packed_1q.c` | Mixed packed 1Q harness |
| `testRotationGate()` | `test/gate/test_gate_pure_1q.c` | Rotation pure harness |
| `testRotationGateMixed()` | `test/gate/test_gate_packed_1q.c` | Rotation mixed harness |
| `testTwoQubitGate()` | `test/gate/test_gate_pure_2q.c` | Two-qubit pure harness |
| `testTwoQubitGateMixed()` | `test/gate/test_gate_packed_2q.c` | Two-qubit mixed harness |
| `testThreeQubitGate()` | `test/gate/test_gate_pure_3q.c` | Three-qubit pure harness |
| `testThreeQubitGateMixed()` | `test/gate/test_gate_packed_3q.c` | Three-qubit mixed harness |
| `mat_single_qubit_gate()` | `test/utility/gatemat.c` | Single-qubit reference matrix |
| `mat_two_qubit_gate()` | `test/utility/gatemat.c` | Two-qubit reference matrix |
| `mat_three_qubit_gate()` | `test/utility/gatemat.c` | Three-qubit reference matrix |

See `docs/GATE_TEST_STRATEGY.md` for comprehensive test counts and implementation plan.

---

## 8. Implementation Phases

### Phase 1: Single-Qubit Gates

Gates: X, Y, Z, H, S, Sdg, T, Tdg, Hy, P, Rx, Ry, Rz.

**Prerequisites:** `GATE_VALIDATE` macro, dispatcher framework, `TRAVERSE_PURE_1Q` and
`TRAVERSE_PACKED_BLOCKS` macros, within-tile and cross-tile traversal for tiled format.

### Phase 2: Two-Qubit Gates

Gates: CX, CY, CZ, SWAP.

**Prerequisites:** `mat_two_qubit_gate()` reference builder,
`testTwoQubitGate()` and `testTwoQubitGateMixed()` harnesses.

**Note:** All four gates are implemented in-place using the two-bit insertion approach (Section 5.2).

### Phase 3: Three-Qubit Gate (Toffoli)

**Prerequisites:** `mat_three_qubit_gate()` reference builder, three-qubit test harness.
