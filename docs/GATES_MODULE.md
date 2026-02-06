# Gates Module: Technical Specification

**Module:** Quantum Gate Operations
**Header:** `include/gate.h`
**Implementation:** `src/gate_pure.c` (pure), `src/gate_packed.c` (mixed packed), `src/gate_tiled.c` (mixed tiled), `src/gate.c` (dispatchers)
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

| Gate | Function | Matrix | Category | Pure | Packed | Tiled |
|------|----------|--------|----------|------|--------|-------|
| Pauli-X | `x(state, target)` | `[[0,1],[1,0]]` | Pauli | TODO | TODO | TODO |
| Pauli-Y | `y(state, target)` | `[[0,-i],[i,0]]` | Pauli | TODO | TODO | TODO |
| Pauli-Z | `z(state, target)` | `[[1,0],[0,-1]]` | Pauli | TODO | TODO | TODO |
| Hadamard | `h(state, target)` | `[[1,1],[1,-1]]/√2` | Clifford | TODO | TODO | TODO |
| Phase S | `s(state, target)` | `[[1,0],[0,i]]` | Clifford | TODO | TODO | TODO |
| S-dagger | `sdg(state, target)` | `[[1,0],[0,-i]]` | Clifford | TODO | TODO | TODO |
| T | `t(state, target)` | `[[1,0],[0,e^(iπ/4)]]` | Non-Clifford | TODO | TODO | TODO |
| T-dagger | `tdg(state, target)` | `[[1,0],[0,e^(-iπ/4)]]` | Non-Clifford | TODO | TODO | TODO |
| Hadamard-Y | `hy(state, target)` | `[[1,-1],[1,1]]/√2` | Clifford | TODO | TODO | TODO |
| P(θ) | `p(state, target, θ)` | `[[1,0],[0,e^(iθ)]]` | Rotation | TODO | TODO | TODO |
| Rx(θ) | `rx(state, target, θ)` | `[[c,-is],[-is,c]]` | Rotation | TODO | TODO | TODO |
| Ry(θ) | `ry(state, target, θ)` | `[[c,-s],[s,c]]` | Rotation | TODO | TODO | TODO |
| Rz(θ) | `rz(state, target, θ)` | `[[e^(-iθ/2),0],[0,e^(iθ/2)]]` | Rotation | TODO | TODO | TODO |

### 2.2 Two-Qubit Gates

Two-qubit gates take `(state, control, target)` for controlled gates or `(state, q1, q2)` for symmetric
gates.

| Gate | Function | Category | Pure | Packed | Tiled |
|------|----------|----------|------|--------|-------|
| CNOT | `cx(state, ctrl, tgt)` | Controlled-Pauli | TODO | TODO | TODO |
| CY | `cy(state, ctrl, tgt)` | Controlled-Pauli | TODO | TODO | TODO |
| CZ | `cz(state, ctrl, tgt)` | Controlled-Pauli | TODO | TODO | TODO |
| CS | `cs(state, ctrl, tgt)` | Controlled-Phase | TODO | TODO | TODO |
| CSdg | `csdg(state, ctrl, tgt)` | Controlled-Phase | TODO | TODO | TODO |
| CH | `ch(state, ctrl, tgt)` | Controlled-Clifford | TODO | TODO | TODO |
| CHy | `chy(state, ctrl, tgt)` | Controlled-Clifford | TODO | TODO | TODO |
| CT | `ct(state, ctrl, tgt)` | Controlled-Non-Clifford | TODO | TODO | TODO |
| CTdg | `ctdg(state, ctrl, tgt)` | Controlled-Non-Clifford | TODO | TODO | TODO |
| CP(θ) | `cp(state, ctrl, tgt, θ)` | Controlled-Rotation | TODO | TODO | TODO |
| CPdg(θ) | `cpdg(state, ctrl, tgt, θ)` | Controlled-Rotation | TODO | TODO | TODO |
| SWAP | `swap(state, q1, q2)` | Permutation | TODO | TODO | TODO |

**CZ symmetry:** CZ(a,b) = CZ(b,a).

**RZZ(θ) — future work:** RZZ(θ) = diag(e^(-iθ/2), e^(iθ/2), e^(iθ/2), e^(-iθ/2)) is relevant for
QAOA circuits. Since Z⊗Z is diagonal, RZZ applies only phases — no off-diagonal mixing. Its packed and
tiled implementations require (row, col)-dependent phase scaling rather than block mixing. Deferred to a
future revision.

### 2.3 Three-Qubit Gates

| Gate | Function | Description | Pure | Packed | Tiled |
|------|----------|-------------|------|--------|-------|
| Toffoli | `ccx(state, ctrl1, ctrl2, tgt)` | Doubly-controlled X | TODO | TODO | TODO |

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
        case MIXED_TILED:  x_tiled(state, target);   break;
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

### 4.2 Two-Qubit Gates: TRAVERSE_PURE_2Q

For a controlled-U gate with control c and target t, enumerate all base indices where both c and t bits
are 0 using a three-level nested loop. Define `lo = min(c, t)`, `hi = max(c, t)`:

```c
#define TRAVERSE_PURE_2Q(state, control, target, PAIR_OP)
    lo = min(control, target);
    hi = max(control, target);
    stride_lo = 1 << lo;
    stride_hi = 1 << hi;
    step_lo   = 1 << (lo + 1);
    step_hi   = 1 << (hi + 1);
    incr_ctrl = 1 << control;
    incr_tgt  = 1 << target;

    for (i = 0; i < dim; i += step_hi) {           // upper bits (above hi)
        for (j = 0; j < stride_hi; j += step_lo) { // middle bits (between lo and hi)
            for (k = 0; k < stride_lo; ++k) {       // lower bits (below lo)
                base = i + j + k;
                a = data + base + incr_ctrl;             // control=1, target=0
                b = data + base + incr_ctrl + incr_tgt;  // control=1, target=1
                PAIR_OP(a, b);
            }
        }
    }
```

OpenMP: `collapse(2)` on the outer two loops.

**SWAP:** Not a controlled-U gate. Uses the same three-level stride loop but applies to a different pair:

```c
for (i, j, k as above) {
    base = i + j + k;
    idx_01 = base | incr_hi;    // lo-bit=0, hi-bit=1
    idx_10 = base | incr_lo;    // lo-bit=1, hi-bit=0
    swap(data[idx_01], data[idx_10]);
}
```

### 4.3 Three-Qubit Gates: Four-Level Stride

For Toffoli (CCX) with controls c1, c2 and target t, the stride approach extends to four nested loops.
Sort the three qubit positions: `lo < mid < hi`.

```c
for (i = 0; i < dim; i += step_hi) {
    for (j = 0; j < stride_hi; j += step_mid) {
        for (k = 0; k < stride_mid; k += step_lo) {
            for (l = 0; l < stride_lo; ++l) {       // stride-1 inner loop
                base = i + j + k + l;
                idx_c1c2_t0 = base | incr_c1 | incr_c2;
                idx_c1c2_t1 = base | incr_c1 | incr_c2 | incr_tgt;
                swap(data[idx_c1c2_t0], data[idx_c1c2_t1]);
            }
        }
    }
}
```

OpenMP: `collapse(3)` on the outer three loops. Higher-order multi-controlled gates belong to a separate
decomposition module.

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
| Rz | diagonal phases | a unchanged, b phase-rotated |

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

OpenMP: `schedule(dynamic)` on the outer `bc` loop. Threshold: `dim >= 64` (n ≥ 6 qubits).

### 5.2 Controlled Gates: One-Bit Insertion with Four-Block Decomposition

A controlled-U gate V = P₀⊗I + P₁⊗U produces four density matrix blocks under conjugation VρV†:

```
H'₀₀ = H₀₀                   (control=0 in both row and col: unchanged)
H'₀₁ = H₀₁ · U†              (control=0 in row, 1 in col: right-multiply by U†)
H'₁₀ = U · H₁₀               (control=1 in row, 0 in col: left-multiply by U)
H'₁₁ = U · H₁₁ · U†          (control=1 in both: full conjugation)
```

**Implementation:** Use the same `TRAVERSE_PACKED_BLOCKS` with `insertBit0` at the **target** position.
The control bits remain as natural values in the loop indices. The BLOCK_OP checks control bit values
and branches:

```c
int row_ctrl = (r0 >> control) & 1;
int col_ctrl = (c0 >> control) & 1;
if (!row_ctrl && !col_ctrl) { /* H_00: skip */ }
else if (row_ctrl && col_ctrl) { /* H_11: full conjugation */ }
else if (row_ctrl) { /* H_10: left-multiply by U */ }
else { /* H_01: right-multiply by U† */ }
```

**One-sided operations for specific gates:**

| Gate | H₁₁ (full conj) | H₁₀ (left U) | H₀₁ (right U†) |
|------|------------------|---------------|-----------------|
| CX | swap 00↔11, swap 01↔10 | swap rows (00↔10, 01↔11) | swap cols (00↔01, 10↔11) |
| CZ | negate off-diag (01, 10) | negate row 1 (10, 11) | negate col 1 (01, 11) |
| CY | swap + phase | swap + phase rows | swap + phase cols |
| CS | rotate off-diag by i | rotate row 1 by i | rotate col 1 by -i |
| CH | butterfly on block | butterfly on rows | butterfly on cols |

**CX and CZ are zero-flop in all four blocks** — only memory moves and sign flips.

**Multi-controlled gates (Toffoli):** The same pattern extends to multiple controls. For CCX with
controls c1, c2:

```c
dim_t ctrl_mask = (1 << ctrl1) | (1 << ctrl2);
int row_ctrl = (r0 & ctrl_mask) == ctrl_mask;
int col_ctrl = (c0 & ctrl_mask) == ctrl_mask;
```

The four-way branch is identical.

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
| `testSingleQubitGate()` | `test/src/test_gate_pure.c` | Pure state harness |
| `testSingleQubitGateMixed()` | `test/src/test_gate_packed.c` | Mixed packed harness |
| `testRotationGate()` | `test/src/test_gate_pure.c` | Rotation pure harness |
| `testRotationGateMixed()` | `test/src/test_gate_packed.c` | Rotation mixed harness |
| `testTwoQubitGate()` | `test/src/test_gate_pure.c` | Two-qubit pure harness |
| `testTwoQubitGateMixed()` | `test/src/test_gate_packed.c` | Two-qubit mixed harness |
| `mat_single_qubit_gate()` | `test/src/gatemat.c` | Single-qubit reference matrix |
| `mat_two_qubit_gate()` | `test/src/gatemat.c` | Two-qubit reference matrix |
| `mat_three_qubit_gate()` | `test/src/gatemat.c` | Three-qubit reference matrix (TODO) |

See `docs/GATE_TEST_STRATEGY.md` for comprehensive test counts and implementation plan.

---

## 8. Implementation Phases

### Phase 1: Single-Qubit Gates

Gates: X, Y, Z, H, S, Sdg, T, Tdg, Hy, P, Rx, Ry, Rz.

**Prerequisites:** `GATE_VALIDATE` macro, dispatcher framework, `TRAVERSE_PURE_1Q` and
`TRAVERSE_PACKED_BLOCKS` macros, within-tile and cross-tile traversal for tiled format.

### Phase 2: Two-Qubit Gates

Gates: CX, CY, CZ, CS, CSdg, CH, CHy, CT, CTdg, CP, CPdg, SWAP.

**Prerequisites:** `TRAVERSE_PURE_2Q` macro, `mat_two_qubit_gate()` reference builder,
`testTwoQubitGate()` and `testTwoQubitGateMixed()` harnesses.

**Action required:** Replace `cx_mixed()` output-buffer approach with in-place four-block decomposition
(Section 5.2) before implementing remaining controlled gates.

### Phase 3: Three-Qubit Gate (Toffoli)

**Prerequisites:** `mat_three_qubit_gate()` reference builder, three-qubit test harness.
