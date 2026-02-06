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
- **Stride-based traversal**: nested loops with power-of-two strides expose contiguous inner-loop access
  to the compiler, enabling auto-vectorization and hardware prefetching
- **Three-way dispatch**: single `x(state, target)` function dispatches to `x_pure()`, `x_packed()`, or
  `x_tiled()` based on `state->type` enum
- **Direct packed operations**: mixed state gates operate directly on lower-triangular packed or tiled
  storage (no unpacking required)
- **Macro-based traversal**: common loop structures factored into traversal macros with gate-specific
  operations passed as macro parameters
- **In-place**: all gates modify state data in-place with O(1) extra memory (exception: SWAP on packed
  storage, see Section 5.3)
- **Test-driven validation**: reference implementations using full BLAS matrix operations verify correctness

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

| Gate | Function | Description | Category | Pure | Packed | Tiled |
|------|----------|-------------|----------|------|--------|-------|
| CNOT | `cx(state, ctrl, tgt)` | Controlled-X | Controlled-Pauli | TODO | TODO | TODO |
| CY | `cy(state, ctrl, tgt)` | Controlled-Y | Controlled-Pauli | TODO | TODO | TODO |
| CZ | `cz(state, ctrl, tgt)` | Controlled-Z | Controlled-Pauli | TODO | TODO | TODO |
| CS | `cs(state, ctrl, tgt)` | Controlled-S | Controlled-Phase | TODO | TODO | TODO |
| CSdg | `csdg(state, ctrl, tgt)` | Controlled-S† | Controlled-Phase | TODO | TODO | TODO |
| CH | `ch(state, ctrl, tgt)` | Controlled-H | Controlled-Clifford | TODO | TODO | TODO |
| CHy | `chy(state, ctrl, tgt)` | Controlled-Hy | Controlled-Clifford | TODO | TODO | TODO |
| CT | `ct(state, ctrl, tgt)` | Controlled-T | Controlled-Non-Clifford | TODO | TODO | TODO |
| CTdg | `ctdg(state, ctrl, tgt)` | Controlled-T† | Controlled-Non-Clifford | TODO | TODO | TODO |
| CP(θ) | `cp(state, ctrl, tgt, θ)` | Controlled-P(θ) | Controlled-Rotation | TODO | TODO | TODO |
| CPdg(θ) | `cpdg(state, ctrl, tgt, θ)` | Controlled-P†(θ) | Controlled-Rotation | TODO | TODO | TODO |
| SWAP | `swap(state, q1, q2)` | Qubit swap | Permutation | TODO | TODO | TODO |

**CZ symmetry:** CZ is symmetric in control and target: CZ(a,b) = CZ(b,a). The implementation may exploit
this but the API takes `(control, target)` for consistency.

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

Internal implementations (`x_pure`, `x_packed`, `x_tiled`) trust validated inputs and perform no
additional checks.

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

**Rationale:** Gate functions may be called directly by user code, where a bare `assert()` would be
disabled by `-DNDEBUG` in release builds, silently proceeding with invalid state. `GATE_VALIDATE` fires
unconditionally and provides a diagnostic message. This is stricter than assert (always active) but more
informative (names the failing condition).

**Validated conditions:**

| Condition | Message |
|-----------|---------|
| `state && state->data` | null state or data pointer |
| `target < state->qubits` | target qubit out of range |
| `control < state->qubits` | control qubit out of range |
| `control != target` | control and target must differ |

Gate functions return `void`. There are no error codes.

### 3.3 Two-Qubit Validation

Two-qubit dispatchers additionally validate distinct qubit indices:

```c
GATE_VALIDATE(control < state->qubits, "cx: control qubit out of range");
GATE_VALIDATE(target < state->qubits, "cx: target qubit out of range");
GATE_VALIDATE(control != target, "cx: control and target must differ");
```

### 3.4 Naming Convention

Public API uses short names: `x()`, `cx()`, `rx()`. Internal implementations use suffixed names:
`x_pure()`, `x_packed()`, `x_tiled()`.

### 3.5 File Organization

| File | Contents |
|------|----------|
| `include/gate.h` | Public API declarations with matrix documentation |
| `src/gate.c` | Dispatchers with input validation |
| `src/gate_pure.c` | Pure state implementations (`*_pure`) |
| `src/gate_packed.c` | Mixed packed implementations (`*_packed`) |
| `src/gate_tiled.c` | Mixed tiled implementations (`*_tiled`) |

---

## 4. Pure State Gate Algorithms

### 4.1 Single-Qubit Gates: TRAVERSE_PURE_1Q

For qubit q in an n-qubit state, pairs (ψ[idx], ψ[idx + 2^q]) are processed independently. The stride-
based nested loop makes the contiguous inner-loop access pattern explicit and compiler-visible, enabling
auto-vectorization and hardware prefetching.

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

The inner loop increments by 1 (stride-1 access), while the outer loop steps over pairs of contiguous
blocks. `collapse(2)` gives dim/2 total iterations regardless of target position:
- target=0 (stride=1): outer loop has dim/2 iterations, inner has 1
- target=n-1 (stride=dim/2): outer loop has 1 iteration, inner has dim/2
- All cases: dim/2 total units of work, evenly distributable by OpenMP

**Why stride-based is better than `insertBit0`:** Both produce identical index sequences. But the stride
loop exposes the contiguous access pattern to the compiler: `data[i + j]` with `j++` is trivially
vectorizable. `insertBit0(k, target)` computes each address via bit manipulation — the compiler cannot
prove consecutive `k` values yield consecutive addresses, so it falls back to scalar code.

OpenMP threshold: `dim >= 4096` (n ≥ 12 qubits).

### 4.2 Two-Qubit Gates: TRAVERSE_PURE_2Q

For a controlled-U gate with control c and target t, enumerate all base indices where both c and t bits
are 0 using a three-level nested loop. The same stride-based approach extends naturally from the 1Q case.

Define `lo = min(c, t)`, `hi = max(c, t)`:

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
            for (k = 0; k < stride_lo; ++k) {       // lower bits (below lo), stride-1
                base = i + j + k;
                a = data + base + incr_ctrl;             // control=1, target=0
                b = data + base + incr_ctrl + incr_tgt;  // control=1, target=1
                PAIR_OP(a, b);
            }
        }
    }
```

**Proof of correctness:** The three loops enumerate exactly dim/4 base values with bits `lo` and `hi`
both equal to 0. The outer loop covers bits above `hi` (2^(n-hi-1) iterations), the middle loop covers
bits between `lo` and `hi` (2^(hi-lo-1) iterations), and the inner loop covers bits below `lo` (2^lo
iterations). Total: 2^(n-2) = dim/4.

**Memory access patterns:**
- The innermost loop is always stride-1, just like `TRAVERSE_PURE_1Q`
- When `target = lo` (target < control): pair elements are separated by `stride_lo` — the close case
- When `target = hi` (target > control): pair elements are separated by `stride_hi` — far apart, but
  each stream is individually contiguous

**OpenMP:** `collapse(2)` on the outer two loops distributes dim/(4·stride_lo) chunks, each running the
stride-1 inner loop of `stride_lo` iterations. This keeps the vectorizable inner loop intact while giving
OpenMP enough parallel units.

For CX (U=X), PAIR_OP is a swap. For CZ (U=Z), PAIR_OP negates `*b`. For CY (U=Y), PAIR_OP applies
the Y gate formula. Any controlled-U gate reuses this traversal with the appropriate single-qubit PAIR_OP.

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

Total iterations: dim/8. OpenMP `collapse(3)` on the outer three loops. Toffoli is the only three-qubit
gate in this module; higher-order multi-controlled gates belong to a separate decomposition module.

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

**Key helper: `insertBit0(val, pos)`** inserts a 0-bit at position `pos`, enabling direct enumeration
of all indices with the target bit = 0. No wasted iterations.

**BLOCK_OP receives:**
- `data`: pointer to state data array
- `idx00, idx01, idx10, idx11`: packed indices for the 4 block elements
- `lower_01`: 1 if (r0, c1) is in the lower triangle, 0 if via conjugate
- `diag_block`: 1 if on the diagonal (bc == br), where idx01 == idx10

OpenMP: `schedule(dynamic)` on the outer `bc` loop. Threshold: `dim >= 64` (n ≥ 6 qubits).

### 5.2 Controlled Gates: One-Bit Insertion with Four-Block Decomposition

A controlled-U gate V = P₀⊗I + P₁⊗U produces four density matrix blocks under conjugation VρV†:

```
H'₀₀ = H₀₀                   (control=0 in both row and col: unchanged)
H'₀₁ = H₀₁ · U†              (control=0 in row, 1 in col: right-multiply by U†)
H'₁₀ = U · H₁₀               (control=1 in row, 0 in col: left-multiply by U)
H'₁₁ = U · H₁₁ · U†          (control=1 in both: full conjugation)
```

**Cross-block modification is required.** Only H₀₀ is left unchanged. The cross-blocks H₀₁ and H₁₀
undergo one-sided multiplication and must not be skipped.

**Left-multiply by U (H₁₀ block):** Each column of the 2×2 sub-block is transformed by U:

```
ρ'₀₀ = U₀₀·ρ₀₀ + U₀₁·ρ₁₀
ρ'₀₁ = U₀₀·ρ₀₁ + U₀₁·ρ₁₁
ρ'₁₀ = U₁₀·ρ₀₀ + U₁₁·ρ₁₀
ρ'₁₁ = U₁₀·ρ₀₁ + U₁₁·ρ₁₁
```

**Right-multiply by U† (H₀₁ block):** Each row of the 2×2 sub-block is transformed by U†:

```
ρ'₀₀ = ρ₀₀·conj(U₀₀) + ρ₀₁·conj(U₁₀)
ρ'₀₁ = ρ₀₀·conj(U₀₁) + ρ₀₁·conj(U₁₁)
ρ'₁₀ = ρ₁₀·conj(U₀₀) + ρ₁₁·conj(U₁₀)
ρ'₁₁ = ρ₁₀·conj(U₀₁) + ρ₁₁·conj(U₁₁)
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

**Why one-bit insertion over two-bit insertion for packed format:** Two-bit insertion would eliminate the
H₀₀ skip iterations (~25% of total) and remove the 4-way branch. However, in the packed format:
- The 4-way branch cost is ~1-3% of total per-iteration cost (2-4 cycles branch vs 50-200 cycles memory)
- The control bit flips slowly as `br` increments, giving the branch predictor long predictable runs
- Two-bit insertion would require three separate traversal passes (H₁₁, H₁₀, H₀₁), each with its own
  OpenMP barrier and triangular constraint logic
- Extending to multi-controlled gates (Toffoli) with two-bit insertion requires combinatorial passes;
  one-bit insertion extends trivially via a mask check

**Overhead:** 4 integer operations per block (shift + AND for row and col). Negligible vs memory access
cost. Compensated by skipping 1/4 of all blocks (H₀₀).

**One-sided operations for specific gates:**

| Gate | H₁₁ (full conj) | H₁₀ (left U) | H₀₁ (right U†) |
|------|------------------|---------------|-----------------|
| CX | swap 00↔11, swap 01↔10 | swap rows (00↔10, 01↔11) | swap cols (00↔01, 10↔11) |
| CZ | negate off-diag (01, 10) | negate row 1 (10, 11) | negate col 1 (01, 11) |
| CY | swap + phase | swap + phase rows | swap + phase cols |
| CS | rotate off-diag by i | rotate row 1 by i | rotate col 1 by -i |
| CH | butterfly on block | butterfly on rows | butterfly on cols |

**CX and CZ are zero-flop in all four blocks** — only memory moves and sign flips.

**Memory:** O(1) extra memory. All operations are in-place on the 4 block elements.

**Multi-controlled gates (Toffoli):** The same pattern extends to multiple controls. For CCX with
controls c1, c2:

```c
dim_t ctrl_mask = (1 << ctrl1) | (1 << ctrl2);
int row_ctrl = (r0 & ctrl_mask) == ctrl_mask;
int col_ctrl = (c0 & ctrl_mask) == ctrl_mask;
```

The four-way branch is identical. CCX is zero-flop in all blocks (same as CX).

### 5.3 SWAP (Packed)

SWAP is a permutation gate: perm(x) = x with bits q1 and q2 exchanged. Like CX, this is zero-flop
(only memory moves). Unlike CX, SWAP is not a controlled-U gate and does not decompose into four blocks.

**The packed format complicates in-place SWAP.** For a 4×4 block of the density matrix (indexed by the
two qubit values in row and column), SWAP permutes rows and columns by exchanging the A↔B indices.
This produces 6 element swaps per block. On diagonal base blocks (r_base == c_base), Hermitian symmetry
causes some swaps to alias the same packed location, reducing to 3 swaps + 1 in-place conjugation. But
on off-diagonal base blocks, the lower-triangle constraint makes the triangle status of each element
data-dependent, introducing up to 12 branches per block and aliasing hazards where different swaps map
to the same packed element through conjugate symmetry.

**Three implementation options:**

1. **CX decomposition (recommended):** SWAP = CX(a,b) · CX(b,a) · CX(a,b). Uses the in-place four-block
   CX implementation three times. O(1) extra memory, correct by construction, reuses existing
   infrastructure. Cost: 3× the memory traffic of a direct implementation.

2. **Output-buffer approach:** Allocate O(dim²) temporary, iterate over all lower-triangle elements
   computing perm(r) and perm(c), gather from source, write to destination. Simple (12 lines of core
   logic), no aliasing, embarrassingly parallel. Cost: O(dim²) temporary allocation.

3. **In-place two-bit insertion:** Enumerate 4×4 blocks via `insertBits2_0`, perform swaps with
   conditional conjugation. Requires loading all block values before writing to avoid aliasing. Complex
   (~70 lines), high bug risk. Not recommended for packed format.

### 5.4 Per-Gate Optimizations (Mixed Packed)

| Gate | Operation | Optimization |
|------|-----------|-------------|
| X | swap blocks | Zero multiplications (swaps only) |
| Y | swap + negate | Zero multiplications |
| Z | negate off-diag | Zero multiplications |
| S | (1,0)*=i, (0,1)*=-i | Zero multiplications (component swap) |
| Sdg | (1,0)*=-i, (0,1)*=i | Zero multiplications (component swap) |
| T | phase rotation | 2 muls per off-diag element |
| Tdg | phase rotation | 2 muls per off-diag element |
| H | butterfly | 4 adds + 2 muls per block |
| Rx | 2×2 block mixing | ics terms via component arithmetic |
| Ry | 2×2 block mixing | Real coefficients cs, c², s² precomputed |
| Rz | off-diag phases | Diagonal unchanged, cos/sin for phase rotation |

---

## 6. Mixed State Gate Algorithms (Tiled)

### 6.1 Tile Structure

Tiles are TILE_DIM × TILE_DIM (default 32×32 = 16 KB). See `docs/STATE_MODULE.md` Section 4.3 for
index formulas. Key properties:

- **Lower-triangular at tile level**: only tiles T(tr, tc) where tr ≥ tc are stored
- **Row-major within tiles**: enables SIMD-friendly sequential access
- **Full storage within tiles**: both upper and lower triangle of the matrix block are stored within
  each tile, unlike packed format

### 6.2 Within-Tile Operations (target < LOG_TILE_DIM)

When the target qubit index is less than LOG_TILE_DIM (default 5), all 4 butterfly elements lie within
the same tile. **No conjugation logic is needed:**

- **Diagonal tiles** (tr == tc): store the full TILE_DIM × TILE_DIM region including both lower and
  upper triangle portions. All 4 butterfly elements are directly accessible.
- **Off-diagonal tiles** (tr > tc): for any butterfly anchored at (r0, c0) within this tile,
  r0 ≥ tr × TILE_DIM ≥ (tc + 1) × TILE_DIM > c1, so all 4 elements are in the lower triangle.

**Cache performance:** 1 tile (16 KB) fits in L1 cache (~32–64 KB). Operations on the rightmost
LOG_TILE_DIM qubits stay entirely within L1, avoiding cache misses.

```c
// Within-tile kernel: NO conjugation branches
static inline void apply_butterfly_within_tile(
    cplx_t *tile, dim_t lr0, dim_t lc0, dim_t stride
) {
    dim_t lr1 = lr0 + stride, lc1 = lc0 + stride;
    cplx_t rho00 = tile[lr0 * TILE_DIM + lc0];
    cplx_t rho01 = tile[lr0 * TILE_DIM + lc1];
    cplx_t rho10 = tile[lr1 * TILE_DIM + lc0];
    cplx_t rho11 = tile[lr1 * TILE_DIM + lc1];
    // Apply gate, write back directly — no lower_01/diag_block flags needed
}
```

**Advantages over packed format:**
- No `lower_01` flag
- No `diag_block` special cases
- No conditional conjugation on read/write
- Cleaner auto-vectorization

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

**Cache performance:** 4 tiles (64 KB) fit in L2 cache (~128–256 KB). Operations on higher qubit
indices access 4 tiles simultaneously but remain within L2.

| Case | Working Set | Cache Level | Conjugation |
|------|-------------|-------------|-------------|
| Within-tile (target < LOG_TILE_DIM) | 16 KB (1 tile) | L1 | None |
| Cross-tile (target ≥ LOG_TILE_DIM) | 64 KB (4 tiles) | L2 | May be needed at tile level |

### 6.4 Controlled Gates on Tiled Storage

Same four-block decomposition as packed (Section 5.2). For within-tile controlled gates (both target and
control < LOG_TILE_DIM), no conjugation is needed — the gate branches only on control bit values.

**Two-bit insertion for tiled controlled gates:** Unlike the packed format, the tiled format eliminates
`lower_01`/`diag_block` complications within tiles. This makes two-bit insertion a viable strategy for
tiled controlled gates: insert 0-bits at both target and control positions, then enumerate the three
non-trivial block types (H₁₁, H₁₀, H₀₁) in separate branch-free passes. Each pass performs a uniform
operation on every iteration, enabling better instruction pipelining than the 4-way branch.

For within-tile operations (both qubits < LOG_TILE_DIM), all four butterfly elements lie within the same
tile with full row-major storage — no triangle constraints, no conjugation, no aliasing. This is exactly
the regime where two-bit insertion pays off.

For cross-tile operations, the analysis depends on how many of the four qubit-indexed tiles lie in the
stored lower triangle of the tile grid. The tile-level triangle constraint is simpler than the packed
element-level constraint: it is determined once per tile group, not per element.

### 6.5 Parallelization (Tiled)

```c
// Within-tile: parallelize over tiles
#pragma omp parallel for schedule(dynamic)
for (dim_t t = 0; t < n_lower_tiles; ++t)

// Cross-tile: parallelize over tile groups
#pragma omp parallel for schedule(dynamic)
for (dim_t g = 0; g < n_tile_groups; ++g)
```

For n=8: 36 tiles → 36 parallel units (within-tile) or 9 groups (cross-tile).

---

## 7. Compiler Optimization and SIMD

Gate operations are memory-bandwidth bound. Performance depends on auto-vectorization.

**Recommended flags:**

| Flag | Purpose |
|------|---------|
| `-O3` | Aggressive optimizations including auto-vectorization |
| `-march=native` | Best SIMD instructions for current CPU (AVX2, AVX-512, etc.) |
| `-ffast-math` | Relax IEEE float semantics (enables reordering for vectorization) |
| `-funroll-loops` | Unroll small loops for better pipelining |

**CMake release configuration:**
```cmake
set(CMAKE_C_FLAGS_RELEASE "-O3 -march=native -ffast-math -DNDEBUG")
```

**Patterns that vectorize well:** The stride-based inner loops (`TRAVERSE_PURE_1Q`, `TRAVERSE_PURE_2Q`)
have stride-1 access, no loop-carried dependencies, and simple operations (add, multiply, negate).

**Tiled row-major rationale:** Elements within tiles are row-major for SIMD-friendly sequential access
during gate kernels. Do not pass tile pointers to LAPACK routines which expect column-major.

---

## 8. OpenMP Parallelization

**Thresholds:** Different thresholds account for different work-per-gate:
- **Pure states**: `OMP_THRESHOLD = 4096` → n ≥ 12 qubits (O(dim) work)
- **Mixed states**: `OMP_THRESHOLD = 64` → n ≥ 6 qubits (O(dim²) work)

**Pure state strategy (1Q):** `collapse(2)` flattens the two-level loop, giving dim/2 total iterations
evenly distributed regardless of target qubit position.

**Pure state strategy (2Q):** `collapse(2)` on the outer two of three loops. The stride-1 inner loop
remains intact as a vectorizable unit within each parallel chunk. Total parallel granularity:
dim/(4·stride_lo) chunks of stride_lo contiguous iterations each.

**Pure state strategy (3Q):** `collapse(3)` on the outer three of four loops. Same principle as 2Q.

**Mixed packed strategy:** `schedule(dynamic)` on the outer `bc` loop of `TRAVERSE_PACKED_BLOCKS`.
Dynamic scheduling balances load when inner loop iteration count varies with `bc`.

**Mixed tiled strategy:** Parallelize over tiles (within-tile case) or tile groups (cross-tile case)
with dynamic scheduling.

**Build configuration:**
```bash
cmake -DENABLE_OPENMP=ON ..   # Enable (default)
cmake -DENABLE_OPENMP=OFF ..  # Disable for single-threaded execution
```

---

## 9. Thread Safety

State objects are not thread-safe for concurrent mutation. Users must create separate `state_t`
objects for parallel evaluation. Gate module achieves within-gate thread safety through disjoint
memory regions in the butterfly decomposition.

---

## 10. Unit Tests

### 10.1 Verification Approach

Gate functions are validated against reference implementations:
- **Pure states**: ψ' = Uψ via matrix-vector multiplication (`zgemv`)
- **Mixed states**: ρ' = UρU† via matrix-matrix multiplication (`zgemm`)

Reference implementations build full matrices via Kronecker products (I⊗...⊗U⊗...⊗I), limited to
n ≤ 4 qubits. Production code uses sparse stride-based indexing validated by testing all qubit positions.

### 10.2 Test State Basis

**Single-qubit informationally complete set (6 states):**
- Computational: |0⟩, |1⟩ — detects amplitude errors
- Hadamard: |+⟩, |-⟩ — detects real phase errors
- Circular: |+i⟩, |-i⟩ — detects complex phase errors

**Multi-qubit states:** Tensor products of the above, plus Bell states (n=2), GHZ and W states (n≥3).

**Mixed states:** All pure states as density matrices, maximally mixed I/2ⁿ, random statistical
mixtures (10-20 per system size, deterministic seed).

### 10.3 Coverage Strategy

| Gate Type | System Sizes | Positions | Tests per Gate |
|-----------|-------------|-----------|----------------|
| Single-qubit | n=1,2,3,4 | All targets | ~11,800 |
| Two-qubit | n=2,3,4 | All (ctrl, tgt) pairs | ~37,700 |
| Three-qubit | n=3,4 | All (c1, c2, tgt) triples | ~40,000 |

Each gate is tested on all three storage formats: PURE, MIXED_PACKED, MIXED_TILED.

**Numerical tolerance:** ε = 1e-12 (double precision with accumulated rounding over O(2ⁿ) operations).

**Runtime target:** Complete all gate tests in < 2 minutes.

### 10.4 Test Infrastructure

| Component | Location | Purpose |
|-----------|----------|---------|
| `testSingleQubitGate()` | `test/src/test_gate_pure.c` | Pure state harness |
| `testSingleQubitGateMixed()` | `test/src/test_gate_packed.c` | Mixed packed harness |
| `testRotationGate()` | `test/src/test_gate_pure.c` | Rotation pure harness |
| `testRotationGateMixed()` | `test/src/test_gate_packed.c` | Rotation mixed harness |
| `testTwoQubitGate()` | `test/src/test_gate_pure.c` | Two-qubit pure harness |
| `testTwoQubitGateMixed()` | `test/src/test_gate_packed.c` | Two-qubit mixed harness |
| `mat_single_qubit_gate()` | `test/src/gatemat.c` | Build I⊗...⊗U⊗...⊗I |
| `mat_two_qubit_gate()` | `test/src/gatemat.c` | Build full two-qubit matrix |
| `mat_three_qubit_gate()` | `test/src/gatemat.c` | Build full three-qubit matrix (TODO) |

See `docs/GATE_TEST_STRATEGY.md` for comprehensive test counts and implementation plan.

---

## 11. Implementation Phases

Each gate is implemented and tested across all three storage formats before proceeding to the next gate.
Within each phase, the order per gate is: pure → packed → tiled, with all tests passing before moving on.

### Phase 1: Single-Qubit Gates

For each gate (X, Y, Z, H, S, Sdg, T, Tdg, Hy, P, Rx, Ry, Rz):
1. Implement `*_pure` in `gate_pure.c`, test against `zgemv` reference
2. Implement `*_packed` in `gate_packed.c`, test against `zgemm` reference
3. Implement `*_tiled` in `gate_tiled.c`, test against `zgemm` reference

**Prerequisites:** `GATE_VALIDATE` macro, dispatcher framework, `TRAVERSE_PURE_1Q` and
`TRAVERSE_PACKED_BLOCKS` macros, within-tile and cross-tile traversal for tiled format.

### Phase 2: Two-Qubit Gates

For each gate (CX, CY, CZ, CS, CSdg, CH, CHy, CT, CTdg, CP, CPdg, SWAP):
1. Implement `*_pure` using `TRAVERSE_PURE_2Q`, test against `zgemv` reference
2. Implement `*_packed` using one-bit insertion with four-block decomposition, test against `zgemm`
3. Implement `*_tiled` using two-bit insertion (within-tile) or appropriate cross-tile strategy, test

**Prerequisites:** `TRAVERSE_PURE_2Q` macro, `mat_two_qubit_gate()` reference builder,
`testTwoQubitGate()` and `testTwoQubitGateMixed()` harnesses.

**Action required:** Replace `cx_mixed()` output-buffer approach with in-place four-block decomposition
(Section 5.2) before implementing remaining controlled gates.

### Phase 3: Three-Qubit Gate (Toffoli)

1. Implement `ccx_pure` using four-level stride, test against `zgemv` reference
2. Implement `ccx_packed` using one-bit insertion with multi-control mask, test against `zgemm`
3. Implement `ccx_tiled`, test against `zgemm`

**Prerequisites:** `mat_three_qubit_gate()` reference builder, three-qubit test harness.
