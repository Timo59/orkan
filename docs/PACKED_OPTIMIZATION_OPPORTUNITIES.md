# Packed Storage Optimization Opportunities

Research notes from a comparative analysis of Qiskit-Aer's flat-array `DensityMatrix<double>`
implementation vs qlib's packed lower-triangular storage. These optimizations target the
small-qubit regime (n ≤ 6) where Aer's simpler indexing gives it a speed advantage for
permutation gates (X, CX, SWAP). They are not thesis-critical — qlib already wins decisively
at n ≥ 5 — but may be worth revisiting if small-n performance becomes important.

**Date:** 2026-03-15
**Status:** Research only — no code changes made.

---

## Background: How the Two Systems Work

**Aer** stores the density matrix as a flat `dim²`-element array. It treats ρ as a `2n`-qubit
state vector via column-stacking: `vec[r + c*dim] = ρ[r,c]`. A gate on qubit `q` becomes an
operation on qubits `{q, q+n}` of this 2n-qubit vector. The iteration is:

```cpp
// Aer's apply_x (densitymatrix.hpp:395-404)
auto lambda = [&](const areg_t<4> &inds) -> void {
    std::swap(data_[inds[0]], data_[inds[3]]);   // swap ρ00 ↔ ρ11
    std::swap(data_[inds[1]], data_[inds[2]]);   // swap ρ01 ↔ ρ10
};
const areg_t<2> qubits = {{qubit, qubit + num_qubits()}};
apply_lambda(lambda, qubits);  // iterates k = 0..dim²/4, computes 4 indices per k
```

**qlib** stores only the lower triangle of ρ (Hermiticity: ρ[r,c] = conj(ρ[c,r])), saving
~50% memory. The cost is that every element access must determine whether it's in the lower
triangle (direct access) or upper triangle (read the transposed element and conjugate):

```c
// qlib's x_packed (gate_packed_1q.c:67-150)
TRAVERSE_PACKED_BLOCKS(state, target, X_BLOCK_OP);
```

which expands to a nested loop with `insertBit0`, `pack_idx`, `lower_01` checks, and
`diag_block` checks.

---

## Finding 1: Precomputed Index Tables for Small n

**Estimated impact:** 30-50% for permutation gates at 2-5 qubits.
**Complexity:** Low.

### The problem

At n=2 (dim=4, packed array = 10 elements), the X gate has only 3 block iterations. Each
iteration computes `insertBit0` twice, `pack_idx` once or twice, evaluates `lower_01`,
evaluates `diag_block`, then does 2 swaps. The index overhead per swap is enormous relative
to the swap itself.

### Aer's approach

Aer's `indexes()` function precomputes all 4 element indices for a 2-qubit block in one shot,
then the lambda just does `swap(data[inds[0]], data[inds[3]])`. The indices are computed via
`index0()` (bit insertion) plus bit-ORs — no division, no triangle checks.

### The idea

For small qubit counts (n ≤ 6), precompute ALL swap pair indices into a stack-allocated array
before entering the hot loop. Then the hot loop becomes a trivial swap sequence with zero
index computation.

### Current code (X gate hot path)

```c
// gate_packed_1q.c — the full TRAVERSE_PACKED_BLOCKS + X_BLOCK_OP
for (gate_idx_t bc = 0; bc < half_dim; ++bc) {
    gate_idx_t c0 = insertBit0(bc, target);
    gate_idx_t c1 = c0 | incr;
    gate_idx_t col1_base = c1 * (2 * dim - c1 + 1) / 2;   // <-- multiply
    gate_idx_t col0_base = c0 * (2 * dim - c0 + 1) / 2;   // <-- multiply

    for (gate_idx_t br = bc; br < half_dim; ++br) {
        gate_idx_t r0 = insertBit0(br, target);
        gate_idx_t r1 = r0 | incr;

        gate_idx_t idx00 = col0_base + (r0 - c0);
        gate_idx_t idx10 = idx00 + incr;
        gate_idx_t idx11 = col1_base + (r1 - c1);

        int lower_01 = (r0 >= c1);                          // <-- branch
        gate_idx_t idx01 = lower_01 ? (col1_base + (r0 - c1))
                                    : pack_idx(dim, c1, r0); // <-- pack_idx if upper

        int diag_block = (bc == br);                         // <-- branch

        // X_BLOCK_OP: 2 swaps with conditional conj()
        cplx_t tmp = data[idx00]; data[idx00] = data[idx11]; data[idx11] = tmp;
        if (diag_block) {
            data[idx10] = conj(data[idx10]);
        } else {
            cplx_t rho01 = lower_01 ? data[idx01] : conj(data[idx01]);
            cplx_t rho10 = data[idx10];
            data[idx10] = rho01;
            data[idx01] = lower_01 ? rho10 : conj(rho10);
        }
    }
}
```

### Proposed code (precomputed index table for small n)

```c
void x_packed(state_t *state, const qubit_t target) {
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    const gate_idx_t half_dim = dim >> 1;

    /* Small-n fast path: precompute all swap operations */
    if (state->qubits <= 6) {
        const gate_idx_t incr = (gate_idx_t)1 << target;
        cplx_t * restrict data = state->data;

        /* Maximum swap pairs for n=6: half_dim=32, ~528 blocks, ~1056 swap ops.
         * Stack-allocate a table of (idxA, idxB, flags) triples. */
        struct swap_op { gate_idx_t a, b; int do_conj; };
        struct swap_op ops[2048];  /* generous upper bound for n<=6 */
        int n_ops = 0;

        for (gate_idx_t bc = 0; bc < half_dim; ++bc) {
            gate_idx_t c0 = insertBit0(bc, target);
            gate_idx_t c1 = c0 | incr;
            gate_idx_t col0_base = c0 * (2 * dim - c0 + 1) / 2;
            gate_idx_t col1_base = c1 * (2 * dim - c1 + 1) / 2;

            for (gate_idx_t br = bc; br < half_dim; ++br) {
                gate_idx_t r0 = insertBit0(br, target);
                gate_idx_t r1 = r0 | incr;

                /* (0,0) <-> (1,1): always a plain swap */
                ops[n_ops++] = (struct swap_op){
                    col0_base + (r0 - c0),    /* idx00 */
                    col1_base + (r1 - c1),    /* idx11 */
                    0                          /* no conj */
                };

                /* (0,1) <-> (1,0): depends on triangle position */
                gate_idx_t idx10 = col0_base + (r0 - c0) + incr;
                int lower_01 = (r0 >= c1);
                gate_idx_t idx01 = lower_01 ? (col1_base + (r0 - c1))
                                            : pack_idx(dim, c1, r0);

                if (bc == br) {
                    /* Diagonal block: idx01 == idx10, just conjugate in place */
                    ops[n_ops++] = (struct swap_op){ idx10, idx10, 2 }; /* flag=2: conj-in-place */
                } else {
                    ops[n_ops++] = (struct swap_op){ idx01, idx10, lower_01 ? 0 : 1 };
                }
            }
        }

        /* Hot loop: execute precomputed swaps — no index math, no branches */
        for (int i = 0; i < n_ops; ++i) {
            gate_idx_t a = ops[i].a, b = ops[i].b;
            if (ops[i].do_conj == 2) {
                data[a] = conj(data[a]);        /* diagonal block conj-in-place */
            } else if (ops[i].do_conj == 0) {
                cplx_t tmp = data[a]; data[a] = data[b]; data[b] = tmp;
            } else {
                /* conj swap: stored value is conj of logical value */
                cplx_t tmp = conj(data[a]);
                data[a] = conj(data[b]);
                data[b] = tmp;
            }
        }
        return;
    }

    /* Large-n path: existing TRAVERSE_PACKED_BLOCKS */
    TRAVERSE_PACKED_BLOCKS(state, target, X_BLOCK_OP);
}
```

### Why this helps

The precomputation loop runs once. The hot loop is a flat array of swap operations with no
`insertBit0`, no `pack_idx`, no `lower_01` branch evaluation. At n=2, you get ~6 precomputed
ops executed in a tight loop. The branch predictor sees a simple sequential pattern. This
mirrors what Aer does — its `indexes()` function precomputes, then the lambda executes.

---

## Finding 2: Specialized X Without TRAVERSE_PACKED_BLOCKS

**Estimated impact:** 20-40% for X at 2-5 qubits.
**Complexity:** Medium.

### The problem

`TRAVERSE_PACKED_BLOCKS` is a generic macro designed to work for ALL 1Q gates (X, H, Z, Rx,
Ry, Rz, S, T, etc.). It always computes all 4 packed indices, always evaluates `lower_01`
and `diag_block`, and passes everything to the BLOCK_OP. For X (and other permutation gates),
most of this information is wasted.

### The key insight

For X on qubit `t`, the permutation `perm(x) = x XOR (1<<t)` means that within a single
packed column `c`, elements at row `r` (bit t=0) pair with elements at row `r | (1<<t)` (bit
t=1). These elements are exactly `incr = 1<<t` positions apart in the packed array within the
same column.

### Proposed approach (column-stride iteration)

```c
void x_packed_fast(state_t *state, const qubit_t target) {
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    const gate_idx_t incr = (gate_idx_t)1 << target;
    const gate_idx_t half_dim = dim >> 1;
    cplx_t * restrict data = state->data;

    /* Phase 1: Within-column swaps.
     * For each column c in 0..dim-1, elements at rows with bit t=0
     * pair with elements at rows with bit t=1.
     * In packed storage, column c starts at col_base = c*(2*dim-c+1)/2
     * and has (dim - c) elements at rows c, c+1, ..., dim-1.
     * For rows r0 (bit t=0) and r1=r0|incr (bit t=1), both in column c,
     * their packed positions are col_base+(r0-c) and col_base+(r1-c),
     * separated by exactly `incr` positions. */

    gate_idx_t col_base = 0;  /* incremental: base(c+1) = base(c) + (dim-c) */

    for (gate_idx_t c = 0; c < dim; ++c) {
        gate_idx_t col_len = dim - c;  /* number of rows in this column */

        for (gate_idx_t local_br = 0; local_br < col_len; ) {
            gate_idx_t r = c + local_br;
            if ((r & incr) == 0 && (r | incr) < dim) {
                /* r has bit t=0; partner is r|incr, offset = incr within column */
                gate_idx_t pos0 = col_base + local_br;
                gate_idx_t pos1 = col_base + local_br + incr;
                /* Both are in column c, so this is a same-column swap.
                 * Same column in lower triangle = no conjugation needed. */
                cplx_t tmp = data[pos0];
                data[pos0] = data[pos1];
                data[pos1] = tmp;
                local_br += 1;
            } else {
                local_br += 1;
            }
        }
        col_base += col_len;  /* incremental update: no multiply, no division */
    }

    /* Phase 2: Cross-column swaps.
     * Element ρ[r,c0] pairs with ρ[r,c1] where c1 = c0 ^ incr.
     * In packed storage, these are in DIFFERENT columns, so we need
     * both column bases. When one is in upper triangle, conjugate. */
    /* ... (more complex, but same idea: iterate over column pairs) ... */
}
```

The key advantage: `col_base` is updated incrementally (`col_base += dim - c`), eliminating
all `pack_idx` multiplications. Within-column swaps are trivially at offset `incr` apart.
Cross-column swaps still need some care but benefit from the incremental column base.

Note: this restructuring may interact with OpenMP parallelization differently. The triangular
iteration in `TRAVERSE_PACKED_BLOCKS` is designed for OpenMP compatibility (the outer loop
parallelizes over `bc`). A column-stride approach would need a different parallelization
strategy.

---

## Finding 3: Branchless Triangle Access

**Estimated impact:** 10-20% at 2-5 qubits.
**Complexity:** Medium.

### The problem

Every off-diagonal block in `TRAVERSE_PACKED_BLOCKS` evaluates multiple branches:

```c
int lower_01 = (r0 >= c1);                           // branch 1
gate_idx_t idx01 = lower_01 ? (col1_base + (r0 - c1))
                            : pack_idx(dim, c1, r0);  // branch 2 (+ pack_idx if upper)

// And in the BLOCK_OP:
if (diag_block) { ... }                               // branch 3
cplx_t rho01 = lower_01 ? data[idx01] : conj(data[idx01]);  // branch 4
data[idx01] = lower_01 ? rho10 : conj(rho10);               // branch 5
```

At small n, these branches have poor prediction because the pattern is irregular.

### Aer's approach

No branches at all. `data[inds[0..3]]` are direct flat-array accesses:
```cpp
std::swap(data_[inds[0]], data_[inds[3]]);  // no branch
std::swap(data_[inds[1]], data_[inds[2]]);  // no branch
```

### Proposed fix

**Step A — Eliminate `pack_idx` call in upper-triangle case:**

Precompute a row-base lookup table:
```c
// Before the loop, precompute row bases for the upper-triangle fallback
gate_idx_t row_bases[dim];  // or use VLA for small dim
for (gate_idx_t r = 0; r < dim; ++r)
    row_bases[r] = r * (2 * dim - r + 1) / 2;
```

Then replace:
```c
// Old:
gate_idx_t idx01 = lower_01 ? (col1_base + (r0 - c1))
                            : pack_idx(dim, c1, r0);  // multiply inside pack_idx

// New:
gate_idx_t idx01 = lower_01 ? (col1_base + (r0 - c1))
                            : (row_bases[r0] + (c1 - r0));  // table lookup, no multiply
```

**Step B — Branchless conjugation:**

```c
// Old (branchy):
cplx_t rho01 = lower_01 ? data[idx01] : conj(data[idx01]);

// New (branchless):
cplx_t raw01 = data[idx01];
double sign = 1.0 - 2.0 * (1 - lower_01);  // +1.0 if lower, -1.0 if upper
cplx_t rho01 = CMPLX(creal(raw01), sign * cimag(raw01));
```

This works because `conj(z) = CMPLX(creal(z), -cimag(z))` — conjugation only flips the
sign of the imaginary part.

**Step C — Eliminate `diag_block` branch:**

Merge the diagonal case into the general path:
```c
cplx_t raw01 = data[idx01];
double sign01 = diag_block ? -1.0 : (lower_01 ? 1.0 : -1.0);
cplx_t rho01 = CMPLX(creal(raw01), sign01 * cimag(raw01));
cplx_t rho10 = data[idx10];

// For X gate: swap (0,1) <-> (1,0)
data[idx10] = rho01;
data[idx01] = CMPLX(creal(rho10), sign01 * cimag(rho10));
```

On diagonal blocks, `idx01 == idx10`, so this reads the same location twice, conjugates,
and writes back — which is exactly what the current `data[idx10] = conj(data[idx10])` does.

At larger n the branch predictor has enough history to predict well, so the benefit
diminishes.

---

## Finding 4: Architectural Lesson — Separate Gate Semantics from Indexing

**Impact:** Code maintainability and extensibility, not direct performance.
**Complexity:** Medium.

### The problem

`cx_packed` (`gate_packed_cx.c`) is 153 lines with deeply nested if/else branches that
interleave "which pairs to swap" (gate semantics) with "where are they in packed storage"
(indexing). Every new permutation gate requires re-deriving this logic from scratch.

Compare Aer's CX — 6 lines:
```cpp
void apply_cnot(const uint_t qctrl, const uint_t qtrgt) {
    std::vector<std::pair<uint_t, uint_t>> pairs = {
        {{1, 3}, {4, 12}, {5, 15}, {6, 14}, {7, 13}, {9, 11}}};
    const reg_t qubits = {{qctrl, qtrgt, qctrl + nq, qtrgt + nq}};
    apply_permutation_matrix(qubits, pairs);  // generic engine handles iteration
}
```

### Proposed abstraction for qlib

```c
typedef struct { int from, to; } perm_pair_t;

/**
 * Apply a permutation gate to a packed density matrix.
 * `qubits` lists the n_q qubits involved (1 for X, 2 for CX/SWAP).
 * `pairs` lists swap pairs in the 2^(2*n_q) truth table.
 * The function handles all packed-storage indexing and conjugation.
 */
void permute_packed(state_t *state, const qubit_t *qubits, int n_q,
                    const perm_pair_t *pairs, int n_pairs);
```

Then gates become trivial:
```c
void x_packed(state_t *state, const qubit_t target) {
    // X truth table on 2 qubits {target, target+n}: swap 0↔3, swap 1↔2
    const perm_pair_t pairs[] = {{0, 3}, {1, 2}};
    permute_packed(state, &target, 1, pairs, 2);
}

void cx_packed(state_t *state, const qubit_t control, const qubit_t target) {
    // CX truth table on 4 qubits: same pairs as Aer
    const perm_pair_t pairs[] = {{1,3}, {4,12}, {5,15}, {6,14}, {7,13}, {9,11}};
    const qubit_t qs[] = {control, target};
    permute_packed(state, qs, 2, pairs, 6);
}

void swap_packed(state_t *state, const qubit_t q1, const qubit_t q2) {
    const perm_pair_t pairs[] = {{1,2}, {4,8}, {5,10}, {6,9}, {7,11}, {13,14}};
    const qubit_t qs[] = {q1, q2};
    permute_packed(state, qs, 2, pairs, 6);
}
```

The `permute_packed()` internals are non-trivial (essentially the current `cx_packed`
generalized), but they are written once and shared by all permutation gates. This:

1. Reduces `cx_packed` from 153 lines to ~6 lines
2. Makes it trivial to add new permutation gates (Toffoli, multi-controlled X, etc.)
3. Concentrates optimization effort: optimize `permute_packed()` once, all gates benefit
4. Enables the precomputed-index-table optimization (Finding 1) to be implemented once

---

## What qlib Should NOT Do

**Do not switch to a flat dim² array.** The 2× memory savings from Hermitian packing is
essential for the larger qubit counts (10-20) that are qlib's primary use case. The small-n
overhead is a fixed cost that can be mitigated with the targeted optimizations above without
abandoning the architectural advantage.

**Do not add a separate flat-array backend just for small n.** The code maintenance burden
of a third storage format would outweigh the benefits. The precomputed-index-table approach
(Finding 1) achieves similar performance gains within the existing packed format.

---

## Key Source Files Referenced

- `extern/aer-dm/include/simulators/statevector/indexes.hpp` — Aer's `index0()` and `apply_lambda`
- `extern/aer-dm/include/simulators/density_matrix/densitymatrix.hpp` — Aer's specialized gates
- `extern/aer-dm/include/simulators/statevector/qubitvector.hpp` — Aer's `apply_permutation_matrix`
- `include/gate.h` — qlib's `insertBit0`, `pack_idx` definitions
- `src/gate/packed/gate_packed_1q.c` — qlib's packed 1Q gate with `TRAVERSE_PACKED_BLOCKS`
- `src/gate/packed/gate_packed_cx.c` — qlib's packed CX with triangle branching
- `src/gate/tiled/gate_tiled_x.c` — qlib's tiled X gate
- `src/gate/tiled/gate_tiled_cx.c` — qlib's tiled CX gate (644 lines)
