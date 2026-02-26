# Circuit Module: Technical Specification

**Module:** Parametrised Quantum Circuits and Pauli Operations
**Status:** Partially implemented — Pauli string construction, phase, and pure-state application complete. Exponential (`pauli_exp.c`) and circuit types (`circuit.c`) pending.
**Depends on:** State module, Gate module
**Consumed by:** MEAS module
**Last updated:** 2026-02-26

---

## Position in Module Stack

```
State ──── Gate
              │
              ▼
           Circuit   ◄── (this module)
              │
              ▼
            MEAS
```

Circuit is the bridging layer between Gate and MEAS. It owns the circuit
representation (`pqc_t`, `lcu_circuit_t`) and the Pauli string machinery
that MEAS consumes for expectation value, gradient, and moment matrix
computation. `pqc_t`, `param_channel_t`, `lcu_layer_t`, and `lcu_circuit_t`
formerly resided in MEAS; they move here.

---

## Module Layout

```
include/
  circuit.h                    # Public API — implemented types and declarations only
                               # (param_channel_t, pqc_t, lcu_* types not yet present)

src/circuit/
  pauli.c                      # ✓ DONE: pauli_new, pauli_from_str, pauli_free,
                               #         pauli_phase, pauli_apply_pure
  pauli_exp.c                  # ✗ PENDING: Pauli exponential PURE/PACKED/TILED
  circuit.c                    # ✗ PENDING: pqc_t, lcu_circuit_t constructors,
                               #            channel constructors

test/circuit/
  test_circuit.c               # ✓ DONE: test runner + main() — 30 tests, 0 failures
  test_pauli.c                 # ✓ DONE: Pauli construction, phase, apply_pure tests
  test_pauli_exp.c             # ✗ PENDING: Pauli exponential tests
```

**Build (circuit tests only):**
```bash
cmake --preset debug
cmake --build --preset debug --target test_circuit
./cmake-build-debug/test/test_circuit
```

---

## Types

### Pauli Label

```c
typedef enum {
    PAULI_I = 0,
    PAULI_X = 1,
    PAULI_Y = 2,
    PAULI_Z = 3
} pauli_label_t;
```

### Pauli String

```c
#include <stdint.h>   // uint64_t

typedef struct {
    pauli_label_t *labels;   // dense array, length n_qubits, heap-allocated
    qubit_t        n_qubits; // matches state_t.qubits (unsigned char, max 255)
    // Precomputed bitmasks (set by pauli_new / pauli_from_str):
    uint64_t  x_mask;  // bit j set  iff  labels[j] ∈ {X, Y}  (flip bit j)
    uint64_t  z_mask;  // bit j set  iff  labels[j] ∈ {Y, Z}  (Z-phase factor)
    uint64_t  y_mask;  // bit j set  iff  labels[j] = Y        (extra i factor)
} pauli_t;
```

**Encoding convention:** labels[j] acts on qubit j — the j-th bit position in the
little-endian computational-basis index (matching STATE_MODULE.md: index 5 = binary
101 = q₀=1, q₁=0, q₂=1).

**Constraint:** n_qubits ≤ 64. Systems requiring more than 64 qubits are out of scope.
`uint64_t` (from `<stdint.h>`) is used for masks rather than `unsigned long` because
`unsigned long` is 32 bits on Windows (LLP64 data model); `uint64_t` is exactly 64 bits
on all platforms per C17. `qubit_t` is used for `n_qubits` to match `state_t.qubits`.

**Derived masks** (all precomputed in pauli_new):
- `x_mask`: bits where the Pauli flips the state. Used to compute the paired index x' = x XOR x_mask.
- `z_mask`: bits that contribute a (−1)^{bit} phase factor. Includes both Y (which has a Z component) and Z labels.
- `y_mask`: bits that contribute an extra factor of i (or −i when the bit is 1). Equals the subset of z_mask where the Pauli is Y.

### Parametrised Channel

```c
typedef struct param_channel {
    void   (*apply)     (state_t *state, double theta, const void *ctx);
    double (*grad_inner)(const state_t *mu, const state_t *tmp,
                         double theta, const void *ctx);
    const void *ctx;      // channel-specific data, owned externally
    int    has_gradient;  // 1 = ideal unitary, grad_inner defined;
                          // 0 = noisy channel, grad_inner is NULL
} param_channel_t;
```

**`apply(state, θ, ctx)`**: Applies the channel in the forward direction.
- Unitary channel: state → U(θ) state (or U(θ)ρU(θ)† for density matrices).
- Noisy channel: state → Σ_a K_a state K_a† (full Kraus map; θ may be used or ignored).

**Adjoint:** For all unitary channel types, U(θ)† = U(−θ) (i.e., exp(+iθ/2 P) =
exp(−i(−θ)/2 P)). The backward pass calls `apply(state, −theta, ctx)` directly; no
separate `apply_adjoint` pointer is needed or stored.

**`grad_inner(mu, tmp, θ, ctx)`**: Returns −2 Im ⟨μ|G|tmp⟩ (pure state) or
−2 Im Tr(H_f G σ) (mixed state), where G is the generator of U(θ). Only meaningful
when `has_gradient = 1`.

**ctx lifetime:** The `ctx` pointer is borrowed — the channel does not own it. The
caller must keep `ctx` alive for at least the lifetime of the channel.

### Flat Parametrised Circuit

```c
typedef struct {
    param_channel_t **channels;  // array of pointers, length n_channels
    int               n_channels;
    double           *params;    // one scalar per channel, length n_channels; owned
    int               is_noisy;  // 0: all channels are ideal; 1: at least one is noisy
} pqc_t;
```

A circuit is **entirely ideal or contains noise**. `is_noisy` is set automatically by
`pqc_new` based on whether any channel has `has_gradient = 0`. When `is_noisy = 1`,
MEAS will refuse to compute gradients.

### LCU Layer

```c
typedef struct {
    void (**unitaries)(state_t *state);  // K unparametrised unitary function pointers
    int    K;
} lcu_layer_t;
```

Each unitary is a fixed gate sequence (no free parameter). The function pointers do
not receive a context argument — the caller captures all required data in a closure
or pre-bound struct. The `K` unitaries are indexed 0 to K−1 and selected by the
multi-index used in moment matrix computation.

### LCU Circuit

```c
typedef struct {
    lcu_layer_t **layers;   // L layer pointers
    int           L;
    int           N;        // total multi-indices: N = ΠK_k, product over k = 0..L−1
} lcu_circuit_t;
```

---

## Pauli Strings

### Construction

```c
// Allocate from a dense label array. Precomputes all three masks.
// Returns NULL on OOM or if n_qubits > 64.
pauli_t *pauli_new(const pauli_label_t *labels, qubit_t n_qubits);

// Parse from a string of characters 'I','X','Y','Z' of length n_qubits.
// Returns NULL on invalid characters, wrong length, or OOM.
pauli_t *pauli_from_str(const char *s, qubit_t n_qubits);

// Free labels array and struct. Safe to call with NULL.
void pauli_free(pauli_t *p);
```

**String ordering in `pauli_from_str`:** Big-endian — `s[0]` maps to qubit `n_qubits−1`
(MSB) and `s[n_qubits−1]` maps to qubit 0 (LSB). This matches standard physics notation
where the leftmost symbol acts on the highest-index qubit. Example for n=4:

    pauli_from_str("IXYZ", 4):
      s[0]='I' → labels[3] = I  (qubit 3)
      s[1]='X' → labels[2] = X  (qubit 2)
      s[2]='Y' → labels[1] = Y  (qubit 1)
      s[3]='Z' → labels[0] = Z  (qubit 0)
      x_mask = 0b0110 = 6,  z_mask = 0b0011 = 3,  y_mask = 0b0010 = 2

`pauli_new` uses the opposite (little-endian) convention: `labels[j]` acts on qubit j.
Use `pauli_from_str("ZYX", 3)` to get the same result as `pauli_new({X,Y,Z}, 3)`.

### Phase Function

For Pauli string P applied to computational basis state |x⟩:

    P |x⟩ = ε(x) · |x XOR x_mask⟩

where the phase ε(x) is one of {+1, −1, +i, −i}. Derivation:
- Each X factor at qubit j: flips bit j, no phase.
- Each Y factor at qubit j: flips bit j, contributes factor i · (−1)^{x_j}.
- Each Z factor at qubit j: no flip, contributes factor (−1)^{x_j}.

Combining all factors:

    ε(x) = i^{popcount(y_mask)} · (−1)^{popcount(x AND z_mask)}

Since ε(x) ∈ {+1, −1, +i, −i}, it can be represented as an integer pair:

```c
typedef struct { int real; int imag; } pauli_phase_t;
// Examples: +1 = {1,0}, -1 = {-1,0}, +i = {0,1}, -i = {0,-1}

pauli_phase_t pauli_phase(const pauli_t *P, uint64_t x);
```

**Key identity:** ε(x') = ε(x) · (−1)^{popcount(y\_mask)}.

Proof: popcount((x XOR x\_mask) AND z\_mask) = popcount((x AND z\_mask) XOR y\_mask)
= popcount(x AND z\_mask) + popcount(y\_mask) (mod 2), using x\_mask AND z\_mask = y\_mask.
Therefore ε(x') = i^{pop(y)} · (−1)^{pop(x AND z)+pop(y)} = ε(x) · (−1)^{pop(y\_mask)}.

Corollaries:
- ε(x) · ε(x') = (−1)^{popcount(y\_mask)} (depends only on the Pauli string, not x).
- ε(x) · ε*(x') = ε(x) · ε(x)* · (−1)^{pop(y)} = (−1)^{pop(y\_mask)} ∈ {+1, −1}.
- If popcount(y\_mask) is even: ε(x') = ε(x). If odd: ε(x') = −ε(x).

This last corollary is useful for the pauli\_exp implementation: when pop(y\_mask) is odd,
the phase factor for the paired element is always the negation of the local phase,
so only one phase evaluation per pair is needed.

**Important identity:** ε(x) · ε*(x') = (−1)^{popcount(y\_mask)} = ±1 (always real).
This follows from P being Hermitian: (PρP)_{xx} is real.

### Pure-State Application

```c
// Apply Pauli string P to a pure statevector in-place: |ψ⟩ → P|ψ⟩
// Requires state->type == PURE and (P->x_mask | P->z_mask) < 2^state->qubits.
// Fails fast via GATE_VALIDATE on null pointers, wrong state type, or mask overflow.
void pauli_apply_pure(state_t *state, const pauli_t *P);
```

**Algorithm — three code paths dispatched by mask structure:**

**Path 1 — Identity** (`x_mask = 0, z_mask = 0`): immediate return. All labels are I;
the Pauli string is a global phase on pure states and requires no work.

**Path 2 — Diagonal** (`x_mask = 0, z_mask ≠ 0`): `P|ψ⟩[i] = ε(i) · ψ[i]` for
all i. No amplitude swapping occurs. Two sub-paths:

- **Pure Z-string** (`y_mask = 0`): ε(i) = (−1)^{popcount(i AND z\_mask)} ∈ {±1}.
  Implemented as a branchless real sign flip `data[i] *= 1.0 - 2.0*parity` — no complex
  arithmetic required.
- **Diagonal with Y-factors** (`y_mask ≠ 0`): unreachable in practice, because
  `y_mask ≠ 0` implies `x_mask ≠ 0`. The branch exists for completeness but is dead
  code for any valid `pauli_t`. See **Known Pitfalls** below.

**Path 3 — Anti-diagonal** (`x_mask ≠ 0`): amplitudes come in pairs (j, j'), where
j' = j XOR x_mask. For each pair:

    tmp     = ψ[j]
    ψ[j]   = ε(j') · ψ[j']      // (P|ψ⟩)[j]  = ε(j') · ψ[j']
    ψ[j']  = ε(j)  · tmp        // (P|ψ⟩)[j'] = ε(j)  · ψ[j]

Since ε ∈ {±1, ±i}, each multiplication is zero-flop: implemented via `apply_ipow`
which uses only component swaps and negations (no floating-point multiplications).

**Pair enumeration** (critical implementation detail): there are exactly dim/2 unique
pairs. To enumerate them without double-counting, insert a 0-bit at the **lowest set
bit of x_mask** (call it `split_bit = __builtin_ctzll(x_mask)`). For k = 0..dim/2−1:

    j  = insertBit0(k, split_bit)   // j has 0 at split_bit
    j' = j XOR x_mask               // j' has 1 at split_bit (and all other X-bits flipped)

This uses a single `insertBit0` call regardless of how many bits x_mask has. The
naive approach of inserting zeros at ALL x_mask bit positions would give only
dim/2^{popcount(x\_mask)} pairs — correct only for single-bit x_masks — and silently
skip all remaining pairs for multi-qubit Pauli strings. This was an actual bug caught
by the test suite during implementation.

**Phase computation:** use the key identity ε(j') = ε(j) · (−1)^{popcount(y\_mask)}
to compute ε(j') from ε(j) with one extra arithmetic operation rather than a second
popcount call:

    parity_j  = popcount(j  AND z_mask) & 1
    exp_j     = (eta_exp + 2 * parity_j) & 3      // ε(j)  = i^{exp_j}
    exp_jp    = (exp_j  + 2 * y_parity) & 3       // ε(j') = i^{exp_jp}

where `eta_exp = popcount(y_mask) & 3` and `y_parity = popcount(y_mask) & 1` are
constants for the Pauli string.

**OpenMP:** `#pragma omp parallel for if(dim >= OMP_THRESHOLD)` where
`OMP_THRESHOLD = 4096`. Pairs are independent — no data races.

**Validation:**
```c
GATE_VALIDATE(state && state->data, "...: null state or data pointer");
GATE_VALIDATE(P,                    "...: null Pauli pointer");
GATE_VALIDATE(state->type == PURE,  "...: expected PURE state");
GATE_VALIDATE(state->qubits >= 64 ||
              (P->x_mask | P->z_mask) >> state->qubits == 0,
              "...: Pauli mask exceeds state dimension");
```

The mask check guards against `>> state->qubits` being UB when `state->qubits ≥ 64`
(shift of a 64-bit type by 64+ is undefined in C17).

### Known Pitfalls

**Dead code branch:** The "Diagonal with Y-factors" sub-path (`x_mask = 0, y_mask ≠ 0`)
is structurally unreachable. By construction, `y_mask ⊆ x_mask` (every Y label sets
both its x_mask and z_mask bit), so `y_mask ≠ 0` implies `x_mask ≠ 0`. The branch is
retained in `pauli.c` as defensive code with a comment noting its unreachability.

**Pair enumeration — single insertBit0 only:** Only insert at `split_bit = ctz(x_mask)`,
not at all x_mask bit positions. Chaining insertBit0 at all positions yields
dim/2^{popcount(x\_mask)} pairs, which is wrong for multi-qubit Pauli strings.

### Identity Pauli Strings

When `x_mask = 0` (all labels are I or Z), the Pauli string is diagonal in the
computational basis: P|x⟩ = ε(x)|x⟩ with no bit flip.

- **Pure state:** exp(−iθ/2 P) applies a state-dependent phase exp(−iθ/2 ε(x)) to
  each amplitude.
- **Density matrix:** ρ'_{xy} = exp(−iθ/2 (ε(x) − ε(y))) ρ_{xy}. The diagonal
  ρ_{xx} is unchanged (global phase cancels). Off-diagonal elements acquire a
  relative phase.

The all-identity string I⊗n (x_mask = 0, z_mask = 0) is a global phase on pure
states and a complete no-op on density matrices.

---

## Pauli Exponential

### Mathematical Definition

    U(θ) = exp(−iθ/2 · P)

Since P is Hermitian and P² = I:

    U(θ) = cos(θ/2) I − i sin(θ/2) P

The adjoint is U(θ)† = exp(+iθ/2 P) = U(−θ).

### Algorithm: PURE State

```
pauli_exp_pure(state, P, theta):
    c = cos(theta/2),  s = sin(theta/2)

    for each x in 0 .. 2^n - 1  with  x < (x XOR x_mask):
        x' = x XOR x_mask
        eps_x  = ε(x)   as (re, im)
        eps_x' = ε(x')  as (re, im)

        a = state->data[x]
        b = state->data[x']

        state->data[x]  = c*a - i*s*eps_x' * b   // ε(x') is the phase for row x
        state->data[x'] = c*b - i*s*eps_x  * a   // ε(x)  is the phase for row x'
```

**Derivation:** P_{row, col} = ε(col) · δ_{row, col XOR x\_mask}.
So (P|ψ⟩)_x = ε(x') · ψ(x')  (the only nonzero column contributing to row x is x').
Therefore ψ'(x) = c·ψ(x) − i·s·ε(x')·ψ(x'), and symmetrically for x'.

The inner update is a 2×2 rotation mixing `a` and `b`. Since ε(x') ∈ {±1, ±i}, each
multiplication by `−i·ε(x')` reduces to a combination of negation, swapping real/imag,
and sign changes — no general complex multiply is needed. When pop(y\_mask) is odd,
eps\_x' = −eps\_x (by the key identity), halving the number of phase evaluations.

**OpenMP:** parallel over pairs, `if (dim >= OMP_THRESHOLD)`.

**Edge case x_mask = 0:** For non-all-identity diagonals (z_mask ≠ 0), iterate over
all x without pairing, applying state->data[x] *= exp(−iθ/2 · ε(x)).

### Algorithm: MIXED_PACKED and MIXED_TILED

For density matrices ρ → U(θ) ρ U(θ)†:

    ρ'_{xy} = c² ρ_{xy} + s² ε(x) ε*(y) ρ_{x',y'}
             − ics ε(x') ρ_{x',y} + ics ε*(y') ρ_{x,y'}

where c = cos(θ/2), s = sin(θ/2), x' = x XOR x_mask, y' = y XOR x_mask.

**Note on the s² coefficient:** ε(x')ε*(y') = ε(x)ε*(y). This follows from the key
identity ε(x') = ε(x)·(−1)^{pop(y\_mask)}: the (−1) factors cancel in the product,
so `ε(x)ε*(y)` and `ε(x')ε*(y')` are interchangeable. The cross terms do NOT share
this cancellation and must use ε(x') and ε*(y') explicitly.

**Four-element coupling:** Each output element ρ'_{xy} is a linear combination of
four input elements: ρ_{xy}, ρ_{x',y'}, ρ_{x',y}, ρ_{x,y'}. These four elements
form one independent mixing group. Elements in a group are processed atomically —
all four inputs must be read before any output is written.

**Group enumeration:** Groups are indexed by the representative pair (x, y) with:
- x < x' = x XOR x\_mask  (pick the lower row of the row-pair), **and**
- y < y' = y XOR x\_mask  (pick the lower column of the column-pair).

The second condition (y < y') is equivalent to: the lowest set bit of x\_mask is 0
in y. Without it, both (x, y) and (x, y') would qualify as representatives for the
same 4-element group, causing each group to be processed twice.

For each representative (x, y), the triangle membership of all four elements
(upper/lower/diagonal) is determined once before processing. Upper-triangle elements
are accessed via Hermitian conjugation of the transposed index.

**Special case y = x (diagonal representative):** When y = x, the group contains
both diagonal elements ρ_{xx} and ρ_{x'x'} along with the off-diagonal pair ρ_{xx'}
and ρ_{x'x}. Since y' = x' > x = y, the condition y < y' is satisfied, so diagonal
representatives are included in the main traversal with no special case needed.

**MIXED_PACKED traversal:**

Iterate over representative pairs (x, y) with x < x' and y < y':

1. Determine: `x' = x XOR x_mask`, `y' = y XOR x_mask`.
2. Load all four elements (accessing upper-triangle via conjugate if needed).
3. Compute the four new values using the formula above.
4. Write back, enforcing Hermitian symmetry for upper-triangle positions.

Diagonal block guard: when x = y (diagonal representative), ρ'_{xx} and ρ'_{x'x'}
are guaranteed real by unitarity (UρU† is Hermitian; its diagonal is real). In
floating-point, rounding may introduce a tiny imaginary part — zero it explicitly.
This is a numerical hygiene step, not a mathematical correction.

**MIXED_TILED traversal:** Same 4-element mixing kernel. The traversal structure
mirrors `gate_tiled_h.c` (Hadamard on tiled storage): iterate over tile pairs
(tr, tc) with tr ≥ tc, then over the within-tile element groups using the
dtile_read/dtile_write helpers from `gate_tiled.h`. Diagonal tiles require the same
conjugation guard as in single-qubit tiled gate implementations.

**OpenMP:** parallel over tile rows (MIXED_TILED) or element rows (MIXED_PACKED),
`if (dim >= OMP_THRESHOLD)`. `OMP_THRESHOLD = 512`.

### Gradient Inner Product

For U(θ) = exp(−iθ/2 P), the generator is G = P/2:

    grad_inner(mu, tmp, θ) = −2 Im ⟨μ|G|tmp⟩ = −Im ⟨μ|P|tmp⟩

#### Pure State — O(2^n)

```
S = Σ_{x=0}^{2^n−1}  conj(mu->data[x]) · ε(x') · tmp->data[x']
    where x' = x XOR x_mask
return −Im S
```

**Derivation:** ⟨μ|P|tmp⟩ = Σ_x μ\*(x) (P|tmp⟩)_x = Σ_x μ\*(x) · ε(x') · tmp(x').
The phase factor is ε(x'), not ε(x): (P|tmp⟩)_x = Σ_z P_{xz} tmp(z) = ε(x') · tmp(x')
(only z = x' contributes, since P_{xz} = ε(z) · δ_{x, z XOR x\_mask} is nonzero
only when z = x' = x XOR x\_mask).

One loop over all amplitudes. Each iteration: one phase lookup (at x'), one multiply-accumulate.

#### Mixed State — O(4^n), Direct Storage Traversal

**Derivation.** Expanding `Tr(H_f · P · σ)` using the matrix elements of P
(P_{ab} = ε(b) δ_{a, b'} where b' = b XOR x_mask):

    Tr(H_f P σ) = Σ_{a,b} (H_f)_{b, a'} · ε(a) · σ_{ab}

where `a' = a XOR x_mask`. This is a bilinear form: each element σ_{ab} is paired
with the element of H_f at coordinates `(b, a')`.

**Reduction to the stored lower triangle.** Since only lower-triangle elements
(row ≥ col) are stored, and both H_f and σ are Hermitian, expand by cases:

```
S = 0

// Diagonal elements of σ  (one contribution each)
for each a = 0..d−1:
    S += lookup_hf(a, a') · ε(a) · σ_{aa}

// Off-diagonal stored elements of σ  (a > b, two contributions each)
for each (a, b) with a > b:
    S += lookup_hf(b, a') · ε(a) · σ_{ab}        // (σ in lower tri)
    S += lookup_hf(a, b') · ε(b) · conj(σ_{ab})  // (σ upper tri by symmetry)

return −Im S
```

where:
```
lookup_hf(r, c):
    if r >= c:  return mu->data[packed_index(r, c)]          // lower tri, direct
    else:       return conj(mu->data[packed_index(c, r)])    // upper tri, conjugate
```

**Correctness check — identity string.** For P = I^⊗n (x_mask = 0, z_mask = 0):
ε(a) = 1 for all a and a' = a. Then S = Σ_{a,b} (H_f)_{ba} σ_{ab} = Tr(H_fᵀ σ).
Since H_f is Hermitian, H_fᵀ = H_f* and Tr(H_f* σ) = conj(Tr(H_f σ)).
Because H_f and σ are both Hermitian, Tr(H_f σ) is real, so −Im S = 0. ✓
(The all-identity Pauli is a no-op on density matrices, so its gradient is zero.)

**MIXED_PACKED:** Direct pointer arithmetic on `mu->data` and `tmp->data`.
`packed_index(r, c)` = `c * d - c*(c+1)/2 + r` (matching STATE_MODULE.md).
Iterate over `a = 0..d−1` for diagonal, then over off-diagonal pairs using the
same outer loop structure as the packed 2Q gate traversal.
**OpenMP:** `reduction(+:S_real, S_imag)` over the outer loop; `if (dim >= OMP_THRESHOLD)`.

**MIXED_TILED:** Replace `packed_index` with `dtile_read(mu->data, r, c, n_tiles)`
from `gate_tiled.h`, which handles tile offset arithmetic and conjugation for
upper-triangle tile blocks. The outer loop iterates over tile pairs (tr, tc) with
tr ≥ tc, then over local indices within each tile — the same structure as the
tiled Pauli exponential traversal.
**OpenMP:** parallel over tile rows; `if (dim >= OMP_THRESHOLD)`.

---

## Built-in Channel Constructors

```c
// Pauli exponential channel. U(θ) = exp(−iθ/2 P).
// Generator G = P/2. has_gradient = 1.
// P is borrowed — caller must keep P alive for the lifetime of the channel.
// Returns NULL on OOM.
param_channel_t *pauli_exp_channel_new(const pauli_t *P);

// Noisy (Kraus) channel. The caller encapsulates the entire Kraus map
// ρ → Σ_a K_a ρ K_a† in one function.
// apply_adjoint = NULL, grad_inner = NULL, has_gradient = 0.
// ctx is borrowed. Returns NULL on OOM.
param_channel_t *kraus_channel_new(
    void       (*apply)(state_t *state, double theta, const void *ctx),
    const void  *ctx
);

// Generic channel with user-supplied function pointers.
// Caller sets has_gradient and may pass NULL for grad_inner if inapplicable.
// Adjoint is always obtained by calling apply with −theta; no separate pointer.
// ctx is borrowed. Returns NULL on OOM.
param_channel_t *channel_new(
    void   (*apply)     (state_t *, double, const void *),
    double (*grad_inner)(const state_t *, const state_t *, double, const void *),
    const void *ctx,
    int    has_gradient
);

// Free the param_channel_t struct. Does NOT free ctx.
void channel_free(param_channel_t *ch);
```

**Pauli exponential channel internals:**

The ctx pointer holds a `const pauli_t *`. The three function pointers are set to
internal static functions:
- `apply`:      dispatches on `state->type` to `pauli_exp_pure` / `pauli_exp_packed` /
               `pauli_exp_tiled`; the adjoint U(θ)† is obtained by the caller passing −theta
- `grad_inner`: computes −Im ⟨μ|P|tmp⟩ as described above

---

## Circuit API

### Flat Parametrised Circuit

```c
// Create a circuit.
// - Copies the channels pointer array (not the channel objects themselves).
// - Allocates params[n_channels], initialised to 0.0.
// - Sets is_noisy = 1 if any channel has has_gradient == 0, else 0.
// Returns NULL on OOM.
pqc_t *pqc_new(param_channel_t **channels, int n_channels);

// Free the circuit struct and params array.
// Does NOT free the channel objects or their ctx pointers.
void pqc_free(pqc_t *circuit);

// Parameter access
void   pqc_set_params(pqc_t *circuit, const double *params); // copies n_channels values
double pqc_get_param(const pqc_t *circuit, int c);
void   pqc_set_param(pqc_t *circuit, int c, double val);
```

### LCU Structures

```c
// LCU layer holding K unparametrised unitary function pointers.
// Copies the fn_arr array into an internally allocated array (layer->unitaries).
// Returns NULL on OOM.
lcu_layer_t *lcu_layer_new(void (**fn_arr)(state_t *), int K);

// Free the lcu_layer_t struct and its internal function pointer array.
void lcu_layer_free(lcu_layer_t *layer);

// LCU circuit holding L layers.
// Copies the layers pointer array (does not copy the lcu_layer_t objects).
// Computes N = ΠK_k and stores it.
// Returns NULL on OOM.
lcu_circuit_t *lcu_circuit_new(lcu_layer_t **layers, int L);

// Free the circuit struct and the internal layers pointer array.
// Does NOT free the lcu_layer_t objects.
void lcu_circuit_free(lcu_circuit_t *circuit);

// Decode a linear multi-index I (0 ≤ I < circuit->N) into per-layer indices.
// Writes layer-k index into idx_out[k] for k = 0..L−1.
// Multi-index ordering: idx_out[0] is the least significant (changes fastest).
void lcu_multi_index(const lcu_circuit_t *circuit, int I, int *idx_out);
```

**Multi-index encoding:** idx_out[0] cycles fastest, idx_out[L−1] slowest:

    I = idx_out[0] + K_0 * (idx_out[1] + K_1 * (idx_out[2] + ...))

Applying multi-index I to an initial state means sequentially applying:
`layers[0].unitaries[idx_out[0]](state)`, then `layers[1].unitaries[idx_out[1]](state)`,
etc. (layer 0 first, layer L−1 last).

---

## Interface with MEAS

The MEAS module is the exclusive consumer of `pqc_t`, `lcu_circuit_t`, and all
circuit-level types. The Circuit module exposes these types and their constructors;
all computation algorithms (expectation value, gradient, moment matrix) live in MEAS.

**MEAS uses from Circuit:**
- `pqc_t` — for `expectation_value` and `gradient_adjoint`
- `lcu_circuit_t` — for `moment_matrix`
- `param_channel_t` function pointers — the 3-pass adjoint algorithm calls `apply`,
  `apply_adjoint`, and `grad_inner` through these pointers
- `lcu_multi_index` — to decode column and row indices during moment matrix computation

**MEAS defines (outside Circuit's scope):**
- `diag_hamiltonian_t` — diagonal cost Hamiltonian and its operations
- QAOA channel constructors — cost channel (diagonal phase multiply using
  `diag_hamiltonian_t`) and mixer channel (n sequential Rx calls); both use
  `channel_new` from Circuit
- `qaoa_circuit_new` — assembles a depth-p QAOA `pqc_t`
- `expectation_value`, `gradient_adjoint`, `moment_matrix`, `moment_matrix_load`

**Consequence for MEAS_MODULE.md:** All type definitions for `param_channel_t`,
`pqc_t`, `lcu_layer_t`, and `lcu_circuit_t` must be removed from MEAS and replaced
with `#include "circuit.h"`. The following specific changes are required in MEAS:

- **`param_channel_t`**: Remove the `apply_adjoint` field (it was never there in MEAS;
  this module does not add it). Pass 2 of the backward algorithm calls
  `apply(state, −theta, ctx)` directly.
- **`pqc_t.params`**: Now struct-owned — `pqc_new` allocates, `pqc_free` frees. Remove
  any caller-side `free(params)` in MEAS. The MEAS spec comment "caller-owned" was
  incorrect.
- **`pqc_free`**: Defined here; remove the `pqc_free` declaration from the MEAS API
  section to avoid duplicate declarations.
- **`lcu_layer_t.unitaries`**: Field is now named `unitaries` (renamed from `apply`
  to avoid collision with `param_channel_t.apply`). MEAS already used `unitaries`.
- **`lcu_circuit_t.layers`**: This module uses `lcu_layer_t **layers` (array of
  pointers; layers are externally owned and borrowed). MEAS used `lcu_layer_t *layers`
  (flat array). MEAS must be updated to `lcu_layer_t **layers`.
- **`unitary_fn_t`**: MEAS Preliminary 2 referenced a shared `unitary_fn_t` typedef.
  This module uses the raw function pointer type `void (**)(state_t *)` directly in the
  struct; `unitary_fn_t` is not defined here. MEAS can define it locally if needed for
  readability, or use the raw type.

---

## Test Coverage

**Last verified:** 2026-02-26 — 30 tests, 0 failures, tolerance ε = 1e−12.

**Reference strategy:** Every `check_pauli_apply` call constructs the full 2^n × 2^n
Pauli matrix (`build_pauli_matrix`), multiplies it against the original statevector
(CBLAS `zgemv` via `zmv`), then compares element-wise against `pauli_apply_pure`
output. Linalg helpers from `test/utility/linalg.c` are linked into `test_circuit`.

### test_pauli.c — Implemented (30 tests)

#### Pauli Construction and Masks (5 tests)

| Test | Description |
|------|-------------|
| `test_pauli_masks_X` | n=1 X: verify x_mask=1, z_mask=0, y_mask=0 |
| `test_pauli_masks_Y` | n=1 Y: verify x_mask=1, z_mask=1, y_mask=1 |
| `test_pauli_masks_Z` | n=1 Z: verify x_mask=0, z_mask=1, y_mask=0 |
| `test_pauli_masks_str` | `pauli_from_str("IXYZ",4)` big-endian: x_mask=6, z_mask=3, y_mask=2; "IIII" (all zero); "XXXX" (x_mask=15) |
| `test_pauli_from_str_negative` | NULL string, invalid char 'A', embedded invalid char, length-too-short, length-too-long — all return NULL |

#### pauli_phase (1 test)

| Test | Description |
|------|-------------|
| `test_pauli_phase` | Spot-checks integer phase components for X(x=0,1), Y(x=0,1), Z(x=0,1) using `TEST_ASSERT_EQUAL_INT` |

#### pauli_apply_pure — Single-qubit (4 tests)

All four single-qubit tests run against all 6 standard states:
{|0⟩, |1⟩, |+⟩, |−⟩, |+y⟩, |−y⟩} (Z, X, and Y eigenbases).

| Test | Pauli | x_mask | z_mask | split_bit | Code path |
|------|-------|--------|--------|-----------|-----------|
| `test_pauli_apply_I_q0` | I | 0 | 0 | — | Identity (early return) |
| `test_pauli_apply_X_q0` | X | 1 | 0 | 0 | Anti-diagonal |
| `test_pauli_apply_Y_q0` | Y | 1 | 1 | 0 | Anti-diagonal + Y-phase |
| `test_pauli_apply_Z_q0` | Z | 0 | 1 | — | Diagonal pure-Z |

#### pauli_apply_pure — Two-qubit strings (7 tests)

All n=2 tests run against all 7 defined states:
{|00⟩, |01⟩, |10⟩, |11⟩, |++⟩, |+0⟩, |1+⟩}.

| Test | Pauli | x_mask | z_mask | split_bit | Notes |
|------|-------|--------|--------|-----------|-------|
| `test_pauli_apply_XX` | XX | 3 | 0 | 0 | Multi-bit x_mask, no phase |
| `test_pauli_apply_ZZ` | ZZ | 0 | 3 | — | Diagonal pure-Z, full mask |
| `test_pauli_apply_XZ` | XZ | 1 | 2 | 0 | Mixed X+Z types |
| `test_pauli_apply_YY` | YY | 3 | 3 | 0 | Multi-bit x_mask, η=2 → −1 phase |
| `test_pauli_apply_IZ` | IZ | 0 | 1 | — | Diagonal, non-trivial z position |
| `test_pauli_apply_IX` | IX | 2 | 0 | **1** | **split_bit=1**: first non-trivial insertBit0 |
| `test_pauli_apply_IY` | IY | 2 | 2 | **1** | **split_bit=1** with Y-phase |

#### pauli_apply_pure — Three-qubit strings (6 tests)

All n=3 tests run against all 11 defined states:
8 computational basis states {|000⟩…|111⟩} plus GHZ = (|000⟩+|111⟩)/√2,
W = (|001⟩+|010⟩+|100⟩)/√3, and |+++⟩.

| Test | Pauli | x_mask | z_mask | split_bit | Notes |
|------|-------|--------|--------|-----------|-------|
| `test_pauli_apply_XYZ` | XYZ | 3 | 6 | 0 | Mixed types; inline spot-check XYZ\|000⟩=i\|011⟩ |
| `test_pauli_apply_ZZZ` | ZZZ | 0 | 7 | — | Full diagonal, dual inline+reference check |
| `test_pauli_apply_IXI` | IXI | 2 | 0 | **1** | **split_bit=1**, 3-qubit |
| `test_pauli_apply_XII` | XII | 4 | 0 | **2** | **split_bit=2** |
| `test_pauli_apply_IZI` | IZI | 0 | 2 | — | Non-uniform diagonal, z_mask=2 |
| `test_pauli_apply_ZIZ` | ZIZ | 0 | 5 | — | Alternating diagonal, z_mask=5 |

#### pauli_apply_pure — Four-qubit strings (4 tests, MAXQUBITS=4)

All n=4 tests run against all 16 computational basis states plus PSI4_GHZ =
(|0000⟩+|1111⟩)/√2.

| Test | Pauli (labels[0..3]) | x_mask | z_mask | split_bit | Notes |
|------|----------------------|--------|--------|-----------|-------|
| `test_pauli_apply_4q` | XZYX | 13 | 6 | 0 | Mixed types, η=1 (i-phase) |
| `test_pauli_apply_4q_split1` | IIXI | 2 | 0 | **1** | **split_bit=1**, 4-qubit |
| `test_pauli_apply_4q_split2` | IXII | 4 | 0 | **2** | **split_bit=2** |
| `test_pauli_apply_4q_split3` | XIII | 8 | 0 | **3** | **split_bit=3**: all {0,1,2,3} covered |

#### Round-trip, involution, and edge cases (3 tests)

| Test | Description |
|------|-------------|
| `test_pauli_from_str_roundtrip` | `pauli_from_str("ZYX",3)` vs `pauli_new({X,Y,Z},3)`: identical masks and application results on all 8 basis states |
| `test_pauli_apply_involution` | P²=I for X, Y, Z (n=1), XX (n=2), XYZ (n=3) applied to \|0…0⟩ |
| `test_pauli_apply_identity_noop` | I⊗³ on \|+++⟩ is an exact no-op |

**split_bit coverage summary:** All four `insertBit0` positions {0, 1, 2, 3} are exercised
at n=4 (MAXQUBITS). Position 0 is the degenerate case (`mask=0`, plain left shift);
positions 1, 2, 3 exercise the nontrivial `(val & mask) | ((val & ~mask) << 1)` path.

### test_pauli_exp.c — Pending

Tests for `pauli_exp_pure`, `pauli_exp_packed`, `pauli_exp_tiled`, and `grad_inner`
will be added when `pauli_exp.c` is implemented. See the **Pauli Exponential** section
above for the planned test design.

---

## Build Integration

The circuit library is a STATIC library in `src/CMakeLists.txt`. It links against
the shared `q` library (which bundles `state` and `gate`) rather than against separate
`state` and `gate` targets, because those are not exposed as standalone CMake targets.

```cmake
# In src/CMakeLists.txt (current implementation — pauli.c only):
add_library(circuit STATIC
    "${CMAKE_CURRENT_SOURCE_DIR}/circuit/pauli.c"
)
target_include_directories(circuit
    PUBLIC  "${CMAKE_SOURCE_DIR}/include"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled"  # for gate_tiled.h (future use)
)
target_compile_definitions(circuit PRIVATE LOG_TILE_DIM=${LOG_TILE_DIM})
if(DEFINED OMP_THRESHOLD)
    target_compile_definitions(circuit PRIVATE OMP_THRESHOLD=${OMP_THRESHOLD})
endif()
target_link_libraries(circuit PUBLIC q)   # q exposes state.h and gate.h transitively
if(ENABLE_OPENMP)
    target_link_libraries(circuit PRIVATE omp_compiler_flags)
endif()

# In test/CMakeLists.txt:
set(TEST_CIRCUIT_SRC
    "${CIRCUIT_DIR}/test_circuit.c"
    "${CIRCUIT_DIR}/test_pauli.c"
    "${UTILITY_DIR}/linalg.c"      # zmv() for reference matrix-vector products
    "${UNITY_SRC}/unity.c"
)
add_executable(test_circuit ${TEST_CIRCUIT_SRC})
target_compile_definitions(test_circuit PRIVATE
    UNITY_INCLUDE_DOUBLE UNITY_SUPPORT_64
    LOG_TILE_DIM=${LOG_TILE_DIM})
if(DEFINED OMP_THRESHOLD)
    target_compile_definitions(test_circuit PRIVATE OMP_THRESHOLD=${OMP_THRESHOLD})
endif()
target_include_directories(test_circuit PRIVATE
    ${UNITY_SRC}
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_link_libraries(test_circuit PRIVATE circuit q)
add_test(NAME test_circuit COMMAND test_circuit)
```

**When `pauli_exp.c` and `circuit.c` are added**, extend the `circuit` STATIC library:
```cmake
add_library(circuit STATIC
    "${CMAKE_CURRENT_SOURCE_DIR}/circuit/pauli.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/circuit/pauli_exp.c"   # add when ready
    "${CMAKE_CURRENT_SOURCE_DIR}/circuit/circuit.c"      # add when ready
)
```
Add `test_pauli_exp.c` to `TEST_CIRCUIT_SRC` at the same time. The `gate/tiled`
PRIVATE include is already wired for `gate_tiled.h` access needed by `pauli_exp.c`.

**Dependency rationale:**
- `q` is a PUBLIC dependency of `circuit` because `circuit.h` exposes `state_t *` in
  function signatures — consumers of `circuit.h` need `state.h` (and `gate.h`)
  transitively.
- `omp_compiler_flags` is PRIVATE because OpenMP is an implementation detail of
  the parallel loops in `pauli.c`; it does not appear in `circuit.h`.
- `linalg.c` is a test-only dependency (provides `zmv` for the reference oracle).
  It is not linked into the `circuit` library itself.

---

## Deferred — Out of Scope

- **General Pauli-string Hamiltonian evolution** (ΣJ_k P_k): no decision on
  Trotterisation vs. direct exponentiation; deferred from MEAS.
- **Pauli strings on n > 64 qubits.**
- **Multi-parameter channels** — each channel owns exactly one scalar parameter.
  Circuits with genuinely independent parameters must use one channel per parameter.
- **Shared-parameter channels** — if two gates share the same θ, the caller bundles
  them into one function pointer.
- **Gradient through Kraus channels** — caller uses parameter-shift or finite
  differences via repeated `expectation_value` in MEAS.
- **Parametrised LCU unitaries** — unitaries in `lcu_layer_t` have no free parameter.
- **Circuit optimisation / compilation passes.**
- **QASM / file-based circuit parsing.**
- **Hessians, natural gradients, quantum geometric tensor.**
