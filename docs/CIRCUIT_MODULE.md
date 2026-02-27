# Circuit Module: Technical Specification

**Module:** Parametrised Quantum Circuits and Pauli Operations
**Status:** Partially implemented — Pauli string construction, phase, and pure-state application complete. Exponential (`pauli_exp.c`) and circuit types (`circuit.c`) pending.
**Depends on:** State module (`state.h`, `state_t`, `qubit_t`), Gate module (`gate.h`, `GATE_VALIDATE`)
**Consumed by:** MEAS module
**Last updated:** 2026-02-27 (audited against source)

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
computation.

---

## Directory Layout

```
include/
  circuit.h                    # Public API — pauli_t, pauli_label_t, pauli_phase_t,
                               #   pauli_new/from_str/free/phase/apply_pure declared.
                               #   param_channel_t, pqc_t, lcu_* types NOT YET present.

src/circuit/
  pauli.c                      # DONE: pauli_new, pauli_from_str, pauli_free,
                               #       pauli_phase, pauli_apply_pure
  pauli_exp.c                  # PENDING: pauli_exp_{pure,packed,tiled}, grad_inner
  circuit.c                    # PENDING: pqc_t, lcu_circuit_t constructors,
                               #          channel constructors

test/circuit/
  test_circuit.c               # DONE: test runner + main() — 30 tests, 0 failures
  test_pauli.c                 # DONE: Pauli construction, phase, apply_pure tests
  test_pauli_exp.c             # PENDING: Pauli exponential and grad_inner tests
```

---

## Current Status and TODO

| Item | Status |
|------|--------|
| `pauli_t` type, masks, construction | Complete |
| `pauli_phase` | Complete |
| `pauli_apply_pure` (PURE state) | Complete |
| `pauli_exp_pure` | Pending — algorithm fully specified |
| `pauli_exp_packed` (MIXED_PACKED) | Pending — algorithm fully specified |
| `pauli_exp_tiled` (MIXED_TILED) | Pending — algorithm fully specified |
| `grad_inner` for all state types | Pending — algorithm fully specified |
| `param_channel_t` type | Pending — type specified, not yet in `circuit.h` |
| `pqc_t` type and constructors | Pending — type specified, not yet in `circuit.h` |
| `lcu_layer_t`, `lcu_circuit_t` | Pending — types specified, not yet in `circuit.h` |
| Channel constructors (`pauli_exp_channel_new`, etc.) | Pending |
| `test_pauli_exp.c` | Pending — blocked on `pauli_exp.c` |
| MEAS migration (remove duplicate type defs) | Pending |

**Build (circuit tests only):**
```bash
cmake --preset debug
cmake --build --preset debug --target test_circuit
./cmake-build-debug/test/test_circuit
```

---

## Inputs — What This Module Consumes

| Source | Items used |
|--------|-----------|
| `state.h` (State module) | `state_t`, `qubit_t`, `state_t.type` (`PURE`/`MIXED_PACKED`/`MIXED_TILED`), `state_t.data`, `state_t.qubits` |
| `gate.h` (Gate module) | `GATE_VALIDATE` macro (fast-fail validation in `pauli_apply_pure`); `insertBit0` (pair enumeration in `pauli_apply_pure`); `gate_idx_t` (loop index type); `gate_tiled.h` (private, for `dtile_read`/`dtile_write` in pending `pauli_exp.c`) |
| `<stdint.h>` | `uint64_t` (bitmasks — exactly 64 bits on all platforms, unlike `unsigned long`) |
| `<math.h>` | `cos`, `sin` (used in pending `pauli_exp.c`) |

---

## Outputs — What This Module Produces

| Artifact | Description |
|----------|-------------|
| `libcircuit.a` | Static library linked by `test_circuit` and (future) the MEAS module |
| `include/circuit.h` | Public header; currently declares Pauli API only; will grow to include channel and circuit types |
| CMake target `circuit` | PUBLIC dependency on `q` (exposes `state.h`/`gate.h` transitively to consumers) |

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

**Note:** `circuit.h` (lines 38–40) has an incorrect comment claiming `Y = 0b11` and `Z = 0b10`; the actual declared values above (`Y=2=0b10`, `Z=3=0b11`) are correct. The mask-construction logic in `pauli_new` uses `switch(labels[j])` and is not affected by the comment error.

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

**Encoding convention:** `labels[j]` acts on qubit j — the j-th bit position in the
little-endian computational-basis index (matching STATE_MODULE.md: index 5 = binary
101 = q₀=1, q₁=0, q₂=1).

**Constraint:** n_qubits ≤ 64. Systems requiring more than 64 qubits are out of scope.

**Derived masks** (all precomputed in `pauli_new`):
- `x_mask`: bits where the Pauli flips the state. Used to compute the paired index x' = x XOR x_mask.
- `z_mask`: bits that contribute a (−1)^{bit} phase factor. Includes both Y and Z labels.
- `y_mask`: bits that contribute an extra factor of i (or −i when the bit is 1). Equals the subset of z_mask where the Pauli is Y. By construction, `y_mask ⊆ x_mask`, so `y_mask ≠ 0` implies `x_mask ≠ 0`.

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

**Adjoint:** For all unitary channel types, U(θ)† = U(−θ). The backward pass calls
`apply(state, −theta, ctx)` directly; no separate adjoint pointer is stored.

**`grad_inner(mu, tmp, θ, ctx)`**: Returns −2 Im ⟨μ|G|tmp⟩ (pure state) or
−2 Im Tr(H_f G σ) (mixed state), where G is the generator of U(θ). Only meaningful
when `has_gradient = 1`.

**ctx lifetime:** Borrowed — the channel does not own it. The caller must keep `ctx`
alive for at least the lifetime of the channel.

### Flat Parametrised Circuit

```c
typedef struct {
    param_channel_t **channels;  // array of pointers, length n_channels
    int               n_channels;
    double           *params;    // one scalar per channel, length n_channels; owned
    int               is_noisy;  // 0: all channels are ideal; 1: at least one is noisy
} pqc_t;
```

`is_noisy` is set automatically by `pqc_new` based on whether any channel has
`has_gradient = 0`. When `is_noisy = 1`, MEAS refuses to compute gradients.

### LCU Layer

```c
typedef struct {
    void (**unitaries)(state_t *state);  // K unparametrised unitary function pointers
    int    K;
} lcu_layer_t;
```

Each unitary is a fixed gate sequence (no free parameter). The `K` unitaries are
indexed 0 to K−1 and selected by the multi-index used in moment matrix computation.

### LCU Circuit

```c
typedef struct {
    lcu_layer_t **layers;   // L layer pointers (borrowed; externally owned)
    int           L;
    int           N;        // total multi-indices: N = ΠK_k, product over k = 0..L−1
} lcu_circuit_t;
```

---

## Pauli Strings

### Construction

```c
// Allocate from a dense label array. Precomputes all three masks.
// Returns NULL on OOM, if n_qubits > 64, or if n_qubits > 0 && labels == NULL.
// n_qubits == 0 is valid: returns a Pauli with NULL labels and all masks zero.
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

**Key identity:** ε(x') = ε(x) · (−1)^{popcount(y\_mask)}, where x' = x XOR x\_mask.

Proof: popcount((x XOR x\_mask) AND z\_mask) = popcount((x AND z\_mask) XOR y\_mask)
= popcount(x AND z\_mask) + popcount(y\_mask) (mod 2), using x\_mask AND z\_mask = y\_mask.

Corollaries:
- ε(x) · ε(x') = (−1)^{popcount(y\_mask)} (depends only on the Pauli string, not x).
- ε(x) · ε*(x') = (−1)^{popcount(y\_mask)} ∈ {+1, −1} (always real — P is Hermitian).
- If popcount(y\_mask) is even: ε(x') = ε(x). If odd: ε(x') = −ε(x).

The last corollary is used by `pauli_exp`: when pop(y\_mask) is odd, only one phase
evaluation per pair is needed.

### Pure-State Application

```c
// Apply Pauli string P to a pure statevector in-place: |ψ⟩ → P|ψ⟩
// Requires state->type == PURE and (P->x_mask | P->z_mask) < 2^state->qubits.
// Fails fast via GATE_VALIDATE on null pointers, wrong state type, or mask overflow.
void pauli_apply_pure(state_t *state, const pauli_t *P);
```

**Algorithm — three code paths dispatched by mask structure:**

**Path 1 — Identity** (`x_mask = 0, z_mask = 0`): immediate return.

**Path 2 — Diagonal** (`x_mask = 0, z_mask ≠ 0`): `P|ψ⟩[i] = ε(i) · ψ[i]` for
all i. No amplitude swapping.

- **Pure Z-string** (`y_mask = 0`): ε(i) = (−1)^{popcount(i AND z\_mask)} ∈ {±1}.
  Implemented as a branchless real sign flip `data[i] *= 1.0 - 2.0*parity`.
- **Diagonal with Y factors** (`y_mask ≠ 0, x_mask = 0`): Structurally unreachable —
  `PAULI_Y` always sets both `x_mask` and `y_mask` bits, so `y_mask ≠ 0` implies
  `x_mask ≠ 0`. The branch is present in `pauli.c` (line 217) as dead code for
  completeness but is never executed.

**Path 3 — Anti-diagonal** (`x_mask ≠ 0`): amplitudes come in pairs (j, j'), where
j' = j XOR x_mask. For each pair:

    tmp     = ψ[j]
    ψ[j]   = ε(j') · ψ[j']
    ψ[j']  = ε(j)  · tmp

Since ε ∈ {±1, ±i}, each multiplication is implemented via the internal helper
`apply_ipow(z, exp)` (static inline, `pauli.c` line 151), which rotates a complex
number by i^exp using only component swaps and sign changes — no floating-point
multiplications.

**Pair enumeration** (critical implementation detail): to enumerate dim/2 unique pairs
without double-counting, insert a 0-bit at the **lowest set bit of x_mask** only
(`split_bit = __builtin_ctzll(x_mask)`). For k = 0..dim/2−1:

    j  = insertBit0(k, split_bit)   // j has 0 at split_bit
    j' = j XOR x_mask               // j' has 1 at split_bit

Chaining `insertBit0` at all x\_mask bit positions yields dim/2^{popcount(x\_mask)}
pairs — correct only for single-bit x\_masks. This was an actual bug caught by the
test suite.

**Phase computation:**

    parity_j  = popcount(j  AND z_mask) & 1
    exp_j     = (eta_exp + 2 * parity_j) & 3      // ε(j)  = i^{exp_j}
    exp_jp    = (exp_j  + 2 * y_parity) & 3       // ε(j') = i^{exp_jp}

where `eta_exp = popcount(y_mask) & 3` and `y_parity = popcount(y_mask) & 1` are
constants for the Pauli string.

**Validation:**
```c
GATE_VALIDATE(state && state->data, "...: null state or data pointer");
GATE_VALIDATE(P,                    "...: null Pauli pointer");
GATE_VALIDATE(state->type == PURE,  "...: expected PURE state");
GATE_VALIDATE(state->qubits >= 64 ||
              (P->x_mask | P->z_mask) >> state->qubits == 0,
              "...: Pauli mask exceeds state dimension");
```

The `state->qubits >= 64` guard avoids UB: shifting a 64-bit type by 64+ is undefined
in C17.

**OpenMP:** `#pragma omp parallel for if(dim >= OMP_THRESHOLD)` where
`OMP_THRESHOLD = 4096`. Pairs are independent — no data races.

### Identity Pauli Strings

When `x_mask = 0` (all labels are I or Z), the Pauli string is diagonal in the
computational basis: P|x⟩ = ε(x)|x⟩ with no bit flip.

- **Pure state:** `exp(−iθ/2 P)` applies a state-dependent phase `exp(−iθ/2 ε(x))` to each amplitude.
- **Density matrix:** ρ'_{xy} = exp(−iθ/2 (ε(x) − ε(y))) ρ_{xy}. The diagonal is unchanged (global phase cancels). Off-diagonal elements acquire a relative phase.

The all-identity string I⊗n (`x_mask = 0, z_mask = 0`) is a global phase on pure
states and a complete no-op on density matrices.

---

## Pauli Exponential

### Mathematical Definition

    U(θ) = exp(−iθ/2 · P) = cos(θ/2) I − i sin(θ/2) P

The adjoint is U(θ)† = U(−θ).

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

**Derivation:** P_{row, col} = ε(col) · δ_{row, col XOR x\_mask}. So
(P|ψ⟩)_x = ε(x') · ψ(x'). The inner update is a 2×2 rotation mixing `a` and `b`.
Since ε(x') ∈ {±1, ±i}, the multiplication by `−i·ε(x')` reduces to component
swaps and sign changes. When pop(y\_mask) is odd, eps\_x' = −eps\_x, halving the
number of phase evaluations.

**Edge case x_mask = 0:** For non-all-identity diagonals (z\_mask ≠ 0), iterate over
all x without pairing, applying `state->data[x] *= exp(−iθ/2 · ε(x))`.

**OpenMP:** parallel over pairs, `if (dim >= OMP_THRESHOLD)`.

### Algorithm: MIXED_PACKED and MIXED_TILED

For density matrices ρ → U(θ) ρ U(θ)†:

    ρ'_{xy} = c² ρ_{xy} + s² ε(x) ε*(y) ρ_{x',y'}
             − ics ε(x') ρ_{x',y} + ics ε*(y') ρ_{x,y'}

where c = cos(θ/2), s = sin(θ/2), x' = x XOR x_mask, y' = y XOR x_mask.

**Note:** ε(x')ε*(y') = ε(x)ε*(y) (the (−1)^{pop(y\_mask)} factors cancel in the
product). The cross terms do not share this cancellation.

**Four-element coupling:** Each output element ρ'_{xy} is a linear combination of
four input elements: ρ_{xy}, ρ_{x',y'}, ρ_{x',y}, ρ_{x,y'}. All four inputs must
be read before any output is written.

**Group enumeration:** Groups are indexed by representative (x, y) with x < x' and
y < y' (i.e., the lowest set bit of x\_mask is 0 in both x and y). The condition
y < y' prevents double-counting. Diagonal representatives (y = x) are included in the
main traversal with no special case.

**MIXED_PACKED:** Iterate over representative pairs, load all four elements (accessing
upper-triangle elements via conjugation), compute new values, write back. Diagonal
block guard: zero the imaginary part of ρ'_{xx} and ρ'_{x'x'} explicitly (numerical
hygiene; unitarity guarantees they are mathematically zero).

**MIXED_TILED:** Same kernel. Traversal mirrors `gate_tiled_h.c`: iterate over tile
pairs (tr, tc) with tr ≥ tc, then over within-tile element groups using
`dtile_read`/`dtile_write` from `gate_tiled.h`.

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

**Derivation:** (P|tmp⟩)_x = ε(x') · tmp(x') (only z = x' contributes to row x).
One loop over all amplitudes; each iteration: one phase lookup at x', one multiply-accumulate.

#### Mixed State — O(4^n), Direct Storage Traversal

Expanding `Tr(H_f · P · σ)` using P_{ab} = ε(b) δ_{a, b'} where b' = b XOR x_mask:

    Tr(H_f P σ) = Σ_{a,b} (H_f)_{b, a'} · ε(a) · σ_{ab}

Since only lower-triangle elements are stored and both H_f and σ are Hermitian:

```
S = 0

// Diagonal elements of σ  (one contribution each)
for each a = 0..d−1:
    S += lookup_hf(a, a') · ε(a) · σ_{aa}

// Off-diagonal stored elements of σ  (a > b, two contributions each)
for each (a, b) with a > b:
    S += lookup_hf(b, a') · ε(a) · σ_{ab}        // σ in lower tri
    S += lookup_hf(a, b') · ε(b) · conj(σ_{ab})  // σ upper tri by symmetry

return −Im S
```

where `lookup_hf(r, c)` returns `mu->data[packed_index(r, c)]` for r ≥ c, or its
conjugate for r < c.

**Correctness check:** For P = I^⊗n (x\_mask = 0, z\_mask = 0), S = Tr(H\_fᵀ σ) = 0
since H\_f and σ are both Hermitian and Tr(H\_f σ) is real. ✓

**MIXED_PACKED:** `packed_index(r, c)` = `c * d - c*(c+1)/2 + r` (matching STATE_MODULE.md).
**OpenMP:** `reduction(+:S_real, S_imag)` over the outer loop.

**MIXED_TILED:** Replace `packed_index` with `dtile_read(mu->data, r, c, n_tiles)`.
**OpenMP:** parallel over tile rows.

---

## Public Interface

### Channel Constructors

```c
// Pauli exponential channel. U(θ) = exp(−iθ/2 P). Generator G = P/2. has_gradient = 1.
// P is borrowed — caller must keep P alive for the lifetime of the channel.
// Returns NULL on OOM.
param_channel_t *pauli_exp_channel_new(const pauli_t *P);

// Noisy (Kraus) channel. grad_inner = NULL, has_gradient = 0.
// ctx is borrowed. Returns NULL on OOM.
param_channel_t *kraus_channel_new(
    void       (*apply)(state_t *state, double theta, const void *ctx),
    const void  *ctx
);

// Generic channel with user-supplied function pointers.
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

**Pauli exponential channel internals:** ctx holds a `const pauli_t *`. `apply`
dispatches on `state->type` to `pauli_exp_pure` / `pauli_exp_packed` / `pauli_exp_tiled`.

### Flat Parametrised Circuit

```c
// Create a circuit. Copies the channels pointer array (not the channel objects).
// Allocates params[n_channels] initialised to 0.0. Sets is_noisy automatically.
// Returns NULL on OOM.
pqc_t *pqc_new(param_channel_t **channels, int n_channels);

// Free the circuit struct and params array. Does NOT free channels or their ctx.
void pqc_free(pqc_t *circuit);

// Parameter access
void   pqc_set_params(pqc_t *circuit, const double *params); // copies n_channels values
double pqc_get_param(const pqc_t *circuit, int c);
void   pqc_set_param(pqc_t *circuit, int c, double val);
```

### LCU Structures

```c
// LCU layer. Copies fn_arr into internal array. Returns NULL on OOM.
lcu_layer_t *lcu_layer_new(void (**fn_arr)(state_t *), int K);
void lcu_layer_free(lcu_layer_t *layer);  // frees internal fn ptr array

// LCU circuit. Copies layers pointer array (layers themselves are borrowed).
// Computes N = ΠK_k. Returns NULL on OOM.
lcu_circuit_t *lcu_circuit_new(lcu_layer_t **layers, int L);
void lcu_circuit_free(lcu_circuit_t *circuit);  // does NOT free lcu_layer_t objects

// Decode linear multi-index I (0 ≤ I < circuit->N) into per-layer indices.
// idx_out[0] is least significant (changes fastest).
// Encoding: I = idx_out[0] + K_0*(idx_out[1] + K_1*(idx_out[2] + ...))
void lcu_multi_index(const lcu_circuit_t *circuit, int I, int *idx_out);
```

Applying multi-index I means sequentially applying
`layers[k].unitaries[idx_out[k]](state)` for k = 0..L−1 (layer 0 first).

---

## Interface with MEAS

MEAS is the exclusive consumer of `pqc_t`, `lcu_circuit_t`, and all circuit-level types.
Circuit exposes the types and constructors; all computation algorithms (expectation
value, gradient, moment matrix) live in MEAS.

**MEAS uses from Circuit:**
- `pqc_t` — for `expectation_value` and `gradient_adjoint`
- `lcu_circuit_t` — for `moment_matrix`
- `param_channel_t` function pointers — the 3-pass adjoint algorithm calls `apply` (with
  `−theta` for the adjoint direction) and `grad_inner` through these pointers
- `lcu_multi_index` — to decode column and row indices during moment matrix computation

**MEAS defines (outside Circuit's scope):**
- `diag_hamiltonian_t` and its operations
- QAOA channel constructors (cost and mixer) using `channel_new` from Circuit
- `qaoa_circuit_new`, `expectation_value`, `gradient_adjoint`, `moment_matrix`

**Required MEAS_MODULE.md updates** (migration tasks, not yet done):
- Replace all `param_channel_t`, `pqc_t`, `lcu_layer_t`, `lcu_circuit_t` definitions with `#include "circuit.h"`.
- `pqc_t.params` is now struct-owned (`pqc_new` allocates, `pqc_free` frees) — remove any caller-side `free(params)`.
- `lcu_circuit_t.layers` is now `lcu_layer_t **` (pointer array, borrowed), not `lcu_layer_t *` (flat array).
- `lcu_layer_t.unitaries` field (renamed from `apply` in earlier MEAS drafts).
- Remove the `pqc_free` declaration from MEAS to avoid duplicate declarations.

---

## Test Coverage

**Last verified:** 2026-02-27 — 30 tests, 0 failures, tolerance ε = 1e−12 (via `PRECISION` macro in `test/include/test.h`).

**Reference strategy:** Every `check_pauli_apply` call constructs the full 2^n × 2^n
Pauli matrix (`build_pauli_matrix`), multiplies it against the original statevector
(CBLAS `zgemv` via `zmv`), then compares element-wise against `pauli_apply_pure`
output. Linalg helpers from `test/utility/linalg.c` are linked into `test_circuit`.

### test_pauli.c — Implemented (30 tests)

#### Pauli Construction and Masks (5 tests)

| Test | Description |
|------|-------------|
| `test_pauli_masks_X` | n=1 X: x_mask=1, z_mask=0, y_mask=0 |
| `test_pauli_masks_Y` | n=1 Y: x_mask=1, z_mask=1, y_mask=1 |
| `test_pauli_masks_Z` | n=1 Z: x_mask=0, z_mask=1, y_mask=0 |
| `test_pauli_masks_str` | `pauli_from_str("IXYZ",4)`: x_mask=6, z_mask=3, y_mask=2; "IIII" (all zero); "XXXX" (x_mask=15) |
| `test_pauli_from_str_negative` | NULL string, invalid char 'A', embedded invalid char, length-too-short, length-too-long — all return NULL |

#### pauli_phase (1 test)

| Test | Description |
|------|-------------|
| `test_pauli_phase` | Spot-checks integer phase components for X(x=0,1), Y(x=0,1), Z(x=0,1) |

#### pauli_apply_pure — Single-qubit (4 tests)

All single-qubit tests run against all 6 standard states:
{|0⟩, |1⟩, |+⟩, |−⟩, |+y⟩, |−y⟩} (named `PSI1_0`, `PSI1_1`, `PSI1_P`, `PSI1_M`,
`PSI1_PI`, `PSI1_MI` in the source).

| Test | Pauli | x_mask | z_mask | split_bit | Code path |
|------|-------|--------|--------|-----------|-----------|
| `test_pauli_apply_I_q0` | I | 0 | 0 | — | Identity (early return) |
| `test_pauli_apply_X_q0` | X | 1 | 0 | 0 | Anti-diagonal |
| `test_pauli_apply_Y_q0` | Y | 1 | 1 | 0 | Anti-diagonal + Y-phase |
| `test_pauli_apply_Z_q0` | Z | 0 | 1 | — | Diagonal pure-Z |

#### pauli_apply_pure — Two-qubit strings (7 tests)

All n=2 tests run against {|00⟩, |01⟩, |10⟩, |11⟩, |++⟩, |+0⟩, |1+⟩}.

| Test | Pauli | x_mask | z_mask | split_bit | Notes |
|------|-------|--------|--------|-----------|-------|
| `test_pauli_apply_XX` | XX | 3 | 0 | 0 | Multi-bit x_mask, no phase |
| `test_pauli_apply_ZZ` | ZZ | 0 | 3 | — | Diagonal pure-Z, full mask |
| `test_pauli_apply_XZ` | XZ | 1 | 2 | 0 | Mixed X+Z types |
| `test_pauli_apply_YY` | YY | 3 | 3 | 0 | Multi-bit x_mask, η=2 → −1 phase |
| `test_pauli_apply_IZ` | IZ | 0 | 1 | — | Diagonal, non-trivial z position |
| `test_pauli_apply_IX` | IX | 2 | 0 | **1** | split_bit=1: first non-trivial insertBit0 |
| `test_pauli_apply_IY` | IY | 2 | 2 | **1** | split_bit=1 with Y-phase |

#### pauli_apply_pure — Three-qubit strings (6 tests)

All n=3 tests run against 8 computational basis states plus GHZ = (|000⟩+|111⟩)/√2,
W = (|001⟩+|010⟩+|100⟩)/√3, and |+++⟩.

| Test | Pauli | x_mask | z_mask | split_bit | Notes |
|------|-------|--------|--------|-----------|-------|
| `test_pauli_apply_XYZ` | XYZ | 3 | 6 | 0 | Mixed types; inline spot-check XYZ\|000⟩=i\|011⟩ |
| `test_pauli_apply_ZZZ` | ZZZ | 0 | 7 | — | Full diagonal |
| `test_pauli_apply_IXI` | IXI | 2 | 0 | **1** | split_bit=1, 3-qubit |
| `test_pauli_apply_XII` | XII | 4 | 0 | **2** | split_bit=2 |
| `test_pauli_apply_IZI` | IZI | 0 | 2 | — | Non-uniform diagonal |
| `test_pauli_apply_ZIZ` | ZIZ | 0 | 5 | — | Alternating diagonal |

#### pauli_apply_pure — Four-qubit strings (4 tests, MAXQUBITS=4)

All n=4 tests run against all 16 computational basis states plus PSI4_GHZ =
(|0000⟩+|1111⟩)/√2.

The "Pauli (labels[0..3])" column lists labels in little-endian order (labels[0] = qubit 0).
This differs from big-endian string notation: e.g., `{X,Z,Y,X}` written big-endian is "XYZX".

| Test | Pauli (labels[0..3]) | x_mask | z_mask | split_bit | Notes |
|------|----------------------|--------|--------|-----------|-------|
| `test_pauli_apply_4q` | {X,Z,Y,X} | 13 | 6 | 0 | Mixed types, η=1 (i-phase); y_mask=4 verified in test |
| `test_pauli_apply_4q_split1` | {I,X,I,I} | 2 | 0 | **1** | split_bit=1, 4-qubit |
| `test_pauli_apply_4q_split2` | {I,I,X,I} | 4 | 0 | **2** | split_bit=2 |
| `test_pauli_apply_4q_split3` | {I,I,I,X} | 8 | 0 | **3** | split_bit=3: all {0,1,2,3} covered |

#### Round-trip, involution, and edge cases (3 tests)

| Test | Description |
|------|-------------|
| `test_pauli_from_str_roundtrip` | `pauli_from_str("ZYX",3)` vs `pauli_new({X,Y,Z},3)`: identical masks and application results on all 8 basis states |
| `test_pauli_apply_involution` | P²=I for X, Y, Z (n=1), XX (n=2), XYZ (n=3) applied to \|0…0⟩ |
| `test_pauli_apply_identity_noop` | I⊗³ on \|+++⟩ is an exact no-op |

**split_bit coverage:** All positions {0, 1, 2, 3} exercised at n=4.

### test_pauli_exp.c — Pending

Tests for `pauli_exp_pure`, `pauli_exp_packed`, `pauli_exp_tiled`, and `grad_inner`
will be added when `pauli_exp.c` is implemented.

---

## Build Integration

The circuit library is a STATIC library in `src/CMakeLists.txt`. It links against the
shared `q` library (which bundles `state` and `gate`) rather than against separate
`state` and `gate` targets, because those are not exposed as standalone CMake targets.

```cmake
# In src/CMakeLists.txt (current — pauli.c only):
add_library(circuit STATIC
    "${CMAKE_CURRENT_SOURCE_DIR}/circuit/pauli.c"
    # pauli_exp.c and circuit.c to be added when implemented
)
target_include_directories(circuit
    PUBLIC  "${CMAKE_SOURCE_DIR}/include"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled"  # gate_tiled.h for pauli_exp.c
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

**Dependency rationale:**
- `q` is PUBLIC because `circuit.h` exposes `state_t *` in function signatures —
  consumers of `circuit.h` need `state.h` (and `gate.h`) transitively.
- `omp_compiler_flags` is PRIVATE — OpenMP is an implementation detail not in `circuit.h`.
- `linalg.c` is test-only (provides `zmv`); not linked into the `circuit` library.

---

## Deferred — Out of Scope

- **General Pauli-string Hamiltonian evolution** (ΣJ_k P_k): no decision on Trotterisation vs. direct exponentiation.
- **Pauli strings on n > 64 qubits.**
- **Multi-parameter channels** — each channel owns exactly one scalar parameter.
- **Shared-parameter channels** — caller bundles shared gates into one function pointer.
- **Gradient through Kraus channels** — caller uses parameter-shift or finite differences via repeated `expectation_value` in MEAS.
- **Parametrised LCU unitaries** — unitaries in `lcu_layer_t` have no free parameter.
- **Circuit optimisation / compilation passes.**
- **QASM / file-based circuit parsing.**
- **Hessians, natural gradients, quantum geometric tensor.**

---

## Open Questions

1. ~~**`circuit.h` pending types:** `param_channel_t`, `pqc_t`, `lcu_layer_t`, and
   `lcu_circuit_t` are specified here but not yet present in `include/circuit.h`. It is
   unspecified whether they will be added to the existing header or declared in a
   separate header (e.g., `pqc.h`). The split matters for MEAS's include directives.~~
   **RESOLVED (2026-02-27):** Audited `include/circuit.h` — none of the pending types
   are present, and no separate header (e.g., `pqc.h`) exists yet. The split decision
   remains open as a design choice, but the current state is confirmed: pending types
   are not declared anywhere in committed source.

2. **`OMP_THRESHOLD` for `pauli_exp`:** The spec uses `OMP_THRESHOLD = 512` for
   MIXED_PACKED/MIXED_TILED and `4096` for PURE. These values are stated as design
   intent but have not been validated by benchmarking. Whether they will share the same
   compile-time constant or be defined separately is unspecified. Note: `pauli.c`
   already defines `OMP_THRESHOLD = 4096` (line 25) for the PURE path; the
   MIXED threshold may need a separate symbol when `pauli_exp.c` is implemented.

3. **`pauli_exp_packed` upper-triangle conjugation:** The spec says to access
   upper-triangle elements of σ by conjugating the transposed index, but the exact
   interaction with `packed_index` when x' < x (i.e., the representative row flips to
   a higher index) is not spelled out. The tiled path similarly defers to `dtile_read`
   without specifying the diagonal-tile conjugation guard in detail.

4. **Error handling policy for channel/circuit constructors:** The spec states
   constructors return NULL on OOM, but does not specify whether they print to stderr,
   call a user-supplied error handler, or remain completely silent. This matters for
   integration with the test suite and MEAS error propagation. Note: `pauli_new` and
   `pauli_from_str` are completely silent on failure (no stderr output).

5. **`lcu_layer_t` function pointer closure model:** The spec states function pointers
   capture context "in a closure or pre-bound struct" but gives no concrete mechanism.
   It is unspecified whether the module will provide helpers for constructing such
   closures or leave that entirely to the caller.

6. **`grad_inner` signature for diagonal Pauli strings:** For x\_mask = 0, the
   gradient formula simplifies significantly, but the spec does not state whether a
   separate fast path will be implemented or whether the general formula handles it
   transparently.

7. **Naming collision risk:** The field rename from `apply` to `unitaries` in
   `lcu_layer_t` is documented, but it is not stated whether the old name ever appeared
   in any committed source file (vs. only in MEAS_MODULE.md). If it did, a grep-based
   migration is needed. Note: audited `include/circuit.h` and `src/circuit/pauli.c` —
   neither contains any `lcu_layer_t` definition. The old name `apply` has not appeared
   in any committed circuit source file. Q7 remains open only for MEAS_MODULE.md.
