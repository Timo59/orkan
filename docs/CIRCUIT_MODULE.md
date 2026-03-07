# Circuit Module: Technical Specification

**Module:** Parametrised Quantum Circuits and Pauli Operations
**Status:** Partially implemented — Pauli string construction, phase, pure-state application, and pure-state Pauli exponential complete. Packed/tiled exponential, `grad_inner`, and circuit types (`circuit.c`) pending.
**Depends on:** State module (`state.h`, `state_t`, `qubit_t`), Gate module (`gate.h`, `GATE_VALIDATE`)
**Consumed by:** MEAS module
**Last updated:** 2026-03-07 (circuit merged into q; labels removed from pauli_t; pauli_to_str added; error handling converted to GATE_VALIDATE/abort; dead branches and branching restructured; test oracle replaced with Kronecker-product reference independent of bitmasks; test states from gate-module test_mk_states_pure; 24 individual apply tests replaced by table-driven test_pauli_apply_all; 59 → 36 tests)

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
                               #   pauli_new/from_str/free/to_str/phase/apply_pure declared.
                               #   param_channel_t, pqc_t, lcu_* types NOT YET present.

src/circuit/
  pauli.c                      # DONE: pauli_new, pauli_from_str, pauli_free, pauli_to_str,
                               #       pauli_phase, pauli_apply_pure, pauli_exp_pure.
                               #       Compiled into libq (not a separate library).
  pauli_exp.c                  # PENDING: pauli_exp_{packed,tiled}, grad_inner
  circuit.c                    # PENDING: pqc_t, lcu_circuit_t constructors,
                               #          channel constructors

test/circuit/
  test_circuit.c               # DONE: test runner + main() — 36 tests, 0 failures
  test_pauli.c                 # DONE: all Pauli tests in one file:
                               #   — construction + phase (13 tests)
                               #   — pauli_apply_pure: single table-driven test_pauli_apply_all
                               #     (24 Paulis × all test_mk_states_pure states)
                               #   — round-trip and involution (3 tests)
                               #   — pauli_exp_pure (19 tests)
                               # PENDING: pauli_exp_{packed,tiled}, grad_inner

test/utility/
  gatemat.c / gatemat.h        # mat_pauli_str: builds full Pauli matrix via Kronecker
                               #   products of static 2×2 Pauli matrices.  Linked into
                               #   test_circuit as an independent oracle (no bitmasks).
  test_pure_states.c / .h      # test_mk_states_pure / test_rm_states_pure: generates
                               #   all cb/xb/yb/GHZ/W/Bell states for a given qubit count.
                               #   Used by every apply and exp test for state iteration.
```

---

## Current Status and TODO

| Item | Status |
|------|--------|
| `pauli_t` type, masks, construction | Complete |
| `pauli_phase` | Complete |
| `pauli_apply_pure` (PURE state) | Complete |
| `pauli_exp_pure` | Complete — implemented in `pauli.c`, declared in `circuit.h` |
| `pauli_exp_packed` (MIXED_PACKED) | Pending — algorithm fully specified |
| `pauli_exp_tiled` (MIXED_TILED) | Pending — algorithm fully specified |
| `grad_inner` for all state types | Pending — algorithm fully specified |
| `param_channel_t` type | Pending — type specified, not yet in `circuit.h` |
| `pqc_t` type and constructors | Pending — type specified, not yet in `circuit.h` |
| `lcu_layer_t`, `lcu_circuit_t` | Pending — types specified, not yet in `circuit.h` |
| Channel constructors (`pauli_exp_channel_new`, etc.) | Pending |
| `pauli_exp_pure` tests | Complete — 18 tests in `test_pauli.c` covering all branches, angles, properties |
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
| `<math.h>` | `cos`, `sin` (used in `pauli_exp_pure` in `pauli.c`; will also be used in pending `pauli_exp.c`) |

---

## Outputs — What This Module Produces

| Artifact | Description |
|----------|-------------|
| `libq.dylib` | Circuit sources compiled into the main `q` shared library |
| `include/circuit.h` | Public header; currently declares Pauli API only; will grow to include channel and circuit types |

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

Integer tags with no bit-field meaning. Mask construction in `pauli_new` uses a `switch` and does not depend on the numeric values.

### Pauli String

```c
#include <stdint.h>   // uint64_t

typedef struct {
    qubit_t  n_qubits; // matches state_t.qubits (unsigned char, max 255)
    // Precomputed bitmasks (set by pauli_new / pauli_from_str):
    uint64_t  x_mask;  // bit j set  iff  qubit j ∈ {X, Y}  (flip bit j)
    uint64_t  z_mask;  // bit j set  iff  qubit j ∈ {Y, Z}  (Z-phase factor)
    uint64_t  y_mask;  // bit j set  iff  qubit j = Y        (extra i factor)
} pauli_t;
```

**No label array.** The labels field has been removed. Use `pauli_to_str(P)` for
human-readable introspection; it reconstructs the per-qubit labels from the masks.

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
// Aborts (GATE_VALIDATE) on n_qubits > 64 or NULL labels with n_qubits > 0.
// n_qubits == 0 is valid: labels may be NULL; returns a Pauli with all masks zero.
pauli_t *pauli_new(const pauli_label_t *labels, qubit_t n_qubits);

// Parse from a string of characters 'I','X','Y','Z' of length n_qubits.
// Aborts (GATE_VALIDATE) on NULL s, length mismatch, or invalid character.
pauli_t *pauli_from_str(const char *s, qubit_t n_qubits);

// Free struct. Safe to call with NULL.
void pauli_free(pauli_t *p);

// Reconstruct big-endian character string from masks.
// s[0] = qubit n-1, s[n-1] = qubit 0. Caller must free() the result.
char *pauli_to_str(const pauli_t *P);
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

**Algorithm — two code paths dispatched by mask structure:**

**Path 1 — Diagonal** (`x_mask = 0`):
- Identity (`z_mask = 0`): immediate return.
- Pure Z-string (`z_mask ≠ 0`): `P|ψ⟩[i] = ε(i) · ψ[i]` for all i. ε(i) =
  (−1)^{popcount(i AND z\_mask)} ∈ {±1}. Implemented as a branchless real sign flip
  `data[i] *= 1.0 - 2.0*parity`. (`y_mask = 0` is guaranteed by construction when
  `x_mask = 0`, so no complex phase is needed here.)

**Path 2 — Anti-diagonal** (`x_mask ≠ 0`): amplitudes come in pairs (j, j'), where
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

**Diagonal case x_mask = 0:**
- Identity (z\_mask = 0): global phase e^{−iθ/2} applied to all amplitudes.
- Pure Z-string (z\_mask ≠ 0): per-amplitude multiply by e^{∓iθ/2} selected by
  parity. Two precomputed scalars avoid per-element `apply_ipow` dispatch.

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

**Last verified:** 2026-03-07 — 36 tests, 0 failures, tolerance ε = 1e−12 (via `PRECISION` macro in `test/include/test.h`).

### Reference Oracle

Both helpers receive a **big-endian Pauli string** and build two completely independent
paths — no shared data between the reference and the code under test:

**`check_pauli_apply(const char *str, unsigned n, const cplx_t *psi)`**

| Path | Steps |
|------|-------|
| Reference | `mat_pauli_str(str, n)` → 2ⁿ×2ⁿ matrix built by iterative Kronecker products of static 2×2 Pauli matrices (no bitmasks) → `zmv` (CBLAS) → M·psi |
| Under test | `pauli_from_str(str, n)` → `pauli_t` (bitmask construction) → `pauli_apply_pure` in-place |

**`check_pauli_exp(const char *str, unsigned n, const cplx_t *psi, double theta)`**

| Path | Steps |
|------|-------|
| Reference | Same Kronecker matrix → P\|psi⟩ = M·psi → ref\[k\] = c·psi\[k\] − i·s·(P\|psi⟩)\[k\] |
| Under test | `pauli_from_str(str, n)` → `pauli_t` → `pauli_exp_pure` in-place |

`mat_pauli_str` is implemented in `test/utility/gatemat.c` using four static 2×2
matrices and `zkron` from `test/utility/linalg.c`. It is linked into `test_circuit`
alongside `test/utility/test_pure_states.c`.

**Test state generation:** Every apply/exp test calls `test_mk_states_pure(n, &nvecs)`
which returns all states for `n` qubits:

| n | State count | Contents |
|---|-------------|---------|
| 1 | 6 | cb(2) + xb(2) + yb(2) |
| 2 | 16 | cb(4) + xb(4) + yb(4) + Bell(4) |
| 3 | 26 | cb(8) + xb(8) + yb(8) + GHZ + W |
| 4 | 50 | cb(16) + xb(16) + yb(16) + GHZ + W |

States are freed via `test_rm_states_pure(n, states)` at the end of each test.

### test_pauli.c — Implemented (36 tests)

#### Pauli Construction and Masks (10 tests)

| Test | Description |
|------|-------------|
| `test_pauli_masks_X` | n=1 X: x_mask=1, z_mask=0, y_mask=0 |
| `test_pauli_masks_Y` | n=1 Y: x_mask=1, z_mask=1, y_mask=1 |
| `test_pauli_masks_Z` | n=1 Z: x_mask=0, z_mask=1, y_mask=0 |
| `test_pauli_masks_str` | `pauli_from_str("IXYZ",4)`: x_mask=6, z_mask=3, y_mask=2; "IIII" (all zero); "XXXX" (x_mask=15) |
| `test_pauli_new_n0` | `pauli_new(NULL, 0)` → valid Pauli, all masks zero; `pauli_free` does not crash |
| `test_pauli_new_n64` | `pauli_new` with n=64, labels[63]=X, labels[62]=Y: x_mask bits 62–63 set; boundary of valid range |
| `test_pauli_free_null` | `pauli_free(NULL)` → no crash |
| `test_pauli_from_str_empty` | `pauli_from_str("", 0)` → valid 0-qubit Pauli, all masks zero |
| `test_pauli_mask_invariant` | For Y, YY, YYY, XYZ, IXYZ: verifies `(y_mask & x_mask) == y_mask` and `(y_mask & z_mask) == y_mask` |
| `test_pauli_to_str` | `pauli_to_str` round-trips: I/X/Y/Z (n=1), "IXYZ" (n=4), {X,Y,Z}→"ZYX" (n=3), "" (n=0) |

#### pauli_phase (3 tests)

| Test | Description |
|------|-------------|
| `test_pauli_phase` | Spot-checks integer phase components for X(x=0,1), Y(x=0,1), Z(x=0,1) — covers η=0 and η=1 |
| `test_pauli_phase_eta2` | YY (η=2, i²=−1): x=0→{−1,0}, x=1→{+1,0}, x=2→{+1,0}, x=3→{−1,0} |
| `test_pauli_phase_eta3` | YYY (η=3, i³=−i): x=0→{0,−1}, x=1→{0,+1}, x=3→{0,−1}, x=7→{0,+1} |

#### pauli_apply_pure — All Paulis (1 test)

**`test_pauli_apply_all`** iterates a static table of 24 Pauli strings and for each
calls `check_pauli_apply` over every state from `test_mk_states_pure`. Coverage:

| n | States per Pauli | Paulis | Code paths exercised |
|---|-----------------|--------|----------------------|
| 1 | 6 | I, X, Y, Z | Identity; diagonal Z; non-diagonal ±1/±i phase |
| 2 | 16 | XX, ZZ, ZX, YY, IZ, XI, YI, YZ | Multi-bit x_mask; non-uniform z_mask; η=1 cross-term; Bell states |
| 3 | 26 | ZYX, ZZZ, IXI, XII, IZI, ZIZ, YYY, III | split_bit ∈ {0,1,2}; η=1,3; GHZ and W states |
| 4 | 50 | XYZX, IIXI, IXII, XIII | split_bit ∈ {0,1,2,3}; n=4 GHZ |

**split_bit coverage:** All positions {0, 1, 2, 3} exercised.

#### Round-trip and involution (3 tests)

| Test | Description |
|------|-------------|
| `test_pauli_from_str_roundtrip` | `pauli_from_str("ZYX",3)` vs `pauli_new({X,Y,Z},3)`: identical masks; `pauli_apply_pure` on both gives same result for all 26 n=3 states |
| `test_pauli_apply_involution` | P²=I verified for {"X",1},{"Y",1},{"Z",1},{"XX",2},{"ZYX",3} over all states from `test_mk_states_pure` per case |
| `test_pauli_apply_involution_Z_nontrivial` | Z²=I verified on all 6 n=1 states (subsumes original 3-state check; confirms sign flip fires at index 1) |

### pauli_exp_pure — Implemented (19 tests)

Convention: `pauli_exp_pure(state, P, theta) = exp(−iθ/2 · P) = cos(θ/2)·I − i·sin(θ/2)·P`

**Reference oracle `check_pauli_exp(str, n, psi, theta)`:**
```
mat    = mat_pauli_str(str, n)          // Kronecker path — independent of bitmasks
Ppsi   = zmv(dim, mat, psi)             // CBLAS matrix-vector product
ref[k] = c·psi[k] − i·s·Ppsi[k]       // Euler decomposition reference
```
Compare `pauli_exp_pure` output element-wise within PRECISION.

#### Branch coverage and special angles (12 tests)

All exp tests use `test_mk_states_pure(n, &nvecs)` to iterate over all states.

| Test | Strings tested | States | Branch / property |
|------|---------------|--------|-------------------|
| `test_pauli_exp_identity_global_phase` | "I" (n=1), "II" (n=2) | all n=1/n=2 | Global phase e^{−iθ/2} |
| `test_pauli_exp_pure_z` | "Z","ZZ","ZZZ","IZI" | all n=1/2/3 | Diagonal branch; all parities |
| `test_pauli_exp_pure_nondiaganal_X` | "X" | all n=1 | η=0, y_parity=0 |
| `test_pauli_exp_pure_nondiaganal_Y` | "Y" | all n=1 | η=1 (real rotation matrix) |
| `test_pauli_exp_pure_nondiaganal_XX` | "XX" | all n=2 | η=0, n=2 non-diagonal |
| `test_pauli_exp_pure_nondiaganal_YY` | "YY" | all n=2 | η=2, y_parity=0 |
| `test_pauli_exp_pure_nondiaganal_YYY` | "YYY" | all n=3 | η=3, y_parity=1 — only i^3 path |
| `test_pauli_exp_pure_nondiaganal_XYZ` | "ZYX" | all n=3 | η=1, mixed z-phase + y-parity |
| `test_pauli_exp_theta_zero_is_identity` | "I","Z","X","Y","YYY" | all per n | θ=0: output = input |
| `test_pauli_exp_theta_pi_is_neg_i_pauli` | "X","Z" | all n=1 | θ=π: output = −i·P\|ψ⟩ via `pauli_apply_pure` reference |
| `test_pauli_exp_theta_2pi_is_neg_identity` | "I","Z","X","YYY" | all per n | θ=2π: output = −input |
| `test_pauli_exp_4q` | "XYZX" | all n=4 | n=4; angles 0, π/2, π, 2π |
| `test_pauli_exp_nondiaganal_split_bit` | "XI","IXI","XII","IYI" | all n=2/3 | split_bit ∈ {1,2}; all 9 angles |

#### Mathematical property tests (6 tests)

| Test | Property | Paulis | States |
|------|----------|--------|--------|
| `test_pauli_exp_unitarity` | ‖U(θ)\|ψ⟩‖² = 1 | "I","Z","ZZ","X","ZYX" | all per n; 5 angles |
| `test_pauli_exp_inverse` | U(θ)·U(−θ) = I | "I","Z","X","YYY","ZYX" | all per n; 4 angles |
| `test_pauli_exp_eigenstate_phases` | P\|v⟩=ε\|v⟩ → e^{−iεθ/2}\|v⟩ | "Z","X","ZZ" | inline-constructed eigenstates at θ=π/3 |
| `test_pauli_exp_theta_4pi_is_identity` | U(4π) = +I | "I","Z","X","YYY" | all per n |
| `test_pauli_exp_consistency_with_apply` | U(π) = −i·P | "I","Z","X","Y","XX","YY","YYY","ZYX" | all per n |
| `test_pauli_exp_non_unit_norm` | Norm² invariant | "X","Z","I" | [2,0] (non-normalised) |

---

## Build Integration

Circuit sources are compiled directly into the main `q` shared library. There is no
separate `circuit` CMake target.

```cmake
# In src/CMakeLists.txt:
set(Q_SOURCES
    ${STATE_SOURCES}
    "${CMAKE_CURRENT_SOURCE_DIR}/gate/gate.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/gate/gate_pure.c"
    ${GATE_PACKED_SOURCES}
    ${GATE_TILED_SOURCES}
    "${CMAKE_CURRENT_SOURCE_DIR}/circuit/pauli.c"
    # pauli_exp.c and circuit.c to be added when implemented
)
# q target gains PRIVATE include for gate/tiled (needed by pauli.c for insertBit0):
target_include_directories(q ... PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled")

# In test/CMakeLists.txt:
set(TEST_CIRCUIT_SRC
    "${CIRCUIT_DIR}/test_circuit.c"
    "${CIRCUIT_DIR}/test_pauli.c"
    "${UTILITY_DIR}/linalg.c"          # zmv() for matrix-vector products in oracle
    "${UTILITY_DIR}/gatemat.c"         # mat_pauli_str() Kronecker oracle
    "${UTILITY_DIR}/test_pure_states.c" # test_mk_states_pure / test_rm_states_pure
    "${UNITY_SRC}/unity.c"
)
target_link_libraries(test_circuit PRIVATE q)
```

**Rationale:** The circuit layer is part of the `q` module, not a separate library.
This keeps consumers to a single link target and avoids the layering complexity of a
`circuit` static lib that PUBLIC-links `q`.

**Test utility dependencies:**
- `linalg.c`: provides `zmv` (CBLAS matrix-vector multiply) and `zkron` (Kronecker product).
- `gatemat.c`: provides `mat_pauli_str` — builds the full Pauli matrix from 2×2 blocks via iterative `zkron`. This is the Kronecker-product reference oracle, completely independent of `pauli_t` bitmasks.
- `test_pure_states.c`: provides `test_mk_states_pure(n, &nvecs)` and `test_rm_states_pure(n, states)`. Generates the complete gate-module state set (cb + xb + yb + Bell/GHZ/W) for a given qubit count. Also used by `test_gate`; sharing it ensures both module test suites cover the same state space.

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

4. ~~**Error handling policy for channel/circuit constructors:** The spec states
   constructors return NULL on OOM, but does not specify whether they print to stderr,
   call a user-supplied error handler, or remain completely silent. This matters for
   integration with the test suite and MEAS error propagation. Note: `pauli_new` and
   `pauli_from_str` are completely silent on failure (no stderr output).~~
   **RESOLVED (2026-03-07):** `pauli_new`, `pauli_from_str`, and `pauli_to_str` now
   use `GATE_VALIDATE` (fprintf+exit) for all error conditions, matching the rest of
   the library. All programmer errors (bad arguments, OOM) are fatal. The test suite
   no longer tests NULL-return error paths for these functions.

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
