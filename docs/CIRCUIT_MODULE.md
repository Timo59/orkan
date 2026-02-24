# Circuit Module: Technical Specification

**Module:** Parametrised Quantum Circuits and Pauli Operations
**Status:** Design complete — implementation pending
**Depends on:** State module, Gate module
**Consumed by:** MEAS module
**Last updated:** 2026-02-23

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
  circuit.h                    # Public API: all types and declarations

src/circuit/
  pauli.c                      # Pauli string construction, phase function, application
  pauli_exp.c                  # Pauli exponential: PURE, MIXED_PACKED, MIXED_TILED
  circuit.c                    # pqc_t, lcu_circuit_t constructors/destructors,
                               # built-in channel constructors

test/circuit/
  test_circuit.c               # Test runner + main()
  test_pauli.c                 # Pauli string unit tests
  test_pauli_exp.c             # Pauli exponential unit tests
```

**Build (circuit tests only):**
```bash
cd cmake-build-debug
cmake --build . --target test_circuit
./test/test_circuit
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

// Free labels array and struct.
void pauli_free(pauli_t *p);
```

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

### test_pauli.c

| Test | Description |
|------|-------------|
| `test_pauli_from_str` | Parse "IXYZ", "XXXX", "IIII"; verify labels and masks |
| `test_pauli_masks` | For n=1..4 and all 4^n Pauli strings: verify x_mask, y_mask, z_mask against brute-force label scan |
| `test_pauli_phase` | For n ≤ 3: enumerate all (P, x) pairs; compare pauli_phase against direct matrix element extraction via cblas_zgemv |
| `test_pauli_invalid` | pauli_from_str with bad characters returns NULL |

### test_pauli_exp.c

Single-qubit identity checks (compare against gate module):
| Test | Description |
|------|-------------|
| `test_pauli_exp_X` | exp(−iθ/2 X) ≡ Rx(θ); all 3 backends, all standard θ values |
| `test_pauli_exp_Y` | exp(−iθ/2 Y) ≡ Ry(θ); all 3 backends |
| `test_pauli_exp_Z` | exp(−iθ/2 Z) ≡ Rz(θ) up to global phase; all 3 backends. For MIXED backends the global phase cancels (ρ' = UρU†) so results are identical. For PURE, verify that `pauli_exp_pure(Z)` produces exp(−iθ/2)·Rz(θ)·ψ elementwise (the factor exp(−iθ/2) is state-independent) |
| `test_pauli_exp_I` | all-identity string: verify no-op on density matrices; global phase on PURE |

Multi-qubit:
| Test | Description |
|------|-------------|
| `test_pauli_exp_XX` | 2-qubit XX string; reference via direct matrix exponentiation |
| `test_pauli_exp_ZZ` | 2-qubit ZZ string (diagonal, no bit flips) |
| `test_pauli_exp_XZ` | Mixed-kind 2-qubit string |
| `test_pauli_exp_random` | Random n=3 Pauli strings and angles; reference via scipy.linalg.expm equivalent (dense matrix multiply) |

Gradient:
| Test | Description |
|------|-------------|
| `test_grad_inner` | Compare grad_inner output to finite-difference ∂⟨H_C⟩/∂θ for 1-qubit and 2-qubit Pauli strings; all 3 backends |

**Test states:** Same set as gate module (Z/X/Y-eigenbasis, random density matrices, fixed seed).
**Tolerance:** ε = 1e−12.
**System sizes:** n = 1, 2, 3; tiled: n = 5, 6 to exercise cross-tile paths.

---

## Build Integration

```cmake
add_library(circuit STATIC
    src/circuit/pauli.c
    src/circuit/pauli_exp.c
    src/circuit/circuit.c
)
target_include_directories(circuit PUBLIC include)
target_link_libraries(circuit
    PUBLIC  state          # state_t visible through circuit.h
    PRIVATE gate           # gate.h (insertBit0, pack_idx) used by pauli_exp.c
)
target_include_directories(circuit PRIVATE src/gate/tiled)
# gate_tiled.h lives in src/gate/tiled/ which is a PRIVATE gate header not exposed
# through gate's PUBLIC include path (include/). The PRIVATE include above gives
# pauli_exp.c direct access without leaking gate_tiled.h to circuit's consumers.

if(ENABLE_OPENMP)
    target_link_libraries(circuit PRIVATE OpenMP::OpenMP_C)
endif()

if(BUILD_TESTS)
    add_executable(test_circuit
        test/circuit/test_circuit.c
        test/circuit/test_pauli.c
        test/circuit/test_pauli_exp.c
    )
    target_link_libraries(test_circuit PRIVATE circuit gate state)
    add_test(NAME circuit COMMAND test_circuit)
endif()
```

`state` is a PUBLIC dependency because `param_channel_t` and `pqc_t` expose `state_t *`
in their function pointer signatures — consumers of `circuit.h` need `state.h` transitively.

`gate` is a PRIVATE dependency because `pauli_exp.c` uses `gate_tiled.h` inline helpers
(`tile_off`, `elem_off`, `dtile_read`, `dtile_write`) for the tiled backend traversal.
These helpers do not appear in `circuit.h`.

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
