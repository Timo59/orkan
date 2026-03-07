# Circuit Module: Technical Specification

**Module:** Parametrised Quantum Circuits and Pauli Operations
**Status:** Partially implemented — Pauli string construction, phase, pure-state application, and pure-state Pauli exponential complete. Packed/tiled exponential, `grad_inner`, and circuit types (`circuit.c`) pending.
**Depends on:** State module (`state.h`, `state_t`, `qubit_t`), Gate module (`gate.h`, `GATE_VALIDATE`)
**Consumed by:** MEAS module

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
  circ.h                    # Public API — pauli_t, pauli_phase_t,
                            #   pauli_from_str/free/to_str/phase declared.
                            #   pauli_apply/pauli_exp dispatchers declared.
                            #   param_channel_t, pqc_t, lcu_* types NOT YET present.

src/circ/
  pauli.c                   # pauli_t constructors, pauli_phase, and
                            #   pauli_apply/pauli_exp dispatchers.
                            #   Compiled into libq (not a separate library).
  pauli_pure.c              # pauli_apply_pure, pauli_exp_pure (PURE backend).
  pauli_exp.c               # PENDING: pauli_exp_packed, pauli_exp_tiled, grad_inner
  circuit.c                 # PENDING: pqc_t, lcu_circuit_t constructors,
                            #          channel constructors
```

---

## API Reference

| Function | Description | Parameters | Status |
|----------|-------------|------------|--------|
| `pauli_from_str(s, n)` | Construct Pauli string from character string (big-endian) | `const char* s`, `qubit_t n_qubits` | Implemented |
| `pauli_free(p)` | Free Pauli string | `pauli_t* p` | Implemented |
| `pauli_to_str(P)` | Serialise Pauli string to character string | `const pauli_t* P` | Implemented |
| `pauli_phase(P, x)` | Compute phase ε(x) for P\|x⟩ = ε(x)\|x⊕x_mask⟩ | `const pauli_t* P`, `uint64_t x` | Implemented |
| `pauli_apply(state, P)` | Apply Pauli string in-place: \|ψ⟩ → P\|ψ⟩ | `state_t* state`, `const pauli_t* P` | PURE: Implemented · PACKED: Pending · TILED: Pending |
| `pauli_exp(state, P, θ)` | Apply Pauli exponential exp(−iθ/2 P) in-place | `state_t* state`, `const pauli_t* P`, `double theta` | PURE: Implemented · PACKED: Pending · TILED: Pending |
| `pauli_exp_channel_new(P)` | Construct parametric channel wrapping `pauli_exp` | `const pauli_t* P` | Pending |
| `channel_new(apply, grad, ctx, has_grad)` | Construct generic parametric channel | fn ptrs + `const void* ctx` | Pending |
| `kraus_channel_new(apply, ctx)` | Construct noisy Kraus channel | `apply` fn ptr, `const void* ctx` | Pending |
| `channel_free(ch)` | Free channel struct (does not free ctx) | `param_channel_t* ch` | Pending |
| `pqc_new(channels, n)` | Construct flat parametric circuit | `param_channel_t** channels`, `int n` | Pending |
| `pqc_free(circuit)` | Free circuit struct and params array | `pqc_t* circuit` | Pending |
| `pqc_set_params(circuit, params)` | Set all parameters | `pqc_t*`, `const double*` | Pending |
| `pqc_set_param(circuit, c, val)` | Set single parameter | `pqc_t*`, `int c`, `double val` | Pending |
| `pqc_get_param(circuit, c)` | Get single parameter | `const pqc_t*`, `int c` | Pending |
| `lcu_layer_new(fn_arr, K)` | Construct LCU layer from unitary fn array | `void(**)(state_t*)`, `int K` | Pending |
| `lcu_layer_free(layer)` | Free LCU layer | `lcu_layer_t* layer` | Pending |
| `lcu_circuit_new(layers, L)` | Construct LCU circuit | `lcu_layer_t** layers`, `int L` | Pending |
| `lcu_circuit_free(circuit)` | Free LCU circuit (does not free layers) | `lcu_circuit_t* circuit` | Pending |
| `lcu_multi_index(circuit, I, idx_out)` | Decode linear multi-index into per-layer indices | `const lcu_circuit_t*`, `int I`, `int*` | Pending |

---

## Inputs — What This Module Consumes

| Source | Items used |
|--------|-----------|
| `state.h` (State module) | `state_t`, `qubit_t`, `state_t.type` (`PURE`/`MIXED_PACKED`/`MIXED_TILED`), `state_t.data`, `state_t.qubits` |
| `gate.h` (Gate module) | `GATE_VALIDATE` macro; `insertBit0` (pair enumeration in `pauli_apply_pure`); `gate_idx_t`; `gate_tiled.h` (private, for `dtile_read`/`dtile_write` in pending `pauli_exp.c`) |
| `<stdint.h>` | `uint64_t` (bitmasks — exactly 64 bits on all platforms) |
| `<math.h>` | `cos`, `sin` (used in `pauli_exp_pure`; will also be used in pending `pauli_exp.c`) |

---

## Outputs — What This Module Produces

| Artifact | Description |
|----------|-------------|
| `libq.dylib` | Circuit sources compiled into the main `q` shared library |
| `include/circ.h` | Public header; currently declares Pauli API only; will grow to include channel and circuit types |

---

## Types

### Pauli String

```c
#include <stdint.h>   // uint64_t

typedef struct {
    qubit_t  n_qubits; // matches state_t.qubits (unsigned char, max 255)
    uint64_t  x_mask;  // bit j set  iff  qubit j ∈ {X, Y}
    uint64_t  z_mask;  // bit j set  iff  qubit j ∈ {Y, Z}
    uint64_t  y_mask;  // bit j set  iff  qubit j = Y
} pauli_t;
```

**Constraint:** `n_qubits ≤ 64`. By construction: `y_mask ⊆ x_mask` and `y_mask ⊆ z_mask`.

### Parametrised Channel

```c
typedef struct param_channel {
    void   (*apply)(state_t *state, double theta, const void *ctx);
    const void *ctx;      // channel-specific data, owned externally
} param_channel_t;
```

**`apply(state, θ, ctx)`**: Applies the channel in the forward direction.
- Unitary channel: state → U(θ) state (or U(θ)ρU(θ)† for density matrices).
- Noisy channel: state → Σ_a K_a state K_a† (full Kraus map; θ may be used or ignored).

**Adjoint:** For unitary channels U(θ)† = U(−θ). The backward pass calls
`apply(state, −theta, ctx)` directly; no separate adjoint pointer is stored.

**ctx lifetime:** Borrowed — the channel does not own it. The caller must keep `ctx`
alive for at least the lifetime of the channel.

### Flat Parametrised Circuit

```c
typedef struct {
    param_channel_t **channels;  // array of pointers, length n_channels
    int               n_channels;
    double           *params;    // one scalar per channel, length n_channels; owned
} pqc_t;
```

### LCU Layer

```c
typedef struct {
    void (**unitaries)(state_t *state);  // K unparametrised unitary function pointers
    int    K;
} lcu_layer_t;
```

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
      s[0]='I' → qubit 3  (I)
      s[1]='X' → qubit 2  (X)
      s[2]='Y' → qubit 1  (Y)
      s[3]='Z' → qubit 0  (Z)
      x_mask = 0b0110 = 6,  z_mask = 0b0011 = 3,  y_mask = 0b0010 = 2

### Phase Function

For Pauli string P applied to computational basis state |x⟩:

    P |x⟩ = ε(x) · |x XOR x_mask⟩

where the phase ε(x) ∈ {+1, −1, +i, −i} is:

    ε(x) = i^{popcount(y_mask)} · (−1)^{popcount(x AND z_mask)}

```c
typedef struct { int real; int imag; } pauli_phase_t;
// Examples: +1 = {1,0}, -1 = {-1,0}, +i = {0,1}, -i = {0,-1}

pauli_phase_t pauli_phase(const pauli_t *P, uint64_t x);
```

**Key identity:** ε(x') = ε(x) · (−1)^{popcount(y\_mask)}, where x' = x XOR x\_mask.

Corollaries:
- ε(x) · ε(x') = (−1)^{popcount(y\_mask)} (depends only on the Pauli string, not x).
- ε(x) · ε*(x') ∈ {+1, −1} (always real — P is Hermitian).
- If popcount(y\_mask) is even: ε(x') = ε(x). If odd: ε(x') = −ε(x).

The last corollary is used by `pauli_exp_pure`: when pop(y\_mask) is odd, only one phase
evaluation per pair is needed.

### Application Dispatcher

```c
// Apply Pauli string P to a quantum state in-place: |ψ⟩ → P|ψ⟩
// Fails fast via GATE_VALIDATE on null pointers or mask overflow.
// Dispatches to pauli_apply_pure / pauli_apply_packed / pauli_apply_tiled.
void pauli_apply(state_t *state, const pauli_t *P);
```

**PURE algorithm — two code paths dispatched by mask structure:**

**Path 1 — Diagonal** (`x_mask = 0`):
- Identity (`z_mask = 0`): immediate return.
- Pure Z-string (`z_mask ≠ 0`): `P|ψ⟩[i] = ε(i) · ψ[i]` for all i. Branchless real sign flip
  `data[i] *= 1.0 - 2.0*parity`.

**Path 2 — Anti-diagonal** (`x_mask ≠ 0`): amplitudes come in pairs (j, j'), where
j' = j XOR x_mask. For each pair:

    tmp     = ψ[j]
    ψ[j]   = ε(j') · ψ[j']
    ψ[j']  = ε(j)  · tmp

**Pair enumeration:** insert a 0-bit at the **lowest set bit of x_mask**
(`split_bit = __builtin_ctzll(x_mask)`). For k = 0..dim/2−1:

    j  = insertBit0(k, split_bit)   // j has 0 at split_bit
    j' = j XOR x_mask               // j' has 1 at split_bit

**OpenMP:** `#pragma omp parallel for if(dim >= OMP_THRESHOLD)` where `OMP_THRESHOLD = 4096`.

### Identity Pauli Strings

When `x_mask = 0` (all labels are I or Z), the Pauli string is diagonal in the
computational basis: P|x⟩ = ε(x)|x⟩ with no bit flip.

- **Pure state:** `exp(−iθ/2 P)` applies a state-dependent phase `exp(−iθ/2 ε(x))` to each amplitude.
- **Density matrix:** ρ'_{xy} = exp(−iθ/2 (ε(x) − ε(y))) ρ_{xy}. The diagonal is unchanged. Off-diagonal elements acquire a relative phase.

The all-identity string I⊗n (`x_mask = 0, z_mask = 0`) is a global phase on pure
states and a complete no-op on density matrices.

---

## Pauli Exponential

### Mathematical Definition

    U(θ) = exp(−iθ/2 · P) = cos(θ/2) I − i sin(θ/2) P

The adjoint is U(θ)† = U(−θ).

### Application Dispatcher

```c
// Apply U(theta) = exp(-i*theta/2 * P) in-place.
// Fails fast via GATE_VALIDATE on null pointers or mask overflow.
// Dispatches to pauli_exp_pure / pauli_exp_packed / pauli_exp_tiled.
void pauli_exp(state_t *state, const pauli_t *P, double theta);
```

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

        state->data[x]  = c*a - i*s*eps_x' * b
        state->data[x'] = c*b - i*s*eps_x  * a
```

**Diagonal case x_mask = 0:**
- Identity (z\_mask = 0): global phase e^{−iθ/2} applied to all amplitudes.
- Pure Z-string (z\_mask ≠ 0): per-amplitude multiply by e^{∓iθ/2} selected by parity.

**OpenMP:** parallel over pairs, `if (dim >= OMP_THRESHOLD)`.

### Algorithm: MIXED_PACKED and MIXED_TILED

For density matrices ρ → U(θ) ρ U(θ)†:

    ρ'_{xy} = c² ρ_{xy} + s² ε(x) ε*(y) ρ_{x',y'}
             − ics ε(x') ρ_{x',y} + ics ε*(y') ρ_{x,y'}

where c = cos(θ/2), s = sin(θ/2), x' = x XOR x_mask, y' = y XOR x_mask.

**Four-element coupling:** Each output element ρ'_{xy} is a linear combination of
four input elements: ρ_{xy}, ρ_{x',y'}, ρ_{x',y}, ρ_{x,y'}. All four must
be read before any output is written.

**Group enumeration:** Groups are indexed by representative (x, y) with x < x' and
y < y'. Diagonal representatives (y = x) are included in the main traversal.

**MIXED_PACKED:** Iterate over representative pairs, load all four elements (accessing
upper-triangle elements via conjugation), compute new values, write back.

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

#### Mixed State — O(4^n)

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
    S += lookup_hf(b, a') · ε(a) · σ_{ab}
    S += lookup_hf(a, b') · ε(b) · conj(σ_{ab})

return −Im S
```

**MIXED_PACKED:** `packed_index(r, c)` = `c * d - c*(c+1)/2 + r` (matching STATE_MODULE.md).
**MIXED_TILED:** Replace `packed_index` with `dtile_read(mu->data, r, c, n_tiles)`.

---

## Public Interface

### Channel Constructors

```c
// Pauli exponential channel. U(θ) = exp(−iθ/2 P). P is borrowed.
// Returns NULL on OOM.
param_channel_t *pauli_exp_channel_new(const pauli_t *P);

// Noisy (Kraus) channel. ctx is borrowed. Returns NULL on OOM.
param_channel_t *kraus_channel_new(
    void       (*apply)(state_t *state, double theta, const void *ctx),
    const void  *ctx
);

// Generic channel with user-supplied function pointers.
// ctx is borrowed. Returns NULL on OOM.
param_channel_t *channel_new(
    void   (*apply)(state_t *, double, const void *),
    const void *ctx
);

// Free the param_channel_t struct. Does NOT free ctx.
void channel_free(param_channel_t *ch);
```

### Flat Parametrised Circuit

```c
// Create a circuit. Copies the channels pointer array (not the channel objects).
// Allocates params[n_channels] initialised to 0.0. Returns NULL on OOM.
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

---

## Interface with MEAS

MEAS is the exclusive consumer of `pqc_t`, `lcu_circuit_t`, and all circuit-level types.
Circuit exposes the types and constructors; all computation algorithms (expectation
value, gradient, moment matrix) live in MEAS.

**MEAS uses from Circuit:**
- `pqc_t` — for `expectation_value` and `gradient_adjoint`
- `lcu_circuit_t` — for `moment_matrix`
- `param_channel_t.apply` function pointer — the 3-pass adjoint algorithm calls `apply` (with
  `−theta` for the adjoint direction) through this pointer
- `lcu_multi_index` — to decode column and row indices during moment matrix computation

**MEAS defines (outside Circuit's scope):**
- `diag_hamiltonian_t` and its operations
- QAOA channel constructors (cost and mixer) using `channel_new` from Circuit
- `qaoa_circuit_new`, `expectation_value`, `gradient_adjoint`, `moment_matrix`

**Required MEAS_MODULE.md updates** (migration tasks, not yet done):
- Replace all `param_channel_t`, `pqc_t`, `lcu_layer_t`, `lcu_circuit_t` definitions with `#include "circ.h"`.
- `pqc_t.params` is now struct-owned (`pqc_new` allocates, `pqc_free` frees) — remove any caller-side `free(params)`.
- `lcu_circuit_t.layers` is now `lcu_layer_t **` (pointer array, borrowed), not `lcu_layer_t *` (flat array).
- `lcu_layer_t.unitaries` field (renamed from `apply` in earlier MEAS drafts).

---

## Test Infrastructure

```
test/circ/
  test_circuit.c         # Test runner + main() — 36 tests, 0 failures
  test_pauli.c           # All Pauli tests: construction, phase, apply, exp

test/utility/
  gatemat.c / gatemat.h  # mat_pauli_str: builds full Pauli matrix via Kronecker
                         #   products of static 2×2 matrices. Independent oracle (no bitmasks).
  test_pure_states.c/.h  # test_mk_states_pure / test_rm_states_pure: generates
                         #   all cb/xb/yb/GHZ/W/Bell states for a given qubit count.
                         #   Shared with test_gate.
  linalg.c / linalg.h    # zmv (CBLAS matrix-vector), zkron (Kronecker product).
                         #   Used by gatemat.c to build the oracle matrix.
```

**Build (circuit tests only):**
```bash
cmake --preset debug
cmake --build --preset debug --target test_circuit
./cmake-build-debug/test/test_circuit
```

---

## Test Coverage

**Last verified:** 2026-03-07 — 35 tests, 0 failures, tolerance ε = 1e−12.

### Reference Oracle

**`check_pauli_apply(const char *str, unsigned n, const cplx_t *psi)`**

| Path | Steps |
|------|-------|
| Reference | `mat_pauli_str(str, n)` → 2ⁿ×2ⁿ matrix (Kronecker products, no bitmasks) → `zmv` → M·psi |
| Under test | `pauli_from_str(str, n)` → `pauli_t` → `pauli_apply` in-place |

**`check_pauli_exp(const char *str, unsigned n, const cplx_t *psi, double theta)`**

| Path | Steps |
|------|-------|
| Reference | Same Kronecker matrix → P\|psi⟩ = M·psi → ref\[k\] = c·psi\[k\] − i·s·(P\|psi⟩)\[k\] |
| Under test | `pauli_from_str(str, n)` → `pauli_t` → `pauli_exp` in-place |

**Test state generation:** Every apply/exp test calls `test_mk_states_pure(n, &nvecs)`:

| n | State count | Contents |
|---|-------------|---------|
| 1 | 6 | cb(2) + xb(2) + yb(2) |
| 2 | 16 | cb(4) + xb(4) + yb(4) + Bell(4) |
| 3 | 26 | cb(8) + xb(8) + yb(8) + GHZ + W |
| 4 | 50 | cb(16) + xb(16) + yb(16) + GHZ + W |

### Pauli Construction and Masks (10 tests)

| Test | Description |
|------|-------------|
| `test_pauli_masks_X` | n=1 X: x_mask=1, z_mask=0, y_mask=0 |
| `test_pauli_masks_Y` | n=1 Y: x_mask=1, z_mask=1, y_mask=1 |
| `test_pauli_masks_Z` | n=1 Z: x_mask=0, z_mask=1, y_mask=0 |
| `test_pauli_masks_str` | `pauli_from_str("IXYZ",4)`: x_mask=6, z_mask=3, y_mask=2; "IIII"; "XXXX" |
| `test_pauli_from_str_n64` | n=64, s[0]='X' (qubit 63), s[1]='Y' (qubit 62): boundary of valid range |
| `test_pauli_free_null` | `pauli_free(NULL)` → no crash |
| `test_pauli_from_str_empty` | `pauli_from_str("", 0)` → valid 0-qubit Pauli, all masks zero |
| `test_pauli_mask_invariant` | For Y, YY, YYY, XYZ, IXYZ: verifies `(y_mask & x_mask) == y_mask` and `(y_mask & z_mask) == y_mask` |
| `test_pauli_to_str` | Round-trips: I/X/Y/Z (n=1), "IXYZ" (n=4), {X,Y,Z}→"ZYX" (n=3), "" (n=0) |

### pauli_phase (3 tests)

| Test | Description |
|------|-------------|
| `test_pauli_phase` | Spot-checks for X(x=0,1), Y(x=0,1), Z(x=0,1) |
| `test_pauli_phase_eta2` | YY (η=2, i²=−1): all four x values |
| `test_pauli_phase_eta3` | YYY (η=3, i³=−i): x=0,1,3,7 |

### pauli_apply (1 test)

**`test_pauli_apply_all`** iterates 24 Pauli strings, for each calling `check_pauli_apply`
over every state from `test_mk_states_pure`. Round-trip and involution tests: 3 tests.

| n | States per Pauli | Paulis | Code paths |
|---|-----------------|--------|------------|
| 1 | 6 | I, X, Y, Z | Identity; diagonal Z; non-diagonal ±1/±i phase |
| 2 | 16 | XX, ZZ, ZX, YY, IZ, XI, YI, YZ | Multi-bit x_mask; Bell states |
| 3 | 26 | ZYX, ZZZ, IXI, XII, IZI, ZIZ, YYY, III | split_bit ∈ {0,1,2}; GHZ and W |
| 4 | 50 | XYZX, IIXI, IXII, XIII | split_bit ∈ {0,1,2,3}; n=4 GHZ |

### pauli_exp (19 tests)

Convention: `pauli_exp(state, P, theta) = exp(−iθ/2 · P) = cos(θ/2)·I − i·sin(θ/2)·P`

#### Branch coverage and special angles (13 tests)

| Test | Strings | Branch / property |
|------|---------|-------------------|
| `test_pauli_exp_identity_global_phase` | "I", "II" | Global phase e^{−iθ/2} |
| `test_pauli_exp_pure_z` | "Z","ZZ","ZZZ","IZI" | Diagonal branch |
| `test_pauli_exp_nondiaganal_X` | "X" | η=0, y_parity=0 |
| `test_pauli_exp_nondiaganal_Y` | "Y" | η=1 (real rotation) |
| `test_pauli_exp_nondiaganal_XX` | "XX" | η=0, n=2 |
| `test_pauli_exp_nondiaganal_YY` | "YY" | η=2, y_parity=0 |
| `test_pauli_exp_nondiaganal_YYY` | "YYY" | η=3, y_parity=1 — only i^3 path |
| `test_pauli_exp_nondiaganal_XYZ` | "ZYX" | η=1, mixed z-phase |
| `test_pauli_exp_theta_zero_is_identity` | "I","Z","X","Y","YYY" | θ=0: output = input |
| `test_pauli_exp_theta_pi_is_neg_i_pauli` | "X","Z" | θ=π: output = −i·P\|ψ⟩ |
| `test_pauli_exp_theta_2pi_is_neg_identity` | "I","Z","X","YYY" | θ=2π: output = −input |
| `test_pauli_exp_4q` | "XYZX" | n=4; angles 0, π/2, π, 2π |
| `test_pauli_exp_nondiaganal_split_bit` | "XI","IXI","XII","IYI" | split_bit ∈ {1,2} |

#### Mathematical property tests (6 tests)

| Test | Property | Paulis |
|------|----------|--------|
| `test_pauli_exp_unitarity` | ‖U(θ)\|ψ⟩‖² = 1 | "I","Z","ZZ","X","ZYX" |
| `test_pauli_exp_inverse` | U(θ)·U(−θ) = I | "I","Z","X","YYY","ZYX" |
| `test_pauli_exp_eigenstate_phases` | P\|v⟩=ε\|v⟩ → e^{−iεθ/2}\|v⟩ | "Z","X","ZZ" |
| `test_pauli_exp_theta_4pi_is_identity` | U(4π) = +I | "I","Z","X","YYY" |
| `test_pauli_exp_consistency_with_apply` | U(π) = −i·P | "I","Z","X","Y","XX","YY","YYY","ZYX" |
| `test_pauli_exp_non_unit_norm` | Norm² invariant | "X","Z","I" |

---

## Build Integration

Circuit sources are compiled directly into the main `q` shared library.

```cmake
# In src/CMakeLists.txt:
set(Q_SOURCES
    ${STATE_SOURCES}
    "${CMAKE_CURRENT_SOURCE_DIR}/gate/gate.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/gate/gate_pure.c"
    ${GATE_PACKED_SOURCES}
    ${GATE_TILED_SOURCES}
    "${CMAKE_CURRENT_SOURCE_DIR}/circ/pauli.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/circ/pauli_pure.c"
    # pauli_exp.c and circuit.c to be added when implemented
)
target_include_directories(q ... PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled")

# In test/CMakeLists.txt:
set(CIRCUIT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/circ")
set(TEST_CIRCUIT_SRC
    "${CIRCUIT_DIR}/test_circuit.c"
    "${CIRCUIT_DIR}/test_pauli.c"
    "${UTILITY_DIR}/linalg.c"
    "${UTILITY_DIR}/gatemat.c"
    "${UTILITY_DIR}/test_pure_states.c"
    "${UNITY_SRC}/unity.c"
)
target_link_libraries(test_circuit PRIVATE q)
```

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

1. **`OMP_THRESHOLD` for `pauli_exp`:** `OMP_THRESHOLD = 512` for MIXED_PACKED/MIXED_TILED
   and `4096` for PURE. These values are design intent, not validated by benchmarking.
   Whether they share one compile-time constant or use separate symbols is unspecified.

2. **`pauli_exp_packed` upper-triangle conjugation:** The exact interaction with
   `packed_index` when x' < x (representative row flips to a higher index) is not
   spelled out. The tiled path similarly defers to `dtile_read` without specifying
   the diagonal-tile conjugation guard in detail.

3. **`lcu_layer_t` function pointer closure model:** The spec states function pointers
   capture context "in a closure or pre-bound struct" but gives no concrete mechanism.

4. **`grad_inner` signature for diagonal Pauli strings:** For x\_mask = 0, the
   gradient formula simplifies significantly; it is unspecified whether a separate
   fast path will be implemented.

5. **`circuit.h` split decision:** Pending types (`param_channel_t`, `pqc_t`,
   `lcu_layer_t`, `lcu_circuit_t`) are not yet in `circ.h`. Whether they will be
   added to the existing header or a separate header is unresolved.
