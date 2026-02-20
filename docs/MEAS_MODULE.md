# MEAS Module: Parametrised Quantum Circuit Simulation & Gradient Computation

**Module:** Expectation Value Evaluation and Adjoint Gradient Computation
**Status:** Design complete — implementation pending
**Last Updated:** 2026-02-20

---

## Scope

This module evaluates expectation values ⟨H_C⟩ = Tr(H_C ρ(θ)) and their gradients ∂⟨H_C⟩/∂θ
for parametrised quantum circuits (PQCs). The framework is designed for arbitrary PQCs but the
**first implementation target is QAOA** (Quantum Approximate Optimization Algorithm).

Parameter optimisation is performed by an external library. This module exposes:
- `expectation_value(circuit, state, H_C, params)` → scalar
- `gradient(circuit, state, H_C, params)` → parameter vector

The module handles both **ideal** (unitary-only) and **noisy** (Kraus channel) circuits.

---

## Cost Hamiltonian Restriction

H_C is restricted to **diagonal Hamiltonians in the Z basis**:

    H_C = Σ_x h_x |x⟩⟨x|

This means:
- Expectation value: ⟨H_C⟩ = Σ_x h_x p_x, a weighted sum of populations
- For pure states: ⟨H_C⟩ = Σ_x h_x |ψ(x)|²  (no matrix multiply, O(2^n))
- For mixed states: ⟨H_C⟩ = Σ_x h_x ρ_{xx}  (diagonal of ρ only)
- Co-state initialisation: H_C|ψ_f⟩ is element-wise multiplication by h_x (free)

General non-diagonal observables are out of scope.

---

## Gradient Method: Adjoint Differentiation

**Chosen method: adjoint differentiation** (not parameter-shift rule).

### Rationale

The parameter-shift rule requires 2 full circuit evaluations per parameter. For a circuit with
2p parameters (QAOA with p layers) this gives O(p) evaluations. Adjoint differentiation requires
O(1) evaluations regardless of circuit depth — at the cost of memory or extra passes.

Additionally, adjoint differentiation places no restriction on the generator's eigenvalue spectrum,
accommodating Hamiltonian evolution gates with arbitrary generators.

### Mathematical Foundation

For a parametrised gate at position k, U_k(θ_k) = exp(−iθ_k G_k), embedded in circuit
U_f · U_k(θ_k) · V:

    ∂⟨H_C⟩/∂θ_k = −2 Im Tr(H_f^(k) · G_k · ρ_k)

where:
- ρ_k = V ρ_0 V†  is the state just before gate k (Schrödinger picture, forward)
- H_f^(k) = U_f† H_C U_f  is H_C pulled back to position k (Heisenberg picture, backward)

Both ρ_k and H_f^(k) are density-matrix-shaped objects (2^n × 2^n Hermitian operators).
H_C starts diagonal but becomes a full Hermitian matrix after mixer layers.
H_f is stored in the same MIXED_PACKED or MIXED_TILED layout as ρ.

---

## Circuit Abstraction: Parametrised Channel

### Definition

The fundamental unit is a **parametrised channel**: a struct grouping all gates that share
one circuit parameter, together with optional Kraus noise operators.

Each channel exposes exactly three operations to the gradient algorithm:

1. `apply(state, θ)`          — forward application (unitary part then noise)
2. `apply_adjoint(state, θ)`  — backward propagation of co-state
3. `grad_inner(μ, tmp, θ)`    — computes −2 Im Tr(H_f · G · σ), the gradient contribution

For the noisy case the forward application is:  σ = U(θ) ρ U(θ)†  then  ρ' = Σ_m K_m σ K_m†
The backward propagation is the adjoint channel:  H_f ← Σ_m K_m† H_f K_m  then  H_f ← U†H_fU
The gradient inner product uses σ = U(θ) ρ_k U(θ)†  (state after unitary, before noise).

### Parameter Sharing Within a Channel

All gates within one channel share one scalar parameter θ. The total gradient contribution
of the channel is the sum over all constituent gates:

    ∂⟨H_C⟩/∂θ = Σ_{gates in channel} (−2 Im Tr(H_f^(k) G_k σ_k))

For QAOA cost layers (commuting RZZ gates), this reduces to one evaluation via the combined
generator G = H_C (see QAOA section below).

Parameters are NOT shared across channels. Each channel owns exactly one parameter.

### Circuit

A circuit is a flat array of N_ch parametrised channel pointers and a flat parameter vector θ.
For QAOA with p layers: N_ch = 2p channels (p cost + p mixer).

Fixed (non-parametrised) operations are represented as channels with no gradient contribution
(θ is ignored, `grad_inner` returns 0).

### LCU Channel (Fixed, No Gradient)

A **linear combination of unitaries** channel: M = Σ_i c_i U_i, applied to pure states only.
It is non-parametrised — its coefficients are determined externally and not subject to gradient
search. It still participates in the backward co-state propagation.

**GAP — see Open Gaps section.**

---

## Gradient Algorithm

### State Type × Circuit Type

| Circuit     | State type   | Gradient algorithm                     | Memory              |
|-------------|--------------|----------------------------------------|---------------------|
| Ideal       | PURE         | 3-pass (see below)                     | 3 state vectors     |
| Ideal       | MIXED        | 3-pass (same logic, density matrix)    | 2 density matrices  |
| Noisy       | MIXED        | Forward (store p DMs) + backward on H_f| O(p) density matrices|

Never: PURE state with Kraus channels. Kraus noise destroys purity; the state must be MIXED.

### 3-Pass Algorithm (Ideal Circuits, Pure and Mixed)

Exploits unitarity U_k U_k† = I to propagate both the state and co-state forward in the same
pass. No intermediate states are stored. Total cost: 3 × L gate applications, where L is the
number of gates in the circuit.

**Pass 1 — forward, compute final state:**
```
|tmp⟩ = |ψ_0⟩
for k = 1 to L:
    |tmp⟩ = U_k(+θ_k)|tmp⟩
|ψ_f⟩ = |tmp⟩
```

**Initialise co-state:**
```
E = Σ_x h_x |tmp(x)|²            # expectation value (free for diagonal H_C)
|μ⟩ = H_C |ψ_f⟩                  # element-wise multiply by h_x
```

**Pass 2 — backward, initialise co-state at circuit entrance:**
```
for k = L, L-1, ..., 1:
    |μ⟩ = U_k†(θ_k)|μ⟩           # adjoint gate = gate with negated angle
|μ_0⟩ = |μ⟩
```

**Pass 3 — forward, accumulate gradients:**
```
|tmp⟩ = |ψ_0⟩
|μ⟩   = |μ_0⟩
for k = 1 to L:
    |tmp⟩ = U_k(+θ_k)|tmp⟩       # same forward gate as Pass 1
    |μ⟩   = U_k(+θ_k)|μ⟩         # same forward gate (unitarity peels off one U_k†)
    grad_k += −2 Im ⟨μ|G_k|tmp⟩  # inner product; G_k applied to |tmp⟩ temporarily
```

For MIXED ideal circuits, replace state vectors with density matrices and apply gates as
U(·)U†. The 3-pass logic is identical; |μ⟩ becomes a density-matrix-shaped H_f.

### Noisy Circuit Algorithm (Mixed States, Kraus Channels)

The 3-pass trick requires invertibility (U_k U_k† = I). Kraus channels are generally not
invertible, so the trick breaks. Instead:

**Pass 1 — forward, store intermediate density matrices:**
```
ρ^(0) = ρ_0
for c = 1 to N_ch:
    σ^(c) = U_c(+θ_c) ρ^(c-1) U_c†(+θ_c)    # unitary part
    ρ^(c) = Σ_m K_m^(c) σ^(c) K_m^(c)†        # noise part
    store ρ^(c-1)                               # store entrance state
ρ_f = ρ^(N_ch)
E = Σ_x h_x (ρ_f)_{xx}
```

**Pass 2 — backward, propagate H_f and accumulate gradients:**
```
H_f = diag(h_x)                   # initialise as diagonal H_C
for c = N_ch, N_ch-1, ..., 1:
    H_f = Σ_m K_m^(c)† H_f K_m^(c)            # adjoint noise channel
    σ^(c) = U_c(θ_c) ρ^(c-1) U_c†(θ_c)        # recompute from stored ρ^(c-1)
    grad_c = −2 Im Tr(H_f · G_c · σ^(c))       # gradient contribution
    H_f = U_c†(θ_c) H_f U_c(θ_c)              # adjoint unitary
```

Memory cost: N_ch density matrices stored at channel entrances.
For QAOA (N_ch = 2p, n qubits): 2p × 4^n × 16 bytes.
At n = 10, p = 5: ~320 MB. At n = 12, p = 5: ~5 GB. Hard limit around n = 11–12.

---

## QAOA: First Implementation Target

### Circuit Structure

Depth-p QAOA on n qubits: 2p channels, 2p parameters (γ_1…γ_p, β_1…β_p).

```
Channel 1: cost layer,  param γ_1
Channel 2: mixer layer, param β_1
Channel 3: cost layer,  param γ_2
Channel 4: mixer layer, param β_2
...
Channel 2p-1: cost layer,  param γ_p
Channel 2p:   mixer layer, param β_p
```

### Cost Channel

Implements U_C(γ) = exp(−iγ H_C) where H_C = Σ_{(i,j)∈E} w_{ij} Z_i Z_j / 2.

Since H_C is diagonal, this is **element-wise phase multiplication** in O(2^n):
    ψ(x) ← ψ(x) · exp(−iγ h_x)   where h_x = Σ_{(i,j)∈E} w_{ij} z_i(x) z_j(x) / 2

Generator: G_C = H_C (diagonal). Gradient inner product:
    −2 Im ⟨μ|H_C|tmp⟩ = −2 Im Σ_x h_x μ*(x) tmp(x)

No RZZ decomposition required. No Trotter approximation.

### Mixer Channel

Implements U_B(β) = Π_j exp(−iβ X_j / 2) = exp(−iβ Σ_j X_j / 2).

All single-qubit Rx(β) gates commute ([X_i, X_j] = 0 for i ≠ j). Applied as n sequential
Rx(β) gate calls on each qubit.

Generator: G_B = Σ_j X_j / 2. Gradient inner product:
    −2 Im ⟨μ|G_B|tmp⟩ = Σ_j (−Im ⟨μ|X_j|tmp⟩)

Each term requires one X gate applied to |tmp⟩, then a dot product with |μ⟩.

### Initial State

|ψ_0⟩ = |+⟩^⊗n (uniform superposition) via `state_plus()` — already implemented.

### Required Gates (Implementation Status)

| Gate/Operation                  | Status in gate module          | Needed for    |
|---------------------------------|-------------------------------|---------------|
| Rx(θ)                           | ✓ All backends                | Mixer channel |
| Rz(θ)                           | ✓ All backends                | Misc          |
| Diagonal phase multiply exp(−iγ h_x) | ✗ New operation required  | Cost channel  |
| RZZ(θ, i, j)                    | ✗ Deferred in GATE_MODULE.md  | (optional if cost channel uses diagonal phase) |

The cost channel should be implemented as a direct diagonal phase multiply, not as
decomposed RZZ gates. RZZ is needed only if individual two-qubit gradients are required
(not the case for QAOA, which treats the whole cost layer as one channel).

---

## Data Structures (Planned)

```c
// Diagonal cost Hamiltonian
typedef struct {
    double   *h;       // h[x] = diagonal element for basis state x, length 2^n
    qubit_t   qubits;
} diag_hamiltonian_t;

// Abstract parametrised channel — function pointer interface
typedef struct param_channel {
    void   (*apply)        (state_t *state, double theta, const void *ctx);
    void   (*apply_adjoint)(state_t *state, double theta, const void *ctx);
    double (*grad_inner)   (const state_t *mu, const state_t *tmp,
                            double theta, const void *ctx);
    const void *ctx;       // channel-specific data (qubit indices, Hamiltonian, Kraus ops…)
    int    has_gradient;   // 0 for fixed channels (LCU, state prep)
} param_channel_t;

// Circuit: flat array of channels + parameter vector
typedef struct {
    param_channel_t **channels;
    int               n_channels;
    double           *params;      // length n_channels (one param per channel)
} pqc_t;
```

---

## Planned API

```c
// Evaluate expectation value ⟨H_C⟩ = Σ_x h_x p_x
double expectation_value(
    const pqc_t            *circuit,
    const state_t          *initial,
    const diag_hamiltonian_t *H_C
);

// Compute full gradient vector via adjoint differentiation
// grad_out[c] = ∂⟨H_C⟩/∂params[c] for each channel c
void gradient_adjoint(
    const pqc_t            *circuit,
    const state_t          *initial,
    const diag_hamiltonian_t *H_C,
    double                 *grad_out   // length circuit->n_channels
);

// QAOA-specific constructors
pqc_t *qaoa_circuit_new(
    int p,                             // number of layers
    const diag_hamiltonian_t *H_C,     // cost Hamiltonian
    int n_qubits
);
void pqc_free(pqc_t *circuit);
```

---

## Open Gaps

The following questions are **unresolved** and must be answered before implementation proceeds.
The next agent should address each gap explicitly.

### GAP 1 — LCU Channel Backward Pass (Blocking)

A linear combination of unitaries channel M = Σ_i c_i U_i is non-parametrised and applies to
pure states. "Norm preservation" is stated, but its exact meaning is unresolved:

**Case A:** M is constrained to be exactly unitary (e.g., by optimising the coefficients c_i
under a unitarity constraint). Then `apply_adjoint` is simply M† = Σ_i c_i* U_i†. The channel
slots cleanly into the existing framework.

**Case B:** M is not unitary; norm is preserved by post-hoc renormalisation:
    |ψ_out⟩ = M|ψ⟩ / ‖M|ψ⟩‖
Then `apply_adjoint` is the adjoint of the tangent map (state-dependent):
    d(N)|_|ψ⟩ (|δ⟩) = P⊥_{M|ψ⟩} M|δ⟩ / ‖M|ψ⟩‖
This requires |ψ⟩ and M|ψ⟩ to be stored during the forward pass and makes the backward pass
nonlinear. The implementation is substantially more complex.

**Resolution needed:** Which case applies?

### GAP 2 — Kraus Noise Model (Blocking for Noisy QAOA)

The adjoint channel backward pass E†(H_f) = Σ_m K_m† H_f K_m is designed, but the concrete
Kraus operators for the noise model are unspecified. Required:
- Which noise channels are needed? (depolarising, amplitude damping, dephasing, phase flip, …)
- Are Kraus operators gate-type-dependent (e.g., two-qubit depolarising on CX), or applied
  uniformly after every gate?
- Are noise parameters fixed inputs or also read from the channel struct?

### GAP 3 — Diagonal Phase Multiply (New Gate Operation)

The cost channel exp(−iγ H_C) for diagonal H_C requires applying element-wise phase factors
ψ(x) ← ψ(x) · exp(−iγ h_x). This operation does not exist in the current gate infrastructure.
It must be added to the gate module before the cost channel can be implemented.

For mixed states: ρ_{xy} ← ρ_{xy} · exp(−iγ(h_x − h_y)). This is a diagonal conjugation on
the density matrix — also new, but structurally similar to existing diagonal gate patterns.

### GAP 4 — Gradient Inner Product for Mixed States

The gradient contribution Tr(H_f · G · σ) for density matrix σ and generator G requires
computing a Hilbert-Schmidt inner product with an inserted single- or two-qubit operator.
For the mixer (G = Σ_j X_j/2): Tr(H_f X_j σ) per qubit j. Efficient computation strategy
(partial trace, direct element access) is not yet specified.

### GAP 5 — General Pauli String Hamiltonian Evolution (Deferred)

The framework is intended to support exp(−iθ Σ_α c_α P_α) for non-diagonal, potentially
non-commuting Pauli string Hamiltonians. QAOA does not require this. Two options:
- Trotterisation: exp(−iθH) ≈ Π_α exp(−iθ c_α P_α) — approximate, uses existing gates
- Direct matrix exponentiation: exact, O(8^n) or O(2^{3n})

No decision made. Deferred until after QAOA implementation.

### GAP 6 — RSWAP Gate

RSWAP is listed as a parametrised gate type in the general framework. Definition (is it
exp(−iθ SWAP/2)?) and role in the first implementation are unspecified. Not required for
standard QAOA. Defer or clarify.

### GAP 7 — Memory Limit for Noisy Simulation

The noisy gradient algorithm stores N_ch = 2p density matrices. At n = 12 qubits, p = 5
layers: ~5 GB. A hard limit around n = 11–12 is implied. The acceptable qubit / layer
range for the target hardware has not been stated, and no checkpoint strategy has been
specified for larger systems.

---

## Deferred (Not in Scope for First Implementation)

- Hessian / second-order gradients
- Natural gradient (quantum geometric tensor / Fisher information metric)
- General non-diagonal observables
- Multi-parameter channels (gates with more than one independent parameter)
- RSWAP and other non-standard parametrised gates
- Trotter-based Hamiltonian simulation
