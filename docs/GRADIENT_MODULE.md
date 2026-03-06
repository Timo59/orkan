# Gradient Computation for Parametrised Quantum Circuits

## Overview

This document summarises the theory and proposed implementation for computing gradients of parametrised quantum circuits, covering both the pure-state (state-vector) and mixed-state (density matrix) settings.

---

## Circuit Structure

- Elementary parametrised gates: $R_x(\theta) = e^{-i\theta X/2}$, $R_y(\theta) = e^{-i\theta Y/2}$, $R_z(\theta) = e^{-i\theta Z/2}$.
- A parameter $\theta_k$ may appear in multiple gates (parameter sharing).
- Gates may be grouped into layers (QAOA-style) where each layer has a unique parameter symbol.
- Cost function: $f(\boldsymbol{\theta}) = \langle O \rangle = \langle\psi(\boldsymbol{\theta})|O|\psi(\boldsymbol{\theta})\rangle$ or $\operatorname{tr}(O\rho(\boldsymbol{\theta}))$.

---

## Adjoint Differentiation (Primary Method)

Reverse-mode AD adapted to quantum circuits. Cost is $O(L \cdot 2^n)$ **independent of the number of parameters**.

### Gradient Identity (Pure States)

Define forward states $|\phi_j\rangle = U_j\cdots U_1|0\rangle$ and adjoint states $|\lambda_j\rangle = U_{j+1}^\dagger\cdots U_L^\dagger\,O\,|\psi\rangle$. For gate $j$ with Pauli generator $G_j$ and parameter $\theta_k$:

$$\frac{\partial f}{\partial \theta_k}\bigg|_{\text{gate }j} = \operatorname{Im}\!\left[\langle\lambda_j|G_j|\phi_j\rangle\right]$$

**Proof sketch.** Differentiate $f$ using $dU_j/d\theta = -\frac{i}{2}G_j U_j$. The two resulting terms are complex conjugates of each other (since $O$ is Hermitian). Their difference is $2i\operatorname{Im}(\cdot)$; the factor $-\frac{i}{2}\cdot 2i = +1$ leaves $\operatorname{Im}[\langle\lambda_j|G_j|\phi_j\rangle]$ with no residual prefactor.

**Requirement:** $O$ must be Hermitian. Projectors, Pauli strings, and Hamiltonians are all valid.

### Algorithm

```
ADJOINT_GRADIENT(circuit, θ, O):

  # Forward pass
  |ψ⟩ ← |0...0⟩
  for j = 1 to L:
      |ψ⟩ ← U_j(θ)|ψ⟩

  # Initialise adjoint state
  |λ⟩ ← O|ψ⟩

  # Backward pass
  grad[1..P] ← 0
  for j = L down to 1:
      if gate j is parametrised (generator G_j, parameter index k):
          grad[k] += Im(⟨λ|G_j|ψ⟩)
      |ψ⟩ ← U_j†(θ)|ψ⟩
      |λ⟩ ← U_j†(θ)|λ⟩

  return grad
```

Parameter sharing is handled automatically: `grad[k]` accumulates from all gates sharing $\theta_k$.

### Complexity

| Resource | Cost |
|---|---|
| Time | $O(L \cdot 2^n)$ — one $L$-gate forward pass + one $2L$-gate backward pass |
| Memory | $2 \times 2^n$ complex amplitudes + $P$ gradient doubles |
| Gate applications | Exactly $3L$ |

### Pauli Inner Product Formulas

Let `mask = 1 << q`. Pairs $(x,\, x\oplus\texttt{mask})$ are processed with $x_q = 0$.

| Generator | Formula |
|---|---|
| $Z_q$ | $\sum_x \operatorname{Im}[\lambda^*[x]\psi[x] - \lambda^*[x\oplus m]\psi[x\oplus m]]$ |
| $X_q$ | $\sum_x \operatorname{Im}[\lambda^*[x]\psi[x\oplus m] + \lambda^*[x\oplus m]\psi[x]]$ |
| $Y_q$ | $\sum_x \operatorname{Re}[-\lambda^*[x]\psi[x\oplus m] + \lambda^*[x\oplus m]\psi[x]]$ |

The $Y$ formula uses $[Y_q|\psi\rangle]_x = -(-1)^{x_q} \cdot i \cdot \psi[x\oplus m]$, verified: $x_q=0 \Rightarrow -i\psi[x\oplus m]$ matches $\langle 0|Y|1\rangle = -i$; $x_q=1 \Rightarrow +i\psi[x\oplus m]$ matches $\langle 1|Y|0\rangle = +i$.

---

## Mixed-State (Density Matrix) Formulation

For $f = \operatorname{tr}(O\rho(\boldsymbol{\theta}))$ with $\rho(\boldsymbol{\theta}) = U_L\cdots U_1\,\rho_{\text{in}}\,U_1^\dagger\cdots U_L^\dagger$, define:

- $\rho^{(j)} = U_j\cdots U_1\,\rho_{\text{in}}\,U_1^\dagger\cdots U_j^\dagger$ (forward-propagated state after gate $j$)
- $O^{(0)} = U_1^\dagger\cdots U_L^\dagger\,O\,U_L\cdots U_1$, and $O^{(j)} = U_j O^{(j-1)} U_j^\dagger$ for $j \geq 1$, so that $O^{(j)} = U_{j+1}^\dagger\cdots U_L^\dagger\,O\,U_L\cdots U_{j+1}$

Then:

$$\frac{\partial f}{\partial \theta_k}\bigg|_{\text{gate }j} = \operatorname{Im}\!\left(\operatorname{tr}(O^{(j)}\,G_j\,\rho^{(j)})\right)$$

**Derivation.** Differentiating and using $dU_j/d\theta = -\frac{i}{2}G_j U_j$, setting $z = \operatorname{tr}(O^{(j)} G_j \rho^{(j)})$, and noting $z^* = \operatorname{tr}(O^{(j)}\rho^{(j)}G_j)$ (cyclic trace + Hermiticity), the result follows from $\frac{i}{2}(z^* - z) = \operatorname{Im}(z)$.

For state-vector simulation this reduces to the pure-state formula via $\operatorname{tr}(O^{(j)} G_j \rho^{(j)}) = \langle\lambda_j|G_j|\phi_j\rangle$.

### Indexing convention

The correct pairing at site $j$ is:

- $\rho^{(j)}$: state propagated **through** gates $1,\ldots,j$ (after $U_j$).
- $O^{(j)}$: observable back-propagated through gates $L,\ldots,j{+}1$ — i.e., **not** through $U_j$.

Gate $j$ appears in $\rho^{(j)}$ but not in $O^{(j)}$, so the commutator $[G_j,\cdot]$ is sandwiched at exactly the right position.

### Density-matrix adjoint algorithm

```
DM_ADJOINT_GRADIENT(gates[1..L], rho_in, O_diag):

  # Note: assumes O is diagonal in the computational basis.
  # For a general Hermitian O, replace diag(O_diag) with the full
  # matrix and dot(O_diag, diag(rho)) with tr(O @ rho).

  # --- Backward pass: build O^(0) ---
  # O^(0) = U_1† ... U_L† O U_L ... U_1
  # Start from O as a diagonal matrix; it becomes dense after first conjugation.
  O_curr = diag(O_diag)            # N×N, initially diagonal (N = 2^n)
  for j = L down to 1:
      O_curr = U_j† O_curr U_j    # conjugation: O <- U_j† O U_j

  # --- Forward pass: propagate rho and O simultaneously ---
  rho_curr = copy(rho_in)          # N×N density matrix
  grad[1..P] = 0

  for j = 1 to L:
      O_curr   = U_j  O_curr  U_j†   # O^(j) = U_j O^(j-1) U_j†
      rho_curr = U_j  rho_curr U_j†   # rho^(j) = U_j rho^(j-1) U_j†
      if gate j is parametrised (generator G_j, parameter index k):
          grad[k] += Im( tr( O_curr @ G_j @ rho_curr ) )

  return grad
```

**Live memory:** 2 dense $N\!\times\!N$ complex matrices (`O_curr`, `rho_curr`) = $32\cdot 4^n$ bytes. No intermediate states stored.

**Gate-application cost** (one conjugation $A \mapsto U A U^\dagger$ where $U$ acts on $k$ qubits): each of the $N^2$ entries requires $O(4^k)$ work via the sparse tensor-product structure → $O(N^2) = O(4^n)$ per conjugation for fixed $k$.

### Computational cost comparison

| Method | Gate-channel calls | Total arithmetic | Memory |
|---|---|---|---|
| Adjoint (density matrix) | $3L$ | $O(L \cdot 4^n)$ | $32\cdot 4^n$ bytes (2 matrices) |
| Parameter shift (density matrix) | $2L^2$ | $O(L^2 \cdot 4^n)$ | $16\cdot 4^n$ bytes (1 matrix) |

The adjoint method is $O(L)$ times cheaper in compute at the cost of one extra matrix. For 100-gate circuits this is a 100× speedup.

---

## HPC Enhancements

### 1. Fused Backward-Step Kernel

At each parametrised gate, fuse the inner product accumulation and both adjoint gate applications into one loop over pairs $(x, x\oplus\texttt{mask})$:

```c
for each pair (x, x_flip):
    phi_0 = ψ[x];  phi_1 = ψ[x_flip]
    lam_0 = λ[x];  lam_1 = λ[x_flip]
    grad  += inner_product_contribution(...)   // read old values
    ψ[x] = ...; ψ[x_flip] = ...               // write updated phi
    λ[x] = ...; λ[x_flip] = ...               // write updated lambda
```

No cross-pair dependencies — correctness guaranteed. Memory traffic: $4\text{R}+4\text{W}$ per pair (fused) vs. $8\text{R}+4\text{W}$ (unfused), a 33% reduction, with additional cache-locality gains for large state vectors.

For non-parametrised gates, fuse the dual $U_j^\dagger$ application into one pass similarly.

### 2. SIMD Vectorisation

Pair loops have no cross-pair dependencies. All three Pauli variants are trivially vectorisable (AVX2: 2 complex doubles/register; AVX-512: 4 complex doubles/register; NEON/SVE on ARM).

### 3. Multi-Observable Batching

For $K$ observables: one forward pass, $K$ backward passes sharing gate-adjoint applications. Batch $K$ adjoint states together per backward pass to amortise $|\psi\rangle$ reads. Arithmetic intensity saturates at $\approx 0.375$ flops/byte (memory-bandwidth-bound). Recommend $K = 2$–$8$; beyond that, L1 working set spills.

Optimal memory layout: interleave $K$ adjoint vectors so reading all $K$ values at one index is a contiguous SIMD load.

### 4. QAOA Layer-Level Commutativity

**Theorem.** If all generators $P_1,\ldots,P_M$ in a layer commute pairwise ($[P_a, P_b] = 0$), then all gradient contributions can be evaluated at the layer **boundary** states:

$$\langle\lambda_a|P_a|\phi_a\rangle = \langle\lambda_{\text{layer}}|P_a|\phi_{\text{layer}}\rangle \quad \forall a$$

**Proof.** Since $[P_a, U_{\text{layer}}] = 0$, the expression telescopes to the boundary. The layer then uncomputes as a single diagonal phase sweep.

---

## QAOA Mixer Layer

The mixer $U_{\text{mix}}(\beta) = \prod_i R_x^{(i)}(\beta)$ has all $n$ gates sharing $\beta$ and commuting pairwise. By commutativity:

$$\frac{\partial f}{\partial \beta} = \operatorname{Im}\!\left[\langle\lambda_{\text{layer}}\Big|\sum_i X_i\Big|\phi_{\text{layer}}\rangle\right]$$

### Hadamard Diagonalisation

$$\sum_i X_i = H^{\otimes n}\!\left(\sum_i Z_i\right)H^{\otimes n}, \quad \left(\sum_i Z_i\right)|x\rangle = (n - 2\operatorname{popcount}(x))\,|x\rangle$$

The full backward step for the mixer layer:

1. WHT both $|\phi\rangle$ and $|\lambda\rangle$
2. Diagonal inner product with weight $n - 2\operatorname{popcount}(x)$ → gradient
3. Diagonal phase: $\tilde\phi[x] \leftarrow e^{+i\beta(n - 2\operatorname{popcount}(x))/2}\tilde\phi[x]$, same for $|\tilde\lambda\rangle$
4. Inverse WHT both vectors

Similarly, the uncomputation $U_{\text{mix}}^\dagger = H^{\otimes n} R_z^{\otimes n}(+\beta)\, H^{\otimes n}$ is a diagonal operation in the Hadamard basis — steps 3–4 above. Since the WHT from step 1 is already computed, the forward states $|\tilde\phi\rangle$ and $|\tilde\lambda\rangle$ can be reused directly for uncomputation without an extra transform.

No asymptotic improvement ($O(n\cdot 2^n)$ either way), but WHT butterfly access patterns are sequential and SIMD-friendly, unlike the $n$ stride-$2^i$ passes of the naive approach.

---

## Scope and Limitations

| Case | Handling |
|---|---|
| Single-qubit Pauli rotations Rx/Ry/Rz | Adjoint method, Pauli inner product formulas |
| Multi-qubit Pauli gadgets $e^{-i\theta P_{\text{string}}/2}$ | Same adjoint identity; diagonal observables stored as arrays |
| Non-Pauli generators (e.g., phase gate) | Generator $G$ must be identified; inner product cost may increase |
| Controlled rotations $CR_z(\theta)$ | Generator $\frac{1}{2}(Z_t - Z_cZ_t)$; adjoint method applies |
| Density matrix / noisy circuits | Adjoint method applies exactly — see §Mixed-State (Density Matrix) Formulation above. PSR also applies and is used as reference. |
| $O$ non-Hermitian | Decompose into Hermitian + anti-Hermitian parts; run two backward passes |
| Unitary circuits only | Adjoint method requires all gates to be unitary |
| Noise-after-gate model | **TODO:** Implement the two-pass noisy adjoint algorithm (forward: store states through noisy channels; backward: propagate observable, apply adjoint noise at each layer). Theory and pseudocode in `~/Projects/thesis/4.VQA/vqa.tex` (Algorithm 3, eq:gradPUQCnoise). Complexity: $O(L \cdot d_{\max} \cdot N^3)$ time, $O(L \cdot N^2)$ space. |

---

## Parameter-Shift as Verification

Retain a parameter-shift implementation for correctness testing. For shared parameters, use per-gate shifts (2 evaluations per gate occurrence). Agreement with adjoint gradients should hold to 10+ significant digits; degradation below 8 digits indicates excessive circuit depth.

### PSR pseudocode — density-matrix simulator

```
DM_PSR_GRADIENT(gates[1..L], rho_in, O_diag, shift=π/2):

  grad[1..L] = 0

  for j = 1 to L:
      # (+) shifted circuit
      rho = copy(rho_in)
      for k = 1 to L:
          theta_k = gates[k].theta + (shift if k==j else 0)
          rho = apply_dm_channel(rho, gates[k].G, theta_k)
      f_plus = dot(O_diag, diag(rho))   # tr(O rho) = Σ_i O[i] rho[i,i]

      # (–) shifted circuit
      rho = copy(rho_in)
      for k = 1 to L:
          theta_k = gates[k].theta - (shift if k==j else 0)
          rho = apply_dm_channel(rho, gates[k].G, theta_k)
      f_minus = dot(O_diag, diag(rho))

      grad[j] = (f_plus - f_minus) / 2

  return grad


# Conjugation channel for G² = I (Pauli-string generators):
#   U = cos(θ/2)·I − i·sin(θ/2)·G,  rho_out = U rho U†
apply_dm_channel(rho, G, theta):
    c = cos(theta/2);  s = sin(theta/2)
    return c²·rho + s²·(G @ rho @ G) − i·c·s·(G @ rho − rho @ G)
```

Total circuit evaluations: $2L$, each of $L$ gate-channel calls → $2L^2$ gate-channel calls total.

---

## Ansatz Literature Survey — PSR Compatibility

Moved to `~/Projects/thesis/4.VQA/PSR_ANSATZ_SURVEY.md`. Surveys 33 papers across 8 ansatz families for PSR compatibility (10 YES, 19 PARTIAL, 4 NO).
