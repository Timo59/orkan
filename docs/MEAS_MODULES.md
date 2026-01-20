# Measurement Module: Gradient Computation for Variational Circuits

**Module:** Gradient Computation & Cost Function Evaluation
**Header:** `include/meas.h` (planned)
**Implementation:** `src/meas.c` (planned)
**Last Updated:** 2026-01-08

---

## Overview

This module addresses gradient computation for variational quantum algorithms (VQE, QAOA, quantum
machine learning). The primary challenge: given a parameterized circuit U(θ) and cost function
C(θ) = Tr(H ρ(θ)), compute the gradient ∇C efficiently.

**Key design considerations:**
- **Memory efficiency**: Preserve mHiPSTER's packed hermitian format (N(N+1)/2 elements)
- **Parallelization**: Gradient computations are embarrassingly parallel
- **Accuracy**: Exact gradients (not finite differences) for reliable optimization
- **Shared parameters**: Handle circuits where multiple gates depend on the same parameter

---

## Background: The Optimization Problem

### Variational Quantum Algorithms

Variational algorithms optimize classical parameters θ = (θ₁, ..., θ_p) controlling a quantum circuit:

```
1. Prepare initial state ρ₀ (often |0...0⟩⟨0...0|)
2. Apply parameterized circuit: ρ(θ) = U(θ) ρ₀ U(θ)†
3. Measure cost function: C(θ) = Tr(H ρ(θ))
4. Update parameters: θ ← θ - η ∇C(θ)
5. Repeat until convergence
```

**Examples:**
- **VQE**: H = molecular Hamiltonian, minimize ground state energy
- **QAOA**: H = problem Hamiltonian (MaxCut, etc.), minimize/maximize objective
- **QML**: H encodes classification loss or generative model likelihood

### The Gradient Challenge

For p parameters, we need ∇C = (∂C/∂θ₁, ..., ∂C/∂θ_p).

| Method | Evaluations | Accuracy | Memory |
|--------|-------------|----------|--------|
| Finite differences | 2p | O(ε) truncation + O(1/ε) rounding | O(N(N+1)/2) |
| Parameter-shift | 2p (or 2m per shared param) | Exact | O(N(N+1)/2) |
| Backpropagation | 1 forward + 1 backward | Exact | O(d · N(N+1)/2) |

---

## Automatic Differentiation (Autodiff)

Autodiff computes **exact derivatives** of functions expressed as computer programs by systematically
applying the chain rule to elementary operations.

### Computational Graphs

Any computation can be represented as a directed acyclic graph (DAG):
- **Nodes**: Variables (inputs, intermediates, outputs)
- **Edges**: Elementary operations (+, -, ×, /, sin, exp, etc.)

**Example:** f(x₁, x₂) = x₁x₂ + sin(x₁)

```
Computation sequence:
    v₁ = x₁               (input)
    v₂ = x₂               (input)
    v₃ = v₁ × v₂          (multiplication)
    v₄ = sin(v₁)          (sine)
    v₅ = v₃ + v₄          (addition) → output f

Graph structure:
    x₁ ──┬──→ [×] ──→ v₃ ──┐
         │                  ├──→ [+] ──→ f
         └──→ [sin] ──→ v₄ ─┘
    x₂ ─────→ [×]
```

### Forward Mode (Tangent Mode)

Propagates derivatives **forward** alongside the primal computation.

For each variable vᵢ, simultaneously compute the **tangent** v̇ᵢ = ∂vᵢ/∂x_seed:

```
Primal                    Tangent (seed: ∂/∂x₁)
──────                    ──────────────────────
v₁ = x₁                   v̇₁ = 1                    (∂x₁/∂x₁ = 1)
v₂ = x₂                   v̇₂ = 0                    (∂x₂/∂x₁ = 0)
v₃ = v₁ × v₂              v̇₃ = v̇₁·v₂ + v₁·v̇₂ = x₂   (product rule)
v₄ = sin(v₁)              v̇₄ = cos(v₁)·v̇₁ = cos(x₁) (chain rule)
v₅ = v₃ + v₄              v̇₅ = v̇₃ + v̇₄ = x₂ + cos(x₁)
```

**Result:** ∂f/∂x₁ = x₂ + cos(x₁)

**Cost:** One forward pass computes derivative w.r.t. ONE input. For p inputs: p passes.

**Use case:** Few inputs, many outputs (e.g., sensitivity analysis, ODEs).

### Reverse Mode (Adjoint Mode) = Backpropagation

Propagates derivatives **backward** from output to inputs.

**Phase 1: Forward pass** — compute and store all intermediate values.

**Phase 2: Backward pass** — for each variable, compute the **adjoint** v̄ᵢ = ∂f/∂vᵢ.

Starting from output (v̄_output = 1), propagate backwards using chain rule:

```
Backward propagation:
──────────────────────
v̄₅ = 1                                    (∂f/∂f = 1, seed)
v̄₄ = v̄₅ · ∂v₅/∂v₄ = 1 · 1 = 1            (v₅ = v₃ + v₄)
v̄₃ = v̄₅ · ∂v₅/∂v₃ = 1 · 1 = 1            (v₅ = v₃ + v₄)
v̄₂ = v̄₃ · ∂v₃/∂v₂ = 1 · v₁ = x₁          (v₃ = v₁ × v₂)
v̄₁ = v̄₃ · ∂v₃/∂v₁ + v̄₄ · ∂v₄/∂v₁        (v₁ feeds both v₃ and v₄)
   = 1 · v₂ + 1 · cos(v₁)
   = x₂ + cos(x₁)
```

**Result:** ∂f/∂x₁ = x₂ + cos(x₁), ∂f/∂x₂ = x₁ — BOTH gradients in ONE backward pass.

**Cost:** One forward + one backward pass gives derivatives w.r.t. ALL inputs.

**Use case:** Many inputs, few outputs (optimization, machine learning, variational algorithms).

---

## Backpropagation: Detailed Algorithm

### The Chain Rule for Computational Graphs

For a node v with children w₁, ..., w_k (nodes that use v as input):

```
v̄ = Σⱼ w̄ⱼ · ∂wⱼ/∂v
```

Each node **accumulates** gradient contributions from ALL paths to the output.

### Memory Requirements

Reverse mode must **store all intermediate values** from the forward pass because local gradients
∂wⱼ/∂v often depend on the primal values.

**Example:** For v₃ = v₁ × v₂, computing ∂v₃/∂v₁ = v₂ requires knowing v₂.

For quantum circuit simulation:
- **Forward pass**: O(d · 2ⁿ) operations for d gates on n qubits
- **Storage**: O(d · N(N+1)/2) if storing all intermediate density matrices
- **Alternative**: O(N(N+1)/2) with recomputation (checkpoint strategies)

### Pseudocode: Reverse-Mode Autodiff

```
FORWARD_PASS(inputs, graph):
    tape = []                           // Record operations
    for node in topological_order(graph):
        node.value = compute(node.op, node.inputs)
        tape.append((node, node.inputs, node.value))
    return tape

BACKWARD_PASS(tape, output_node):
    adjoints = {output_node: 1.0}       // Seed with ∂f/∂f = 1
    for (node, inputs, value) in reverse(tape):
        if node not in adjoints:
            continue
        v_bar = adjoints[node]
        for (input_node, local_grad) in compute_local_gradients(node, inputs, value):
            if input_node not in adjoints:
                adjoints[input_node] = 0.0
            adjoints[input_node] += v_bar * local_grad    // Accumulate!
    return adjoints
```

### Key Insight: Gradient Accumulation

When a variable is used multiple times (shared parameters, reused computations), the backward
pass **automatically sums** contributions from all uses:

```
θ used in gates U₁(θ), U₂(θ), U₃(θ):

Forward:  θ → U₁ → ... → U₂ → ... → U₃ → C

Backward: θ̄ = ∂C/∂U₁ · ∂U₁/∂θ + ∂C/∂U₂ · ∂U₂/∂θ + ∂C/∂U₃ · ∂U₃/∂θ
              \_____________/   \_____________/   \_____________/
                 path 1            path 2            path 3
```

This is precisely what neural networks with weight sharing (e.g., CNNs, RNNs) rely on.

---

## Application to Quantum Circuit Simulation

### State Vector Formulation

For |ψ(θ)⟩ = U_d(θ_d) ··· U₁(θ₁) |0⟩ and C(θ) = ⟨ψ(θ)|H|ψ(θ)⟩:

```
Forward pass:
    |ψ₀⟩ = |0...0⟩
    |ψ₁⟩ = U₁(θ₁)|ψ₀⟩     → store |ψ₀⟩
    |ψ₂⟩ = U₂(θ₂)|ψ₁⟩     → store |ψ₁⟩
    ...
    |ψ_d⟩ = U_d(θ_d)|ψ_{d-1}⟩
    C = ⟨ψ_d|H|ψ_d⟩

Backward pass:
    Compute ∂C/∂θ_d, ∂C/∂θ_{d-1}, ..., ∂C/∂θ₁
    using stored intermediate states
```

### Density Matrix / Hermitian Operator Formulation (mHiPSTER)

For ρ(θ) = U(θ) ρ₀ U(θ)† and C(θ) = Tr(H ρ(θ)):

```
Forward pass:
    ρ₀ = |0...0⟩⟨0...0|
    ρ₁ = U₁(θ₁) ρ₀ U₁(θ₁)†     → store ρ₀
    ρ₂ = U₂(θ₂) ρ₁ U₂(θ₂)†     → store ρ₁
    ...
    ρ_d = U_d(θ_d) ρ_{d-1} U_d(θ_d)†
    C = Tr(H ρ_d)

Backward pass:
    Need gradients through conjugation: ∂(UρU†)/∂θ
```

**Challenge for mHiPSTER:** The packed column-major format requires custom gradient rules for the
conjugation operation. Standard autodiff frameworks (TensorFlow, PyTorch, JAX) expect dense arrays.

---

## Parameter-Shift Rule

For gates of the form U(θ) = exp(-iθG/2) where generator G has eigenvalues ±1:

```
∂C/∂θ = ½[C(θ + π/2) - C(θ - π/2)]
```

**Proof sketch:** The expectation value C(θ) = ⟨ψ|U†(θ)HU(θ)|ψ⟩ can be shown to have sinusoidal
dependence on θ with period 2π. The parameter-shift formula follows from the derivative of
sin/cos at phase shifts of π/2.

### Advantages

| Property | Benefit |
|----------|---------|
| Exact | No truncation or rounding errors |
| Memory-efficient | Only need to evaluate C(θ), not store intermediates |
| Hardware-compatible | Works on real quantum devices (no backprop possible) |
| Parallelizable | All 2p shifted evaluations are independent |

### Gate Compatibility

| Gate Type | Generator | Eigenvalues | Shift |
|-----------|-----------|-------------|-------|
| RX(θ), RY(θ), RZ(θ) | X/2, Y/2, Z/2 | ±1/2 | π/2 |
| CRX, CRY, CRZ | Controlled Pauli | ±1/2, 0 | π/2 (with correction) |
| General rotation | Arbitrary | Varies | Generalized formula |

### Generalized Parameter-Shift

For generators with eigenvalues {λ₁, ..., λ_r}, the gradient requires more shifted evaluations.
For two distinct eigenvalues ±λ:

```
∂C/∂θ = λ · [C(θ + π/(4λ)) - C(θ - π/(4λ))] / sin(π/2)
```

---

## Handling Shared Parameters

When parameter θ appears in m different gates U₁(θ), ..., U_m(θ):

### Method 1: Individual Shifts (Exact)

Apply parameter-shift to each gate individually, sum the results:

```
∂C/∂θ = Σᵢ₌₁ᵐ ½[C(θ⁺ⁱ) - C(θ⁻ⁱ)]
```

Where θ⁺ⁱ means: shift gate i by +π/2, keep all other gates at θ.

**Cost:** 2m circuit evaluations per shared parameter.

**Example:** QAOA with 10 edges sharing γ requires 20 evaluations for ∂C/∂γ.

### Method 2: Stochastic Sampling

For large m, randomly sample which gate to shift:

```
∂C/∂θ ≈ m · ½[C(θ⁺ʲ) - C(θ⁻ʲ)]    where j ~ Uniform(1, ..., m)
```

**Properties:**
- Unbiased estimator of true gradient
- Variance scales as O(m²/k) for k samples
- Cost: 2k evaluations (you choose k)

### Method 3: Combined Generator (Special Cases)

If all gates U₁(θ), ..., U_m(θ) can be written as exp(-iθG_total) where G_total = Σᵢ Gᵢ and all
Gᵢ commute with compatible eigenvalue structure, a single shift suffices.

**Example:** QAOA cost layer exp(-iγ Σ_edges Z_i Z_j) implemented as single rotation.

### Method 4: Autodiff (Automatic Accumulation)

Reverse-mode autodiff handles shared parameters automatically—the backward pass accumulates
gradients from all uses without explicit summation.

**Cost:** 1 forward + 1 backward (independent of m).

**Trade-off:** Requires O(d · N(N+1)/2) memory for intermediate storage.

---

## Comparison for mHiPSTER

| Method | Circuit Evals | Memory | Implementation | When to Use |
|--------|---------------|--------|----------------|-------------|
| Parameter-shift (full) | 2Σmᵢ | O(N(N+1)/2) | Simple | Few shared params |
| Stochastic shift | 2k (tunable) | O(N(N+1)/2) | Simple | Many shared params |
| Combined generator | 2p | O(N(N+1)/2) | Gate restructuring | QAOA, commuting gates |
| Autodiff (reverse) | 2 (fwd+bwd) | O(d·N(N+1)/2) | Custom gradients | Memory available |

### Recommended Approach for mHiPSTER

Given mHiPSTER's focus on memory efficiency:

**Primary:** Parameter-shift rule with parallelization
- Preserves packed hermitian format
- Embarrassingly parallel (ideal for pool-based execution)
- No framework dependencies

**Hybrid strategy:**
1. Group parameters by sharing pattern
2. Small groups (m ≤ 5): Full parameter-shift
3. Large groups (m > 5): Stochastic sampling or combined generator
4. Parallelize all independent evaluations

---

## Implementation Design

### Data Structures

```c
typedef struct {
    int* gate_indices;      // Which gates use this parameter
    int num_gates;          // m = number of gates sharing this parameter
    double shift;           // Shift amount (π/2 for standard gates)
} param_group_t;

typedef struct {
    param_group_t* groups;  // Array of parameter groups
    int num_params;         // p = total number of parameters
    int num_groups;         // Number of distinct parameters
} gradient_spec_t;
```

### Core Functions (Planned)

| Function | Description | Status |
|----------|-------------|--------|
| `cost_eval()` | Evaluate C(θ) = Tr(H ρ(θ)) | 🔲 TODO |
| `gradient_parameter_shift()` | Full parameter-shift gradient | 🔲 TODO |
| `gradient_stochastic()` | Stochastic parameter-shift | 🔲 TODO |
| `gradient_parallel()` | Parallel gradient computation | 🔲 TODO |

### API Design

```c
// Evaluate cost function
double cost_eval(
    const state_t* rho_init,      // Initial state
    const circuit_t* circuit,     // Parameterized circuit
    const double* params,         // Parameter values
    const hermitian_t* H          // Observable/Hamiltonian
);

// Compute gradient via parameter-shift
void gradient_parameter_shift(
    const state_t* rho_init,
    const circuit_t* circuit,
    const double* params,
    const hermitian_t* H,
    const gradient_spec_t* spec,
    double* grad_out              // Output: p-dimensional gradient
);

// Parallel gradient computation
void gradient_parallel(
    const state_t* rho_init,
    const circuit_t* circuit,
    const double* params,
    const hermitian_t* H,
    const gradient_spec_t* spec,
    int num_threads,
    double* grad_out
);
```

---

## Performance Considerations

### Parallelization Strategy

Parameter-shift evaluations are independent → ideal for parallel execution.

```
For p parameters with sharing pattern (m₁, m₂, ..., m_p):
    Total evaluations = 2 Σᵢ mᵢ

With T threads:
    Runtime ≈ (2 Σᵢ mᵢ) / T × time_per_circuit
```

### Memory Access Patterns

Each circuit evaluation uses mHiPSTER's O(N² · 2^k) conjugation algorithm:
- Memory-bandwidth bound for large N
- Cache-friendly for packed column-major format
- Independent evaluations don't share state (no false sharing)

### Precision Considerations

| Precision | Memory | Gradient Accuracy | Use Case |
|-----------|--------|-------------------|----------|
| double (complex128) | 16 bytes/element | ~15 digits | Default |
| float (complex64) | 8 bytes/element | ~7 digits | Memory-constrained |

Parameter-shift provides exact gradients regardless of precision (no finite-difference errors).

---

## Integration with Optimization

### Gradient Descent Variants

```c
// Standard gradient descent
for (int iter = 0; iter < max_iter; iter++) {
    gradient_parallel(..., params, ..., grad);
    for (int i = 0; i < num_params; i++) {
        params[i] -= learning_rate * grad[i];
    }
}

// Adam optimizer (recommended for VQE/QAOA)
// Requires first and second moment estimates
```

### Convergence Monitoring

| Metric | Formula | Threshold |
|--------|---------|-----------|
| Gradient norm | ‖∇C‖ | < 1e-6 |
| Cost change | \|C_{t} - C_{t-1}\| | < 1e-8 |
| Parameter change | ‖θ_t - θ_{t-1}‖ | < 1e-6 |

---

## References

1. **Mitarai et al. (2018)**: "Quantum circuit learning" — Original parameter-shift rule
2. **Schuld et al. (2019)**: "Evaluating analytic gradients on quantum hardware"
3. **Mari et al. (2021)**: "Estimating the gradient and higher-order derivatives on quantum hardware"
4. **Bergholm et al. (2018)**: PennyLane — Autodiff framework for quantum computing
5. **Jones et al. (2023)**: QuEST distributed simulation — HPC considerations applicable to gradient computation

---

## Future Extensions

- **Hessian computation**: Second-order optimization (Newton methods)
- **Natural gradient**: Fisher information metric for better convergence
- **Noise-aware gradients**: Gradient estimation under realistic noise models
- **Checkpointing**: Trade compute for memory in autodiff approach
