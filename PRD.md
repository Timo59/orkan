# Product Requirements Document: qSim

**Author:** [Your Name]
**Status:** Draft
**Last Updated:** 2025-12-24

## Motivation

Existing quantum frameworks (Qiskit, Cirq, PennyLane) simulate hardware abstractions—measurement
statistics, post-selection for linear combinations of unitaries (LCU), parameter shift rules for
gradients. These constraints make sense for hardware-focused development but introduce unnecessary
overhead for fundamental algorithm research.

qSim provides **direct mathematical access** to quantum state evolution:
- LCU implemented as complex-valued array operations (not post-selection)
- Gradients computed analytically from state derivatives (not parameter shift)
- Mixed states (density matrices) with efficient packed storage
- Direct state access for algorithm analysis and debugging

**Primary use case**: Variational algorithm research using LCU and mixed-state dynamics, where
mathematical operations on states matter more than hardware-realistic execution.

## Purpose

qSim is a quantum simulation library providing **exact analytical access** to quantum state
evolution. Users prepare states, apply circuits, and retrieve quantities (expectation values,
gradients, matrix elements, probability amplitudes) through direct computation rather than
statistical sampling. This enables rapid algorithm prototyping without hardware constraints like
repeated circuit execution or measurement limitations.

qSim provides a **complementary computational model** to quantum hardware for algorithm development.
Within simulator qubit limits (n≤32), it enables exact analytical access to quantum states—computing
gradients, expectation values, and matrix elements directly rather than through measurement
statistics. This eliminates **prototyping overhead** (block encoding for linear combinations,
repeated circuit execution for sampling) at the cost of classical exponential scaling.

Built on the qHiPSTER architecture with extensions for density matrix simulation, qSim exploits the
sparse, local structure of quantum operations to:
- **Avoid exponential memory overhead**: No Kronecker product construction (infeasible beyond n≈15)
- **Minimize error accumulation**: O(2^n) local operations vs O(2^3n) dense matrix multiplication
- **Enable exact computation**: Direct analytical access to amplitudes, gradients, and observables

External optimization algorithms integrate via clean C API.

## Target Users

**Primary**: Algorithm researchers developing variational quantum algorithms requiring direct state
manipulation, LCU operations, and mixed-state dynamics.

**Typical workflow**: Algorithm prototyping at scales (n≤32 qubits) where exact simulation enables
analytical validation before hardware deployment.

**Not intended for**: Hardware-realistic simulation, large-scale distributed computing, or
production quantum applications.

## Core Principles

1. **Exact analytical computation**: Direct access to state vectors and density matrices for
   computing observables, gradients, and matrix elements without statistical sampling
2. **Dual state representation**: Supports pure states (state vectors) and mixed states (density
   matrices with packed lower-triangular storage for memory efficiency)
3. **Direct sparse operations**: Quantum gates implemented as functions acting directly on state
   arrays, exploiting locality to avoid constructing full gate matrices
4. **Library not application**: Provides primitives (gates, observables, gradients, measurements)
   for building variational algorithms; optimization layer is external

## Scope

### In Scope
- Pure state representation (state vectors) and manipulation
- Mixed state representation (density matrices with packed lower-triangular storage)
- Fundamental single-qubit gates (Pauli-X/Y/Z, Hadamard, phase gates, rotations)
- Fundamental two-qubit gates (CNOT, CZ, SWAP, controlled rotations)
- Three-qubit gates (Toffoli)
- Exact computation of observable expectation values
- Gradient computation for parameterized gates
- Matrix element calculation for observables
- Measurement outcome sampling from state amplitudes
- C API with BLAS/LAPACK backend for performance
- Comprehensive unit tests for all modules
- API documentation and usage examples

### Out of Scope
- Optimization algorithms (external to library)
- Circuit optimization and compilation
- Hardware backend integration
- Noise models (planned for future after noise-free version is complete)
- Error mitigation techniques
- Quantum error correction
- Multi-node distributed simulation
- Visualization tools
- Algorithm-level parallelization (users implement externally)

## Scalability Targets

**Qubit capacity**:
- Pure states: n=32 qubits (64GB RAM required)
- Mixed states: n=16 qubits (32GB packed storage, 64-128GB workstation recommended)

**Performance envelope**:
- Typical use case: Overnight optimization campaigns (<48h runtime)
- Target: Algorithm prototyping and hyperparameter exploration on workstation hardware

**Design constraint**: RAM-limited rather than compute-limited. Packed storage and direct sparse
operations avoid constructing full gate matrices via Kronecker products, which becomes infeasible
beyond n≈15 qubits (2GB+ per matrix). This approach pushes simulator capacity to hardware memory
ceiling.

## Success Criteria

1. **Correctness**: All unit tests pass with numerical tolerance ≤1e-12
   - Pure states: Norm preservation, gate unitarity
   - Mixed states: Hermiticity, trace preservation
   - Observable computation matches analytical results

2. **Scalability**: Successfully simulate within memory constraints
   - Pure states: n=32 qubits (64GB RAM)
   - Mixed states: n=16 qubits (32GB packed storage)
   - Runtime: <48h for typical optimization campaigns
   - Test suite: <2 minutes for all gate tests

3. **Reliability**: Error handling
   - All functions return status codes (0=success, negative=error)
   - Memory allocation failures handled gracefully
   - No undefined behavior or crashes on valid inputs

4. **Documentation**:
   - API fully documented (function signatures, parameters, return codes)
   - Usage examples for core workflows (state preparation, circuit application, observable
     computation)
   - Technical specifications for each module

5. **Performance baseline**:
   - Functional correctness prioritized over optimization
   - Competitive with research-grade simulators for n≤20
   - (Formal benchmarks deferred to post-v1.0)

## Technical Architecture

### Implementation

Implementation in C with BLAS/LAPACK for numerical operations. State representations use
complex-valued arrays (`cplx_t`) with platform-specific optimizations (Apple Accelerate on macOS,
OpenBLAS on Linux). All operations leverage 64-bit integer support (ILP64) for large state spaces.

Built on qHiPSTER architecture (cite: Intel-QS) with extensions for density matrix simulation.

### Sparse Operations Design

Gates operate directly on state arrays by exploiting locality—a single-qubit gate affects only
2^(n-1) amplitude pairs regardless of position. This avoids constructing full 2^n × 2^n gate
matrices via Kronecker products, which becomes infeasible beyond n≈15 qubits (2GB+ per matrix).

**For n=20 qubits:**
- Sparse approach: ~1M operations, ~16MB state
- Dense approach: Build 1M×1M matrix = 1TB memory + 1 trillion operations (impossible)

**For n=32 qubits:**
- Sparse approach: ~4B operations, ~64GB state
- Dense approach: Requires >16 exabytes (physically impossible)

### Packed Density Matrix Storage

Mixed states use packed lower-triangular storage for Hermitian density matrices (50% memory
reduction). Gates operate directly on packed format via generalized qHiPSTER algorithms—no unpacking
required in production code.

**Memory savings critical at target scales:**
- n=16 mixed state: Full matrix 64GB → Packed 32GB
- Enables n=16 on 64GB workstations; full storage would require 128GB+

### Numerical Precision

Double-precision complex arithmetic (IEEE 754) via BLAS. Sparse gate operations with local structure
provide superior numerical stability compared to dense matrix methods—errors accumulate through
O(2^n) operations per gate rather than O(2^3n) for full matrix multiplication.

Test tolerance: ~1e-12 (BLAS precision with margin for accumulated errors). Formal error analysis
deferred pending empirical validation.

### Error Handling

Functions return integer status codes (0 = success, negative = error). Minimal input validation for
performance—callers responsible for providing valid parameters. No assertions or program termination
on errors.

Error code definitions maintained in technical specifications.

### Parallelization Characteristics

Gate operations are memory-bandwidth limited rather than compute-limited. Empirical tests show
gate-level parallelization (OpenMP) achieves ~3× speedup on 10 cores—limited by memory access
patterns in O(2^n) streaming operations.

**RAM constraints:**
- Large problems (n=32 pure, n=16 mixed) consume most/all available memory
- Multiple parallel instances infeasible at target qubit counts
- Smaller problems allow users to instantiate multiple state objects for parallel evaluation

**Library scope**: Provides thread-safe state object creation. Algorithm-level parallelization
(running multiple circuit evaluations concurrently) implemented by users external to library.

**Future consideration**: GPU acceleration may improve memory-bandwidth-limited operations.

---

Detailed technical specifications, architecture decisions, and API documentation are maintained in
separate technical documents.
