# Product Requirements Document: qSim v1.0

**Document Version:** 1.2
**Last Updated:** 2025-12-22
**Status:** Draft
**Owner:** Timo Ziegler

---

## Executive Summary

qSim is a high-performance quantum computing simulation library written in C, designed to accelerate quantum algorithm 
conception by enabling quantum information theorists and algorithm designers to work at the ideal operation 
level—deferring hardware implementation complexity until algorithmic benefits are proven.

**The Problem:** Existing simulators (Qiskit, Cirq, etc.) mimic quantum hardware behavior, forcing algorithm designers
to prematurely address implementation overhead. This creates a conception bottleneck: Obtaining N samples requires N
circuit executions, performing non-physical operations, despite being implementable on real devices, like Linear 
Combination of Unitaries (LCU), and computing advanced quantities like the gradient of a parametrised quantum circuit
require awkward workarounds.

**The Solution:** qSim provides simulator-native capabilities that exploit linear algebra uniquely available in
simulation: (1) dual state representations (pure/mixed) with optimized storage, (2) direct non-physical operations
(3) exact computation of observable quantities, and (4) one-time execution with unlimited sampling.

**Target:** Workstation-scale systems (up to 30 qubits pure states, 15 qubits mixed states) with cross-platform support
(macOS, Linux, x86, ARM).

---

## 1. Product Vision & Differentiation

### 1.1 Key Differentiator

qSim explicitly distinguishes between capabilities that are:
1. **Simulator-native** (direct, efficient): Fundamental operations act on the quantum state representation, i.e., a 
complex-valued array. As the native gate set of most quantum computers mainly comprises local operations, these can be
implemented using a form of sparse matrix multiplication. Composite operations like sums or products can be established
using basic linear algebra. Measurements of a quantum state are samples drawn from a probability distribution given
by the absolute squares of the state's amplitudes in the measurement basis.
2. **Hardware-implementable with overhead**: All simulator-native operations *can* be implemented on real quantum 
computers, but require additional resources (ancilla qubits, repeated/mid-circuit measurements, circuit decomposition).
N measurements of the output state of a quantum circuits requires N circuit executions.


This enables a two-phase workflow:
- **Phase 1 (qSim)**: Validate algorithm at ideal operation level, prove value
- **Phase 2 (Hardware)**: Add implementation overhead, test robustness

**Impact**: Algorithm designers discard non-promising quantum algorithms 3-5x faster.

---

## 2. Goals & Success Metrics

### 2.1 Primary Goals

1. Enable quantum algorithm designers to work at the ideal operation level without premature hardware constraints
2. Provide simulator-native capabilities for hardware-implementable operations which are not hardware-native
3. Achieve competitive performance (workstation-scale, up to 30 qubits)
4. Ensure cross-platform portability (macOS, Linux, x86, ARM)

### 2.2 Success Metrics

| Metric | Target |
|--------|--------|
| Sampling efficiency | 1000× faster than circuit re-execution |
| Memory efficiency | ≤ 2× theoretical minimum |
| Execution speed | Competitive with Qiskit Aer (within 3×) |
| Test coverage | ≥ 90% |
| Cross-platform support | Pass test suite on macOS/Linux, x86/ARM |

---

## 3. Target Audience

### 3.1 Primary Audience

**Quantum Algorithm Researchers**
- PhD/Postdoc researchers developing and optimizing quantum algorithms (VQE, QAOA, novel approaches)
- Graduate students in quantum computing research
- **Key Needs**: Rapid prototyping without hardware overhead, direct access to key quantities and operations, fast
iteration cycles
- **Pain Points**: Existing simulators force premature hardware constraints and inefficient sampling workflow

### 3.2 Secondary Audience

**Research Software Engineers**
- Software developers integrating quantum algorithms into research workflows
- **Key Needs**: Reliable C API, good documentation, cross-platform compatibility
- **Pain Points**: Language barriers, dependency management

---

## 4. User Stories

**As a** quantum algorithm theorist,
**I want** to apply linear combinations of unitaries directly to quantum states,
**So that** I can prototype LCU-based algorithms naturally without manual decomposition into physically implementable 
gates.

**As a** variational algorithm researcher,
**I want** exact gradients of parametrized circuits computed automatically,
**So that** I can optimize ansätze efficiently without parameter-shift rule overhead or numerical differentiation.

**As a** quantum information scientist,
**I want** to compute transition matrix elements ⟨ψᵢ|O|ψⱼ⟩ directly,
**So that** I can analyze algorithm dynamics without designing complex ancilla-based measurement schemes.

**As an** algorithm designer,
**I want** to execute circuits once and sample indefinitely from the result,
**So that** I can efficiently analyze measurement statistics and shot-noise effects without repeated circuit execution.

**As a** quantum algorithms engineer,
**I want** direct access to state amplitudes and phases,
**So that** I can debug algorithms by inspecting quantum interference and identifying sources of error.

**As a** research scientist,
**I want** both exact (shot-noise free) and sampled measurement modes,
**So that** I can validate ideal algorithm performance before analyzing realistic finite-sampling effects.

**As a** research software engineer,
**I want** a well-documented C API with clear examples,
**So that** I can integrate quantum simulations into research pipelines without language binding complexity.

---

## 5. Key Feature Requirements

### 5.1 Must-Have Features (P0): Core Algorithm Design Enablers

- **REQ-001**: Support pure state (statevector) and mixed state (density matrix) representation
- **REQ-002**: Implement native set of quantum gates (single- and two-qubit) efficiently
- **REQ-003**: Apply compositions directly, e.g., products and linear combinations
- **REQ-004**: Compute quantities subject to measurements exactly, e.g., expectation values and gradients of PQCs
- **REQ-005**: Sample from output quantum state efficiently by accessing the amplitude's squares

### 5.2 Should-Have Features (P1)
- **REQ-080**: OpenMP parallelization for large circuits
- **REQ-081**: Optimized BLAS/LAPACK integration
- **REQ-082**: Comprehensive API documentation with examples
- **REQ-083**: Input validation and error handling

### 5.3 Nice-to-Have Features (P2)
- **REQ-090**: Support for arbitrary unitary matrices
- **REQ-091**: Additional gates (iSWAP, CZ, CY, etc.)
- **REQ-092**: State vector inspection and debugging utilities

---

## 6. Non-Functional Requirements

### 6.1 Performance Requirements

- **NFR-001**: Simulate 30-qubit circuit with 100 gates in < 60 seconds (on reference hardware)
- **NFR-002**: Memory usage ≤ 2× theoretical minimum (2^n complex numbers for n qubits)
- **NFR-003**: Support parallel execution via OpenMP with near-linear scaling (up to 8 cores)

**Reference Hardware**: x86_64 workstation, 64GB RAM, 8-core CPU @ 3.0GHz

### 6.2 Scalability Requirements

- **NFR-010**: Support up to 30 qubits for pure states
- **NFR-011**: Support up to 15 qubits for mixed states
- **NFR-012**: Graceful degradation or clear error for out-of-memory conditions

### 6.3 Portability Requirements

- **NFR-020**: Compile and run on Linux (kernel 4.x+)
- **NFR-021**: Compile and run on macOS (11.0+)
- **NFR-022**: Support x86_64 and ARM64 architectures
- **NFR-023**: No platform-specific code outside cmake configuration

### 6.4 Reliability Requirements

- **NFR-030**: Test coverage ≥ 90% (line coverage)
- **NFR-031**: Zero memory leaks under Valgrind/AddressSanitizer
- **NFR-032**: All numerical results accurate to machine precision (≈10^-15 for double)
- **NFR-033**: Deterministic behavior (same input → same output)

### 6.5 Quality & Build Requirements

- **NFR-040**: Modular architecture with comprehensive test suite (unit and integration tests)
- **NFR-041**: CMake-based build system (≥ 3.27) with automated dependency resolution
- **NFR-042**: All public APIs documented

## 7. Technology Stack

- **Language**: C17
- **Build System**: CMake 3.27+
- **Linear Algebra**: CBLAS/LAPACK with ILP64 support (Accelerate on macOS, OpenBLAS on Linux)
- **Parallelization**: OpenMP (optional)
- **Testing**: Unity framework

**Note**: Detailed technical specifications, compiler requirements, and architectural constraints are documented in `docs/TECHNICAL_SPEC.md`.

---

## 8. Out of Scope (v1.0)

**Deferred to Future Versions:**
- Noise models, GPU acceleration, distributed simulation, circuit optimization
- Language bindings (Python, Julia, Rust)
- QASM import/export
- Advanced simulation backends (tensor networks, stabilizer, MPS)

---

## 9. Key Assumptions

- Users have basic understanding of quantum computing concepts
- Target platforms have sufficient memory for intended qubit counts (64GB+ RAM for 30 qubits)
- Developers can use C APIs and understand manual memory management

---

## 10. Key Risks

| Risk | Mitigation |
|------|------------|
| BLAS library incompatibility across platforms | Bundle OpenBLAS, provide verification tools |
| Performance doesn't meet competitive benchmarks | Early performance testing, profiling, optimization |
| API complexity hinders adoption | User testing, comprehensive documentation, examples |

---

## 11. Release Criteria (v1.0)

Version 1.0 will be released when:
- All P0 features implemented and tested
- Test suite passes on macOS/Linux, x86/ARM
- Test coverage ≥ 90%, zero memory leaks
- Performance benchmarks meet targets
- API reference and user guide complete
- 30-qubit Grover and VQE algorithms run successfully

---

## 12. Future Roadmap (Post-v1.0)

### v1.1 (Planned)
- QASM 2.0 import/export
- Enhanced error handling and debugging features
- Noise model foundation

### v1.2 (Planned)
- Basic circuit optimization passes
- Stabilizer simulator backend
- Extended gate library

### v2.0 (Vision)
- GPU acceleration (CUDA/OpenCL)
- Distributed simulation support
- Advanced noise models

---

## 13. Appendix

### 13.1 Glossary

**Core Quantum Concepts:**
- **Statevector**: Pure quantum state |ψ⟩ = ∑ᵢ αᵢ|i⟩ (2^n complex amplitudes for n qubits)
- **Density Matrix**: Mixed quantum state ρ (Hermitian matrix, stored as lower-triangle for efficiency)
- **Expectation Value**: ⟨O⟩ = ⟨ψ|O|ψ⟩ (observable measurement average)

**Algorithm-Specific:**
- **LCU (Linear Combination of Unitaries)**: Operation ∑ᵢ αᵢUᵢ, natural in algorithms but requiring decomposition on hardware
- **Gradient**: Derivative ∂⟨H⟩/∂θ for circuit parameter optimization
- **Transition Matrix Element**: Quantum amplitude ⟨ψᵢ|O|ψⱼ⟩ between states
- **VQE**: Variational Quantum Eigensolver (hybrid algorithm for ground state finding)
- **QAOA**: Quantum Approximate Optimization Algorithm

**Simulator Terminology:**
- **Simulator-Native**: Operations efficient in simulation (gradients, matrix elements, unlimited sampling) but requiring overhead on hardware
- **Shot Noise**: Statistical uncertainty from finite measurement samples
- **ILP64**: 64-bit BLAS interface for large state vectors

### 13.2 References

- qSim Repository: https://github.com/[org]/qsim (TBD)
- Technical Specification: See `docs/TECHNICAL_SPEC.md`
- API Reference: See `docs/API.md`
- Architecture Overview: See `docs/ARCHITECTURE.md`

### 13.3 Document History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2025-12-22 | [Author] | Initial PRD creation |
| 1.1 | 2025-12-22 | [Author] | Major repositioning: Focus on algorithm design enablers |
| 1.2 | 2025-12-22 | [Author] | Consolidated duplicate content, removed technical details better suited for TECHNICAL_SPEC.md, simplified personas and requirements |

---

**Document Classification**: Internal
**Review Cycle**: Quarterly or as needed for major changes
