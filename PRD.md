# Product Requirements Document: qSim

**Author:** [Your Name]
**Status:** Draft
**Last Updated:** 2025-12-22

## Purpose

qSim is a quantum state simulator designed to evaluate variational quantum algorithms for combinatorial optimization research. It provides the computational feedback a noise-free quantum computer would deliver—state evolution under quantum gates and measurement outcomes—without implementing optimization routines. All optimization algorithms operate externally to the library, using qSim solely as an oracle for quantum state queries.

## Core Principles

1. **Evaluation-focused**: Simulator serves primarily for evaluating the effect of variational quantum algorithms
2. **Dual state representation**: Supports pure states (state vectors) and mixed states (density matrices; lower triangle storage for memory efficiency)
3. **Direct gate implementation**: Fundamental quantum gates implemented as functions acting directly on state arrays, exploiting the sparse structure of local quantum operations
4. **Exact computation**: Quantities like expectation values, gradients, and matrix elements computed directly and exactly from complex-valued arrays
5. **Measurement simulation**: Output state amplitudes define probability distributions for measurement sampling

## Scope

### In Scope
- Pure state representation (state vectors) and manipulation
- Mixed state representation (density matrices with lower-triangular storage)
- Fundamental single-qubit gates (Pauli-X/Y/Z, Hadamard, phase gates, rotations)
- Fundamental two-qubit gates (CNOT, CZ, SWAP, controlled rotations)
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

## Success Criteria

1. **Correctness**: All unit tests pass for each module (state representation, gates, observables, measurements)
2. **Documentation**: API fully documented with usage examples
3. **Performance**: Efficient memory usage via lower-triangular density matrix storage and BLAS-optimized operations
4. **Usability**: Clean C API suitable for integration with external optimization routines

## Technical Architecture

Implementation in C with BLAS/LAPACK for numerical operations. State representations use complex-valued arrays (`cplx_t`) with platform-specific optimizations (Apple Accelerate on macOS, OpenBLAS on Linux). All operations leverage 64-bit integer support (ILP64) for large state spaces.

Detailed technical specifications, architecture decisions, and API documentation are maintained in separate technical documents.
