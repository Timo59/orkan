# qSim

A quantum computer simulator library written in C for representing, manipulating, and working with quantum states.

## Description

qSim is a quantum computing simulation library focused on efficient quantum state representation and manipulation using industry-standard numerical libraries. The library supports both pure states (represented as statevectors) and mixed states (represented as density matrices).

## Features

- **Quantum State Representation**
  - Pure states as statevectors
  - Mixed states as density matrices (lower triangle representation)

- **State Manipulation Operations**
  - Initialize quantum states
  - Create uniform superposition (plus state)
  - Deep copy states
  - Efficient memory management

- **High-Performance Computing**
  - CBLAS (Basic Linear Algebra Subprograms) integration for optimized operations
  - OpenMP support for parallel computation
  - Cross-platform support (macOS and Linux)

- **Platform-Specific Optimizations**
  - macOS: Apple Accelerate framework (vecLib/LAPACK)
  - Linux: OpenBLAS for BLAS operations

## Requirements

- CMake 3.27 or higher
- C17-compatible compiler
- OpenMP (optional, for parallel computation)
- BLAS library with **64-bit integer support (ILP64)**:
  - macOS: Accelerate framework (built-in, ILP64 native)
  - Linux: **Bundled OpenBLAS** (default, built with ILP64 automatically)

### BLAS Configuration

qSim requires BLAS/LAPACK libraries with **64-bit integer (ILP64)** support. By default:

- **macOS**: Uses Apple Accelerate framework with native ILP64 support
- **Linux**: Automatically downloads and builds OpenBLAS 0.3.27 with ILP64 support

**Note**: The first build on Linux takes 10-30 minutes to compile OpenBLAS, but subsequent builds are fast as it's cached.

## Building

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

### Build Options

- `CMAKE_BUILD_TYPE`: Set to `Debug` or `Release`
- `USE_SYSTEM_OPENBLAS`: Use system OpenBLAS instead of bundled version (default: `OFF`)
- `ENABLE_ASAN`: Enable AddressSanitizer for memory debugging (default: `OFF`)

### Using System OpenBLAS (Advanced)

If you have OpenBLAS installed with ILP64 support, you can use it instead of the bundled version:

```bash
cmake -DUSE_SYSTEM_OPENBLAS=ON ..
cmake --build .
```

**Warning**: Most system OpenBLAS installations use 32-bit integers (LP64) and are **not compatible**. Only enable this if you're certain your OpenBLAS was built with `INTERFACE64=1`.

To verify your system BLAS has ILP64 support:

```bash
# Compile the verification tool
cd tools
gcc -o verify_ilp64 verify_ilp64.c -lopenblas -DOPENBLAS_USE64BITINT -I/usr/include

# Run verification
./verify_ilp64
```

Expected output: `✓ ILP64 support detected (sizeof(blasint) = 8)`

If verification fails, use the bundled OpenBLAS (default behavior).

## Project Structure

```
qsim/
├── CMakeLists.txt          # Main build configuration
├── README.md               # Project documentation
├── src/                    # Source implementation
│   ├── state.c            # State manipulation implementation
│   └── CMakeLists.txt     # Source build rules
├── include/               # Public headers
│   ├── q_types.h          # Type definitions
│   ├── state.h            # State structures and functions
│   └── qlib.h             # Main library header
├── docs/                  # Technical documentation
│   ├── STATE_MODULE.md    # State representation specification
│   ├── GATES_MODULE.md    # Gate operations specification
│   └── MEAS_MODULES.md    # Gradient computation for variational circuits
├── test/                  # Test suite
│   ├── test_state.c       # State tests
│   ├── test.h             # Test utilities
│   ├── test_qhipster.h    # QHipster tests
│   └── CMakeLists.txt     # Test build rules
├── cmake/                 # CMake modules
│   ├── PlatformConfig.cmake    # OS-specific configuration
│   ├── CompilerFlags.cmake     # Compiler and linker flags
│   └── Dependencies.cmake      # External dependencies
└── extern/                # External dependencies
    └── Unity/             # Testing framework (submodule)
```

## Core Components

### q_types.h
Fundamental type definitions for the quantum simulator:
- `qubit_t`: Qubit identifier
- `cplx_t`: Complex number representation
- `dim_t`: Dimension type
- Utility macros for common operations

### state.h
Quantum state structures and manipulation functions:
- `state_type_t`: Enum for PURE or MIXED state types
- `state_t`: Main state structure containing type, data array, and qubit count
- State manipulation functions:
  - `state_free()`: Free allocated memory
  - `state_len()`: Get array size for state representation
  - `state_init()`: Initialize a quantum state
  - `state_plus()`: Initialize uniform superposition
  - `state_cp()`: Create a deep copy of a state

### qlib.h
Main library header that includes all core components.

## Testing

The project uses the Unity testing framework (included as a git submodule).

```bash
# Run tests from build directory
ctest
```

## Project Status

The project is in active development. Current focus is on implementing core quantum state representation and manipulation. Additional modules are planned:
- Circuit operations
- Gate implementations (see `docs/GATES_MODULE.md`)
- Measurement operations
- Gradient computation for variational algorithms (see `docs/MEAS_MODULES.md`)

## Documentation

Technical specifications are available in the `docs/` directory:
- **STATE_MODULE.md**: State representation (pure/mixed), packed storage format, BLAS integration
- **GATES_MODULE.md**: Gate operations, test infrastructure, implementation phases
- **MEAS_MODULES.md**: Gradient computation methods (parameter-shift, backpropagation) for VQE/QAOA

## Theoretical Foundation

This implementation is based on the algorithms and complexity analysis from **PhD Thesis Chapter 2**:

**Location**: `~/Projects/thesis/2.Simulation/`

### Code-to-Theory Mapping

| qlib Component | Thesis Section | Theory Provided |
|----------------|----------------|-----------------|
| Packed column major format | sec4.tex | Memory layout, offset/stride formulas (Eq 6-10) |
| `src/qhipster.c` | sec3.tex | Pure state algorithm (Alg 4), HPC optimizations |
| `src/mhipster.c` | sec2.tex, sec5.tex | Bit-insertion (Alg 1), hermitian conjugation (Alg 6-8) |
| Gate dispatch strategy | sec1.tex | k-local complexity O(r·N²·2^k), Propositions 1-3 |
| Zero-flop Pauli gates | sec1.tex Ex 1, sec2.tex Alg 3 | Memory-only operations for X, CNOT |

### Key Theoretical Results

- **Proposition 3**: k-local channels achieve O(r·N²·2^k) complexity (vs O(r·N³) dense)
- **Algorithm 1**: Bit-insertion reconstructs global indices from block coordinates
- **Equation 12**: Update formula for single-qubit hermitian conjugation
- **Table 1**: 10³×–10⁷× speedup for 10-25 qubits with native gates

For full theoretical background, proofs, and complexity analysis: `cat ~/Projects/thesis/2.Simulation/NOTES.md`
