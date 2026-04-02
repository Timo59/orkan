# Orkan

A fast quantum state simulator written in C17. Named after the *Orkan* — the strongest wind force classification in German meteorology.

Supports pure states (statevectors) and mixed states (density matrices) with three storage backends optimised for different use cases.

## Features

- **Pure states** as statevectors
- **Mixed states** as density matrices in two formats:
  - Packed lower-triangular (LAPACK-compatible layout)
  - Tiled/blocked (cache-optimised for gate operations)
- **18 quantum gates**: Pauli (X, Y, Z), Clifford (H, S, T, ...), rotations (Rx, Ry, Rz, P), controlled (CX, CY, CZ, SWAP), Toffoli (CCX)
- **Quantum channels** via superoperator formalism (Kraus representation)
- **OpenMP parallelisation** across all backends
- **Cross-platform**: macOS (Accelerate) and Linux (OpenBLAS)

## Requirements

- CMake >= 3.27
- C17-compatible compiler
- OpenMP (optional): Homebrew `libomp` on macOS, system package on Linux

## Building

```bash
cmake -B build
cmake --build build
```

Or using presets:

```bash
cmake --preset release
cmake --build --preset release
```

## Testing

```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug
```

Tests require BLAS (Accelerate on macOS, auto-downloaded OpenBLAS on Linux). The Debug preset uses small tiles (8x8) to exercise multi-tile code paths at low qubit counts.

## Installation

```bash
cmake --preset release
cmake --build --preset release
cmake --install cmake-build-release --prefix /usr/local
```

## Using Orkan in your project

After installation, add to your `CMakeLists.txt`:

```cmake
find_package(Orkan 0.1 REQUIRED)
add_executable(myapp main.c)
target_link_libraries(myapp Orkan::orkan)
```

Minimal example:

```c
#include <orkan/qlib.h>
#include <stdio.h>

int main(void) {
    // Create a 3-qubit pure state initialised to |+>^3
    state_t state = {.type = PURE, .data = NULL, .qubits = 0};
    state_plus(&state, 3);

    // Apply gates
    h(&state, 0);          // Hadamard on qubit 0
    cx(&state, 0, 1);      // CNOT: control=0, target=1
    rz(&state, 2, 0.5);    // Rz(0.5) on qubit 2

    // Read amplitudes
    for (idx_t i = 0; i < state_len(&state); i++) {
        cplx_t amp = state_get(&state, i, 0);
        printf("|%llu> : %.4f + %.4fi\n",
               (unsigned long long)i, creal(amp), cimag(amp));
    }

    state_free(&state);
    return 0;
}
```

## Project Structure

```
include/          Public headers (installed)
src/              Library sources (state, gate, channel modules)
  internal/       Internal headers (not installed)
test/             CTest-based test suite
benchmark/        Performance benchmarks
cmake/            CMake modules and package config template
```

## Benchmarks

Benchmarks compare Orkan against QuEST, Qulacs, and Qiskit-Aer on single-gate and Kraus channel throughput. Build with:

```bash
cmake --preset release -DBUILD_BENCHMARKS=ON
cmake --build --preset release --target bench_gate bench_kraus
./cmake-build-release/benchmark/bench_gate --help
```

## License

Copyright (C) 2025–2026 Timo Ziegler

This program is free software: you can redistribute it and/or modify it under the terms of the [GNU Affero General Public License v3.0](LICENSE).
