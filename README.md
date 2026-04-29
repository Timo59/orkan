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
- **Diagonal-observable expectation values** (`mean`, `matel_diag`)
- **Diagonal-unitary evolution** for parameterised circuits (`exp_diag`)
- **OpenMP parallelisation** across all backends
- **macOS and Linux** supported

## Requirements

| Component | When needed | Source |
|---|---|---|
| CMake >= 3.27 | Always | system package |
| C17 / C++17 compiler | Always | system |
| OpenMP | Library (optional) | Homebrew `libomp` on macOS, system package on Linux |
| BLAS/LAPACK (ILP64) | Tests & benchmarks only | Apple Accelerate (macOS), bundled OpenBLAS via FetchContent (Linux) |

The library itself does **not** require BLAS. BLAS is only linked into test and benchmark targets.

## Integration

Two ways to consume Orkan from another CMake project. Both pull only the library; BLAS, QuEST, Qulacs, and Aer are *not* required for downstream consumers.

### Option A — System install + `find_package`

```bash
cmake --preset release
cmake --build --preset release
cmake --install cmake-build-release --prefix /usr/local
```

In your project's `CMakeLists.txt`:

```cmake
find_package(Orkan 0.1 REQUIRED)
add_executable(myapp main.c)
target_link_libraries(myapp PRIVATE Orkan::orkan)
```

### Option B — Submodule + `add_subdirectory`

```bash
git submodule add https://github.com/Timo59/orkan.git extern/orkan
```

No `--recursive` is needed: Orkan's own submodules (QuEST, Qulacs) are only used to build *Orkan's* benchmarks and are not pulled in by `add_subdirectory`.

In your project's `CMakeLists.txt`:

```cmake
add_subdirectory(extern/orkan)
add_executable(myapp main.c)
target_link_libraries(myapp PRIVATE Orkan::orkan)
```

This path needs no system-wide install, no BLAS on the consumer machine, and propagates OpenMP transitively.

## Quick start

```c
#include <orkan/orkan.h>
#include <stdio.h>

int main(void) {
    // 3-qubit pure state initialised to |+>^3
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

Build with either integration path above.

## Building from source

```bash
cmake --preset release
cmake --build --preset release
```

`debug` and `release` presets differ in tile dimension (8x8 vs 32x32). The Debug preset is required for the test suite — see `docs/BUILD_SYSTEM.md` for details.

## Testing

```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug
```

Five test executables (`test_state`, `test_gate`, `test_channel`, `test_meas`, `test_circ`) verify each backend against a BLAS-based reference at tolerance 1e-12. Tests require BLAS (Accelerate on macOS, auto-fetched OpenBLAS on Linux).

## Benchmarks

Orkan's tiled and packed backends are benchmarked against [QuEST](https://github.com/QuEST-Kit/QuEST), [Qulacs](https://github.com/qulacs/qulacs), and the density-matrix kernel from [Qiskit Aer](https://github.com/Qiskit/qiskit-aer) on single-gate throughput and Kraus channel application.

### Reproducing the paper's results

Reference: [arXiv:2604.15765](https://arxiv.org/abs/2604.15765). Hardware and software environment used for the published numbers are specified in the paper itself; follow that environment for binary-comparable timings. The repository's job is the build/run procedure:

```bash
# 1. Clone with submodules (QuEST, Qulacs)
git clone --recurse-submodules https://github.com/Timo59/orkan.git
cd orkan
# Or, if already cloned without --recurse-submodules:
# git submodule update --init --recursive

# 2. Build benchmarks (Release required for representative timings)
cmake --preset release -DBUILD_BENCHMARKS=ON
cmake --build --preset release --target bench_gate bench_kraus

# 3. Run a benchmark and capture CSV
./cmake-build-release/benchmark/bench_gate H \
    --min-qubits 2 --max-qubits 12 --output csv > h_gate.csv
./cmake-build-release/benchmark/bench_kraus --channel depolarize \
    --output csv > kraus_depolarize.csv
```

The Qiskit Aer backend is vendored at `extern/aer-dm/` and may not be present in your clone. If absent, the Aer column reports `nan`; the other backends still run.

Full CLI surface, output formats (`console`, `csv`, `pgfplots`), bandwidth-estimation methodology, and timing-loop details are in `docs/BENCHMARK.md`.

## Project Structure

```
include/         Public headers: orkan.h, state.h, gate.h, channel.h, meas.h, circ.h
src/             Library sources (state, gate, channel, meas, circ modules)
  internal/      Internal headers (not installed)
test/            CTest-based test suite (state, gate, channel, meas, circ)
benchmark/       Performance benchmarks (bench_gate, bench_kraus)
profile/         Micro-profiling harnesses for individual gate kernels
cmake/           Platform detection, dependency resolution, package config template
extern/          Third-party sources for benchmarks: QuEST, qulacs, aer-dm
```

## Documentation

- Per-module specifications: `docs/<MODULE>_MODULE.md` (state, gate, channel, meas, circ)
- Build system, presets, install layout: `docs/BUILD_SYSTEM.md`
- Test framework, fixtures, reference strategy: `docs/TEST_SUITE.md`
- Benchmark CLI, methodology, output formats: `docs/BENCHMARK.md`
- Profiling framework: `docs/PROFILE.md`

## Citation

If you use Orkan in academic work, please cite:

> Ziegler, T. *Orkan: a fast quantum state simulator.* arXiv:2604.15765 (2026). https://arxiv.org/abs/2604.15765

## License

Copyright (C) 2025–2026 Timo Ziegler

This program is free software: you can redistribute it and/or modify it under the terms of the [GNU Affero General Public License v3.0](LICENSE).
