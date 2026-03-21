# CLAUDE.md

This file provides guidance to Claude Code when working with this repository.

## Build & Test

```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug
```

> **Warning:** Always use `--preset debug` (or explicitly pass `-DCMAKE_BUILD_TYPE=Debug`).
> The Debug build sets `LOG_TILE_DIM=3` (8×8 tiles), which is required for tests to exercise
> multi-tile code paths at small qubit counts. A Release build (`LOG_TILE_DIM=5`, 32×32 tiles)
> will silently skip those paths, causing tests to pass despite untested code.
> `CMakeLists.txt` enforces this: configuring a `*debug*` directory with a non-Debug build type
> (or vice versa) is a `FATAL_ERROR`.

## Repository Layout

```
QSim/                         # C17/C++17 quantum simulator
├── src/                      # Main library source (state/ and gate/ subdirs)
├── include/                  # Public headers (qlib.h, state.h, gate.h, …)
├── test/                     # CTest-based integration test suite
├── benchmark/                # Performance benchmarks
├── profile/                  # Profiling harness
├── extern/                   # Third-party deps: QuEST, qulacs (benchmark only)
├── cmake/
│   ├── PlatformConfig.cmake  # OS/architecture detection
│   ├── Dependencies.cmake    # find_package calls and extern wiring
│   └── CompilerFlags.cmake   # Interface library for compiler flags
├── docs/                     # Per-module technical specs (see Module table below)
├── tools/                    # Dev utilities (e.g. verify_ilp64.c)
├── archive/                  # Archived/experimental code – do not modify
├── CMakeLists.txt            # Root build definition
└── CMakePresets.json         # Defines debug and release presets
```

## Build Options

| Option | Default | Description |
|---|---|---|
| `BUILD_TESTS` | `ON` | Build CTest test suite |
| `BUILD_BENCHMARKS` | `OFF` | Build benchmark targets |
| `ENABLE_OPENMP` | `ON` | Enable OpenMP parallelization |
| `ENABLE_ASAN` | `OFF` | Enable AddressSanitizer + UBSan |
| `LOG_TILE_DIM` | `3` (Debug) / `5` (Release) | Log₂ of tile dimension — **do not override manually** |

## Modules

**Read only the Technical Spec for your module. Do not read any other Technical Spec. Do not read Thesis Notes unless explicitly asked to add or update them. Each Technical Spec must include a Build Integration section documenting its `CMakeLists.txt`.**

**Documentation update policy: After completing any task — and always before committing changes — update the Technical Spec for the affected module(s) to reflect the current state of the code. The spec must remain the authoritative, accurate description of what is actually implemented.**

| Module      | Status      | Technical Spec          | Thesis Notes |
|-------------|-------------|-------------------------|--------------|
| Build System | Complete   | `docs/BUILD_SYSTEM.md`  | —            |
| State       | Complete    | `docs/STATE_MODULE.md`  | `~/Projects/thesis/2.Simulation/notes/STATE_THESIS_NOTES.md` |
| Gate        | In Progress | `docs/GATE_MODULE.md`   | `~/Projects/thesis/2.Simulation/notes/GATES_THESIS_NOTES.md` |
| Channel     | In Progress | `docs/CHANNEL_MODULE.md`| `~/Projects/thesis/2.Simulation/notes/CHANNEL_THESIS_NOTES.md` |
| Measurement | Pending     | `docs/MEAS_MODULE.md`   | —            |
| Tests       | In Progress | `docs/TEST_SUITE.md`    | —            |
| Benchmark   | In Progress | `docs/BENCHMARK.md`     | `~/Projects/thesis/2.Simulation/notes/BENCHMARK_THESIS_NOTES.md` |
| Profile     | In Progress | `docs/PROFILE.md`       | —            |
| Circuit     | In Progress | `docs/CIRCUIT_MODULE.md`| —            |
