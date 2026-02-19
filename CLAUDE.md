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

## Modules

**Read only the Technical Spec for your module. Do not read any other Technical Spec. Do not read Thesis Notes unless explicitly asked to add or update them.**

| Module      | Status      | Technical Spec          | Thesis Notes |
|-------------|-------------|-------------------------|--------------|
| State       | Complete    | `docs/STATE_MODULE.md`  | `~/Projects/thesis/2.Simulation/notes/STATE_THESIS_NOTES.md` |
| Gate        | In Progress | `docs/GATE_MODULE.md`   | `~/Projects/thesis/2.Simulation/notes/GATES_THESIS_NOTES.md` |
| Measurement | Pending     | `docs/MEAS_MODULES.md`  | —            |
| Tests       | In Progress | `docs/TEST_SUITE.md`    | —            |
| Benchmark   | In Progress | `docs/BENCHMARK.md`     | `~/Projects/thesis/2.Simulation/notes/BENCHMARK_THESIS_NOTES.md` |
| Profile     | In Progress | `docs/PROFILE.md`       | —            |
