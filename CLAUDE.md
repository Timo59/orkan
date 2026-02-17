# CLAUDE.md

This file provides guidance to Claude Code when working with this repository.

## Build & Test

```bash
mkdir cmake-build-debug && cd cmake-build-debug
cmake ..
cmake --build .
ctest
```

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
