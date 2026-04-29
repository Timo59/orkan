# CLAUDE.md

Instructions for Claude Code agents working on Orkan.

---

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

---

## Repository Layout

```
orkan/                        # C17/C++17 quantum simulator
├── src/                      # Main library source (state/ and gate/ subdirs)
├── include/                  # Public headers (orkan.h, state.h, gate.h, …)
├── test/                     # CTest-based integration test suite
├── benchmark/                # Performance benchmarks
├── profile/                  # Profiling harness
├── extern/                   # Third-party deps: QuEST, qulacs (benchmark only)
├── cmake/
│   ├── PlatformConfig.cmake  # OS/architecture detection
│   ├── Dependencies.cmake    # find_package calls and extern wiring
│   └── CompilerFlags.cmake   # Interface library for compiler flags
├── docs/                     # Per-module technical specs (see Module Routing below)
├── tools/                    # Dev utilities (e.g. verify_ilp64.c)
├── archive/                  # Archived/experimental code – do not modify
├── CMakeLists.txt            # Root build definition
└── CMakePresets.json         # Defines debug and release presets
```

---

## Build Options

| Option | Default | Description |
|---|---|---|
| `BUILD_TESTS` | `ON` | Build CTest test suite |
| `BUILD_BENCHMARKS` | `OFF` | Build benchmark targets |
| `ENABLE_OPENMP` | `ON` | Enable OpenMP parallelization |
| `ENABLE_ASAN` | `OFF` | Enable AddressSanitizer + UBSan |
| `LOG_TILE_DIM` | `3` (Debug) / `5` (Release) | Log₂ of tile dimension — **do not override manually** |

---

## Module Routing

Before starting any task, identify which module it belongs to and read **only** that module's Technical Spec.

| Module       | Status      | Technical Spec           |
|--------------|-------------|--------------------------|
| Build System | Complete    | `docs/BUILD_SYSTEM.md`   |
| State        | Complete    | `docs/STATE_MODULE.md`   |
| Gate         | In Progress | `docs/GATE_MODULE.md`    |
| Channel      | In Progress | `docs/CHANNEL_MODULE.md` |
| Measurement  | In Progress | `docs/MEAS_MODULE.md`    |
| Circuit      | In Progress | `docs/CIRCUIT_MODULE.md` |
| Tests        | In Progress | `docs/TEST_SUITE.md`     |
| Benchmark    | In Progress | `docs/BENCHMARK.md`      |
| Profile      | In Progress | `docs/PROFILE.md`        |

If a task touches multiple modules, read each relevant spec but finish one module's changes before starting the next.

---

## Reading Rules

**Always read before any task:**
1. This file (`CLAUDE.md`)
2. The Technical Spec for the module(s) your task belongs to
3. `README.md` for build and install instructions

**Do NOT read unless explicitly instructed:**
- Other modules' Technical Specs
- Files under `archive/` (legacy / experimental)

---

## Behavioral Rules

1. **Never edit, build, or commit on `main` unless the user explicitly asks for it.** Before making any change to a tracked file, run `git rev-parse --abbrev-ref HEAD`. If the result is `main` and the user has not explicitly authorised work on `main`, stop and create or switch to a feature branch (`git checkout -b feat/<short-name>`) before proceeding. `main` is reserved for fast-forward merges of completed work.
2. After completing any task — and **always before committing changes** — update the Technical Spec for the affected module(s) to reflect the current state of the code. The spec must remain the authoritative description of what is actually implemented.
3. Each Technical Spec must include a **Build Integration** section documenting its `CMakeLists.txt`.
4. Do not extend code in `archive/`. Treat it as read-only.
5. Do not add historical data, changelogs, or update timestamps to documentation files. Documentation describes the current state only — git history serves that purpose.
