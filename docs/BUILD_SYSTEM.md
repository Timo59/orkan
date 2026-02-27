# Build System

This document covers the shared build infrastructure for QSim. Per-module `CMakeLists.txt` files are documented in their respective Technical Specs under a Build Integration section.

---

## File Hierarchy

```
CMakeLists.txt            ← root orchestrator
CMakePresets.json         ← preset definitions (debug / release)
cmake/
├── PlatformConfig.cmake  ← OS detection and platform-specific settings
├── Dependencies.cmake    ← BLAS/LAPACK resolution and FetchContent wiring
└── CompilerFlags.cmake   ← interface libraries for compiler flags
src/CMakeLists.txt        ← defines library targets q (shared) and circuit (static)
test/CMakeLists.txt       ← fetches Unity v2.6.1; defines test_state, test_gate, test_circuit
benchmark/CMakeLists.txt  ← defines bench_mixed; auto-detects QuEST and Qulacs from extern/
profile/CMakeLists.txt    ← defines profile_x_tiled, profile_h_tiled
tools/
└── verify_ilp64.c        ← standalone utility to validate ILP64 BLAS ABI at runtime
extern/                   ← third-party sources: QuEST and qulacs (benchmark only)
```

---

## Inputs and Dependencies

| Dependency | Required | Resolution |
|---|---|---|
| C17 / C++17 compiler | Always | Set via `CMAKE_C_STANDARD 17` / `CMAKE_CXX_STANDARD 17` after `project()` |
| CMake ≥ 3.27 | Always | Declared via `cmake_minimum_required(VERSION 3.27)` |
| BLAS/LAPACK ILP64 | Always | Accelerate (macOS), FetchContent OpenBLAS v0.3.27 (Linux default), or system OpenBLAS |
| OpenMP | Optional (`ENABLE_OPENMP=ON`) | Homebrew `libomp` on macOS (brew must be present), system package on Linux |
| Unity v2.6.1 | Tests only | FetchContent (git tag v2.6.1) — requires network or pre-populated cache |
| QuEST | Benchmark, optional | Auto-detected from `extern/QuEST` or system paths |
| Qulacs (csim only) | Benchmark, optional | Auto-detected from `extern/qulacs/src/csim`; built in-tree as `csim_static` |
| Eigen3 | Benchmark (Qulacs), optional | Searched in Homebrew and system paths |

---

## Outputs and Artifacts

| Artifact | Location | Notes |
|---|---|---|
| `q` shared library | `cmake-build-{debug,release}/src/` | Defined in `src/CMakeLists.txt` |
| `circuit` static library | `cmake-build-{debug,release}/src/` | Defined in `src/CMakeLists.txt`; links `q` PUBLIC |
| Test executables | `cmake-build-debug/test/` | Only when `BUILD_TESTS=ON` and `test/` exists |
| `bench_mixed` executable | `cmake-build-release/benchmark/` | Only when `BUILD_BENCHMARKS=ON` and Release build |
| Profile executables | `cmake-build-{debug,release}/profile/` | Included whenever `profile/` directory exists |
| `QSim_compiler_flags` | CMake interface target | Consumed by all QSim targets via `target_link_libraries` |
| `blas_compiler_flags` | CMake interface target | Consumed by all QSim targets via `target_link_libraries` |
| `omp_compiler_flags` | CMake interface target | Consumed by all QSim targets via `target_link_libraries` |

The three interface libraries (`QSim_compiler_flags`, `blas_compiler_flags`, `omp_compiler_flags`) are created in `cmake/CompilerFlags.cmake`. The `q` target links all three as **PUBLIC**, so downstream targets that link `q` inherit them transitively.

---

## CMakeLists.txt

Root orchestrator. Responsible for:

- **Project declaration**: name `QSim`, C17/C++17, version 1.0.0
- **Default build type**: falls back to Debug when no build type or preset is specified (single-config generators only)
- **Directory/type guard**: a `*debug*` or `*Debug*` binary directory paired with a non-Debug build type (or vice versa for Release) is a `FATAL_ERROR` — prevents silent Release builds inside debug directories
- **`LOG_TILE_DIM`**: set to `3` in Debug (8×8 tiles) and `5` in Release (32×32 tiles); validated to be in [3, 8]; **do not override manually**
- **`OMP_THRESHOLD`**: set to `4` in Debug builds only — overrides the per-file default minimum qubit count for OpenMP parallelization; not set in Release (per-file defaults apply)
- **`ENABLE_ASAN`**: when ON, appends `-fsanitize=address,undefined -fno-omit-frame-pointer` to C, CXX, EXE linker, and shared linker flags globally
- **Subdirectories**: includes `src` unconditionally; `test` when `BUILD_TESTS=ON` AND `test/` exists; `benchmark` when `BUILD_BENCHMARKS=ON` AND `benchmark/` exists; `profile` when `profile/` exists (no CMake option)
- **`verified_install` custom target**: runs CTest then `cmake --install`; depends on `test_state` and `test_gate`

### Build Options

| Option | Default | Description |
|---|---|---|
| `BUILD_TESTS` | `ON` | Build CTest test suite |
| `BUILD_BENCHMARKS` | `OFF` | Build benchmark targets |
| `ENABLE_OPENMP` | `ON` | Enable OpenMP parallelization |
| `ENABLE_ASAN` | `OFF` | Enable AddressSanitizer + UBSan |
| `USE_SYSTEM_OPENBLAS` | `OFF` | Linux only: declared in `Dependencies.cmake`; use system OpenBLAS instead of FetchContent build; must be ILP64 |
| `LOG_TILE_DIM` | `3` / `5` | Log₂ of tile dimension — **do not override manually** |

---

## CMakePresets.json

Defines three preset types (schema version 3):

| Type | Names | Effect |
|---|---|---|
| `configurePresets` | `debug`, `release` | Sets `binaryDir` (`cmake-build-debug` / `cmake-build-release`) and `CMAKE_BUILD_TYPE` — this is what the root guard validates against |
| `buildPresets` | `debug`, `release` | References matching configure preset |
| `testPresets` | `debug` only | References debug configure preset; enables `outputOnFailure` |

No `release` test preset exists by design — tests are only run against Debug builds.

---

## cmake/PlatformConfig.cmake

Applies platform-specific settings early in the configure step:

| Platform | Settings |
|---|---|
| macOS (`APPLE`) | `CMAKE_MACOSX_RPATH=1` — enables `@rpath`-based dynamic library lookup |
| Linux (`UNIX` and not `APPLE`) | `CMAKE_POSITION_INDEPENDENT_CODE=ON` — required for shared library targets |
| Other | `FATAL_ERROR "Unknown system: ${CMAKE_SYSTEM_NAME}"` |

---

## cmake/Dependencies.cmake

Resolves the BLAS/LAPACK dependency and exports three variables consumed by `CompilerFlags.cmake`:

| Variable | Purpose |
|---|---|
| `QSIM_BLAS_COMPILE_DEFINITIONS` | Preprocessor defines for ILP64 ABI |
| `QSIM_BLAS_LIBRARIES` | Library or framework to link against |
| `QSIM_BLAS_INCLUDE_DIRS` | Include paths |

ILP64 (64-bit integer BLAS interface) is required throughout. Resolution strategy by platform:

| Platform | Strategy |
|---|---|
| macOS | Apple Accelerate framework — ILP64 via `ACCELERATE_NEW_LAPACK` + `ACCELERATE_LAPACK_ILP64` defines |
| Linux (default) | Bundled OpenBLAS v0.3.27 via FetchContent (tarball from GitHub); built with `INTERFACE64=1` and `DYNAMIC_ARCH=ON`; requires network on first configure |
| Linux (`USE_SYSTEM_OPENBLAS=ON`) | System OpenBLAS searched in `/usr/lib`, `/usr/local/lib`, `/opt/OpenBLAS/lib`; FATAL_ERROR if not found; user is responsible for ILP64 build |

`USE_SYSTEM_OPENBLAS` is declared in this file (not the root `CMakeLists.txt`) as it is a dependency-specific concern.

> **Offline note:** When `USE_SYSTEM_OPENBLAS=OFF` on Linux, FetchContent downloads a tarball from GitHub. An offline configure will fail unless the CMake FetchContent cache (`_deps/` in the build directory) is pre-populated.

---

## cmake/CompilerFlags.cmake

Creates three interface libraries linked by all targets via `target_link_libraries(...PUBLIC ...)` on the `q` target.

### `QSim_compiler_flags`

General project-wide flags:

| Flag | Condition |
|---|---|
| `-fno-strict-aliasing` | Always — required for correct complex number aliasing |
| `-O0` | Debug only — suppresses implicit optimisation |
| `-march=native` | Release only — enables CPU-specific optimisations |
| `-lm` | Linux only (`UNIX AND NOT APPLE`) — math library (implicit on macOS) |

### `blas_compiler_flags`

Forwards the three `QSIM_BLAS_*` variables from `Dependencies.cmake` as compile definitions, include directories, and link libraries.

### `omp_compiler_flags`

Wires OpenMP when `ENABLE_OPENMP=ON`:

- **macOS**: runs `brew --prefix libomp` to locate Homebrew's `libomp`; if `brew` is not found, issues `FATAL_ERROR`. Sets `OpenMP_C_FLAGS`, `OpenMP_CXX_FLAGS`, `OpenMP_C_LIB_NAMES`, `OpenMP_CXX_LIB_NAMES`, and `OpenMP_omp_LIBRARY` as hints for `find_package(OpenMP REQUIRED)`. Also explicitly adds `${LIBOMP_PREFIX}/include` as an include directory (FindOpenMP does not always propagate it on macOS).
- **Linux**: `find_package(OpenMP REQUIRED COMPONENTS C CXX)` resolves without hints.
- Exposes `OpenMP::OpenMP_C` and `OpenMP::OpenMP_CXX` imported targets.

When `ENABLE_OPENMP=OFF`, the target is created but empty — module targets link it unconditionally without needing conditionals.

---

## src/CMakeLists.txt

Defines two library targets. See **STATE_MODULE.md** and **GATE_MODULE.md** for full source lists and gate organisation.

### Target: `q` (shared library)

- Sources: `state/state.c`, `state/state_pure.c`, `state/state_packed.c`, `state/state_tiled.c`, `gate/gate.c`, `gate/gate_pure.c`, all `gate/packed/gate_packed_*.c`, all `gate/tiled/gate_tiled_*.c`
- `VERSION` and `SOVERSION` set from `PROJECT_VERSION`
- Compile definitions (PRIVATE): `LOG_TILE_DIM`, `OMP_THRESHOLD` (Debug only)
- Links (PUBLIC): `QSim_compiler_flags`, `blas_compiler_flags`, `omp_compiler_flags`
- Install rules: library to `${CMAKE_INSTALL_LIBDIR}`; headers (`qlib.h`, `q_types.h`, `state.h`, `gate.h`, `circuit.h`) to `${CMAKE_INSTALL_INCLUDEDIR}`

### Target: `circuit` (static library)

- Sources: `circuit/pauli.c`
- Private include: `gate/tiled/` (for tiled gate internals)
- Compile definitions (PRIVATE): `LOG_TILE_DIM`, `OMP_THRESHOLD` (Debug only)
- Links (PUBLIC): `q`; (PRIVATE): `omp_compiler_flags` when `ENABLE_OPENMP=ON`

---

## test/CMakeLists.txt

Fetches Unity v2.6.1 (git) via FetchContent as the test framework. Defines three CTest executables:

| Target | Sources (abbreviated) | Links |
|---|---|---|
| `test_state` | `state/test_state*.c` | `q` (PRIVATE) |
| `test_gate` | `gate/test_gate_*.c`, `utility/` helpers | `q` (PRIVATE) |
| `test_circuit` | `circuit/test_circuit.c`, `test_pauli.c`, `utility/linalg.c` | `circuit`, `q` (PRIVATE) |

All targets define `UNITY_INCLUDE_DOUBLE`, `UNITY_SUPPORT_64`, and `LOG_TILE_DIM`. `test_circuit` also defines `OMP_THRESHOLD` when set.

---

## benchmark/CMakeLists.txt

Issues an early `return()` if the build type is not Release (benchmarks require optimised builds).

Auto-detects optional third-party frameworks:

| Framework | Detection | How built |
|---|---|---|
| QuEST | `extern/QuEST` or system paths | Pre-built; linked as `${QuEST_LIBRARY}` |
| Qulacs (csim) | `extern/qulacs/src/csim` | Built in-tree as `csim_static` (STATIC); needs Eigen3 |
| OpenMP | Homebrew hints on macOS, then `find_package(OpenMP)` | Used by `csim_static` and `bench_mixed` |

Defines `bench_mixed` executable (output: `${CMAKE_BINARY_DIR}/benchmark/`). Conditionally compiled with `WITH_QUEST` and/or `WITH_QULACS` defines depending on what is found.

Also defines a `bench` custom target (when `BUILD_TESTS=ON`) that runs CTest then `bench_mixed`.

---

## profile/CMakeLists.txt

Defines two micro-profiling executables (output: `${CMAKE_BINARY_DIR}/profile/`):

| Target | Source |
|---|---|
| `profile_x_tiled` | `src/profile_x_tiled.c` |
| `profile_h_tiled` | `src/profile_h_tiled.c` |

Both link `q` PRIVATE. No build-type guard — can build in Debug or Release.

---

## Install Support

`cmake --install` is supported. Rules are defined in `src/CMakeLists.txt`:
- Shared library `q` → `${CMAKE_INSTALL_LIBDIR}`
- Public headers → `${CMAKE_INSTALL_INCLUDEDIR}`

The `verified_install` custom target in the root runs CTest first, then performs the install.

No CPack configuration exists.

---

## Status

**Complete.** The build system is stable and fully described by this spec. No known missing functionality.

---

## Open Questions

1. ~~**Minimum CMake version**~~ — **Resolved**: `cmake_minimum_required(VERSION 3.27)` in root `CMakeLists.txt`.
2. ~~**Linkage model for interface libraries**~~ — **Resolved**: `q` links all three interface libraries as PUBLIC; downstream targets that link `q` inherit them transitively. The `circuit` target additionally links `omp_compiler_flags` PRIVATE.
3. ~~**Install / packaging support**~~ — **Resolved**: `cmake --install` is supported via rules in `src/CMakeLists.txt`. No CPack configuration exists.
4. **`profile/` subdirectory**: included conditionally "when the directory exists" rather than behind a CMake option. The rationale for directory-existence detection instead of a `BUILD_PROFILE` option (consistent with `BUILD_TESTS` / `BUILD_BENCHMARKS`) is not stated.
5. ~~**FetchContent network behaviour**~~ — **Resolved**: FetchContent downloads a tarball from GitHub on first configure. An offline configure will fail unless the `_deps/` cache is pre-populated. No explicit offline fallback is implemented.
6. **`verify_ilp64.c` integration**: the file exists under `tools/` but it is not invoked from any CMakeLists.txt. It appears to be a standalone manual-run utility; confirm whether it is used in CI or only ad hoc.
