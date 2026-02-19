# Build System

This document covers the shared build infrastructure for QSim. Per-module `CMakeLists.txt` files are documented in their respective Technical Specs under a Build Integration section.

---

## File Hierarchy

```
CMakeLists.txt            ← root orchestrator
CMakePresets.json         ← preset definitions (debug / release)
cmake/
├── PlatformConfig.cmake  ← OS detection and platform-specific settings
├── Dependencies.cmake    ← BLAS/LAPACK resolution
└── CompilerFlags.cmake   ← interface libraries for compiler flags
```

---

## CMakeLists.txt

Root orchestrator. Responsible for:

- **Project declaration**: name `QSim`, C17/C++17, version 1.0.0
- **Default build type**: falls back to Debug when no build type or preset is specified (single-config generators only)
- **Directory/type guard**: a `*debug*` binary directory paired with a non-Debug build type (or vice versa) is a `FATAL_ERROR` — prevents silent Release builds inside debug directories
- **`LOG_TILE_DIM`**: set to `3` in Debug (8×8 tiles) and `5` in Release (32×32 tiles); validated to be in [3, 8]; **do not override manually**
- **`ENABLE_ASAN`**: when ON, appends `-fsanitize=address,undefined -fno-omit-frame-pointer` to C, CXX, and linker flags
- **Subdirectories**: includes `src` unconditionally; `test` when `BUILD_TESTS=ON`; `benchmark` when `BUILD_BENCHMARKS=ON`; `profile` when the directory exists

### Build Options

| Option | Default | Description |
|---|---|---|
| `BUILD_TESTS` | `ON` | Build CTest test suite |
| `BUILD_BENCHMARKS` | `OFF` | Build benchmark targets |
| `ENABLE_OPENMP` | `ON` | Enable OpenMP parallelization |
| `ENABLE_ASAN` | `OFF` | Enable AddressSanitizer + UBSan |
| `LOG_TILE_DIM` | `3` / `5` | Log₂ of tile dimension — do not override manually |

---

## CMakePresets.json

Defines three preset types:

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
| macOS | `CMAKE_MACOSX_RPATH=1` — enables `@rpath`-based dynamic library lookup |
| Linux | `CMAKE_POSITION_INDEPENDENT_CODE=ON` — required for shared library targets |
| Other | `FATAL_ERROR` — Windows is not supported |

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
| macOS | Apple Accelerate framework — ILP64 guaranteed via `ACCELERATE_LAPACK_ILP64` |
| Linux (default) | Bundled OpenBLAS v0.3.27 via FetchContent, built with `INTERFACE64=1` |
| Linux (`USE_SYSTEM_OPENBLAS=ON`) | System OpenBLAS — must be built with `INTERFACE64=1`; validated at configure time |

> **Note:** `USE_SYSTEM_OPENBLAS` is declared in this file rather than the root `CMakeLists.txt` as it is a dependency-specific concern.

---

## cmake/CompilerFlags.cmake

Creates three interface libraries linked by all targets:

### `QSim_compiler_flags`

General project-wide flags:

| Flag | Condition |
|---|---|
| `-fno-strict-aliasing` | Always — required for correct complex number aliasing |
| `-O0` | Debug only — suppresses implicit optimisation |
| `-march=native` | Release only — enables CPU-specific optimisations |
| `-lm` | Linux only — math library (implicit on macOS) |

### `blas_compiler_flags`

Forwards the three `QSIM_BLAS_*` variables from `Dependencies.cmake` as compile definitions, include directories, and link libraries.

### `omp_compiler_flags`

Wires OpenMP when `ENABLE_OPENMP=ON`:

- **macOS**: uses `brew --prefix libomp` to locate Homebrew's `libomp` dynamically (works on both Apple Silicon and Intel), then sets hints for `find_package(OpenMP)`
- **Linux**: `find_package(OpenMP)` resolves without hints
- Exposes `OpenMP::OpenMP_C` and `OpenMP::OpenMP_CXX` imported targets
- When `ENABLE_OPENMP=OFF`: library is created but empty — targets link it unconditionally without needing conditionals
