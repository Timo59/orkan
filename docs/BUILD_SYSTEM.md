# Build System

Shared build infrastructure for QSim. Per-module build details are in their respective Technical Specs.

---

## File Hierarchy

```
CMakeLists.txt            Root orchestrator
CMakePresets.json         Preset definitions (debug / release)
cmake/
├── PlatformConfig.cmake  OS detection and platform-specific settings
├── Dependencies.cmake    BLAS/LAPACK resolution (tests/benchmarks only)
└── CompilerFlags.cmake   Interface libraries for compiler flags
src/
├── CMakeLists.txt        Library target (shared)
└── internal/             Internal headers (index.h, utils.h) — not installed
test/CMakeLists.txt       Unity test framework; test executables
benchmark/CMakeLists.txt  Benchmark executables; auto-detects QuEST, Qulacs, Aer
profile/CMakeLists.txt    Micro-profiling executables
tools/
└── verify_ilp64.c        Standalone utility to validate ILP64 BLAS ABI
extern/                   Third-party sources (benchmark only): QuEST, qulacs, aer-dm
```

---

## Dependencies

| Dependency | When needed | Resolution |
|---|---|---|
| C17 / C++17 compiler | Always | `CMAKE_C_STANDARD 17` / `CMAKE_CXX_STANDARD 17` |
| CMake >= 3.27 | Always | `cmake_minimum_required(VERSION 3.27)` |
| OpenMP | Library (optional, `ENABLE_OPENMP=ON`) | Homebrew `libomp` on macOS; system package on Linux |
| BLAS/LAPACK ILP64 | Tests and benchmarks only | Accelerate (macOS) or FetchContent OpenBLAS (Linux) |
| Unity v2.6.1 | Tests only | FetchContent (git, `EXCLUDE_FROM_ALL`) |
| QuEST, Qulacs, Aer | Benchmark (optional) | Auto-detected from `extern/` or system paths |

BLAS is **not** a dependency of the library itself. It is only linked by test and benchmark targets. `Dependencies.cmake` is conditionally included when `BUILD_TESTS` or `BUILD_BENCHMARKS` is ON.

---

## Build Options

| Option | Default | Description |
|---|---|---|
| `BUILD_TESTS` | `ON` | Build CTest test suite |
| `BUILD_BENCHMARKS` | `OFF` | Build benchmark targets |
| `ENABLE_OPENMP` | `ON` | Enable OpenMP parallelization |
| `ENABLE_ASAN` | `OFF` | Enable AddressSanitizer + UBSan |
| `LOG_TILE_DIM` | `3` (Debug) / `5` (Release) | Log2 of tile dimension — do not override manually |

---

## Presets (CMakePresets.json)

| Preset | Build type | Binary dir | Tile dim |
|---|---|---|---|
| `debug` | Debug | `cmake-build-debug` | 8x8 |
| `release` | Release | `cmake-build-release` | 32x32 |

A directory/type mismatch (e.g. `cmake-build-debug` with Release) is a `FATAL_ERROR`.

---

## Targets

All target names are derived from the project name (`PROJECT_NAME_LOWER`).

| Target | Type | Output | Notes |
|---|---|---|---|
| `${PROJECT_NAME_LOWER}` | Shared library | `libqsim.dylib` / `libqsim.so` | State, gate, channel, circuit modules |
| `test_state` | Executable | `test/test_state` | State module tests |
| `test_gate` | Executable | `test/test_gate` | Gate module tests |
| `test_circuit` | Executable | `test/test_circuit` | Circuit module tests |
| `test_channel` | Executable | `test/test_channel` | Channel module tests |
| `bench_gate` | Executable | `benchmark/bench_gate` | Gate throughput benchmark (Release only) |
| `bench_kraus` | Executable | `benchmark/bench_kraus` | Kraus channel benchmark (Release only) |
| `profile_*` | Executables | `profile/` | Micro-profiling harnesses |
| `verified_install` | Custom | — | Runs tests, then installs |

---

## Interface Libraries (cmake/CompilerFlags.cmake)

| Library | Purpose | Linked by |
|---|---|---|
| `QSim_compiler_flags` | `-fno-strict-aliasing`, `-O0` (Debug), `-march=native` (Release), `-lm` (Linux) | Library (PUBLIC) |
| `blas_compiler_flags` | BLAS compile definitions, include dirs, link libraries | Tests, benchmarks (PRIVATE) |
| `omp_compiler_flags` | OpenMP flags and libraries | Library (PUBLIC) |

---

## Install

```bash
cmake --install cmake-build-release --prefix /usr/local
```

Installs:
- `lib/libqsim.{so,dylib}` with versioning (SOVERSION from `PROJECT_VERSION_MAJOR`)
- `include/${PROJECT_NAME_LOWER}/` containing: `qlib.h`, `q_types.h`, `state.h`, `gate.h`, `channel.h`

Users include headers as `#include <qsim/qlib.h>`. The subdirectory name is derived from the project name.

The `verified_install` target runs all tests before installing:
```bash
cmake --build cmake-build-debug --target verified_install
```

---

## Platform Notes

| Platform | BLAS | OpenMP | Position-independent code |
|---|---|---|---|
| macOS | Apple Accelerate (ILP64 native) | Homebrew `libomp` required | `CMAKE_MACOSX_RPATH=1` |
| Linux | FetchContent OpenBLAS with ILP64 | System package | `CMAKE_POSITION_INDEPENDENT_CODE=ON` |
