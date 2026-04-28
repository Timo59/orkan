# Test Suite

Cross-module test framework: shared fixtures, reference computation, build integration, and per-module test-architecture conventions.

Per-module test inventories are in each module's spec (`docs/<MODULE>_MODULE.md`). This file describes what's shared.

---

## Framework

[Unity](https://github.com/throwtheswitch/unity) v2.6.1, fetched via CMake `FetchContent` at configure time. `unity.c` is compiled directly into each test executable (not a shared library). Configured with `UNITY_INCLUDE_DOUBLE` and `UNITY_SUPPORT_64` to enable double-precision and 64-bit integer assertions.

---

## Build Type Requirement

Tests **must** be built and run with `CMAKE_BUILD_TYPE=Debug`. The project sets `LOG_TILE_DIM` based on build type:

| Build type | `LOG_TILE_DIM` | Tile size | Effect on tests |
|---|---|---|---|
| **Debug** | 3 | 8x8 | Multi-tile code paths exercised at small qubit counts |
| **Release** | 5 | 32x32 | Multi-tile paths silently skipped — false passes |

`CMakeLists.txt` enforces this: configuring a `*debug*` directory with a non-Debug build type (or vice versa) is a `FATAL_ERROR`.

---

## Directory Structure

```
test/
  CMakeLists.txt           Build definitions for all test executables
  include/                 Shared test headers
    test.h                 Constants (MAXQUBITS=4, PRECISION=1e-12, SQRT2/INVSQRTn);
                           custom assertion macros (TEST_ASSERT_COMPLEX_WITHIN,
                           TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL,
                           TEST_ASSERT_EQUAL_DOUBLE_ARRAY_TOL,
                           TEST_ASSERT_NOT_EQUAL_PTR)
    test_gate.h            Gate harness declarations; MAXQUBITS_TILED=6;
                           gate-constant arrays (XMAT, HMAT, CXMAT, CCXMAT, ISWAPMAT,
                           HHMAT); TEST_THETAS / NUM_THETAS;
                           testTwoFromMat / testTwoFromMatDouble / testTwoFromMatSymmetry
    test_pure_states.h     Pure state fixture declarations
    test_mixed_states.h    Mixed state fixture declarations
    test_mixed_utils.h     Pack/unpack helpers, tiled-state assertion utilities
    gatemat.h              Full N-qubit matrix builder declarations: mat_id,
                           mat_single_qubit_gate, mat_two_qubit_gate,
                           mat_three_qubit_gate; rotation builders mat_rx/ry/rz/p
    linalg.h               Linear algebra helpers: zmv, zumu, zkron
  state/                   State module tests (-> test_state)
  gate/                    Gate module tests (-> test_gate)
  channel/                 Channel module tests (-> test_channel)
  meas/                    Measurement module tests (-> test_meas)
  circ/                    Circuit module tests (-> test_circ)
  utility/                 Test helpers and fixtures
    gatemat.c              Full N-qubit matrix builders via Kronecker products
    linalg.c               zmv (zgemv), zumu (zgemm), zkron (column-major Kronecker);
                           CBLAS header selected by platform
    test_pure_states.c     Pure state fixtures: cb, xb, yb, Bell, GHZ, W
    test_mixed_states.c    Mixed state fixtures: cb, xb, yb, Bell, GHZ, W,
                           maximally-mixed, random mixture
    test_mixed_utils.c     density_unpack/pack, tiled_state_from_full,
                           assert_tiled_equals_full
```

---

## Dependencies

| Dependency | Use | Executables |
|---|---|---|
| State module (`include/state.h`) | `state_t`, allocation/init/copy/get/set | all |
| Gate module (`include/gate.h`) | All gate application functions | `test_gate` |
| Channel module (`include/channel.h`) | `kraus_t`, `superop_t`, `channel_1q` | `test_channel` |
| Measurement module (`include/meas.h`) | `mean`, `matel_diag` | `test_meas` |
| Circuit module (`include/circ.h`) | `exp_diag` | `test_circ` |
| CBLAS (platform) | `zgemv`, `zgemm` in `utility/linalg.c` for reference matrix products | `test_gate`, `test_channel`, `test_meas`, `test_circ` |
| Unity v2.6.1 (fetched) | Assertion macros, `RUN_TEST`, runner | all |

`test_state` does **not** link against CBLAS — `linalg.c` is excluded from that target.

`utility/linalg.c` selects the CBLAS header by platform: `<vecLib/cblas_new.h>` on Apple, `<cblas.h>` on Linux.

Internal types (`cplx_t`, fixture state types) are defined in test headers and not exported to production code.

---

## Build Integration

Configured in `test/CMakeLists.txt`. All executables are built when `BUILD_TESTS=ON` (default).

```cmake
add_executable(test_state    state/test_state.c    state/...)
add_executable(test_gate     gate/test_gate.c      gate/... utility/...)
add_executable(test_channel  channel/test_channel.c channel/... utility/...)
add_executable(test_meas     meas/test_meas.c      meas/...)
add_executable(test_circ     circ/test_circ.c      circ/... utility/...)
```

Unity source (`unity.c`) is compiled directly into each executable. All targets link the `${PROJECT_NAME_LOWER}` library. CBLAS is linked transitively via the library or directly via `linalg.c` where used.

CTest registers each executable as a single test entry:

```cmake
add_test(NAME test_state   COMMAND test_state)
add_test(NAME test_gate    COMMAND test_gate)
add_test(NAME test_channel COMMAND test_channel)
add_test(NAME test_meas    COMMAND test_meas)
add_test(NAME test_circ    COMMAND test_circ)
```

No per-test-function CTest labels are defined; all Unity tests within an executable run as one CTest entry.

`LOG_TILE_DIM` is passed as a compile definition to all targets. Modules that depend on it also receive `OMP_THRESHOLD`.

---

## Gate Test Architecture

Two-layer harness pattern:

- **`test_gate.c`** — single entry point. Contains `main()`, all named `test_<Gate>_<representation>()` wrappers (e.g. `test_X_pure`, `test_CX_tiled`), and all `RUN_TEST()` registrations.

- **`test_gate_<repr>_<Nq>.c`** — one harness per (state representation, gate arity). Each implements a parameterized driver (e.g. `testTwoQubitGateTiled`) that applies a gate, builds the reference full-matrix via `gatemat.c`, and asserts element-wise equality across all fixture states.

- **`test_gate_tiled_2q_mat.c`** — separate harness for the internal `two_from_mat` kernel, which accepts an arbitrary 4x4 matrix rather than a named gate function pointer. Implements `testTwoFromMat`, `testTwoFromMatDouble` (catches stale mirror-slot bugs), `testTwoFromMatSymmetry` (argument-order invariance). Defines `ISWAPMAT` and `HHMAT` constants used via `extern` in `test_gate.h`. Forward-declares `two_from_mat`, bypassing the public gate header.

Rotation matrix builders (`mat_rx`, `mat_ry`, `mat_rz`, `mat_p`, `mat_rx_pi3`) live in `utility/gatemat.c` and are declared in `include/gatemat.h`. They fill a 4-element column-major `cplx_t[4]` with the 2x2 gate matrix for a given angle, shared by pure, packed, and tiled 1Q rotation harnesses.

---

## Naming Conventions

- Unity test functions: `test_<Gate>_<repr>` for gate tests, `test_<operation>_<repr>` for state tests, `test_<feature>_<description>` for module-specific tests.
- Harness driver functions (not directly registered with `RUN_TEST`): upper-camel-case — `testSingleQubitGate`, `testTwoQubitGateTiled`, `testTwoFromMat`.

---

## Reference Computation

Tests verify backends against a BLAS-based ground truth at tolerance 1e-12:

- Gate tests: `zgemv` (pure), `zgemm` (mixed) on full N-qubit matrices built via Kronecker products in `gatemat.c`.
- Channel tests: explicit `sum_k K_k rho K_k-dagger` with full N-qubit Kraus matrices.
- Measurement tests: element-wise reduction on full state vector.
- Circuit tests: element-wise complex exponential.

---

## Running

```bash
cmake --preset debug          # Configure (enforces Debug build type)
cmake --build --preset debug  # Build
ctest --preset debug          # Run all tests

# Or run executables directly after building:
./cmake-build-debug/test/test_state
./cmake-build-debug/test/test_gate
./cmake-build-debug/test/test_channel
./cmake-build-debug/test/test_meas
./cmake-build-debug/test/test_circ
```
