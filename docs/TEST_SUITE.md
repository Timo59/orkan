# Test Suite

Cross-module test infrastructure: framework, build-type requirement, shared fixtures and reference computation, CTest registration. **Per-module test inventories live in each module's spec** (`docs/<MODULE>_MODULE.md`).

---

## Framework

[Unity](https://github.com/throwtheswitch/unity) v2.6.1, fetched via CMake `FetchContent` at configure time. `unity.c` is compiled directly into each test executable (not a shared library). Configured with `UNITY_INCLUDE_DOUBLE` and `UNITY_SUPPORT_64` to enable double-precision and 64-bit integer assertions.

---

## Build Type Requirement

Tests **must** be built and run with `CMAKE_BUILD_TYPE=Debug`. `LOG_TILE_DIM` is set by build type:

| Build type | `LOG_TILE_DIM` | Tile size | Effect on tests |
|---|---|---|---|
| **Debug** | 3 | 8x8 | Multi-tile code paths exercised at small qubit counts |
| **Release** | 5 | 32x32 | Multi-tile paths silently skipped — false passes |

`CMakeLists.txt` enforces this: configuring a `*debug*` directory with a non-Debug build type (or vice versa) is a `FATAL_ERROR`.

---

## Shared Test Infrastructure

Each module's test directory is named after the module (`test/<module>/`) and is documented in that module's spec. The files below are shared across all five test executables.

```
test/include/                Shared headers
  test.h                     Constants (MAXQUBITS=4, PRECISION=1e-12, SQRT2/INVSQRTn);
                             custom assertion macros (TEST_ASSERT_COMPLEX_WITHIN,
                             TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL,
                             TEST_ASSERT_EQUAL_DOUBLE_ARRAY_TOL,
                             TEST_ASSERT_NOT_EQUAL_PTR)
  test_gate.h                Gate harness declarations; MAXQUBITS_TILED=6;
                             gate-constant arrays (XMAT, HMAT, CXMAT, CCXMAT, ISWAPMAT,
                             HHMAT); TEST_THETAS / NUM_THETAS
  test_pure_states.h         Pure state fixture declarations
  test_mixed_states.h        Mixed state fixture declarations
  test_mixed_utils.h         Pack/unpack helpers, tiled-state assertion utilities
  gatemat.h                  Full N-qubit matrix builder declarations
  linalg.h                   Linear algebra helpers: zmv, zumu, zkron

test/utility/                Shared implementations
  gatemat.c                  Full N-qubit matrix builders via Kronecker products;
                             rotation builders mat_rx/ry/rz/p, mat_rx_pi3
  linalg.c                   zmv (zgemv), zumu (zgemm), zkron (column-major Kronecker);
                             CBLAS header selected by platform (vecLib on Apple,
                             cblas.h on Linux)
  test_pure_states.c         Pure state fixtures: cb, xb, yb, Bell, GHZ, W
  test_mixed_states.c        Mixed state fixtures: cb, xb, yb, Bell, GHZ, W,
                             maximally-mixed, random mixture
  test_mixed_utils.c         density_unpack/pack, tiled_state_from_full,
                             assert_tiled_equals_full
```

Internal types (`cplx_t`, fixture state types) are defined in test headers and not exported to production code.

---

## Reference Computation

All tests verify backends against a BLAS-based ground truth at tolerance **1e-12**:

| Module | Reference strategy |
|---|---|
| State | No BLAS; tests check storage layout directly |
| Gate | `zgemv` (pure), `zgemm` (mixed) on full N-qubit matrices built by Kronecker products in `gatemat.c` |
| Channel | Explicit `sum_k K_k rho K_k-dagger` with full N-qubit Kraus matrices |
| Measurement | Element-wise reduction on full state vector |
| Circuit | Element-wise complex exponential |

`test_state` does **not** link against CBLAS — `linalg.c` is excluded from that target.

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

`LOG_TILE_DIM` is passed as a compile definition to all targets. Modules that depend on it also receive `OMP_THRESHOLD`.

CTest registers each executable as a single test entry — no per-test-function CTest labels:

```cmake
add_test(NAME test_state   COMMAND test_state)
add_test(NAME test_gate    COMMAND test_gate)
add_test(NAME test_channel COMMAND test_channel)
add_test(NAME test_meas    COMMAND test_meas)
add_test(NAME test_circ    COMMAND test_circ)
```

All Unity tests within an executable run as one CTest pass/fail.

---

## Naming Conventions

- Unity test functions: `test_<Gate>_<repr>` for gate tests, `test_<operation>_<repr>` for state tests, `test_<feature>_<description>` for module-specific tests.
- Harness driver functions (not directly registered with `RUN_TEST`): upper-camel-case — `testSingleQubitGate`, `testTwoQubitGateTiled`, `testTwoFromMat`.

Module-specific harness architectures (e.g. the gate module's two-layer harness) are documented in the respective module spec.

---

## Running

```bash
cmake --preset debug          # Configure (enforces Debug build type)
cmake --build --preset debug  # Build
ctest --preset debug          # Run all tests

# Or run executables directly:
./cmake-build-debug/test/test_state
./cmake-build-debug/test/test_gate
./cmake-build-debug/test/test_channel
./cmake-build-debug/test/test_meas
./cmake-build-debug/test/test_circ
```
