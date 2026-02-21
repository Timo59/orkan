# Test Suite

## Framework

[Unity](https://github.com/throwtheswitch/unity) v2.6.1 (C unit testing), fetched via CMake
`FetchContent` at configure time. `unity.c` is compiled directly into each test executable (not
a shared library). Unity is configured with `UNITY_INCLUDE_DOUBLE` and `UNITY_SUPPORT_64` to
enable double-precision and 64-bit integer assertions.

## Directory Structure

```
test/
├── CMakeLists.txt              Build definitions for both test executables (85 lines)
├── include/                    Shared test headers
│   ├── test.h                  Constants (MAXQUBITS=4, PRECISION=1e-12, SQRT2/INVSQRTn),
│   │                           custom assertion macros (COMPLEX_WITHIN, EQUAL_COMPLEX_ARRAY_TOL,
│   │                           EQUAL_DOUBLE_ARRAY_TOL, NOT_EQUAL_PTR) (154 lines)
│   ├── test_gate.h             Declarations for all gate harness functions and rotation matrix
│   │                           builders (mat_rx/ry/rz/p) (224 lines)
│   ├── test_pure_states.h      Pure state fixture declarations (114 lines)
│   ├── test_mixed_states.h     Mixed state fixture declarations (145 lines)
│   ├── test_mixed_utils.h      Pack/unpack helpers, tiled-state assertion utilities,
│   │                           MAXQUBITS_TILED=6 (85 lines)
│   ├── gatemat.h               Full N-qubit matrix builder declarations: mat_id,
│   │                           mat_single_qubit_gate, mat_two_qubit_gate,
│   │                           mat_three_qubit_gate (86 lines)
│   └── linalg.h                Linear algebra helper declarations: zmv, zumu, zkron (73 lines)
├── state/                      State module tests (→ test_state executable)
│   ├── test_state.c            Runner (main), dispatch tests; 41 RUN_TEST registrations (230 lines)
│   ├── test_state_pure.c       10 pure statevector tests (335 lines)
│   ├── test_state_packed.c     13 packed density matrix tests (448 lines)
│   └── test_state_tiled.c      16 tiled density matrix tests (532 lines)
├── gate/                       Gate module tests (→ test_gate executable)
│   ├── test_gate.c             Runner (main) + all named test_<Gate>_<repr>() functions
│   │                           + RUN_TEST registrations for 62 tests total (210 lines)
│   ├── test_gate_pure_1q.c     Pure 1Q harness: testSingleQubitGate, testRotationGate,
│   │                           mat_rx/ry/rz/p builders (270 lines)
│   ├── test_gate_pure_2q.c     Pure 2Q harness: testTwoQubitGate (123 lines)
│   ├── test_gate_pure_3q.c     Pure 3Q harness: testThreeQubitGate (127 lines)
│   ├── test_gate_packed_1q.c   Packed 1Q harness: testSingleQubitGateMixed,
│   │                           testRotationGateMixed (276 lines)
│   ├── test_gate_packed_2q.c   Packed 2Q harness: testTwoQubitGateMixed (150 lines)
│   ├── test_gate_packed_3q.c   Packed 3Q harness: testThreeQubitGateMixed (155 lines)
│   ├── test_gate_tiled_1q.c    Tiled 1Q harness: testSingleQubitGateTiled,
│   │                           testRotationGateTiled, testSingleQubitMatGateTiled,
│   │                           testSingleQubitMatGateTiledDouble (424 lines)
│   ├── test_gate_tiled_2q.c    Tiled 2Q harness: testTwoQubitGateTiled (128 lines)
│   └── test_gate_tiled_3q.c    Tiled 3Q harness: testThreeQubitGateTiled (147 lines)
└── utility/                    Test helpers and fixtures
    ├── gatemat.c               Full N-qubit matrix builders via Kronecker products and
    │                           basis-state expansion (194 lines)
    ├── linalg.c                zmv (M*v via zgemv), zumu (U*M*U† via zgemm),
    │                           zkron (A⊗B, column-major) (111 lines)
    ├── test_pure_states.c      Pure state fixtures: cb, xb, yb, Bell, GHZ, W states (394 lines)
    ├── test_mixed_states.c     Mixed state fixtures: cb, xb, yb, Bell, GHZ, W,
    │                           maximally-mixed, random mixture (678 lines)
    └── test_mixed_utils.c      density_unpack/pack, tiled_state_from_full,
                                assert_tiled_equals_full (109 lines)
```

## Gate Test Architecture

The gate test suite uses a two-layer harness pattern:

- **`test_gate.c`** — the single entry point. Contains `main()`, all named
  `test_<Gate>_<representation>()` wrapper functions (e.g. `test_X_pure`, `test_CX_tiled`,
  `test_CCX_packed`), and all `RUN_TEST()` registrations.

- **`test_gate_<repr>_<Nq>.c`** — one harness per (state representation, gate arity). Each file
  implements a single parameterized driver (e.g. `testTwoQubitGateTiled`) that applies a gate,
  constructs the reference full-matrix via `gatemat.c`, and asserts element-wise equality across
  all fixture states. The named wrappers in `test_gate.c` call these drivers with concrete gate
  function pointers and reference matrices.

This separation keeps per-gate boilerplate in one place while the harness is exercised uniformly
across all fixture states.

## Test Inventory

### `test_state` — 41 tests

| Group | Count | File |
|-------|------:|------|
| PURE state | 10 | `state/test_state_pure.c` |
| MIXED_PACKED state | 13 | `state/test_state_packed.c` |
| MIXED_TILED state | 16 | `state/test_state_tiled.c` |
| Dispatch layer | 2 | `state/test_state.c` |
| **Total** | **41** | |

### `test_gate` — 62 tests

| Group | Count | Representation | Gates |
|-------|------:|----------------|-------|
| Pure 1Q | 13 | pure | X, Y, Z, H, Hy, S, Sdg, T, Tdg, Rx, Ry, Rz, P |
| Pure 2Q | 4 | pure | CX, CY, CZ, SWAP |
| Pure 3Q | 1 | pure | CCX |
| Packed 1Q | 13 | packed | X, Y, Z, H, Hy, S, Sdg, T, Tdg, Rx, Ry, Rz, P |
| Packed 2Q | 4 | packed | CX, CY, CZ, SWAP |
| Packed 3Q | 1 | packed | CCX |
| Tiled 1Q | 13 | tiled | X, Y, Z, H, Hy, S, Sdg, T, Tdg, Rx, Ry, Rz, P |
| Tiled 2Q | 4 | tiled | CX, CY, CZ, SWAP |
| Tiled 3Q | 1 | tiled | CCX |
| Tiled `single_from_mat` | 8 | tiled | ID, X, Y, H, T, Rx(π/3), X·H, H·T |
| **Total** | **62** | | |

## Executables

| Executable   | Sources                        | Description |
|------------- |--------------------------------|-------------|
| `test_state` | `state/*` + Unity              | State allocation, initialization, copy, get/set |
| `test_gate`  | `gate/*` + `utility/*` + Unity | Gate application across all representations |

## Build Type Requirement

Tests **must** be built and run with `CMAKE_BUILD_TYPE=Debug`. The project sets `LOG_TILE_DIM`
based on build type:

| Build type | `LOG_TILE_DIM` | Tile size | Effect on tests |
|---|---|---|---|
| **Debug** | 3 | 8×8 | Multi-tile code paths exercised at small qubit counts ✅ |
| **Release** | 5 | 32×32 | Multi-tile paths silently skipped → false passes ❌ |

A Release build inside `cmake-build-debug` is possible if the CMake cache was previously
configured with `Release`. Use the preset commands below to guarantee the correct build type.

## Running

```bash
cmake --preset debug          # Configure (enforces Debug build type)
cmake --build --preset debug  # Build
ctest --preset debug          # Run all tests

# Or run executables directly after building:
./cmake-build-debug/test/test_state
./cmake-build-debug/test/test_gate
```
