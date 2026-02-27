# Test Suite

## Overview

Three CTest-registered executables (`test_state`, `test_gate`, `test_circuit`) validate the
State, Gate, and Circuit modules across all three state representations (pure statevector, packed
density matrix, tiled density matrix). Tests are written in C17 using the Unity framework and
must always be built and run with `--preset debug`.

**Status: In Progress** — State, Gate, and Circuit module tests are active. Measurement module
tests do not yet exist.

## Framework

[Unity](https://github.com/throwtheswitch/unity) v2.6.1 (C unit testing), fetched via CMake
`FetchContent` at configure time. `unity.c` is compiled directly into each test executable (not
a shared library). Unity is configured with `UNITY_INCLUDE_DOUBLE` and `UNITY_SUPPORT_64` to
enable double-precision and 64-bit integer assertions.

## Directory Structure

```
test/
├── CMakeLists.txt              Build definitions for all three test executables (121 lines)
├── include/                    Shared test headers
│   ├── test.h                  Constants (MAXQUBITS=4, PRECISION=1e-12, SQRT2/INVSQRTn),
│   │                           custom assertion macros (TEST_ASSERT_COMPLEX_WITHIN,
│   │                           TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL,
│   │                           TEST_ASSERT_EQUAL_DOUBLE_ARRAY_TOL,
│   │                           TEST_ASSERT_NOT_EQUAL_PTR) (154 lines)
│   ├── test_gate.h             Declarations for all gate harness functions; MAXQUBITS_TILED=6;
│   │                           gate-constant arrays (XMAT, HMAT, CXMAT, CCXMAT, ISWAPMAT extern,
│   │                           HHMAT extern, …); TEST_THETAS/NUM_THETAS;
│   │                           testTwoFromMat/testTwoFromMatDouble/testTwoFromMatSymmetry
│   │                           declarations (228 lines)
│   ├── test_pure_states.h      Pure state fixture declarations (114 lines)
│   ├── test_mixed_states.h     Mixed state fixture declarations (145 lines)
│   ├── test_mixed_utils.h      Pack/unpack helpers, tiled-state assertion utilities;
│   │                           packed_index() static inline (106 lines)
│   ├── gatemat.h               Full N-qubit matrix builder declarations: mat_id,
│   │                           mat_single_qubit_gate, mat_two_qubit_gate,
│   │                           mat_three_qubit_gate; rotation matrix builders:
│   │                           mat_rx/ry/rz/p, mat_rx_pi3 (95 lines)
│   └── linalg.h                Linear algebra helper declarations: zmv, zumu, zkron (73 lines)
├── state/                      State module tests (→ test_state executable)
│   ├── test_state.c            Runner (main), dispatch tests; 41 RUN_TEST registrations (230 lines)
│   ├── test_state_pure.c       10 pure statevector tests (335 lines)
│   ├── test_state_packed.c     13 packed density matrix tests (448 lines)
│   └── test_state_tiled.c      16 tiled density matrix tests (532 lines)
├── gate/                       Gate module tests (→ test_gate executable)
│   ├── test_gate.c             Runner (main) + all named test_<Gate>_<representation>() wrappers
│   │                           + all RUN_TEST registrations for 70 tests total (231 lines)
│   ├── test_gate_pure_1q.c     Pure 1Q harness: testSingleQubitGate, testRotationGate (221 lines)
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
│   ├── test_gate_tiled_2q_mat.c Tiled two_from_mat harnesses: testTwoFromMat,
│   │                           testTwoFromMatDouble, testTwoFromMatSymmetry;
│   │                           also defines ISWAPMAT[16] and HHMAT[16] (408 lines)
│   └── test_gate_tiled_3q.c    Tiled 3Q harness: testThreeQubitGateTiled (130 lines)
├── circuit/                    Circuit module tests (→ test_circuit executable)
│   ├── test_circuit.c          Runner (main); 30 RUN_TEST registrations; setUp/tearDown stubs
│   │                           (168 lines)
│   └── test_pauli.c            All test function implementations: pauli construction/masks,
│   │                           pauli_phase, pauli_apply_pure (n=1..4), from_str round-trip,
│   │                           involution, identity edge case (1310 lines)
└── utility/                    Test helpers and fixtures
    ├── gatemat.c               Full N-qubit matrix builders via Kronecker products and
    │                           basis-state expansion; rotation matrix builders
    │                           mat_rx/ry/rz/p, mat_rx_pi3 (251 lines)
    ├── linalg.c                zmv (M*v via zgemv), zumu (U*M*U† via zgemm),
    │                           zkron (A⊗B, column-major); CBLAS header selected by
    │                           platform: <vecLib/cblas_new.h> (Apple),
    │                           <cblas.h> (Linux) (111 lines)
    ├── test_pure_states.c      Pure state fixtures: cb, xb, yb, Bell, GHZ, W states (394 lines)
    ├── test_mixed_states.c     Mixed state fixtures: cb, xb, yb, Bell, GHZ, W,
    │                           maximally-mixed, random mixture (653 lines)
    └── test_mixed_utils.c      density_unpack/pack (use packed_index), tiled_state_from_full,
                                assert_tiled_equals_full (104 lines)
```

## Dependencies

### Consumed from other modules

| Dependency | What is used | Executables |
|---|---|---|
| State module (`include/state.h`) | `State`, `StateType`, allocation/init/copy/get/set API | all three |
| Gate module (`include/gate.h`) | All gate application functions under test | `test_gate` |
| Circuit module (`include/circuit.h`) | `pauli_t`, `pauli_new`, `pauli_from_str`, `pauli_phase`, `pauli_apply_pure`, `pauli_free` | `test_circuit` |
| CBLAS (platform) | `zgemv`, `zgemm` used by `utility/linalg.c` for reference matrix products | `test_gate`, `test_circuit` |
| Unity v2.6.1 (fetched) | Assertion macros, `RUN_TEST`, test runner infrastructure | all three |

`test_state` does **not** link against CBLAS; `linalg.c` is not compiled into it.

### Internal test-only types

`cplx_t` (complex double) and fixture state types are defined in test headers and are not exported
to production code.

## Build Integration

Configured in `test/CMakeLists.txt`. All three executables are built when `BUILD_TESTS=ON`
(default).

```
add_executable(test_state   state/test_state.c  state/test_state_pure.c  …)
add_executable(test_gate    gate/test_gate.c    gate/…  utility/…)
add_executable(test_circuit circuit/test_circuit.c  circuit/test_pauli.c
                            utility/linalg.c)
```

Unity source (`unity.c`) is compiled directly into each executable. `test_state` and `test_gate`
link against the `q` library target only. `test_circuit` links against both the `circuit` and
`q` library targets. All three link the platform CBLAS library transitively via `q` or directly
through `linalg.c`.

CTest registers each executable as a single test entry:

```
add_test(NAME test_state   COMMAND test_state)
add_test(NAME test_gate    COMMAND test_gate)
add_test(NAME test_circuit COMMAND test_circuit)
```

No per-test-function CTest labels are defined; all Unity tests run as one CTest entry per
executable.

`LOG_TILE_DIM` is passed as a compile definition to all three targets. `test_circuit` also
receives `OMP_THRESHOLD` if defined.

## Gate Test Architecture

The gate test suite uses a two-layer harness pattern:

- **`test_gate.c`** — the single entry point. Contains `main()`, all named
  `test_<Gate>_<representation>()` wrapper functions (e.g. `test_X_pure`, `test_CX_tiled`,
  `test_CCX_packed`, `test_two_from_mat_cx`), and all `RUN_TEST()` registrations.

- **`test_gate_<repr>_<Nq>.c`** — one harness per (state representation, gate arity). Each file
  implements a single parameterized driver (e.g. `testTwoQubitGateTiled`) that applies a gate,
  constructs the reference full-matrix via `gatemat.c`, and asserts element-wise equality across
  all fixture states.

- **`test_gate_tiled_2q_mat.c`** — a separate harness for the internal `two_from_mat` kernel,
  which accepts an arbitrary 4×4 matrix rather than a named gate function pointer. Implements
  `testTwoFromMat` (single application over all qubit pairs in [2, MAXQUBITS_TILED]),
  `testTwoFromMatDouble` (double application, catches stale mirror-slot bugs), and
  `testTwoFromMatSymmetry` (verifies argument-order invariance). Also defines `ISWAPMAT` and
  `HHMAT` constants used via `extern` declarations in `test_gate.h`. Contains a forward
  declaration of `two_from_mat` that bypasses the public gate header.

**Rotation matrix builders** (`mat_rx`, `mat_ry`, `mat_rz`, `mat_p`, `mat_rx_pi3`) are defined
in `utility/gatemat.c` and declared in `include/gatemat.h`. They fill a 4-element column-major
`cplx_t[4]` with the 2×2 gate matrix for a given angle, and are shared by the pure, packed, and
tiled 1Q rotation harnesses.

## Circuit Test Architecture

`test_circuit.c` is the Unity runner (entry point, `main()`, `setUp`/`tearDown` stubs). All test
function implementations live in `test_pauli.c`, which is compiled as a separate translation unit.
`test_circuit.c` contains only forward declarations and `RUN_TEST()` registrations.

The Pauli tests use a `check_pauli_apply` helper (static, defined in `test_pauli.c`) that:
1. Builds the full 2^n × 2^n Pauli matrix via `build_pauli_matrix` (internal helper).
2. Computes the reference result `ref = M * psi` via `zmv` (CBLAS).
3. Constructs a `state_t` (ownership-transfer via `state_init`).
4. Calls `pauli_apply_pure` in-place.
5. Asserts element-wise equality with `TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL`.
6. Follows the goto-cleanup pattern for leak-free failure handling.

## Test Naming Convention

Unity test functions follow the pattern `test_<Gate>_<repr>` for gate tests,
`test_<operation>_<repr>` for state tests, and `test_pauli_<description>` for circuit/Pauli
tests. Harness driver functions (not directly registered with `RUN_TEST`) use upper-camel-case:
`testSingleQubitGate`, `testTwoQubitGateTiled`, `testTwoFromMat`, etc.

## Test Inventory

### `test_state` — 41 tests

| Group | Count | File |
|-------|------:|------|
| PURE state | 10 | `state/test_state_pure.c` |
| MIXED_PACKED state | 13 | `state/test_state_packed.c` |
| MIXED_TILED state | 16 | `state/test_state_tiled.c` |
| Dispatch layer | 2 | `state/test_state.c` |
| **Total** | **41** | |

### `test_gate` — 70 tests

| Group | Count | Representation | Gates / harness |
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
| Tiled `two_from_mat` | 8 | tiled | CX, SWAP, CZ, iSWAP, H⊗H (single); CX double, iSWAP double, iSWAP symmetry |
| **Total** | **70** | | |

### `test_circuit` — 30 tests

| Group | Count | Description |
|-------|------:|-------------|
| Pauli construction / masks | 5 | `test_pauli_masks_X/Y/Z`, `test_pauli_masks_str`, `test_pauli_from_str_negative` |
| `pauli_phase` | 1 | Spot-checks for X, Y, Z on col=0,1 |
| Single-qubit application (n=1) | 4 | I, X, Y, Z on 6 standard states each |
| Two-qubit strings (n=2) | 7 | XX, ZZ, XZ, YY, IZ, IX, IY |
| Three-qubit strings (n=3) | 6 | XYZ, ZZZ, IXI, XII, IZI, ZIZ |
| Four-qubit strings (n=4) | 4 | XZYX, IIXI (split_bit=1), IXII (split_bit=2), XIII (split_bit=3) |
| Round-trip and involution | 2 | `pauli_from_str` round-trip, P²=I involution |
| Identity edge case | 1 | III on uniform superposition |
| **Total** | **30** | |

## Build Type Requirement

Tests **must** be built and run with `CMAKE_BUILD_TYPE=Debug`. The project sets `LOG_TILE_DIM`
based on build type:

| Build type | `LOG_TILE_DIM` | Tile size | Effect on tests |
|---|---|---|---|
| **Debug** | 3 | 8×8 | Multi-tile code paths exercised at small qubit counts |
| **Release** | 5 | 32×32 | Multi-tile paths silently skipped — false passes |

## Running

```bash
cmake --preset debug          # Configure (enforces Debug build type)
cmake --build --preset debug  # Build
ctest --preset debug          # Run all tests

# Or run executables directly after building:
./cmake-build-debug/test/test_state
./cmake-build-debug/test/test_gate
./cmake-build-debug/test/test_circuit
```

## TODO

- Add tests for the Measurement module once it is implemented.
- Evaluate adding per-test-function CTest labels to enable `ctest -R` filtering by gate or
  representation without rebuilding.
- Document `MAXQUBITS` (currently 4) and `MAXQUBITS_TILED` (currently 6) rationale — why these
  bounds are sufficient to cover multi-tile paths at `LOG_TILE_DIM=3`.

## Open Questions

3. **Fixture completeness**: The state fixtures (cb, xb, yb, Bell, GHZ, W, maximally-mixed,
   random mixture) are listed but the rationale for this set is not documented. Are there known
   coverage gaps (e.g. no product states with relative phases, no separable mixed states)?

4. **`MAXQUBITS` bound justification**: `test.h` sets `MAXQUBITS=4`. It is unstated whether this
   is a memory constraint, a test-speed constraint, or a deliberate coverage choice, and whether
   it is sufficient to catch all off-by-one errors in qubit-index handling.

---

**Resolved questions** (retained for audit trail):

1. **BLAS dependency** (resolved): `utility/linalg.c` selects the CBLAS header by platform
   preprocessor guard: `<vecLib/cblas_new.h>` on Apple, `<cblas.h>` on Linux. There is no
   separate CMake option; the BLAS library is linked unconditionally via the platform-selected
   header. `test_state` does not link BLAS (`linalg.c` is not compiled into it).

2. **CTest registration** (resolved): `test_state`, `test_gate`, and `test_circuit` are each
   registered as a single CTest entry (`add_test(NAME … COMMAND …)`). No Unity CMake integration
   for per-function splitting exists. All Unity tests within an executable run as one CTest
   pass/fail.

5. **Circuit module tests** (resolved): A `test_circuit` executable exists at
   `test/circuit/test_circuit.c` (runner) + `test/circuit/test_pauli.c` (implementations).
   It is fully registered in `test/CMakeLists.txt` and adds 30 tests covering `pauli_apply_pure`
   and associated construction/phase/mask API.

6. **Test isolation** (resolved): `setUp()` and `tearDown()` are empty stubs in both
   `test_state.c` and `test_gate.c`. Fixture states are constructed and freed within each
   individual test function body. There is no shared mutable state across `RUN_TEST` calls;
   order-dependent failures are not a risk from fixture sharing.
