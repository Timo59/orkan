# Gate Module

Quantum gate operations for pure and mixed states across three backends.

---

## Directory Layout

```
include/
  gate.h                         Public API: 18 gate declarations

src/gate/
  gate.c                         Dispatch layer with input validation
  gate_pure.c                    Pure-state implementations (statevector)
  packed/
    gate_packed_1q.c             All 13 single-qubit packed implementations
    gate_packed_cx.c             CNOT packed
    gate_packed_cy.c             CY packed
    gate_packed_cz.c             CZ packed
    gate_packed_swap.c           SWAP packed
    gate_packed_ccx.c            Toffoli packed
  tiled/
    gate_tiled.h                 Shared declarations, OMP threshold, function prototypes
    gate_tiled_1q.c              single_from_mat() — generic 2x2 unitary helper
    gate_tiled_2q.c              two_from_mat() — generic 4x4 unitary helper
    gate_tiled_{x,y,z,...}.c     13 single-qubit tiled implementations (one per file)
    gate_tiled_{cx,cy,cz,...}.c  4 two-qubit + 1 three-qubit tiled implementations

src/internal/
  index.h                        Shared index helpers: insertBit0, pack_idx, tile_off, elem_off

test/gate/
  test_gate.c                    Test runner; named per-gate wrappers; RUN_TEST registrations
  test_gate_pure_{1q,2q,3q}.c    Pure-state harnesses
  test_gate_packed_{1q,2q,3q}.c  Packed harnesses
  test_gate_tiled_{1q,2q,2q_mat,3q}.c  Tiled harnesses
```

---

## Public API

Declared in `include/gate.h`. Every public function validates inputs (null pointer, qubit range, duplicate qubits) before dispatching.

### Single-Qubit Gates

| Gate | Signature |
|------|-----------|
| Pauli-X | `x(state, target)` |
| Pauli-Y | `y(state, target)` |
| Pauli-Z | `z(state, target)` |
| Hadamard | `h(state, target)` |
| Hadamard-Y | `hy(state, target)` |
| Phase S | `s(state, target)` |
| S-dagger | `sdg(state, target)` |
| T | `t(state, target)` |
| T-dagger | `tdg(state, target)` |
| P(theta) | `p(state, target, theta)` |
| Rx(theta) | `rx(state, target, theta)` |
| Ry(theta) | `ry(state, target, theta)` |
| Rz(theta) | `rz(state, target, theta)` |

### Two-Qubit Gates

| Gate | Signature |
|------|-----------|
| CNOT | `cx(state, control, target)` |
| CY | `cy(state, control, target)` |
| CZ | `cz(state, control, target)` |
| SWAP | `swap_gate(state, q1, q2)` |

`swap_gate` avoids a name clash with the C standard library.

### Three-Qubit Gates

| Gate | Signature |
|------|-----------|
| Toffoli | `ccx(state, ctrl1, ctrl2, target)` |

Input validation uses `GATE_VALIDATE` (defined in `gate.c`), which prints to stderr and calls `exit()` on failure. This is a deliberate fail-fast choice: silent corruption is worse than crashing.

---

## Architecture

### Dispatch

Every public gate function in `gate.c` validates inputs, then dispatches on `state->type`:

```c
switch (state->type) {
    case PURE:         gate_pure(state, ...);   break;
    case MIXED_PACKED: gate_packed(state, ...);  break;
    case MIXED_TILED:  gate_tiled(state, ...);   break;
}
```

Dispatch macros: `DISPATCH_1Q` (9 gates), `DISPATCH_ROT` (4 rotation gates), `DISPATCH_2Q` (3 controlled gates). `swap_gate` and `ccx` are dispatched inline due to non-standard naming or arity.

### Pure-State Backend

Operates on the statevector. Single-qubit gates process amplitude pairs `(psi[i], psi[i + 2^target])` via a traversal macro. Multi-qubit gates enumerate base indices with `insertBit0` / `insertBits2_0`.

### Mixed-Packed Backend

Implements rho -> U rho U-dagger directly on lower-triangular packed storage. Single-qubit gates iterate block pairs; two-qubit gates work on 4x4 sub-blocks; three-qubit gates on 8x8 sub-blocks. Operation classes:

- **Zero-flop** (X, CX, SWAP): amplitude swaps and conjugations only
- **Negation permutation** (Y, CY): swaps with sign negation
- **Diagonal** (Z, S, Sdg, T, Tdg, P, Rz, CZ): only off-diagonal elements change
- **Mixing** (H, Hy, Rx, Ry): all block elements participate

### Mixed-Tiled Backend

Operates on lower-triangular tile storage. Single-qubit gates split into within-tile (target < LOG_TILE_DIM) and cross-tile (target >= LOG_TILE_DIM) cases. Multi-qubit gates have additional cases depending on whether each qubit is tile-local or spans tile boundaries.

For density matrices, Rz and P differ only by a global phase that cancels in U rho U-dagger. The tiled backend implements `rz_tiled` as the primary and `p_tiled` delegates to it; the packed backend reverses this. In the pure backend both are distinct (global phase is observable).

### OpenMP

Each backend defines its own threshold below which OpenMP overhead exceeds computation cost:

| Backend | Threshold | Rationale |
|---------|-----------|-----------|
| Pure | 4096 | O(dim) work per gate |
| Packed | 512 | O(dim^2) work per gate |
| Tiled | 512 | O(dim^2) work per gate |

In Debug builds, CMake overrides `OMP_THRESHOLD=4` so multi-threaded paths are exercised at small qubit counts.

---

## Build Integration

Gate sources are compiled into `${PROJECT_NAME_LOWER}` via `src/CMakeLists.txt`. Packed and tiled sources are collected by glob; `gate.c` and `gate_pure.c` are listed explicitly.

Test target `test_gate` in `test/CMakeLists.txt` links the library, BLAS (for reference computations in `test/utility/linalg.c`), and Unity v2.6.1.

Installed header: `gate.h` → `${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME_LOWER}/`.

```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug -R test_gate
```

---

## Tests

Test framework, build-type requirement, and shared fixtures: see `docs/TEST_SUITE.md`. Gate matrices are defined in `test/include/test_gate.h`; reference matrix builders are in `test/utility/gatemat.c`.

### Test Architecture

Two-layer harness pattern:

- **`test_gate.c`** — single entry point. Contains `main()`, all named `test_<Gate>_<representation>()` wrappers (e.g. `test_X_pure`, `test_CX_tiled`), and all `RUN_TEST()` registrations.
- **`test_gate_<repr>_<Nq>.c`** — one harness per (state representation, gate arity). Each implements a parameterized driver (e.g. `testTwoQubitGateTiled`) that applies a gate, builds the reference full-matrix via `gatemat.c`, and asserts element-wise equality across all fixture states.
- **`test_gate_tiled_2q_mat.c`** — separate harness for the internal `two_from_mat` kernel, which accepts an arbitrary 4x4 matrix rather than a named gate. Implements `testTwoFromMat`, `testTwoFromMatDouble` (catches stale mirror-slot bugs), `testTwoFromMatSymmetry` (argument-order invariance). Defines `ISWAPMAT` and `HHMAT` constants used via `extern` in `test_gate.h`. Forward-declares `two_from_mat`, bypassing the public gate header.

Rotation matrix builders (`mat_rx`, `mat_ry`, `mat_rz`, `mat_p`, `mat_rx_pi3`) live in `utility/gatemat.c` (declared in `gatemat.h`). They fill a 4-element column-major `cplx_t[4]` with the 2x2 gate matrix for a given angle, shared by pure, packed, and tiled 1Q rotation harnesses.

Each harness sweeps all qubit positions across multiple system sizes and input states. Rotation gates are tested at multiple angles via `TEST_THETAS`. Correctness verified against a BLAS-based reference (`zgemv` for pure, `zgemm` for mixed) at tolerance 1e-12.

### Per-gate sweep tests

Each registered gate is verified independently on each storage representation. All three tables apply the same coverage strategy described above.

| Representation | Tests |
|---|---|
| PURE (`test_gate_pure_{1q,2q,3q}.c`) | `test_X_pure`, `test_Y_pure`, `test_Z_pure`, `test_H_pure`, `test_Hy_pure`, `test_S_pure`, `test_Sdg_pure`, `test_T_pure`, `test_Tdg_pure`, `test_Rx_pure`, `test_Ry_pure`, `test_Rz_pure`, `test_P_pure`, `test_CX_pure`, `test_CY_pure`, `test_CZ_pure`, `test_SWAP_pure`, `test_CCX_pure` |
| MIXED_PACKED (`test_gate_packed_{1q,2q,3q}.c`) | `test_X_packed`, `test_Y_packed`, `test_Z_packed`, `test_H_packed`, `test_Hy_packed`, `test_S_packed`, `test_Sdg_packed`, `test_T_packed`, `test_Tdg_packed`, `test_Rx_packed`, `test_Ry_packed`, `test_Rz_packed`, `test_P_packed`, `test_CX_packed`, `test_CY_packed`, `test_CZ_packed`, `test_SWAP_packed`, `test_CCX_packed` |
| MIXED_TILED (`test_gate_tiled_{1q,2q,3q}.c`) | `test_X_tiled`, `test_Y_tiled`, `test_Z_tiled`, `test_H_tiled`, `test_Hy_tiled`, `test_S_tiled`, `test_Sdg_tiled`, `test_T_tiled`, `test_Tdg_tiled`, `test_Rx_tiled`, `test_Ry_tiled`, `test_Rz_tiled`, `test_P_tiled`, `test_CX_tiled`, `test_CY_tiled`, `test_CZ_tiled`, `test_SWAP_tiled`, `test_CCX_tiled` |

### Internal `single_from_mat` kernel (`test_gate_tiled_1q.c`)

Verifies the 1Q tiled kernel on arbitrary 2x2 unitaries.

| Test | Matrix |
|---|---|
| `test_single_from_mat_id_tiled` | Identity |
| `test_single_from_mat_x_tiled` | Pauli-X |
| `test_single_from_mat_y_tiled` | Pauli-Y |
| `test_single_from_mat_h_tiled` | Hadamard |
| `test_single_from_mat_t_tiled` | T gate |
| `test_single_from_mat_rx_tiled` | Rx(π/3) |
| `test_single_from_mat_xh_tiled` | X · H composition |
| `test_single_from_mat_ht_tiled` | H · T composition |

### Internal `two_from_mat` kernel (`test_gate_tiled_2q_mat.c`)

Verifies the 2Q tiled kernel on arbitrary 4x4 unitaries.

| Test | Verifies |
|---|---|
| `test_two_from_mat_cx` | CNOT applied via `two_from_mat` matches named `cx` |
| `test_two_from_mat_swap` | SWAP applied via `two_from_mat` |
| `test_two_from_mat_cz` | CZ applied via `two_from_mat` |
| `test_two_from_mat_iswap` | iSWAP (not in the named-gate set) |
| `test_two_from_mat_hh` | H ⊗ H |
| `test_two_from_mat_cx_double` | CX applied twice — catches stale mirror-slot state between applications |
| `test_two_from_mat_iswap_double` | iSWAP applied twice |
| `test_two_from_mat_symmetry` | Argument-order invariance for symmetric matrices |
