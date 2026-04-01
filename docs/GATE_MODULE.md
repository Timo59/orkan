# Gate Module

Quantum gate operations for pure and mixed states across three backends.

---

## Overview

The gate module implements 18 quantum gates, each with three backends: pure state (statevector), packed mixed state (density matrix), and tiled mixed state (blocked density matrix). All public functions dispatch at runtime on `state->type`.

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
  test_gate.c                    Test runner (56 tests)
  test_gate_pure_{1q,2q,3q}.c   Pure-state test harnesses
  test_gate_packed_{1q,2q,3q}.c  Packed test harnesses
  test_gate_tiled_{1q,2q,2q_mat,3q}.c  Tiled test harnesses
```

---

## Public API

Declared in `include/gate.h`. Every function validates inputs (null pointer, qubit range, duplicate qubits) before dispatching.

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

---

## Dispatch Architecture

Every public gate function in `gate.c` validates inputs, then dispatches on `state->type`:

```c
switch (state->type) {
    case PURE:         gate_pure(state, ...);   break;
    case MIXED_PACKED: gate_packed(state, ...);  break;
    case MIXED_TILED:  gate_tiled(state, ...);   break;
}
```

Dispatch macros: `DISPATCH_1Q` (9 gates), `DISPATCH_ROT` (4 rotation gates), `DISPATCH_2Q` (3 controlled gates). `swap_gate` and `ccx` are dispatched inline due to non-standard naming or arity.

Input validation uses `GATE_VALIDATE` (defined in `gate.c`), which prints to stderr and calls `exit()` on failure. This is a deliberate fail-fast choice for a scientific computing library where silent corruption is worse than crashing.

---

## Backend Algorithms

### Pure State

Operates on the statevector. Single-qubit gates process amplitude pairs `(psi[i], psi[i + 2^target])` via a traversal macro. Multi-qubit gates enumerate base indices with `insertBit0` / `insertBits2_0`.

### Mixed Packed

Implements rho -> U rho U-dagger directly on lower-triangular packed storage. Single-qubit gates iterate block pairs; two-qubit gates work on 4x4 sub-blocks; three-qubit gates on 8x8 sub-blocks. Gate operation classes:

- **Zero-flop** (X, CX, SWAP): amplitude swaps and conjugations only
- **Negation permutation** (Y, CY): swaps with sign negation
- **Diagonal** (Z, S, Sdg, T, Tdg, P, Rz, CZ): only off-diagonal elements change
- **Mixing** (H, Hy, Rx, Ry): all block elements participate

### Mixed Tiled

Operates on lower-triangular tile storage. Single-qubit gates split into within-tile (target < LOG_TILE_DIM) and cross-tile (target >= LOG_TILE_DIM) cases. Multi-qubit gates have additional cases depending on whether each qubit is tile-local or spans tile boundaries.

For density matrices, Rz and P differ only by a global phase that cancels in U rho U-dagger. The tiled backend implements `rz_tiled` as the primary and `p_tiled` delegates to it; the packed backend reverses this. In the pure backend both are distinct (global phase is observable).

---

## OpenMP Parallelization

Each backend defines its own `OMP_THRESHOLD` — the minimum dimension below which OpenMP overhead exceeds the computation cost:

| Backend | Threshold | Rationale |
|---------|-----------|-----------|
| Pure | 4096 | O(dim) work per gate |
| Packed | 512 | O(dim^2) work per gate |
| Tiled | 512 | O(dim^2) work per gate |

In Debug builds, CMake overrides `OMP_THRESHOLD=4` so multi-threaded paths are exercised at small qubit counts.

---

## Build Integration

Gate sources are compiled into the shared library target `${PROJECT_NAME_LOWER}`, defined in `src/CMakeLists.txt`. Packed and tiled sources are collected by glob; `gate.c` and `gate_pure.c` are listed explicitly.

Test target: `test_gate`, defined in `test/CMakeLists.txt`. Links the library, BLAS (for reference computations in `test/utility/linalg.c`), and Unity v2.6.1.

```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug --tests-regex test_gate
```

---

## Test Infrastructure

Each backend has dedicated test harnesses that verify all qubit positions across multiple system sizes and input states. Rotation gates are tested at multiple angles. Correctness is verified against BLAS-based reference computation (`zgemv` for pure, `zgemm` for mixed) at tolerance 1e-12.

| Harness | File | Backend |
|---------|------|---------|
| `testSingleQubitGate` | `test_gate_pure_1q.c` | PURE |
| `testRotationGate` | `test_gate_pure_1q.c` | PURE |
| `testTwoQubitGate` | `test_gate_pure_2q.c` | PURE |
| `testThreeQubitGate` | `test_gate_pure_3q.c` | PURE |
| `testSingleQubitGateMixed` | `test_gate_packed_1q.c` | MIXED_PACKED |
| `testRotationGateMixed` | `test_gate_packed_1q.c` | MIXED_PACKED |
| `testTwoQubitGateMixed` | `test_gate_packed_2q.c` | MIXED_PACKED |
| `testThreeQubitGateMixed` | `test_gate_packed_3q.c` | MIXED_PACKED |
| `testSingleQubitGateTiled` | `test_gate_tiled_1q.c` | MIXED_TILED |
| `testRotationGateTiled` | `test_gate_tiled_1q.c` | MIXED_TILED |
| `testTwoQubitGateTiled` | `test_gate_tiled_2q.c` | MIXED_TILED |
| `testTwoQubitGateTiledMat` | `test_gate_tiled_2q_mat.c` | MIXED_TILED |
| `testThreeQubitGateTiled` | `test_gate_tiled_3q.c` | MIXED_TILED |

Gate matrices are defined in `test/include/test_gate.h`. Reference matrix builders are in `test/utility/gatemat.c`.
