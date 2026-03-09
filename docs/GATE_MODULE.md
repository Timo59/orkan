# Gate Module: Technical Specification

**Module:** Quantum Gate Operations
**Status:** Nearly complete — CCX tiled implementation exists but its test harness is unimplemented (see TODO)
**Last verified build:** 51 gate test functions, 0 failures (tolerance 1e-12)
**Last updated:** 2026-03-08 (spec cleanup)

---

## TODO

- [ ] Implement `testThreeQubitGateTiled()` in `test/gate/test_gate_tiled_2q.c` (declared, not implemented)
- [ ] Add `RUN_TEST(test_CCX_tiled)`, `RUN_TEST(test_CY_tiled)`, `RUN_TEST(test_CZ_tiled)` in `test/gate/test_gate.c`
- [ ] Delete stale files: `gate_tiled_cx.c.old`, `gate_tiled_swap.c.old`

---

## Directory Layout

```
include/
  gate.h                         # Public API: 18 gate declarations, GATE_VALIDATE macro,
                                 #   gate_idx_t typedef, insertBit0/insertBits2_0/pack_idx
                                 #   static inlines (shared across all backends)

src/gate/
  gate.c                         # Central dispatch: DISPATCH_1Q, DISPATCH_ROT, DISPATCH_2Q;
                                 #   CCX and swap_gate dispatched inline
  gate_pure.c                    # All pure-state implementations (13×1Q + 4×2Q + 1×3Q)
  packed/
    gate_packed_1q.c             # All 13 single-qubit packed implementations
    gate_packed_cx.c
    gate_packed_cy.c
    gate_packed_cz.c
    gate_packed_swap.c
    gate_packed_ccx.c
  tiled/
    gate_tiled.h                 # Shared inline helpers: tile_off, elem_off, dtile_read/write,
                                 #   OMP_THRESHOLD; also declares all tiled gate prototypes.
                                 #   No gate_tiled.c exists.
    gate_tiled_x.c               # ]
    gate_tiled_y.c               # ]
    gate_tiled_z.c               # ] 13 single-qubit tiled implementations
    gate_tiled_h.c               # ]  (one gate per file)
    gate_tiled_hy.c              # ]
    gate_tiled_s.c               # ]
    gate_tiled_sdg.c             # ]
    gate_tiled_t.c               # ]
    gate_tiled_tdg.c             # ]
    gate_tiled_rx.c              # ]
    gate_tiled_ry.c              # ]
    gate_tiled_rz.c              # ]
    gate_tiled_p.c               # ] (delegates to rz_tiled)
    gate_tiled_cx.c              # CX tiled — reference implementation for 2Q tiled pattern
    gate_tiled_cy.c              # CY tiled
    gate_tiled_cz.c              # CZ tiled
    gate_tiled_swap.c            # SWAP tiled
    gate_tiled_ccx.c             # CCX tiled
    gate_tiled_cx.c.old          # Stale — superseded; should be deleted (excluded from build)
    gate_tiled_swap.c.old        # Stale — superseded; should be deleted (excluded from build)

test/gate/
  test_gate.c                    # Test registry: function definitions + main()
                                 #   with 51 RUN_TESTs (18 pure + 18 packed + 15 tiled)
  test_gate_pure_1q.c            # Pure 1Q harnesses + rotation matrix builders
  test_gate_pure_2q.c            # Pure 2Q harness
  test_gate_pure_3q.c            # Pure 3Q harness
  test_gate_packed_1q.c          # Packed 1Q harnesses
  test_gate_packed_2q.c          # Packed 2Q harness
  test_gate_packed_3q.c          # Packed 3Q harness
  test_gate_tiled_1q.c           # Tiled 1Q harnesses
  test_gate_tiled_2q.c           # Tiled 2Q+3Q harness (testThreeQubitGateTiled declared, not impl)

test/include/
  test_gate.h                    # Typedefs, gate matrix constants, harness declarations

test/utility/
  gatemat.c                      # mat_single_qubit_gate(), mat_two_qubit_gate(),
                                 #   mat_three_qubit_gate()
```

---

## Dependencies

| Dependency | Source | Purpose |
|---|---|---|
| `state_t` and `StateType` enum | `include/state.h` | State representation consumed by every gate |
| `qlib.h` | `include/qlib.h` | Common macros, error handling primitives |
| `gate_tiled.h` | `src/gate/tiled/gate_tiled.h` | Element-access helpers used by all tiled implementations |
| OpenMP | System / compiler | Optional parallelization in packed and tiled backends |

---

## Outputs

The gate module compiles into the `libq` shared library. It does not define any persistent data structures of its own. All gate functions operate in-place on the state passed by pointer.

---

## Build Integration

The gate module is compiled as part of `src/CMakeLists.txt` into the `libq` shared library. There is no separate `CMakeLists.txt` inside `src/gate/`. Packed and tiled sources are collected by glob; `gate.c` and `gate_pure.c` are listed explicitly.

---

## Public Interface

All functions declared in `include/gate.h`. Every function validates its inputs with `GATE_VALIDATE` (fail-fast: prints to stderr and calls `exit(EXIT_FAILURE)`) before dispatching.

### Single-Qubit Gates (all complete)

| Gate | Signature | Pure | Packed | Tiled |
|---|---|---|---|---|
| Pauli-X | `x(state, tgt)` | ✓ | ✓ | ✓ |
| Pauli-Y | `y(state, tgt)` | ✓ | ✓ | ✓ |
| Pauli-Z | `z(state, tgt)` | ✓ | ✓ | ✓ |
| Hadamard | `h(state, tgt)` | ✓ | ✓ | ✓ |
| Hadamard-Y | `hy(state, tgt)` | ✓ | ✓ | ✓ |
| Phase S | `s(state, tgt)` | ✓ | ✓ | ✓ |
| S† | `sdg(state, tgt)` | ✓ | ✓ | ✓ |
| T | `t(state, tgt)` | ✓ | ✓ | ✓ |
| T† | `tdg(state, tgt)` | ✓ | ✓ | ✓ |
| P(θ) | `p(state, tgt, θ)` | ✓ | ✓ | ✓ |
| Rx(θ) | `rx(state, tgt, θ)` | ✓ | ✓ | ✓ |
| Ry(θ) | `ry(state, tgt, θ)` | ✓ | ✓ | ✓ |
| Rz(θ) | `rz(state, tgt, θ)` | ✓ | ✓ | ✓ |

### Two-Qubit Gates (all complete)

| Gate | Signature | Pure | Packed | Tiled |
|---|---|---|---|---|
| CNOT | `cx(state, ctrl, tgt)` | ✓ | ✓ | ✓ |
| CY | `cy(state, ctrl, tgt)` | ✓ | ✓ | ✓ |
| CZ | `cz(state, ctrl, tgt)` | ✓ | ✓ | ✓ |
| SWAP | `swap_gate(state, q1, q2)` | ✓ | ✓ | ✓ |

`swap` is named `swap_gate` to avoid a clash with the C standard library.

### Three-Qubit Gates

| Gate | Signature | Pure | Packed | Tiled |
|---|---|---|---|---|
| Toffoli | `ccx(state, c1, c2, tgt)` | ✓ | ✓ | ✓ |

All three backends are implemented. The tiled test harness (`testThreeQubitGateTiled`) and its `RUN_TEST` registration remain outstanding (see TODO).

---

## Dispatch Architecture

Every public gate validates inputs with `GATE_VALIDATE`, then dispatches on `state->type`:

```c
switch (state->type) {
    case PURE:         gate_pure(state, ...);   break;
    case MIXED_PACKED: gate_packed(state, ...);  break;
    case MIXED_TILED:  gate_tiled(state, ...);   break;
    default: GATE_VALIDATE(0, "gate: unknown state type");
}
```

Dispatch is handled via macros in `gate.c`: `DISPATCH_1Q` and `DISPATCH_ROT` for single-qubit gates, `DISPATCH_2Q` for two-qubit gates. CCX and `swap_gate` are dispatched inline — CCX because there is no `DISPATCH_3Q` macro, and `swap_gate` because the public name does not match the internal `swap_` prefix expected by `DISPATCH_2Q`.

---

## Algorithm Summaries

### Pure State (PURE)

Operates directly on the state vector ψ. Single-qubit gates process amplitude pairs `(ψ[i], ψ[i + 2^target])` via a traversal macro. Two- and three-qubit gates enumerate base indices with the relevant qubit bits set to 0 using `insertBits2_0` / chained `insertBit0`, then apply the gate to the corresponding amplitude tuple.

### Mixed Packed (MIXED_PACKED)

Implements ρ → UρU† directly on lower-triangular packed storage, without unpacking. Single-qubit gates iterate lower-triangle block pairs; two-qubit gates work on 4×4 sub-blocks; three-qubit gates on 8×8 sub-blocks. Gates fall into four operation classes:

- **Zero-flop** (X, CX, SWAP): amplitude swaps and conjugations only — no arithmetic.
- **Negation permutation** (Y, CY): swaps with sign negation; the ±i factors combine to −1 in UρU†.
- **Diagonal** (Z, S, Sdg, T, Tdg, P, Rz, CZ): only off-diagonal elements change.
- **Mixing** (H, Hy, Rx, Ry): all elements of each block participate; Hermitian symmetry is enforced on diagonal blocks.

### Mixed Tiled (MIXED_TILED)

Operates on lower-triangular tile storage (only tiles with tr ≥ tc are stored). Single-qubit gates split into two cases based on whether the target qubit falls within a tile (within-tile butterfly) or spans tile boundaries (cross-tile butterfly). Two-qubit gates have four cases depending on whether each qubit is local or at the tile level. The same structure extends to CCX. Triangle status (lower/upper/diagonal) is determined once per tile group.

**P/Rz:** For density matrices, Rz and P differ only by a global phase that cancels in UρU†. In the packed backend `p_packed` is the full implementation and `rz_packed` delegates to it; in the tiled backend the roles are reversed. In the pure backend, `p_pure` and `rz_pure` are distinct (global phase is physically observable).

---

## Test Infrastructure

### Test Harnesses

| Harness | File | State Type |
|---|---|---|
| `testSingleQubitGate(fn, mat)` | `test_gate_pure_1q.c` | PURE |
| `testRotationGate(fn, mat_fn)` | `test_gate_pure_1q.c` | PURE |
| `testTwoQubitGate(fn, mat)` | `test_gate_pure_2q.c` | PURE |
| `testThreeQubitGate(fn, mat)` | `test_gate_pure_3q.c` | PURE |
| `testSingleQubitGateMixed(fn, mat)` | `test_gate_packed_1q.c` | MIXED_PACKED |
| `testRotationGateMixed(fn, mat_fn)` | `test_gate_packed_1q.c` | MIXED_PACKED |
| `testTwoQubitGateMixed(fn, mat)` | `test_gate_packed_2q.c` | MIXED_PACKED |
| `testThreeQubitGateMixed(fn, mat)` | `test_gate_packed_3q.c` | MIXED_PACKED |
| `testSingleQubitGateTiled(fn, mat)` | `test_gate_tiled_1q.c` | MIXED_TILED |
| `testRotationGateTiled(fn, mat_fn)` | `test_gate_tiled_1q.c` | MIXED_TILED |
| `testTwoQubitGateTiled(fn, mat)` | `test_gate_tiled_2q.c` | MIXED_TILED |
| `testThreeQubitGateTiled(fn, mat)` | `test_gate_tiled_2q.c` | DECLARED in `test_gate.h`, NOT IMPLEMENTED |

Each harness tests all qubit positions across multiple system sizes, using a variety of input states and — for rotation gates — multiple angles. Correctness is verified against a BLAS-based reference (pure states via `zgemv`; mixed states via `zgemm`) at tolerance ε = 1e-12. Gate matrix constants are in `test/include/test_gate.h`; reference matrix builders are in `test/utility/gatemat.c`.

---

## Known Issues / Deferred

- **RZZ gate:** Diagonal two-qubit gate useful for QAOA. Deferred to a future revision.
- **Higher-order controlled gates** (CS, CT, CH, etc.): Out of scope for this module.
