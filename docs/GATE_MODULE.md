# Gate Module: Technical Specification

**Module:** Quantum Gate Operations
**Status:** Nearly complete — CCX tiled implementation exists but its test harness is unimplemented (see TODO)
**Last verified build:** 51 gate test functions, 0 failures (tolerance 1e-12)
**Last updated:** 2026-02-27 (full source audit)

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
  gate.c                         # Central dispatch: DISPATCH_1Q, DISPATCH_ROT, DISPATCH_2Q,
                                 #   DISPATCH_2Q_NO_TILED (defined but unused in current code);
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
                                 #   OMP_THRESHOLD (512); also declares all tiled gate prototypes.
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
    gate_tiled_cy.c              # CY tiled — full implementation (4 cases, same structure as cx)
    gate_tiled_cz.c              # CZ tiled — full implementation (4 cases, negate-only)
    gate_tiled_swap.c            # SWAP tiled
    gate_tiled_ccx.c             # CCX tiled — full implementation (dispatched from gate.c)
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

## Dependencies (Inputs)

| Dependency | Source | Purpose |
|---|---|---|
| `state_t` and `StateType` enum | `include/state.h` | State representation consumed by every gate |
| `qlib.h` | `include/qlib.h` | Common macros, error handling primitives |
| `gate_tiled.h` | `src/gate/tiled/gate_tiled.h` | Element-access helpers used by all tiled implementations |
| BLAS (`zgemv`, `zgemm`) | Platform BLAS / OpenBLAS | Test reference computations only |
| OpenMP | System / compiler | Optional parallelization in packed and tiled backends |

---

## Outputs / Artifacts

| Artifact | Description |
|---|---|
| Object files in `src/gate/` | Compiled into the main `qsim` library target |
| `test/test_gate` binary | CTest executable: 51 test functions |

The gate module is a consumer of `state_t`; it does not define any persistent data structures of its own. All gate functions operate in-place on the state passed by pointer.

---

## Build Integration

The gate module is compiled as part of `src/CMakeLists.txt` as part of the `q` shared library. There is no separate `CMakeLists.txt` inside `src/gate/`. Key points:

- Packed gate sources are collected via `file(GLOB GATE_PACKED_SOURCES "gate/packed/gate_packed_*.c")`.
- Tiled gate sources are collected via `file(GLOB GATE_TILED_SOURCES "gate/tiled/gate_tiled_*.c")`. This glob matches only files ending in `.c`, so `gate_tiled_cx.c.old` and `gate_tiled_swap.c.old` are excluded from compilation.
- `gate.c` and `gate_pure.c` are listed explicitly.
- Test target `test_gate` is defined under `test/` and linked against the library and BLAS.
- `LOG_TILE_DIM` is set to `3` (Debug) or `5` (Release) by the build system. Do not override manually.
- In Debug builds, `OMP_THRESHOLD=4` is injected as a compile definition, overriding the per-file `#ifndef OMP_THRESHOLD` defaults. This forces OpenMP to activate at tiny state sizes during testing.

**Build gate tests only:**
```bash
cmake --preset debug
cmake --build --preset debug --target test_gate
./cmake-build-debug/test/test_gate
```

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

**Dispatch macros in `gate.c`:** `DISPATCH_1Q`, `DISPATCH_ROT`, `DISPATCH_2Q`, `DISPATCH_2Q_NO_TILED`. All 13 single-qubit gates use `DISPATCH_1Q` or `DISPATCH_ROT`. All four two-qubit gates use `DISPATCH_2Q`. There is no `DISPATCH_3Q`; CCX and `swap_gate` are dispatched inline in `gate.c`.

`DISPATCH_2Q_NO_TILED` is defined in `gate.c` but is not currently used by any gate — it was intended for CY/CZ before their tiled implementations were written.

**`swap_gate` inline dispatch rationale:** `DISPATCH_2Q(name, ...)` generates calls to `name_pure`, `name_packed`, `name_tiled`. The internal implementations are named `swap_pure/packed/tiled` (not `swap_gate_pure/packed/tiled`), so the macro cannot be used with the public name `swap_gate`. The inline dispatch calls the correct internal names directly.

---

## Algorithm Summaries

### Pure State (PURE)

**Single-qubit** — `TRAVERSE_PURE_1Q(state, target, PAIR_OP)`:
Processes pairs `(ψ[i], ψ[i + 2^target])` where the target bit differs. Nested loop: outer stride `2×2^target`, inner `2^target`. OpenMP: `if(dim >= OMP_THRESHOLD)` where `OMP_THRESHOLD=4096` by default (overridden to `4` in Debug builds by CMake).

**Two-qubit** — flat loop with `insertBits2_0(k, lo, hi)`:
Enumerates `dim/4` base indices with both qubit bits = 0. No traversal macro; each gate inlines the loop. SWAP swaps the `(01, 10)` amplitude pair; CX/CY/CZ swap/phase the `(ctrl=1,tgt=0)` and `(ctrl=1,tgt=1)` amplitudes.

**Three-qubit** — flat loop with chained `insertBit0()`:
Enumerates `dim/8` base indices with all three qubit bits = 0. CCX swaps amplitudes at `(c1=1, c2=1, tgt=0)` and `(c1=1, c2=1, tgt=1)`.

### Mixed Packed (MIXED_PACKED)

Implements ρ → UρU† directly on lower-triangular packed storage. No unpacking required.

**Single-qubit** — `TRAVERSE_PACKED_BLOCKS(state, target, BLOCK_OP)`:
Iterates lower-triangle pairs `(bc ≤ br)` of column/row base indices, each producing 4 packed indices `(idx00, idx01, idx10, idx11)`. `lower_01` and `diag_block` flags control conjugation and skipping. OpenMP `schedule(static)`.

**Two-qubit** — `insertBits2_0()` with 4×4 block pairs:
Two nested loops `(bc ≤ br)`, each producing 4 row and 4 column indices — a 4×4 sub-block of the density matrix. Fast path when `r00 > c11` (all in lower triangle); slow path for upper-triangle crossings. OpenMP `schedule(static, 1)`.

**Three-qubit** — same pattern with `insertBit0()` chained three times, producing 8×8 blocks. OpenMP `schedule(static, 1)`.

**Per-gate operation types in packed:**
- **Zero-flop** (X, CX, SWAP): only swap pairs and conjugations — no multiply/add arithmetic.
- **Negation permutation** (Y, CY): swaps with sign negation. The ±i entries of Y/CY and their adjoints combine to −1 in UρU†, so no complex multiplication — only negation and conjugation.
- **Diagonal** (Z, S, Sdg, T, Tdg, P, Rz, CZ): only off-diagonal elements change; diagonal elements unchanged. Guard `idx01` write with `if (!diag_block)` to avoid double-processing the aliased element.
- **Mixing** (H, Hy, Rx, Ry): all 4 elements of each 2×2 block participate; Hermitian symmetry enforced by setting `a01 = conj(a10)` on diagonal blocks.

### Mixed Tiled (MIXED_TILED)

`TILE_DIM` and `LOG_TILE_DIM` are build-preset-controlled: `LOG_TILE_DIM=3` (TILE_DIM=8) in Debug, `LOG_TILE_DIM=5` (TILE_DIM=32) in Release. Only lower-triangular tiles (tr ≥ tc) are stored. Within diagonal tiles, only lr ≥ lc elements are populated.

**Single-qubit** — two traversal modes:

1. **Within-tile** (`target < LOG_TILE_DIM`): all 4 butterfly elements are in the same tile. Diagonal tiles require special handling — the (0,1) element at `(lr0, lc1)` may be in the upper triangle and must be read as `conj(elem[lc1, lr0])`.

2. **Cross-tile** (`target ≥ LOG_TILE_DIM`): butterfly spans 4 tiles T00/T01/T10/T11, grouped by tile-level bit `tile_bit = target - LOG_TILE_DIM`. T(tr0,tc1) may be an upper-triangular tile — access as conjugate of T(tc1,tr0).

**No shared traversal macros exist for tiled gates.** Each per-gate `.c` file inlines the full traversal logic independently for both the within-tile and cross-tile paths, including all triangle sub-cases. `gate_tiled.h` provides the four element-access helpers (`tile_off`, `elem_off`, `dtile_read`, `dtile_write`) and `OMP_THRESHOLD`, plus all tiled gate function declarations. See `gate_tiled_x.c` as the simplest reference implementation and `gate_tiled_cx.c` for the two-qubit cross-tile pattern.

**P/Rz delegation:** The global phase from Rz cancels in UρU†, making Rz and P identical for density matrices, so one delegates to the other. The delegation direction differs by backend:
- **Packed:** `p_packed()` is the full implementation; `rz_packed()` delegates to it.
- **Tiled:** `rz_tiled()` is the full implementation; `p_tiled()` delegates to it.
- **Pure:** `p_pure()` is distinct — pure states retain the global phase difference.

**Key pitfall — diag_block aliasing:** When `lr0 == lc0` (within-tile diagonal), `a01` and `a10` alias the same memory location. Diagonal gates must guard with `if (!diag_block)`. Mixing gates must enforce Hermitian symmetry: after computing `a10`, set `a01 = conj(a10)`.

**Two-qubit tiled:** Same four-block decomposition as packed, with an additional case split based on qubit position relative to `LOG_TILE_DIM`. Each of `cx_tiled`, `cy_tiled`, `cz_tiled`, `swap_tiled`, and `ccx_tiled` follows the same 4-case structure (both-local / target-tile-level / control-tile-level / both-tile-level). Triangle status (lower/upper/diagonal) is determined once per tile group before processing elements. See `gate_tiled_cx.c` as the reference implementation for the 2Q cross-tile pattern; `gate_tiled_cy.c` adds per-element phase factors; `gate_tiled_cz.c` is negate-only (no permutation).

**OMP_THRESHOLD:** `512` by default for packed and tiled backends (defined via `#ifndef` guards in `gate_packed_1q.c` and `gate_tiled.h` respectively; each file defines independently). Pure state backend uses `4096` by default (defined in `gate_pure.c`). In Debug builds, CMake injects `OMP_THRESHOLD=4` as a compile definition, overriding all per-file defaults for the `q` library target.

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

### Adding the Missing Tiled Tests

The tiled section of `test_gate.c` currently has 15 `RUN_TEST` entries (all 13 single-qubit gates + CX + SWAP). Missing: `test_CY_tiled`, `test_CZ_tiled`, `test_CCX_tiled`. The implementations exist; only the test wiring is absent.

For `ccx_tiled`, first implement `testThreeQubitGateTiled()` in `test/gate/test_gate_tiled_2q.c`, then in `test/gate/test_gate.c`:

```c
void test_CY_tiled(void)  { testTwoQubitGateTiled(cy, CYMAT); }
void test_CZ_tiled(void)  { testTwoQubitGateTiled(cz, CZMAT); }
void test_CCX_tiled(void) { testThreeQubitGateTiled(ccx, CCXMAT); }

// In main():
RUN_TEST(test_CY_tiled);
RUN_TEST(test_CZ_tiled);
RUN_TEST(test_CCX_tiled);
```

Gate matrix constants (`CCXMAT`, `CXMAT`, `CYMAT`, `CZMAT`, `SWAPMAT`, etc.) are in `test/include/test_gate.h`. Reference matrix builders (`mat_three_qubit_gate`, `mat_two_qubit_gate`, `mat_single_qubit_gate`) are in `test/utility/gatemat.c`.

### Test Coverage Design

- **System sizes:** 1Q gates: n=1,2,3,4; 2Q gates: n=2,3,4; 3Q gates: n=3,4
- **All qubit positions** tested (all targets, all ctrl/tgt pairs, all triples)
- **Test states (pure):** All 2^n Z-eigenbasis, X-eigenbasis, and Y-eigenbasis states (3×2^n vectors total); Bell states (n=2) or GHZ+W states (n≥3)
- **Test states (mixed):** Same as pure converted to density matrices, plus maximally mixed state (I/2^n) and 10 random mixed states (fixed seed: 1000 + call_count)
- **Rotation angles:** θ ∈ {0, π/4, π/2, π, 3π/2, −π/3}
- **Tiled edge case:** `MAXQUBITS_TILED=6` to exercise cross-tile paths (n=5 is single-tile boundary, n=6 requires multiple tiles)
- **Tolerance:** ε = 1e-12
- **Reference:** pure states via `zgemv`; mixed states via `zgemm` (UρU†)

---

## Known Issues / Deferred

- **RZZ gate:** Diagonal two-qubit gate useful for QAOA. Deferred to a future revision.
- **Higher-order controlled gates** (CS, CT, CH, etc.): Out of scope for this module.

---

## Open Questions

5. **`swap_gate` inline dispatch rationale — RESOLVED.** `DISPATCH_2Q(name, ...)` expands to calls to `name_pure`, `name_packed`, `name_tiled`. The internal implementations are named `swap_pure/packed/tiled` (not `swap_gate_pure/packed/tiled`), so the macro cannot bridge the public name `swap_gate` to the internal `swap` prefix. Inline dispatch is the correct final pattern.

6. **OpenMP threshold unification — RESOLVED.** Each translation unit defines `OMP_THRESHOLD` independently via an `#ifndef` guard: `gate_pure.c` defaults to `4096`; `gate_packed_1q.c`, `gate_packed_cx.c` (and all other packed files), and `gate_tiled.h` default to `512`. The values differ by design (pure state work is O(dim); packed/tiled work is O(dim²), justifying a lower threshold). In Debug builds, CMake overrides all three to `4` via a compile definition that takes precedence over the `#ifndef` guards. No silent divergence risk.
