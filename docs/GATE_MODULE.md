# Gate Module: Consolidated Technical Specification

**Module:** Quantum Gate Operations
**Status:** ⚠️ Nearly complete — one implementation remaining (see TODO)
**Last verified build:** 36 gate test functions, 0 failures (tolerance 1e-12)
**Last updated:** 2026-02-17

---

## TODO

**One remaining item:**

- [ ] Implement `ccx_tiled` in `src/gate/tiled/gate_tiled_ccx.c` (a `.new` stub exists)
- [ ] Implement `testThreeQubitGateTiled()` in `test/gate/test_gate_tiled_2q.c` (declared, not implemented)
- [ ] Add `RUN_TEST(test_CCX_tiled)` in `test/gate/test_gate.c`
- [ ] Wire up dispatch in `src/gate/gate.c` (`DISPATCH_3Q` or inline)

All other gates are implemented and tested for all three memory layouts (pure, packed, tiled).

---

## Module Layout

```
include/
  gate.h                         # Public API (18 gate declarations), GATE_VALIDATE macro

src/gate/
  gate.c                         # Three-way dispatch for all gates, extern declarations
  gate_pure.c                    # All pure implementations (13×1Q + 4×2Q + 1×3Q)
  packed/
    gate_packed_1q.c             # All 13 single-qubit packed implementations
    gate_packed_cx.c
    gate_packed_cy.c
    gate_packed_cz.c
    gate_packed_swap.c
    gate_packed_ccx.c
  tiled/
    gate_tiled.h                 # Shared inline helpers: tile_off, elem_off, dtile_read/write
    gate_tiled_x.c               # ]
    gate_tiled_y.c               # ]
    gate_tiled_z.c               # ] 13 single-qubit tiled implementations
    gate_tiled_h.c               # ]  (each file = one gate)
    gate_tiled_hy.c              # ]
    gate_tiled_s.c               # ]
    gate_tiled_sdg.c             # ]
    gate_tiled_t.c               # ]
    gate_tiled_tdg.c             # ]
    gate_tiled_rx.c              # ]
    gate_tiled_ry.c              # ]
    gate_tiled_rz.c              # ]
    gate_tiled_p.c               # ] (delegates to rz_tiled)
    gate_tiled_cx.c              # CX tiled
    gate_tiled_swap.c            # SWAP tiled
    gate_tiled_ccx.c.new         # ← WORK IN PROGRESS (rename to .c when done)

test/gate/
  test_gate.c                    # Test registry: function defs + main() with 36 RUN_TESTs
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
  gatemat.c                      # mat_single_qubit_gate(), mat_two_qubit_gate(), mat_three_qubit_gate()
```

**Build (gate tests only — full build is broken due to stale benchmark headers):**
```bash
cd cmake-build-debug
cmake --build . --target test_gate
./test/test_gate
```

---

## Gate Completion Status

### Single-Qubit Gates (all complete)

| Gate | API | Pure | Packed | Tiled |
|------|-----|------|--------|-------|
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

| Gate | API | Pure | Packed | Tiled |
|------|-----|------|--------|-------|
| CNOT | `cx(state, ctrl, tgt)` | ✓ | ✓ | ✓ |
| CY | `cy(state, ctrl, tgt)` | ✓ | ✓ | ✓ |
| CZ | `cz(state, ctrl, tgt)` | ✓ | ✓ | ✓ |
| SWAP | `swap_gate(state, q1, q2)` | ✓ | ✓ | ✓ |

Note: `swap` is named `swap_gate` in the public API to avoid the C standard library clash.

### Three-Qubit Gates

| Gate | API | Pure | Packed | Tiled |
|------|-----|------|--------|-------|
| Toffoli | `ccx(state, c1, c2, tgt)` | ✓ | ✓ | **missing** |

---

## Dispatch Architecture

Every public gate validates inputs with `GATE_VALIDATE`, then dispatches:

```c
switch (state->type) {
    case PURE:         gate_pure(state, ...);   break;
    case MIXED_PACKED: gate_packed(state, ...);  break;
    case MIXED_TILED:  gate_tiled(state, ...);   break;
    default: GATE_VALIDATE(0, "gate: unknown state type");
}
```

`GATE_VALIDATE(cond, msg)` prints to stderr and calls `exit(EXIT_FAILURE)` — intentional fail-fast for a library that never recovers from invalid input.

**Dispatch macros in `gate.c`:** `DISPATCH_1Q`, `DISPATCH_ROT`, `DISPATCH_2Q`, `DISPATCH_2Q_NO_TILED`. For CCX, use `DISPATCH_3Q` (to be added) or inline the switch.

---

## Algorithm Summaries

### Pure State (state_t type PURE)

**Single-qubit** — `TRAVERSE_PURE_1Q(state, target, PAIR_OP)`:
Processes pairs `(ψ[i], ψ[i + 2^target])` where the target bit differs. Nested loop: outer stride `2×2^target`, inner `2^target`. OpenMP: `if(dim >= 4096)`.

**Two-qubit** — flat loop with `insertBits2_0(k, lo, hi)`:
Enumerates `dim/4` base indices with both qubit bits = 0. No traversal macro; each gate inlines the loop. SWAP swaps the `(01, 10)` amplitude pair; CX/CY/CZ swap/phase the `(ctrl=1,tgt=0)` and `(ctrl=1,tgt=1)` amplitudes.

**Three-qubit** — flat loop with chained `insertBit0()`:
Enumerates `dim/8` base indices with all three qubit bits = 0. CCX swaps amplitudes at `(c1=1, c2=1, tgt=0)` and `(c1=1, c2=1, tgt=1)`.

### Mixed Packed (state_t type MIXED_PACKED)

Implements ρ → UρU† directly on the lower-triangular packed storage. No unpacking required.

**Single-qubit** — `TRAVERSE_PACKED_BLOCKS(state, target, BLOCK_OP)`:
Iterates lower-triangle pairs `(bc ≤ br)` of column/row base indices, each producing 4 packed indices `(idx00, idx01, idx10, idx11)`. `lower_01` and `diag_block` flags control conjugation and skipping. OpenMP `schedule(static)`, `if(dim >= 64)`.

**Two-qubit** — `insertBits2_0()` with 4×4 block pairs:
Two nested loops `(bc ≤ br)`, each producing 4 row and 4 column indices — a 4×4 sub-block of the density matrix. Fast path when `r00 > c11` (all in lower triangle); slow path for upper-triangle crossings. OpenMP `schedule(static, 1)`.

**Three-qubit** — same pattern with `insertBit0()` chained three times, producing 8×8 blocks. OpenMP `schedule(static, 1)`.

**Per-gate operation types in packed:**
- **Zero-flop** (X, CX, SWAP): only swap pairs — no arithmetic.
- **Phased permutation** (Y, CY): swaps with ±i phase factors.
- **Diagonal** (Z, S, Sdg, T, Tdg, P, Rz, CZ): only off-diagonal elements change; diagonal elements unchanged. Guard with `if (!diag_block)`.
- **Mixing** (H, Hy, Rx, Ry): all 4 elements of each 2×2 block participate.

### Mixed Tiled (state_t type MIXED_TILED)

TILE_DIM = 32, LOG_TILE_DIM = 5. Only lower-triangular tiles (tr ≥ tc) stored. Within diagonal tiles, only lr ≥ lc elements are populated.

**Single-qubit** — two traversal modes:

1. **Within-tile** (`target < LOG_TILE_DIM`): all 4 butterfly elements are in the same tile. Diagonal tiles require special handling — the (0,1) element at `(lr0, lc1)` may be in the upper triangle and must be read as `conj(elem[lc1, lr0])`.

2. **Cross-tile** (`target ≥ LOG_TILE_DIM`): butterfly spans 4 tiles T00/T01/T10/T11, grouped by tile-level bit `tile_bit = target - LOG_TILE_DIM`. T(tr0,tc1) may be an upper-triangular tile — access as conjugate of T(tc1,tr0).

**Shared macros in `gate_tiled.c`** (used by single-qubit non-rotation gates):
```c
TRAVERSE_TILED_WITHIN(state, target, OFFDIAG_OP, CONJ_OP)
  // OFFDIAG_OP(a00, a01, a10, a11)          — 4 args, off-diagonal blocks
  // CONJ_OP(a00, a01, a10, a11, lower_01, diag_block) — 6 args, diagonal blocks

TRAVERSE_TILED_CROSS(state, target, BLOCK_OP)
  // BLOCK_OP(a00, a01, a10, a11, lower_01, diag_block) — 6 args

DISPATCH_TILED_1Q(state, target, OFFDIAG_OP, CROSS_OP)
  // Dispatches within vs cross based on target < LOG_TILE_DIM
```

Rotation gates (Rx, Ry, Rz) cannot use these macros (need angle parameter in closure); they inline the traversal logic directly using `TRAVERSE_TILED_WITHIN_ROT` / `TRAVERSE_TILED_CROSS_ROT`.

**P/Rz delegation:** `p_packed()` and `p_tiled()` contain the actual implementation. `rz_packed()` and `rz_tiled()` are one-line wrappers — the global phase from Rz cancels in UρU†, making Rz and P identical for density matrices. `p_pure()` is distinct.

**Key pitfall — diag_block aliasing:** When `lr0 == lc0` (within-tile) or the two element pointers land on the same memory location (diagonal of diagonal tile), `a01` and `a10` alias the same location. Diagonal gates must guard with `if (!diag_block)`. Mixing gates must enforce Hermitian symmetry: after computing `a10`, set `a01 = conj(a10)`.

**Two-qubit tiled:** Same four-block decomposition as packed (Section 5.2 of GATES_MODULE.md). Triangle status determined once per tile group. See `gate_tiled_cx.c` as the reference implementation.

**OMP_THRESHOLD:** 512 (n < 9 qubits), defined in `gate_tiled.h`.

---

## Implementing `ccx_tiled` — Quick Reference

The `ccx_packed` implementation in `gate_packed_ccx.c` is the closest structural analog. The tiled version follows the same control-flow pattern as `cx_tiled` but with three-bit insertion.

**Algorithm sketch:**
1. Sort `(c1, c2, tgt)` into `lo < mid < hi`; compute `incr_c1`, `incr_c2`, `incr_tgt` as `1 << position`.
2. Outer loops over tile pairs `(btr, btc)` with `btr ≥ btc` (lower triangle), skipping tiles where neither `c1=c2=1` quadrant exists.
3. For each tile group, load 8 row indices and 8 column indices using chained `insertBit0()` at the tile level.
4. Inner loop processes the 8×8 block: when both controls are 1, swap the target-bit-0 and target-bit-1 element pair. Apply standard lower-triangle / conjugation logic for off-diagonal and diagonal tiles.

See the "Controlled Gates on Tiled Storage" section above and `gate_tiled_cx.c` for the concrete cross-tile traversal pattern.

---

## Test Infrastructure

### Running Tests

```bash
cd cmake-build-debug
cmake --build . --target test_gate && ./test/test_gate
```

**Note:** Do not use `cmake --build .` (full build) — `benchmark/src/bench_mixed.c` includes stale headers that no longer exist.

### Test Harnesses

| Harness | File | State Type |
|---------|------|------------|
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
| `testThreeQubitGateTiled(fn, mat)` | `test_gate_tiled_2q.c` | **DECLARED, NOT IMPLEMENTED** |

### Adding a New Test (e.g., `ccx_tiled`)

In `test/gate/test_gate.c`:
```c
void test_CCX_tiled(void) { testThreeQubitGateTiled(ccx, CCXMAT); }

// In main():
RUN_TEST(test_CCX_tiled);
```

Gate matrix constants (`CCXMAT`, `CXMAT`, `CYMAT`, `CZMAT`, `SWAPMAT`, etc.) are in `test/include/test_gate.h`. Reference matrix builders (`mat_three_qubit_gate`, `mat_two_qubit_gate`, `mat_single_qubit_gate`) are in `test/utility/gatemat.c`.

### Test Coverage Design

- **System sizes:** 1Q gates: n=1,2,3,4; 2Q gates: n=2,3,4; 3Q gates: n=3,4
- **All qubit positions** tested (all targets, all ctrl/tgt pairs, all triples)
- **Test states:** |0⟩, |1⟩, |+⟩, |−⟩, |+i⟩, |−i⟩ plus Bell, GHZ, W states for pure; same as density matrices plus maximally mixed and 10 random mixtures (fixed seed) for mixed
- **Rotation angles:** θ ∈ {0, π/4, π/2, π, 3π/2, −π/3}
- **Tiled edge case:** `MAXQUBITS_TILED=6` to exercise cross-tile paths (n=5 is single-tile boundary, n=6 requires multiple tiles)
- **Tolerance:** ε = 1e-12
- **Reference:** pure states via `zgemv`; mixed states via `zgemm` (UρU†)

---

## Known Issues / Deferred

- **Benchmark build broken:** `benchmark/src/bench_mixed.c` imports `state_blocked.h` and `mhipster_block.h` which no longer exist. Fix by updating the benchmark source.
- **RZZ gate:** Diagonal two-qubit gate useful for QAOA. Deferred to a future revision.
- **Higher-order controlled gates** (CS, CT, CH, etc.): Out of scope for this module.
