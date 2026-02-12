# Gates Module Implementation Status

**Last updated:** 2026-02-12
**Last verified build:** 36 gate tests, 0 failures (all at 1e-12 tolerance)

## How to Build and Test

```bash
cd /Users/timo/Code/qlib/cmake-build-debug

# Build only the test targets (benchmark target is broken, see Known Issues)
cmake --build . --target test_state --target test_gate

# Run tests
./test/test_gate
```

**Known issue:** The full `cmake --build .` fails because `benchmark/src/bench_mixed.c`
includes stale headers (`state_blocked.h`, `mhipster_block.h`) that no longer exist. Use
`--target test_gate` to build only the test targets.

---

## Completed Work

### Phase 0: Infrastructure (DONE)

All infrastructure refactoring is complete:

1. **File renames** (`git mv`): `qhipster.c` -> `gate_pure.c`, `mhipster.c` -> `gate_packed.c`,
   test files similarly renamed.
2. **Suffix rename**: `_mixed` -> `_packed` in all internal function names.
3. **GATE_VALIDATE macro** added to `include/gate.h` with `exit(EXIT_FAILURE)` semantics.
   All gate signatures changed from `qs_error_t` to `void`.
4. **Three-way dispatch** in `src/gate/gate.c`: `switch(state->type)` dispatches to `*_pure()`,
   `*_packed()`, `*_tiled()` for all 18 public gates. Macros: `DISPATCH_1Q`, `DISPATCH_ROT`,
   `DISPATCH_2Q`, `DISPATCH_2Q_NO_TILED`.
5. **Tolerance tightened** to `1e-12` in `test/include/test.h`.
6. **Tiled test harnesses** created in `test/gate/test_gate_tiled_1q.c` and `test/gate/test_gate_tiled_2q.c`.
7. **Tiled gate files created**: `src/gate/tiled/gate_tiled.c` (1Q + rotation), `src/gate/tiled/gate_tiled_cx.c`, `src/gate/tiled/gate_tiled_swap.c`. Shared inline helpers in `src/gate/tiled/gate_tiled.h`.
8. **CMakeLists.txt updated** for all renames and new files.

### Phase 1: Single-Qubit Gates (DONE)

All 13 single-qubit gates implemented across all 3 memory layouts (pure, packed, tiled):

| Gate | Pure | Packed | Tiled | Status |
|------|------|--------|-------|--------|
| X | pre-existing | pre-existing | new | TESTED |
| Y | pre-existing | pre-existing | new | TESTED |
| Z | pre-existing | pre-existing | new | TESTED |
| H | pre-existing | pre-existing | new | TESTED |
| S | pre-existing | pre-existing | new | TESTED |
| Sdg | pre-existing | pre-existing | new | TESTED |
| T | pre-existing | pre-existing | new | TESTED |
| Tdg | pre-existing | pre-existing | new | TESTED |
| Hy | **new** | **new** | **new** | TESTED |
| Rx | pre-existing | pre-existing | new | TESTED |
| Ry | pre-existing | pre-existing | new | TESTED |
| Rz | pre-existing | pre-existing | new | TESTED |
| P | **new** | **new** | **new** | TESTED |

**P/Rz delegation note:** `rz_packed()` and `rz_tiled()` delegate to `p_packed()` and
`p_tiled()` respectively, since `P(θ)=diag(1, e^(iθ))` and `Rz(θ)=diag(e^(-iθ/2), e^(iθ/2))`
produce identical density matrix conjugations (the global phase cancels in ρ'=UρU†).
`p_packed()` contains the actual implementation; `rz_packed()` is a one-line wrapper.
`p_pure()` is distinct — it only phases the |1⟩ amplitude.

---

## Phase 2: Two-Qubit Gates (DONE for pure + packed)

4 two-qubit gates implemented across pure and packed layouts. Tiled dispatch temporarily
disabled (under construction).

| Gate | Pure | Packed | Tiled | Test matrix constant |
|------|------|--------|-------|---------------------|
| CX | done | done | done* | `CXMAT` |
| CY | done | done | — | `CYMAT` |
| CZ | done | done | — | `CZMAT` |
| SWAP | done | done | done* | `SWAPMAT` |

*Tiled implementations exist but dispatch is temporarily disabled.

**Implementation locations:**
- CX: `cx_pure` in `gate_pure.c`, `cx_packed` in `packed/gate_packed_cx.c`, `cx_tiled` in `tiled/gate_tiled_cx.c`
- CY: `cy_pure` in `gate_pure.c`, `cy_packed` in `packed/gate_packed_cy.c`
- CZ: `cz_pure` in `gate_pure.c`, `cz_packed` in `packed/gate_packed_cz.c`
- SWAP: `swap_pure` in `gate_pure.c`, `swap_packed` in `packed/gate_packed_swap.c`, `swap_tiled` in `tiled/gate_tiled_swap.c`

### Test Harness Status

Test harnesses that **exist and work:**
- `testTwoQubitGate()` — pure
- `testTwoQubitGateMixed()` — packed
- `testTwoQubitGateTiled()` — tiled
- `testThreeQubitGate()` — pure
- `testThreeQubitGateMixed()` — packed
- `mat_three_qubit_gate()` — reference matrix builder for 3Q gates

Test harnesses that are **declared but NOT implemented:**
- `testThreeQubitGateTiled()` — tiled (needed when tiled CCX is implemented)

### Implementation Pattern

All two-qubit gates use flat loops with `insertBits2_0()` for index enumeration (no traversal macro):
- **Pure:** iterate over base indices; when control=1, apply single-qubit operation
- **Packed:** two-bit insertion with `insertBits2_0()`, processing swap/transform pairs per 4×4 block
- **Tiled:** same two-bit insertion logic but with tile coordinate arithmetic

**SWAP is different** — it's a permutation gate, not a controlled gate. It swaps basis states
|...01...> <-> |...10...> where the two qubit positions differ.

### Adding a New Test

For each new gate, add to `test/gate/test_gate.c`:

```c
// Function definition
void test_CY_pure(void) {testTwoQubitGate(cy, CYMAT);}
void test_CY_packed(void) {testTwoQubitGateMixed(cy, CYMAT);}
void test_CY_tiled(void) {testTwoQubitGateTiled(cy, CYMAT);}

// In main(), add:
RUN_TEST(test_CY_pure);
RUN_TEST(test_CY_packed);
RUN_TEST(test_CY_tiled);
```

---

## Phase 3: Three-Qubit Gate — Toffoli (DONE for pure + packed)

| Component | Status |
|-----------|--------|
| `mat_three_qubit_gate()` reference builder | done |
| `CCXMAT` 8x8 test constant | done |
| `testThreeQubitGate()` (pure) | done |
| `testThreeQubitGateMixed()` (packed) | done |
| `testThreeQubitGateTiled()` (tiled) | declared, not implemented |
| `ccx_pure` | done |
| `ccx_packed` | done |
| `ccx_tiled` | not implemented |

**Implementation locations:**
- CCX: `ccx_pure` in `gate_pure.c`, `ccx_packed` in `packed/gate_packed_ccx.c`

---

## Architecture Notes for Implementers

### Tiled Storage Model — Critical Details

The tiled storage is the most complex layout. Key facts discovered during Phase 1:

1. **Tile coordinates:** A density matrix of dimension `dim = 2^n` is partitioned into
   tiles of size `TILE_DIM x TILE_DIM` (TILE_DIM=32, LOG_TILE_DIM=5). Only lower-triangular
   tiles (tr >= tc) are stored.

2. **Within tiles:** Elements are stored **row-major**. For off-diagonal tiles (tr > tc),
   the full tile is stored. For diagonal tiles (tr == tc), **only the lower triangle
   (lr >= lc) is populated** — upper-triangle positions contain zero/garbage.

3. **Hermitian symmetry:** `state_get(row, col)` automatically returns `conj(state_get(col, row))`
   for upper-triangle accesses. `state_set(row, col, val)` stores only if row >= col.

4. **Two traversal modes for 1Q gates:**
   - **Within-tile** (target < LOG_TILE_DIM): butterfly elements are in the same tile.
     Diagonal tiles need special handling — the (0,1) element at `(lr0, lc1)` may be
     in the upper triangle of the tile, requiring conjugate read from `(lc1, lr0)`.
   - **Cross-tile** (target >= LOG_TILE_DIM): butterfly spans 4 tiles.
     T(tr0,tc1) may be an upper-triangular tile, stored as conj of T(tc1,tr0).

5. **diag_block aliasing:** When `lr0==lc0` (within-tile) or `lr==lc` (cross-tile),
   the (0,1) and (1,0) pointers can alias the same memory location (since row0==col0
   means the conjugate pair shares storage). For diagonal gates (Z, S, T, etc.),
   this means you must **not** operate on both — use `if (!diag_block)` guards.
   For mixing gates (H, Hy, Rx, Ry), the diag_block case must enforce Hermitian
   symmetry: set a01 = conj(a10) after computing a10.

### Key Macros in `src/gate/tiled/gate_tiled.c`

```
TRAVERSE_TILED_WITHIN(state, target, OFFDIAG_OP, DIAG_OP)
  - OFFDIAG_OP(a00, a01, a10, a11) — 4 args, for off-diagonal tiles
  - DIAG_OP(a00, a01, a10, a11, lower_01, diag_block) — 6 args, for diagonal tiles

TRAVERSE_TILED_CROSS(state, target, BLOCK_OP)
  - BLOCK_OP(a00, a01, a10, a11, lower_01, diag_block) — 6 args

DISPATCH_TILED_1Q(state, target, OFFDIAG_OP, CROSS_OP)
  - Dispatches within vs cross based on target < LOG_TILE_DIM
```

For rotation gates (Rx, Ry, Rz), the angle parameter prevents using macros; the traversal
logic is inlined directly in the function body.

### Key Macros in `src/gate/packed/gate_packed_1q.c`

```
TRAVERSE_PACKED_BLOCKS(state, target, BLOCK_OP)
  - BLOCK_OP(data, idx00, idx01, idx10, idx11, lower_01, diag_block)
  - data is the packed array, idx* are packed indices
```

Multi-qubit packed gates (`src/gate/packed/gate_packed_<name>.c`, one file per gate) do not use a
traversal macro; each gate inlines its own loop with gate-specific pair handling.

### Key Macros in `src/gate/gate_pure.c`

```
TRAVERSE_PURE_1Q(state, target, PAIR_OP)
  - PAIR_OP(a, b) — pointers to amplitudes with target bit 0 and 1
```

No multi-qubit traversal macro exists. All two-qubit and three-qubit pure gates inline their flat loops.

### Public API Naming

- All gates: `void name(state_t*, ...)` in `include/gate.h`
- Internal: `name_pure()`, `name_packed()`, `name_tiled()` — declared `extern` in `src/gate/gate.c`
- Exception: SWAP is `swap_gate()` in the public API (avoids C keyword clash), but
  internal functions are `swap_pure()`, `swap_packed()`, `swap_tiled()`

### CX Implementation (reference for 2Q pattern)

**`cx_pure` (gate_pure.c:282):** Uses `insertBits2_0()` to enumerate base indices with
both control and target bits = 0. Then for each base, swaps amplitudes where control=1.
Simple loop: `n_base = dim/4` iterations.

**`cx_packed` (gate_packed_cx.c):** In-place permutation using `insertBits2_0()` to
enumerate quarter-dim blocks. Each block produces 4 row × 4 column indices; 6 swap
pairs are processed per block. Fast path (r00 > c11) handles all pairs in lower triangle;
slow path handles upper-triangle crossings with conjugation. CY, CZ, and SWAP follow the
same structural pattern.

---

## File Reference

| File | Lines | Purpose |
|------|-------|---------|
| `include/gate.h` | 170 | GATE_VALIDATE, 18 gate declarations |
| `src/gate/gate.c` | ~210 | Dispatchers with validation, extern declarations |
| `src/gate/gate_pure.c` | ~435 | 18 pure implementations (13 1Q + 4 2Q + 1 3Q) |
| `src/gate/packed/gate_packed_1q.c` | ~700 | 13 single-qubit packed implementations |
| `src/gate/packed/gate_packed_<name>.c` | varies | One file per multi-qubit packed gate (CX, CY, CZ, SWAP, CCX) |
| `src/gate/tiled/gate_tiled.c` | ~845 | 13 tiled 1Q + rotation implementations |
| `src/gate/tiled/gate_tiled_cx.c` | 721 | cx_tiled implementation |
| `src/gate/tiled/gate_tiled_swap.c` | 503 | swap_tiled implementation |
| `src/gate/tiled/gate_tiled.h` | 70 | Shared inline helpers (tile_off, elem_off, dtile_read/write) |
| `test/gate/test_gate.c` | ~165 | Test function defs + main() with 36 RUN_TESTs |
| `test/gate/test_gate_pure_1q.c` | ~270 | Pure 1Q harnesses + rotation matrix builders |
| `test/gate/test_gate_pure_2q.c` | ~120 | Pure 2Q harness |
| `test/gate/test_gate_pure_3q.c` | ~125 | Pure 3Q harness |
| `test/gate/test_gate_packed_1q.c` | ~275 | Packed 1Q harnesses |
| `test/gate/test_gate_packed_2q.c` | ~150 | Packed 2Q harness |
| `test/gate/test_gate_packed_3q.c` | ~155 | Packed 3Q harness |
| `test/gate/test_gate_tiled_1q.c` | ~225 | Tiled 1Q harnesses |
| `test/gate/test_gate_tiled_2q.c` | ~130 | Tiled 2Q harness |
| `test/include/test_gate.h` | ~220 | Typedefs, matrices, harness declarations |
| `test/utility/gatemat.c` | ~195 | Reference matrix builders (1Q, 2Q, 3Q) |
| `test/CMakeLists.txt` | ~105 | Build config for test_state and test_gate |

---

## Full Plan

The original multi-phase plan is at `/Users/timo/.claude/plans/silly-squishing-harbor.md`.
It contains detailed derivations and implementation notes for each gate. Refer to it for
the conjugation formulas, especially for the controlled-gate four-block decomposition pattern.

The authoritative gate spec is at `docs/GATES_MODULE.md`.
