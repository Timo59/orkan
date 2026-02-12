# Gates Module Implementation Status

**Last updated:** 2026-02-06
**Last verified build:** 41 tests, 0 failures (all at 1e-12 tolerance)

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
4. **Three-way dispatch** in `src/gate.c`: `switch(state->type)` dispatches to `*_pure()`,
   `*_packed()`, `*_tiled()` for all 26 gates. Macros: `DISPATCH_1Q`, `DISPATCH_ROT`,
   `DISPATCH_2Q`, `DISPATCH_2Q_ROT`, `DISPATCH_2Q_NO_TILED`, `DISPATCH_2Q_ROT_NO_TILED`.
5. **Tolerance tightened** to `1e-12` in `test/include/test.h`.
6. **Tiled test harnesses** created in `test/src/test_gate_tiled.c` with all needed helpers.
7. **Tiled gate files created**: `src/gate_tiled.c` (1Q + rotation), `src/gate_tiled_cx.c`, `src/gate_tiled_swap.c`. Shared inline helpers in `src/gate_tiled.h`.
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

## Current Work: Phase 2 — Two-Qubit Gates (NOT STARTED)

### What Needs to Be Done

12 two-qubit gates need implementation across 3 layouts = 36 functions.

**Implemented two-qubit gates:**
- CX (CNOT): `cx_pure` in `src/gate_pure.c`, `cx_packed` in `src/gate_packed_cx.c`, `cx_tiled` in `src/gate_tiled_cx.c`
- SWAP: `swap_packed` in `src/gate_packed_swap.c`, `swap_tiled` in `src/gate_tiled_swap.c`

**Remaining stubs** (pure/packed call `GATE_VALIDATE(0, ...)`; tiled dispatchers reject with "not implemented for tiled"):

| Gate | Pure | Packed | Tiled | Test matrix constant |
|------|------|--------|-------|---------------------|
| CX | done | done | done | `CXMAT` (exists) |
| CY | STUB | STUB | — | `CYMAT` (exists) |
| CZ | STUB | STUB | — | `CZMAT` (exists) |
| CS | STUB | STUB | — | `CSMAT` (exists) |
| CSdg | STUB | STUB | — | `CSDGMAT` (exists) |
| CH | STUB | STUB | — | declared `mat_ch_build()` **NOT IMPLEMENTED** |
| CHy | STUB | STUB | — | declared `mat_chy_build()` **NOT IMPLEMENTED** |
| CT | STUB | STUB | — | `CTMAT` (exists) |
| CTdg | STUB | STUB | — | `CTDGMAT` (exists) |
| CP(θ) | STUB | STUB | — | declared `mat_cp()` **NOT IMPLEMENTED** |
| CPdg(θ) | STUB | STUB | — | declared `mat_cpdg()` **NOT IMPLEMENTED** |
| SWAP | STUB | STUB | done | `SWAPMAT` (exists) |

### Test Harness Status

Test harnesses that **already exist and work:**
- `testTwoQubitGate()` — pure (in `test/src/test_gate_pure.c:186`)
- `testTwoQubitGateMixed()` — packed (in `test/src/test_gate_packed.c:105`)
- `testTwoQubitGateTiled()` — tiled (in `test/src/test_gate_tiled.c:310`)
- `testTwoQubitRotationGateTiled()` — tiled parameterized 2Q (in `test/src/test_gate_tiled.c:403`)

Test harnesses that are **declared but NOT implemented:**
- `testTwoQubitRotationGate()` — pure parameterized 2Q (needed for CP, CPdg)
- `testTwoQubitRotationGateMixed()` — packed parameterized 2Q (needed for CP, CPdg)
- `testThreeQubitGate()` — pure (needed for Phase 3)
- `testThreeQubitGateMixed()` — packed (needed for Phase 3)
- `testThreeQubitGateTiled()` — tiled (needed for Phase 3)

Matrix builder functions that are **declared but NOT implemented:**
- `mat_ch_build()` — builds 4x4 CH matrix (contains irrational √2 entries)
- `mat_chy_build()` — builds 4x4 CHy matrix
- `mat_cp(double theta, cplx_t *mat)` — builds 4x4 CP(θ) matrix
- `mat_cpdg(double theta, cplx_t *mat)` — builds 4x4 CPdg(θ) matrix

### Recommended Implementation Order

The original plan (in `/Users/timo/.claude/plans/silly-squishing-harbor.md`) recommended:

1. **CX tiled** first (establishes tiled 2Q traversal pattern)
2. **TRAVERSE_PURE_2Q macro** (optional: refactor cx_pure to use it, then all 2Q pure gates reuse)
3. Then gate-by-gate: CY, CZ, CS, CSdg, CH, CHy, CT, CTdg, CP, CPdg, SWAP

For each gate (except SWAP): the pattern is a **controlled-U** gate where:
- Pure: iterate over (control,target) pairs; when control=1, apply single-qubit U_OP
- Packed: two-bit insertion with `insertBits2_0()`, processing swap/transform pairs per 4×4 block
- Tiled: same two-bit insertion logic but with tile coordinate arithmetic

**SWAP is different** — it's a permutation gate, not a controlled gate. It swaps basis states
|...01...> <-> |...10...> where the two qubit positions differ.

### Adding a New Test

For each new gate, add to `test/src/test_gate.c`:

```c
// Function definition (near line 80)
void test_CY_pure(void) {testTwoQubitGate(cy, CYMAT);}
void test_CY_packed(void) {testTwoQubitGateMixed(cy, CYMAT);}
void test_CY_tiled(void) {testTwoQubitGateTiled(cy, CYMAT);}

// In main(), add:
RUN_TEST(test_CY_pure);
RUN_TEST(test_CY_packed);
RUN_TEST(test_CY_tiled);
```

For parameterized gates (CP, CPdg), use `testTwoQubitRotationGate*` harnesses instead.

---

## Phase 3: Three-Qubit Gate — Toffoli (NOT STARTED)

Single gate (CCX) across 3 layouts. Requires:
1. `mat_three_qubit_gate()` in `test/src/gatemat.c` (builds 8x8 Kronecker product)
2. `CCXMAT` 8x8 constant in test header
3. `testThreeQubitGate()`, `testThreeQubitGateMixed()`, `testThreeQubitGateTiled()` harness implementations
4. CCX implementations in all 3 layout files

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

### Key Macros in `src/gate_tiled.c`

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

### Key Macros in `src/gate_packed_1q.c`

```
TRAVERSE_PACKED_BLOCKS(state, target, BLOCK_OP)
  - BLOCK_OP(data, idx00, idx01, idx10, idx11, lower_01, diag_block)
  - data is the packed array, idx* are packed indices
```

Two-qubit packed gates (`src/gate_packed_<name>.c`, one file per gate) do not use a
traversal macro; each gate inlines the `insertBits2_0()` loop with gate-specific pair handling.

### Key Macros in `src/gate_pure.c`

```
TRAVERSE_PURE_1Q(state, target, PAIR_OP)
  - PAIR_OP(a, b) — pointers to amplitudes with target bit 0 and 1
```

No 2Q traversal macro exists yet. `cx_pure` uses `insertBits2_0()` manually.

### Public API Naming

- All gates: `void name(state_t*, ...)` in `include/gate.h`
- Internal: `name_pure()`, `name_packed()`, `name_tiled()` — declared `extern` in `src/gate.c`
- Exception: SWAP is `swap_gate()` in the public API (avoids C keyword clash), but
  internal functions are `swap_pure()`, `swap_packed()`, `swap_tiled()`

### Existing CX Implementations (reference for 2Q pattern)

**`cx_pure` (gate_pure.c:189):** Uses `insertBits2_0()` to enumerate base indices with
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
| `include/gate.h` | 119 | GATE_VALIDATE, all 26 gate declarations |
| `src/gate.c` | 268 | Dispatchers with validation, extern declarations |
| `src/gate_pure.c` | 395 | 14 pure implementations + 12 stubs |
| `src/gate_packed_1q.c` | ~710 | 13 single-qubit packed implementations |
| `src/gate_packed_<name>.c` | varies | One file per two-qubit packed gate (CX, CY, CZ, SWAP + stubs) |
| `src/gate_packed_3q.c` | ~20 | Three-qubit packed gate stubs (CCX) |
| `src/gate_tiled.c` | 791 | 13 tiled 1Q + rotation implementations |
| `src/gate_tiled_cx.c` | 721 | cx_tiled implementation |
| `src/gate_tiled_swap.c` | 503 | swap_tiled implementation |
| `src/gate_tiled.h` | 70 | Shared inline helpers (tile_off, elem_off, dtile_read/write) |
| `test/src/test_gate.c` | ~150 | Test function defs + main() with 41 RUN_TESTs |
| `test/src/test_gate_pure.c` | ~370 | Pure harnesses + rotation matrix builders |
| `test/src/test_gate_packed.c` | ~180 | Packed harnesses |
| `test/src/test_gate_tiled.c` | ~485 | Tiled harnesses (1Q, 2Q, 2Q-rotation all done) |
| `test/include/test_gate.h` | 241 | Typedefs, matrices, harness declarations |
| `test/CMakeLists.txt` | 117 | Build config for test_state and test_gate |

---

## Full Plan

The original multi-phase plan is at `/Users/timo/.claude/plans/silly-squishing-harbor.md`.
It contains detailed derivations and implementation notes for each gate. Refer to it for
the conjugation formulas, especially for the controlled-gate four-block decomposition pattern.

The authoritative gate spec is at `docs/GATES_MODULE.md`.
