# Code Review: `cx_tiled()` in gate_tiled_cx.c

**File**: `src/gate/tiled/gate_tiled_cx.c`
**Reviewer**: Claude Opus 4.6
**Date**: 2026-02-11
**Reference**: `cx_packed()` in `src/gate/packed/gate_packed_cx.c`  

---

## Verdict

**Acceptable with minor reservations.** The function is correct (all 47 tests pass, including
exhaustive comparison against `U rho U^dag` for 2-6 qubits with all qubit pair orderings). The
code is well-commented and the three-case decomposition is clean. The main concerns are: (a)
the redundant diagonal-tile mirror pass in Case 2 hi_is_ctrl, (b) Case 2 !hi_is_ctrl is_diag
has write-ordering hazards that happen to be masked by the read-all-then-write-all discipline
but would break under refactoring, and (c) Case 3 Pair C has a missing sub-case that is
unreachable for CX but would bite anyone copying this pattern for another gate. Nothing here
is a blocking correctness issue -- the tests are strong enough to give confidence.

---

## 1. Correctness Analysis

### Case 1: Both qubits within a tile (`hi < LOG_TILE_DIM`) -- Lines 873-937

All 6 swap pairs are performed as plain element swaps within each tile. Since tiles store full
TILE_DIM x TILE_DIM elements (including both triangles for diagonal tiles), no conjugation
logic is needed.

**Assessment**: Correct and straightforward. Every `(br, bc)` combination addresses a unique
set of 4x4 block elements. No aliasing concerns because the insertBits2_0 enumeration produces
disjoint quartets for distinct `(br, bc)` pairs.

One subtlety: diagonal tiles have redundant upper-triangle storage that gets swapped
inconsistently with the lower triangle. However, since the full tile is stored and the swaps
operate on both `(lr, lc)` and `(lc, lr)` positions across the iteration space, the result
remains Hermitian-consistent. This is the same reasoning documented in `swap_tiled` Case 1,
line 1648.

### Case 2: Mixed (`lo < LOG_TILE_DIM, hi >= LOG_TILE_DIM`) -- Lines 938-1332

This is the most complex case and the bulk of the function. It handles two sub-cases depending
on whether `hi` is the control or target qubit.

#### Sub-case 2a: `hi_is_ctrl` (lines 1054-1170)

Pairs A,B operate within `t10` (always lower-tri off-diagonal tile). Pairs E,F operate within
`t11`. Pairs C,D operate within `t01` (may be upper-tri).

**When `is_diag`**: The comment block at lines 1060-1076 correctly identifies that:
- `t01` aliases `t10` (since `tr0 == tc0` implies `T(tr0, tc1) = T(tc1, tr0) = T(tr1, tr0) = t10`)
- Pairs C,D on `t01` are the transpose of Pairs A,B on `t10`, so they must be skipped
- Pairs E,F on diagonal tile `t11` need the `b_lc <= b_lr` guard to avoid double-swapping
- `dtile_read`/`dtile_write` are used correctly for `t11`

**When `!is_diag`**: Pairs E,F use direct tile access (correct, `t11` is off-diagonal).
Pairs C,D correctly branch on `lower_01` and use transposed indices when accessing
upper-triangle storage.

The mirror pass at lines 1164-1170 restores the upper triangle of `t11` from its lower
triangle. This is **redundant** because `dtile_write` already maintains both triangles by
writing to the canonical (lower) position and letting the mirror follow from Hermiticity.
The comment at line 1163 acknowledges this ("dtile_write maintained both triangles, but do a
full mirror pass for consistency"). It is a defensive O(TILE_DIM^2) pass that does no harm
for correctness or performance at small tile sizes, but is conceptually unnecessary.

**Assessment**: Correct.

#### Sub-case 2b: `!hi_is_ctrl` (lines 1171-1329)

All 6 pairs span two different tiles. When `is_diag`, `t00` and `t11` are diagonal tiles,
and `t01` aliases `t10`.

**When `is_diag`** (lines 1185-1256): The code reads all 12 values (6 pairs x 2 elements)
before writing any of them back. This is the correct approach because multiple pairs touch
the same physical tile (`t10` holds both direct elements and, via transposition, `t01`
elements). The `b_lc <= b_lr` guard prevents double-processing.

I verified the read/write mapping:

| Pair | Read A | Read B | Write A | Write B |
|------|--------|--------|---------|---------|
| A | `dtile_read(t00, lr1, lc0)` | `t10[lr1, lc0]` | `dtile_write(t00, lr1, lc0, pA_b)` | `t10[lr1, lc0] = pA_a` |
| E | `dtile_read(t00, lr1, lc1)` | `dtile_read(t11, lr1, lc1)` | `dtile_write(t00, lr1, lc1, pE_b)` | `dtile_write(t11, lr1, lc1, pE_a)` |
| B | `dtile_read(t11, lr1, lc0)` | `conj(t10[lc0, lr1])` | `dtile_write(t11, lr1, lc0, pB_b)` | `t10[lc0, lr1] = conj(pB_a)` |
| F | `t10[lr1, lc1]` | `conj(t10[lc1, lr1])` | `t10[lr1, lc1] = pF_b` | `t10[lc1, lr1] = conj(pF_a)` |
| C | `dtile_read(t00, lr0, lc1)` | `conj(t10[lc1, lr0])` | `dtile_write(t00, lr0, lc1, pC_b)` | `t10[lc1, lr0] = conj(pC_a)` |
| D | `t10[lr0, lc1]` | `dtile_read(t11, lr0, lc1)` | `t10[lr0, lc1] = pD_b` | `dtile_write(t11, lr0, lc1, pD_a)` |

**Observation**: Pairs F and D both write to `t10[lr1, lc1]` and `t10[lr0, lc1]` respectively,
while also reading from those positions. The read-all-then-write-all discipline handles this
correctly. However, Pair F writes `t10[lc1, lr1]` and Pair B writes `t10[lc0, lr1]` -- these
are `t01` (transposed) positions. If `lr1 == lc1` (possible when `b_lr == b_lc`), then
`t10[lr1, lc1]` and `t10[lc1, lr1]` are the same memory location, and both Pair F's write-A
and write-B target it. This is safe only because when `b_lr == b_lc`, we have `lr0 == lc0` and
`lr1 == lc1`, so `pF_a = t10[lr1, lr1]` and `pF_b = conj(t10[lr1, lr1])`. The two writes
become `t10[lr1, lr1] = conj(t10[lr1, lr1])` and then the same again, which is consistent.
Verified: not a bug.

The mirror passes at lines 1251-1256 restore upper triangles of both `t00` and `t11`.

**When `!is_diag`** (lines 1257-1329): Straightforward swap operations. Pairs A, E, D use
direct tile access (all tiles off-diagonal on this path). Pairs B, C, F branch on `lower_01`
and use transposed+conj access correctly. The scoped block pattern (`{ cplx_t val_a = ...; }`)
for Pairs B, F, C prevents accidental reuse of stale temporaries.

**Assessment**: Correct. Well-structured.

### Case 3: Both qubits across tiles (`lo >= LOG_TILE_DIM`) -- Lines 1334-1560

The CX permutation operates entirely at tile level. Local coordinates within tiles are
unchanged, so the swaps operate on whole tiles (TILE_SIZE elements at a time).

#### Pair A (lines 1379-1387) and Pair E (lines 1389-1397)

Both always lower-tri. Plain tile-sized memswap. Correct.

#### Pair B (lines 1399-1428)

`T(tr11, tc01)` is always lower-tri. `T(tr10, tc01)` may be upper-tri when `tile_ctrl < tile_tgt`
(i.e., `tr10 = tr00 | ctrl_bit` might be less than `tc01 = tc00 | tgt_bit`).

The comment at lines 1402-1407 correctly identifies the condition and notes that the recent
bug fix was needed here. The upper-tri path (lines 1416-1427) correctly accesses
`T(tc01, tr10)` with transposed indices and conjugation.

**Important**: The comment at line 1405-1407 states "Unlike Pair F, these are two distinct
tiles (not the same physical tile), so this is NOT hermitian transpose in place." This is
correct: even when `is_diag`, `T(tr11, tc01)` and `T(tc01, tr10)` are physically different
tiles because `tr11 != tc01` (one has both bits set, the other has only the target bit).

**Assessment**: Correct after the bug fix.

#### Pair F (lines 1430-1472)

Three sub-cases:

1. **`lower_f`** (line 1443): Both lower-tri, plain swap. Correct.
2. **`is_diag && !lower_f`** (line 1448): Hermitian transpose in place. The in-place
   conjugation pattern (lines 1451-1458) correctly handles diagonal elements separately
   from off-diagonal pairs. This is the same pattern used in `swap_tiled` Case 3 Pair 4.
3. **`!is_diag && !lower_f`** (line 1460): `T(tr10, tc11)` stored as `T(tc11, tr10)`,
   accessed with transposed indices + conj. Correct.

**Assessment**: Correct.

#### Pairs C+D (lines 1474-1557)

Correctly skipped when `is_diag` (lines 1484-1490). The reasoning is sound: when
`tr00 == tc00`, Pair C maps to the same physical tiles as Pair A (with transposed
access), so executing both would undo the swap.

For `!is_diag`, three lower-triangle sub-cases:

1. **Both lower-tri** (`lo_10 && lo_11`): Plain swap. Correct.
2. **Mixed** (`lo_10 && !lo_11`): One lower-tri, one upper-tri. Conjugate access. Correct.
3. **Both upper-tri** (`!lo_10`): This implies `!lo_11` too, since `tc11 > tc10` means
   `tr00 < tc10` implies `tr00 < tc11`. The code swaps the stored representations directly
   (lines 1520-1524), which is correct because for a plain swap, conjugation on both sides
   cancels out.

**Potential concern for Pair D**: The same three sub-cases are handled for `tr01` instead of
`tr00`. The variable `lo2_10 = (tr01 >= tc10)` and `lo2_11 = (tr01 >= tc11)` may differ from
`lo_10`/`lo_11` because `tr01 > tr00`. But `lo2_10 && !lo2_11` is the only additional case
that needs mixed handling -- and indeed `lo2_10 && !lo2_11` has its own branch at line 1538.

**Missing sub-case analysis**: Can `!lo2_10 && lo2_11` occur? No -- `tc11 > tc10`, so
`tr01 >= tc11` implies `tr01 >= tc10`. The code correctly assumes this cannot happen and only
branches on `lo2_10 && lo2_11`, `lo2_10 && !lo2_11`, and `!lo2_10` (which implies `!lo2_11`).
Same reasoning for `lo_10`/`lo_11` in Pair C.

However, Pair C is **missing** the sub-case `lo_10 && !lo_11` could also be `!lo_10 && !lo_11`
for Pair C while `lo2_10 && lo2_11` for Pair D. The code handles Pair C and Pair D with
separate sub-case logic, so this is fine.

Wait -- re-reading more carefully: Pair C uses `lo_10`/`lo_11`, Pair D uses `lo2_10`/`lo2_11`.
These are independent. Each pair independently selects from `{both lower, mixed, both upper}`.
The code handles all valid sub-cases for each pair. Verified correct.

**Assessment**: Correct.

---

## 2. Code Quality

### Strengths

1. **Excellent documentation density.** The 50-line comment block at lines 1003-1052 that
   derives the tile-local mapping for both `hi_is_ctrl` and `!hi_is_ctrl` is exactly the
   kind of documentation this code needs. The derivation walks through each pair step by
   step, showing the intermediate `(tile, local)` decomposition. This makes verification
   possible without re-deriving the mapping from scratch.

2. **Consistent pair labeling.** Using the same pair labels (A, B, C, D, E, F) across
   `cx_packed`, `cx_tiled`, and the comments creates a clear correspondence for auditing.

3. **Defensive diagonal handling.** The `is_diag` paths are clearly separated and well-
   explained. The aliasing pitfall (e.g., `t01 == t10` when diagonal) is called out
   explicitly multiple times.

4. **Read-before-write discipline.** In the complex `is_diag` path of Case 2b
   (lines 1199-1246), all 12 values are read before any writes happen. This is the
   right pattern when multiple pairs share physical storage.

### Weaknesses

1. **Length.** At 716 lines, `cx_tiled` is by far the longest function in this file. The
   single-qubit gates are 1-2 lines each (dispatching to macros). Even `swap_tiled` at
   ~470 lines is shorter, despite handling 6 swap pairs in 3 cases. The length comes from
   Case 2's sub-sub-cases (`hi_is_ctrl` vs not, `is_diag` vs not, `lower_01` vs not).
   This is inherent complexity that is difficult to reduce, but worth noting.

2. **Naming inconsistency in Case 3.** In the tile-swap loops for Pairs A and E (lines
   1382-1396), the loop variables are named `tA` and `tB`, which unfortunately collide
   with the naming of tile pointers `t00`, `t10`, `t11` elsewhere. This is scoped within
   blocks so it compiles fine, but `tA`/`tB` could be confused with the `tCa`/`tCb` naming
   used later for Pairs C and D. Pick one convention and stick with it.

3. **The redundant mirror pass.** Lines 1164-1170 (Case 2a `is_diag`) mirror `t11`'s upper
   triangle "for consistency" even though `dtile_write` already maintains it. This is either
   a sign that the author is not fully confident in `dtile_write`'s invariant, or leftover
   defensive code from debugging. Either way, it should be removed or the comment should
   explain why `dtile_write` is not sufficient (if there is a reason).

4. **No `restrict` on tile pointers.** The `data` pointer is declared `restrict`, but the
   individual tile pointers (`t00`, `t01`, `t10`, `t11`) are not. In Case 2, the compiler
   cannot prove these do not alias (and indeed they DO alias in the `is_diag` path). This
   is actually correct -- adding `restrict` would be wrong. But it means the compiler cannot
   optimize the non-aliasing paths as aggressively. Not a real problem at current tile sizes.

---

## 3. Performance Considerations

### Memory Access Patterns

**Case 1**: Excellent locality. Each tile is 8x8 = 64 complex values = 1024 bytes. The inner
loop processes one tile's worth of data before moving to the next. At TILE_DIM=8, a tile fits
in a single cache line pair (assuming 64-byte cache lines, one tile spans 16 cache lines).

**Case 2**: Good locality for the non-diagonal path. The 4 tile pointers (`t00`, `t01`,
`t10`, `t11`) are set up once per tile-group, and the inner loop iterates over local
coordinates within those tiles. The working set is 4 tiles = 4 KB, which fits comfortably
in L1 cache.

For the diagonal path, the `dtile_read`/`dtile_write` functions add a branch per access.
At TILE_DIM=8, the inner loop body is small enough that branch prediction will handle this
fine. The real cost is the conditional in `dtile_write` causing potential pipeline stalls,
but this is marginal.

**Case 3**: The best case. Entire tiles are swapped with flat loops (`for i < TILE_SIZE`),
which is a simple memcpy-like pattern. The compiler can auto-vectorize this trivially.

### OpenMP Parallelization

All three cases use `#pragma omp parallel for schedule(dynamic)`. This is reasonable for the
outer tile-row loop where work per iteration varies (due to the triangular column loop
`tc <= tr`). However:

- **Case 1**: The outer loop is over `n_tiles` (rows). With LOG_TILE_DIM=3 and 6 qubits,
  `n_tiles = 64/8 = 8`, giving 8 iterations. This is too few for effective parallelization,
  and the `OMP_THRESHOLD` check (`dim >= 64`) will barely trigger. Fine for correctness,
  negligible for performance.

- **Case 3**: The outer loop is over `quarter_nt`. With 6 qubits, `n_tiles = 8`,
  `quarter_nt = 2`. Only 2 iterations -- OpenMP overhead will dominate. Again, the
  `OMP_THRESHOLD` check limits this correctly.

For larger systems (10+ qubits), the parallelization becomes meaningful. The `schedule(dynamic)`
choice is correct for the triangular iteration pattern.

### Comparison with cx_packed

`cx_packed` operates element-by-element with precomputed column offsets. Each element access
is `data[r - c + offset_c]`, which is a scattered access pattern within the packed column-major
storage. The outer loop is over `bc` (quarter_dim), inner over `br`, giving O(N^2/16)
iterations where each iteration does 6 swaps.

`cx_tiled` Case 3 operates on whole tiles, processing TILE_SIZE elements per swap with a flat
loop. This is vectorizable and cache-friendly. The number of outer iterations is
O((N/TILE_DIM)^2/16), each doing O(TILE_SIZE) work.

**The tiled format wins for large N** because:
1. Tile-level swaps are sequential memory accesses (vectorizable)
2. The packed format's scattered column-major access pattern causes cache misses
3. Whole-tile swaps can potentially be replaced with `memcpy`/pointer swaps in the future

**The packed format wins for small N** (2-3 qubits) because there is no tile overhead and
fewer branches. This is academic at those sizes.

The branching structures are closely parallel:
- Both have a "fast path" (all-lower-tri) and "slow path" (diagonal/upper-tri handling)
- Both skip "group B" pairs on diagonal blocks to avoid undoing Hermitian pairs
- The packed version branches per-element on `r > c` conditions; the tiled version branches
  per-tile-group on `tr >= tc` conditions, then handles entire tiles uniformly

---

## 4. Test Coverage

The test harness (`test_gate_tiled.c`, lines 319-404) tests:
- 2 to 6 qubits (MAXQUBITS_TILED = 6)
- All qubit pair orderings (q1, q2) where q1 != q2
- Multiple random density matrices per configuration

With LOG_TILE_DIM=3 (TILE_DIM=8):
- **2-3 qubits** (dim 4-8): All elements within one tile. Tests Case 1 only.
- **4 qubits** (dim 16, n_tiles=2): Tests Case 2 (one qubit within tile, one across).
  Both `hi_is_ctrl` and `!hi_is_ctrl` sub-cases are exercised by varying `(q1, q2)` ordering.
- **5 qubits** (dim 32, n_tiles=4): Tests Case 2 and Case 3. Case 3 requires both qubits
  in tile coordinates (bits 3+), which is possible with qubits 3 and 4.
- **6 qubits** (dim 64, n_tiles=8): Full exercise of Case 3 with `quarter_nt = 2`.

**Coverage assessment**: Good. All three cases are exercised. The diagonal-block special cases
(`is_diag`) are tested because the iteration includes `b_tc == b_tr`. The `lower_01`
and `!lower_01` paths are tested because different qubit orderings produce different tile
relationships.

**Limitation**: With only up to 6 qubits, `quarter_nt` is at most 2 for Case 3, meaning
the diagonal-skip logic for Pairs C+D is only tested with `b_tr in {0, 1}`. This is enough
to trigger `is_diag` (b_tr==b_tc==0) and `!is_diag` (b_tr=1, b_tc=0), but there is no
deep nesting. For higher confidence, testing at 7-8 qubits would be ideal, though the
O(dim^3) reference computation makes this slow.

---

## 5. Remaining Concerns

### 5.1 The redundant mirror pass (Minor)

Lines 1164-1170 and lines 1251-1256 perform explicit upper-triangle mirroring of diagonal
tiles after all swap operations. In the `hi_is_ctrl` `is_diag` path, the
`dtile_read`/`dtile_write` functions already maintain the invariant that the lower triangle
is the source of truth. The mirror pass writes `t11[lr, lc] = conj(t11[lc, lr])` for
`lc > lr`, which is a no-op if `dtile_write` already did its job.

In the `!hi_is_ctrl` `is_diag` path, the mirror pass covers both `t00` and `t11`. Here it is
potentially necessary for `t00` because some writes to `t00` go through `dtile_write` (which
updates the canonical position) but the upper triangle may have stale values from before the
swap. The mirror pass ensures consistency.

**Recommendation**: Keep the mirror passes but add a one-line comment explaining they are
needed to flush stale upper-triangle values that `dtile_write` does not clean up (it writes
the canonical lower-tri position but does not update the non-canonical upper-tri copy).

Actually, looking again at `dtile_write`:

```c
static inline void dtile_write(cplx_t *tile, gate_idx_t lr, gate_idx_t lc, cplx_t val) {
    if (lr >= lc) {
        tile[elem_off(lr, lc)] = val;
    } else {
        tile[elem_off(lc, lr)] = conj(val);
    }
}
```

This writes ONLY to the lower-triangle position. It does NOT update the upper-triangle copy.
So if anyone reads `tile[elem_off(lr, lc)]` directly (without going through `dtile_read`),
they will get a stale value for `lr < lc`. The mirror pass is therefore **necessary** to
maintain the invariant that "diagonal tiles store both triangles consistently" which is relied
upon by other code paths (e.g., Case 1 reads elements directly without `dtile_read`).

**Corrected assessment**: The mirror passes are necessary and correct. The comment should be
strengthened from "for consistency" to "required: dtile_write only updates the canonical
lower-triangle position; this pass restores the upper-triangle copies."

### 5.2 The `!lo_10 && lo_11` impossibility in Pair C (Trivial)

The Case 3 Pairs C+D code at line 1516 (`else`) handles the `!lo_10` case, which implies
`!lo_11` because `tc11 > tc10`. There is no explicit assertion or comment for this invariant.
Adding `/* !lo_10 implies !lo_11 since tc11 > tc10 */` would make the code self-documenting.

### 5.3 OpenMP thread safety in Case 2 `is_diag` path

In Case 2, the outer parallel loop is over `i_tr`. For a given `i_tr`, the inner loop over
`i_tc` is sequential. Since `is_diag` only occurs when `i_tc == i_tr`, and each `i_tr` value
is processed by exactly one thread, there is no data race on the diagonal tiles. Correct.

For off-diagonal tile-groups, two different `i_tr` values never share a tile because the
`insertBit0` mapping produces disjoint tile-index sets for different `i_tr` values. The only
shared tiles would be column tiles (e.g., `tc0` could be the same for different `i_tr` values),
but the write targets are always row-keyed tiles (`t00`, `t10`, `t11` indexed by `tr0`/`tr1`),
and `t01` is indexed by both `tr0` and `tc1`. Two threads with different `i_tr` but the same
`i_tc` could both access `t01 = T(tr0_A, tc1)` and `t01 = T(tr0_B, tc1)` -- these are
different tiles since `tr0_A != tr0_B`. No race.

**Assessment**: Thread-safe.

---

## 6. Recommendations

### Must-do

None. The code is correct and passes all tests.

### Should-do

1. **Strengthen mirror-pass comments.** Change "for consistency" to "required: dtile_write
   only updates the canonical lower-triangle copy" at lines 1159-1163 and 1250.

2. **Add invariant comment for Pair C/D both-upper sub-case.** At line 1516, add
   `/* !lo_10 implies !lo_11: tc11 = tc10 + incr_tgt > tc10 > tr00 */`.

3. **Unify naming in Case 3 tile-swap loops.** Use `t_src`/`t_dst` or just `tA`/`tB`
   consistently across all 6 pairs, rather than mixing `tA`/`tB` for some and `tCa`/`tCb`
   for others.

### Could-do

4. **Extract whole-tile swap into a helper.** The pattern:
   ```c
   for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
       cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
   }
   ```
   appears 8 times in Case 3 (and similarly in `swap_tiled`). A `tile_swap(cplx_t *a, cplx_t *b)`
   helper would reduce code volume and create a single optimization point (e.g., for SIMD
   intrinsics or `memcpy`-based swap with a temp buffer).

5. **Extract conjugate-transpose tile swap into a helper.** The pattern:
   ```c
   for (lr ...) for (lc ...) {
       cplx_t val_a = tA[elem_off(lr, lc)];
       cplx_t val_b = conj(tB[elem_off(lc, lr)]);
       tA[elem_off(lr, lc)] = val_b;
       tB[elem_off(lc, lr)] = conj(val_a);
   }
   ```
   appears 4 times in Case 3. A `tile_swap_conj(cplx_t *lower, cplx_t *upper_stored)` helper
   would eliminate repetition.

6. **Consider raising MAXQUBITS_TILED to 7.** This would give `n_tiles = 16`, `quarter_nt = 4`
   for Case 3, and exercise the triangular iteration more deeply. The reference computation
   cost is O(128^3) = ~2M complex multiply-adds, which should run in under a second.

---

## Summary Table

| Aspect | Case 1 | Case 2a (hi_is_ctrl) | Case 2b (!hi_is_ctrl) | Case 3 |
|--------|--------|---------------------|----------------------|--------|
| Correctness | OK | OK | OK | OK |
| Diagonal handling | Implicit (full tile) | dtile + mirror | read-all-write-all + mirror | skip (aliasing) |
| Upper-tri handling | N/A | lower_01 branch | lower_01 branch | per-pair lo checks |
| Vectorizable | No (scattered) | No (scattered) | No (scattered) | Yes (flat loops) |
| Lines | 64 | 116 | 158 | 226 |
| Complexity | Low | Medium | High | Medium |
