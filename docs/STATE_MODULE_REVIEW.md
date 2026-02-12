# State Module Code Review

**Reviewer:** Claude (Torvalds-style)
**Date:** 2026-02-06
**Files reviewed:**
- `/Users/timo/Code/qlib/include/state.h`
- `/Users/timo/Code/qlib/include/q_types.h`
- `/Users/timo/Code/qlib/include/utils.h`
- `/Users/timo/Code/qlib/src/state/state.c`
- `/Users/timo/Code/qlib/src/state/state_pure.c`
- `/Users/timo/Code/qlib/src/state/state_packed.c`
- `/Users/timo/Code/qlib/src/state/state_tiled.c`
- `/Users/timo/Code/qlib/test/src/test_state.c`
- `/Users/timo/Code/qlib/test/src/test_state_pure.c`
- `/Users/timo/Code/qlib/test/src/test_state_packed.c`
- `/Users/timo/Code/qlib/test/src/test_state_tiled.c`

**Test results:** 41/41 PASS

---

## Verdict

This is solid, well-structured code. The architecture is clean, the dispatch pattern is simple and
appropriate, the math is correct, and it shows evidence of careful iteration (overflow fixes, alignment
contracts, padding handling). Acceptable as-is for a research library. What follows are the things
that would make it better, and a few things that could bite you later.

---

## Critical Issues

### 1. `state_tiled_len()` overflow in tile count arithmetic

In `/Users/timo/Code/qlib/src/state/state_tiled.c`, line 104:

```c
return (n_tiles * (n_tiles + 1) >> 1) * TILE_SIZE;
```

Unlike `state_packed_len()` which carefully uses `(dim >> 1) * (dim + 1)` to avoid intermediate
overflow, here `n_tiles * (n_tiles + 1)` can overflow `dim_t` before the right-shift. With
`dim_t = long` (64-bit), this is not a practical concern today since you would run out of RAM
long before `n_tiles` approaches 2^31, but `state_packed_len()` went to the trouble of preventing
this, and `state_tiled_len()` did not. The inconsistency is a code smell. Either both functions
should handle it, or neither should (and the assert is sufficient). Pick one policy.

### 2. `tile_offset()` same overflow risk

In `/Users/timo/Code/qlib/src/state/state_tiled.c`, line 52:

```c
return tr * (tr + 1) / 2 + tc;
```

No `size_t` cast. Again, this only matters for absurdly large tile counts that are impossible in
practice, but `packed_idx()` has these casts and `tile_offset()` does not. Be consistent.

### 3. `SETREAL`/`SETIMAG` are strict aliasing violations

In `/Users/timo/Code/qlib/include/q_types.h`, lines 77-78:

```c
#define SETREAL(a, b)  (*((double*) &(a)) = (b))
#define SETIMAG(a, b)  (*(((double*) &(a)) + 1) = (b))
```

This assumes `cplx_t` is laid out as `{double real, double imag}` and accesses it through a
`double*` cast. This works in practice on every compiler you will ever use because `cplx_t` is
a C99 `_Complex double` or equivalent struct, and C guarantees this layout for complex types
(C11 6.2.5p13). However, the cast from `__LAPACK_double_complex*` to `double*` technically
depends on `__LAPACK_double_complex` being layout-compatible with `double _Complex`. If Apple
ever changes their opaque type, this breaks silently. A safer approach would be:

```c
#define SETREAL(a, b)  (((double _Complex*) &(a))[0] = (b) + cimag(a) * I)
```

Or just use C99 `__real__` and `__imag__` GCC extensions which both Clang and GCC support. This
is not urgent, but it is the kind of thing that bites you during a platform migration and you
spend three days tracking down.

### 4. `state_get` asserts `col == 0` for PURE but the API accepts it

In `/Users/timo/Code/qlib/src/state/state.c`, line 183:

```c
assert(col == 0 && "col must be 0 for PURE states");
```

This is a debug-only check. In release builds (`-DNDEBUG`), `col` is silently ignored for PURE
states, which is the documented behavior. But any caller who passes a non-zero `col` to a PURE
state has a bug -- they are confused about what type they are working with. The assert is the
right call here. No change needed, just flagging that this is a conscious design decision.

---

## Improvements

### 5. `state_alloc_aligned()` uses `memset` on complex type

In `/Users/timo/Code/qlib/src/state/state.c`, lines 31-35:

```c
cplx_t *ptr = aligned_alloc(STATE_ALIGNMENT, size);
if (ptr) {
    memset(ptr, 0, size);
}
```

`memset(ptr, 0, size)` sets all bytes to zero. For IEEE 754 doubles, all-zero-bytes is
`+0.0`, and for a complex type stored as two doubles, this gives `0.0 + 0.0i`. This is
correct on every IEEE 754 platform. But if you want to be pedantic about it, you could
document *why* this is safe. A `calloc`-like call would also work here but would not give
you alignment guarantees. What you have is fine.

### 6. `cblas_zcopy` in `state_cp` -- overkill for a memcpy

In `/Users/timo/Code/qlib/src/state/state.c`, line 153:

```c
cblas_zcopy(len, state->data, 1, out.data, 1);
```

With stride 1, `cblas_zcopy` is just a memcpy with extra function-call overhead (BLAS dispatch,
type checking). For a simple contiguous copy, `memcpy(out.data, state->data, len * sizeof(cplx_t))`
would be faster and clearer. The BLAS call only buys you something when you need non-unit strides,
which you do not here.

That said, this is a one-time allocation cost, not a hot path, so it does not matter for
performance. It is a clarity issue, not a performance issue.

### 7. `state_pure_plus()` loop should be a `memset`-style fill or BLAS call

In `/Users/timo/Code/qlib/src/state/state_pure.c`, lines 40-42:

```c
for (dim_t i = 0; i < len; ++i) {
    state->data[i] = amplitude;
}
```

For a pure state with 2^32 elements, this is a scalar loop filling ~64 GB. The compiler *may*
auto-vectorize this, but you are leaving performance on the table. A `memset`-style approach
or BLAS `zset` (if available) would be better for large states. Same applies to
`state_packed_plus()` and the bulk branch of `state_tiled_plus()`.

For the target qubit counts in the PRD (max 32 qubits for pure, 16 for mixed), this
initialization cost is dwarfed by gate application time, so this is a minor point. But it is
the kind of thing that separates "works" from "works fast."

### 8. `state_init` reinit path calls `free()` then `state_len()`

In `/Users/timo/Code/qlib/src/state/state.c`, lines 71-76:

```c
if (state->data) {
    free(state->data);
    state->data = NULL;
}
state->qubits = qubits;
```

When reinitializing, you free the old data, update `qubits`, and then compute `state_len()`
with the new qubit count. This is correct. But there is a subtle thing: if the new allocation
fails, the old data is gone and you have a state with `data == NULL` and `qubits == 0`. This
is documented and handled correctly. Good.

### 9. `#pragma once` AND `#ifndef` guard in `q_types.h`

In `/Users/timo/Code/qlib/include/q_types.h`, lines 1-4:

```c
#pragma once

#ifndef QTYPES_H
#define QTYPES_H
```

Pick one. `#pragma once` is sufficient and supported by every compiler that matters (GCC, Clang,
MSVC). The `#ifndef` guard is redundant. It is not a bug, just noise. Every other header in the
project uses only `#ifndef`, so either drop `#pragma once` here for consistency, or add it
everywhere. Consistency matters more than the choice itself.

### 10. `test.h` uses `index_t` which is not defined

In `/Users/timo/Code/qlib/test/include/test.h`, lines 59 and 76:

```c
for (index_t _i = 0; _i < (num_elements); ++_i) { \
```

`index_t` is not defined in `test.h`, `q_types.h`, or any included header. This compiles because
`TEST_ASSERT_EQUAL_DOUBLE_ARRAY_TOL` and `TEST_ASSERT_EQUAL_DOUBLE_ARRAY_ABS` are never
instantiated in the state tests. But this is a landmine -- the first person to use these macros
will get a cryptic compiler error. Meanwhile, the complex array macro on line 117 correctly uses
`size_t`. Fix `index_t` to `size_t` or `dim_t`.

### 11. `znprint` output format quirk

In `/Users/timo/Code/qlib/include/utils.h`, line 39:

```c
cimag(x) >= 0 ? printf("+i%.7f", cimag(x)) : printf("-i%.7f", fabs(cimag(x)));
```

This prints `-i0.0000000` when `cimag(x)` is negative zero (`-0.0`), which can happen with
complex arithmetic. Negative zero compares as `>= 0` per IEEE 754, so this is actually fine --
you get `+i0.0000000` for negative zero. Never mind, this is correct. But the `+i`/`-i` notation
is unconventional (most tools print `0.5+0.3i` not `0.5+i0.3`). This is a cosmetic choice, not a
bug, but may confuse someone comparing output against numpy or Mathematica.

### 12. Print macros expand to enormous code

`mprint_packed` and `mprint_tiled` in `/Users/timo/Code/qlib/include/utils.h` are macros that
expand to ~20 lines of code at every call site. These are only used in print functions, so code
size is not a practical concern, but they should be actual functions. The only reason they are
macros is to support `_Generic` dispatch via `nprint`, but you could just make them functions
that take `cplx_t*` and call `znprint` directly.

### 13. Missing test: `state_pure_len(0)`

The packed tests include `test_packed_zero_qubit` which verifies `state_packed_len(0) == 1`.
But there is no equivalent test for `state_pure_len(0)`. This should return 1 (a 0-qubit system
has a single amplitude), and `POW2(0, dim_t)` does give 1, so it works. But test it.

### 14. Missing test: allocation failure paths

None of the tests verify behavior when `malloc`/`aligned_alloc` fails. The code handles this
(sets `data = NULL`, prints to stderr), but it is never tested. This is understandable -- mocking
`malloc` in C is painful -- but worth noting.

### 15. The spec documents conversion functions that do not exist

Section 11 of `/Users/timo/Code/qlib/docs/STATE_MODULE.md` lists `packed_to_full`,
`full_to_packed`, `tiled_to_full`, `full_to_tiled`. These functions do not exist in the codebase.
Either implement them or remove them from the spec. A spec that documents nonexistent functions
is worse than no spec at all, because it wastes the time of anyone trying to use them.

---

## Specific Feedback

### Architecture: Dispatch Pattern

The `state.c` dispatch layer using `switch(state->type)` is clean and appropriate for 3 types.
If the type count grows beyond ~5, consider a vtable (array of function pointers), but for now
this is the right level of abstraction. Do not over-engineer it.

The decision to store only `{type, data, qubits}` and compute everything else (dim, n_tiles,
len) on demand is correct. No stale cached values, no synchronization issues. The cost is a
shift and some arithmetic, which is negligible compared to any actual gate operation.

### `state_alloc_aligned`: Good

```c
static cplx_t *state_alloc_aligned(dim_t len) {
    size_t size = len * sizeof(cplx_t);
    size = (size + STATE_ALIGNMENT - 1) & ~(size_t)(STATE_ALIGNMENT - 1);
    cplx_t *ptr = aligned_alloc(STATE_ALIGNMENT, size);
    if (ptr) {
        memset(ptr, 0, size);
    }
    return ptr;
}
```

Correctly rounds up to alignment boundary (required by `aligned_alloc`'s contract). 64-byte
alignment is correct for AVX-512 and cache line alignment. The rounding ensures the `size`
argument is a multiple of `STATE_ALIGNMENT`. This is correct.

One edge case: `len == 0` gives `size == 0`, and `aligned_alloc(64, 0)` is implementation-defined
(C17 7.22.3.1). On macOS, this returns a valid pointer that must be freed; on Linux, it may
return NULL. Since `state_packed_len(0)` returns 1, this path is not currently reachable with 0
elements, but it is worth being aware of.

### `state_init` alignment contract: Good

```c
if (((uintptr_t)*data & (STATE_ALIGNMENT - 1)) == 0) {
    state->data = *data;
} else {
    // copy into aligned storage
}
```

This is exactly right. Take aligned buffers directly, copy unaligned ones. The ownership transfer
semantics (nullifying `*data`) are clean and prevent use-after-transfer bugs. The fact that
unaligned buffers get copied and freed is transparent to the caller.

### `packed_idx` overflow prevention: Good

```c
static inline dim_t packed_idx(dim_t dim, dim_t row, dim_t col) {
    return (dim_t)((size_t)col * dim - (size_t)col * (col + 1) / 2 + row);
}
```

The `size_t` casts prevent intermediate overflow of `col * dim` when `dim_t` is 32-bit. With
64-bit `dim_t` (ILP64), this is belt-and-suspenders. Fine.

### `state_packed_len` special case: Correct but fragile

```c
if (qubits == 0) {
    return 1;
}
return (dim >> 1) * (dim + 1);
```

The special case for `qubits == 0` is necessary because `dim = 1` is odd, making
`(dim >> 1) * (dim + 1) = 0 * 2 = 0`, which is wrong. The comment explains this clearly.
The formula `(dim >> 1) * (dim + 1)` avoids overflow by exploiting the fact that
`dim = 2^n` is always even for `n >= 1`. Correct.

### `state_tiled_plus` padding handling: Good

```c
if (state->qubits < LOG_TILE_DIM) {
    for (dim_t row = 0; row < dim; ++row) {
        for (dim_t col = 0; col < dim; ++col) {
            state->data[row * TILE_DIM + col] = prefactor;
        }
    }
} else {
    for (dim_t i = 0; i < len; ++i) {
        state->data[i] = prefactor;
    }
}
```

For small systems, this correctly writes only the physical `dim x dim` corner and leaves
padding as zero (from `state_alloc_aligned`'s `memset`). For large systems where `dim` is a
multiple of `TILE_DIM`, there is no padding and bulk fill is correct. This handles the edge
case that would otherwise put nonzero values in padding positions, which could cause subtle
bugs in gate operations that iterate over full tiles.

Note: the small-system branch writes to positions `row * TILE_DIM + col`, which is the
row-major intra-tile layout (correct for the single tile at offset 0). Good.

But there is a missing case: what if `qubits >= LOG_TILE_DIM` but `dim` is NOT a multiple
of `TILE_DIM`? This cannot happen because `dim = 2^qubits` and `TILE_DIM = 2^LOG_TILE_DIM`,
so `qubits >= LOG_TILE_DIM` implies `dim >= TILE_DIM` and `dim` is a power of 2 that is a
multiple of `TILE_DIM` (also a power of 2). The logic is sound -- but deserves a comment
explaining *why* the bulk fill is safe (no padding when `qubits >= LOG_TILE_DIM`).

### `state_free` NULL safety: Good

```c
void state_free(state_t *state) {
    if (!state) return;
    free(state->data);
    state->data = NULL;
    state->qubits = 0;
}
```

`free(NULL)` is a no-op per C standard, so the double-free safety comes from setting
`state->data = NULL` after free. Correct. Setting `qubits = 0` is a nice touch for
debugging.

### Test quality

The tests are thorough. The 41-test suite covers:
- Happy paths for all three storage types
- Edge cases (0 qubits, 1 qubit, tile boundaries)
- Hermitian symmetry (upper triangle returns conjugate)
- Ownership transfer semantics
- Double-free safety
- Deep copy independence
- Cross-tile access (8 qubits = 256x256 = 8 tiles per dimension)
- Minimal multi-tile case (6 qubits = 2 tiles per dimension)

Missing coverage:
- `state_pure_len(0)` (trivial but should be tested for completeness)
- Behavior when `state->type` is set to a garbage value (would trigger assert)
- Large qubit counts that stress the overflow guards
- Interaction between `state_init` reinit and ownership transfer (reinit with
  non-NULL data after previous init)

---

## Performance Analysis

### Hot path identification

The state module is not the hot path. Gate operations (in `qhipster.c`/`mhipster.c`) are where
99%+ of compute time is spent. The state module's job is to allocate correctly, provide correct
indexing, and get out of the way. It does this well.

### Allocation strategy

64-byte aligned allocation is correct for:
- AVX-512 (512-bit = 64-byte vectors)
- Cache line alignment (most x86 CPUs use 64-byte cache lines)
- Apple Silicon (128-byte cache lines on M1+ would benefit from 128-byte alignment, but 64
  is sufficient for correctness)

The `aligned_alloc` + `memset` pattern is one extra pass over memory compared to a hypothetical
`aligned_calloc`, but since this is initialization cost, it is irrelevant.

### Indexing overhead

`tiled_index()` computes: 2 shifts, 2 masks, 1 multiply, 2 adds, 1 shift-right (for
`tr*(tr+1)/2`). This is ~5 cycles on modern x86. For element-wise access via `state_get/set`,
this adds ~5 cycles per element. But `state_get/set` should NOT be used in gate inner loops --
the gate implementations should compute tile pointers once and then iterate within tiles using
simple pointer arithmetic. The per-element accessor is for debugging and testing only.

`packed_idx()` is simpler: 1 multiply, 1 multiply, 1 shift-right, 1 subtract, 1 add. Also ~5
cycles. Same observation applies.

If these accessors are being called in inner loops, that is a design bug in the caller, not in
the state module.

### Memory overhead of tiled format

For small systems, `MIXED_TILED` wastes significant memory:
- 1 qubit (2x2 matrix): stores 32x32 = 1024 elements instead of 3. That is 341x overhead.
- 4 qubits (16x16 matrix): stores 1024 elements instead of 136. That is 7.5x overhead.
- 5 qubits (32x32 matrix): stores 1024 elements instead of 528. That is 1.9x overhead.

At 8 qubits (256x256 matrix, 8 tiles per dimension), the overhead drops to about 9% (36 full
tiles vs. 36 tiles of 528 elements each in an ideal packed format). This is the right tradeoff
for the target workloads (8-16 qubits for mixed states), where cache locality during gate
operations far outweighs the memory overhead.

### `cblas_zcopy` vs `memcpy`

As noted above, `cblas_zcopy` with unit stride is functionally identical to `memcpy` but goes
through BLAS dispatch overhead. For a one-time copy during `state_cp`, this does not matter.
But for code clarity, `memcpy` is better.

---

## Summary of Action Items

| Priority | Item | File | Lines |
|----------|------|------|-------|
| Low | Fix `index_t` to `size_t` in test macros | `test/include/test.h` | 59, 76 |
| Low | Remove `#pragma once` from `q_types.h` (or add to all headers) | `include/q_types.h` | 1 |
| Low | Add overflow protection to `tile_offset()` and `state_tiled_len()` or document why it is not needed | `src/state/state_tiled.c` | 52, 104 |
| Low | Document why bulk fill is safe in `state_tiled_plus()` | `src/state/state_tiled.c` | 133 |
| Low | Remove or implement conversion functions listed in spec | `docs/STATE_MODULE.md` | 369-372 |
| Low | Add `state_pure_len(0)` test | `test/src/test_state_pure.c` | -- |
| Cosmetic | Replace `cblas_zcopy` with `memcpy` in `state_cp` | `src/state/state.c` | 153 |
| Cosmetic | Consider making print macros into functions | `include/utils.h` | 77-131 |

Nothing here is blocking. The code is correct, the tests pass, the architecture is sensible.
Ship it and move on to the gate module.

---

## Final Assessment

This is good, disciplined C code. It is simple where it should be simple (`state_pure.c`),
careful where it needs to be careful (`packed_idx` overflow, tiled padding), and the dispatch
layer is clean boilerplate. The test suite is thorough and covers the cases that matter.

The main risk in this codebase is not in the state module itself -- it is in whether the gate
modules use `state_get/set` in inner loops instead of computing tile/packed pointers directly.
If the gate implementations are as well-structured as the state module, this library should
perform well.

What I respect about this code:
- No unnecessary abstractions
- No "enterprise" patterns (no factories, no builders, no visitor pattern)
- The struct is 3 fields. Not 30.
- `state_free` is 4 lines. Not 40.
- The math is correct and the comments explain the non-obvious parts

What I would push back on in a review:
- The undeclared `index_t` in test.h is sloppy
- The spec documenting nonexistent functions is misleading
- The inconsistent overflow handling between packed and tiled is a code smell

But overall: this is code that knows what it is and does not try to be more. That is exactly
what good infrastructure code should look like.
