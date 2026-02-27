# State Module — Technical Specification

**Status:** Complete

---

## Overview

The state module represents n-qubit quantum systems in three storage formats:

| Type | Mathematical object | Use case |
|------|---------------------|----------|
| `PURE` | State vector ∈ ℂ^(2^n) | Pure quantum states |
| `MIXED_PACKED` | Density matrix, Hermitian packed | Mixed states; LAPACK-compatible |
| `MIXED_TILED` | Density matrix, blocked layout | Mixed states; cache-optimised gate operations |

All public functions dispatch at runtime on `state->type`. Type-specific implementations (`state_pure_*`, `state_packed_*`, `state_tiled_*`) are declared in `state.h` and called only through the dispatch layer or the type-specific files.

---

## Directory Layout

```
include/
  qlib.h              Top-level convenience header; includes state.h and gate.h
  state.h             Public API, type definitions, tile constants, type-specific function declarations
  q_types.h           Foundational types: cplx_t, dim_t, qubit_t, error codes, utility macros
  utils.h             Print macros: vprint, mprint_packed, mprint_tiled
src/
  state/
    state.c           Dispatch layer — routes all public API calls by state->type
    state_pure.c      PURE state implementations
    state_packed.c    MIXED_PACKED state implementations
    state_tiled.c     MIXED_TILED state implementations
test/
  state/
    test_state.c        Test runner + dispatch layer tests (2 tests)
    test_state_pure.c   PURE unit tests (10 tests)
    test_state_packed.c MIXED_PACKED unit tests (13 tests)
    test_state_tiled.c  MIXED_TILED unit tests (16 tests)
```

---

## Dependencies (Inputs)

| Dependency | Source | Notes |
|------------|--------|-------|
| CBLAS (`cblas_zcopy`) | vecLib (macOS) / OpenBLAS ILP64 (Linux) | Used in `state_cp`; OpenBLAS first build takes 10–30 min |
| LAPACK packed format | Same BLAS library | MIXED_PACKED layout is directly usable with `zhpmv`, `zhpr`, `zhpevd` |
| `q_types.h` | `include/` | Provides `cplx_t`, `dim_t`, `qubit_t` |

---

## Build Integration

The state module sources are compiled into the shared library target **`q`**, defined in `src/CMakeLists.txt`.

**Library target:** `q` (shared, `libq.so` / `libq.dylib`)

**State sources compiled into `q`:**
```cmake
src/state/state.c
src/state/state_pure.c
src/state/state_packed.c
src/state/state_tiled.c
```

**Key compile definitions passed to `q`:**
- `LOG_TILE_DIM=${LOG_TILE_DIM}` — set by CMake from build type; 3 (Debug) or 5 (Release).
- `OMP_THRESHOLD=4` — Debug-only override that lowers the minimum dimension for OpenMP parallelization so multi-threaded paths are exercised at small qubit counts.

**Test target:** `test_state` (executable), defined in `test/CMakeLists.txt`. Fetches Unity v2.6.1 via `FetchContent`. Links against `q`. Compile definitions: `UNITY_INCLUDE_DOUBLE`, `UNITY_SUPPORT_64`, `LOG_TILE_DIM=${LOG_TILE_DIM}`.

**Installed artifacts** (via `cmake --install`):
- `libq.so` / `libq.dylib` → `${CMAKE_INSTALL_LIBDIR}`
- `include/state.h`, `include/q_types.h`, `include/qlib.h` → `${CMAKE_INSTALL_INCLUDEDIR}`

**Run state tests:**
```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug --tests-regex test_state
```

---

## Types

### `state_t`

```c
typedef struct state {
    state_type_t type;   // Storage format; must be set before any other call
    cplx_t       *data;  // Heap-allocated, 64-byte-aligned; owned by struct
    qubit_t       qubits;
} state_t;
```

Invariants: `data != NULL` iff storage is allocated and valid. `data == NULL` means uninitialized, freed, or OOM.

### `state_type_t`

```c
typedef enum { PURE, MIXED_PACKED, MIXED_TILED } state_type_t;
```

`MIXED` is a deprecated alias for `MIXED_PACKED` (defined in `state.h`).

### Types from `q_types.h`

| Type | Underlying type | Notes |
|------|-----------------|-------|
| `cplx_t` | `__LAPACK_double_complex` (macOS) / `openblas_complex_double` (Linux) | 16 bytes |
| `dim_t` | `__LAPACK_int` (macOS) / `blasint` (Linux) | 64-bit on both platforms (ILP64); overflow-safe for any memory-feasible qubit count |
| `qubit_t` | `unsigned char` | Max 255 qubits |

`q_types.h` also defines:
- `POW2(a, T)`, `SETREAL(a, b)`, `SETIMAG(a, b)`, `SWAP(a, b, T)` — used throughout gate implementations.
- `MIN(a, b)`, `MAX(a, b)` — guarded by `#ifndef` to avoid conflicts with system headers.

### Error codes (`qs_error_t`)

Defined in `q_types.h` as a `typedef enum`:

| Code | Value | Meaning |
|------|-------|---------|
| `QS_OK` | 0 | Success |
| `QS_ERR_NULL` | -1 | Null pointer argument |
| `QS_ERR_OOM` | -2 | Out of memory |
| `QS_ERR_QUBIT` | -3 | Invalid qubit index |
| `QS_ERR_TYPE` | -4 | Invalid state type for operation |
| `QS_ERR_FILE` | -5 | File I/O error |
| `QS_ERR_FORMAT` | -6 | Invalid file format |
| `QS_ERR_PARAM` | -7 | Invalid parameter value |

Note: the public API currently uses `assert` and OOM-via-NULL rather than returning `qs_error_t`. The enum is defined for future use or by other modules.

---

## Constants

Defined in `state.h` and set by CMake at compile time via `-DLOG_TILE_DIM=<n>`:

```c
#ifndef LOG_TILE_DIM
#define LOG_TILE_DIM 5          // Fallback only; value is injected by CMake
#endif
#define TILE_DIM  (1 << LOG_TILE_DIM)   // 8 (Debug) or 32 (Release)
#define TILE_SIZE (TILE_DIM * TILE_DIM) // 64 (Debug) or 1024 (Release)
```

**Do not set `LOG_TILE_DIM` manually.** CMake enforces correct values per build type (3 for Debug, 5 for Release). Valid range enforced by CMake: 3–8.

**Rationale for release default (LOG_TILE_DIM = 5, TILE_DIM = 32):**
- One tile = 32 × 32 × 16 bytes = 16 KB → fits in L1 cache (~32–64 KB)
- Four tiles = 64 KB → fits in L2 cache (~256 KB), enabling cross-tile operations without thrashing

---

## Memory Layouts

### PURE — State Vector

Stores the 2^n complex amplitudes of |ψ⟩ = Σᵢ cᵢ|i⟩ as a flat array.

**Basis encoding:** little-endian binary. Index i = 5 (binary `101`) encodes |q₀=1, q₁=0, q₂=1⟩.

**Storage:**
```
data[i] = cᵢ      i ∈ {0, ..., 2^n − 1}
```

**Array length:** `2^n`

---

### MIXED_PACKED — Hermitian Packed Lower Triangle

A density matrix ρ ∈ ℂ^(d×d) (d = 2^n) is Hermitian: ρ† = ρ. Only the lower triangle is stored, in column-major order. This matches LAPACK's packed Hermitian format (`uplo='L'`) and is directly usable with `zhpmv`, `zhpr`, `zhpevd`.

**Layout for d = 4:**
```
ρ = [ρ₀₀  ·    ·    ·  ]     Packed array (column-major, lower triangle):
    [ρ₁₀  ρ₁₁  ·    ·  ]     index:  0     1     2     3     4     5     6     7     8     9
    [ρ₂₀  ρ₂₁  ρ₂₂  ·  ]     value: ρ₀₀  ρ₁₀  ρ₂₀  ρ₃₀  ρ₁₁  ρ₂₁  ρ₃₁  ρ₂₂  ρ₃₂  ρ₃₃
    [ρ₃₀  ρ₃₁  ρ₃₂  ρ₃₃]
```

**Index formula** for element (row, col) with row ≥ col:
```c
idx = col * d - col * (col + 1) / 2 + row
```

**Upper triangle:** `ρ[row, col] = conj(ρ[col, row])` for row < col — never stored, computed on access.

**Array length:** `d * (d + 1) / 2`

---

### MIXED_TILED — Blocked Layout

The density matrix is partitioned into TILE_DIM × TILE_DIM tiles. Hermiticity is exploited **at the tile level only**: upper-triangular tile blocks (tile_row < tile_col) are not stored. The lower-triangular tiles (tile_row ≥ tile_col) are stored in row-major tile order.

**Every stored tile, including diagonal tiles, is a full TILE_DIM × TILE_DIM block of TILE_SIZE elements.** There is no secondary triangular packing within a tile. A diagonal tile T(tr, tr) stores all TILE_DIM² elements — both the portion below and the portion above the matrix diagonal that falls within that tile's row/column range.

This layout keeps gate operands contiguous in memory. The gate module accesses tiles directly by pointer arithmetic rather than element-by-element indexing.

**Tile grid for d = 2·TILE_DIM (n = LOG_TILE_DIM + 1):**
```
Tile grid (2×2, lower triangle):
  T(0,0)  ·
  T(1,0)  T(1,1)

Storage order: T(0,0), T(1,0), T(1,1)
```

**Derived quantities** (computed, never stored):
```c
dim_t n_tiles = (dim + TILE_DIM - 1) / TILE_DIM;  // ceil(2^n / TILE_DIM)
```

**Index formula** for element (row, col) with row ≥ col:
```c
tile_row    = row / TILE_DIM;
tile_col    = col / TILE_DIM;
tile_offset = tile_row * (tile_row + 1) / 2 + tile_col;
local_row   = row % TILE_DIM;
local_col   = col % TILE_DIM;
idx         = tile_offset * TILE_SIZE + local_row * TILE_DIM + local_col;
```

**Upper triangle:** same as packed — access via conjugate of transposed coordinates.

**Array length:** `n_tiles * (n_tiles + 1) / 2 * TILE_SIZE`

**Padding:** When `dim < TILE_DIM` (n < LOG_TILE_DIM), a single full tile is allocated; only the `dim × dim` corner is physical. Padding elements must remain zero.

**LAPACK incompatibility warning:** Intra-tile layout is row-major. Do not pass tile pointers to LAPACK routines that expect column-major packed storage.

---

## Public API

All functions in `state.h` dispatch on `state->type`. Initialize `type` before any call.

### `state_len`
```c
dim_t state_len(const state_t *state);
```
Returns the number of `cplx_t` elements in `data` for the given type and qubit count.

### `state_init`
```c
void state_init(state_t *state, qubit_t qubits, cplx_t **data);
```
Allocates storage. Frees any previous `state->data`. If `data == NULL`, allocates zero-initialized, 64-byte-aligned storage. If `data != NULL`, transfers ownership (`*data` set to NULL); copies to aligned storage first if necessary.

### `state_free`
```c
void state_free(state_t *state);
```
Frees `data`, sets `data = NULL` and `qubits = 0`. Safe to call on NULL or already-freed state.

### `state_plus`
```c
void state_plus(state_t *state, qubit_t qubits);
```
Initializes to the uniform superposition:
- **PURE:** |+⟩^⊗n, all amplitudes = 1/√d
- **MIXED_PACKED / MIXED_TILED:** ρ = |+⟩⟨+|^⊗n, all matrix elements = 1/d (including off-diagonals)

Note: this is a rank-1 pure state embedded in density matrix form, **not** the maximally mixed state I/d.

### `state_cp`
```c
state_t state_cp(const state_t *state);
```
Returns a deep copy with an independent allocation. Allocates 64-byte-aligned storage, then copies with `cblas_zcopy`. On OOM, the returned state has `data == NULL`. Requires `state->data != NULL`; calling on an uninitialized or freed state triggers an `assert`.

### `state_get` / `state_set`
```c
cplx_t state_get(const state_t *state, dim_t row, dim_t col);
void   state_set(state_t *state, dim_t row, dim_t col, cplx_t val);
```
Element access with Hermitian symmetry handling. For PURE, `col` must be 0; passing a non-zero `col` triggers an `assert` (crash in debug builds, undefined behavior if asserts are disabled). For mixed types, accessing the upper triangle (row < col) returns/stores the conjugate of the transposed element.

**MIXED_TILED diagonal tile mirroring:** `state_tiled_set` applies additional mirroring for diagonal tiles (where `row / TILE_DIM == col / TILE_DIM`): setting a lower-triangle element also writes the conjugate into the corresponding upper-triangle position within the same tile, so the full tile stays consistent for gate code that reads both sides directly.

### `state_print`
```c
void state_print(const state_t *state);
```
Prints to stdout. Delegates to `mprint_packed` / `mprint_tiled` macros in `utils.h`.

---

## Memory Ownership

- `state_t` owns its `data` pointer after initialization.
- `state_init(..., &data)` transfers ownership; caller's pointer is set to NULL.
- No shared ownership. Each buffer has exactly one owner.
- Calling `state_init()` on an already-initialized state frees the previous allocation.

---

## Error Handling

| Error class | Handling |
|-------------|----------|
| Programmer error (NULL state, invalid type) | `assert()` — crash immediately |
| Out-of-memory | `state->data = NULL`, message to stderr |

After `state_init()` or `state_plus()`, check `state->data != NULL` before use.

---

## Memory Footprint

Each `cplx_t` is 16 bytes. Release TILE_DIM = 32; Debug TILE_DIM = 8.

| Qubits | dim | PURE | MIXED_PACKED | MIXED_TILED |
|--------|-----|------|--------------|-------------|
| 1 | 2 | 32 B | 48 B | 16 KB |
| 2 | 4 | 64 B | 160 B | 16 KB |
| 5 | 32 | 512 B | 8.4 KB | 16 KB |
| 8 | 256 | 4 KB | 526 KB | 576 KB |
| 10 | 1024 | 16 KB | 8.4 MB | 8.6 MB |
| 12 | 4096 | 64 KB | 134 MB | 136 MB |

MIXED_TILED figures assume Release (TILE_DIM = 32). In Debug (TILE_DIM = 8), the minimum allocation is 1 KB (one tile) and overhead differs significantly for small qubit counts.

MIXED_TILED wastes space for small systems (n < LOG_TILE_DIM) but converges to MIXED_PACKED overhead for n ≥ 8.

---

## Platform Notes

`dim_t` is 64-bit on both platforms. On Linux, OpenBLAS must be compiled with ILP64 (handled by CMake; first build takes 10–30 min). On macOS, vecLib provides ILP64 natively.

---

## Open Questions

1. ~~**`state_get`/`state_set` for PURE with col ≠ 0:**~~ **RESOLVED.** Both `state_get` (state.c:183) and `state_set` (state.c:201) contain `assert(col == 0 && "col must be 0 for PURE states")`. A non-zero `col` triggers an assert crash in debug builds; behavior with asserts disabled is undefined. Note: the `state.h` comment says "col is ignored" for PURE — this is misleading and contradicts the implementation. The assert is the authoritative behavior.

2. ~~**Alignment contract for `state_init` ownership-transfer path:**~~ **RESOLVED.** The alignment threshold is 64 bytes (`STATE_ALIGNMENT = 64`, state.c:19), checked as `((uintptr_t)*data & 63) == 0`. If the buffer is not 64-byte aligned: new aligned storage is allocated, the data is `memcpy`'d, and the original buffer is `free()`'d. The caller's pointer is set to NULL in both aligned and non-aligned paths (state.c:96).

3. **Thread safety:** No statement is made about whether any `state_t` operation is safe to call concurrently on distinct state objects, or on the same state object from multiple threads. No synchronization is present in any implementation file.

4. ~~**`MIXED` deprecated alias:**~~ **RESOLVED.** `MIXED` is a preprocessor `#define` to `MIXED_PACKED` (state.h:50). It is resolved at compile time and never appears as a separate `case` in any switch statement. No removal timeline is documented.

5. ~~**Padding element invariant enforcement:**~~ **RESOLVED.** `state_init` calls `state_alloc_aligned` which uses `memset(ptr, 0, size)` (state.c:33) — this zero-initializes all allocated bytes including padding. In the ownership-transfer path, when a non-aligned buffer is copied into new aligned storage, the new allocation is also zero-initialized before the `memcpy`, so padding bytes are zero. If the caller passes an already-aligned buffer, it is taken directly with no zeroing — the caller is responsible for ensuring padding is zero in that path.

6. ~~**`state_plus` for MIXED_TILED with padding:**~~ **RESOLVED.** `state_tiled_plus` (state_tiled.c:124–137) explicitly checks `if (state->qubits < LOG_TILE_DIM)` and loops only over the `dim × dim` physical corner. Padding elements remain zero from `state_init`'s `memset`.

7. **Memory footprint table accuracy for Debug builds:** The MIXED_TILED column was computed assuming TILE_DIM = 32 (Release). The table should either be annotated for which build type it applies to, or split by build type for clarity. (Partially addressed in the footnote, but the raw numbers in the table are potentially misleading for Debug users.)
