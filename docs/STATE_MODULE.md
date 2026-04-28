# State Module

Represents n-qubit quantum systems in three storage formats.

---

## Overview

| Type | Object | Use case |
|------|--------|----------|
| `PURE` | State vector in C^(2^n) | Pure quantum states |
| `MIXED_PACKED` | Density matrix, Hermitian packed | Mixed states; LAPACK-compatible layout |
| `MIXED_TILED` | Density matrix, blocked layout | Mixed states; cache-optimised gate operations |

All public functions dispatch at runtime on `state->type`.

---

## Directory Layout

```
include/
  qlib.h              Umbrella header; includes state.h, gate.h, channel.h
  state.h             Public API and type definitions
  q_types.h           Foundational types: cplx_t, idx_t, qubit_t, tile macros
src/
  state/
    state_internal.h  Type-specific function declarations (not public)
    state.c           Dispatch layer
    state_pure.c      PURE implementations
    state_packed.c    MIXED_PACKED implementations
    state_tiled.c     MIXED_TILED implementations
  internal/
    index.h           Shared index computation helpers
    utils.h           Debug print macros
test/
  state/
    test_state.c        Test runner + dispatch tests (2 tests)
    test_state_pure.c   PURE unit tests (10 tests)
    test_state_packed.c MIXED_PACKED unit tests (13 tests)
    test_state_tiled.c  MIXED_TILED unit tests (16 tests)
```

---

## Types

### `state_t`

```c
typedef struct state {
    state_type_t type;   // Storage format; must be set before any call
    cplx_t       *data;  // Heap-allocated, 64-byte-aligned; owned by struct
    qubit_t       qubits;
} state_t;
```

`data != NULL` iff storage is allocated and valid.

### `state_type_t`

```c
typedef enum { PURE, MIXED_PACKED, MIXED_TILED } state_type_t;
```

### `qop_t`

```c
typedef void (*qop_t)(state_t *state);
```

Function pointer type for quantum operations that modify a state in-place.
Operations with additional parameters (target qubit, rotation angle, ...)
must be wrapped in a function matching this signature.

### Types from `q_types.h`

| Type | Underlying | Notes |
|------|------------|-------|
| `cplx_t` | `double _Complex` | C99, 16 bytes |
| `idx_t` | `uint64_t` | Index/dimension type; unsigned for bitwise ops |
| `qubit_t` | `unsigned char` | Max 255 qubits |

---

## Tile Configuration

Defined in `q_types.h`, set by CMake at compile time via `-DLOG_TILE_DIM=N`:

```c
#define TILE_DIM  (1 << LOG_TILE_DIM)   // 8 (Debug) or 32 (Release)
#define TILE_SIZE (TILE_DIM * TILE_DIM)  // 64 (Debug) or 1024 (Release)
```

Do not set `LOG_TILE_DIM` manually. CMake sets it to 3 (Debug) or 5 (Release). Valid range: 3-8.

Release default (32x32): one tile = 16 KB, fits in L1 cache.

---

## Memory Layouts

### PURE — State Vector

Flat array of 2^n complex amplitudes. Index i encodes basis state |i> in little-endian binary.

```
data[i] = <i|psi>      i in {0, ..., 2^n - 1}
```

Length: `2^n`

### MIXED_PACKED — Hermitian Packed Lower Triangle

Column-major packed lower triangle, matching LAPACK `uplo='L'` format.

```
Index formula (row >= col):
  idx = col * dim - col * (col + 1) / 2 + row
```

Upper triangle: `rho[row, col] = conj(rho[col, row])` — computed on access, never stored.

Length: `dim * (dim + 1) / 2`

### MIXED_TILED — Blocked Layout

Density matrix partitioned into TILE_DIM x TILE_DIM tiles. Only lower-triangular tile blocks are stored, in row-major tile order. Every stored tile is a full TILE_DIM x TILE_DIM block — diagonal tiles store both triangles, with no secondary packing within a tile.

```
Index formula (row >= col):
  tile_row    = row / TILE_DIM
  tile_col    = col / TILE_DIM
  tile_offset = tile_row * (tile_row + 1) / 2 + tile_col
  local_row   = row % TILE_DIM
  local_col   = col % TILE_DIM
  idx         = tile_offset * TILE_SIZE + local_row * TILE_DIM + local_col
```

Intra-tile layout is **row-major**. Do not pass tile pointers to LAPACK routines.

Length: `n_tiles * (n_tiles + 1) / 2 * TILE_SIZE`

Padding: when dim < TILE_DIM, a single full tile is allocated; only the dim x dim corner is physical. Padding elements are zero.

---

## Public API

Declared in `include/state.h`. All functions dispatch on `state->type`.

| Function | Signature | Description |
|----------|-----------|-------------|
| `state_len` | `idx_t state_len(const state_t *state)` | Number of elements in storage |
| `state_init` | `void state_init(state_t *state, qubit_t qubits, cplx_t **data)` | Allocate or take ownership of storage |
| `state_free` | `void state_free(state_t *state)` | Free storage; safe to call on NULL or freed state |
| `state_plus` | `void state_plus(state_t *state, qubit_t qubits)` | Initialize to \|+>^n (pure) or \|+><+\|^n (mixed) |
| `state_cp` | `state_t state_cp(const state_t *state)` | Deep copy with independent allocation |
| `state_get` | `cplx_t state_get(const state_t *state, idx_t row, idx_t col)` | Element access with Hermitian symmetry |
| `state_set` | `void state_set(state_t *state, idx_t row, idx_t col, cplx_t val)` | Element write with Hermitian symmetry |
| `state_print` | `void state_print(const state_t *state)` | Print to stdout |

For `PURE` states, `col` must be 0 in `state_get`/`state_set` (enforced by assert).

For `MIXED_TILED`, `state_set` mirrors writes in diagonal tiles to maintain both triangles.

---

## Memory Ownership

- `state_t` owns its `data` pointer after initialization.
- `state_init(..., &data)` transfers ownership; caller's pointer is set to NULL.
- No shared ownership. Calling `state_init` on an already-initialized state frees the previous allocation.

---

## Error Handling

| Error | Handling |
|-------|----------|
| Programmer error (NULL state, invalid type, col != 0 for PURE) | `assert()` — crash in debug builds |
| Out-of-memory | `state->data = NULL`, message to stderr |

Check `state->data != NULL` after `state_init` or `state_plus`.

---

## Build Integration

State sources are compiled into the shared library target `${PROJECT_NAME_LOWER}` (e.g. `libqsim`), defined in `src/CMakeLists.txt`.

Sources: `src/state/state.c`, `state_pure.c`, `state_packed.c`, `state_tiled.c`

Test target: `test_state`, defined in `test/CMakeLists.txt`. Links the library and Unity v2.6.1.

Installed headers: `state.h`, `q_types.h`, `qlib.h` to `${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME_LOWER}/`

```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug --tests-regex test_state
```

---

## Memory Footprint (Release, TILE_DIM = 32)

| Qubits | dim | PURE | MIXED_PACKED | MIXED_TILED |
|--------|-----|------|--------------|-------------|
| 1 | 2 | 32 B | 48 B | 16 KB |
| 5 | 32 | 512 B | 8.4 KB | 16 KB |
| 8 | 256 | 4 KB | 526 KB | 576 KB |
| 10 | 1024 | 16 KB | 8.4 MB | 8.6 MB |
| 12 | 4096 | 64 KB | 134 MB | 136 MB |

MIXED_TILED overhead is significant for small systems (n < 5) but converges to MIXED_PACKED for n >= 8.
