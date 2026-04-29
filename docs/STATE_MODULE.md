# State Module

Represents n-qubit quantum systems in three storage formats. All public functions dispatch at runtime on `state->type`.

| Type | Object | Use case |
|------|--------|----------|
| `PURE` | State vector in C^(2^n) | Pure quantum states |
| `MIXED_PACKED` | Density matrix, Hermitian packed | Mixed states; LAPACK-compatible layout |
| `MIXED_TILED` | Density matrix, blocked layout | Mixed states; cache-optimised gate operations |

---

## Directory Layout

```
include/
  orkan.h             Umbrella header; includes state.h, gate.h, channel.h, meas.h, circ.h
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
test/state/
  test_state.c        Test runner + dispatch tests
  test_state_pure.c   PURE unit tests
  test_state_packed.c MIXED_PACKED unit tests
  test_state_tiled.c  MIXED_TILED unit tests
```

---

## Public API

Declared in `include/state.h`. All functions dispatch on `state->type`.

### Types

```c
typedef struct state {
    state_type_t type;   // Storage format; must be set before any call
    cplx_t       *data;  // Heap-allocated, 64-byte-aligned; owned by struct
    qubit_t       qubits;
} state_t;

typedef enum { PURE, MIXED_PACKED, MIXED_TILED } state_type_t;

typedef void (*qop_t)(state_t *state);
```

`data != NULL` iff storage is allocated and valid.

`qop_t` is the function-pointer type for quantum operations that modify a state in-place. Operations with extra parameters (target qubit, rotation angle, ...) must be wrapped to match this signature.

### Types from `q_types.h`

| Type | Underlying | Notes |
|------|------------|-------|
| `cplx_t` | `double _Complex` | C99, 16 bytes |
| `idx_t` | `uint64_t` | Index/dimension type; unsigned for bitwise ops |
| `qubit_t` | `unsigned char` | Max 255 qubits |

### Functions

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

For `PURE` states, `col` must be 0 in `state_get`/`state_set` (enforced by assert). For `MIXED_TILED`, `state_set` mirrors writes in diagonal tiles to maintain both triangles.

### Memory Ownership

- `state_t` owns its `data` pointer after initialization.
- `state_init(..., &data)` transfers ownership; caller's pointer is set to NULL.
- No shared ownership. Calling `state_init` on an already-initialized state frees the previous allocation.

### Error Handling

| Error | Handling |
|-------|----------|
| Programmer error (NULL state, invalid type, col != 0 for PURE) | `assert()` — crash in debug builds |
| Out-of-memory | `state->data = NULL`, message to stderr |

Check `state->data != NULL` after `state_init` or `state_plus`.

---

## Architecture

### Tile Configuration

Defined in `q_types.h`, set by CMake at compile time via `-DLOG_TILE_DIM=N`:

```c
#define TILE_DIM  (1 << LOG_TILE_DIM)   // 8 (Debug) or 32 (Release)
#define TILE_SIZE (TILE_DIM * TILE_DIM)  // 64 (Debug) or 1024 (Release)
```

CMake sets `LOG_TILE_DIM` to 3 (Debug) or 5 (Release); valid range 3-8. Do not override manually. Release default (32x32): one tile = 16 KB, fits in L1 cache.

### Memory Layouts

#### PURE — State Vector

Flat array of 2^n complex amplitudes. Index i encodes basis state |i> in little-endian binary.

```
data[i] = <i|psi>      i in {0, ..., 2^n - 1}
```

Length: `2^n`

#### MIXED_PACKED — Hermitian Packed Lower Triangle

Column-major packed lower triangle, matching LAPACK `uplo='L'` format.

```
Index formula (row >= col):
  idx = col * dim - col * (col + 1) / 2 + row
```

Upper triangle: `rho[row, col] = conj(rho[col, row])` — computed on access, never stored.

Length: `dim * (dim + 1) / 2`

#### MIXED_TILED — Blocked Layout

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

### Memory Footprint (Release, TILE_DIM = 32)

| Qubits | dim | PURE | MIXED_PACKED | MIXED_TILED |
|--------|-----|------|--------------|-------------|
| 1 | 2 | 32 B | 48 B | 16 KB |
| 5 | 32 | 512 B | 8.4 KB | 16 KB |
| 8 | 256 | 4 KB | 526 KB | 576 KB |
| 10 | 1024 | 16 KB | 8.4 MB | 8.6 MB |
| 12 | 4096 | 64 KB | 134 MB | 136 MB |

MIXED_TILED overhead is significant for small systems (n < 5) but converges to MIXED_PACKED for n >= 8.

---

## Build Integration

State sources are compiled into `${PROJECT_NAME_LOWER}` via `src/CMakeLists.txt`.

Sources: `src/state/state.c`, `state_pure.c`, `state_packed.c`, `state_tiled.c`

Test target `test_state` in `test/CMakeLists.txt` links the library and Unity v2.6.1.

Installed headers: `state.h`, `q_types.h`, `orkan.h` → `${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME_LOWER}/`.

```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug -R test_state
```

---

## Tests

`test_state` does not link BLAS — `linalg.c` is excluded from this target. Test framework, build-type requirement, and shared fixtures: see `docs/TEST_SUITE.md`.

### PURE (`test/state/test_state_pure.c`)

| Test | Verifies |
|---|---|
| `test_pure_len` | `state_len()` returns `2^qubits` |
| `test_pure_init_null_data` | `state_init` allocates fresh storage when `*data == NULL` |
| `test_pure_init_ownership` | `state_init` transfers ownership; caller's pointer is set to NULL |
| `test_pure_free_double_free` | `state_free` is safe on a NULL or already-freed state |
| `test_pure_plus_values` | `state_plus` initialises to \|+>^n (uniform amplitude 1/sqrt(2^n)) |
| `test_pure_get_set_roundtrip` | `state_set` followed by `state_get` returns the value written |
| `test_pure_cp_independence` | `state_cp` produces an independent allocation |
| `test_pure_single_qubit` | Single-qubit state behaviour |
| `test_pure_four_qubits` | Four-qubit state behaviour |
| `test_pure_reinit` | Re-initialising an already-initialised state frees the previous allocation |

### MIXED_PACKED (`test/state/test_state_packed.c`)

| Test | Verifies |
|---|---|
| `test_packed_len` | `state_len()` returns `dim*(dim+1)/2` |
| `test_packed_init_null_data` | `state_init` allocates packed storage |
| `test_packed_init_ownership` | Ownership transfer to packed state |
| `test_packed_free_double_free` | Safe double-free |
| `test_packed_plus_values` | `state_plus` initialises to \|+><+\|^n |
| `test_packed_get_lower` | Direct read on lower triangle |
| `test_packed_get_upper` | Upper-triangle read returns `conj` of lower |
| `test_packed_set_lower` | Direct write on lower triangle |
| `test_packed_set_upper` | Upper-triangle write stores `conj` into lower |
| `test_packed_cp_independence` | Deep copy independence |
| `test_packed_hermitian_symmetry` | `rho[r,c] = conj(rho[c,r])` invariant maintained |
| `test_packed_single_qubit` | Single-qubit packed state |
| `test_packed_zero_qubit` | 0-qubit edge case |

### MIXED_TILED (`test/state/test_state_tiled.c`)

| Test | Verifies |
|---|---|
| `test_tiled_len_small` | Storage length when `dim < TILE_DIM` (single-tile case) |
| `test_tiled_len_multi_tile` | Storage length when `dim >= TILE_DIM` |
| `test_tiled_init_null_data` | Tiled allocation |
| `test_tiled_init_ownership` | Ownership transfer to tiled state |
| `test_tiled_free_double_free` | Safe double-free |
| `test_tiled_plus_values` | `state_plus` on tiled storage |
| `test_tiled_get_within_tile` | Read when row and col live in the same tile |
| `test_tiled_get_cross_tile` | Read when row and col live in different tiles |
| `test_tiled_get_upper` | Upper-triangle read returns `conj` of mirror |
| `test_tiled_set_roundtrip` | Set/get round-trip across tile boundaries |
| `test_tiled_cp_independence` | Deep copy of tiled storage |
| `test_tiled_hermitian_symmetry` | Symmetry maintained on tiled layout |
| `test_tiled_single_qubit` | Single-qubit (single-tile) state |
| `test_tiled_tile_boundaries` | Index correctness at tile-boundary rows/cols |
| `test_tiled_set_upper` | Upper-triangle write mirrors into lower (and into both halves of diagonal tiles) |
| `test_tiled_multi_tile` | Multi-tile state behaviour |

### Dispatch (`test/state/test_state.c`)

| Test | Verifies |
|---|---|
| `test_dispatch_roundtrip` | Dispatcher selects the correct backend for each `state_type_t` |
| `test_dispatch_print` | `state_print` works for all three storage types |
