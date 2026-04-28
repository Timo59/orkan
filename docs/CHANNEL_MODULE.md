# Channel Module

Simulates local quantum channels on mixed states via the superoperator formalism.

A channel E acts on a density matrix as E(rho) = sum_k K_k rho K_k-dagger, where {K_k} are the Kraus operators. The module converts the Kraus representation to a 4x4 superoperator matrix S, then applies S element-wise to 2x2 blocks of the density matrix across packed and tiled backends. Only `MIXED_PACKED` and `MIXED_TILED` states are supported; `PURE` states are rejected.

---

## Directory Layout

```
include/
  channel.h                    Public API: kraus_t, superop_t, kraus_to_superop(), channel_1q()

src/channel/
  channel.c                    Superoperator construction + channel_1q dispatcher
  channel_packed.c             MIXED_PACKED backend
  channel_tiled.c              MIXED_TILED backend

src/internal/
  index.h                      Shared index helpers (insertBit0, col_off, tile_off, elem_off)

test/channel/
  test_channel.c               Test runner
  test_channel_superop.c       kraus_to_superop unit tests
  test_channel_packed.c        Packed sweep + analytical + trace tests
  test_channel_tiled.c         Tiled sweep + analytical + trace + cross-backend tests
```

---

## Public API

Declared in `include/channel.h`.

### Types

| Type | Description |
|------|-------------|
| `kraus_t` | Kraus representation: `n_qubits`, `n_terms`, `data` (concatenated Kraus matrices, row-major) |
| `superop_t` | Superoperator matrix: `n_qubits`, `data` (N^2 x N^2 matrix, row-major) |

### Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `kraus_to_superop` | `superop_t kraus_to_superop(const kraus_t *kraus)` | Convert Kraus to superoperator. Caller must `free(sop.data)`. |
| `channel_1q` | `void channel_1q(state_t *state, const superop_t *sop, qubit_t target)` | Apply 1-qubit channel in-place. Dispatches to packed or tiled backend. |

Input validation uses `CHANNEL_VALIDATE` (defined in `channel.c`), which prints to stderr and calls `exit()` on failure.

---

## Architecture

### Dispatch

`channel_1q` validates inputs (null pointers, qubit range, superoperator locality), then dispatches:

```c
switch (state->type) {
    case PURE:         error — not supported
    case MIXED_PACKED: channel_packed_1q(state, sop->data, target);  break;
    case MIXED_TILED:  channel_tiled_1q(state, sop->data, target);   break;
    default:           error — unknown state type
}
```

### Kraus-to-Superoperator Conversion

For N = 2^n and superoperator dimension M = N^2:

```
S_{(i,j),(k,l)} = sum_a K_a[i,k] * conj(K_a[j,l])
```

For 1-qubit channels (N=2), this produces a 4x4 matrix.

### SOP_PRELOAD / SOP_APPLY Macros

Both backends define identical file-local macros:

- **SOP_PRELOAD**: Decomposes 16 complex superoperator entries into 32 `const double` locals for register promotion.
- **SOP_APPLY**: Computes the 4x4 product in real arithmetic (64 multiplications, 56 additions per block), enabling SIMD auto-vectorisation.

Conjugation for Hermitian symmetry is handled at load/store sites by sign-flipping imaginary parts, keeping SOP_APPLY itself pure.

### Packed Backend

Outer loop over column-base (OpenMP-parallelised), inner loop over row-base. Three symmetry branches handle lower triangle, upper-lower mixed, and diagonal blocks.

### Tiled Backend

Two cases based on target qubit position:

- **Case 1** (target < LOG_TILE_DIM): All four block elements within the same tile.
- **Case 2** (target >= LOG_TILE_DIM): Elements span four tiles. Three sub-cases for tile triangle structure: all-lower (vectorisable flat loop), mixed (conjugate-transpose access), and diagonal.

### OpenMP

`OMP_THRESHOLD` is defined locally in `channel_packed.c` and `channel_tiled.c` with a default of 512. In Debug builds, CMake overrides it to 4.

---

## Build Integration

Channel sources are compiled into `${PROJECT_NAME_LOWER}` via `src/CMakeLists.txt`.

Sources: `src/channel/channel.c`, `channel_packed.c`, `channel_tiled.c`

Test target `test_channel` in `test/CMakeLists.txt` links the library, BLAS, and Unity v2.6.1.

Installed header: `channel.h` → `${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME_LOWER}/`.

```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug -R test_channel
```

---

## Tests

| File | Coverage |
|---|---|
| `test/channel/test_channel.c` | Runner |
| `test/channel/test_channel_superop.c` | Analytical superoperator verification |
| `test/channel/test_channel_packed.c` | Packed sweep, analytical, trace preservation |
| `test/channel/test_channel_tiled.c` | Tiled sweep, analytical, trace preservation, cross-backend agreement |

Sweeps cover qubit counts, all target positions, and multiple input states (identity, bit-flip, phase-flip, depolarising, amplitude damping). Reference is computed via BLAS `sum_k K_k rho K_k-dagger` using full n-qubit Kraus matrices, tolerance 1e-12.
