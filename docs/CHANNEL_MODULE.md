# Channel Module

Simulates local quantum channels on mixed states via the superoperator formalism.

---

## Overview

A quantum channel E acts on a density matrix as E(rho) = sum_k K_k rho K_k-dagger, where {K_k} are the Kraus operators. The module converts the Kraus representation to a 4x4 superoperator matrix S, then applies S element-wise to 2x2 blocks of the density matrix across packed and tiled backends. Only `MIXED_PACKED` and `MIXED_TILED` states are supported; `PURE` states are rejected.

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
  test_channel.c               Test runner (35 tests)
  test_channel_superop.c       kraus_to_superop unit tests (7 tests)
  test_channel_packed.c        Packed sweep + analytical + trace tests (12 tests)
  test_channel_tiled.c         Tiled sweep + analytical + trace + cross-backend tests (16 tests)
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

## Dispatch Architecture

`channel_1q` validates inputs (null pointers, qubit range, superoperator locality), then dispatches:

```c
switch (state->type) {
    case PURE:         error — not supported
    case MIXED_PACKED: channel_packed_1q(state, sop->data, target);  break;
    case MIXED_TILED:  channel_tiled_1q(state, sop->data, target);   break;
    default:           error — unknown state type
}
```

---

## Algorithm

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

---

## OpenMP Parallelization

`OMP_THRESHOLD` is defined locally in `channel_packed.c` and `channel_tiled.c` with a default of 512. In Debug builds, CMake overrides it to 4 so multi-threaded paths are exercised at small qubit counts.

---

## Build Integration

Channel sources are compiled into the shared library target `${PROJECT_NAME_LOWER}`, defined in `src/CMakeLists.txt`.

Sources: `src/channel/channel.c`, `channel_packed.c`, `channel_tiled.c`

Test target: `test_channel`, defined in `test/CMakeLists.txt`. Links the library, BLAS, and Unity v2.6.1.

Installed header: `channel.h` to `${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME_LOWER}/`

```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug -R test_channel
```

---

## Test Infrastructure

### Test Harnesses

| Harness | File | Backend |
|---------|------|---------|
| `testChannel1qPacked` | `test_channel_packed.c` | MIXED_PACKED |
| `testChannel1qTiled` | `test_channel_tiled.c` | MIXED_TILED |

Both sweep qubit counts, all target positions, and multiple input states. Reference is computed via BLAS-based sum_k K_k rho K_k-dagger using full n-qubit Kraus matrices.

### Test Groups

| Group | Count | Description |
|-------|-------|-------------|
| kraus_to_superop | 7 | Analytical superoperator verification |
| Packed sweep | 5 | Full sweep for identity, bit-flip, phase-flip, depolarising, amplitude damping |
| Tiled sweep | 5 | Same channels on tiled states |
| Packed analytical | 4 | Boundary cases (full decay, full depolarisation, involution) |
| Tiled analytical | 4 | Same on tiled states |
| Trace preservation | 6 | Tr(E(rho)) = Tr(rho) for 3 channels x 2 backends |
| Cross-backend | 4 | Element-wise packed vs tiled agreement |
| **Total** | **35** | |
