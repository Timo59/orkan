# MEAS Module: Diagonal Observable Expectation Value

**Module:** Measurement
**Status:** In progress
**Last Updated:** 2026-04-24

---

## Scope

This module provides measurement operations on quantum states. Currently
implemented:

| Function        | Pure | Mixed (packed) | Mixed (tiled) |
|-----------------|:----:|:--------------:|:-------------:|
| `mean`          | yes  | yes            | yes           |
| `matel_diag`    | yes  | —              | —             |
| `sample_matrix` | yes  | —              | —             |

---

## Directory Layout

```
QSim/
├── include/
│   └── meas.h                  # Public API
├── src/
│   └── meas/
│       ├── meas.c              # Dispatch layer
│       ├── meas_pure.c         # Pure state backend
│       ├── meas_packed.c       # Packed mixed state backend
│       └── meas_tiled.c        # Tiled mixed state backend
├── test/
│   ├── include/
│   │   └── test_meas.h         # Shared test helpers
│   └── meas/
│       ├── test_meas.c         # Test runner (22 tests)
│       ├── test_meas_pure.c    # Pure state tests
│       ├── test_meas_packed.c  # Packed state tests
│       └── test_meas_tiled.c   # Tiled state tests
└── docs/
    └── MEAS_MODULE.md          # This file
```

---

## API

### Types

```c
typedef void (*qop_t)(state_t *state);
```

Function pointer type for quantum operations that modify a state in-place.
Operations with additional parameters (target qubit, rotation angle, ...) must
be wrapped in a function matching this signature.

### `mean`

```c
double mean(const state_t *state, const double *obs);
```

Computes the expectation value of a diagonal observable:

- **Pure:**  `<H> = sum_x obs[x] |psi(x)|^2`
- **Mixed:** `<H> = sum_x obs[x] rho_{xx}`

**Parameters:**

| Parameter | Type              | Description                                      |
|-----------|-------------------|--------------------------------------------------|
| `state`   | `const state_t *` | Quantum state in any representation               |
| `obs`     | `const double *`  | Diagonal coefficients, length `2^(state->qubits)` |
| return    | `double`          | Real expectation value `Tr(rho H)`                |

**Preconditions** (checked via `MEAS_VALIDATE`; calls `fprintf(stderr) + exit()`
on failure, surviving Release builds): `state`, `state->data`, and `obs` must
be non-NULL; `state->type` must be a valid `state_type_t` variant.

### `matel_diag`

```c
cplx_t matel_diag(const state_t *bra, const state_t *ket, const double *obs);
```

Matrix element of a diagonal observable between two pure states:

`<bra| H |ket> = sum_x obs[x] conj(bra[x]) ket[x]`

When `bra == ket` this reduces to `mean()` (imaginary part zero to machine
precision).

**Parameters:**

| Parameter | Type              | Description                                      |
|-----------|-------------------|--------------------------------------------------|
| `bra`     | `const state_t *` | Pure quantum state (bra vector)                   |
| `ket`     | `const state_t *` | Pure quantum state (ket vector)                   |
| `obs`     | `const double *`  | Diagonal coefficients, length `2^(bra->qubits)`   |
| return    | `cplx_t`          | Complex matrix element `<bra\| H \|ket>`          |

**Preconditions:** Both states must be `PURE` with matching qubit counts.
Mixed states cause program termination with a diagnostic message.

### `sample_matrix`

```c
void sample_matrix(const state_t *state, const double *obs,
                   qop_t *const *ops, const idx_t *sizes, idx_t n_layers,
                   cplx_t *out);
```

Computes the Hermitian sample matrix `S` with entries

`S_{i,j} = <iota| V_i^+ H V_j |iota>`

where `V_k = U_{k_L}^{(L)} ... U_{k_1}^{(1)}` is the composite operation for
multi-index `k`, and `H = diag(obs)`.

The flat index `k` maps to multi-index `(k_1, ..., k_L)` via mixed-radix
decomposition with `k_1` least significant (fastest-varying):

`k = k_L * (n_{L-1} * ... * n_1) + ... + k_2 * n_1 + k_1`

Output is the packed lower triangle of `S` in column-major order, `K(K+1)/2`
elements where `K = sizes[0] * ... * sizes[n_layers-1]`.

**Parameters:**

| Parameter   | Type                | Description                                      |
|-------------|---------------------|--------------------------------------------------|
| `state`     | `const state_t *`   | Pure quantum state \|iota>                       |
| `obs`       | `const double *`    | Diagonal coefficients, length `2^(state->qubits)` |
| `ops`       | `qop_t *const *`    | Layered operations: `ops[l][k]` is the k-th op in layer l |
| `sizes`     | `const idx_t *`     | Number of operations per layer, length `n_layers` |
| `n_layers`  | `idx_t`             | Number of operation layers (L)                    |
| `out`       | `cplx_t *`          | Packed lower triangle, caller-allocated, length `K*(K+1)/2` |

**Preconditions:** State must be `PURE`. Mixed states cause program termination.
`n_layers` must be positive. Memory: holds at most two state-vector copies
simultaneously.

---

## Backend Implementations

### Pure (`meas_pure.c`)

**`mean_pure`:** Iterates all `2^n` amplitudes, accumulates
`obs[x] * |psi[x]|^2` where `|psi|^2 = creal^2 + cimag^2` (avoids `cabs`
sqrt).  OpenMP: `#pragma omp parallel for reduction(+:sum)` over the full
state vector.

**`matel_diag_pure`:** Same loop structure as `mean_pure` but accumulates two
`double` reductions (`sum_re`, `sum_im`) for the real and imaginary parts of
`obs[x] * conj(bra[x]) * ket[x]`.  The loop body expands complex arithmetic
into real operations (`ar*br + ai*bi`, `ar*bi - ai*br`) for auto-vectorization.

**`apply_composite`:** Static inline helper.  Decomposes flat index `k` into
multi-index via repeated `% sizes[l]` / `/ sizes[l]` (k_1 least significant)
and applies `ops[l][k_l]` sequentially.

**`sample_matrix_pure`:** Iterates the packed lower triangle column by column.
For each column `j`, builds `phi_j = V_j |iota>` (the "ket") via `state_cp` +
`apply_composite`.  The diagonal entry uses `phi_j` as both bra and ket.
Sub-diagonal entries build `phi_i` from the input state each time.  No nested
OpenMP — gate/channel operations already parallelise internally.

### Packed (`meas_packed.c`)

Diagonal elements of packed lower-triangular column-major storage sit at
indices `0, dim, 2*dim-1, 3*dim-3, ...` with strides `dim, dim-1, dim-2, ...`

Each OpenMP thread computes its starting packed index once via
`pack_idx(dim, x_start, x_start)` from `index.h`, then uses the
stride-decrement pattern sequentially within its chunk.

Only `creal` is needed: diagonal elements of a Hermitian matrix are real.

### Tiled (`meas_tiled.c`)

Diagonal elements live on diagonal tiles `(k, k)`.  Within each tile the
diagonal is strided by `TILE_DIM + 1` (row-major intra-tile layout).

Diagonal tile positions (in `TILE_SIZE` units): `0, 2, 5, 9, 14, ...`
Stride between consecutive diagonal tiles: `2, 3, 4, 5, ...` (starts at 2,
increments by 1 each iteration).

Each OpenMP thread processes a contiguous chunk of diagonal tiles.  Within
each tile, diagonal elements are gathered to a contiguous stack array
`double diag[TILE_DIM]` to enable auto-vectorisation of the dot product.

`len_tile` is computed once: either `dim` (single tile, `qubits < LOG_TILE_DIM`)
or `TILE_DIM` (all diagonal tiles are full since both `dim` and `TILE_DIM` are
powers of 2).

---

## Build Integration

The measurement module sources are compiled into the `orkan` shared library
alongside the state, gate, and channel modules.

**Sources** added to `Q_SOURCES` in `src/CMakeLists.txt`:

```cmake
set(MEAS_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/meas/meas.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/meas/meas_pure.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/meas/meas_packed.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/meas/meas_tiled.c"
)
```

**Public header** `include/meas.h` is installed alongside other headers.

**Umbrella header** `include/qlib.h` includes `meas.h`.

**Test target** `test_meas` is registered in `test/CMakeLists.txt` and added
to `verified_install` DEPENDS in the root `CMakeLists.txt`.

---

## Test Suite

33 tests across seven categories, run at qubit counts 1 through `MAXQUBITS` (4).
Shared helpers `fill_identity_obs` and `fill_uniform_obs` live in
`test/include/test_meas.h`.

| Test | State | Observable | Expected |
|------|-------|-----------|----------|
| `test_mean_pure_uniform` | \|+⟩^n | obs[x]=1 | 1.0 |
| `test_mean_pure_identity` | \|+⟩^n | obs[x]=x | (2^n−1)/2 |
| `test_mean_pure_basis` | \|0⟩ | obs[x]=x | 0.0 |
| `test_mean_pure_basis_k` | \|k⟩ | obs[x]=x | k |
| `test_mean_pure_imaginary` | (i/√2, i/√2) | obs={0,1} | 0.5 |
| `test_mean_pure_negative_obs` | \|+⟩^n | obs[x]=−1 | −1.0 |
| `test_mean_pure_mixed_sign_obs` | \|+⟩^2 | {−3,+1,−1,+3} | 0.0 |
| `test_matel_diag_pure_self` | \|+⟩^n, bra=ket | obs[x]=x | mean(s, obs) |
| `test_matel_diag_pure_orthogonal` | \|0⟩, \|1⟩ | obs[x]=x | 0 |
| `test_matel_diag_pure_conjugate_symmetry` | \|+⟩^n, \|0⟩ | obs[x]=x | conj(swap) |
| `test_matel_diag_pure_known` | (1,i)/√2, \|+⟩ | {0,1} | −i/2 |
| `test_sample_matrix_pure_trivial` | \|+⟩^n, 1 id op | obs[x]=x | [mean] |
| `test_sample_matrix_pure_two_ops` | \|0⟩, {I,X} | {0,1} | {0,0,1} |
| `test_sample_matrix_pure_known_values` | \|+⟩, {I,H} | {0,1} | {0.5, 0, 0} |
| `test_sample_matrix_pure_multi_layer` | \|0⟩, 2 layers | {0,1} | known 4×4 |
| `test_sample_matrix_pure_asymmetric_layers` | \|0⟩, sizes={3,2} | {0,1} | diag discriminates ordering |
| `test_sample_matrix_pure_state_unchanged` | \|+⟩^2, {I,X} | obs[x]=x | input state preserved |
| `test_matel_diag_pure_complex_2q` | \|+⟩^2, i\|+⟩^2 | {0,1,2,3} | 1.5i |
| `test_mean_packed_uniform` | \|+⟩⟨+\|^⊗n | obs[x]=1 | 1.0 |
| `test_mean_packed_identity` | \|+⟩⟨+\|^⊗n | obs[x]=x | (2^n−1)/2 |
| `test_mean_packed_pure_embed` | \|0⟩⟨0\| | obs[x]=x | 0.0 |
| `test_mean_packed_basis_k` | \|k⟩⟨k\| | obs[x]=x | k |
| `test_mean_packed_negative_obs` | \|+⟩⟨+\|^⊗n | obs[x]=−1 | −1.0 |
| `test_mean_packed_mixed_sign_obs` | \|+⟩⟨+\|^2 | {−3,+1,−1,+3} | 0.0 |
| `test_mean_tiled_uniform` | \|+⟩⟨+\|^⊗n | obs[x]=1 | 1.0 |
| `test_mean_tiled_identity` | \|+⟩⟨+\|^⊗n | obs[x]=x | (2^n−1)/2 |
| `test_mean_tiled_pure_embed` | \|0⟩⟨0\| | obs[x]=x | 0.0 |
| `test_mean_tiled_small` | n=1 (< LOG_TILE_DIM) | obs[x]=1 | 1.0 |
| `test_mean_tiled_multi` | n=LOG_TILE_DIM+1 | obs[x]=x | (2^n−1)/2 |
| `test_mean_tiled_basis_k` | \|k⟩⟨k\| | obs[x]=x | k |
| `test_mean_tiled_negative_obs` | \|+⟩⟨+\|^⊗n | obs[x]=−1 | −1.0 |
| `test_mean_tiled_mixed_sign_obs` | \|+⟩⟨+\|^2 | {−3,+1,−1,+3} | 0.0 |
| `test_mean_consistency` | same state, 3 types | same obs | equal |
