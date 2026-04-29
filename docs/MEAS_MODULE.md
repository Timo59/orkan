# Measurement Module

Computes expectation values and matrix elements of diagonal observables on quantum states.

| Function | Pure | MIXED_PACKED | MIXED_TILED |
|---|:---:|:---:|:---:|
| `mean` | yes | yes | yes |
| `matel_diag` | yes | — | — |

A `sample_matrix` operation is planned but not currently implemented; see `docs/SAMPLE_MATRIX.md` for the implementation brief.

---

## Directory Layout

```
include/
  meas.h                  Public API

src/meas/
  meas.c                  Dispatch + input validation
  meas_pure.c             PURE backend
  meas_packed.c           MIXED_PACKED backend
  meas_tiled.c            MIXED_TILED backend

test/meas/
  test_meas.c             Test runner
  test_meas_pure.c        PURE tests
  test_meas_packed.c      MIXED_PACKED tests
  test_meas_tiled.c       MIXED_TILED tests
```

---

## Public API

Declared in `include/meas.h`.

### `mean`

```c
double mean(const state_t *state, const double *obs);
```

Computes the expectation value of a diagonal observable:

- **Pure:**  `<H> = sum_x obs[x] |psi(x)|^2`
- **Mixed:** `<H> = sum_x obs[x] rho_{xx}`

| Parameter | Type | Description |
|---|---|---|
| `state` | `const state_t *` | Quantum state in any representation |
| `obs` | `const double *` | Diagonal coefficients, length `2^(state->qubits)` |
| return | `double` | Real expectation value `Tr(rho H)` |

### `matel_diag`

```c
cplx_t matel_diag(const state_t *bra, const state_t *ket, const double *obs);
```

Matrix element of a diagonal observable between two pure states:

`<bra| H |ket> = sum_x obs[x] conj(bra[x]) ket[x]`

When `bra == ket`, reduces to `mean()` (imaginary part zero to machine precision). Both states must be `PURE` with matching qubit counts.

| Parameter | Type | Description |
|---|---|---|
| `bra` | `const state_t *` | Pure quantum state (bra vector) |
| `ket` | `const state_t *` | Pure quantum state (ket vector) |
| `obs` | `const double *` | Diagonal coefficients, length `2^(bra->qubits)` |
| return | `cplx_t` | Complex matrix element `<bra\| H \|ket>` |

Input validation uses `MEAS_VALIDATE`: non-NULL `state`/`state->data`/`obs`; valid `state->type`. Failures print to stderr and call `exit()`.

---

## Architecture

### Pure — `mean_pure`

Iterates all `2^n` amplitudes and accumulates `obs[x] * |psi[x]|^2` where `|psi|^2 = creal^2 + cimag^2` (avoiding the `cabs` sqrt). OpenMP: `#pragma omp parallel for reduction(+:sum)` over the full state vector.

### Pure — `matel_diag_pure`

Same loop structure as `mean_pure` but accumulates two `double` reductions (`sum_re`, `sum_im`) for the real and imaginary parts of `obs[x] * conj(bra[x]) * ket[x]`. The loop body expands complex arithmetic into real operations (`ar*br + ai*bi`, `ar*bi - ai*br`) for auto-vectorization.

### Packed — `mean_packed`

Diagonal elements of packed lower-triangular column-major storage sit at indices `0, dim, 2*dim-1, 3*dim-3, ...` with strides `dim, dim-1, dim-2, ...`.

Each OpenMP thread computes its starting packed index once via `pack_idx(dim, x_start, x_start)` from `index.h`, then uses the stride-decrement pattern sequentially within its chunk.

Only `creal` is needed: diagonal elements of a Hermitian matrix are real.

### Tiled — `mean_tiled`

Diagonal elements live on diagonal tiles `(k, k)`. Within each tile the diagonal is strided by `TILE_DIM + 1` (row-major intra-tile layout).

Diagonal tile positions (in `TILE_SIZE` units): `0, 2, 5, 9, 14, ...`. Stride between consecutive diagonal tiles: starts at 2, increments by 1.

Each OpenMP thread processes a contiguous chunk of diagonal tiles. Within each tile, diagonal elements are gathered to a contiguous stack array `double diag[TILE_DIM]` to enable auto-vectorisation of the dot product.

`len_tile` is computed once: either `dim` (single tile, `qubits < LOG_TILE_DIM`) or `TILE_DIM` (all diagonal tiles full, since both `dim` and `TILE_DIM` are powers of 2).

---

## Build Integration

Sources compiled into `${PROJECT_NAME_LOWER}` alongside the state/gate/channel/circ modules:

```cmake
set(MEAS_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/meas/meas.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/meas/meas_pure.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/meas/meas_packed.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/meas/meas_tiled.c"
)
```

The current `meas` module is BLAS-free. A future `sample_matrix` implementation will introduce an ILP64 BLAS dependency confined to its own translation unit.

Public header `include/meas.h` is installed alongside other module headers. Umbrella header `include/orkan.h` includes `meas.h`.

Test target `test_meas` is registered in `test/CMakeLists.txt` and added to `verified_install` DEPENDS in the root `CMakeLists.txt`.

```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug -R test_meas
```

---

## Tests

Test framework, build-type requirement, and shared fixtures: see `docs/TEST_SUITE.md`. Shared helpers `fill_identity_obs`, `fill_uniform_obs` in `test/include/test_meas.h`. Tests run at qubit counts 1 through `MAXQUBITS` (4). Tolerance 1e-12.

### PURE — `mean` (`test/meas/test_meas_pure.c`)

| Test | State | Observable | Expected |
|---|---|---|---|
| `test_mean_pure_uniform` | \|+>^n | `obs[x] = 1` | 1.0 |
| `test_mean_pure_identity` | \|+>^n | `obs[x] = x` | (2^n − 1)/2 |
| `test_mean_pure_basis` | \|0> | `obs[x] = x` | 0.0 |
| `test_mean_pure_basis_k` | \|k> | `obs[x] = x` | k |
| `test_mean_pure_imaginary` | (i/√2, i/√2) | {0, 1} | 0.5 |
| `test_mean_pure_negative_obs` | \|+>^n | `obs[x] = −1` | −1.0 |
| `test_mean_pure_mixed_sign_obs` | \|+>^2 | {−3, +1, −1, +3} | 0.0 |

### PURE — `matel_diag` (`test/meas/test_meas_pure.c`)

| Test | Verifies |
|---|---|
| `test_matel_diag_pure_self` | bra=ket reduces to `mean(state, obs)` |
| `test_matel_diag_pure_orthogonal` | <0\|H\|1> = 0 for any diagonal H |
| `test_matel_diag_pure_conjugate_symmetry` | <a\|H\|b> = conj(<b\|H\|a>) on real states |
| `test_matel_diag_pure_conjugate_symmetry_complex` | Same on complex states; ab ≠ ba in general |
| `test_matel_diag_pure_known` | Hand-computed value for (1,i)/√2 vs \|+>, obs={0,1} → −i/2 |
| `test_matel_diag_pure_complex_2q` | 2-qubit hand-computed value: \|+>^2 vs i\|+>^2, obs={0,1,2,3} → 1.5i |

### MIXED_PACKED — `mean` (`test/meas/test_meas_packed.c`)

| Test | State | Observable | Expected |
|---|---|---|---|
| `test_mean_packed_uniform` | \|+><+\|^n | `obs[x] = 1` | 1.0 |
| `test_mean_packed_identity` | \|+><+\|^n | `obs[x] = x` | (2^n − 1)/2 |
| `test_mean_packed_pure_embed` | \|0><0\| | `obs[x] = x` | 0.0 |
| `test_mean_packed_basis_k` | \|k><k\| | `obs[x] = x` | k |
| `test_mean_packed_negative_obs` | \|+><+\|^n | `obs[x] = −1` | −1.0 |
| `test_mean_packed_mixed_sign_obs` | \|+><+\|^2 | {−3, +1, −1, +3} | 0.0 |

### MIXED_TILED — `mean` (`test/meas/test_meas_tiled.c`)

| Test | State | Observable | Expected |
|---|---|---|---|
| `test_mean_tiled_uniform` | \|+><+\|^n | `obs[x] = 1` | 1.0 |
| `test_mean_tiled_identity` | \|+><+\|^n | `obs[x] = x` | (2^n − 1)/2 |
| `test_mean_tiled_pure_embed` | \|0><0\| | `obs[x] = x` | 0.0 |
| `test_mean_tiled_small` | n=1 (single tile) | `obs[x] = 1` | 1.0 |
| `test_mean_tiled_multi` | n = LOG_TILE_DIM+1 | `obs[x] = x` | (2^n − 1)/2 |
| `test_mean_tiled_basis_k` | \|k><k\| | `obs[x] = x` | k |
| `test_mean_tiled_negative_obs` | \|+><+\|^n | `obs[x] = −1` | −1.0 |
| `test_mean_tiled_mixed_sign_obs` | \|+><+\|^2 | {−3, +1, −1, +3} | 0.0 |

### Cross-backend (`test/meas/test_meas.c`)

| Test | Verifies |
|---|---|
| `test_mean_consistency` | `mean()` returns identical values on PURE / MIXED_PACKED / MIXED_TILED for the same logical state and observable |
