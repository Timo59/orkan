# Circuit Module: Diagonal Unitary Evolution

**Module:** Circuit
**Status:** In progress
**Last Updated:** 2026-04-15

---

## Scope

This module provides parameterised circuit operations on quantum states.
Currently implemented:

| Function   | Pure | Mixed (packed) | Mixed (tiled) |
|------------|:----:|:--------------:|:-------------:|
| `exp_diag` | yes  | yes            | yes           |

---

## Directory Layout

```
QSim/
├── include/
│   └── circ.h                  # Public API
├── src/
│   └── circ/
│       ├── circ.c              # Dispatch layer
│       ├── circ_pure.c         # Pure state backend
│       ├── circ_packed.c       # Packed mixed state backend
│       └── circ_tiled.c        # Tiled mixed state backend
├── test/
│   ├── include/
│   │   └── test_circ.h         # Shared test helpers
│   └── circ/
│       ├── test_circ.c         # Test runner (18 tests)
│       ├── test_circ_pure.c    # Pure state tests
│       ├── test_circ_packed.c  # Packed state tests
│       └── test_circ_tiled.c   # Tiled state tests
└── docs/
    └── CIRCUIT_MODULE.md       # This file
```

---

## API

```c
void exp_diag(state_t *state, const double *diag, double t);
```

Applies the diagonal unitary U = exp(-i D t) in-place:

- **Pure:**  `psi[x] *= exp(-i diag[x] t)`
- **Mixed:** `rho[r][c] *= exp(-i (diag[r] - diag[c]) t)`

Diagonal density-matrix elements are unchanged (phase difference is zero).

No implicit factor of 1/2.  To match the QAOA convention e^{-i gamma C},
pass gamma directly.  For the mixer layer with RX gates, pass 2*beta.

**Parameters:**

| Parameter | Type              | Description                                      |
|-----------|-------------------|--------------------------------------------------|
| `state`   | `state_t *`       | Quantum state (any representation), modified in-place |
| `diag`    | `const double *`  | Diagonal coefficients, length `2^(state->qubits)` |
| `t`       | `double`          | Evolution parameter                               |

**Preconditions** (checked via `CIRC_VALIDATE`; calls `fprintf(stderr) + exit()`
on failure): `state`, `state->data`, and `diag` must be non-NULL; `state->type`
must be a valid `state_type_t` variant.

---

## Backend Implementations

### Pure (`circ_pure.c`)

Element-wise phase rotation: for each basis state x, multiply `psi[x]` by
`e^{-i diag[x] t} = cos(diag[x]*t) - i sin(diag[x]*t)`.

Complex multiply expanded as real arithmetic:
`(c re + s im) + i(c im - s re)` where `c = cos(angle)`, `s = sin(angle)`.

The loop is trig-bound; `#pragma omp parallel for simd` permits the compiler
to use vector math (libmvec / SVML) when available.  `restrict` qualifiers
on `data` and `diag` enable alias analysis.

OpenMP threshold: 4096.

### Packed (`circ_packed.c`)

Phase precomputation: `cos_d[x] = cos(diag[x]*t)`, `sin_d[x] = sin(diag[x]*t)`
for all `x = 0..dim-1`.  This reduces trig calls from O(dim^2) to O(dim).
The precomputation loop is itself OpenMP-parallelised with `schedule(static)`.

Column iteration: outer loop over columns `c`, inner loop over rows `r = c+1`
to `dim-1`.  Elements within a column are contiguous in packed storage
(index = `col_off(dim, c) + (r - c)`).

Phase factor for element (r, c) reconstructed via the product identity:
`e^{-i(d[r]-d[c])t} = (cos_d[r] - i sin_d[r])(cos_d[c] + i sin_d[c])`.

`#pragma omp simd` on the inner loop enables auto-vectorisation of the
arithmetic (no trig calls in the inner loop).  `schedule(dynamic, 1)` on
the outer loop balances work (column c has dim-c-1 elements).

OpenMP threshold: 512.

### Tiled (`circ_tiled.c`)

Same precomputation strategy as packed.  Iteration follows the tiled
lower-triangular layout: tile pairs (tr, tc) with tr >= tc.

- **Diagonal tiles (tr == tc):** strict lower triangle only, skipping
  diagonal elements (lr > lc, lr != lc).
- **Off-diagonal tiles (tr > tc):** all elements.

`cplx_t * restrict tile = data + tile_off(tr, tc)` eliminates redundant
base-offset additions and enables compiler alias analysis within each tile.

`#pragma omp simd` on inner loops.  `schedule(dynamic, 1)` on tile rows
(row tr processes tr+1 tiles, so later rows have more work).

OpenMP threshold: 512.

---

## Build Integration

The circuit module sources are compiled into the `orkan` shared library
alongside the state, gate, channel, and meas modules.

**Sources** added to `Q_SOURCES` in `src/CMakeLists.txt`:

```cmake
set(CIRC_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/circ/circ.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/circ/circ_pure.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/circ/circ_packed.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/circ/circ_tiled.c"
)
```

**Public header** `include/circ.h` is installed alongside other headers.

**Umbrella header** `include/qlib.h` includes `circ.h`.

**Test target** `test_circ` is registered in `test/CMakeLists.txt` and added
to `verified_install` DEPENDS in the root `CMakeLists.txt`.

---

## Test Suite

18 tests across four categories, run at qubit counts 1 through `MAXQUBITS` (4).
Shared helpers `fill_ramp_diag`, `fill_zero_diag`, and `fill_uniform_diag` live
in `test/include/test_circ.h`.

| Test | State | Observable | Expected |
|------|-------|-----------|----------|
| `test_exp_diag_pure_basis_k` | \|k⟩ | diag[x]=x, t=1 | amp[k] = e^{-ik}, others 0 |
| `test_exp_diag_pure_plus` | \|+⟩^n | diag[x]=x, t=pi/4 | amp[x] = (1/sqrt(dim)) e^{-ix pi/4} |
| `test_exp_diag_pure_zero_obs` | \|+⟩^n | diag[x]=0 | unchanged |
| `test_exp_diag_pure_global_phase` | \|+⟩^n | diag[x]=1, t=pi/3 | all amps *= e^{-i pi/3} |
| `test_exp_diag_pure_unitarity` | \|+⟩^n | diag[x]=x, t=1.23 | sum \|psi[x]\|^2 = 1 |
| `test_exp_diag_packed_basis_k` | \|k⟩⟨k\| | diag[x]=x | unchanged (diagonal state) |
| `test_exp_diag_packed_plus` | \|+⟩⟨+\| | diag[x]=x, t=pi/4 | rho[r,c] = (1/dim) e^{-i(r-c)pi/4} |
| `test_exp_diag_packed_double_apply` | \|+⟩⟨+\| | diag[x]=x, t=pi/4, applied twice | rho[r,c] = (1/dim) e^{-i(r-c)pi/2} |
| `test_exp_diag_packed_zero_obs` | \|+⟩⟨+\| | diag[x]=0 | unchanged |
| `test_exp_diag_packed_identity_obs` | \|+⟩⟨+\| | diag[x]=1 | unchanged (phase diffs zero) |
| `test_exp_diag_tiled_basis_k` | \|k⟩⟨k\| | diag[x]=x | unchanged (diagonal state) |
| `test_exp_diag_tiled_plus` | \|+⟩⟨+\| | diag[x]=x, t=pi/4 | rho[r,c] = (1/dim) e^{-i(r-c)pi/4} |
| `test_exp_diag_tiled_double_apply` | \|+⟩⟨+\| | diag[x]=x, t=pi/4, applied twice | rho[r,c] = (1/dim) e^{-i(r-c)pi/2} |
| `test_exp_diag_tiled_zero_obs` | \|+⟩⟨+\| | diag[x]=0 | unchanged |
| `test_exp_diag_tiled_identity_obs` | \|+⟩⟨+\| | diag[x]=1 | unchanged (phase diffs zero) |
| `test_exp_diag_tiled_small` | n=1 (single tile) | {0,1}, t=pi/2 | rho[1,0] = -0.5i |
| `test_exp_diag_tiled_multi` | n=LOG_TILE_DIM+1 | diag[x]=x, t=pi/4 | element-wise verify |
| `test_exp_diag_consistency` | same state, 3 types | same diag | equal across types |
