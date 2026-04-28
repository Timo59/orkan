# Circuit Module

Parameterised circuit operations on quantum states. Currently provides `exp_diag` (diagonal unitary evolution) across PURE, MIXED_PACKED, and MIXED_TILED.

---

## Directory Layout

```
include/
  circ.h                  Public API

src/circ/
  circ.c                  Dispatch layer
  circ_pure.c             PURE backend
  circ_packed.c           MIXED_PACKED backend
  circ_tiled.c            MIXED_TILED backend

test/circ/
  test_circ.c             Test runner
  test_circ_pure.c        PURE tests
  test_circ_packed.c      MIXED_PACKED tests
  test_circ_tiled.c       MIXED_TILED tests
```

---

## Public API

Declared in `include/circ.h`.

### `exp_diag`

```c
void exp_diag(state_t *state, const double *diag, double t);
```

Applies the diagonal unitary `U = exp(-i D t)` in-place:

- **Pure:**  `psi[x] *= exp(-i diag[x] t)`
- **Mixed:** `rho[r][c] *= exp(-i (diag[r] - diag[c]) t)`

Diagonal density-matrix elements are unchanged (phase difference is zero).

No implicit factor of 1/2. To match the QAOA convention `e^{-i gamma C}`, pass `gamma` directly. For the mixer layer with RX gates, pass `2*beta`.

| Parameter | Type | Description |
|---|---|---|
| `state` | `state_t *` | Quantum state (any representation), modified in-place |
| `diag` | `const double *` | Diagonal coefficients, length `2^(state->qubits)` |
| `t` | `double` | Evolution parameter |

Input validation uses `CIRC_VALIDATE` (in `circ.c`): non-NULL `state`, `state->data`, `diag`; valid `state->type`. Failures print to stderr and call `exit()`.

---

## Architecture

### Pure — `exp_diag_pure`

Element-wise phase rotation: for each basis state `x`, multiply `psi[x]` by `e^{-i diag[x] t} = cos(diag[x]*t) - i sin(diag[x]*t)`.

Complex multiply expanded as real arithmetic: `(c·re + s·im) + i(c·im - s·re)` where `c = cos(angle)`, `s = sin(angle)`. The loop is trig-bound; `#pragma omp parallel for simd` permits vector-math libraries (libmvec / SVML). `restrict` qualifiers on `data` and `diag` enable alias analysis.

OpenMP threshold: 4096.

### Packed — `exp_diag_packed`

Phase precomputation: `cos_d[x] = cos(diag[x]*t)`, `sin_d[x] = sin(diag[x]*t)` for all `x = 0..dim-1`. Reduces trig calls from O(dim^2) to O(dim). The precomputation loop is OpenMP-parallelised with `schedule(static)`.

Column iteration: outer over columns `c`, inner over rows `r = c+1..dim-1`. Elements within a column are contiguous (`index = col_off(dim, c) + (r - c)`).

Phase factor for element `(r, c)` reconstructed via `e^{-i(d[r]-d[c])t} = (cos_d[r] - i sin_d[r])(cos_d[c] + i sin_d[c])`.

`#pragma omp simd` on the inner loop. `schedule(dynamic, 1)` on the outer loop balances column work (column c has dim-c-1 elements).

OpenMP threshold: 512.

### Tiled — `exp_diag_tiled`

Same precomputation strategy. Iteration follows the tiled lower-triangular layout: tile pairs `(tr, tc)` with `tr >= tc`.

- **Diagonal tiles** (`tr == tc`): strict lower triangle only, skipping the diagonal (lr > lc, lr != lc).
- **Off-diagonal tiles** (`tr > tc`): all elements.

`cplx_t * restrict tile = data + tile_off(tr, tc)` eliminates redundant base-offset additions and enables compiler alias analysis within each tile.

`#pragma omp simd` on inner loops. `schedule(dynamic, 1)` on tile rows (row `tr` processes `tr+1` tiles).

OpenMP threshold: 512.

---

## Build Integration

Sources compiled into `${PROJECT_NAME_LOWER}` via `src/CMakeLists.txt`:

```cmake
set(CIRC_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/circ/circ.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/circ/circ_pure.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/circ/circ_packed.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/circ/circ_tiled.c"
)
```

Public header `include/circ.h` is installed alongside other module headers. Umbrella header `include/orkan.h` includes `circ.h`.

Test target `test_circ` is registered in `test/CMakeLists.txt` and added to `verified_install` DEPENDS in the root `CMakeLists.txt`.

```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug -R test_circ
```

---

## Tests

| File | Coverage |
|---|---|
| `test/circ/test_circ.c` | Runner |
| `test/circ/test_circ_pure.c` | PURE: basis states, \|+>^n, zero/identity diag, unitarity |
| `test/circ/test_circ_packed.c` | MIXED_PACKED: \|+><+\|, double-apply, zero/identity diag |
| `test/circ/test_circ_tiled.c` | MIXED_TILED: same coverage + small-state (single-tile) + multi-tile |

Shared helpers `fill_ramp_diag`, `fill_zero_diag`, `fill_uniform_diag` live in `test/include/test_circ.h`. Tests run at qubit counts 1 through `MAXQUBITS` (4). Reference is element-wise complex exponential at tolerance 1e-12.
