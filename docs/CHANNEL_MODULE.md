# Channel Module — Technical Specification

**Status:** In Progress (1-qubit channels complete; multi-qubit channels not yet implemented)
**Last verified build:** 35 test functions, 0 failures (tolerance 1e-12)

---

## Overview

The channel module simulates local quantum channels given by their Kraus representation on mixed quantum states. A quantum channel E acts on a density matrix as:

E(ρ) = Σ_k K_k ρ K_k†

where {K_k} are the Kraus operators satisfying the completeness relation Σ_k K_k† K_k = I.

The module converts the Kraus representation to a superoperator matrix S such that vec(E(ρ)) = S · vec(ρ), then applies the superoperator element-wise to 2×2 blocks of the density matrix across both packed and tiled storage formats. Only `MIXED_PACKED` and `MIXED_TILED` state types are supported; `PURE` states are rejected at the dispatcher.

---

## Directory Layout

```
include/
  channel.h                    # Public API: kraus_t, superop_t, idx_t typedef,
                               #   kraus_to_superop(), channel_1q(),
                               #   OMP_THRESHOLD, insertBit0, col_off, tile_off, elem_off
                               #   static inlines (shared across backends)

src/channel/
  channel.c                    # kraus_to_superop conversion + channel_1q dispatcher
  channel_packed.c             # MIXED_PACKED backend: channel_packed_1q()
  channel_tiled.c              # MIXED_TILED backend: channel_tiled_1q()

test/channel/
  test_channel.c               # Test runner: main() + 35 RUN_TESTs
  test_channel_superop.c       # Standalone kraus_to_superop unit tests (7 tests)
  test_channel_packed.c        # Packed sweep harness + analytical + trace tests (12 tests)
  test_channel_tiled.c         # Tiled sweep harness + analytical + trace + cross-backend (16 tests)
```

---

## Dependencies

| Dependency | Source | Purpose |
|---|---|---|
| `state_t` and `state_type_t` enum | `include/state.h` | State representation consumed by channel functions |
| `qlib.h` | `include/qlib.h` | `TILE_DIM`, `TILE_SIZE`, `LOG_TILE_DIM`, `POW2` macro |
| OpenMP | System / compiler | Optional parallelisation in packed and tiled backends |

---

## Build Integration

The channel module is compiled as part of `src/CMakeLists.txt` into the `libq` shared library. There is no separate `CMakeLists.txt` inside `src/channel/`.

**Channel sources compiled into `q`:**
```cmake
src/channel/channel.c
src/channel/channel_packed.c
src/channel/channel_tiled.c
```

**Installed header:** `include/channel.h` → `${CMAKE_INSTALL_INCLUDEDIR}`

**Test target:** `test_channel` (executable), defined in `test/CMakeLists.txt`. Links against `q`. Compile definitions: `UNITY_INCLUDE_DOUBLE`, `UNITY_SUPPORT_64`, `LOG_TILE_DIM=${LOG_TILE_DIM}`.

**Run channel tests:**
```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug -R test_channel -V
```

---

## Public Interface

All functions declared in `include/channel.h`. Input validation uses `fprintf(stderr, ...)` + `exit(EXIT_FAILURE)` (consistent with gate module's `GATE_VALIDATE` pattern).

### Types

| Type | Description |
|---|---|
| `idx_t` | `uint64_t` — used for all index computations in the channel module |
| `kraus_t` | Kraus representation: `n_qubits`, `n_terms`, `data` (concatenation of Kraus matrices, row-major) |
| `superop_t` | Superoperator matrix: `n_qubits`, `data` (N²×N² matrix, row-major) |

### Functions

| Function | Signature | Description |
|---|---|---|
| `kraus_to_superop` | `superop_t kraus_to_superop(const kraus_t *kraus)` | Convert Kraus representation to superoperator matrix. Caller must `free(sop.data)`. |
| `channel_1q` | `void channel_1q(state_t *state, const superop_t *sop, qubit_t target)` | Apply 1-qubit channel in-place. Dispatches to packed or tiled backend. |

### Static Inline Helpers (in `channel.h`)

| Helper | Purpose |
|---|---|
| `insertBit0(val, pos)` | Insert a 0 bit at position `pos` — enumerates indices with target bit = 0 |
| `col_off(dim, col)` | Packed column-major offset: index of first element in column `col` |
| `tile_off(tr, tc)` | Tiled lower-triangle offset: index of first element in tile `(tr, tc)` |
| `elem_off(lr, lc)` | Intra-tile offset: element `(lr, lc)` within a `TILE_DIM × TILE_DIM` tile |

---

## Algorithm

### Kraus-to-Superoperator Conversion

For an n-qubit channel with Hilbert space dimension N = 2^n and superoperator dimension M = N²:

```
S_{(i,j),(k,l)} = Σ_a K_a[i,k] · conj(K_a[j,l])
```

The implementation uses four nested loops over `block_row`, `row`, `block_col`, `col` to fill the M×M superoperator matrix. For 1-qubit channels (N=2), this produces a 4×4 matrix. Allocation uses `calloc` for zero-initialised accumulation.

### Shared Macros: `SOP_PRELOAD` and `SOP_APPLY`

Both `channel_tiled.c` and `channel_packed.c` define identical copies of these two macros (file-local, not in the header):

**`SOP_PRELOAD(sop)`**: Decomposes the 16 complex superoperator entries (`sop[0]`..`sop[15]`) into 32 `const double` locals named `_s0r`, `_s0i`, ..., `_s15r`, `_s15i` via `creal()`/`cimag()`. Called once at function entry before the hot loop. The heap pointer `sop` is never touched again after preloading. Named `const double` locals guarantee register promotion; the compiler's register allocator handles spill decisions on individual doubles rather than a monolithic 256-byte heap block.

**`SOP_APPLY(r0,i0, r1,i1, r2,i2, r3,i3, o0r,o0i, o1r,o1i, o2r,o2i, o3r,o3i)`**: Computes the 4×4 superoperator-vector product entirely in real arithmetic. Takes 8 input doubles (real and imaginary parts of 4 density-matrix elements) and writes 8 output doubles. Each of the 8 outputs is a sum of 4 terms of the form `_sXr*(rY) - _sXi*(iY)` (for real parts) or `_sXr*(iY) + _sXi*(rY)` (for imaginary parts), expanding the complex product `(a+bi)(c+di)` into independent real multiply-adds. Total: 64 multiplications and 56 additions per block. Requires `SOP_PRELOAD` in the enclosing scope.

This real-arithmetic decomposition (a) guarantees register promotion of superoperator coefficients, (b) eliminates the cross-lane dependencies inherent in complex-complex multiplication (`ac-bd`, `ad+bc`), and (c) enables SIMD auto-vectorization where all operations are independent per-lane multiply-adds.

At load and store sites, conjugation for Hermitian symmetry is handled by sign-flipping the imaginary part (`-cimag()` at load, `CMPLX(re, -im)` at store), keeping `SOP_APPLY` itself pure (no conjugation logic).

### Packed Backend (`channel_packed_1q`)

Operates directly on lower-triangular packed column-major storage without unpacking. The target qubit divides each density matrix element into a 2×2 block indexed by the target qubit's value in the row and column indices.

**Traversal:** Outer loop over column-base `bc` (OpenMP-parallelised), inner loop over row-base `br ≥ bc`. For each `(br, bc)` pair, four packed indices are computed:
- `idx00 = col_off(dim, c0) + (r0 - c0)` — ρ(r0, c0)
- `idx10 = idx00 + incr` — ρ(r1, c0)
- `idx11 = col_off(dim, c1) + (r1 - c1)` — ρ(r1, c1)
- `idx01` — ρ(r0, c1): branch-dependent (see below)

Three of the four elements (idx00, idx10, idx11) are always in the lower triangle and are read unconditionally before branching. The fourth element (idx01) requires branch-dependent handling:

**Three symmetry branches:**

1. **Branch 1: Lower triangle** (`r0 > c1`): Element ρ(r0, c1) is in the lower triangle. `idx01 = col1_off + (r0 - c1)`. Four reads, `SOP_APPLY`, four writes.
2. **Branch 2: Upper-lower mixed** (`br ≠ bc`, `r0 ≤ c1`): Element ρ(r0, c1) is in the upper triangle. Read via conjugation: `idx01_conj = col_off(dim, r0) + (c1 - r0)`, with `ci1 = -cimag(data[idx01_conj])`. After `SOP_APPLY`, the write-back conjugates the result: `data[idx01_conj] = CMPLX(o1r, -o1i)`.
3. **Branch 3: Diagonal block** (`br == bc`): Element ρ(r0, c1) = conj(ρ(r1, c0)), so `cr1 = cr2, ci1 = -ci2` (reuses already-loaded values with sign flip). Only three elements are written (idx00, idx10, idx11); the fourth (ρ(r0, c1)) is the Hermitian conjugate of the written ρ(r1, c0).

### Tiled Backend (`channel_tiled_1q`)

Operates on lower-triangular tile storage (tiles with `tr ≥ tc` stored, each tile is `TILE_DIM × TILE_DIM` row-major).

**Case 1** (`target < LOG_TILE_DIM`): All four elements of each 2×2 block reside within the same tile. Iterates over all stored tiles; within each tile, the `insertBit0`-based indexing produces non-overlapping 2×2 blocks.

**Case 2** (`target ≥ LOG_TILE_DIM`): The four elements span four different tiles. The target bit is shifted to tile-level coordinates (`target_tile = target - LOG_TILE_DIM`). Three sub-cases:

- **2a** (`tr0 > tc1`): All four tiles (`t00`, `t01`, `t10`, `t11`) are in the stored lower triangle. Flat loop over `TILE_SIZE` elements with `#pragma clang loop vectorize(enable)`. This is the only vectorizable case: the flat loop over contiguous memory with pure-real `SOP_APPLY` arithmetic enables the compiler to emit SIMD instructions (NEON `vld2`/`fmla`/`vst2` on ARM, `vfmadd` on x86).
- **2b** (`btc ≠ btr`, `tr0 ≤ tc1`): Tile `(tr0, tc1)` is above the diagonal; accessed via conjugate-transpose of stored tile `(tc1, tr0)`. Element-wise `i,j` loop with transposed indexing `elem_off(j, i)`. Conjugation at load: `i1 = -cimag(t01[elem_off(j,i)])`. Conjugation at store: `t01[elem_off(j,i)] = CMPLX(o1r, -o1i)`.
- **2c** (`btc == btr`): Diagonal tile pair. Loop `j ≤ i`. Tile `(tr0, tc1)` is the conjugate-transpose of `(tr1, tc0)` (since `tc1 = tr1` and `tc0 = tr0` when `btc == btr`), so reads come from `t10[elem_off(j,i)]` with conjugation: `i1 = -cimag(t10[elem_off(j,i)])`. Write-backs enforce Hermitian symmetry: `t00[elem_off(j,i)] = CMPLX(o0r, -o0i)`, `t10[elem_off(j,i)] = CMPLX(o1r, -o1i)`, `t11[elem_off(j,i)] = CMPLX(o3r, -o3i)`.

---

## Test Infrastructure

### Test Harnesses

| Harness | File | State Type |
|---|---|---|
| `testChannel1qPacked(kraus)` | `test_channel_packed.c` | MIXED_PACKED |
| `testChannel1qTiled(kraus)` | `test_channel_tiled.c` | MIXED_TILED |

Both harnesses sweep `n_qubits` (packed: 1–4, tiled: 1–6), all target positions, and all states from `test_mk_states_mixed`. The reference is computed via `zsumu` (BLAS-based Σ K_k ρ K_k†) using full n-qubit Kraus matrices built by `mat_single_qubit_gate`. The 2×2 Kraus operators are transposed from row-major to column-major before passing to `mat_single_qubit_gate`.

### Test Groups

| Group | Count | Description |
|---|---|---|
| `kraus_to_superop` standalone | 7 | Analytical superoperator verification (identity, X, dephasing, boundary values) |
| Packed sweep | 5 | Full (n, target, state) sweep for identity, bit-flip, phase-flip, depolarising, amplitude damping |
| Tiled sweep | 5 | Same channels on tiled states, n up to 6 |
| Packed analytical / boundary | 4 | Hardcoded: γ=1 decay, p=1 depolarisation, p=0.5 decoherence, p=1 involution |
| Tiled analytical / boundary | 4 | Same tests on tiled states |
| Trace preservation | 6 | Tr(E(ρ)) = Tr(ρ) verified for 3 channels × 2 backends |
| Cross-backend consistency | 4 | Element-wise packed vs tiled agreement for all 4 channels |
| **Total** | **35** | |

### Code Path Coverage (Debug: `LOG_TILE_DIM=3`, `TILE_DIM=8`)

| Code path | First exercised at |
|---|---|
| Packed Branch 1 (`r0 > c1`) | n=3, target=0 |
| Packed Branch 2 (`br≠bc`, `r0≤c1`) | n=2, target=1 |
| Packed Branch 3 (`br==bc`) | n=1 (always) |
| Tiled Case 1 (`target < 3`) | n=1..6, target ∈ {0,1,2} |
| Tiled Case 2a (`tr0 > tc1`) | n=5, target=3 |
| Tiled Case 2b (`btc≠btr`, `tr0≤tc1`) | n=5, target=4 |
| Tiled Case 2c (`btc==btr`) | n=4, target=3 |

### Channel Constructors (test-only, static in each test file)

| Constructor | Kraus Operators |
|---|---|
| `make_bitflip(p)` | K0=√(1-p)·I, K1=√p·X |
| `make_phaseflip(p)` | K0=√(1-p)·I, K1=√p·Z |
| `make_depolarizing(p)` | K0=√(1-3p/4)·I, K1=√(p/4)·X, K2=√(p/4)·Y, K3=√(p/4)·Z |
| `make_ampdamp(γ)` | K0=[[1,0],[0,√(1-γ)]], K1=[[0,√γ],[0,0]] |

---

## Known Issues / Deferred

- **Multi-qubit channels:** Only 1-qubit channels are implemented (`channel_1q`). The `kraus_to_superop` function supports arbitrary `n_qubits`, but no `channel_2q` dispatcher or backend exists.
- **Error handling:** Uses `exit(EXIT_FAILURE)` for all errors, consistent with the gate module. A project-wide transition to `qs_error_t` return codes is deferred.
- **Internal helpers in public header:** `insertBit0`, `col_off`, `tile_off`, `elem_off` are exposed in `channel.h`. These are implementation details that could be moved to an internal header in a future cleanup.
- **`idx_t` naming:** The generic name `idx_t` will be unified with the gate module's `gate_idx_t` in a future refactor.
- **BLAS acceleration:** The 4×4 superoperator application uses real-arithmetic decomposition with `SOP_PRELOAD`/`SOP_APPLY` macros. For multi-qubit channels (larger superoperators), a BLAS-based path (`dgemv`/`dgemm` on the real-decomposed form, or `zgemv`/`zgemm`) would be beneficial.
- **Benchmark integration:** The benchmark suite (`bench_kraus_*`) includes qlib packed and tiled backends (`bench_kraus_qlib.c`), which pre-build the superoperator at init and call `channel_1q()` in the hot path.
