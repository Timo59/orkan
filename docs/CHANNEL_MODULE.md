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

Test framework, build-type requirement, and shared fixtures: see `docs/TEST_SUITE.md`. Sweeps cover qubit counts, all target positions, and multiple input states. Reference is computed via BLAS `sum_k K_k rho K_k-dagger` using full n-qubit Kraus matrices, tolerance 1e-12.

### `kraus_to_superop` (`test/channel/test_channel_superop.c`)

| Test | Verifies |
|---|---|
| `test_superop_identity` | Identity channel maps to identity superoperator |
| `test_superop_x_gate` | Single-Kraus X gate matches direct conjugation |
| `test_superop_bitflip_p0` | Bit-flip with p=0 reduces to identity |
| `test_superop_ampdamp_g0` | Amplitude damping with γ=0 reduces to identity |
| `test_superop_ampdamp_g1` | Amplitude damping with γ=1 (full decay) |
| `test_superop_depolarizing_full` | Depolarising channel at full mixing |
| `test_superop_dephasing` | Dephasing channel preserves diagonal, kills off-diagonal |

### MIXED_PACKED (`test/channel/test_channel_packed.c`)

| Test | Verifies |
|---|---|
| `test_channel_1q_identity_packed` | Identity channel leaves state unchanged |
| `test_channel_1q_bitflip_packed` | Bit-flip channel sweep across qubit positions |
| `test_channel_1q_phaseflip_packed` | Phase-flip channel sweep |
| `test_channel_1q_depolar_packed` | Depolarising channel sweep |
| `test_channel_1q_ampdamp_packed` | Amplitude-damping channel sweep |
| `test_ampdamp_g1_decay_packed` | Full-decay (γ=1) drives state to \|0><0\| |
| `test_depolar_full_packed` | Full depolarisation drives state to maximally mixed |
| `test_phaseflip_coherence_packed` | Phase-flip kills off-diagonal coherences |
| `test_bitflip_p1_involution_packed` | Bit-flip at p=1 applied twice = identity |
| `test_trace_preserved_bitflip_packed` | Tr(E(rho)) = Tr(rho) for bit-flip |
| `test_trace_preserved_depolar_packed` | Trace preservation for depolarising |
| `test_trace_preserved_ampdamp_packed` | Trace preservation for amplitude damping |

### MIXED_TILED (`test/channel/test_channel_tiled.c`)

| Test | Verifies |
|---|---|
| `test_channel_1q_identity_tiled` | Identity channel on tiled storage |
| `test_channel_1q_bitflip_tiled` | Bit-flip sweep on tiled |
| `test_channel_1q_phaseflip_tiled` | Phase-flip sweep on tiled |
| `test_channel_1q_depolar_tiled` | Depolarising sweep on tiled |
| `test_channel_1q_ampdamp_tiled` | Amplitude-damping sweep on tiled |
| `test_ampdamp_g1_decay_tiled` | Full decay on tiled |
| `test_depolar_full_tiled` | Full depolarisation on tiled |
| `test_phaseflip_coherence_tiled` | Phase-flip coherence kill on tiled |
| `test_bitflip_p1_involution_tiled` | Bit-flip involution on tiled |
| `test_trace_preserved_bitflip_tiled` | Trace preservation on tiled (bit-flip) |
| `test_trace_preserved_depolar_tiled` | Trace preservation on tiled (depolar) |
| `test_trace_preserved_ampdamp_tiled` | Trace preservation on tiled (ampdamp) |

### Cross-backend (`test/channel/test_channel_tiled.c`)

| Test | Verifies |
|---|---|
| `test_packed_tiled_agree_bitflip` | Element-wise agreement between packed and tiled outputs (bit-flip) |
| `test_packed_tiled_agree_phaseflip` | Cross-backend agreement (phase-flip) |
| `test_packed_tiled_agree_depolar` | Cross-backend agreement (depolarising) |
| `test_packed_tiled_agree_ampdamp` | Cross-backend agreement (amplitude damping) |
