# Persistent Agent Memory — qlib HPC Engineer

## Project Quick Reference

- **Build:** `cmake --preset debug && cmake --build --preset debug && ctest --preset debug`
- **Debug preset mandatory:** sets `LOG_TILE_DIM=3` (8x8 tiles) for multi-tile test coverage
- **Key headers:** `include/state.h`, `include/gate.h`, `include/circuit.h`
- **Key sources:** `src/circuit/pauli.c` (Pauli strings + apply + exp_pure)

## Module Status (as of 2026-02-27)

| Module | Status |
|--------|--------|
| State | Complete |
| Gate | In Progress |
| Circuit: pauli construction/phase/apply | Complete |
| Circuit: `pauli_exp_pure` | **Complete** (added to `pauli.c`) |
| Circuit: `pauli_exp_packed/tiled`, `grad_inner` | Pending |
| Circuit: channel/pqc/lcu types | Pending |

## Pauli String Conventions

- `y_mask != 0` structurally implies `x_mask != 0` (Y sets both bits in `pauli_new`)
- Phase formula: `eps(x) = i^{popcount(y_mask)} * (-1)^{popcount(x & z_mask)}`
- Exponent: `exp_j = (eta_exp + 2*parity_j) & 3`, `exp_jp = (exp_j + 2*y_parity) & 3`
- Pair enumeration: `split_bit = __builtin_ctzll(x_mask)`, ONE insertBit0 only — NOT at every x_mask bit

## pauli_exp_pure Algorithm (module spec convention)

`U(theta) = exp(-i*theta/2 * P) = cos(theta/2)*I - i*sin(theta/2)*P`

Key trick: `-i * i^k = i^{(k+3)&3}` — avoids all complex multiplication.

- **Identity** (x_mask==0, z_mask==0): global phase `CMPLX(c, -s)` on all amplitudes
- **Pure-Z diagonal** (y_mask==0): two precomputed scalars `CMPLX(c,-s)` and `CMPLX(c,+s)`
- **General diagonal** (dead code): `c*data[i] + s*apply_ipow(data[i], (exp_i+3)&3)`
- **Non-diagonal**: read both old values first, then:
  ```c
  data[j]  = c*old_j  + s*apply_ipow(old_jp, (exp_jp+3)&3)
  data[jp] = c*old_jp + s*apply_ipow(old_j,  (exp_j +3)&3)
  ```

See `docs/pauli_exp_notes.md` for full derivation details.
