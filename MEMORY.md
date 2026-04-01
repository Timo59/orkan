# Agent Memory - qlib Project

- [feedback_approve_edits.md](memory/feedback_approve_edits.md) — Always wait for user approval before writing/editing any file

## Key Architecture Facts

### Tiled Gate Layout (Debug build: LOG_TILE_DIM=3, TILE_DIM=8)
- Tiles are 8x8, TILE_SIZE=64 elements stored row-major (NOT column-major like LAPACK)
- Only lower-triangular tiles are stored at tile level
- With MAXQUBITS_TILED=6, dim_tile ranges from 1 (n=2,3) to 8 (n=6)
- OMP_THRESHOLD=512: OpenMP only fires for dim>=512, i.e., n>=9 for state vectors

### Case Routing in Two-Qubit Tiled Gates (CX/CY/CZ/SWAP)
- Case 1: both q1,q2 < LOG_TILE_DIM (within-tile)
- Case 2: control < LOG_TILE_DIM, target >= LOG_TILE_DIM (mixed)
- Case 3: control >= LOG_TILE_DIM, target < LOG_TILE_DIM (mixed, reversed)
- Case 4: both q1,q2 >= LOG_TILE_DIM (tile-level)

### DISPATCH_2Q vs DISPATCH_2Q_NO_TILED
- cy() and cz() use DISPATCH_2Q (full tiled dispatch) — tiled implementations EXIST
- The description "DISPATCH_2Q_NO_TILED" in the analysis prompt was incorrect
- gate.c line 179-180 confirms: cy and cz both use DISPATCH_2Q

### test_mk_states_mixed() Contents
- 3*dim computational/X/Y basis states + 1 maximally mixed + 10 random + entangled
- For n=2: additionally 4 Bell states; for n>2: GHZ + W states

### GATE_VALIDATE fires exit() — not an assertion
- If CY/CZ were unimplemented for tiled, calling them would kill the test process
- Since they ARE implemented, test_CY_tiled and test_CZ_tiled are valid tests

## Circuit Module — pauli.c Test Suite (as of 2026-02-27)

- **62 tests total** in `test/circuit/test_pauli.c` / `test_circuit.c`, 0 failures
- `pauli_new` / `pauli_from_str` / `pauli_free`: 13 tests — includes n=0, NULL labels, n>64, n=64 boundary, free(NULL), case-sensitivity, empty string, mask invariant
- `pauli_phase`: 3 tests — all 4 eta_exp values (0,1,2,3) covered via X/Y/YY/YYY
- `pauli_apply_pure`: 43 tests — n=1 to n=4, all split_bit values {0,1,2,3}, YYY (eta_exp=3), YZ (mixed), Z-involution on non-trivial states
- `pauli_exp_pure` convention: `exp(-i*theta/2 * P)` = `cos(theta/2)*I - i*sin(theta/2)*P`
  - theta=0: identity; theta=pi: -i*P; theta=2*pi: -I; theta=4*pi: +I (true periodicity)
  - XX gives `-i*s|11>`, YY gives `+i*s|11>` (different sign due to eta_exp difference)
  - 19 tests covering all branches, all eta_exp values, split_bit {0,1,2}, unitarity, inverse, eigenstate phases, non-unit-norm

## Known Coverage Gaps in testTwoQubitGateTiled (see full analysis)
- OMP_THRESHOLD (dim>=512) never reached: MAXQUBITS_TILED=6 gives dim=64 max
- Involution property (gate applied twice = identity) not tested
- Cleanup path (goto cleanup) is only reachable on malloc failure — not tested
- Debug diagnostic block (goto done_check) is redundant: it prints first mismatch
  but does NOT call TEST_FAIL; assert_tiled_equals_full does the actual assertion
- State count computation in test_rm_states_mixed duplicates logic from test_mk_states_mixed
  (fragile, not a test gap but a maintenance hazard)
