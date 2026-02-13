# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Related Project: Thesis Chapter 2 (Theory)

**Location**: `~/Projects/thesis/2.Simulation/`

This codebase implements the algorithms described in the thesis chapter. When implementing or debugging, agents should reference the theoretical foundations:

| qlib File | Thesis Reference | Theory Provided |
|-----------|------------------|-----------------|
| `src/state/state.c` | sec4.tex (Eq 6-10, Alg 5) | Packed column major format, offset/stride formulas |
| `src/gate/gate_pure.c` | sec3.tex (Alg 4) | Pure state gate algorithm, invariant subspace decomposition |
| `src/gate/packed/gate_packed_*.c` | sec2.tex (Alg 1-3), sec5.tex (Alg 6-8) | Bit-insertion, hermitian conjugation, Pauli-X optimization |
| `src/gate/gate.c` | sec1.tex (Prop 1-3) | k-local complexity bounds, dispatch strategy |
| `include/state.h` | sec4.tex | Memory layout specification |

**Key thesis results for implementers**:
- **Proposition 3** (sec1.tex): O(r·N²·2^k) complexity bound justifies per-gate optimization
- **Algorithm 1** (sec2.tex): `insertBits()` implementation spec
- **Equation 12** (sec5.tex): Update formula for single-qubit conjugation
- **Zero-flop gates**: Pauli-X, CNOT require only memory swaps (sec1.tex Example 1, sec2.tex Alg 3)

**To work on thesis**: `cd ~/Projects/thesis/2.Simulation && cat CLAUDE.md`

## Build & Test

```bash
mkdir cmake-build-debug && cd cmake-build-debug
cmake ..
cmake --build .
ctest                    # Run all tests
./test/test_state        # State module tests (15 tests)
./test/test_gate         # Gate module tests (~11,800 per gate)
```

First Linux build takes 10-30 min (compiles OpenBLAS with ILP64).

## Key Source Files

| File | Purpose |
|------|---------|
| `src/state/` | State allocation, initialization, copy |
| `src/gate/gate.c` | Gate dispatchers with input validation |
| `src/gate/gate_pure.c` | Pure state gate implementations |
| `src/gate/packed/gate_packed_1q.c` | Mixed state packed 1-qubit gate implementations |
| `src/gate/packed/gate_packed_<name>.c` | Mixed state packed 2-qubit gates (one file per gate) |
| `src/gate/tiled/gate_tiled.c` | Mixed state tiled 1Q + rotation gate implementations |
| `src/gate/tiled/gate_tiled_<name>.c` | Mixed state tiled 2Q gates (one file per gate) |

## Current Status

- State module: complete
- Single-qubit gates: X, Y, Z, H, Hy, S, Sdg, T, Tdg, P, Rx, Ry, Rz (pure + packed; tiled under construction)
- Two-qubit gates: CX, CY, CZ, SWAP (pure + packed; tiled under construction for CX, SWAP)
- Three-qubit gate: CCX (pure + packed)

## Documentation

- `PRD.md` - Architecture, data structures, API design
- `docs/GATES_MODULE.md` - Gate implementation details and test strategy
- `docs/STATE_MODULE.md` - State module specification
- `docs/TEST_SUITE.md` - Test directory structure and how to run tests (read only when working on tests)
- `~/Projects/thesis/2.Simulation/NOTES.md` - Full theoretical background
