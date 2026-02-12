# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Related Project: Thesis Chapter 2 (Theory)

**Location**: `~/Projects/thesis/2.Simulation/`

This codebase implements the algorithms described in the thesis chapter. When implementing or debugging, agents should reference the theoretical foundations:

| qlib File | Thesis Reference | Theory Provided |
|-----------|------------------|-----------------|
| `src/state.c` | sec4.tex (Eq 6-10, Alg 5) | Packed column major format, offset/stride formulas |
| `src/gate_pure.c` | sec3.tex (Alg 4) | Pure state gate algorithm, invariant subspace decomposition |
| `src/gate_packed_*.c` | sec2.tex (Alg 1-3), sec5.tex (Alg 6-8) | Bit-insertion, hermitian conjugation, Pauli-X optimization |
| `src/gate.c` | sec1.tex (Prop 1-3) | k-local complexity bounds, dispatch strategy |
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
| `src/state.c` | State allocation, initialization, copy |
| `src/gate.c` | Gate dispatchers with input validation |
| `src/gate_pure.c` | Pure state gate implementations |
| `src/gate_packed_1q.c` | Mixed state packed 1-qubit gate implementations |
| `src/gate_packed_<name>.c` | Mixed state packed 2-qubit gates (one file per gate) |
| `src/gate_packed_3q.c` | Mixed state packed 3-qubit gate stubs |
| `src/gate_tiled.c` | Mixed state tiled gate implementations |
| `test/src/test_gate_pure.c` | Pure state test harness |
| `test/src/test_gate_packed.c` | Mixed state test harness |
| `test/src/gatemat.c` | Reference implementations (Kronecker products) |

## Current Status

- State module: complete
- Gates: X, Y, Z, H, S, Sdg, T, Tdg implemented (pure + mixed)
- Rotation gates: Rx, Ry, Rz implemented (pure + mixed)
- Two-qubit gates: CNOT (cx) implemented (pure + mixed)

## Documentation

- `PRD.md` - Architecture, data structures, API design
- `docs/GATES_MODULE.md` - Gate implementation details and test strategy
- `docs/STATE_MODULE.md` - State module specification
- `~/Projects/thesis/2.Simulation/NOTES.md` - Full theoretical background
