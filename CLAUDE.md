# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

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
| `src/qhipster.c` | Pure state gate implementations |
| `src/mhipster.c` | Mixed state gate implementations |
| `test/src/test_qhipster.c` | Pure state test harness |
| `test/src/test_mhipster.c` | Mixed state test harness |
| `test/src/gatemat.c` | Reference implementations (Kronecker products) |

## Current Status

- State module: complete
- Gates: X, Y, Z implemented (pure + mixed)
- In progress: H, S, T, rotation gates

## Documentation

- `PRD.md` - Architecture, data structures, API design
- `docs/GATES_MODULE.md` - Gate implementation details and test strategy
- `docs/STATE_MODULE.md` - State module specification
