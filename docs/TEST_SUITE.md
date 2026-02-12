# Test Suite

## Framework

[Unity](https://github.com/throwtheswitch/unity) (C unit testing), included as a git submodule in `extern/Unity/`.

## Directory Structure

```
test/
├── CMakeLists.txt          Build definitions for both test executables
├── include/                Shared test headers
│   ├── test.h              Common test macros and includes
│   ├── test_gate.h         Gate test declarations
│   ├── test_pure_states.h  Pure state fixture declarations
│   ├── test_mixed_states.h Mixed state fixture declarations
│   ├── test_mixed_utils.h  Packed density matrix conversion utilities
│   ├── gatemat.h           Reference gate matrix declarations
│   └── linalg.h            Linear algebra helper declarations
├── state/                  State module tests (→ test_state executable)
│   ├── test_state.c        Runner (main)
│   ├── test_state_pure.c   Pure statevector tests
│   ├── test_state_packed.c Packed density matrix tests
│   └── test_state_tiled.c  Tiled density matrix tests
├── gate/                   Gate module tests (→ test_gate executable)
│   ├── test_gate.c         Runner (main)
│   ├── test_gate_pure_1q.c Pure state single-qubit gate tests
│   ├── test_gate_pure_2q.c Pure state two-qubit gate tests
│   ├── test_gate_packed_1q.c Packed mixed state single-qubit gate tests
│   ├── test_gate_packed_2q.c Packed mixed state two-qubit gate tests
│   ├── test_gate_packed_3q.c Packed mixed state three-qubit gate tests
│   └── test_gate_tiled.c   Tiled mixed state gate tests
└── utility/                Test helpers and fixtures
    ├── gatemat.c           Reference gate matrices (Kronecker products)
    ├── linalg.c            Linear algebra helpers
    ├── test_pure_states.c  Pre-built pure state fixtures
    ├── test_mixed_states.c Pre-built mixed state fixtures
    └── test_mixed_utils.c  Packed density matrix conversion utilities
```

## Executables

| Executable   | Sources                        | Description |
|------------- |--------------------------------|-------------|
| `test_state` | `state/*` + Unity              | State allocation, initialization, copy, get/set |
| `test_gate`  | `gate/*` + `utility/*` + Unity | Gate application across all representations |

## Running

```bash
cd cmake-build-debug
ctest                    # Run all tests
./test/test_state        # State module only
./test/test_gate         # Gate module only
```
