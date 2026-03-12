# Qiskit-Aer DensityMatrix API Surface

Discovered by iterative compilation of a probe program against the vendored headers.
Commit: c5841316078b2f141418459faf12e67070888fd7

## Constructor Convention

```cpp
AER::QV::DensityMatrix<double> dm(n);
```

`n` is the number of **logical qubits**. The class stores a 2^n x 2^n density
matrix as a column-stacked 2^(2n)-element complex vector. Internally it
delegates to `UnitaryMatrix<double>(n)`, which calls
`BaseVector::set_num_qubits(2 * n)`.

The constructor allocates memory but does **not** initialize to any state.
You must call `dm.initialize()` explicitly to set the |0><0| state.

## Data Accessor

```cpp
std::complex<double> *p = dm.data();   // non-const
const std::complex<double> *p = dm.data();  // const overload
```

- Returns a raw `std::complex<double>*` (i.e., `std::complex<data_t>*`).
- Pointer is valid until the object is destroyed or `set_num_qubits` is called.
- Element count: `dm.size()` returns `1ULL << (2 * n)`.
  For n=4 qubits: 256 elements.
- Layout: column-major vectorization of the density matrix.
  Element `(row, col)` is at `data_[row + col * 2^n]`.
- `data()[0]` is `rho[0,0]` = 1.0+0i after `dm.initialize()` for |0><0|.

## sizeof

```
sizeof(DensityMatrix<double>) = 144   (on arm64 macOS)
```

This is the struct size only; the actual matrix data lives on the heap.

## Gate Methods

### Methods with dedicated implementations on DensityMatrix

| Gate   | Method signature |
|--------|-----------------|
| X      | `dm.apply_x(uint_t qubit)` |
| Y      | `dm.apply_y(uint_t qubit)` |
| CNOT   | `dm.apply_cnot(uint_t qctrl, uint_t qtrgt)` |
| CY     | `dm.apply_cy(uint_t qctrl, uint_t qtrgt)` |
| SWAP   | `dm.apply_swap(uint_t q0, uint_t q1)` |
| ECR    | `dm.apply_ecr(uint_t q0, uint_t q1)` |
| Phase  | `dm.apply_phase(uint_t q, const complex_t &phase)` |
| CPhase | `dm.apply_cphase(uint_t q0, uint_t q1, const complex_t &phase)` |
| Toffoli| `dm.apply_toffoli(uint_t qctrl0, uint_t qctrl1, uint_t qtrgt)` |
| Reset  | `dm.apply_reset(const reg_t &qubits)` |

### H, Z and general single/multi-qubit unitaries: use apply_unitary_matrix

H and Z do NOT have dedicated `apply_h` / `apply_z` methods.
Apply them via the general unitary interface:

```cpp
#include "framework/linalg/matrix_utils.hpp"   // for AER::Linalg::Matrix::H etc.

// H gate on qubit 0
auto h_vec = AER::Utils::vectorize_matrix(AER::Linalg::Matrix::H);
dm.apply_unitary_matrix({0}, h_vec);  // reg_t = std::vector<uint_t>

// Z gate on qubit 1
auto z_vec = AER::Utils::vectorize_matrix(AER::Linalg::Matrix::Z);
dm.apply_unitary_matrix({1}, z_vec);
```

`AER::Linalg::Matrix` provides pre-built gate matrices (type `AER::cmatrix_t`):
`I, X, Y, Z, H, S, SDG, T, TDG, CX, CY, CZ, SWAP, ...`

`AER::Utils::vectorize_matrix` column-stacks the matrix into `cvector_t<double>`.

`apply_unitary_matrix` signature:
```cpp
void apply_unitary_matrix(const reg_t &qubits, const cvector_t<double> &mat);
```
where `reg_t = std::vector<uint_t>` and `cvector_t<double> = std::vector<std::complex<double>>`.

### Superoperator (Kraus/channel) interface

```cpp
void apply_superop_matrix(const reg_t &qubits, const cvector_t<double> &mat);
void apply_diagonal_unitary_matrix(const reg_t &qubits, const cvector_t<double> &diag);
void apply_diagonal_superop_matrix(const reg_t &qubits, const cvector_t<double> &diag);
```

## Required Compilation Flags

```
-std=c++17
-I extern/aer-dm/include
-I <nlohmann-json-prefix>/include    # e.g. $(brew --prefix nlohmann-json)/include
-I <openmp-include>                  # for omp.h (from clifford.hpp via operations.hpp)
-L <openmp-lib-dir> -lomp            # or -fopenmp on GCC/Linux
```

On macOS (arm64):
```
-I /opt/homebrew/opt/nlohmann-json/include
-I /opt/homebrew/Caskroom/miniconda/base/include
-L /opt/homebrew/Cellar/libomp/21.1.8/lib -lomp
```

Note: OpenMP is required at link time because `clifford.hpp` (pulled in
transitively via `operations.hpp` -> `qubitvector.hpp` -> `densitymatrix.hpp`)
has inline implementations that call `omp_get_num_threads`.

## Vendor Patch Applied

`densitymatrix.hpp` was patched to add an explicit include of
`framework/linalg/matrix_utils.hpp` (provides `Linalg::SMatrix` used by
`apply_reset`). The upstream file is missing this include and relies on
`densitymatrix_state.hpp` pulling it in transitively — that include chain is
not present when using the DensityMatrix header standalone.

Diff:
```
--- densitymatrix.hpp (upstream)
+++ densitymatrix.hpp (vendored)
@@ -18,6 +18,7 @@
 #include "framework/utils.hpp"
+#include "framework/linalg/matrix_utils.hpp"
 #include "simulators/unitary/unitarymatrix.hpp"
```

## Headers Copied

All are from `src/` in the Qiskit-Aer repository at the commit above.

```
framework/avx2_detect.hpp
framework/blas_protos.hpp
framework/circuit.hpp
framework/config.hpp
framework/creg.hpp
framework/json.hpp
framework/json_parser.hpp
framework/lapack_protos.hpp
framework/linalg/almost_equal.hpp
framework/linalg/eigensystem.hpp
framework/linalg/enable_if_numeric.hpp
framework/linalg/linalg.hpp
framework/linalg/linops/linops_aer_vector.hpp
framework/linalg/linops/linops_array.hpp
framework/linalg/linops/linops_generic.hpp
framework/linalg/linops/linops_json.hpp
framework/linalg/linops/linops_map.hpp
framework/linalg/linops/linops_matrix.hpp
framework/linalg/linops/linops_unordered_map.hpp
framework/linalg/linops/linops_vector.hpp
framework/linalg/matrix_utils.hpp
framework/linalg/matrix_utils/matrix_defs.hpp
framework/linalg/matrix_utils/smatrix_defs.hpp
framework/linalg/matrix_utils/vmatrix_defs.hpp
framework/linalg/square.hpp
framework/linalg/vector.hpp
framework/linalg/vector_json.hpp
framework/matrix.hpp
framework/noise_utils.hpp
framework/operations.hpp
framework/opset.hpp
framework/rng.hpp
framework/stl_ostream.hpp
framework/types.hpp
framework/utils.hpp
misc/common_macros.hpp
misc/hacks.hpp
misc/warnings.hpp
simulators/density_matrix/densitymatrix.hpp       (patched — see above)
simulators/stabilizer/binary_vector.hpp
simulators/stabilizer/clifford.hpp
simulators/stabilizer/pauli.hpp
simulators/statevector/indexes.hpp
simulators/statevector/qubitvector.hpp
simulators/statevector/qv_avx2.hpp
simulators/statevector/transformer.hpp
simulators/statevector/transformer_avx2.hpp
simulators/unitary/unitarymatrix.hpp
```

## External Dependencies NOT Vendored

- **nlohmann/json** (3.x): system package or `brew install nlohmann-json`
- **omp.h / libomp**: system OpenMP; on macOS `brew install libomp`
