# `sample_matrix` — implementation brief

This document specifies the `sample_matrix` operation that an Orkan
contributor is expected to implement. The contributor has full freedom
in API surface and internal data structures provided the requirements
below are met.

---

## 1. Mathematical definition

Consider a variational ansatz built from `L` channels:

```
|psi> = K_L … K_2 K_1 |iota>,         K_l = sum_m  w_l[m] · M_l[m]
```

Each channel `K_l` is a linear combination of operations `M_l[m]` with
complex weights `w_l[m]`. The weights are the variational parameters;
the operations `M_l[m]` are fixed.

For a chosen channel index `l ∈ [1, L]`, the **sample matrix** at that
channel is the `K_l × K_l` Hermitian matrix

```
S^(l)_{a,b}  =  <psi^(l)_a|  H  |psi^(l)_b>,

|psi^(l)_a>  =  K_L … K_{l+1} · M_l[a] · K_{l-1} … K_1 |iota>,
```

where `H = diag(obs)` is a real diagonal observable in the
computational basis. `S^(l)` is the central object used to update the
weights `w_l` in a generalised eigenproblem step (see §6 for refs).

`H` is real diagonal so `S^(l)` is Hermitian and the diagonal entries
are real.

The implementation must accommodate **arbitrarily many channels per
ansatz**. Initial usage will compute `S` at one or two channels per
optimiser step, but the implementation must not bake in either limit.

---

## 2. Inputs

The function consumes:

- **`|iota>`** — the starting pure state (length `2^n` complex
  amplitudes, `n` qubits).
- **`obs[0..2^n)`** — the diagonal of `H` (real-valued).
- **The ansatz layer specification**, one entry per channel `l`:
  - `M_l[0..K_l)` — array of operations. The representation of an
    operation is the implementer's choice; the existing convention in
    the codebase is `qop_t = void (*)(state_t *)` (see
    `include/state.h`).
  - `w_l[0..K_l)` — complex weight vector for channel `l`. May be
    omitted (NULL or sentinel) at the channel index `l*` for which the
    sample matrix is being computed, since the weights at that channel
    are not used by the formula.
- **The channel index `l*`** at which to compute `S^(l*)`. (If a
  batched-multi-channel API is preferred, accept a list.)
- **A workspace budget** in bytes (a configurable cap on the total
  intermediate state-vector storage the function may hold). May be
  surfaced as a CMake-time constant (see §4) or a runtime parameter.

---

## 3. Output

A column-major `K_{l*} × K_{l*}` complex matrix in a caller-allocated
buffer:

```
S[i + j * K_{l*}] = S^(l*)_{i, j}  (i,j ∈ [0, K_{l*}))
```

Both triangles must be filled (caller convenience). Diagonal entries
must have their imaginary part set to exactly zero (achievable by
construction; assert in debug builds).

---

## 4. Memory-optimality requirements

These are **hard requirements**, not suggestions:

1. **Bounded peak workspace.** The number of full-N (= `2^n`)
   complex-double buffers held simultaneously by the function must be
   bounded by a value derived from the workspace budget. The function
   must not allocate `O(K_{l*})` full state vectors when the budget
   does not permit it.

2. **Prefix sharing.** The state `chi_pre := K_{l*-1} … K_1 |iota>` is
   identical for every column of `S`. It must be computed at most once
   per call. Equivalently, layers below the differentiated channel
   must be applied in a single sweep, not per-leaf.

3. **Streaming when `K_{l*} · 2^n · 16` exceeds the budget.** Process
   the leaves in column blocks of size `B` chosen from the budget.
   Materialise no more than ~`2B` leaf vectors at any moment.

4. **In-place where possible.** The input state buffer may be consumed
   (mutated) by the implementation. Document this clearly in the API
   docstring; callers who must preserve `|iota>` pass a copy.

5. **No per-call allocation in inner loops.** All scratch allocations
   happen at function entry; freed at exit. Inner loops touch only
   already-allocated buffers.

---

## 5. Constraints

- **No `apply_lcu` primitive.** The library no longer exposes a
  generic linear-combination-of-operations primitive. The
  channel-application logic (state copy, op application, weighted
  accumulation) must be inlined inside the sample-matrix
  implementation, sharing the function's scratch buffers.

- **BLAS dependency must remain optional for the rest of the
  library.** Orkan core currently has no BLAS dependency, and this
  must continue to hold for consumers who do not call `sample_matrix`.
  Two acceptable strategies:
  - Confine BLAS-using code to a single translation unit, gated by a
    CMake option (e.g. `BUILD_SAMPLE_MATRIX`, default OFF). Consumers
    who do not enable it never link BLAS.
  - Make BLAS discovery a soft requirement of `sample_matrix` only;
    skip building the symbol if BLAS not found.

  Either way, building Orkan as a submodule of a project that does not
  use `sample_matrix` MUST NOT require an ILP64 BLAS install on the
  consumer's machine.

- **ILP64 when BLAS is used.** When the BLAS-using code is built, the
  BLAS implementation must be ILP64 (64-bit integer ABI). Apple
  Accelerate provides this natively via the
  `ACCELERATE_LAPACK_ILP64` macro on `vecLib/cblas_new.h`. On Linux
  the user supplies an ILP64 OpenBLAS or Intel MKL; symbol-suffix
  variants (Debian's `_64_` suffix) must be handled, either by
  requiring the user to provide a no-suffix build or by a
  configure-time symbol-suffix probe and a generated routing header.
  See `tools/verify_ilp64.c` for an existing standalone ABI test.

- **No bundling.** Do not pull in BLAS via `FetchContent`. The
  consumer is responsible for installing BLAS.

---

## 6. Recommended algorithm

The previous implementation explored, in summary form for the next
implementer:

```
1. Identify channel positions (no special "fixed/branching" flag —
   l* is the only channel that branches; all others apply with their
   given w_l).
2. Apply layers 1..l*-1 to a working copy of |iota> in place.
   The result is chi_pre, kept in the input buffer (slot saving).
3. For each j-block jb in 0, B, 2B, …, K_{l*}:
   a. For each column c in [0, B_j):
        - copy chi_pre into Phi_j[:, c]
        - apply M_{l*}[jb + c] to Phi_j[:, c]
        - apply layers l*+1 … L to Phi_j[:, c] (each as the LCU
          sum_m w_{l'}[m] · M_{l'}[m], inlined using the function's
          two scratch buffers)
   b. Y_j := obs[:, None] * Phi_j  (element-wise scale, parallel)
   c. cblas_zgemm(ConjTrans, NoTrans, B_j, B_j, N, 1, Phi_j, N, Y_j, N,
                  0, &S[jb + jb*K], K)        # diagonal block
   d. for ib in jb+B, jb+2B, …, K_{l*}:
        build Phi_i (B_i columns) into a reused buffer
        cblas_zgemm(ConjTrans, NoTrans, B_i, B_j, N, 1, Phi_i, N, Y_j,
                    N, 0, &S[ib + jb*K], K)   # off-diagonal block
        mirror conj(*) into S[jb + ib*K]
```

Concurrent peak buffers: `2B + small constant` of size `N`.

This is one viable design; alternatives (e.g. Hermitian-rank-K
updates, batched zgemm, GPU offload) are acceptable provided §4 is
respected.

---

## 7. Verification

The implementation MUST ship with tests covering:

- **Correctness against brute force.** For small fixtures
  (`K_{l*} ≤ 8`, `n ≤ 4`), build each `|psi^(l*)_a>` independently
  using only `state_cp`, `qop_t` calls, and a hand-rolled length-N
  reduction for `<a| H |b>`. Compare to `sample_matrix` output to
  machine precision.

- **Hermiticity.** `|S_{i,j} - conj(S_{j,i})| < 1e-12` for every
  off-diagonal entry on a non-trivial fixture.

- **Real diagonal.** `|cimag(S_{i,i})| < 1e-12` for every `i`.

- **Streaming consistency.** A test that runs once with workspace
  budget large enough that `B = K_{l*}` (no streaming) and again with
  a budget forcing `B < K_{l*}`; outputs must be bit-equal.

- **Multi-channel correctness.** A two-channel ansatz with branching
  at position `l*` and non-trivial weights at every other channel;
  compared against brute force.

- **ASan clean.** `cmake -DENABLE_ASAN=ON` build runs the test suite
  with no UBSan/ASan errors.

---

## 8. Documentation deliverables

Alongside the implementation, update:

- **`docs/MEAS_MODULE.md`** (or a new `docs/SAMPLE_MATRIX.md` if the
  function lives in its own module): API table, algorithm sketch,
  parallelism model, build integration (BLAS gating).
- **`docs/BUILD_SYSTEM.md`**: any new CMake options.
- **`CLAUDE.md`** module table: status entry for the new feature.

---

## 9. Out of scope

- Mixed-state generalisation. Pure-state sample matrix is sufficient
  for the immediate use case.
- Non-diagonal observables. `H` is given as `obs[x]` only.
- Resurrection of `apply_lcu` as a public primitive. Its previous
  inclusion was tied to the prior sample-matrix design and has been
  intentionally removed.
