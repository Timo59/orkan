# Moment Matrix: Implementation Plan

**Module:** MEAS — `moment_matrix()`
**Status:** Design complete — implementation pending
**Last updated:** 2026-02-28
**Supersedes:** The column-by-column MPI broadcast scheme described in MEAS_MODULE.md §"Computation Algorithm" is replaced by the cyclic distribution scheme in §7 below.

---

## 1. Mathematical Findings

### On the user's notation `tr(E_{j_p}(…E_{j_1}(ρ)))`

This formula produces only the **diagonal** of the moment matrix (M_{II}), not the full
bilinear object M_{IJ}. If every `E_j` is a unitary channel `E_j(ρ) = U_j ρ U_j†`, the
plain trace equals `tr(ρ) = 1` identically. At least the last "channel" must incorporate
`H_C`, giving `M_{II} = Tr(H_C · U_I ρ U_I†)`. The full off-diagonal matrix

    M_{IJ} = Tr(H_C · U_I ρ_0 U_J†)

requires a **bilinear** structure that no single sequential-channel composition can express.

**Conclusion:** The moment matrix computation reduces to the Gram matrix identity

    M_{IJ} = Σ_x h_x φ_I(x) φ_J*(x)   where |φ_I⟩ = U_I|ψ_0⟩

which is a weighted inner product over pure state vectors.

---

## 2. Core Design Decision: Option C (Pure-State Only)

The existing `ρ → U ρ U†` (packed/tiled backends) is **not needed and not used**.
The existing `|ψ⟩ → U|ψ⟩` (PURE state gate operations) is **sufficient**.

| Option | Verdict |
|--------|---------|
| **A: Reuse `ρ → UρU†` only** | ✗ Cannot produce off-diagonal M_{IJ} without the polarization identity (4× cost, requires non-unitary W — incompatible with gate primitives). |
| **B: Add `ρ → Uρ` and `ρ → ρU†`** | ✗ Costs O(N·L·n·D²) flops — a factor of D = 2^n worse than Option C. Only justified for features (quantum instruments, measurement backaction) that are explicitly out of scope. |
| **C: Pure-state only** | ✓ O(N·L·n·D + N²·D) flops, memory O(min(N,B)·D), maps cleanly to BLAS3 ZGEMM. |

**For mixed initial states:** purification `ρ_0 = Σ_k p_k|ψ_k⟩⟨ψ_k|` reduces the
mixed-state problem to `r` independent pure-state runs (r = rank(ρ_0)), each costing the
same as Option C. Since mixed initial states are explicitly deferred in MEAS_MODULE.md,
this path is not implemented now.

---

## 3. Algorithm

The core identity is:

    M = F^H diag(h) F

where F is the D×N matrix with column I equal to `|φ_I⟩ = U_I|ψ_0⟩`.

Two modes are selected at runtime by whether `N · D · 16 bytes` fits in memory.

### Mode A — BLAS (small systems: N·D ≤ 64M complex numbers ≈ 1 GB)

```
1. Allocate F ∈ ℂ^{D×N}
2. Compute all N columns via depth-first tree traversal (§3.1)
3. Form T = diag(h) · F        [element-wise row scaling, O(N·D)]
4. M = ZGEMM(F^H, T)           [N×N result, O(N²·D)]
5. Write lower triangle to file
```

**BLAS call:**
```c
cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans,
            N, N, D, &one, F, D, T, D, &zero, M_dense, N);
```

Note: `h` can be negative (Ising-type cost functions), so ZHERK (which requires a
positive semi-definite kernel) is **not** the right tool — use ZGEMM.

### Mode B — Streaming (large systems: N·D > 1 GB)

```
For each column J = 0..N-1:
  1. Compute |φ_J⟩ via depth-first tree traversal (§3.1)
  2. Compute g_J[x] = h[x] · conj(φ_J[x])        [O(D)]
  3. For each row I = J..N-1:
       a. Compute |φ_I⟩ via tree traversal         [L gate applications]
       b. M[I][J] = Σ_x φ_I(x) · g_J[x]           [O(D) inner product]
       c. M[J][I] = conj(M[I][J])
```

**Recompute cost (without tree):** N·(N+1)/2 full state evolutions × L layers = O(N²·L·D).
**With depth-first tree per column:** O(N²·K/(K−1)·D) gate flops.

*Worked example:* n=20, N=100, L=10, K=2 (binary per layer):
- Tree gate savings: ~1.6×10¹¹ flops
- Inner product cost: ~5×10¹⁰ flops
- **Net overall speedup: ~2.8× vs. naive column-by-column recomputation**

### 3.1 Depth-First Tree Traversal

Multi-indices sharing the same first k layers share the intermediate state
`U_{k;i_k}…U_{1;i_1}|ψ_0⟩`. Depth-first traversal with a stack of L+1 saved state
vectors avoids recomputing these shared prefixes.

**Memory:** (L+1)·D·16 bytes. At n=20, L=10: 176 MB.

**Gate applications** (downward edges only):

    total_edges = Σ_{k=0}^{L-1} (Π_{j=0}^{k} K_j)  =  K·(K^L − 1)/(K − 1)   [uniform K]

vs. naive N·L. Ratio = K/(L·(K−1)). At K=2, L=10: **5× fewer gate applications**.

```c
// Depth-first traversal skeleton
static void dfs(
        const lcu_circuit_t *circuit,
        int level, int prefix_I,
        state_t *stack,          // stack[0..L], pre-allocated
        leaf_fn_t leaf_cb, void *ctx)
{
    if (level == circuit->L) {
        leaf_cb(prefix_I, &stack[level], ctx);
        return;
    }
    lcu_layer_t *layer = circuit->layers[level];
    int stride = circuit->stride[level]; // Π_{j=level+1}^{L-1} K_j, precomputed
    for (int j = 0; j < layer->K; j++) {
        state_cp_inplace(&stack[level + 1], &stack[level]); // save parent
        layer->unitaries[j](&stack[level + 1]);              // apply U
        dfs(circuit, level + 1, prefix_I + j * stride, stack, leaf_cb, ctx);
    }
}
```

`circuit->stride[level]` is precomputed in `lcu_circuit_new` as `Π_{j=level+1}^{L-1} K_j`
(so `stride[L−1] = 1`, `stride[0] = N/K_0`). This matches the least-significant-first
encoding of `lcu_multi_index`: `I = Σ_k i_k · stride[k]`.

### 3.2 Flop / Memory Summary

| | Mode A (BLAS) | Mode B (streaming) |
|---|---|---|
| **Memory for columns** | O(N·D) — all stored | O(L·D) — stack only |
| **Gate flops** | O(N·K/(K-1)·D) — one tree | O(N²·K/(K-1)·D) — tree per column |
| **Assembly flops** | O(N²·D) — ZGEMM | O(N²·D) — inner products |
| **Feasible when** | N·D·16 ≤ 1 GB | Always |

---

## 4. File Structure

```
Prerequisites (Circuit module — must complete first):
  include/circuit.h          ← Add lcu_layer_t, lcu_circuit_t, lcu_multi_index decl
                               Add stride[] field to lcu_circuit_t
  src/circuit/circuit.c      ← Implement lcu_layer_new/free, lcu_circuit_new/free,
                                 lcu_multi_index, precompute stride[] in constructor

New files (MEAS module):
  include/meas.h             ← diag_hamiltonian_t, moment_matrix(), moment_matrix_load()
  src/meas/moment_matrix.c   ← Full implementation: dfs_tree, Mode A, Mode B, file I/O
  test/meas/test_moment_matrix.c

CMake:
  src/CMakeLists.txt         ← add_library(meas STATIC ...), link circuit q BLAS
  test/CMakeLists.txt        ← add_executable(test_moment_matrix ...)
```

---

## 5. Public Interface

```c
// include/meas.h

typedef struct {
    double  *h;       // diagonal of H_C, length 2^n_qubits, caller-owned
    int      n_qubits;
} diag_hamiltonian_t;

// Compute and write the moment matrix (pure initial state only).
// initial->type must be PURE. Writes lower triangle to output_path (MMAT format).
// Returns 0 on success, -1 on I/O error (errno set).
// Programming errors (NULL pointers, wrong state type) → assert.
int moment_matrix(
    const lcu_circuit_t      *circuit,
    const state_t            *initial,
    const diag_hamiltonian_t *H_C,
    const char               *output_path
);

// Load a previously written moment matrix. Allocates *M_out (caller must free).
// Returns N (number of multi-indices) on success, -1 on error.
int moment_matrix_load(
    const char     *path,
    int            *N_out,
    double complex **M_out   // lower triangle, N*(N+1)/2 elements
);
```

---

## 6. Internal Structure of `moment_matrix.c`

```c
#define MMAT_STORE_THRESHOLD  (1LL << 26)   // 64M complex numbers ≈ 1 GB

// --- Tree traversal helpers ---
typedef void (*leaf_fn_t)(int I, const state_t *phi, void *ctx);

static void dfs(const lcu_circuit_t*, int level, int prefix_I,
                state_t *stack, leaf_fn_t cb, void *ctx);

// --- Mode A helpers ---
// Leaf callback: writes phi into column I of F
static void cb_store_column(int I, const state_t *phi, void *ctx);

static int moment_matrix_blas(const lcu_circuit_t*, const state_t *initial,
                               const diag_hamiltonian_t*, const char *path);

// --- Mode B helpers ---
static double complex inner_h(const state_t *phi_I, const state_t *phi_J,
                               const double *h, dim_t D);

static int moment_matrix_stream(const lcu_circuit_t*, const state_t *initial,
                                const diag_hamiltonian_t*, const char *path);

// --- Dispatch ---
int moment_matrix(const lcu_circuit_t *circuit, const state_t *initial,
                  const diag_hamiltonian_t *H_C, const char *output_path)
{
    assert(initial && initial->type == PURE && circuit && H_C && output_path);
    dim_t D = 1LL << H_C->n_qubits;
    long long N = circuit->N;
    if (N * D <= MMAT_STORE_THRESHOLD)
        return moment_matrix_blas(circuit, initial, H_C, output_path);
    else
        return moment_matrix_stream(circuit, initial, H_C, output_path);
}
```

---

## 7. MPI Distribution Strategy (ENABLE_MPI=ON)

**Replaces** the column-by-column broadcast scheme in MEAS_MODULE.md.

**Use cyclic column distribution:**

- Rank r owns columns J where `J % n_ranks == r`
- Each rank independently computes its columns and all corresponding rows (I ≥ J)
- No `MPI_Bcast` — each rank traverses the tree autonomously
- Gather M entries to rank 0 via `MPI_Gatherv` (data volume: N²·16 bytes, negligible)
- **Load balance:** cyclic assignment distributes row counts uniformly (rank r gets a mix
  of small and large J values), avoiding the ~(2P−1)× imbalance of block column assignment

**Why not the broadcast scheme from MEAS_MODULE.md:**
At n=20, N=100, 10 GB/s interconnect: N broadcasts × 16 MB = 160 ms of pure
communication overhead added to ~660 ms of compute — a 24% tax with no compensating
benefit. At n=25, each broadcast is 512 MB; the tax grows to 5 s per run. Cyclic
distribution eliminates this entirely.

---

## 8. Implementation Order

```
Step 1  Add stride[] to lcu_circuit_t; implement lcu_layer_new/free,
        lcu_circuit_new/free, lcu_multi_index in src/circuit/circuit.c.
        Add all types to include/circuit.h.

Step 2  Create include/meas.h (diag_hamiltonian_t, declarations only).

Step 3  Implement src/meas/moment_matrix.c:
          3a. dfs() traversal helper (no BLAS dependency, unit-testable)
          3b. Mode B (streaming) — validate on small n immediately
          3c. Mode A (BLAS ZGEMM path) — add once Mode B passes

Step 4  Wire CMake: meas static library, test_moment_matrix executable.

Step 5  Write test/meas/test_moment_matrix.c:
          • n=2, N=4 (2 layers × 2 unitaries): brute-force reference comparison
          • Hermitian symmetry: M[I][J] == conj(M[J][I]) for all I,J
          • Diagonal: M[II] equals direct expectation value Tr(H_C · U_I ρ U_I†)
          • Mode A vs Mode B agreement at small n
          • File round-trip: write then moment_matrix_load, compare

Step 6  (deferred) MPI cyclic distribution path.
```

---

## 9. What is NOT Needed

| Item | Decision | Rationale |
|------|----------|-----------|
| New `ρ → Uρ` primitive | ✗ Skip | Pure states only; off-diagonal M_{IJ} handled by inner products over state vectors |
| New `ρ → ρU†` primitive | ✗ Skip | Same reason |
| Polarization identity via `ρ → UρU†` | ✗ Skip | 4× cost, requires non-unitary W, incompatible with gate module's unitary assumption |
| Mixed initial state support | ✗ Deferred | Explicitly out of scope per MEAS_MODULE.md |
| Column-by-column MPI broadcast | ✗ Replaced | Cyclic column distribution is simpler, faster, and better load-balanced (§7) |
| ZHERK instead of ZGEMM | ✗ Use ZGEMM | h can be negative; ZHERK assumes positive semi-definite kernel |

The existing `ρ → UρU†` implementations (packed and tiled backends) are **not modified**.
All moment matrix computation operates on PURE state vectors exclusively.
