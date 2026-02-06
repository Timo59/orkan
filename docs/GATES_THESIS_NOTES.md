# Gate Algorithms: Thesis-Level Analysis

Supporting material for the thesis chapter on classical simulation of quantum circuits.
Central question: **how do algorithmic design choices interact with hardware characteristics
to generate performance in quantum gate simulation?**

Organized to follow the chapter's logical arc: index arithmetic and vectorization (how
performance arises from the pure-state algorithm), parallelization (how it scales), then
mixed-state complications (how storage format determines algorithmic structure and cost).

---

## 1. Index Arithmetic: `insertBit0` and Stride-Based Enumeration

A single-qubit gate U acting on qubit q of an n-qubit state pairs amplitudes whose indices
differ only in bit q. Enumerating these pairs efficiently is the fundamental operation of the
simulator. Two approaches exist.

**`insertBit0(val, pos)`** inserts a 0-bit at position `pos` in the binary representation of
`val`, shifting all higher bits left by one:

```
insertBit0(val, pos) = (val >> pos) << (pos + 1) | (val & ((1 << pos) - 1))
```

This maps a dense loop counter `val ∈ [0, dim/2)` to the sparse set of array indices with
bit q = 0, without wasted iterations. For the mixed-state packed traversal
(`TRAVERSE_PACKED_BLOCKS`), `insertBit0` is the natural primitive: it maps a dense counter
to sparse row/column indices in triangular storage where no regular stride pattern exists.

**Stride-based enumeration** achieves the same index sequence through nested loops with
power-of-two strides:

```
stride = 2^q
for (i = 0; i < dim; i += 2·stride)
    for (j = 0; j < stride; ++j)
        pair: (i+j, i+j+stride)
```

Both produce identical index pairs. The choice between them is not mathematical but
microarchitectural, and has significant performance consequences (Section 2).

**Gate-agnostic traversal.** The traversal macros (`TRAVERSE_PURE_1Q`, `TRAVERSE_PURE_2Q`)
accept a `PAIR_OP` parameter — the gate-specific operation applied to each element pair. The
same traversal serves every gate: for CX, PAIR_OP is a swap; for CZ, PAIR_OP negates the
target-bit-1 amplitude; for CY, PAIR_OP applies the Y rotation to the pair. Any
controlled-U gate reuses `TRAVERSE_PURE_2Q` with the appropriate single-qubit PAIR_OP. This
separation of traversal from arithmetic is what makes the per-gate optimizations in
Section 10 composable with the parallelization and vectorization strategies.

---

## 2. Auto-Vectorization: Why Stride-Based Loops Outperform Bit-Manipulation

The stride loop exposes the contiguous access pattern to the compiler: `data[i + j]` with
`j++` is trivially vectorizable. `insertBit0(k, target)` computes each address via bit
manipulation — the compiler cannot prove that consecutive `k` values yield consecutive
addresses, so it falls back to scalar code.

The stride-based inner loops have three properties that enable auto-vectorization:

1. **Stride-1 access**: the innermost loop increments `j` by 1
2. **No loop-carried dependencies**: each pair (a, b) is processed independently
3. **Simple arithmetic**: add, multiply, negate — operations with direct SIMD equivalents

These are exactly the conditions under which compilers emit vector instructions (SSE/AVX on
x86, NEON on ARM). When `insertBit0` is used instead, the compiler sees an opaque integer
function in the address calculation and cannot prove the sequential access pattern, even
though it exists.

**Verification.** Inspecting compiler output (`-S -O3 -march=native`) should confirm that
stride-based inner loops produce packed vector instructions (`vaddpd`/`vmulpd` on 256-bit
YMM registers with AVX2), while `insertBit0`-based loops produce scalar `addsd`/`mulsd`
on individual XMM registers. This is a concrete, reproducible demonstration that two
mathematically equivalent index computations produce different machine code — a result that
connects the algorithmic design directly to hardware performance.

---

## 3. Compiler Optimization and SIMD Interaction

Gate operations are memory-bandwidth bound at scale. Performance depends on the compiler's
ability to translate the algorithmic access patterns into efficient machine code.

**Optimization flags and their roles:**

| Flag | Effect on Gate Kernels |
|------|----------------------|
| `-O3` | Enables auto-vectorization, function inlining, loop unrolling |
| `-march=native` | Selects best available SIMD width (AVX2: 4 doubles, AVX-512: 8 doubles) |
| `-ffast-math` | Permits floating-point reassociation, enabling reductions to vectorize |
| `-funroll-loops` | Unrolls short inner loops for better instruction pipelining |

**Why `-ffast-math` matters for quantum simulation.** Strict IEEE 754 semantics prevent
reordering of floating-point additions, which blocks vectorization of reduction patterns
(e.g., summing contributions in a Hadamard butterfly `a' = (a+b)/√2`). Since gate
operations accumulate O(2^n) roundoff regardless, the additional reassociation error from
`-ffast-math` is negligible compared to the inherent numerical noise of the simulation.
The tradeoff is clear: measurably faster gate application for immeasurably different
numerical results.

**Row-major within tiles, column-major for BLAS.** Tiled storage uses row-major layout
within each tile for SIMD-friendly sequential access during gate kernels — the inner loop
sweeps consecutive memory addresses. However, tile data pointers must not be passed to
LAPACK/BLAS routines that expect column-major (Fortran) layout without transposition. The
reference implementation (Kronecker products via `zgemm`) uses separate full-format matrices
for this reason.

---

## 4. Loop Count Derivations and Load Balancing

### TRAVERSE_PURE_1Q: dim/2 iterations

The inner loop increments by 1 (stride-1 access), while the outer loop steps over pairs of
contiguous blocks. `collapse(2)` gives dim/2 total iterations regardless of target position:
- target = 0 (stride = 1): outer loop has dim/2 iterations, inner has 1
- target = n−1 (stride = dim/2): outer loop has 1 iteration, inner has dim/2
- All cases: dim/2 total units of work, evenly distributable by OpenMP

This target-position invariance is essential: it guarantees uniform parallel efficiency
regardless of which qubit the gate acts on. No qubit is "harder" than another.

### TRAVERSE_PURE_2Q: dim/4 iterations

The three loops enumerate exactly dim/4 base values with bits `lo` and `hi` both equal to 0.
The outer loop covers bits above `hi` (2^(n−hi−1) iterations), the middle loop covers bits
between `lo` and `hi` (2^(hi−lo−1) iterations), and the inner loop covers bits below `lo`
(2^lo iterations). Total: 2^(n−2) = dim/4.

### TRAVERSE_PURE_3Q: dim/8 iterations

Three qubit positions remove three bits from the index space. The four nested loops yield
2^(n−3) = dim/8 iterations.

### Connection to the complexity bound

These iteration counts confirm Proposition 3: an m-qubit gate on an n-qubit system requires
O(2^n / 2^m) pair operations, each involving a 2^m × 2^m block. For m = 1: dim/2 pairs of 2
elements. For m = 2: dim/4 groups of 4 elements. The total work scales as O(dim) regardless
of m (for fixed m), with the constant determined by the gate arithmetic (Section 10).

---

## 5. Memory Access Patterns and Bandwidth Bound

### Pure state: single-qubit gates

The innermost loop accesses `data[i+j]` and `data[i+j+stride]` with `j` incrementing by 1.
The first stream is perfectly sequential. The second stream is sequential at an offset of
`stride = 2^q` elements. For low-index qubits (small q), both streams fit in the same cache
line. For high-index qubits (large q), the streams are far apart but individually
sequential — hardware prefetchers handle both patterns efficiently.

### Pure state: two-qubit gates

- When `target = lo` (target < control): pair elements are separated by `stride_lo` — the
  close case. Both elements likely share a cache line for small `lo`.
- When `target = hi` (target > control): pair elements are separated by `stride_hi` — far
  apart, but each stream is individually contiguous.

The innermost loop is always stride-1, preserving vectorizability regardless of which qubit
is control and which is target.

### Bandwidth bound

At scale (n ≥ 12), single-qubit gate operations are limited by memory bandwidth, not
arithmetic throughput. A single-qubit gate touches every element exactly once: 16 bytes
read + 16 bytes written per `double complex` element, for a total of 32·2^n bytes of memory
traffic. The achievable gate rate is therefore bounded by

    gates/second ≤ bandwidth / (32 · 2^n)

This bound is **independent of gate type** — an X gate (zero arithmetic) and an Rx gate
(trig-based mixing) move the same number of bytes. The zero-flop gates (Section 10) are
therefore useful benchmarks: they measure the pure overhead of the traversal and memory
subsystem, providing a ceiling for all other gates at the same system size.

---

## 6. OpenMP Parallelization Strategies

**Thresholds.** Different thresholds reflect different work-per-gate scaling:
- **Pure states**: `OMP_THRESHOLD = 4096` → n ≥ 12 qubits (O(dim) work per gate)
- **Mixed states**: `OMP_THRESHOLD = 64` → n ≥ 6 qubits (O(dim²) work per gate)

Below these thresholds, the overhead of thread creation and synchronization exceeds the
parallel benefit.

**Pure state strategy (1Q):** `collapse(2)` flattens the two-level loop, giving dim/2 total
iterations evenly distributed regardless of target qubit position.

**Pure state strategy (2Q):** `collapse(2)` on the outer two of three loops. The stride-1
inner loop remains intact as a vectorizable unit within each parallel chunk. Total parallel
granularity: dim/(4·stride_lo) chunks of stride_lo contiguous iterations each.

**Pure state strategy (3Q):** `collapse(3)` on the outer three of four loops. Same principle
as 2Q.

**Mixed packed strategy:** `schedule(dynamic)` on the outer `bc` loop of
`TRAVERSE_PACKED_BLOCKS`. Dynamic scheduling balances load because the inner loop iteration
count varies with `bc` (triangular iteration space: br ranges from bc to half_dim).

**Mixed tiled strategy:** Parallelize over tiles (within-tile case) or tile groups
(cross-tile case) with dynamic scheduling.

**Concrete example (n = 8):** dim = 256, giving 36 lower-triangular tiles → 36 independent
parallel units (within-tile operations) or 9 tile groups (cross-tile operations).

**Thread safety within gates.** State objects are not thread-safe for concurrent mutation —
users must create separate `state_t` objects for parallel circuit evaluation. Within a
single gate application, thread safety is guaranteed by the butterfly decomposition: each
thread operates on a disjoint set of element pairs, determined by the loop partitioning.
No synchronization is needed beyond the implicit OpenMP barrier at loop exit.

---

## 7. Mixed State Packed Traversal: The BLOCK_OP Contract

For ρ → UρU†, the conjugation mixes 2×2 blocks of the density matrix indexed by the target
qubit. The traversal macro `TRAVERSE_PACKED_BLOCKS` uses `insertBit0` at the target position
to enumerate all (row, column) base pairs in the lower triangle.

**BLOCK_OP interface.** The macro invokes `BLOCK_OP` with:

- `data`: pointer to the state data array
- `idx00`, `idx01`, `idx10`, `idx11`: packed storage indices for the four block elements
  ρ(r0,c0), ρ(r0,c1), ρ(r1,c0), ρ(r1,c1)
- `lower_01`: 1 if element (r0, c1) is in the lower triangle of the full matrix; 0 if it
  must be accessed via its conjugate at position (c1, r0)
- `diag_block`: 1 if the block lies on the matrix diagonal (bc == br), where the (0,1) and
  (1,0) sub-block elements are conjugates of each other

The `lower_01` and `diag_block` flags encode the consequences of Hermitian symmetry in
packed storage. Every BLOCK_OP must handle both cases correctly — reading from the conjugate
position when `lower_01 = 0`, and avoiding double-counting on diagonal blocks. This
per-element branching is the principal computational overhead of the packed format compared
to full storage, and the fundamental motivation for the tiled format (Section 9).

---

## 8. Controlled Gate Conjugation: Four-Block Decomposition

A controlled-U gate V = |0⟩⟨0| ⊗ I + |1⟩⟨1| ⊗ U produces four density matrix blocks under
conjugation VρV†:

```
H'₀₀ = H₀₀                   (control=0 in both row and col: unchanged)
H'₀₁ = H₀₁ · U†              (control=0 in row, 1 in col: right-multiply by U†)
H'₁₀ = U · H₁₀               (control=1 in row, 0 in col: left-multiply by U)
H'₁₁ = U · H₁₁ · U†          (control=1 in both: full conjugation)
```

**Cross-block modification is required.** Only H₀₀ is left unchanged. The cross-blocks H₀₁
and H₁₀ undergo one-sided multiplication and must not be skipped — an implementation that
only applies conjugation to the H₁₁ block produces incorrect mixed state evolution.

### Left-multiply by U (H₁₀ block)

Each column of the 2×2 sub-block is transformed by U:

```
ρ'₀₀ = U₀₀·ρ₀₀ + U₀₁·ρ₁₀
ρ'₀₁ = U₀₀·ρ₀₁ + U₀₁·ρ₁₁
ρ'₁₀ = U₁₀·ρ₀₀ + U₁₁·ρ₁₀
ρ'₁₁ = U₁₀·ρ₀₁ + U₁₁·ρ₁₁
```

### Right-multiply by U† (H₀₁ block)

Each row of the 2×2 sub-block is transformed by U†:

```
ρ'₀₀ = ρ₀₀·conj(U₀₀) + ρ₀₁·conj(U₁₀)
ρ'₀₁ = ρ₀₀·conj(U₀₁) + ρ₀₁·conj(U₁₁)
ρ'₁₀ = ρ₁₀·conj(U₀₀) + ρ₁₁·conj(U₁₀)
ρ'₁₁ = ρ₁₀·conj(U₀₁) + ρ₁₁·conj(U₁₁)
```

These are the explicit component forms of the matrix products UH and HU† respectively,
applied block-wise to 2×2 sub-matrices of the density matrix. For specific gates, many terms
vanish or simplify: CX and CZ require zero multiplications in all four blocks.

---

## 9. Design Trade-off: One-Bit vs Two-Bit Insertion

The choice of traversal strategy for controlled gates on packed vs tiled storage illustrates
how **data layout determines algorithmic structure** — the same mathematical operation
(VρV†) leads to qualitatively different implementations depending on the storage format.

### Packed format: one-bit insertion preferred

One-bit insertion (insert 0-bit at the target position only) keeps control bits as natural
values in the loop indices and uses a 4-way branch to dispatch to the appropriate block
operation (H₀₀ skip, H₁₁ conjugation, H₁₀ left-multiply, H₀₁ right-multiply).

Two-bit insertion would eliminate the H₀₀ skip iterations (~25% of total) and remove the
4-way branch. However, in the packed format:

- The 4-way branch cost is ~1–3% of total per-iteration cost (2–4 cycles branch vs 50–200
  cycles memory access)
- The control bit flips slowly as `br` increments, giving the branch predictor long
  predictable runs (misprediction rate < 1%)
- Two-bit insertion would require three separate traversal passes (H₁₁, H₁₀, H₀₁), each
  with its own OpenMP barrier and its own triangular constraint logic
  (`lower_01`/`diag_block` must be recomputed per pass)
- Extending to multi-controlled gates (Toffoli) with two-bit insertion requires a
  combinatorial number of passes; one-bit insertion extends trivially via a mask check

Overhead of one-bit insertion: 4 integer operations per block (shift + AND for row and
column control bits). Negligible compared to memory access cost. Compensated by skipping
1/4 of all blocks (H₀₀).

### Tiled format: two-bit insertion becomes viable

Unlike the packed format, tiled storage uses full (non-triangular) storage within each tile.
This eliminates the `lower_01` and `diag_block` complications that make multi-pass traversal
expensive in packed format.

Two-bit insertion for tiled controlled gates: insert 0-bits at both target and control
positions, then enumerate the three non-trivial block types (H₁₁, H₁₀, H₀₁) in separate
branch-free passes. Each pass performs a uniform operation on every iteration, enabling
better instruction pipelining than the 4-way branch.

For within-tile operations (both qubits < LOG_TILE_DIM), all four butterfly elements lie
within the same tile with full row-major storage — no triangle constraints, no conjugation,
no aliasing. This is exactly the regime where two-bit insertion pays off: the branch-free
inner loop auto-vectorizes cleanly.

For cross-tile operations, the triangle constraint reappears at the tile level: some of the
four tiles in a group may lie in the upper triangle of the tile grid and must be accessed via
their conjugate. But this constraint is determined once per tile group (not per element),
making it far cheaper than the packed element-level branching.

### Summary

| Property | Packed (one-bit) | Tiled (two-bit) |
|----------|-----------------|-----------------|
| Passes | 1 | 3 |
| Branching | 4-way per block | None within pass |
| Triangle logic | Per element | Per tile group |
| Multi-control extension | Mask check | Combinatorial passes |
| Inner loop vectorization | Branch inhibits | Branch-free enables |

The packed format optimizes for memory (half the storage of full matrices) at the cost of
per-element branching. The tiled format trades ~2× memory for branch-free inner loops. The
choice between them is a memory–compute trade-off governed by the target system size: for
small n where the full state fits in cache, packed format's memory savings dominate; for
large n where bandwidth is saturated, the tiled format's vectorization advantage becomes
decisive.

---

## 10. Per-Gate Optimizations: Zero-Flop Gates and Operation Counts

Proposition 3 establishes O(r·N²·2^k) as the complexity bound for a k-local gate on an
n-qubit mixed state. The constant factor r depends on the gate's arithmetic structure.
Several common gates require **zero floating-point multiplications**, performing only memory
moves and sign flips — they are pure memory operations whose execution time is determined
entirely by the memory subsystem.

### Pure state

| Gate | Operation per pair | Multiplications | Dominant cost |
|------|-------------------|-----------------|--------------|
| X | swap(a, b) | 0 | Memory (2 loads + 2 stores) |
| Y | a' = −ib, b' = ia | 0 | Component swap + negate |
| Z | b' = −b | 0 | Single negation |
| S | b' = ib | 0 | Component swap (real↔imag) |
| Sdg | b' = −ib | 0 | Component swap |
| T | b' = e^{iπ/4}b | 2 | Via (re−im, re+im)/√2 |
| Tdg | b' = e^{−iπ/4}b | 2 | Via (re+im, im−re)/√2 |
| H | butterfly(a, b) | 2 | 2 adds + 2 muls (by 1/√2) |
| Rx(θ) | 2×2 unitary mixing | 4 | cos/sin precomputed |
| Ry(θ) | real rotation | 2 | Real coefficients only |
| Rz(θ) | diagonal phase | 2 | a unchanged, b rotated |

### Mixed state (packed)

| Gate | Block operation | Multiplications per block | Note |
|------|----------------|--------------------------|------|
| X | swap 4 elements | 0 | Pure memory operation |
| Y | swap + negate | 0 | Pure memory operation |
| Z | negate off-diagonal pair | 0 | Pure memory operation |
| S | multiply by ±i | 0 | Component swap (real↔imag) |
| Sdg | multiply by ∓i | 0 | Component swap |
| T | phase rotation | 4 | 2 muls per off-diagonal element |
| Tdg | phase rotation | 4 | 2 muls per off-diagonal element |
| H | butterfly | 6 | 4 adds + 2 muls per 2×2 block |
| Rx(θ) | 2×2 block mixing | 8 | Via component arithmetic |
| Ry(θ) | 2×2 block mixing | 6 | Real coefficients c², s², cs |
| Rz(θ) | off-diagonal phases | 4 | Diagonal unchanged |

The zero-flop gates (X, Y, Z, S, Sdg, and their controlled variants CX, CZ) serve as
**bandwidth benchmarks**: they isolate the cost of traversal and indexing from arithmetic,
establishing the floor for all gate execution times at a given system size. Any gap between
a non-trivial gate (e.g., H) and the zero-flop floor measures the true cost of the gate's
arithmetic, independent of memory system overhead.

---

## 11. Tiled Format: Conjugation-Free Operations and Cache Performance

### Within-tile conjugation proof

**No conjugation logic is needed for within-tile operations** (target qubit < LOG_TILE_DIM):

- **Diagonal tiles** (tr == tc): store the full TILE_DIM × TILE_DIM region including both
  lower and upper triangle portions of the density matrix. All 4 butterfly elements are
  directly accessible at their row-major offsets.
- **Off-diagonal tiles** (tr > tc): for any butterfly anchored at (r0, c0) within this tile,
  r0 ≥ tr · TILE_DIM ≥ (tc + 1) · TILE_DIM > c0 + stride ≥ c1, so all 4 elements satisfy
  r ≥ c and lie in the lower triangle. No conjugation is needed.

This is the theoretical justification for the tiled format's existence: it converts a
per-element runtime decision (which triangle?) into a per-tile structural guarantee (always
lower triangle within tile). The consequence is that within-tile gate kernels have no
data-dependent branches — all elements are accessed at compile-time-predictable offsets,
enabling auto-vectorization.

### Cache performance

| Case | Working Set | Cache Level | Conjugation Overhead |
|------|-------------|-------------|---------------------|
| Within-tile (target < LOG_TILE_DIM) | 16 KB (1 tile) | L1 | None |
| Cross-tile (target ≥ LOG_TILE_DIM) | 64 KB (4 tiles) | L2 | Per tile group |

With TILE_DIM = 32, each tile occupies 32 × 32 × 16 bytes = 16,384 bytes. Typical L1 data
caches are 32–64 KB, so a single tile fits entirely in L1 with room for other working data.
Operations on the rightmost LOG_TILE_DIM = 5 qubits (qubits 0–4) stay within a single tile
and execute entirely from L1.

Cross-tile operations load 4 tiles (64 KB), which fits in L2 (typically 256 KB–1 MB). The
tile size is therefore chosen to match the two natural qubit index ranges to the two
innermost cache levels: low-index qubits → L1 resident, high-index qubits → L2 resident.

### Within-tile kernel

```c
static inline void apply_butterfly_within_tile(
    cplx_t *tile, dim_t lr0, dim_t lc0, dim_t stride
) {
    dim_t lr1 = lr0 + stride, lc1 = lc0 + stride;
    cplx_t rho00 = tile[lr0 * TILE_DIM + lc0];
    cplx_t rho01 = tile[lr0 * TILE_DIM + lc1];
    cplx_t rho10 = tile[lr1 * TILE_DIM + lc0];
    cplx_t rho11 = tile[lr1 * TILE_DIM + lc1];
    // Apply gate to (rho00, rho01, rho10, rho11), write back directly
    // No lower_01 or diag_block flags — all four elements are valid
}
```

---

## 12. SWAP on Packed Storage: A Case Study

SWAP is a permutation gate: perm(x) = x with bits q1 and q2 exchanged. Like CX, it is
zero-flop (only memory moves). Unlike CX, SWAP is not a controlled-U gate and does not
decompose into the four blocks of Section 8.

**The packed format complicates in-place SWAP.** For a 4×4 block of the density matrix
(indexed by the two qubit values in row and column), SWAP permutes rows and columns by
exchanging the q1↔q2 bit positions. This produces 6 element swaps per block. On diagonal
base blocks (r_base == c_base), Hermitian symmetry causes some swaps to alias the same
packed location, reducing to 3 swaps + 1 in-place conjugation. On off-diagonal base blocks,
the lower-triangle constraint makes the triangle status of each element data-dependent,
introducing up to 12 branches per block and aliasing hazards where different swaps map to
the same packed element through conjugate symmetry.

This is an instance of a general problem: **permutation gates on symmetry-reduced storage
require alias analysis that does not arise with full storage.** The packed format's memory
savings come at the cost of making even zero-arithmetic operations algorithmically
non-trivial.

**Three implementation strategies:**

1. **CX decomposition (recommended):** SWAP = CX(a,b) · CX(b,a) · CX(a,b). Reuses the
   in-place four-block CX implementation three times. O(1) extra memory, correct by
   construction. Cost: 3× the memory traffic of a direct implementation. Since SWAP is
   memory-bound, this means ~3× wall-clock time — but with zero implementation complexity.

2. **Output-buffer approach:** Allocate O(dim²) temporary, iterate over all lower-triangle
   elements computing perm(r) and perm(c), gather from source, write to destination. Simple
   (~12 lines of core logic), no aliasing, embarrassingly parallel. Cost: O(dim²) temporary
   allocation, which may be prohibitive at scale.

3. **In-place two-bit insertion:** Enumerate 4×4 blocks via `insertBits2_0`, perform swaps
   with conditional conjugation. Requires loading all 16 block values before writing to
   avoid aliasing. Complex (~70 lines), high bug risk. Not recommended.

The CX decomposition is preferred: it trades a 3× constant-factor slowdown for correctness
by construction and zero additional code. For a gate that appears infrequently in typical
circuits (SWAP is often decomposed at the circuit level), this trade-off is clearly
favorable.
