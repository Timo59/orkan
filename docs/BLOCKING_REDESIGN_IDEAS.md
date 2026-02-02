# Blocking Redesign Ideas

Captured from design discussion, January 2025.

## Problem Statement

The current blocked implementation (`mhipster_block.c`) uses 64×64 tiles for cache locality, but the parallelization strategy doesn't align well with quantum gate operations.

### Current Issues

1. **Tiles don't capture quantum structure**: A 64×64 tile captures indices where the lower 6 bits can vary freely. This has no relationship to Hamming distance between basis states or the structure of single-qubit gates.

2. **Cross-tile gates are awkward**: For gates on qubit k >= 6, the 4 elements of each butterfly span 4 different tiles. The current code iterates all tiles and skips half, rather than directly enumerating tile groups.

3. **Parallelism at wrong level**: We parallelize over tiles, but the natural unit of work for a single-qubit gate is the **butterfly** (a 2×2 transformation of 4 density matrix elements).

## Key Insight

For a single-qubit gate U on target qubit k, each element ρ[i,j] participates in exactly one butterfly:

```
Butterfly defined by (i₀, j₀) where both have bit k = 0:
  [ρ[i₀, j₀],       ρ[i₀, j₀⊕2^k]      ]
  [ρ[i₀⊕2^k, j₀],   ρ[i₀⊕2^k, j₀⊕2^k] ]
```

There are `dim²/4` butterflies, reduced to ~`dim²/8` by Hermitian symmetry.

**The natural parallelization unit is the butterfly, not the tile.**

## Proposed Redesign

### Option A: Gate-Centric Parallelism

Instead of:
```c
parallel_for tile in tiles:
    for butterfly in tile:
        process(butterfly)
```

Do:
```c
parallel_for butterfly_group in butterfly_groups:
    prefetch relevant tiles/cache lines
    for butterfly in group:
        process(butterfly)
```

Where `butterfly_group` is chosen to maximize cache reuse.

### Option B: Tile-Group Enumeration for Cross-Tile

For gates where `k >= log₂(BLOCK_DIM)`:

```
tile_bit = k - log₂(BLOCK_DIM)
tile_incr = 1 << tile_bit
```

Tiles partition into groups of 4:
```
Group (base_row, base_col) where both have tile_bit = 0:
  - (base_row,             base_col)              // (0,0)
  - (base_row,             base_col | tile_incr)  // (0,1)
  - (base_row | tile_incr, base_col)              // (1,0)
  - (base_row | tile_incr, base_col | tile_incr)  // (1,1)
```

Number of groups = `(n_blocks / 2)² / 2` for lower triangle.

Directly enumerate and process these groups rather than iterating all tiles and skipping.

### Option C: Adaptive Block Size

Choose block size based on target qubit:
- For small k: Use small blocks (fit more in cache during butterfly)
- For large k: Use larger blocks (amortize cross-tile overhead)

Or: hierarchical blocking with multiple levels.

## Questions to Resolve

1. What's the actual memory bandwidth vs. compute ratio for these operations?
2. Does SIMD within butterflies matter more than parallelism across butterflies?
3. Should we optimize for specific qubit counts (e.g., 10-14 qubits)?
4. Is the Hermitian storage worth the complexity, or should we store full matrices for small n?

## Benchmark Data Needed

- Profile cross-tile H gate to identify bottleneck (memory or compute)
- Compare butterfly-centric vs tile-centric parallelism
- Measure cache miss rates for different access patterns
