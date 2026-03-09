# Benchmark Methodology Note

**Date:** 2026-03-09
**Concerns:** Apparent CX/SWAP performance regression between `bench.dat` (~2026-02-23)
and `bench_2_to_8_qubits.dat` / `bench_9_to_12_qubits.dat` (~2026-03-08)

---

## Finding: No performance regression

The 4–6× increase in reported CX and SWAP times is not a gate implementation
regression. It is a benchmark methodology change introduced in commit `20d88a7`
("Refactor benchmarks, strip BLAS from library types…", 2026-03-07).

### What changed

| | Old (`bench.dat`) | New (split `.dat` files) |
|---|---|---|
| 2Q sweep | Adjacent pairs `(t, t+1)`, **n−1** pairs | All ordered pairs `(q1 < q2)`, **n(n−1)/2** pairs |
| Reported `time_ms` | Total for `iterations × (n−1)` applications | Total for `iterations × n(n−1)/2` applications |

The ratio of work done is `[n(n−1)/2] / (n−1) = n/2`, which exactly matches the
observed slowdown at every qubit count:

| n | Predicted (n/2) | qlib CX | QuEST CX | Qulacs CX |
|---|---|---|---|---|
| 2 | 1.0× | 0.96× | 1.00× | 0.99× |
| 5 | 2.5× | 2.57× | 2.59× | 2.67× |
| 8 | 4.0× | 4.17× | 4.01× | 4.40× |
| 10 | 5.0× | 5.23× | 5.01× | 5.03× |
| 12 | 6.0× | 6.32× | 6.33× | 6.23× |

Because all three libraries (qlib, QuEST, Qulacs) degrade by the same factor,
the gate implementations themselves are unchanged.

The underlying per-gate cost is identical. For example, qlib_packed CX at n=5:
- Old: 0.795 ms / 4 adjacent pairs ≈ **0.199 ms/pair**
- New: 2.046 ms / 10 all-pairs ≈ **0.205 ms/pair**

---

## Design issue: `time_ms` is aggregate, not per-gate

Currently `time_ms` in the `.dat` output is the total wall time for a full
sweep times `iterations`. Concretely:

- 1Q gates: `time_ms` covers `iterations × n` gate applications (all qubit targets)
- 2Q gates: `time_ms` covers `iterations × n(n−1)/2` gate applications (all pairs)

This means columns are not comparable across gate arity or across runs with
different sweep sizes. At n=8, the X column covers 8 applications while the CX
column covers 28. A reader of the `.dat` files — or a paper — would have no way
to reconstruct the per-gate cost without digging into the source.

The `ops_per_sec` field already encodes the correct normalized value, but it is
not written to the `.dat` output.

---

## Recommendation

**Report time per gate application in all output formats.**

The all-pairs sweep is the right *methodology* — it exercises within-tile,
cross-tile, and tile-boundary code paths and reduces variance. But the *reported
metric* should be normalized:

```
time_per_gate_ms = time_ms / (iterations × sweep_size)
```

where `sweep_size = n` for 1Q gates and `n(n−1)/2` for 2Q gates. Equivalently,
this is `1000 / ops_per_sec`.

### Concrete changes

1. **Add a `time_per_gate_ms` column** to the pgfplots `.dat` output (or replace
   the current `time_ms` mean/std columns with normalized versions). This is the
   number to plot and to cite in a paper.

2. **Normalize std dev the same way** so error bars remain valid:
   `std_per_gate_ms = time_ms_std / (iterations × sweep_size)`.

3. **Keep `ops_per_sec`** — it is already correct and useful as a secondary metric.

4. **Keep the all-pairs sweep** — do not revert to adjacent-only. The broader
   sweep is strictly better for exercising the implementation; only the output
   needs to change.

### Why this matters for the paper

- A reader seeing "CX at n=8 takes X ms" must be able to interpret that as one
  gate application, not an implementation-defined sweep.
- Normalization makes results reproducible: another group benchmarking with a
  different number of pairs or iterations will report the same per-gate number.
- Cross-gate comparison (X vs CX at the same qubit count) becomes valid.
