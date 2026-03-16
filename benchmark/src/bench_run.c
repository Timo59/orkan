/**
 * @file bench_run.c
 * @brief Two-level timing loop (samples x iterations) and position generation.
 */

#include "bench.h"

/* =====================================================================
 * Position generation
 * ===================================================================== */

int bench_gen_positions(qubit_t n_qubits, int gate_qubits, qubit_t *out) {
    int count = 0;
    int n = (int)n_qubits;

    if (gate_qubits == 1) {
        for (int t = 0; t < n && count < BENCH_MAX_POSITIONS; ++t)
            out[count++] = (qubit_t)t;
    } else if (gate_qubits == 2) {
        for (int c = 0; c < n; ++c)
            for (int t = 0; t < n; ++t)
                if (c != t && count < BENCH_MAX_POSITIONS) {
                    out[count * 2]     = (qubit_t)c;
                    out[count * 2 + 1] = (qubit_t)t;
                    count++;
                }
    } else if (gate_qubits == 3) {
        for (int c1 = 0; c1 < n; ++c1)
            for (int c2 = 0; c2 < n; ++c2)
                for (int t = 0; t < n; ++t)
                    if (c1 != c2 && c1 != t && c2 != t
                        && count < BENCH_MAX_POSITIONS) {
                        out[count * 3]     = (qubit_t)c1;
                        out[count * 3 + 1] = (qubit_t)c2;
                        out[count * 3 + 2] = (qubit_t)t;
                        count++;
                    }
    }

    if (gate_qubits == 3 && count == BENCH_MAX_POSITIONS) {
        int full = n * (n - 1) * (n - 2);
        if (full > BENCH_MAX_POSITIONS)
            fprintf(stderr, "Warning: 3Q positions truncated to %d of %d "
                    "for %d qubits\n", BENCH_MAX_POSITIONS, full, n);
    }
    return count;
}

/* =====================================================================
 * Aggregate mode
 * ===================================================================== */

void bench_run_aggregate(const bench_backend_t *be, void *ctx,
                         qubit_t qubits,
                         const qubit_t *positions, int n_pos, int gate_qubits,
                         const bench_options_t *opts,
                         bench_result_t *out) {
    /* Warm-up: execute full sweeps (not timed) */
    for (int w = 0; w < opts->warm_up; ++w)
        for (int p = 0; p < n_pos; ++p)
            be->apply(ctx, &positions[p * gate_qubits]);

    /* Two-level timing: samples x iterations */
    double samples[BENCH_MAX_SAMPLES];
    int n_samples = opts->samples;
    if (n_samples > BENCH_MAX_SAMPLES) n_samples = BENCH_MAX_SAMPLES;

    for (int s = 0; s < n_samples; ++s) {
        bench_timing_barrier();
        uint64_t t0 = bench_time_ns();
        bench_timing_barrier();

        for (int i = 0; i < opts->iterations; ++i)
            for (int p = 0; p < n_pos; ++p)
                be->apply(ctx, &positions[p * gate_qubits]);

        bench_timing_barrier();
        uint64_t t1 = bench_time_ns();

        double total_ns = (double)(t1 - t0);
        samples[s] = total_ns / ((double)opts->iterations * n_pos);
    }

    *out = bench_compute_result(samples, n_samples, be->touch_bytes(qubits));
}

/* =====================================================================
 * Per-qubit mode
 * ===================================================================== */

void bench_run_per_qubit(const bench_backend_t *be, void *ctx,
                         qubit_t qubits,
                         const qubit_t *positions, int n_pos, int gate_qubits,
                         const bench_options_t *opts,
                         bench_pos_result_t *out) {
    double samples[BENCH_MAX_SAMPLES];
    int n_samples = opts->samples;
    if (n_samples > BENCH_MAX_SAMPLES) n_samples = BENCH_MAX_SAMPLES;
    size_t touch = be->touch_bytes(qubits);

    for (int p = 0; p < n_pos; ++p) {
        const qubit_t *pos = &positions[p * gate_qubits];

        /* Warm-up for this position */
        for (int w = 0; w < opts->warm_up; ++w)
            be->apply(ctx, pos);

        /* Two-level timing */
        for (int s = 0; s < n_samples; ++s) {
            bench_timing_barrier();
            uint64_t t0 = bench_time_ns();
            bench_timing_barrier();

            for (int i = 0; i < opts->iterations; ++i)
                be->apply(ctx, pos);

            bench_timing_barrier();
            uint64_t t1 = bench_time_ns();

            samples[s] = (double)(t1 - t0) / opts->iterations;
        }

        out[p].result = bench_compute_result(samples, n_samples, touch);
        for (int q = 0; q < gate_qubits; ++q)
            out[p].pos[q] = pos[q];
        for (int q = gate_qubits; q < 3; ++q)
            out[p].pos[q] = 0;
    }
}
