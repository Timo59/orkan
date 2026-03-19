/**
 * @file bench_kraus_run.c
 * @brief Two-level timing loop and position generation for the Kraus benchmark.
 */

#include "bench_kraus.h"

/* =====================================================================
 * Position generation (same combinatorial logic as gate benchmark)
 * ===================================================================== */

int kraus_gen_positions(qubit_t n_qubits, int tgt_qubits, qubit_t *out) {
    int count = 0;
    int n = (int)n_qubits;

    if (tgt_qubits == 1) {
        for (int t = 0; t < n && count < BENCH_MAX_POSITIONS; ++t)
            out[count++] = (qubit_t)t;
    } else if (tgt_qubits == 2) {
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b)
                if (a != b && count < BENCH_MAX_POSITIONS) {
                    out[count * 2]     = (qubit_t)a;
                    out[count * 2 + 1] = (qubit_t)b;
                    count++;
                }
    } else if (tgt_qubits == 3) {
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b)
                for (int c = 0; c < n; ++c)
                    if (a != b && a != c && b != c
                        && count < BENCH_MAX_POSITIONS) {
                        out[count * 3]     = (qubit_t)a;
                        out[count * 3 + 1] = (qubit_t)b;
                        out[count * 3 + 2] = (qubit_t)c;
                        count++;
                    }
    }
    return count;
}

/* =====================================================================
 * Aggregate mode
 * ===================================================================== */

void kraus_run_aggregate(const bench_kraus_backend_t *be, void *ctx,
                         qubit_t qubits,
                         const qubit_t *positions, int n_pos, int tgt_qubits,
                         const bench_options_t *opts,
                         bench_result_t *out) {
    /* Warm-up */
    for (int w = 0; w < opts->warm_up; ++w)
        for (int p = 0; p < n_pos; ++p)
            be->apply(ctx, &positions[p * tgt_qubits]);

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
                be->apply(ctx, &positions[p * tgt_qubits]);

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

void kraus_run_per_qubit(const bench_kraus_backend_t *be, void *ctx,
                         qubit_t qubits,
                         const qubit_t *positions, int n_pos, int tgt_qubits,
                         const bench_options_t *opts,
                         bench_pos_result_t *out) {
    double samples[BENCH_MAX_SAMPLES];
    int n_samples = opts->samples;
    if (n_samples > BENCH_MAX_SAMPLES) n_samples = BENCH_MAX_SAMPLES;
    size_t touch = be->touch_bytes(qubits);

    for (int p = 0; p < n_pos; ++p) {
        const qubit_t *pos = &positions[p * tgt_qubits];

        for (int w = 0; w < opts->warm_up; ++w)
            be->apply(ctx, pos);

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
        for (int q = 0; q < tgt_qubits; ++q)
            out[p].pos[q] = pos[q];
        for (int q = tgt_qubits; q < 3; ++q)
            out[p].pos[q] = 0;
    }
}
