/**
 * @file bench_circuit_qulacs.cpp
 * @brief Qulacs density matrix circuit-level benchmark runner.
 *
 * Implements bench_circuit_qulacs() and the context alloc/free pair for
 * calibrate_iterations compatibility. Uses Qulacs csim low-level DM gate
 * functions (dm_H_gate, dm_X_gate, dm_RX_gate, dm_CNOT_gate, etc.) directly,
 * bypassing the higher-level cppsim gate object layer for minimal overhead.
 *
 * State is re-initialized to the maximally mixed state (rho = I/dim) before
 * each of the `runs` timed runs. Within a single run, state evolves across
 * `iterations` full circuit applications -- this measures sustained throughput
 * rather than single-shot latency.
 *
 * Guarded by WITH_QULACS; compiled as C++ with extern "C" linkage for the
 * three exported symbols.
 */

#ifdef WITH_QULACS

#include <csim/update_ops_dm.hpp>
#include <csim/type.hpp>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cstdio>

#include "bench_circuit.h"

/*
 * =====================================================================================================================
 * Gate dispatch
 * =====================================================================================================================
 */

/**
 * @brief Apply a single circuit operation to the Qulacs density matrix.
 *
 * Dispatches to the appropriate Qulacs csim dm_*_gate function.
 * Rotation gate signatures: (UINT target, double angle, CTYPE *state, ITYPE dim).
 * CNOT signature: (UINT control, UINT target, CTYPE *state, ITYPE dim).
 */
static void apply_circuit_op_qulacs(const circuit_op_t *op, CTYPE *rho, ITYPE dim) {
    switch (op->gate_type) {
        case GATE_H:   dm_H_gate(static_cast<UINT>(op->q0), rho, dim); break;
        case GATE_X:   dm_X_gate(static_cast<UINT>(op->q0), rho, dim); break;
        case GATE_Z:   dm_Z_gate(static_cast<UINT>(op->q0), rho, dim); break;
        case GATE_RX:  dm_RX_gate(static_cast<UINT>(op->q0), op->angle, rho, dim); break;
        case GATE_RY:  dm_RY_gate(static_cast<UINT>(op->q0), op->angle, rho, dim); break;
        case GATE_RZ:  dm_RZ_gate(static_cast<UINT>(op->q0), op->angle, rho, dim); break;
        case GATE_CX:  dm_CNOT_gate(static_cast<UINT>(op->q0), static_cast<UINT>(op->q1),
                                     rho, dim); break;
        default:
            fprintf(stderr, "bench_circuit_qulacs: unimplemented gate type %d\n",
                    op->gate_type);
    }
}

/*
 * =====================================================================================================================
 * Maximally mixed state initialization
 * =====================================================================================================================
 */

/**
 * @brief Set rho = I/dim: zero all elements, then set diagonal to 1/dim.
 */
static void init_maximally_mixed_qulacs(CTYPE *rho, ITYPE dim) {
    size_t n_elems = static_cast<size_t>(dim) * static_cast<size_t>(dim);
    std::memset(rho, 0, n_elems * sizeof(CTYPE));
    double diag_val = 1.0 / static_cast<double>(dim);
    for (ITYPE i = 0; i < dim; ++i) {
        rho[i * dim + i] = CTYPE(diag_val, 0.0);
    }
}

/*
 * =====================================================================================================================
 * Exported functions (extern "C")
 * =====================================================================================================================
 */

extern "C" {

/**
 * @brief Run circuit benchmark using Qulacs density matrix.
 *
 * Allocates a dim x dim CTYPE density matrix, re-initializes to maximally
 * mixed state before each timed run, and applies the full circuit for
 * `iterations` repetitions per run. Reports statistics across `runs` runs.
 */
bench_result_t bench_circuit_qulacs(const circuit_t *c, int iterations,
                                     int warmup, int runs) {
    (void)warmup; /* warmup is time-based, not count-based */

    bench_result_t result;
    std::memset(&result, 0, sizeof(result));
    result.qubits     = static_cast<qubit_t>(c->n_qubits);
    result.dim        = static_cast<dim_t>(1) << c->n_qubits;
    result.gate_name  = c->name;
    result.method     = "qulacs";
    result.iterations = iterations;
    result.sweep_size = c->n_ops;
    result.runs       = runs;

    ITYPE dim = static_cast<ITYPE>(1) << c->n_qubits;
    size_t matrix_size = static_cast<size_t>(dim) * static_cast<size_t>(dim) * sizeof(CTYPE);
    CTYPE *rho = static_cast<CTYPE *>(std::malloc(matrix_size));
    if (!rho) {
        return result;
    }
    result.memory_bytes = matrix_size;

    /* Time-based warmup: run full circuits until BENCH_PQ_WARMUP_MS elapsed */
    init_maximally_mixed_qulacs(rho, dim);
    uint64_t wt0 = bench_time_ns();
    do {
        for (int op = 0; op < c->n_ops; ++op) {
            apply_circuit_op_qulacs(&c->ops[op], rho, dim);
        }
    } while (bench_ns_to_ms(bench_time_ns() - wt0) < BENCH_PQ_WARMUP_MS);

    /* Timed runs */
    int cnt = (runs < BENCH_MAX_RUNS) ? runs : BENCH_MAX_RUNS;
    double run_times[BENCH_MAX_RUNS];

    for (int r = 0; r < cnt; ++r) {
        /* Re-initialize to maximally mixed state before each timed run */
        init_maximally_mixed_qulacs(rho, dim);

        bench_timing_barrier();
        uint64_t t0 = bench_time_ns();
        bench_timing_barrier();

        for (int i = 0; i < iterations; ++i) {
            for (int op = 0; op < c->n_ops; ++op) {
                apply_circuit_op_qulacs(&c->ops[op], rho, dim);
            }
        }

        bench_timing_barrier();
        run_times[r] = bench_ns_to_ms(bench_time_ns() - t0);
    }

    /* Compute statistics */
    bench_run_stats_t stats = bench_compute_stats(run_times, cnt);
    result.time_ms        = stats.mean;
    result.time_ms_std    = stats.std_dev;
    result.time_ms_min    = stats.min;
    result.time_ms_median = stats.median;
    result.time_ms_cv     = stats.cv;

    double total_ops = static_cast<double>(iterations) * c->n_ops;
    result.ops_per_sec     = (stats.mean > 0.0) ? total_ops / (stats.mean / 1000.0) : 0.0;
    result.time_per_gate_ms = (total_ops > 0.0) ? stats.mean / total_ops : 0.0;

    std::free(rho);
    return result;
}

/**
 * @brief Allocate a Qulacs circuit context for calibrate_iterations.
 *
 * Allocates the density matrix and initializes to maximally mixed state.
 */
qulacs_circuit_ctx_t *qulacs_circuit_ctx_alloc(int n_qubits) {
    qulacs_circuit_ctx_t *ctx = static_cast<qulacs_circuit_ctx_t *>(
        std::malloc(sizeof(qulacs_circuit_ctx_t)));
    if (!ctx) return nullptr;

    ITYPE dim = static_cast<ITYPE>(1) << n_qubits;
    size_t matrix_size = static_cast<size_t>(dim) * static_cast<size_t>(dim) * sizeof(CTYPE);
    CTYPE *rho = static_cast<CTYPE *>(std::malloc(matrix_size));
    if (!rho) {
        std::free(ctx);
        return nullptr;
    }

    init_maximally_mixed_qulacs(rho, dim);

    ctx->rho = rho;
    ctx->dim = static_cast<size_t>(dim);
    return ctx;
}

/**
 * @brief Free a Qulacs circuit context.
 */
void qulacs_circuit_ctx_free(qulacs_circuit_ctx_t *ctx) {
    if (!ctx) return;
    std::free(ctx->rho);
    std::free(ctx);
}

} /* extern "C" */

#endif /* WITH_QULACS */
