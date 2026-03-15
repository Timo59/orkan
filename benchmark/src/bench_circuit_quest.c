/**
 * @file bench_circuit_quest.c
 * @brief QuEST density matrix circuit-level benchmark runner.
 *
 * Implements bench_circuit_quest() and the associated context alloc/free
 * functions for circuit-level benchmarks using QuEST v4's density matrix
 * backend.
 *
 * QuEST v4 API differences from v3 (documented for maintainability):
 *   - Rotation gates renamed: rotateX -> applyRotateX, rotateY -> applyRotateY,
 *     rotateZ -> applyRotateZ (v3 names exist as deprecated macros).
 *   - createDensityQureg() takes only numQubits (no QuESTEnv argument).
 *   - destroyQureg() takes only Qureg (no QuESTEnv argument).
 *   - No initMixedState() for maximally mixed rho = I/dim. We implement this
 *     manually via initBlankState() + direct cpuAmps diagonal writes.
 *
 * Requires bench_quest_init() to have been called before any use.
 */

#ifdef WITH_QUEST

/* quest.h must precede bench_circuit.h: the latter includes bench.h whose
 * WITH_QUEST block declares types using Qureg. */
#include "quest/include/quest.h"
#include "bench_circuit.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * =====================================================================================================================
 * Maximally mixed state initialisation
 * =====================================================================================================================
 */

/**
 * @brief Initialise a QuEST density matrix to the maximally mixed state rho = I/dim.
 *
 * QuEST v4 has no built-in function for this. We zero the state via
 * initBlankState(), then set diagonal elements cpuAmps[i * dim + i] = 1/dim.
 * The Qureg.cpuAmps array stores the density matrix in column-major flat layout:
 * cpuAmps[col * dim + row] = rho(row, col), with dim = 2^numQubits.
 * Diagonal indices cpuAmps[i * dim + i] are invariant under transposition.
 */
static void init_maximally_mixed(Qureg rho) {
    qindex dim = (qindex)1 << rho.numQubits;
    qreal inv_dim = 1.0 / (qreal)dim;

    initBlankState(rho);

    for (qindex i = 0; i < dim; ++i) {
        rho.cpuAmps[i * dim + i] = inv_dim;
    }
}

/*
 * =====================================================================================================================
 * Gate dispatch
 * =====================================================================================================================
 */

/**
 * @brief Apply a single circuit operation to a QuEST density matrix.
 *
 * Gate name mapping (circuit_gate_t -> QuEST v4 function):
 *   CIRCUIT_GATE_H  -> applyHadamard(qureg, target)
 *   GATE_X  -> applyPauliX(qureg, target)
 *   GATE_Z  -> applyPauliZ(qureg, target)
 *   GATE_RX -> applyRotateX(qureg, target, angle)
 *   GATE_RY -> applyRotateY(qureg, target, angle)
 *   GATE_RZ -> applyRotateZ(qureg, target, angle)
 *   GATE_CX -> applyControlledPauliX(qureg, control, target)
 */
static void apply_circuit_op_quest(Qureg rho, const circuit_op_t *op) {
    switch (op->gate_type) {
        case CIRCUIT_GATE_H:   applyHadamard(rho, (int)op->q0);                         break;
        case GATE_X:   applyPauliX(rho, (int)op->q0);                           break;
        case GATE_Z:   applyPauliZ(rho, (int)op->q0);                           break;
        case GATE_RX:  applyRotateX(rho, (int)op->q0, op->angle);               break;
        case GATE_RY:  applyRotateY(rho, (int)op->q0, op->angle);               break;
        case GATE_RZ:  applyRotateZ(rho, (int)op->q0, op->angle);               break;
        case GATE_CX:  applyControlledPauliX(rho, (int)op->q0, (int)op->q1);    break;
        default:
            fprintf(stderr, "bench_circuit_quest: unimplemented gate type %d\n",
                    op->gate_type);
            break;
    }
}

/*
 * =====================================================================================================================
 * Circuit runner helper
 * =====================================================================================================================
 */

/**
 * @brief Execute the full circuit once on a QuEST density matrix.
 */
static void run_circuit_quest(Qureg rho, const circuit_t *c) {
    for (int i = 0; i < c->n_ops; ++i) {
        apply_circuit_op_quest(rho, &c->ops[i]);
    }
}

/*
 * =====================================================================================================================
 * QuEST environment lifecycle
 * =====================================================================================================================
 */

/**
 * @brief Initialise the QuEST environment.  Must be called once before any Qureg is created.
 *
 * Defined here (rather than in bench_quest.c) so that bench_circuit does not need to link
 * the full gate-level benchmark infrastructure from bench_quest.c.
 */
void bench_quest_init(void) {
    initQuESTEnv();
}

/**
 * @brief Finalise the QuEST environment.  Must be called once at program exit after all
 * Quregs have been destroyed.
 */
void bench_quest_cleanup(void) {
    finalizeQuESTEnv();
}

/*
 * =====================================================================================================================
 * Context alloc / free (for calibrate_iterations compatibility)
 * =====================================================================================================================
 */

quest_circuit_ctx_t *quest_circuit_ctx_alloc(int n_qubits) {
    quest_circuit_ctx_t *ctx = malloc(sizeof(*ctx));
    if (!ctx) return NULL;

    /* Heap-allocate a Qureg so we can store it as void* in the opaque context.
     * QuEST returns Qureg by value; we copy it into our allocation. */
    Qureg *qr = malloc(sizeof(Qureg));
    if (!qr) {
        free(ctx);
        return NULL;
    }
    *qr = createDensityQureg(n_qubits);
    initDebugState(*qr);  /* Touch all pages to avoid lazy-allocation artefacts */

    ctx->qureg = qr;
    return ctx;
}

void quest_circuit_ctx_free(quest_circuit_ctx_t *ctx) {
    if (!ctx) return;
    if (ctx->qureg) {
        Qureg *qr = (Qureg *)ctx->qureg;
        destroyQureg(*qr);
        free(qr);
    }
    free(ctx);
}

/*
 * =====================================================================================================================
 * Main benchmark entry point
 * =====================================================================================================================
 */

/**
 * @brief Run a circuit-level benchmark using QuEST density matrix backend.
 *
 * Structure:
 *   1. Create density matrix Qureg once.
 *   2. Time-based warmup: repeatedly execute the circuit until BENCH_PQ_WARMUP_MS elapses.
 *   3. For each of `runs` timed rounds:
 *      a. Re-initialise to maximally mixed state rho = I/dim.
 *      b. Barrier -> start timer -> barrier.
 *      c. Execute the circuit `iterations` times.
 *      d. Barrier -> stop timer.
 *   4. Compute statistics and return.
 *
 * Re-initialisation happens before each run (not between iterations within a run),
 * matching the contract specified in the plan.
 */
bench_result_t bench_circuit_quest(const circuit_t *c, int iterations,
                                   int warmup, int runs)
{
    bench_result_t result = {0};
    result.qubits     = (qubit_t)c->n_qubits;
    result.dim        = (dim_t)1 << c->n_qubits;
    result.gate_name  = c->name;
    result.method     = "quest";
    result.iterations = iterations;
    result.runs       = runs;
    result.sweep_size = c->n_ops;

    dim_t dim = (dim_t)1 << c->n_qubits;
    result.memory_bytes = (size_t)dim * (size_t)dim * sizeof(double) * 2; /* complex<double> */

    Qureg rho = createDensityQureg(c->n_qubits);
    initDebugState(rho);  /* Touch all pages */

    /* --- Time-based warmup ------------------------------------------------ */
    (void)warmup;  /* Use time-based warmup instead of count-based */
    init_maximally_mixed(rho);
    {
        uint64_t wt0 = bench_time_ns();
        do {
            run_circuit_quest(rho, c);
        } while (bench_ns_to_ms(bench_time_ns() - wt0) < BENCH_PQ_WARMUP_MS);
    }

    /* --- Timed runs ------------------------------------------------------- */
    double run_times[BENCH_MAX_RUNS];
    int cnt = (runs < BENCH_MAX_RUNS) ? runs : BENCH_MAX_RUNS;

    for (int r = 0; r < cnt; ++r) {
        /* Re-initialise to maximally mixed state before each run */
        init_maximally_mixed(rho);

        bench_timing_barrier();
        uint64_t t0 = bench_time_ns();
        bench_timing_barrier();

        for (int i = 0; i < iterations; ++i) {
            run_circuit_quest(rho, c);
        }

        bench_timing_barrier();
        run_times[r] = bench_ns_to_ms(bench_time_ns() - t0);
    }

    /* --- Statistics ------------------------------------------------------- */
    bench_run_stats_t stats = bench_compute_stats(run_times, cnt);
    result.time_ms        = stats.mean;
    result.time_ms_std    = stats.std_dev;
    result.time_ms_min    = stats.min;
    result.time_ms_median = stats.median;
    result.time_ms_cv     = stats.cv;

    double total_gate_apps = (double)iterations * (double)c->n_ops;
    result.ops_per_sec = (stats.mean > 0.0)
        ? total_gate_apps / (stats.mean / 1000.0)
        : 0.0;
    result.time_per_gate_ms = (total_gate_apps > 0.0)
        ? stats.mean / total_gate_apps
        : 0.0;

    destroyQureg(rho);
    return result;
}

#endif /* WITH_QUEST */
