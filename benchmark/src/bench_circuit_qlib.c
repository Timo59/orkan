/**
 * @file bench_circuit_qlib.c
 * @brief qlib packed and tiled circuit runners for circuit-level benchmarks.
 *
 * Provides bench_circuit_qlib_packed() and bench_circuit_qlib_tiled(), which
 * execute a circuit_t gate sequence on a qlib density matrix and return timing
 * statistics. Also provides context alloc/free for calibrate_iterations
 * compatibility via the runner_fn function pointer.
 *
 * State is initialized to the maximally mixed state rho = I/dim before each
 * timed run. Within a run, the state evolves across `iterations` circuit
 * applications -- this measures sustained throughput, the correct metric for
 * a VQO inner-loop context.
 */

#include "bench_circuit.h"

#include <assert.h>
#include <string.h>

/*
 * =====================================================================================================================
 * Gate dispatch
 * =====================================================================================================================
 */

/**
 * @brief Apply a single circuit operation to the state via switch dispatch.
 *
 * Covers all 7 gate types used by the circuit constructors (QAOA, VQE-HEA, QV).
 */
static inline void apply_op(state_t *s, const circuit_op_t *op) {
    switch (op->gate_type) {
        case GATE_H:  h(s, op->q0);              break;
        case GATE_X:  x(s, op->q0);              break;
        case GATE_Z:  z(s, op->q0);              break;
        case GATE_RX: rx(s, op->q0, op->angle);  break;
        case GATE_RY: ry(s, op->q0, op->angle);  break;
        case GATE_RZ: rz(s, op->q0, op->angle);  break;
        case GATE_CX: cx(s, op->q0, op->q1);     break;
        default:
            fprintf(stderr, "bench_circuit_qlib: unimplemented gate type %d\n",
                    (int)op->gate_type);
            break;
    }
}

/*
 * =====================================================================================================================
 * Maximally mixed state initialization
 * =====================================================================================================================
 */

/**
 * @brief Initialize the state buffer to rho = I/dim (maximally mixed).
 *
 * Zeroes the entire buffer, then sets diagonal elements rho(i,i) = 1/dim.
 * This is NOT |+><+|/dim (which has all elements = 1/dim); it is the identity
 * matrix scaled by 1/dim. Uses bench_packed_size / bench_tiled_size for the
 * buffer size to ensure the full storage footprint (including padding) is zeroed.
 *
 * @param s      State with data already allocated and type/qubits set.
 * @param qubits Number of qubits.
 * @param type   MIXED_PACKED or MIXED_TILED.
 */
static void init_maximally_mixed(state_t *s, qubit_t qubits, state_type_t type) {
    size_t nbytes = (type == MIXED_PACKED)
                  ? bench_packed_size(qubits)
                  : bench_tiled_size(qubits);
    memset(s->data, 0, nbytes);

    dim_t dim = (dim_t)1 << qubits;
    cplx_t diag_val = 1.0 / (double)dim;
    for (dim_t i = 0; i < dim; ++i) {
        state_set(s, i, i, diag_val);
    }
}

/*
 * =====================================================================================================================
 * Shared circuit runner (parameterized by state type)
 * =====================================================================================================================
 */

/**
 * @brief Run the circuit benchmark for a given state type.
 *
 * Allocates the state once, then for each of `runs` timed runs: re-initializes
 * to maximally mixed, applies the circuit `iterations` times, and records the
 * elapsed wall-clock time. A time-based warmup precedes the timed runs.
 */
static bench_result_t run_circuit(const circuit_t *c, state_type_t type,
                                  int iterations, int warmup, int runs) {
    (void)warmup; /* warmup is time-based, not count-based */

    const char *method = (type == MIXED_PACKED) ? "qlib_packed" : "qlib_tiled";
    size_t mem = (type == MIXED_PACKED)
               ? bench_packed_size((qubit_t)c->n_qubits)
               : bench_tiled_size((qubit_t)c->n_qubits);

    bench_result_t result = {0};
    result.qubits       = (qubit_t)c->n_qubits;
    result.dim          = (dim_t)1 << c->n_qubits;
    result.gate_name    = c->name;
    result.method       = method;
    result.iterations   = iterations;
    result.sweep_size   = c->n_ops;
    result.memory_bytes = mem;

    /* Allocate state */
    state_t state = {0};
    state.type = type;
    state_init(&state, (qubit_t)c->n_qubits, NULL);
    if (!state.data) {
        fprintf(stderr, "bench_circuit_qlib: state allocation failed for %d qubits (%s)\n",
                c->n_qubits, method);
        return result;
    }

    /* Time-based warmup: run full circuit passes until BENCH_PQ_WARMUP_MS elapsed */
    init_maximally_mixed(&state, (qubit_t)c->n_qubits, type);
    uint64_t wt0 = bench_time_ns();
    do {
        for (int op = 0; op < c->n_ops; ++op)
            apply_op(&state, &c->ops[op]);
    } while (bench_ns_to_ms(bench_time_ns() - wt0) < BENCH_PQ_WARMUP_MS);

    /* Timed runs */
    assert(runs >= 1 && runs <= BENCH_MAX_RUNS);
    double run_times[BENCH_MAX_RUNS];

    for (int r = 0; r < runs; ++r) {
        /* Re-initialize to maximally mixed before each run */
        init_maximally_mixed(&state, (qubit_t)c->n_qubits, type);

        bench_timing_barrier();
        uint64_t t0 = bench_time_ns();
        bench_timing_barrier();

        for (int i = 0; i < iterations; ++i) {
            for (int op = 0; op < c->n_ops; ++op) {
                apply_op(&state, &c->ops[op]);
            }
        }

        bench_timing_barrier();
        run_times[r] = bench_ns_to_ms(bench_time_ns() - t0);
    }

    /* Compute statistics */
    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    result.time_ms        = stats.mean;
    result.time_ms_std    = stats.std_dev;
    result.time_ms_min    = stats.min;
    result.time_ms_median = stats.median;
    result.time_ms_cv     = stats.cv;
    result.runs           = runs;

    double total_ops = (double)iterations * c->n_ops;
    result.ops_per_sec    = (stats.mean > 0.0)
                          ? total_ops / (stats.mean / 1000.0)
                          : 0.0;
    result.time_per_gate_ms = (total_ops > 0.0)
                            ? stats.mean / total_ops
                            : 0.0;

    state_free(&state);
    return result;
}

/*
 * =====================================================================================================================
 * Public runner functions
 * =====================================================================================================================
 */

bench_result_t bench_circuit_qlib_packed(const circuit_t *c, int iterations,
                                         int warmup, int runs) {
    return run_circuit(c, MIXED_PACKED, iterations, warmup, runs);
}

bench_result_t bench_circuit_qlib_tiled(const circuit_t *c, int iterations,
                                        int warmup, int runs) {
    return run_circuit(c, MIXED_TILED, iterations, warmup, runs);
}

/*
 * =====================================================================================================================
 * Context alloc/free for calibrate_iterations compatibility
 * =====================================================================================================================
 */

qlib_circuit_ctx_t *qlib_circuit_ctx_alloc(int n_qubits, state_type_t type) {
    qlib_circuit_ctx_t *ctx = malloc(sizeof(*ctx));
    if (!ctx) return NULL;

    ctx->type = type;
    ctx->s = malloc(sizeof(state_t));
    if (!ctx->s) {
        free(ctx);
        return NULL;
    }

    *ctx->s = (state_t){0};
    ctx->s->type = type;
    state_init(ctx->s, (qubit_t)n_qubits, NULL);
    if (!ctx->s->data) {
        free(ctx->s);
        free(ctx);
        return NULL;
    }

    /* Initialize to maximally mixed */
    init_maximally_mixed(ctx->s, (qubit_t)n_qubits, type);
    return ctx;
}

void qlib_circuit_ctx_free(qlib_circuit_ctx_t *ctx) {
    if (!ctx) return;
    if (ctx->s) {
        state_free(ctx->s);
        free(ctx->s);
    }
    free(ctx);
}
