#ifdef WITH_AER_DM

/*
 * bench_circuit_aer_dm.cpp — Qiskit-Aer density matrix circuit runner
 *
 * Implements bench_circuit_aer_dm() and context alloc/free for the Aer-DM
 * backend. Uses the same AER::QV::DensityMatrix<double> API as bench_aer_dm.cpp.
 *
 * Rotation gates (Rx, Ry, Rz) require per-op matrix construction since angles
 * vary (especially in QV). This overhead is inherent to the Aer-DM API and is
 * part of what we measure — it is not a benchmark artifact.
 *
 * Fixed gates (H, X, Z, CX) use precomputed vectorized matrices or direct
 * apply methods where available. Z uses apply_phase(q, -1).
 */

#include "simulators/density_matrix/densitymatrix.hpp"
#include "framework/linalg/matrix_utils.hpp"
#include <complex>
#include <cstdio>
#include <cmath>
#include <cassert>

#include "bench_circuit.h"

/* ------------------------------------------------------------------------------------------------------------------ *
 * Precomputed gate matrices — built once per bench_circuit_aer_dm call
 * ------------------------------------------------------------------------------------------------------------------ */

struct aer_dm_precomp {
    AER::cvector_t h_vec;
    AER::cvector_t x_vec;
    AER::cvector_t cx_vec;
};

static aer_dm_precomp build_precomp(void) {
    aer_dm_precomp p;

    /* H gate */
    p.h_vec = AER::Utils::vectorize_matrix(AER::Linalg::Matrix::H);

    /* X gate: [[0,1],[1,0]] */
    AER::cmatrix_t x_mat(2, 2);
    x_mat(0, 0) = {0, 0}; x_mat(0, 1) = {1, 0};
    x_mat(1, 0) = {1, 0}; x_mat(1, 1) = {0, 0};
    p.x_vec = AER::Utils::vectorize_matrix(x_mat);

    /* Z gate uses apply_phase directly; no precomp needed */

    /* CX gate: 4x4 matrix
     * [[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]] */
    AER::cmatrix_t cx_mat(4, 4);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            cx_mat(i, j) = {0, 0};
    cx_mat(0, 0) = {1, 0};
    cx_mat(1, 3) = {1, 0};
    cx_mat(2, 2) = {1, 0};
    cx_mat(3, 1) = {1, 0};
    p.cx_vec = AER::Utils::vectorize_matrix(cx_mat);

    return p;
}

/* ------------------------------------------------------------------------------------------------------------------ *
 * Per-op dispatch
 * ------------------------------------------------------------------------------------------------------------------ */

static inline void apply_circuit_op_aer_dm(AER::QV::DensityMatrix<double> &dm,
                                            const circuit_op_t *op,
                                            const aer_dm_precomp &pre) {
    switch (op->gate_type) {
    case CIRCUIT_GATE_H: {
        const AER::reg_t reg = {static_cast<AER::uint_t>(op->q0)};
        dm.apply_unitary_matrix(reg, pre.h_vec);
        break;
    }
    case GATE_X: {
        dm.apply_x(static_cast<AER::uint_t>(op->q0));
        break;
    }
    case GATE_Z: {
        dm.apply_phase(static_cast<AER::uint_t>(op->q0), {-1.0, 0.0});
        break;
    }
    case GATE_RX: {
        double c = std::cos(op->angle / 2.0);
        double s = std::sin(op->angle / 2.0);
        AER::cmatrix_t m(2, 2);
        m(0, 0) = {c, 0};      m(0, 1) = {0, -s};
        m(1, 0) = {0, -s};     m(1, 1) = {c, 0};
        const AER::reg_t reg = {static_cast<AER::uint_t>(op->q0)};
        dm.apply_unitary_matrix(reg, AER::Utils::vectorize_matrix(m));
        break;
    }
    case GATE_RY: {
        double c = std::cos(op->angle / 2.0);
        double s = std::sin(op->angle / 2.0);
        AER::cmatrix_t m(2, 2);
        m(0, 0) = {c, 0};      m(0, 1) = {-s, 0};
        m(1, 0) = {s, 0};      m(1, 1) = {c, 0};
        const AER::reg_t reg = {static_cast<AER::uint_t>(op->q0)};
        dm.apply_unitary_matrix(reg, AER::Utils::vectorize_matrix(m));
        break;
    }
    case GATE_RZ: {
        double h = op->angle / 2.0;
        AER::cmatrix_t m(2, 2);
        m(0, 0) = {std::cos(-h), std::sin(-h)};  m(0, 1) = {0, 0};
        m(1, 0) = {0, 0};                         m(1, 1) = {std::cos(h), std::sin(h)};
        const AER::reg_t reg = {static_cast<AER::uint_t>(op->q0)};
        dm.apply_unitary_matrix(reg, AER::Utils::vectorize_matrix(m));
        break;
    }
    case GATE_CX: {
        dm.apply_cnot(static_cast<AER::uint_t>(op->q0),
                      static_cast<AER::uint_t>(op->q1));
        break;
    }
    default:
        fprintf(stderr, "bench_circuit_aer_dm: unimplemented gate type %d\n",
                (int)op->gate_type);
        break;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ *
 * Run a full circuit on dm for `iterations` repetitions (no reinit between iterations)
 * ------------------------------------------------------------------------------------------------------------------ */

static void run_circuit_aer_dm(AER::QV::DensityMatrix<double> &dm,
                                const circuit_t *c,
                                const aer_dm_precomp &pre,
                                int iterations) {
    for (int iter = 0; iter < iterations; ++iter) {
        for (int g = 0; g < c->n_ops; ++g) {
            apply_circuit_op_aer_dm(dm, &c->ops[g], pre);
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ *
 * Maximally mixed state initialization: rho = I / dim
 * ------------------------------------------------------------------------------------------------------------------ */

static void init_maximally_mixed(AER::QV::DensityMatrix<double> &dm, int n_qubits) {
    dm.initialize();  /* zeros out the internal buffer */
    size_t dim = (size_t)1 << n_qubits;
    double diag_val = 1.0 / (double)dim;
    auto *buf = dm.data();
    /* Density matrix is stored as dim x dim row-major complex array */
    for (size_t i = 0; i < dim; ++i) {
        buf[i * dim + i] = std::complex<double>(diag_val, 0.0);
    }
}

/* ------------------------------------------------------------------------------------------------------------------ *
 * Extern "C" exports
 * ------------------------------------------------------------------------------------------------------------------ */

extern "C" {

aer_dm_circuit_ctx_t *aer_dm_circuit_ctx_alloc(int n_qubits) {
    try {
        auto *ctx = new aer_dm_circuit_ctx_t;
        auto *dm = new AER::QV::DensityMatrix<double>(static_cast<unsigned>(n_qubits));
        init_maximally_mixed(*dm, n_qubits);
        ctx->dm       = dm;
        ctx->n_qubits = n_qubits;
        return ctx;
    } catch (...) {
        return nullptr;
    }
}

void aer_dm_circuit_ctx_free(aer_dm_circuit_ctx_t *ctx) {
    if (!ctx) return;
    delete static_cast<AER::QV::DensityMatrix<double> *>(ctx->dm);
    delete ctx;
}

bench_result_t bench_circuit_aer_dm(const circuit_t *c, int iterations,
                                     int warmup, int runs) {
    bench_result_t result = {};
    result.qubits     = (qubit_t)c->n_qubits;
    result.dim        = (dim_t)1 << c->n_qubits;
    result.gate_name  = c->name;
    result.method     = "qiskit_aer_dm";
    result.iterations = iterations;
    result.runs       = runs;
    result.sweep_size = c->n_ops;

    try {
        if (c->n_qubits == 0 || c->n_ops == 0) return result;

        int clamped_runs = (runs > BENCH_MAX_RUNS) ? BENCH_MAX_RUNS : runs;

        AER::QV::DensityMatrix<double> dm(static_cast<unsigned>(c->n_qubits));

        /* Precompute fixed gate matrices once */
        aer_dm_precomp pre = build_precomp();

        size_t dim = (size_t)1 << c->n_qubits;
        result.memory_bytes = dim * dim * sizeof(std::complex<double>);

        /* Time-based warmup */
        init_maximally_mixed(dm, c->n_qubits);
        uint64_t wt0 = bench_time_ns();
        do {
            run_circuit_aer_dm(dm, c, pre, 1);
        } while (bench_ns_to_ms(bench_time_ns() - wt0) < BENCH_PQ_WARMUP_MS);

        /* Timed runs — reinit before each run, not between iterations */
        double run_times[BENCH_MAX_RUNS];
        assert(clamped_runs <= BENCH_MAX_RUNS);

        for (int r = 0; r < clamped_runs; ++r) {
            /* Re-initialize to maximally mixed state before each timed run */
            init_maximally_mixed(dm, c->n_qubits);

            bench_timing_barrier();
            uint64_t t0 = bench_time_ns();
            bench_timing_barrier();

            run_circuit_aer_dm(dm, c, pre, iterations);

            bench_timing_barrier();
            run_times[r] = bench_ns_to_ms(bench_time_ns() - t0);
        }

        bench_run_stats_t stats = bench_compute_stats(run_times, clamped_runs);
        result.time_ms        = stats.mean;
        result.time_ms_std    = stats.std_dev;
        result.time_ms_min    = stats.min;
        result.time_ms_median = stats.median;
        result.time_ms_cv     = stats.cv;
        result.runs           = clamped_runs;
        result.ops_per_sec    = (stats.mean > 0.0)
            ? static_cast<double>(iterations) * c->n_ops / (stats.mean / 1000.0)
            : 0.0;
        result.time_per_gate_ms = (iterations > 0 && c->n_ops > 0)
            ? stats.mean / ((double)iterations * c->n_ops)
            : 0.0;

    } catch (const std::exception &e) {
        fprintf(stderr, "bench_circuit_aer_dm: exception for %d qubits, circuit %s: %s\n",
                c->n_qubits, c->name ? c->name : "(null)", e.what());
        return result;
    } catch (...) {
        fprintf(stderr, "bench_circuit_aer_dm: unknown exception for %d qubits, circuit %s\n",
                c->n_qubits, c->name ? c->name : "(null)");
        return result;
    }

    return result;
}

} /* extern "C" */

#endif /* WITH_AER_DM */
