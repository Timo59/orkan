#ifdef WITH_AER_DM

#include "simulators/density_matrix/densitymatrix.hpp"
#include "framework/linalg/matrix_utils.hpp"
#include <complex>
#include <cstdio>
#include <cstring>
#include <random>
#include <cassert>

#include "bench.h"

/*
 * =====================================================================================================================
 * Qiskit-Aer DensityMatrix benchmark implementation
 * =====================================================================================================================
 */

/* File-scope gate kind enum — hoisted so it is visible to all functions in this TU.
 * Prefix BK_ (bench kind) to avoid colliding with gate.h include guard macro GATE_H. */
enum gate_kind_t { BK_X, BK_H, BK_Z, BK_CX, BK_SWAP, BK_UNKNOWN };

/* Resolve a gate name string to gate_kind_t once, outside any hot path. */
static gate_kind_t resolve_gate_kind(const char *gate_name) {
    if      (std::strcmp(gate_name, "X")    == 0) return BK_X;
    else if (std::strcmp(gate_name, "H")    == 0) return BK_H;
    else if (std::strcmp(gate_name, "Z")    == 0) return BK_Z;
    else if (std::strcmp(gate_name, "CX")   == 0) return BK_CX;
    else if (std::strcmp(gate_name, "SWAP") == 0) return BK_SWAP;
    else                                           return BK_UNKNOWN;
}

/* ------------------------------------------------------------------------------------------------------------------ *
 * Callback context and callback for per-qubit harness
 * ------------------------------------------------------------------------------------------------------------------ */

struct aer_dm_cb_ctx {
    AER::QV::DensityMatrix<double> *dm;
    gate_kind_t                     gate_kind;
    int                             q1;
    int                             q2;
    /* Pre-built vectorized matrices for gates that lack direct named methods */
    AER::cvector_t                  h_vec;
    AER::cvector_t                  z_vec;
    /* Cached reg_t objects to avoid per-call heap allocation */
    AER::reg_t                      q1_reg;
    AER::reg_t                      q2_reg;
};

static void aer_dm_cb(void *ctx) {
    auto *c = static_cast<aer_dm_cb_ctx *>(ctx);
    switch (c->gate_kind) {
    case BK_X:
        c->dm->apply_x(static_cast<AER::uint_t>(c->q1));
        break;
    case BK_H:
        c->dm->apply_unitary_matrix(c->q1_reg, c->h_vec);
        break;
    case BK_Z:
        c->dm->apply_unitary_matrix(c->q1_reg, c->z_vec);
        break;
    case BK_CX:
        c->dm->apply_cnot(static_cast<AER::uint_t>(c->q1),
                          static_cast<AER::uint_t>(c->q2));
        break;
    case BK_SWAP:
        c->dm->apply_swap(static_cast<AER::uint_t>(c->q1),
                          static_cast<AER::uint_t>(c->q2));
        break;
    default:
        break;
    }
}

extern "C" {

/* ------------------------------------------------------------------------------------------------------------------ *
 * Original sweep benchmark — unchanged
 * ------------------------------------------------------------------------------------------------------------------ */

bench_result_t bench_aer_dm(qubit_t qubits, const char *gate_name,
                             int iterations, int warmup, int runs) {
    bench_result_t result = {0};
    result.qubits     = qubits;
    result.dim        = (dim_t)1 << qubits;
    result.gate_name  = gate_name;
    result.method     = "qiskit_aer_dm";
    result.iterations = iterations;
    result.runs       = runs;

    try {
        if (qubits == 0) return result;

        /* Fix 3: clamp runs to avoid overflowing run_times stack buffer */
        if (runs > BENCH_MAX_RUNS) runs = BENCH_MAX_RUNS;

        AER::QV::DensityMatrix<double> dm(qubits);
        dm.initialize();  /* allocates internal buffer; overwritten below */

        /* Random init — faults all pages and avoids trivial no-op gate paths */
        std::mt19937_64 rng(static_cast<uint64_t>(qubits) * 6364136223846793005ULL
                            + 1442695040888963407ULL);
        std::uniform_real_distribution<double> dist(-0.5, 0.5);
        size_t n_elems = (size_t)result.dim * (size_t)result.dim;
        auto *buf = dm.data();
        for (size_t i = 0; i < n_elems; ++i)
            buf[i] = std::complex<double>(dist(rng), dist(rng));
        /* Intentionally invalid density matrix (not Hermitian/trace-1): ensures
         * all memory pages are faulted and gate fast-paths for special states are bypassed. */
        result.memory_bytes = n_elems * sizeof(std::complex<double>);

        /* Precompute vectorized matrices for gates that lack direct methods */
        AER::cvector_t h_vec;
        AER::cvector_t z_vec;

        /* Fix 1: resolve gate type once before warmup/timing loops — no strcmp in hot path */
        gate_kind_t gate_kind = resolve_gate_kind(gate_name);

        if (gate_kind == BK_UNKNOWN)
            return result;

        if (gate_kind == BK_H) {
            h_vec = AER::Utils::vectorize_matrix(AER::Linalg::Matrix::H);
        }
        if (gate_kind == BK_Z) {
            AER::cmatrix_t z_mat(2, 2);
            z_mat(0, 0) = {1, 0}; z_mat(0, 1) = {0, 0};
            z_mat(1, 0) = {0, 0}; z_mat(1, 1) = {-1, 0};
            z_vec = AER::Utils::vectorize_matrix(z_mat);
        }

        bool is_1q = (gate_kind == BK_X || gate_kind == BK_H || gate_kind == BK_Z);
        bool is_2q = (gate_kind == BK_CX || gate_kind == BK_SWAP);

        double run_times[BENCH_MAX_RUNS];

        if (is_1q) {
            /* Single-qubit gate: sweep all qubits per iteration */

            /* Pre-build qubit register objects for apply_unitary_matrix calls (avoids per-call heap allocation) */
            std::vector<AER::reg_t> qubit_regs(qubits);
            for (qubit_t t = 0; t < qubits; ++t)
                qubit_regs[t] = AER::reg_t{static_cast<AER::uint_t>(t)};

            for (int i = 0; i < warmup; ++i) {
                qubit_t t = static_cast<qubit_t>(i % qubits);
                switch (gate_kind) {
                case BK_X:
                    dm.apply_x(t);
                    break;
                case BK_H:
                    dm.apply_unitary_matrix(qubit_regs[t], h_vec);
                    break;
                case BK_Z:
                    dm.apply_unitary_matrix(qubit_regs[t], z_vec);
                    break;
                default:
                    break;
                }
            }

            for (int r = 0; r < runs; ++r) {
                uint64_t start = bench_time_ns();
                for (int i = 0; i < iterations; ++i) {
                    for (qubit_t t = 0; t < qubits; ++t) {
                        switch (gate_kind) {
                        case BK_X:
                            dm.apply_x(t);
                            break;
                        case BK_H:
                            dm.apply_unitary_matrix(qubit_regs[t], h_vec);
                            break;
                        case BK_Z:
                            dm.apply_unitary_matrix(qubit_regs[t], z_vec);
                            break;
                        default:
                            break;
                        }
                    }
                }
                run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
            }

            bench_run_stats_t stats = bench_compute_stats(run_times, runs);
            result.time_ms        = stats.mean;
            result.time_ms_std    = stats.std_dev;
            result.time_ms_min    = stats.min;
            result.time_ms_median = stats.median;
            result.time_ms_cv     = stats.cv;
            result.ops_per_sec    = (stats.mean > 0.0)
                ? static_cast<double>(iterations * qubits) / (stats.mean / 1000.0)
                : 0.0;
        } else if (is_2q) {
            /* Two-qubit gate: all ordered pairs (q1 < q2) */
            qubit_t q1s[MAX_PAIRS], q2s[MAX_PAIRS];
            int num_pairs = build_all_pairs(qubits, q1s, q2s);

            if (num_pairs == 0)
                return result;

            for (int i = 0; i < warmup; ++i) {
                int p = i % num_pairs;
                switch (gate_kind) {
                case BK_CX:
                    dm.apply_cnot(q1s[p], q2s[p]);
                    break;
                case BK_SWAP:
                    dm.apply_swap(q1s[p], q2s[p]);
                    break;
                default:
                    break;
                }
            }

            for (int r = 0; r < runs; ++r) {
                uint64_t start = bench_time_ns();
                for (int i = 0; i < iterations; ++i) {
                    for (int p = 0; p < num_pairs; ++p) {
                        switch (gate_kind) {
                        case BK_CX:
                            dm.apply_cnot(q1s[p], q2s[p]);
                            break;
                        case BK_SWAP:
                            dm.apply_swap(q1s[p], q2s[p]);
                            break;
                        default:
                            break;
                        }
                    }
                }
                run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
            }

            bench_run_stats_t stats = bench_compute_stats(run_times, runs);
            result.time_ms        = stats.mean;
            result.time_ms_std    = stats.std_dev;
            result.time_ms_min    = stats.min;
            result.time_ms_median = stats.median;
            result.time_ms_cv     = stats.cv;
            result.ops_per_sec    = (stats.mean > 0.0)
                ? static_cast<double>(iterations * num_pairs) / (stats.mean / 1000.0)
                : 0.0;
        }
    } catch (const std::exception &e) {
        fprintf(stderr, "bench_aer_dm: exception for %u qubits, gate %s: %s\n",
                (unsigned)qubits, gate_name, e.what());
        return result;
    } catch (...) {
        fprintf(stderr, "bench_aer_dm: unknown exception for %u qubits, gate %s\n",
                (unsigned)qubits, gate_name);
        return result;
    }
    return result;
}

/* ------------------------------------------------------------------------------------------------------------------ *
 * Allocator / reinit / free helpers
 * ------------------------------------------------------------------------------------------------------------------ */

void *bench_aer_dm_alloc(qubit_t qubits) {
    try {
        auto *dm = new AER::QV::DensityMatrix<double>(qubits);
        dm->initialize();
        return dm;
    } catch (...) { return nullptr; }
}

void bench_aer_dm_reinit(void *dm, qubit_t qubits) {
    auto *d = static_cast<AER::QV::DensityMatrix<double>*>(dm);
    d->initialize();
    /* Fill with pseudo-random bytes using bench_fill_random for consistency */
    size_t n_elements = (size_t)1 << (2 * qubits);  /* dim * dim complex elements */
    bench_fill_random(d->data(), n_elements * sizeof(std::complex<double>));
}

void bench_aer_dm_free(void *dm) {
    delete static_cast<AER::QV::DensityMatrix<double> *>(dm);
}

/* ------------------------------------------------------------------------------------------------------------------ *
 * Per-qubit benchmark — single-qubit gate at a fixed target
 * ------------------------------------------------------------------------------------------------------------------ */

bench_result_perq_t bench_aer_dm_at(qubit_t qubits, const char *gate_name,
                                     qubit_t target, int iterations, int runs,
                                     void *dm) {
    bench_result_perq_t result = {};
    result.qubits    = qubits;
    result.dim       = (dim_t)1 << qubits;
    result.gate_name = gate_name;
    result.method    = "qiskit_aer_dm";
    result.target    = target;
    result.target2   = QUBIT_NONE;
    result.is_2q     = 0;

    try {
        gate_kind_t gate_kind = resolve_gate_kind(gate_name);
        if (gate_kind == BK_UNKNOWN || gate_kind == BK_CX || gate_kind == BK_SWAP)
            return result;

        /* Re-init state before measurement to fault all pages */
        bench_aer_dm_reinit(dm, qubits);

        /* Build callback context */
        aer_dm_cb_ctx ctx;
        ctx.dm        = static_cast<AER::QV::DensityMatrix<double> *>(dm);
        ctx.gate_kind = gate_kind;
        ctx.q1        = static_cast<int>(target);
        ctx.q2        = 0;  /* unused for 1Q */
        ctx.q1_reg    = AER::reg_t{static_cast<AER::uint_t>(target)};
        ctx.q2_reg    = AER::reg_t{};  /* unused for 1Q */

        if (gate_kind == BK_H)
            ctx.h_vec = AER::Utils::vectorize_matrix(AER::Linalg::Matrix::H);

        if (gate_kind == BK_Z) {
            AER::cmatrix_t z_mat(2, 2);
            z_mat(0, 0) = {1, 0}; z_mat(0, 1) = {0, 0};
            z_mat(1, 0) = {0, 0}; z_mat(1, 1) = {-1, 0};
            ctx.z_vec = AER::Utils::vectorize_matrix(z_mat);
        }

        int clamped_runs = (runs > BENCH_MAX_RUNS) ? BENCH_MAX_RUNS : runs;

        bench_harness_t h;
        h.call       = aer_dm_cb;
        h.ctx        = &ctx;
        h.iterations = iterations;
        h.runs       = clamped_runs;

        double run_times[BENCH_MAX_RUNS];
        assert(clamped_runs <= BENCH_MAX_RUNS);
        bench_run_timed(&h, run_times, (int)qubits);

        bench_run_stats_t stats = bench_compute_stats(run_times, clamped_runs);
        bench_fill_perq_stats(&result, &stats, iterations);
        result.runs         = clamped_runs;
        result.memory_bytes = ((size_t)1 << qubits) * ((size_t)1 << qubits)
                              * sizeof(std::complex<double>);
    } catch (const std::exception &e) {
        fprintf(stderr, "bench_aer_dm_at: exception: %s\n", e.what());
        return {};
    } catch (...) {
        fprintf(stderr, "bench_aer_dm_at: unknown exception\n");
        return {};
    }

    return result;
}

/* ------------------------------------------------------------------------------------------------------------------ *
 * Per-qubit benchmark — two-qubit gate at fixed (q1, q2)
 * ------------------------------------------------------------------------------------------------------------------ */

bench_result_perq_t bench_aer_dm_2q_at(qubit_t qubits, const char *gate_name,
                                        qubit_t q1, qubit_t q2,
                                        int iterations, int runs,
                                        void *dm) {
    bench_result_perq_t result = {};
    result.qubits    = qubits;
    result.dim       = (dim_t)1 << qubits;
    result.gate_name = gate_name;
    result.method    = "qiskit_aer_dm";
    result.target    = q1;
    result.target2   = q2;
    result.is_2q     = 1;

    try {
        gate_kind_t gate_kind = resolve_gate_kind(gate_name);
        if (gate_kind != BK_CX && gate_kind != BK_SWAP)
            return result;

        /* Re-init state before measurement */
        bench_aer_dm_reinit(dm, qubits);

        /* Build callback context */
        aer_dm_cb_ctx ctx;
        ctx.dm        = static_cast<AER::QV::DensityMatrix<double> *>(dm);
        ctx.gate_kind = gate_kind;
        ctx.q1        = static_cast<int>(q1);
        ctx.q2        = static_cast<int>(q2);
        ctx.q1_reg    = AER::reg_t{static_cast<AER::uint_t>(q1)};
        ctx.q2_reg    = AER::reg_t{static_cast<AER::uint_t>(q2)};
        /* h_vec / z_vec unused for 2Q gates */

        int clamped_runs = (runs > BENCH_MAX_RUNS) ? BENCH_MAX_RUNS : runs;

        bench_harness_t h;
        h.call       = aer_dm_cb;
        h.ctx        = &ctx;
        h.iterations = iterations;
        h.runs       = clamped_runs;

        double run_times[BENCH_MAX_RUNS];
        assert(clamped_runs <= BENCH_MAX_RUNS);
        bench_run_timed(&h, run_times, (int)qubits);

        bench_run_stats_t stats = bench_compute_stats(run_times, clamped_runs);
        bench_fill_perq_stats(&result, &stats, iterations);
        result.runs         = clamped_runs;
        result.memory_bytes = ((size_t)1 << qubits) * ((size_t)1 << qubits)
                              * sizeof(std::complex<double>);
    } catch (const std::exception &e) {
        fprintf(stderr, "bench_aer_dm_2q_at: exception: %s\n", e.what());
        return {};
    } catch (...) {
        fprintf(stderr, "bench_aer_dm_2q_at: unknown exception\n");
        return {};
    }

    return result;
}

} /* extern "C" */

#endif /* WITH_AER_DM */
