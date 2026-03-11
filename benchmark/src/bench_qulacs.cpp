/**
 * @file bench_qulacs.cpp
 * @brief Qulacs framework benchmark integration
 *
 * Qulacs is a fast quantum circuit simulator developed by QunaSys.
 * It provides optimized density matrix operations with SIMD and OpenMP support.
 *
 * Homepage: https://github.com/qulacs/qulacs
 *
 * To enable: cmake -DWITH_QULACS=ON ..
 * Requires: Boost 1.71+, Eigen3 (auto-fetched by Qulacs)
 *
 * Qulacs csim has no global state — no init/finalize, no singleton.
 * The benchmark runs directly in-process: malloc, compute, free, return.
 */

#ifdef WITH_QULACS

#include <csim/update_ops_dm.hpp>
#include <csim/type.hpp>
#include <cstdlib>
#include <cstring>

#include "bench.h"

/*
 * =====================================================================================================================
 * Qulacs benchmark implementation
 * =====================================================================================================================
 */

extern "C" {

/**
 * @brief Run Qulacs density matrix benchmark in-process.
 *
 * Qulacs csim has no global state; this is a direct malloc/compute/free call.
 * Performs `runs` independent timed rounds and returns the full statistical summary.
 */
bench_result_t bench_qulacs(qubit_t qubits, const char *gate_name,
                            int iterations, int warmup, int runs) {
    bench_result_t result = {0};
    result.qubits     = qubits;
    result.dim        = (dim_t)1 << qubits;
    result.gate_name  = gate_name;
    result.method     = "qulacs";
    result.iterations = iterations;
    result.runs       = runs;

    /* Allocate full density matrix (dim x dim complex doubles) */
    ITYPE dim = static_cast<ITYPE>(1) << qubits;
    size_t matrix_size = static_cast<size_t>(dim) * static_cast<size_t>(dim) * sizeof(CTYPE);
    CTYPE* rho = static_cast<CTYPE*>(std::malloc(matrix_size));
    if (!rho) {
        return result;
    }

    /*
     * Fill with random values to match qlib and BLAS baselines.
     * A computational basis state (e.g. |0><0|) leaves Z and CX effectively
     * as no-ops throughout the benchmark, making those timings meaningless.
     * Random initialisation also faults all pages before timing begins.
     */
    size_t n_elems = static_cast<size_t>(dim) * static_cast<size_t>(dim);
    for (size_t i = 0; i < n_elems; ++i) {
        double re = (double)std::rand() / RAND_MAX - 0.5;
        double im = (double)std::rand() / RAND_MAX - 0.5;
        rho[i] = CTYPE(re, im);
    }

    /* Memory: theoretical density matrix size. */
    result.memory_bytes = matrix_size;

    /* Dispatch: single-qubit gates */
    void (*gate_fn_1q)(UINT, CTYPE*, ITYPE) = nullptr;
    if (std::strcmp(gate_name, "X") == 0) {
        gate_fn_1q = dm_X_gate;
    } else if (std::strcmp(gate_name, "H") == 0) {
        gate_fn_1q = dm_H_gate;
    } else if (std::strcmp(gate_name, "Z") == 0) {
        gate_fn_1q = dm_Z_gate;
    }

    /* Dispatch: two-qubit gates */
    void (*gate_fn_2q)(UINT, UINT, CTYPE*, ITYPE) = nullptr;
    if (std::strcmp(gate_name, "CX") == 0) {
        gate_fn_2q = dm_CNOT_gate;
    } else if (std::strcmp(gate_name, "SWAP") == 0) {
        gate_fn_2q = dm_SWAP_gate;
    }

    if (!gate_fn_1q && !gate_fn_2q) {
        std::free(rho);
        return result;
    }

    int cnt = (runs < 200) ? runs : 200;
    double run_times[200];

    if (gate_fn_1q) {
        /* Single-qubit gate: `runs` independent timed rounds */
        for (int i = 0; i < warmup; ++i)
            gate_fn_1q(static_cast<UINT>(i % qubits), rho, dim);

        for (int r = 0; r < cnt; ++r) {
            uint64_t start = bench_time_ns();
            for (int i = 0; i < iterations; ++i)
                for (qubit_t t = 0; t < qubits; ++t)
                    gate_fn_1q(static_cast<UINT>(t), rho, dim);
            run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
        }

        bench_run_stats_t stats = bench_compute_stats(run_times, cnt);
        result.time_ms        = stats.mean;
        result.time_ms_std    = stats.std_dev;
        result.time_ms_min    = stats.min;
        result.time_ms_median = stats.median;
        result.time_ms_cv     = stats.cv;
        result.ops_per_sec    = (stats.mean > 0.0) ? static_cast<double>(iterations * qubits) / (stats.mean / 1000.0) : 0.0;
    } else {
        /* Two-qubit gate: all ordered pairs (q1 < q2) */
        qubit_t q1s[MAX_PAIRS], q2s[MAX_PAIRS];
        int num_pairs = build_all_pairs(qubits, q1s, q2s);

        if (num_pairs == 0) {
            std::free(rho);
            return result;
        }

        for (int i = 0; i < warmup; ++i)
            gate_fn_2q(static_cast<UINT>(q1s[i % num_pairs]),
                       static_cast<UINT>(q2s[i % num_pairs]), rho, dim);

        for (int r = 0; r < cnt; ++r) {
            uint64_t start = bench_time_ns();
            for (int i = 0; i < iterations; ++i)
                for (int p = 0; p < num_pairs; ++p)
                    gate_fn_2q(static_cast<UINT>(q1s[p]), static_cast<UINT>(q2s[p]), rho, dim);
            run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
        }

        bench_run_stats_t stats = bench_compute_stats(run_times, cnt);
        result.time_ms        = stats.mean;
        result.time_ms_std    = stats.std_dev;
        result.time_ms_min    = stats.min;
        result.time_ms_median = stats.median;
        result.time_ms_cv     = stats.cv;
        result.ops_per_sec    = (stats.mean > 0.0) ? static_cast<double>(iterations * num_pairs) / (stats.mean / 1000.0) : 0.0;
    }

    std::free(rho);
    return result;
}

} /* extern "C" */

#endif /* WITH_QULACS */
