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
 * Memory measurement: We fork a child process for each benchmark to get
 * accurate memory measurements, avoiding memory reuse between benchmarks.
 *
 * Statistics: the child performs `runs` independent timed rounds and returns
 * the full bench_run_stats_t summary through the pipe.
 */

#ifdef WITH_QULACS

#include <csim/update_ops_dm.hpp>
#include <csim/type.hpp>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <sys/wait.h>

extern "C" {
#include "bench.h"
}

/*
 * =====================================================================================================================
 * Qulacs benchmark implementation
 * =====================================================================================================================
 */

/**
 * @brief Run Qulacs benchmark in child process (internal)
 *
 * Performs `runs` independent timed rounds and fills all statistical fields
 * in result before writing it to the pipe.
 */
static void bench_qulacs_child(int write_fd, qubit_t qubits, const char *gate_name,
                               int iterations, int warmup, int runs) {
    bench_result_t result = {0};
    result.qubits     = qubits;
    result.dim        = (dim_t)1 << qubits;
    result.iterations = iterations;
    result.runs       = runs;

    /* Allocate full density matrix (dim x dim complex doubles) */
    ITYPE dim = static_cast<ITYPE>(1) << qubits;
    size_t matrix_size = static_cast<size_t>(dim) * static_cast<size_t>(dim) * sizeof(CTYPE);
    CTYPE* rho = static_cast<CTYPE*>(std::malloc(matrix_size));
    if (!rho) {
        write(write_fd, &result, sizeof(result));
        close(write_fd);
        _exit(1);
    }

    /*
     * Initialize to |0><0| state and touch all memory pages.
     * The density matrix is stored in row-major order: rho[i*dim + j]
     */
    std::memset(rho, 0, matrix_size);
    rho[0] = CTYPE(1.0, 0.0);

    /*
     * Memory measurement: Use theoretical size since RSS measurement after fork()
     * is unreliable on macOS (child inherits parent's COW pages, deltas are wrong).
     */
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
        bench_pipe_result_t pr = {0};
        write(write_fd, &pr, sizeof(pr));
        close(write_fd);
        _exit(1);
    }

    int cnt = (runs < 200) ? runs : 200;
    double run_times[200];

    if (gate_fn_1q) {
        /* Single-qubit gate: `runs` independent timed rounds */
        for (int i = 0; i < warmup; ++i)
            gate_fn_1q(static_cast<UINT>(i % qubits), rho, dim);

        for (int r = 0; r < cnt; ++r) {
            auto start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < iterations; ++i)
                for (qubit_t t = 0; t < qubits; ++t)
                    gate_fn_1q(static_cast<UINT>(t), rho, dim);
            auto end = std::chrono::high_resolution_clock::now();
            run_times[r] = static_cast<double>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count())
                / 1000000.0;
        }

        bench_run_stats_t stats = bench_compute_stats(run_times, cnt);
        result.time_ms        = stats.mean;
        result.time_ms_std    = stats.std_dev;
        result.time_ms_min    = stats.min;
        result.time_ms_median = stats.median;
        result.time_ms_cv     = stats.cv;
        result.ops_per_sec    = static_cast<double>(iterations * qubits) / (stats.mean / 1000.0);
    } else {
        /* Two-qubit gate: all ordered pairs (q1 < q2) */
        qubit_t q1s[128], q2s[128];
        int num_pairs = 0;
        for (qubit_t q1 = 0; q1 < qubits; ++q1)
            for (qubit_t q2 = q1 + 1; q2 < qubits; ++q2) {
                if (num_pairs >= 128) break;
                q1s[num_pairs] = q1;
                q2s[num_pairs] = q2;
                num_pairs++;
            }

        if (num_pairs == 0) {
            std::free(rho);
            bench_pipe_result_t pr = {0};
            write(write_fd, &pr, sizeof(pr));
            close(write_fd);
            _exit(0);
        }

        for (int i = 0; i < warmup; ++i)
            gate_fn_2q(static_cast<UINT>(q1s[i % num_pairs]),
                       static_cast<UINT>(q2s[i % num_pairs]), rho, dim);

        for (int r = 0; r < cnt; ++r) {
            auto start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < iterations; ++i)
                for (int p = 0; p < num_pairs; ++p)
                    gate_fn_2q(static_cast<UINT>(q1s[p]), static_cast<UINT>(q2s[p]), rho, dim);
            auto end = std::chrono::high_resolution_clock::now();
            run_times[r] = static_cast<double>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count())
                / 1000000.0;
        }

        bench_run_stats_t stats = bench_compute_stats(run_times, cnt);
        result.time_ms        = stats.mean;
        result.time_ms_std    = stats.std_dev;
        result.time_ms_min    = stats.min;
        result.time_ms_median = stats.median;
        result.time_ms_cv     = stats.cv;
        result.ops_per_sec    = static_cast<double>(iterations * num_pairs) / (stats.mean / 1000.0);
    }

    std::free(rho);

    /* Write scalar fields to parent via a pipe-safe struct.
     * bench_result_t contains pointer fields (gate_name, method) that are
     * meaningless in the parent address space and must not be transmitted. */
    bench_pipe_result_t pr;
    pr.time_ms        = result.time_ms;
    pr.time_ms_std    = result.time_ms_std;
    pr.time_ms_min    = result.time_ms_min;
    pr.time_ms_median = result.time_ms_median;
    pr.time_ms_cv     = result.time_ms_cv;
    pr.ops_per_sec    = result.ops_per_sec;
    pr.memory_bytes   = result.memory_bytes;
    write(write_fd, &pr, sizeof(pr));
    close(write_fd);
    _exit(0);
}

extern "C" {

/**
 * @brief Benchmark Qulacs density matrix gate operations.
 *
 * Forks a child process; the child performs `runs` independent timed rounds
 * and returns the full statistical summary via pipe.
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

    int pipefd[2];
    if (pipe(pipefd) == -1) {
        perror("bench_qulacs: pipe failed");
        return result;
    }

    pid_t pid = fork();
    if (pid == -1) {
        perror("bench_qulacs: fork failed");
        close(pipefd[0]);
        close(pipefd[1]);
        return result;
    }

    if (pid == 0) {
        close(pipefd[0]);
        bench_qulacs_child(pipefd[1], qubits, gate_name, iterations, warmup, runs);
        _exit(0);  /* unreachable */
    }

    close(pipefd[1]);

    bench_pipe_result_t pr;
    ssize_t n = read(pipefd[0], &pr, sizeof(pr));
    close(pipefd[0]);

    int status;
    waitpid(pid, &status, 0);

    if (n == static_cast<ssize_t>(sizeof(pr))) {
        result.time_ms        = pr.time_ms;
        result.time_ms_std    = pr.time_ms_std;
        result.time_ms_min    = pr.time_ms_min;
        result.time_ms_median = pr.time_ms_median;
        result.time_ms_cv     = pr.time_ms_cv;
        result.ops_per_sec    = pr.ops_per_sec;
        result.memory_bytes   = pr.memory_bytes;
    }

    return result;
}

} /* extern "C" */

#endif /* WITH_QULACS */
