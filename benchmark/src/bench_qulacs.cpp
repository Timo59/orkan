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
 * This runs in a forked child and writes results to the pipe.
 */
static void bench_qulacs_child(int write_fd, qubit_t qubits, const char *gate_name,
                               int iterations, int warmup) {
    bench_result_t result = {0};
    result.qubits = qubits;
    result.dim = (dim_t)1 << qubits;
    result.iterations = iterations;

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
     * We set all elements to zero first, then set rho[0,0] = 1.
     * The density matrix is stored in row-major order: rho[i*dim + j]
     */
    std::memset(rho, 0, matrix_size);
    rho[0] = CTYPE(1.0, 0.0);

    /*
     * Memory measurement: Use theoretical size since RSS measurement after fork()
     * is unreliable on macOS (child inherits parent's COW pages, deltas are wrong).
     * The density matrix is the dominant allocation: dim^2 * sizeof(complex<double>).
     */
    result.memory_bytes = matrix_size;

    /* Select gate function based on name */
    void (*gate_fn)(UINT, CTYPE*, ITYPE) = nullptr;
    if (gate_name[0] == 'X') {
        gate_fn = dm_X_gate;
    } else if (gate_name[0] == 'H') {
        gate_fn = dm_H_gate;
    } else if (gate_name[0] == 'Z') {
        gate_fn = dm_Z_gate;
    } else {
        /* Unknown gate */
        std::free(rho);
        write(write_fd, &result, sizeof(result));
        close(write_fd);
        _exit(1);
    }

    /* Warm-up: apply gate to qubit 0 */
    for (int i = 0; i < warmup; ++i) {
        gate_fn(0, rho, dim);
    }

    /* Timed iterations: apply gate to each qubit */
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        for (qubit_t t = 0; t < qubits; ++t) {
            gate_fn(static_cast<UINT>(t), rho, dim);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    result.time_ms = static_cast<double>(duration.count()) / 1000000.0;
    result.ops_per_sec = static_cast<double>(iterations * qubits) / (result.time_ms / 1000.0);

    std::free(rho);

    /* Write result to parent */
    write(write_fd, &result, sizeof(result));
    close(write_fd);
    _exit(0);
}

extern "C" {

/**
 * @brief Benchmark Qulacs density matrix gate operations
 *
 * Forks a child process for each benchmark to ensure accurate memory
 * measurement, avoiding memory reuse between benchmarks.
 *
 * @param qubits Number of qubits
 * @param gate_name Name of the gate ("X", "H", "Z")
 * @param iterations Number of iterations
 * @param warmup Number of warm-up iterations
 * @return bench_result_t Benchmark results
 */
bench_result_t bench_qulacs(qubit_t qubits, const char *gate_name,
                            int iterations, int warmup) {
    bench_result_t result = {0};
    result.qubits = qubits;
    result.dim = (dim_t)1 << qubits;
    result.gate_name = gate_name;
    result.method = "qulacs";
    result.iterations = iterations;

    /* Create pipe for child to send results back */
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
        /* Child process */
        close(pipefd[0]);  /* Close read end */
        bench_qulacs_child(pipefd[1], qubits, gate_name, iterations, warmup);
        _exit(0);  /* Should not reach here */
    }

    /* Parent process */
    close(pipefd[1]);  /* Close write end */

    /* Read result from child */
    bench_result_t child_result;
    ssize_t n = read(pipefd[0], &child_result, sizeof(child_result));
    close(pipefd[0]);

    /* Wait for child to finish */
    int status;
    waitpid(pid, &status, 0);

    if (n == static_cast<ssize_t>(sizeof(child_result))) {
        result.time_ms = child_result.time_ms;
        result.ops_per_sec = child_result.ops_per_sec;
        result.memory_bytes = child_result.memory_bytes;
    }

    return result;
}

} /* extern "C" */

#endif /* WITH_QULACS */
