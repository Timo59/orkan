/**
 * @file bench_qpp.cpp
 * @brief Quantum++ framework benchmark integration
 *
 * Quantum++ is a modern C++17 general-purpose quantum computing library.
 * It uses Eigen3 for linear algebra and supports density matrix simulation.
 *
 * Homepage: https://github.com/softwareQinc/qpp
 *
 * To enable: cmake -DWITH_QPP=ON ..
 * Requires: Eigen3 library
 *
 * Memory measurement: We fork a child process for each benchmark to get
 * accurate memory measurements, avoiding memory reuse between benchmarks.
 */

#ifdef WITH_QPP

#include <qpp/qpp.hpp>
#include <chrono>
#include <unistd.h>
#include <sys/wait.h>

extern "C" {
#include "bench.h"
}

using namespace qpp;

/*
 * =====================================================================================================================
 * Quantum++ benchmark implementation
 * =====================================================================================================================
 */

/**
 * @brief Run Quantum++ benchmark in child process (internal)
 *
 * This runs in a forked child and writes results to the pipe.
 */
static void bench_qpp_child(int write_fd, qubit_t qubits, const char *gate_name,
                            int iterations, int warmup) {
    bench_result_t result = {0};
    result.qubits = qubits;
    result.dim = (dim_t)1 << qubits;
    result.iterations = iterations;

    /* Measure actual memory usage in fresh process */
    size_t mem_before = bench_get_rss();

    /* Create list of qubit dimensions (all qubits have dimension 2) */
    std::vector<idx> dims(qubits, 2);

    /*
     * Create initial density matrix |0...0><0...0|.
     * Use Ones() first to ensure all memory pages are touched by the OS,
     * avoiding lazy allocation that would cause RSS to underreport memory.
     */
    idx dim = 1ULL << qubits;
    cmat rho = cmat::Ones(dim, dim);
    rho.setZero();
    rho(0, 0) = 1.0;

    size_t mem_after = bench_get_rss();
    result.memory_bytes = (mem_after > mem_before) ? (mem_after - mem_before) : 0;

    /* Select gate matrix based on name */
    cmat U;
    if (gate_name[0] == 'X') {
        U = Gates::get_instance().X;
    } else if (gate_name[0] == 'H') {
        U = Gates::get_instance().H;
    } else if (gate_name[0] == 'Z') {
        U = Gates::get_instance().Z;
    } else {
        /* Write empty result and exit */
        write(write_fd, &result, sizeof(result));
        close(write_fd);
        _exit(1);
    }

    /* Warm-up: apply U*rho*U† for each target qubit */
    for (int i = 0; i < warmup; ++i) {
        rho = apply(rho, U, {0}, dims);
    }

    /* Timed iterations */
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        for (qubit_t t = 0; t < qubits; ++t) {
            rho = apply(rho, U, {static_cast<idx>(t)}, dims);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    result.time_ms = static_cast<double>(duration.count()) / 1000000.0;
    result.ops_per_sec = static_cast<double>(iterations * qubits) / (result.time_ms / 1000.0);

    /* Write result to parent */
    write(write_fd, &result, sizeof(result));
    close(write_fd);
    _exit(0);
}

extern "C" {

/**
 * @brief Benchmark Quantum++ density matrix gate operations
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
bench_result_t bench_qpp(qubit_t qubits, const char *gate_name,
                         int iterations, int warmup) {
    bench_result_t result = {0};
    result.qubits = qubits;
    result.dim = (dim_t)1 << qubits;
    result.gate_name = gate_name;
    result.method = "qpp";
    result.iterations = iterations;

    /* Create pipe for child to send results back */
    int pipefd[2];
    if (pipe(pipefd) == -1) {
        perror("bench_qpp: pipe failed");
        return result;
    }

    pid_t pid = fork();
    if (pid == -1) {
        perror("bench_qpp: fork failed");
        close(pipefd[0]);
        close(pipefd[1]);
        return result;
    }

    if (pid == 0) {
        /* Child process */
        close(pipefd[0]);  /* Close read end */
        bench_qpp_child(pipefd[1], qubits, gate_name, iterations, warmup);
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

#endif /* WITH_QPP */
