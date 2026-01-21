/**
 * @file bench_quest.c
 * @brief QuEST framework benchmark integration
 *
 * QuEST (Quantum Exact Simulation Toolkit) is a high-performance,
 * multithreaded, distributed quantum simulator with density matrix support.
 *
 * Homepage: https://quest.qtechtheory.org/
 * GitHub: https://github.com/QuEST-Kit/QuEST
 *
 * To enable: cmake -DWITH_QUEST=ON -DQuEST_DIR=/path/to/quest ..
 *
 * Memory measurement: QuEST cannot be re-initialized after finalization,
 * so we fork a child process for each benchmark to get accurate memory
 * measurements that include the full environment setup cost.
 */

#ifdef WITH_QUEST

#include "bench.h"
#include "quest/include/quest.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>

void bench_quest_init(void) {
    /* No-op: environment initialized in forked child */
}

void bench_quest_cleanup(void) {
    /* No-op: environment finalized in forked child */
}

/**
 * @brief Run QuEST benchmark in child process (internal)
 *
 * This runs in a forked child and writes results to the pipe.
 */
static void bench_quest_child(int write_fd, qubit_t qubits, const char *gate_name,
                               int iterations, int warmup) {
    bench_result_t result = {0};
    result.qubits = qubits;
    result.dim = (dim_t)1 << qubits;
    result.iterations = iterations;

    /*
     * Measure total memory: environment + Qureg allocation.
     * Fresh process ensures we capture the true memory cost.
     */
    size_t mem_before = bench_get_rss();

    initQuESTEnv();
    Qureg rho = createDensityQureg((int)qubits);

    /*
     * Use initDebugState() to ensure all memory pages are touched.
     * initZeroState() may leave pages unmapped due to lazy allocation.
     */
    initDebugState(rho);

    size_t mem_after = bench_get_rss();
    result.memory_bytes = (mem_after > mem_before) ? (mem_after - mem_before) : 0;

    /* Select gate function based on name */
    void (*gate_fn)(Qureg, int) = NULL;
    if (gate_name[0] == 'X') {
        gate_fn = applyPauliX;
    } else if (gate_name[0] == 'H') {
        gate_fn = applyHadamard;
    } else if (gate_name[0] == 'Z') {
        gate_fn = applyPauliZ;
    } else {
        /* Write empty result and exit */
        write(write_fd, &result, sizeof(result));
        close(write_fd);
        destroyQureg(rho);
        finalizeQuESTEnv();
        _exit(1);
    }

    /* Warm-up */
    for (int i = 0; i < warmup; ++i) {
        gate_fn(rho, 0);
    }

    /* Timed iterations */
    uint64_t start = bench_time_ns();
    for (int i = 0; i < iterations; ++i) {
        for (qubit_t t = 0; t < qubits; ++t) {
            gate_fn(rho, (int)t);
        }
    }
    uint64_t end = bench_time_ns();

    result.time_ms = bench_ns_to_ms(end - start);
    result.ops_per_sec = (double)(iterations * qubits) / (result.time_ms / 1000.0);

    /* Write result to parent */
    write(write_fd, &result, sizeof(result));
    close(write_fd);

    /* Cleanup and exit */
    destroyQureg(rho);
    finalizeQuESTEnv();
    _exit(0);
}

/**
 * @brief Benchmark QuEST density matrix gate operations
 *
 * Forks a child process for each benchmark to ensure accurate memory
 * measurement. QuEST cannot be re-initialized after finalization.
 *
 * @param qubits Number of qubits
 * @param gate_name Name of the gate ("X", "H", "Z")
 * @param iterations Number of iterations
 * @param warmup Number of warm-up iterations
 * @return bench_result_t Benchmark results
 */
bench_result_t bench_quest(qubit_t qubits, const char *gate_name,
                           int iterations, int warmup) {
    bench_result_t result = {0};
    result.qubits = qubits;
    result.dim = (dim_t)1 << qubits;
    result.gate_name = gate_name;
    result.method = "quest";
    result.iterations = iterations;

    /* Create pipe for child to send results back */
    int pipefd[2];
    if (pipe(pipefd) == -1) {
        perror("bench_quest: pipe failed");
        return result;
    }

    pid_t pid = fork();
    if (pid == -1) {
        perror("bench_quest: fork failed");
        close(pipefd[0]);
        close(pipefd[1]);
        return result;
    }

    if (pid == 0) {
        /* Child process */
        close(pipefd[0]);  /* Close read end */
        bench_quest_child(pipefd[1], qubits, gate_name, iterations, warmup);
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

    if (n == sizeof(child_result)) {
        result.time_ms = child_result.time_ms;
        result.ops_per_sec = child_result.ops_per_sec;
        result.memory_bytes = child_result.memory_bytes;
    }

    return result;
}

#endif /* WITH_QUEST */
