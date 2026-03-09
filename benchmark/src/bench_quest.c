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
 *
 * Statistics: the child performs `runs` independent timed rounds and returns
 * the full bench_run_stats_t summary through the pipe.
 */

#ifdef WITH_QUEST

#include "bench.h"
#include "quest/include/quest.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>

extern int build_all_pairs(qubit_t qubits, qubit_t *q1_out, qubit_t *q2_out);

void bench_quest_init(void) {
    /* No-op: environment initialized in forked child */
}

void bench_quest_cleanup(void) {
    /* No-op: environment finalized in forked child */
}

/**
 * @brief Run QuEST benchmark in child process (internal)
 *
 * Performs `runs` independent timed rounds and fills all statistical fields
 * in result before writing it to the pipe.
 */
static void bench_quest_child(int write_fd, qubit_t qubits, const char *gate_name,
                               int iterations, int warmup, int runs) {
    bench_result_t result = {0};
    result.qubits     = qubits;
    result.dim        = (dim_t)1 << qubits;
    result.iterations = iterations;
    result.runs       = runs;

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

    /* Dispatch: single-qubit gates */
    void (*gate_fn_1q)(Qureg, int) = NULL;
    if (strcmp(gate_name, "X") == 0) {
        gate_fn_1q = applyPauliX;
    } else if (strcmp(gate_name, "H") == 0) {
        gate_fn_1q = applyHadamard;
    } else if (strcmp(gate_name, "Z") == 0) {
        gate_fn_1q = applyPauliZ;
    }

    /* Dispatch: two-qubit gates */
    void (*gate_fn_2q)(Qureg, int, int) = NULL;
    if (strcmp(gate_name, "CX") == 0) {
        gate_fn_2q = applyControlledPauliX;
    } else if (strcmp(gate_name, "SWAP") == 0) {
        gate_fn_2q = applySwap;
    }

    if (!gate_fn_1q && !gate_fn_2q) {
        /* Unknown gate - write empty (zeroed) pipe result and exit */
        bench_pipe_result_t pr = {0};
        if (write(write_fd, &pr, sizeof(pr)) != (ssize_t)sizeof(pr)) {
            close(write_fd);
            destroyQureg(rho);
            finalizeQuESTEnv();
            _exit(1);
        }
        close(write_fd);
        destroyQureg(rho);
        finalizeQuESTEnv();
        _exit(1);
    }

    double run_times[200];  /* bounded by --runs validation (max 200) */
    int cnt = (runs < 200) ? runs : 200;

    if (gate_fn_1q) {
        /* Single-qubit gate: `runs` independent timed rounds */
        for (int i = 0; i < warmup; ++i)
            gate_fn_1q(rho, (int)(i % qubits));

        for (int r = 0; r < cnt; ++r) {
            uint64_t start = bench_time_ns();
            for (int i = 0; i < iterations; ++i)
                for (qubit_t t = 0; t < qubits; ++t)
                    gate_fn_1q(rho, (int)t);
            run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
        }

        bench_run_stats_t stats = bench_compute_stats(run_times, cnt);
        result.time_ms        = stats.mean;
        result.time_ms_std    = stats.std_dev;
        result.time_ms_min    = stats.min;
        result.time_ms_median = stats.median;
        result.time_ms_cv     = stats.cv;
        result.ops_per_sec    = (double)(iterations * qubits) / (stats.mean / 1000.0);
    } else {
        /* Two-qubit gate: all ordered pairs (q1 < q2) */
        qubit_t q1s[128], q2s[128];
        int num_pairs = build_all_pairs(qubits, q1s, q2s);

        if (num_pairs == 0) {
            bench_pipe_result_t pr = {0};
            if (write(write_fd, &pr, sizeof(pr)) != (ssize_t)sizeof(pr)) {
                close(write_fd);
                destroyQureg(rho);
                finalizeQuESTEnv();
                _exit(1);
            }
            close(write_fd);
            destroyQureg(rho);
            finalizeQuESTEnv();
            _exit(0);
        }

        for (int i = 0; i < warmup; ++i)
            gate_fn_2q(rho, (int)q1s[i % num_pairs], (int)q2s[i % num_pairs]);

        for (int r = 0; r < cnt; ++r) {
            uint64_t start = bench_time_ns();
            for (int i = 0; i < iterations; ++i)
                for (int p = 0; p < num_pairs; ++p)
                    gate_fn_2q(rho, (int)q1s[p], (int)q2s[p]);
            run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
        }

        bench_run_stats_t stats = bench_compute_stats(run_times, cnt);
        result.time_ms        = stats.mean;
        result.time_ms_std    = stats.std_dev;
        result.time_ms_min    = stats.min;
        result.time_ms_median = stats.median;
        result.time_ms_cv     = stats.cv;
        result.ops_per_sec    = (double)(iterations * num_pairs) / (stats.mean / 1000.0);
    }

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
    if (write(write_fd, &pr, sizeof(pr)) != (ssize_t)sizeof(pr)) {
        close(write_fd);
        destroyQureg(rho);
        finalizeQuESTEnv();
        _exit(1);
    }
    close(write_fd);

    destroyQureg(rho);
    finalizeQuESTEnv();
    _exit(0);
}

/**
 * @brief Benchmark QuEST density matrix gate operations.
 *
 * Forks a child process; the child performs `runs` independent timed rounds
 * and returns the full statistical summary via pipe.
 */
bench_result_t bench_quest(qubit_t qubits, const char *gate_name,
                           int iterations, int warmup, int runs) {
    bench_result_t result = {0};
    result.qubits     = qubits;
    result.dim        = (dim_t)1 << qubits;
    result.gate_name  = gate_name;
    result.method     = "quest";
    result.iterations = iterations;
    result.runs       = runs;

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
        close(pipefd[0]);
        bench_quest_child(pipefd[1], qubits, gate_name, iterations, warmup, runs);
        _exit(0);  /* unreachable */
    }

    close(pipefd[1]);

    bench_pipe_result_t pr;
    ssize_t n = read(pipefd[0], &pr, sizeof(pr));
    close(pipefd[0]);

    int status;
    waitpid(pid, &status, 0);

    if (n == (ssize_t)sizeof(pr)) {
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

#endif /* WITH_QUEST */
