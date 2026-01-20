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
 */

#ifdef WITH_QUEST

#include "bench.h"
#include "quest/include/quest.h"
#include <stdio.h>
#include <stdlib.h>

static int quest_initialized = 0;

void bench_quest_init(void) {
    if (!quest_initialized) {
        initQuESTEnv();
        quest_initialized = 1;
    }
}

void bench_quest_cleanup(void) {
    if (quest_initialized) {
        finalizeQuESTEnv();
        quest_initialized = 0;
    }
}

/**
 * @brief Benchmark QuEST density matrix gate operations
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
    /* QuEST stores full density matrix internally */
    result.memory_bytes = bench_dense_size(qubits);
    result.iterations = iterations;

    /* Ensure QuEST is initialized */
    bench_quest_init();

    /* Create density matrix */
    Qureg rho = createDensityQureg((int)qubits);
    initZeroState(rho);

    /* Select gate function based on name */
    void (*gate_fn)(Qureg, int) = NULL;
    if (gate_name[0] == 'X') {
        gate_fn = applyPauliX;
    } else if (gate_name[0] == 'H') {
        gate_fn = applyHadamard;
    } else if (gate_name[0] == 'Z') {
        gate_fn = applyPauliZ;
    } else {
        fprintf(stderr, "bench_quest: unknown gate '%s'\n", gate_name);
        destroyQureg(rho);
        return result;
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

    /* Cleanup Qureg only (not the environment) */
    destroyQureg(rho);

    return result;
}

#endif /* WITH_QUEST */
