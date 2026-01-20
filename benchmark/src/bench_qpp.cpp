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
 */

#ifdef WITH_QPP

#include <qpp/qpp.hpp>
#include <chrono>

extern "C" {
#include "bench.h"
}

using namespace qpp;

/*
 * =====================================================================================================================
 * Quantum++ benchmark implementation
 * =====================================================================================================================
 */

extern "C" {

/**
 * @brief Benchmark Quantum++ density matrix gate operations
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
    /* Quantum++ stores full density matrix */
    result.memory_bytes = bench_dense_size(qubits);
    result.iterations = iterations;

    /* Create list of qubit dimensions (all qubits have dimension 2) */
    std::vector<idx> dims(qubits, 2);

    /* Create initial density matrix |0...0><0...0| */
    idx dim = 1ULL << qubits;
    cmat rho = cmat::Zero(dim, dim);
    rho(0, 0) = 1.0;

    /* Select gate matrix based on name */
    cmat U;
    if (gate_name[0] == 'X') {
        U = Gates::get_instance().X;
    } else if (gate_name[0] == 'H') {
        U = Gates::get_instance().H;
    } else if (gate_name[0] == 'Z') {
        U = Gates::get_instance().Z;
    } else {
        fprintf(stderr, "bench_qpp: unknown gate '%s'\n", gate_name);
        return result;
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

    return result;
}

} /* extern "C" */

#endif /* WITH_QPP */
