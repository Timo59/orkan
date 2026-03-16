/**
 * @file bench_csv.c
 * @brief CSV output formatter for the gate benchmark (stub — Phase 3).
 */

#include "bench.h"

void bench_print_csv_header(const bench_gate_info_t *gate, int per_qubit) {
    (void)gate; (void)per_qubit;
}

void bench_print_csv_agg(qubit_t qubits,
                         const bench_result_t results[BACKEND_COUNT]) {
    (void)qubits; (void)results;
}

void bench_print_csv_perq(qubit_t qubits, bench_backend_id_t backend,
                          const bench_pos_result_t *results, int n_pos,
                          int gate_qubits) {
    (void)qubits; (void)backend; (void)results; (void)n_pos; (void)gate_qubits;
}
